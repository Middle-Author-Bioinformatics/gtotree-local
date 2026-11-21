#!/usr/bin/env python3
"""
GToTree S3 worker for the Amplify React app.

This is the GToTree counterpart of the FeGenie worker: same S3 polling model,
same status.json state machine, same lock / --once cron design. It

1. Scans an input S3 bucket for uploaded GToTree jobs:
     s3://<input-bucket>/<slug>/form-data.txt
2. Downloads the job folder to a local work directory.
3. Parses form-data.txt into a job manifest (SCG set, input files, flags).
4. Builds the input list files GToTree needs (-A / -g / -f) and the command.
5. Runs GToTree.
6. Normalizes the primary outputs to stable, frontend-known names.
7. Bakes a self-contained interactive HTML dashboard (gtotree_report.py).
8. Publishes frontend-ready results to:
     s3://<results-bucket>/<slug>/

Frontend-required result names (published at <slug>/output/ and <slug>/)
-----------------------------------------------------------------------
  - output/GToTree_output.tre           (final Newick tree)
  - output/Aligned_SCGs.faa             (concatenated SCG alignment)
  - output/Genomes_summary_info.tsv     (per-genome QC summary)
  - output/SCG_hit_counts.tsv           (per-genome SCG hit matrix)
  - output/gtotree-runlog.txt           (run log)
  - gtotree-report.html                 (self-contained interactive dashboard)
  - status.json                         (job state)
  - result.json                         (viewer manifest)
  - raw-results.tar.gz                  (full GToTree output dir)

Typical cron usage
------------------
  /usr/bin/python3 /opt/gtotree/gtotree_worker.py --once \
    --input-bucket midauthorbio-gtotree-input \
    --results-bucket midauthorbio-gtotree-results \
    --work-root /data/gtotree-worker \
    --command-prefix "conda run -n gtotree"

If GToTree is already on PATH in the active environment, omit --command-prefix.
"""

from __future__ import annotations

import argparse
import json
import fcntl
import signal
import threading
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from job_validation import validate_manifest, bound_resource_flags

import boto3
from botocore.config import Config
from botocore.exceptions import ClientError


FASTA_EXTS = {".fa", ".fasta", ".fna"}
PROTEIN_EXTS = {".faa"}
GENBANK_EXTS = {".gb", ".gbk", ".gbff", ".gbf"}
SKIP_INPUT_NAMES = {"form-data.txt"}
SLUG_RE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")

# GToTree's advanced flags are user-controlled text spliced from the web app.
# Re-validate every token here (same allow-list philosophy as the frontend /
# gtotree_processing.py). value: number of value tokens the flag consumes.
ALLOWED_GTOTREE_FLAGS = {
    "-t": 0, "-D": 0, "-B": 0, "-N": 0, "-k": 0, "-X": 0, "-P": 0,
    "-L": 1, "-c": 1, "-G": 1, "-T": 1, "-n": 1, "-M": 1, "-j": 1,
}
_TREE_PROGRAMS = {"FastTree", "FastTreeMP", "IQ-TREE"}


@dataclass
class JobManifest:
    slug: str
    job_name: str = "GToTree_output"
    gene_set: str = ""
    accessions: str | None = None
    metadata_file: str | None = None
    proteins: list[str] = field(default_factory=list)
    genbanks: list[str] = field(default_factory=list)
    genomes: list[str] = field(default_factory=list)
    ko_file: str | None = None
    pfam_file: str | None = None
    extra_flags: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def configure_logging(log_path: Path, verbose: bool = False) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handlers: list[logging.Handler] = [
        logging.FileHandler(log_path),
        logging.StreamHandler(sys.stdout),
    ]
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=handlers,
    )


def make_s3(region: str):
    # Regional endpoint avoids browser presigned-PUT CORS failures on redirect.
    return boto3.client(
        "s3",
        region_name=region,
        endpoint_url=f"https://s3.{region}.amazonaws.com",
        config=Config(signature_version="s3v4", s3={"addressing_style": "virtual"}),
    )


def s3_key_exists(s3, bucket: str, key: str) -> bool:
    try:
        s3.head_object(Bucket=bucket, Key=key)
        return True
    except ClientError as e:
        if e.response.get("Error", {}).get("Code") in {"404", "NoSuchKey", "NotFound"}:
            return False
        raise


def put_json(s3, bucket: str, key: str, payload: dict) -> None:
    s3.put_object(
        Bucket=bucket,
        Key=key,
        Body=json.dumps(payload, indent=2, sort_keys=True).encode("utf-8"),
        ContentType="application/json",
        ContentDisposition="inline",
    )


def upload_file(s3, bucket: str, key: str, path: Path, content_type: str | None = None) -> None:
    extra = {"ContentDisposition": "inline"}
    if content_type:
        extra["ContentType"] = content_type
    s3.upload_file(str(path), bucket, key, ExtraArgs=extra)


# ---------------------------------------------------------------------------
# job discovery
# ---------------------------------------------------------------------------
def list_job_slugs(s3, input_bucket: str) -> list[str]:
    """Find run folders by locating <slug>/form-data.txt at the bucket root."""
    paginator = s3.get_paginator("list_objects_v2")
    slugs: set[str] = set()
    for page in paginator.paginate(Bucket=input_bucket):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith("/form-data.txt"):
                parts = key.split("/")
                if len(parts) == 2 and SLUG_RE.match(parts[0]):
                    slugs.add(parts[0])
    return sorted(slugs)


def download_prefix(s3, bucket: str, prefix: str, dest: Path, max_bytes: int = 100 * 1024 * 1024) -> list[Path]:
    dest.mkdir(parents=True, exist_ok=True)
    downloaded: list[Path] = []
    paginator = s3.get_paginator("list_objects_v2")
    objects, total_bytes = [], 0
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            objects.append(obj)
            total_bytes += obj.get("Size", 0)
            if len(objects) > 250 or total_bytes > max_bytes:
                raise ValueError("Submission exceeds file-count or upload-size limits")
    for obj in objects:
        key = obj["Key"]
        if key.endswith("/"):
            continue
        rel = Path(key).relative_to(prefix)
        if len(rel.parts) != 1 or rel.name in (".", "..") or "\\" in rel.name:
            raise ValueError("Job objects must use flat, safe input filenames")
        local = (dest / rel).resolve()
        if local.parent != dest.resolve():
            raise ValueError("Input path escapes the job directory")
        local.parent.mkdir(parents=True, exist_ok=True)
        s3.download_file(bucket, key, str(local))
        downloaded.append(local)
        logging.info("Downloaded s3://%s/%s -> %s", bucket, key, local)
    return downloaded


# ---------------------------------------------------------------------------
# manifest parsing (matches the web app's buildFormData in gtotree-config.ts)
# ---------------------------------------------------------------------------
def _kv(line: str):
    if ":" not in line:
        return None, None
    k, v = line.split(":", 1)
    return k.strip(), v.strip()


def sanitize_flags(raw: str) -> list[str]:
    tokens = shlex.split(raw) if raw else []
    clean, i = [], 0
    while i < len(tokens):
        tok = tokens[i]
        if tok not in ALLOWED_GTOTREE_FLAGS:
            logging.warning("[flags] dropping unrecognized token: %r", tok)
            i += 1
            continue
        if ALLOWED_GTOTREE_FLAGS[tok]:
            if i + 1 >= len(tokens):
                logging.warning("[flags] %s expects a value; dropping", tok)
                break
            val = tokens[i + 1]
            if tok == "-T" and val not in _TREE_PROGRAMS:
                logging.warning("[flags] bad tree program %r; dropping", val)
                i += 2
                continue
            if tok in ("-c", "-G") and not re.match(r"^\d*\.?\d+$", val):
                logging.warning("[flags] %s non-numeric %r; dropping", tok, val)
                i += 2
                continue
            if tok == "-L" and not re.match(r"^[A-Za-z,]+$", val):
                logging.warning("[flags] -L bad ranks %r; dropping", val)
                i += 2
                continue
            if tok in ("-n", "-M", "-j") and not re.match(r"^\d+$", val):
                logging.warning("[flags] %s non-integer %r; dropping", tok, val)
                i += 2
                continue
            clean.extend([tok, val])
            i += 2
        else:
            clean.append(tok)
            i += 1
    return clean


def parse_manifest(path: Path, fallback_slug: str) -> JobManifest:
    m = JobManifest(slug=fallback_slug)
    if not path.exists():
        return m
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            key, val = _kv(line.rstrip("\n"))
            if key is None:
                continue
            if key == "Name":
                m.job_name = val or m.job_name
            elif key == "SCG set":
                m.gene_set = val
            elif key == "Accessions File":
                m.accessions = None if val in ("", "N/A") else val
            elif key == "Genome Metadata File":
                m.metadata_file = None if val in ("", "N/A") else val
            elif key == "Protein File":
                m.proteins.append(val)
            elif key == "GenBank File":
                m.genbanks.append(val)
            elif key == "Genome File":
                m.genomes.append(val)
            elif key == "KOFile":
                m.ko_file = None if val in ("", "N/A") else val
            elif key == "PfamFile":
                m.pfam_file = None if val in ("", "N/A") else val
            elif key == "GtotreeFlags":
                m.extra_flags = sanitize_flags(val)
    return m


# ---------------------------------------------------------------------------
# command construction
# ---------------------------------------------------------------------------
def write_list_file(job_dir: Path, filenames: list[str], list_name: str) -> Path:
    path = job_dir / list_name
    with path.open("w") as fh:
        for n in filenames:
            fh.write(str((job_dir / n).resolve()) + "\n")
    return path


def build_gtotree_command(args, m: JobManifest, job_dir: Path, out_dir: Path) -> list[str]:
    cmd: list[str] = []
    if args.command_prefix:
        cmd.extend(shlex.split(args.command_prefix))
    cmd.append(args.gtotree_bin)

    if m.accessions:
        cmd += ["-a", str((job_dir / m.accessions).resolve())]
    if m.proteins:
        cmd += ["-A", str(write_list_file(job_dir, m.proteins, "protein_files.txt"))]
    if m.genbanks:
        cmd += ["-g", str(write_list_file(job_dir, m.genbanks, "genbank_files.txt"))]
    if m.genomes:
        cmd += ["-f", str(write_list_file(job_dir, m.genomes, "fasta_files.txt"))]

    cmd += ["-H", m.gene_set]  # required

    if m.ko_file:
        cmd += ["-K", str((job_dir / m.ko_file).resolve())]
    if m.pfam_file:
        cmd += ["-p", str((job_dir / m.pfam_file).resolve())]

    cmd += bound_resource_flags(m.extra_flags, args.max_threads)
    cmd += ["-o", str(out_dir), "-F"]
    return cmd


def run_command(cmd: list[str], log_handle, timeout: int = 21600) -> None:
    logging.info("Running command: %s", shlex.join(cmd))
    proc = subprocess.Popen(cmd, stdout=log_handle, stderr=subprocess.STDOUT, text=True, start_new_session=True)
    try:
        result = proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
        except ProcessLookupError:
            pass
        raise RuntimeError("Analysis exceeded its %s-second runtime limit" % timeout)
    if result != 0:
        raise RuntimeError("GToTree/report command failed with exit code %s" % result)


# ---------------------------------------------------------------------------
# output collection
# ---------------------------------------------------------------------------
def find_first(out_dir: Path, patterns: Iterable[str]) -> Path | None:
    """Return the first file whose name matches any regex (case-insensitive)."""
    compiled = [re.compile(p, re.I) for p in patterns]
    paths = sorted(path for path in out_dir.rglob("*") if path.is_file())
    for pattern in compiled:
        for path in paths:
            if pattern.search(path.name):
                return path
    return None


def collect_outputs(out_dir: Path) -> dict[str, Path | None]:
    """Locate GToTree's primary outputs regardless of the -o dir name."""
    return {
        # Prefer the *_mod_names variant (readable taxonomy labels from -t/-D)
        # so tips carry names by default; fall back to the accession-labeled tree.
        "tree": find_first(out_dir, [r"mod_names.*\.tre$", r"\.tre$", r"\.tree$", r"\.nwk$", r"\.newick$"]),
        "alignment": find_first(out_dir, [r"mod_names.*\.faa$", r"^Aligned_SCGs.*\.faa$", r"aligned.*\.faa$"]),
        "summary": find_first(out_dir, [r"Genomes?_summary_info\.tsv$", r"summary_info\.tsv$"]),
        "hits": find_first(out_dir, [r"SCG_hit_counts\.tsv$"]),
        "runlog": find_first(out_dir, [r"gtotree-?runlog\.txt$", r"runlog\.txt$"]),
        "partitions": find_first(out_dir, [r"Partitions\.txt$"]),
        "citations": find_first(out_dir, [r"citations\.txt$"]),
        "itol_colors": find_first(out_dir, [r"iToL-colors\.txt$"]),
    }


def make_tarball(source_dir: Path, tar_gz: Path) -> None:
    with tarfile.open(tar_gz, "w:gz") as tar:
        tar.add(source_dir, arcname=source_dir.name)


def generate_report(args, report_dest: Path, outputs: dict, manifest: JobManifest, log_handle) -> bool:
    """Bake the self-contained interactive HTML dashboard via gtotree_report.py."""
    candidates = []
    if args.report_script:
        candidates.append(Path(args.report_script))
    candidates.append(Path(__file__).resolve().parent / "gtotree_report.py")
    candidates.append(Path.cwd() / "gtotree_report.py")

    seen: set[Path] = set()
    for script in candidates:
        script = script.expanduser().resolve()
        if script in seen or not script.exists():
            continue
        seen.add(script)
        cmd: list[str] = []
        if args.command_prefix:
            cmd.extend(shlex.split(args.command_prefix))
        cmd += ["python", str(script), "-o", str(report_dest), "--title", manifest.job_name,
                "--gene-set", manifest.gene_set, "--report-code", manifest.slug,
                "--contact-email", args.contact_email]
        if manifest.metadata_file:
            cmd += ["--metadata", str(report_dest.parents[1] / "input" / manifest.metadata_file)]
        if outputs.get("tree"):
            cmd += ["--tree", str(outputs["tree"])]
        if outputs.get("alignment"):
            cmd += ["--alignment", str(outputs["alignment"])]
        if outputs.get("summary"):
            cmd += ["--summary", str(outputs["summary"])]
        if outputs.get("hits"):
            cmd += ["--hits", str(outputs["hits"])]
        logging.info("Generating GToTree report with: %s", script)
        try:
            run_command(cmd, log_handle, min(args.job_timeout, 600))
        except Exception:
            logging.exception("Report generation failed with %s", script)
            continue
        if report_dest.exists():
            return True
    logging.warning("gtotree_report.py not found or failed; no HTML dashboard produced.")
    return False


# ---------------------------------------------------------------------------
# per-job pipeline
# ---------------------------------------------------------------------------
def process_job(args, s3, slug: str) -> None:
    result_prefix = f"{slug}/"
    status_key = f"{result_prefix}status.json"

    if not SLUG_RE.fullmatch(slug):
        raise ValueError("Invalid job code")
    if s3_key_exists(s3, args.results_bucket, status_key) and not args.force:
        existing = json.loads(s3.get_object(Bucket=args.results_bucket, Key=status_key)["Body"].read())
        if not (args.retry_failed and existing.get("state") in ("failed", "running")):
            logging.info("Skipping %s; result status already exists.", slug)
            return

    input_prefix = f"{slug}/"
    if not s3_key_exists(s3, args.input_bucket, f"{input_prefix}form-data.txt"):
        logging.info("Skipping %s; no form-data.txt yet.", slug)
        return

    job_root = Path(args.work_root) / slug
    if job_root.exists():
        shutil.rmtree(job_root)
    job_dir = job_root / "input"
    out_dir = job_root / "GToTree_output"
    final_dir = job_root / "frontend_results" / "output"
    log_path = job_root / "run.log"
    job_root.mkdir(parents=True, exist_ok=True)

    put_json(s3, args.results_bucket, status_key,
             {"slug": slug, "state": "running", "started_at": utc_now(), "input_prefix": input_prefix})

    success = False
    manifest = JobManifest(slug=slug)
    try:
        download_prefix(s3, args.input_bucket, input_prefix, job_dir, args.max_upload_bytes)
        manifest = parse_manifest(job_dir / "form-data.txt", fallback_slug=slug)
        if not manifest.gene_set:
            raise RuntimeError("No 'SCG set' in form-data.txt (required for -H).")
        if not (manifest.accessions or manifest.proteins or manifest.genbanks or manifest.genomes):
            raise RuntimeError("No genome inputs (accessions / protein / GenBank / genome).")

        input_summary = validate_manifest(manifest, job_dir, args)
        final_dir.mkdir(parents=True, exist_ok=True)
        with log_path.open("w", encoding="utf-8") as log_handle:
            log_handle.write(f"GToTree worker started {utc_now()}\nSlug: {slug}\n")
            cmd = build_gtotree_command(args, manifest, job_dir, out_dir)
            run_command(cmd, log_handle, args.job_timeout)

            outputs = collect_outputs(out_dir)
            if not (outputs["tree"] or outputs["alignment"]):
                raise RuntimeError("GToTree produced neither a tree nor an alignment.")

            if manifest.metadata_file:
                shutil.copy2(job_dir / manifest.metadata_file, final_dir / "selected-genomes.json")
            # normalize to stable, frontend-known names
            stable = {}
            if outputs["tree"]:
                dest = final_dir / "GToTree_output.tre"
                shutil.copy2(outputs["tree"], dest); stable["tree"] = dest
            if outputs["alignment"]:
                dest = final_dir / "Aligned_SCGs.faa"
                shutil.copy2(outputs["alignment"], dest); stable["alignment"] = dest
            if outputs["summary"]:
                dest = final_dir / "Genomes_summary_info.tsv"
                shutil.copy2(outputs["summary"], dest); stable["summary"] = dest
            if outputs["hits"]:
                dest = final_dir / "SCG_hit_counts.tsv"
                shutil.copy2(outputs["hits"], dest); stable["hits"] = dest
            if outputs["runlog"]:
                shutil.copy2(outputs["runlog"], final_dir / "gtotree-runlog.txt")
            if outputs["partitions"]:
                shutil.copy2(outputs["partitions"], final_dir / "Partitions.txt")
            if outputs["citations"]:
                shutil.copy2(outputs["citations"], final_dir / "citations.txt")
            if outputs["itol_colors"]:
                shutil.copy2(outputs["itol_colors"], final_dir / "iToL-colors.txt")

            report_dest = job_root / "frontend_results" / "gtotree-report.html"
            report_ok = generate_report(args, report_dest, stable, manifest, log_handle)
            if not report_ok:
                raise RuntimeError("Tree analysis finished, but the interactive report could not be generated. See run.log.")

        # tarball of the full raw output dir
        tar_path = job_root / "raw-results.tar.gz"
        make_tarball(out_dir, tar_path)

        # ---- upload everything ----
        content_types = {
            ".tre": "text/plain", ".faa": "text/plain", ".tsv": "text/tab-separated-values",
            ".txt": "text/plain", ".html": "text/html", ".json": "application/json",
        }
        uploaded_files = []
        for f in sorted(final_dir.rglob("*")):
            if f.is_file():
                key = f"{result_prefix}output/{f.relative_to(final_dir)}"
                upload_file(s3, args.results_bucket, key, f, content_types.get(f.suffix, "text/plain"))
                uploaded_files.append(f"output/{f.relative_to(final_dir)}")
        if report_ok and report_dest.exists():
            upload_file(s3, args.results_bucket, f"{result_prefix}gtotree-report.html", report_dest, "text/html")
            uploaded_files.append("gtotree-report.html")
        upload_file(s3, args.results_bucket, f"{result_prefix}raw-results.tar.gz", tar_path, "application/gzip")
        uploaded_files.append("raw-results.tar.gz")
        upload_file(s3, args.results_bucket, f"{result_prefix}run.log", log_path, "text/plain")
        uploaded_files.append("run.log")

        # viewer manifest (what the React ResultsView reads)
        put_json(s3, args.results_bucket, f"{result_prefix}result.json", {
            "slug": slug, "run_name": manifest.job_name, "gene_set": manifest.gene_set,
            "has_tree": bool(stable.get("tree")), "has_alignment": bool(stable.get("alignment")),
            "has_report": report_ok, "files": uploaded_files, "input_summary": input_summary,
        })
        # completion status
        put_json(s3, args.results_bucket, status_key, {
            "slug": slug, "state": "complete", "completed_at": utc_now(),
            "result_prefix": result_prefix, "has_report": report_ok, "files": uploaded_files,
        })
        success = True
        logging.info("Completed job %s -> s3://%s/%s", slug, args.results_bucket, result_prefix)

    except Exception as e:
        logging.exception("Job %s failed", slug)
        put_json(s3, args.results_bucket, status_key,
                 {"slug": slug, "state": "failed", "failed_at": utc_now(), "error": str(e)})
        if log_path.exists():
            upload_file(s3, args.results_bucket, f"{result_prefix}run.log", log_path, "text/plain")
        if not args.continue_on_error:
            raise
    finally:
        if args.clean and success and job_root.exists():
            shutil.rmtree(job_root)


# ---------------------------------------------------------------------------
# lock + main loop
# ---------------------------------------------------------------------------
def acquire_lock(lock_path: Path):
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    handle = lock_path.open("a+")
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        handle.close()
        raise FileExistsError(str(lock_path))
    handle.seek(0); handle.truncate(); handle.write(str(os.getpid())); handle.flush()
    return handle


def release_lock(handle):
    # Do not unlink: another process may be waiting on this same inode.
    fcntl.flock(handle, fcntl.LOCK_UN)
    handle.close()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Poll S3 for GToTree jobs, run GToTree, publish frontend-ready results.")
    p.add_argument("--input-bucket", default=os.getenv("GTOTREE_INPUT_BUCKET", "midauthorbio-gtotree-input"))
    p.add_argument("--results-bucket", default=os.getenv("GTOTREE_RESULTS_BUCKET", "midauthorbio-gtotree-results"))
    p.add_argument("--region", default=os.getenv("AWS_REGION", "us-east-2"))
    p.add_argument("--work-root", default=os.getenv("GTOTREE_WORK_ROOT", "/tmp/gtotree-worker"))
    p.add_argument("--gtotree-bin", default=os.getenv("GTOTREE_BIN", "GToTree"))
    p.add_argument("--command-prefix", default=os.getenv("GTOTREE_COMMAND_PREFIX", ""),
                   help='Optional prefix, e.g. "conda run -n gtotree"')
    p.add_argument("--report-script", default=os.getenv("GTOTREE_REPORT_SCRIPT", ""),
                   help="Optional path to gtotree_report.py.")
    p.add_argument("--once", action="store_true", help="Run one scan and exit (cron mode).")
    p.add_argument("--interval", type=int, default=300, help="Polling interval seconds when not using --once.")
    p.add_argument("--force", action="store_true", help="Reprocess jobs even if status.json exists.")
    p.add_argument("--clean", action="store_true", help="Delete local work dir after each job.")
    p.add_argument("--continue-on-error", action="store_true", default=True, help="Continue after a job failure (default).")
    p.add_argument("--stop-on-error", action="store_false", dest="continue_on_error", help="Exit after a job failure.")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--lock-file", default=os.getenv("GTOTREE_LOCK_FILE", "/tmp/gtotree-worker.lock"))
    p.add_argument("--log-file", default=os.getenv("GTOTREE_WORKER_LOG", "/tmp/gtotree-worker.log"))
    p.add_argument("--job", action="append", help="Process only this result code; repeat to select several jobs.")
    p.add_argument("--retry-failed", action="store_true", help="Retry selected failed/interrupted jobs; requires --job.")
    p.add_argument("--max-references", type=int, default=100)
    p.add_argument("--max-uploaded-genomes", type=int, default=100)
    p.add_argument("--max-total-genomes", type=int, default=200)
    p.add_argument("--max-upload-bytes", type=int, default=100 * 1024 * 1024)
    p.add_argument("--max-threads", type=int, default=16)
    p.add_argument("--job-timeout", type=int, default=21600, help="Maximum analysis runtime in seconds (default 6 hours).")
    p.add_argument("--contact-email", default=os.getenv("MAB_DASHBOARD_EMAIL", "ark@midauthorbio.com"))
    args = p.parse_args()
    if (args.force or args.retry_failed) and not args.job:
        p.error("--force and --retry-failed require --job; they must not rerun the entire queue")
    for option in ("max_references", "max_uploaded_genomes", "max_total_genomes", "max_upload_bytes", "max_threads", "job_timeout", "interval"):
        if getattr(args, option) < 1:
            p.error(option.replace("_", "-") + " must be positive")
    if args.max_references > 100:
        p.error("max-references cannot exceed the 100-reference selection contract")
    if args.job:
        args.once = True
    if args.job and any(not SLUG_RE.fullmatch(code) for code in args.job):
        p.error("Invalid --job result code")
    return args


def main() -> int:
    args = parse_args()
    configure_logging(Path(args.log_file), verbose=args.verbose)
    lock_path = Path(args.lock_file)

    try:
        lock_handle = acquire_lock(lock_path)
    except FileExistsError:
        logging.warning("Lock already exists: %s (another worker running?)", lock_path)
        return 0

    stop = threading.Event()
    def request_stop(_signum, _frame):
        logging.info("Shutdown requested; finishing the current job before exiting.")
        stop.set()
    signal.signal(signal.SIGTERM, request_stop)
    signal.signal(signal.SIGINT, request_stop)
    try:
        s3 = make_s3(args.region)
        while True:
            slugs = args.job or list_job_slugs(s3, args.input_bucket)
            logging.info("Found %d GToTree job(s).", len(slugs))
            for slug in slugs:
                if stop.is_set():
                    break
                try:
                    process_job(args, s3, slug)
                except Exception:
                    if not args.continue_on_error:
                        raise
                    logging.exception("Continuing after failure on %s", slug)
            if args.once or stop.wait(args.interval):
                break
    finally:
        release_lock(lock_handle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
