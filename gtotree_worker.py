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


def download_prefix(s3, bucket: str, prefix: str, dest: Path) -> list[Path]:
    dest.mkdir(parents=True, exist_ok=True)
    downloaded: list[Path] = []
    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith("/"):
                continue
            rel = Path(key).relative_to(prefix)
            local = dest / rel
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
        cmd.extend(args.command_prefix.split())
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

    cmd += m.extra_flags
    cmd += ["-o", str(out_dir), "-F"]
    return cmd


def run_command(cmd: list[str], log_handle) -> None:
    logging.info("Running command: %s", " ".join(cmd))
    proc = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"GToTree failed with exit code {proc.returncode}")


# ---------------------------------------------------------------------------
# output collection
# ---------------------------------------------------------------------------
def find_first(out_dir: Path, patterns: Iterable[str]) -> Path | None:
    """Return the first file whose name matches any regex (case-insensitive)."""
    compiled = [re.compile(p, re.I) for p in patterns]
    for path in sorted(out_dir.rglob("*")):
        if path.is_file() and any(c.search(path.name) for c in compiled):
            return path
    return None


def collect_outputs(out_dir: Path) -> dict[str, Path | None]:
    """Locate GToTree's primary outputs regardless of the -o dir name."""
    return {
        # <out>.tre  or  <out>_aligned_SCGs_mod_names.tre
        "tree": find_first(out_dir, [r"\.tre$", r"\.tree$", r"\.nwk$", r"\.newick$"]),
        # Aligned_SCGs.faa or Aligned_SCGs_mod_names.faa
        "alignment": find_first(out_dir, [r"^Aligned_SCGs.*\.faa$", r"aligned.*\.faa$"]),
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
            cmd.extend(args.command_prefix.split())
        cmd += ["python", str(script), "-o", str(report_dest), "--title", manifest.job_name]
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
            run_command(cmd, log_handle)
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

    if s3_key_exists(s3, args.results_bucket, status_key) and not args.force:
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

    manifest = JobManifest(slug=slug)
    try:
        download_prefix(s3, args.input_bucket, input_prefix, job_dir)
        manifest = parse_manifest(job_dir / "form-data.txt", fallback_slug=slug)
        if not manifest.gene_set:
            raise RuntimeError("No 'SCG set' in form-data.txt (required for -H).")
        if not (manifest.accessions or manifest.proteins or manifest.genbanks or manifest.genomes):
            raise RuntimeError("No genome inputs (accessions / protein / GenBank / genome).")

        final_dir.mkdir(parents=True, exist_ok=True)
        with log_path.open("w", encoding="utf-8") as log_handle:
            log_handle.write(f"GToTree worker started {utc_now()}\nSlug: {slug}\n")
            cmd = build_gtotree_command(args, manifest, job_dir, out_dir)
            run_command(cmd, log_handle)

            outputs = collect_outputs(out_dir)
            if not (outputs["tree"] or outputs["alignment"]):
                raise RuntimeError("GToTree produced neither a tree nor an alignment.")

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

        # tarball of the full raw output dir
        tar_path = job_root / "raw-results.tar.gz"
        make_tarball(out_dir, tar_path)

        # ---- upload everything ----
        content_types = {
            ".tre": "text/plain", ".faa": "text/plain", ".tsv": "text/tab-separated-values",
            ".txt": "text/plain", ".html": "text/html",
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
            "has_report": report_ok, "files": uploaded_files,
        })
        # completion status
        put_json(s3, args.results_bucket, status_key, {
            "slug": slug, "state": "complete", "completed_at": utc_now(),
            "result_prefix": result_prefix, "has_report": report_ok, "files": uploaded_files,
        })
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
        if args.clean and job_root.exists():
            shutil.rmtree(job_root)


# ---------------------------------------------------------------------------
# lock + main loop
# ---------------------------------------------------------------------------
def acquire_lock(lock_path: Path):
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, str(os.getpid()).encode())
    os.close(fd)


def release_lock(lock_path: Path):
    try:
        lock_path.unlink()
    except FileNotFoundError:
        pass


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
    p.add_argument("--continue-on-error", action="store_true", help="Continue after a job failure.")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--lock-file", default=os.getenv("GTOTREE_LOCK_FILE", "/tmp/gtotree-worker.lock"))
    p.add_argument("--log-file", default=os.getenv("GTOTREE_WORKER_LOG", "/tmp/gtotree-worker.log"))
    return p.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(Path(args.log_file), verbose=args.verbose)
    lock_path = Path(args.lock_file)

    try:
        acquire_lock(lock_path)
    except FileExistsError:
        logging.warning("Lock already exists: %s (another worker running?)", lock_path)
        return 0

    s3 = make_s3(args.region)
    try:
        while True:
            slugs = list_job_slugs(s3, args.input_bucket)
            logging.info("Found %d GToTree job(s).", len(slugs))
            for slug in slugs:
                try:
                    process_job(args, s3, slug)
                except Exception:
                    if not args.continue_on_error:
                        raise
                    logging.exception("Continuing after failure on %s", slug)
            if args.once:
                break
            time.sleep(args.interval)
    finally:
        release_lock(lock_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
