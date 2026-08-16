import boto3
import os
import subprocess
import json
import re
import shlex
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

# GToTree "tower": polls the INPUT bucket for new submissions, runs GToTree per
# submission, and writes results to the RESULTS bucket under the same flat slug.
#
# There is NO emailing — runs are addressed purely by their result code (slug),
# and the web app's results viewer polls the results bucket.

s3_client = boto3.client("s3")

# ---- buckets (override via env) -------------------------------------------
INPUT_BUCKET = os.environ.get("GTOTREE_INPUT_BUCKET", "midauthorbio-gtotree-input")
RESULTS_BUCKET = os.environ.get("GTOTREE_RESULTS_BUCKET", "midauthorbio-gtotree-results")

# ---- local working dirs ----------------------------------------------------
LOCAL_BASE_DIR = os.environ.get("GTOTREE_DATA_DIR", "/home/ark/MAB/gtotree/data")
BASE_OUTPUT_DIR = os.environ.get("GTOTREE_RESULTS_DIR", "/home/ark/MAB/gtotree/results")
CONDA = os.environ.get("CONDA_BIN", "/home/ark/miniconda3/bin/conda")
CONDA_ENV = os.environ.get("GTOTREE_CONDA_ENV", "gtotree")

log_file_path = os.environ.get("GTOTREE_PROCESSED_LOG", "/home/ark/MAB/gtotree/processed_folders.log")
failed_log_file_path = os.environ.get("GTOTREE_FAILED_LOG", "/home/ark/MAB/gtotree/failed_folders.log")

MAX_PARALLEL_JOBS = int(os.environ.get("GTOTREE_MAX_PARALLEL_JOBS", "4"))
log_lock = threading.Lock()

# ---------------------------------------------------------------------------
# Advanced-flag allow-list. The web app's GtotreeFlags string is user-controlled
# text spliced into the command, so we re-validate every token against this set
# and drop anything unexpected. Value-taking flags are followed by one token.
# ---------------------------------------------------------------------------
ALLOWED_GTOTREE_FLAGS = {
    "-t": 0, "-D": 0, "-B": 0, "-N": 0, "-k": 0, "-X": 0, "-P": 0,
    "-L": 1, "-c": 1, "-G": 1, "-T": 1, "-n": 1, "-M": 1, "-j": 1,
}
_TREE_PROGRAMS = {"FastTree", "FastTreeMP", "IQ-TREE"}


def load_seen_folders(path):
    if os.path.exists(path):
        with open(path) as f:
            return set(line.strip() for line in f)
    return set()


def append_seen_folder(path, folder):
    with log_lock:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "a") as f:
            f.write(folder + "\n")


def _kv(line):
    if ":" not in line:
        return None, None
    k, v = line.split(":", 1)
    return k.strip(), v.strip()


def _sanitize_flags(raw):
    """Validate the web app's GtotreeFlags against the allow-list."""
    tokens = shlex.split(raw) if raw else []
    clean, i = [], 0
    while i < len(tokens):
        tok = tokens[i]
        if tok not in ALLOWED_GTOTREE_FLAGS:
            print(f"[flags] dropping unrecognized token: {tok!r}")
            i += 1
            continue
        if ALLOWED_GTOTREE_FLAGS[tok]:
            if i + 1 >= len(tokens):
                print(f"[flags] {tok} expects a value; dropping")
                break
            val = tokens[i + 1]
            if tok == "-T" and val not in _TREE_PROGRAMS:
                print(f"[flags] bad tree program {val!r}; dropping")
                i += 2
                continue
            if tok in ("-c", "-G") and not re.match(r"^\d*\.?\d+$", val):
                print(f"[flags] {tok} non-numeric {val!r}; dropping")
                i += 2
                continue
            if tok in ("-L",) and not re.match(r"^[A-Za-z,]+$", val):
                print(f"[flags] {tok} bad ranks {val!r}; dropping")
                i += 2
                continue
            if tok in ("-n", "-M", "-j") and not re.match(r"^\d+$", val):
                print(f"[flags] {tok} non-integer {val!r}; dropping")
                i += 2
                continue
            clean.extend([tok, val])
            i += 2
        else:
            clean.append(tok)
            i += 1
    return clean


def list_folders_in_bucket(bucket_name):
    """Find run folders by locating <slug>/form-data.txt at the bucket root."""
    paginator = s3_client.get_paginator("list_objects_v2")
    folders = set()
    for page in paginator.paginate(Bucket=bucket_name):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith("form-data.txt"):
                parts = key.split("/")
                if len(parts) == 2 and parts[0]:
                    folders.add(f"{parts[0]}/")
    return sorted(folders)


def download_s3_folder(bucket_name, s3_folder, local_folder):
    os.makedirs(local_folder, exist_ok=True)
    paginator = s3_client.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket_name, Prefix=s3_folder):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith("/"):
                continue
            rel = os.path.relpath(key, s3_folder)
            dest = os.path.join(local_folder, rel)
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            print(f"Downloading s3://{bucket_name}/{key} -> {dest}")
            s3_client.download_file(bucket_name, key, dest)


def extract_form_data(folder_path):
    """Parse the web app's form-data.txt.

        Name: <run name>
        SCG set: <gene set>                 -> -H target (required)
        Accessions File: <file|N/A>         -> -a
        Protein File: <file>   (repeated)   -> -A
        GenBank File: <file>   (repeated)   -> -g
        Genome File: <file>    (repeated)   -> -f
        KOFile: <file|N/A>                  -> -K
        PfamFile: <file|N/A>                -> -p
        GtotreeFlags: <extra flags>

    Returns a dict of everything needed to build the command.
    """
    form = os.path.join(folder_path, "form-data.txt")
    name = "GToTree_output"
    gene_set = None
    accessions = None
    proteins, genbanks, genomes = [], [], []
    ko_file = pfam_file = None
    raw_flags = ""

    if os.path.exists(form):
        with open(form) as f:
            for line in f:
                key, val = _kv(line.rstrip("\n"))
                if key is None:
                    continue
                if key == "Name":
                    name = val or name
                elif key == "SCG set":
                    gene_set = val
                elif key == "Accessions File":
                    accessions = None if val == "N/A" else val
                elif key == "Protein File":
                    proteins.append(val)
                elif key == "GenBank File":
                    genbanks.append(val)
                elif key == "Genome File":
                    genomes.append(val)
                elif key == "KOFile":
                    ko_file = None if val == "N/A" else val
                elif key == "PfamFile":
                    pfam_file = None if val == "N/A" else val
                elif key == "GtotreeFlags":
                    raw_flags = val

    return {
        "name": name,
        "gene_set": gene_set,
        "accessions": accessions,
        "proteins": proteins,
        "genbanks": genbanks,
        "genomes": genomes,
        "ko_file": ko_file,
        "pfam_file": pfam_file,
        "extra_flags": _sanitize_flags(raw_flags),
    }


def _write_list_file(folder_path, filenames, list_name):
    """Write a single-column file of absolute paths GToTree can consume (-A/-g/-f)."""
    path = os.path.join(folder_path, list_name)
    with open(path, "w") as f:
        for n in filenames:
            f.write(os.path.join(folder_path, n) + "\n")
    return path


def upload_file_to_s3(bucket, prefix, local_file):
    key = f"{prefix}{os.path.basename(local_file)}" if prefix.endswith("/") \
        else f"{prefix}/{os.path.basename(local_file)}"
    s3_client.upload_file(local_file, bucket, key)
    print(f"Uploaded {local_file} -> s3://{bucket}/{key}")


def upload_directory(local_dir, bucket, prefix):
    import mimetypes
    if not os.path.isdir(local_dir):
        return
    for root, _, files in os.walk(local_dir):
        for fn in files:
            fp = os.path.join(root, fn)
            rel = os.path.relpath(fp, local_dir)
            key = f"{prefix}/{rel}".replace("\\", "/")
            ctype, _ = mimetypes.guess_type(fp)
            if ctype is None:
                # Newick / alignment / log render nicely as text in the browser
                ctype = "text/plain" if re.search(r"\.(tre|tree|nwk|newick|faa|fa|fasta|aln|log|txt|tsv)$", fn, re.I) \
                    else "binary/octet-stream"
            s3_client.upload_file(fp, bucket, key,
                                  ExtraArgs={"ContentType": ctype, "ContentDisposition": "inline"})
            print(f"Uploaded {fp} -> s3://{bucket}/{key}")


def s3_key_exists(bucket, key):
    try:
        s3_client.head_object(Bucket=bucket, Key=key)
        return True
    except s3_client.exceptions.ClientError as e:
        if e.response.get("Error", {}).get("Code", "") in ("404", "NoSuchKey", "NotFound"):
            return False
        raise


def build_gtotree_command(fd, folder_path, output_dir):
    """Assemble the GToTree argument list from parsed form data."""
    cmd = ["GToTree"]

    if fd["accessions"]:
        cmd += ["-a", os.path.join(folder_path, fd["accessions"])]
    if fd["proteins"]:
        cmd += ["-A", _write_list_file(folder_path, fd["proteins"], "protein_files.txt")]
    if fd["genbanks"]:
        cmd += ["-g", _write_list_file(folder_path, fd["genbanks"], "genbank_files.txt")]
    if fd["genomes"]:
        cmd += ["-f", _write_list_file(folder_path, fd["genomes"], "fasta_files.txt")]

    # required HMM / SCG set
    cmd += ["-H", fd["gene_set"]]

    # KO / Pfam target files (their flags are NOT in extra_flags)
    if fd["ko_file"]:
        cmd += ["-K", os.path.join(folder_path, fd["ko_file"])]
    if fd["pfam_file"]:
        cmd += ["-p", os.path.join(folder_path, fd["pfam_file"])]

    # validated advanced flags
    cmd += fd["extra_flags"]

    # output dir + force overwrite (idempotent re-runs)
    cmd += ["-o", output_dir, "-F"]
    return cmd


def process_s3_folder(s3_folder, input_bucket, results_bucket):
    """Run one submission end-to-end and push results to the results bucket."""
    try:
        print(f"Processing {s3_folder}")
        local_folder = os.path.join(LOCAL_BASE_DIR, s3_folder)
        download_s3_folder(input_bucket, s3_folder, local_folder)

        fd = extract_form_data(local_folder)
        if not fd["gene_set"]:
            print(f"[FAILED] {s3_folder}: no SCG set in form-data.txt")
            append_seen_folder(failed_log_file_path, s3_folder)
            return False
        if not (fd["accessions"] or fd["proteins"] or fd["genbanks"] or fd["genomes"]):
            print(f"[FAILED] {s3_folder}: no genome inputs")
            append_seen_folder(failed_log_file_path, s3_folder)
            return False

        output_dir = os.path.join(BASE_OUTPUT_DIR, s3_folder.rstrip("/"), "output")
        # GToTree wants its -o dir to not exist (or use -F to overwrite)
        os.makedirs(os.path.dirname(output_dir), exist_ok=True)

        cmd = build_gtotree_command(fd, local_folder, output_dir)
        full = [CONDA, "run", "-n", CONDA_ENV, *cmd]
        print("[gtotree] Running:\n " + " ".join(cmd))

        result = subprocess.run(full, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[gtotree] ERROR in {s3_folder}:\n{result.stderr[-4000:]}")
            append_seen_folder(failed_log_file_path, s3_folder)
            return False

        # Verify a tree (or, in -N mode, an alignment) was produced.
        produced_tree = any(re.search(r"\.(tre|tree|nwk|newick)$", f, re.I)
                            for f in os.listdir(output_dir)) if os.path.isdir(output_dir) else False
        produced_aln = any(re.search(r"aligned.*\.(faa|fa|fasta)$", f, re.I)
                           for f in os.listdir(output_dir)) if os.path.isdir(output_dir) else False
        if not (produced_tree or produced_aln):
            print(f"[gtotree] FAILED — no tree/alignment in {output_dir}")
            append_seen_folder(failed_log_file_path, s3_folder)
            return False

        # ---- push results to the RESULTS bucket under <slug>/output/... ----
        upload_directory(output_dir, results_bucket, f"{s3_folder}output")

        # a small manifest for the viewer / future use
        manifest = {"slug": s3_folder.rstrip("/"), "run_name": fd["name"],
                    "gene_set": fd["gene_set"], "has_tree": produced_tree}
        man_path = os.path.join(BASE_OUTPUT_DIR, s3_folder.rstrip("/"), "result.json")
        with open(man_path, "w") as mf:
            json.dump(manifest, mf, indent=2)
        s3_client.upload_file(man_path, results_bucket, f"{s3_folder}result.json",
                              ExtraArgs={"ContentType": "application/json", "ContentDisposition": "inline"})

        append_seen_folder(log_file_path, s3_folder)
        print(f"Completed {s3_folder}")
        return True

    except Exception as exc:
        print(f"[FAILED] {s3_folder}: {exc}")
        append_seen_folder(failed_log_file_path, s3_folder)
        return False


if __name__ == "__main__":
    print(f"[CONFIG] input bucket   = {INPUT_BUCKET}")
    print(f"[CONFIG] results bucket = {RESULTS_BUCKET}")

    folders = list_folders_in_bucket(INPUT_BUCKET)
    seen = load_seen_folders(log_file_path)
    failed = load_seen_folders(failed_log_file_path)

    new_folders = []
    for f in folders:
        if f in seen:
            continue
        if f in failed:
            print(f"[SKIP] previously failed: {f}")
            continue
        # done once the manifest exists in the results bucket
        if s3_key_exists(RESULTS_BUCKET, f"{f}result.json"):
            print(f"[SKIP] existing results for {f}")
            append_seen_folder(log_file_path, f)
            continue
        new_folders.append(f)

    print(f"Folders found: {folders}")
    if not new_folders:
        print("[SCAN] no new folders")
    else:
        print(f"[PARALLEL] processing {len(new_folders)} folder(s), up to {MAX_PARALLEL_JOBS} at once")
        done = failed_ct = 0
        with ThreadPoolExecutor(max_workers=MAX_PARALLEL_JOBS) as ex:
            futs = {ex.submit(process_s3_folder, f, INPUT_BUCKET, RESULTS_BUCKET): f
                    for f in new_folders}
            for fut in as_completed(futs):
                f = futs[fut]
                try:
                    ok = fut.result()
                except Exception as exc:
                    print(f"[FAILED] worker error {f}: {exc}")
                    append_seen_folder(failed_log_file_path, f)
                    ok = False
                done += 1 if ok else 0
                failed_ct += 0 if ok else 1
        print(f"[PARALLEL] completed {done}, failed {failed_ct}")
    print("All folders processed.")
