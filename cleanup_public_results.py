#!/usr/bin/env python3
"""Audit public results; --apply backs up originals locally before sanitizing.

The apply command holds the cron's outer flock; stop the persistent worker first.
No analysis is rerun. Only recognized GToTree result-code folders are processed.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path

import boto3
from gtotree_report import build_html
from public_results import is_public_file, public_data, public_manifest, public_status, public_text

CODE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
REPORT_DATA = re.compile(r"<script>\s*window\.GT_DATA\s*=\s*(.*?);\s*</script>", re.S)


def plan_job(s3, bucket, slug, objects, reports_only=False):
    """Prepare all safe replacements before changing any remote object."""
    def read(name):
        return s3.get_object(Bucket=bucket, Key=slug + "/" + name)["Body"].read()
    updates, removals = {}, []
    available = sorted(name for name in objects if is_public_file(name))
    for name in objects:
        if reports_only and not (is_public_file(name) and name.endswith(".html")):
            continue
        if name in ("status.json", "result.json"):
            data = json.loads(read(name))
            if name == "status.json":
                safe = public_status(data, slug)
            else:
                safe = public_manifest(data, slug)
            safe["files"] = available
            updates[name] = (json.dumps(safe, indent=2) + "\n").encode()
        elif is_public_file(name):
            raw = read(name).decode("utf-8", errors="replace")
            if name.endswith(".html"):
                match = REPORT_DATA.search(raw)
                if not match:
                    raise ValueError("Report data could not be read; this result was left unchanged")
                data = json.loads(match.group(1))
                safe = build_html(data.get("title", slug), data.get("gene_set", ""),
                                  data.get("newick", ""), data.get("alignment"), data.get("summary"),
                                  data.get("hits"), data.get("reference_metadata"), slug)
            elif name.endswith(".json"):
                safe = json.dumps(public_data(json.loads(raw)), ensure_ascii=False, indent=2) + "\n"
            else:
                safe = public_text(raw)
            updates[name] = safe.encode("utf-8")
        else:
            removals.append(name)
    return updates, removals


def apply_job(s3, bucket, slug, objects, updates, removals, backup_root):
    # Back up every affected object before the first remote write/delete.
    backup = backup_root / slug
    backup.mkdir(parents=True, exist_ok=False, mode=0o700)
    os.chmod(backup, 0o700)
    manifest = {}
    for name in sorted(set(updates) | set(removals)):
        key = slug + "/" + name
        filename = hashlib.sha256(key.encode()).hexdigest()
        target = backup / filename
        s3.download_file(bucket, key, str(target))
        os.chmod(target, 0o600)
        if "Size" in objects[name] and target.stat().st_size != objects[name]["Size"]:
            raise RuntimeError("Backup size mismatch; no remote objects were changed")
        with target.open("rb") as handle:
            os.fsync(handle.fileno())
        manifest[key] = filename
    with (backup / "index.json").open("x", encoding="utf-8") as handle:
        os.chmod(backup / "index.json", 0o600)
        json.dump(manifest, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    # Scientific data and safe metadata first; final status last.
    def publish(name):
        suffix = Path(name).suffix
        content_type = {".html": "text/html", ".json": "application/json",
                        ".tsv": "text/tab-separated-values"}.get(suffix, "text/plain")
        s3.put_object(Bucket=bucket, Key=slug + "/" + name, Body=updates[name],
                      ContentType=content_type, CacheControl="no-cache")
    for name in sorted(updates):
        if name != "status.json":
            publish(name)
    for name in removals:
        s3.delete_object(Bucket=bucket, Key=slug + "/" + name)
    if "status.json" in updates:
        publish("status.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bucket", default="midauthorbio-gtotree-results")
    parser.add_argument("--region", default="us-east-2")
    parser.add_argument("--job", action="append", help="Only this result code; repeat for several.")
    parser.add_argument("--apply", action="store_true", help="Back up and apply; default is a read-only audit.")
    parser.add_argument("--reports-only", action="store_true", help="Rebuild existing HTML reports without changing any other objects or rerunning analysis.")
    parser.add_argument("--lock-file", type=Path, default=Path("/tmp/gtotree-scraper.lock"))
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    parser.add_argument("--backup-dir", type=Path, default=Path.home() / ".mab-private" / "gtotree-results" / stamp)
    args = parser.parse_args()
    if args.job and any(not CODE.fullmatch(code) for code in args.job):
        parser.error("Invalid result code")
    os.umask(0o077)
    lock = None
    if args.apply:
        lock = args.lock_file.open("a+")
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            lock.close()
            parser.error("The worker is running. Stop it gracefully and let its current job finish before applying cleanup.")
    s3 = boto3.client("s3", region_name=args.region)
    jobs = {}
    for page in s3.get_paginator("list_objects_v2").paginate(Bucket=args.bucket):
        for obj in page.get("Contents", []):
            parts = obj["Key"].split("/", 1)
            if len(parts) == 2 and CODE.fullmatch(parts[0]) and (not args.job or parts[0] in args.job):
                jobs.setdefault(parts[0], {})[parts[1]] = obj
    failures = 0
    for slug, objects in sorted(jobs.items()):
        if not any(n in objects for n in ("status.json", "result.json", "gtotree-report.html", "report.html")):
            continue
        try:
            plan = plan_job(s3, args.bucket, slug, objects, reports_only=args.reports_only)
            updates, removals = plan
            action = f"refresh {len(updates)} HTML report(s); preserve all other files" if args.reports_only else f"sanitize {len(updates)} scientific/status files; remove {len(removals)} diagnostic/other files"
            print(f"{slug}: {action}", flush=True)
            if args.apply:
                apply_job(s3, args.bucket, slug, objects, updates, removals, args.backup_dir)
        except Exception as exc:
            failures += 1
            print(f"ERROR {slug}: {exc}", flush=True)
    if args.apply:
        print(f"Private originals: {args.backup_dir}")
    else:
        print("Read-only audit. Add --apply to back up and update the public results.")
    if lock:
        lock.close()
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
