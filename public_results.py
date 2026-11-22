"""Public scientific output policy. Server diagnostics stay in the work directory."""
from __future__ import annotations

import json
import re
from pathlib import Path

PUBLIC_NAMES = frozenset({
    "GToTree_output.tre", "Aligned_SCGs.faa", "Genomes_summary_info.tsv",
    "SCG_hit_counts.tsv", "selected-genomes.json", "Partitions.txt",
    "citations.txt", "iToL-colors.txt",
})
PUBLIC_REPORTS = frozenset({"gtotree-report.html", "report.html"})
PUBLIC_ERROR = "Analysis failed. Contact MAB with the result code."
_POSIX = re.compile(r'''(?<![\w:/<])/(?!/)[^\s/:"'<>;(),\[\]{}\\]+(?:/[^\s/:"'<>;(),\[\]{}\\]+)*/?''')
_WINDOWS = re.compile(r'''\b[A-Za-z]:[\\/][^\s:"'<>;(),\[\]{}]+''')
_FILE_URI = re.compile(r'''file://[^\s"'<>;(),\[\]{}]+''', re.I)


def public_text(text: str) -> str:
    """Replace absolute POSIX/Windows/file-URI paths with their leaf names.

    Relative biological identifiers, ratios and HTTP(S) citation URLs survive.
    Apply to data, never to executable JavaScript/HTML/CSS source.
    """
    def leaf(match):
        return match.group(0).replace("\\", "/").rstrip("/").rsplit("/", 1)[-1] or "[path]"
    return _POSIX.sub(leaf, _WINDOWS.sub(leaf, _FILE_URI.sub(leaf, text)))


def public_data(value):
    if isinstance(value, str):
        return public_text(value)
    if isinstance(value, list):
        return [public_data(v) for v in value]
    if isinstance(value, dict):
        return {public_text(k): public_data(v) for k, v in value.items()}
    return value


def is_public_file(relative: str) -> bool:
    if relative in PUBLIC_REPORTS:
        return True
    parts = relative.split("/")
    return (len(parts) == 1 and parts[0] in PUBLIC_NAMES) or (
        len(parts) == 2 and parts[0] == "output" and parts[1] in PUBLIC_NAMES)


def sanitize_file(path: Path) -> None:
    text = path.read_text(encoding="utf-8", errors="replace")
    if path.suffix == ".json":
        text = json.dumps(public_data(json.loads(text)), ensure_ascii=False, indent=2) + "\n"
    else:
        text = public_text(text)
    path.write_text(text, encoding="utf-8")


def public_status(data: dict, slug: str) -> dict:
    if data.get("state") not in {"running", "complete", "failed"}:
        raise ValueError("Invalid job status")
    result = {"slug": slug, "state": data["state"]}
    for key in ("started_at", "completed_at", "failed_at"):
        if isinstance(data.get(key), str):
            result[key] = public_text(data[key])[:64]
    if isinstance(data.get("has_report"), bool):
        result["has_report"] = data["has_report"]
    if isinstance(data.get("files"), list):
        result["files"] = [f for f in data["files"] if isinstance(f, str) and is_public_file(f)]
    if data["state"] == "failed":
        result["error"] = PUBLIC_ERROR
    return result


def public_manifest(data: dict, slug: str) -> dict:
    result = {"slug": slug}
    for key in ("run_name", "gene_set", "has_tree", "has_alignment", "has_report", "input_summary"):
        if key in data:
            result[key] = public_data(data[key])
    result["files"] = [f for f in data.get("files", []) if isinstance(f, str) and is_public_file(f)]
    return result
