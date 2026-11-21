"""Validation for public GToTree submissions; no AWS dependency."""
import json
import re
from pathlib import Path

ACCESSION = re.compile(r"^GC[AF]_\d+(?:\.\d+)?$")
METADATA_KEYS = ("assembly_accession", "organism_name", "infraspecific_name", "biosample", "bioproject",
                 "assembly_level", "group", "host", "isolation_source", "geo_loc_name", "country", "dataset_id")


def input_path(root, name):
    if not name or name in (".", "..") or not re.fullmatch(r"[A-Za-z0-9._-]{1,255}", name):
        raise ValueError("Inputs must have plain filenames without directories: %r" % name)
    path = (Path(root) / name).resolve()
    if path.parent != Path(root).resolve() or not path.is_file():
        raise ValueError("Missing or invalid input file: %s" % name)
    return path


def read_accessions(path, limit=100):
    if path.stat().st_size > 1024 * 1024:
        raise ValueError("Accession list exceeds 1 MB")
    lines = re.split(r"[\s,;]+", path.read_text(encoding="utf-8-sig").strip())
    accessions = list(dict.fromkeys(line.upper() for line in lines if line))
    if not accessions or any(not ACCESSION.fullmatch(a) for a in accessions):
        raise ValueError("Use a list of valid GCA_/GCF_ assembly accessions")
    if len(accessions) > limit:
        raise ValueError("At most %s reference genomes are allowed" % limit)
    return accessions


def read_metadata(path, allowed_accessions=None, limit=100):
    if path.stat().st_size > 1024 * 1024:
        raise ValueError("Reference metadata exceeds 1 MB")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("version") != 1 or data.get("source") != "assembly-grasp":
        raise ValueError("Unsupported reference-metadata format")
    rows = data.get("genomes")
    if not isinstance(rows, list) or not 1 <= len(rows) <= limit:
        raise ValueError("Invalid number of reference-metadata rows")
    clean, seen = [], set()
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("Invalid reference-metadata row")
        record = {}
        for key in METADATA_KEYS:
            v = row.get(key, "")
            if not isinstance(v, str) or len(v) > 1000 or any(c in v for c in "\r\n\x00"):
                raise ValueError("Invalid reference metadata: %s" % key)
            if v.strip():
                record[key] = v.strip()
        accession = record.get("assembly_accession", "").upper()
        if not ACCESSION.fullmatch(accession) or accession in seen:
            raise ValueError("Metadata accessions must be valid and unique")
        if allowed_accessions is not None and accession not in allowed_accessions:
            raise ValueError("Metadata accession is absent from the submitted reference list: %s" % accession)
        record["assembly_accession"] = accession
        seen.add(accession)
        clean.append(record)
    return {"version": 1, "source": "assembly-grasp", "created_at": str(data.get("created_at", ""))[:64], "genomes": clean}


def validate_manifest(m, job_dir, args):
    references = []
    if m.accessions:
        path = input_path(job_dir, m.accessions)
        references = read_accessions(path, args.max_references)
        # Normalize only the local working copy; the S3 source stays unchanged.
        path.write_text("\n".join(references) + "\n", encoding="utf-8")
    groups = [(m.proteins, {".faa", ".fa", ".fasta", ".fna"}),
              (m.genbanks, {".gb", ".gbk", ".gbff", ".gbf"}),
              (m.genomes, {".fa", ".fasta", ".fna"})]
    uploads = [name for names, _ in groups for name in names]
    declared = uploads + [x for x in (m.accessions, m.ko_file, m.pfam_file, m.metadata_file) if x]
    reserved = {"form-data.txt", "protein_files.txt", "genbank_files.txt", "fasta_files.txt"}
    if len(declared) != len(set(declared)) or any(n in reserved for n in declared):
        raise ValueError("Declared input filenames overlap or use a reserved name")
    if len(uploads) > args.max_uploaded_genomes or len(references) + len(uploads) > args.max_total_genomes:
        raise ValueError("Submission exceeds the uploaded-genome or total-genome limit")
    for names, extensions in groups:
        for name in names:
            path = input_path(job_dir, name)
            if path.suffix.lower() not in extensions or path.stat().st_size == 0:
                raise ValueError("Invalid or empty genome input: %s" % name)
            with path.open(encoding="utf-8", errors="replace") as handle:
                first = next((line.strip() for line in handle if line.strip()), "")
            if not (first.startswith("LOCUS") if extensions == {".gb", ".gbk", ".gbff", ".gbf"} else first.startswith(">")):
                raise ValueError("Input is not the expected FASTA/GenBank format: %s" % name)
    for name in (m.ko_file, m.pfam_file):
        if name:
            input_path(job_dir, name)
    if not references and not uploads:
        raise ValueError("No genome inputs supplied")
    if not re.fullmatch(r"[A-Za-z0-9_.-]{1,100}", m.gene_set):
        raise ValueError("Invalid single-copy-gene set")
    metadata = read_metadata(input_path(job_dir, m.metadata_file), set(references), args.max_references) if m.metadata_file else None
    if metadata:
        input_path(job_dir, m.metadata_file).write_text(json.dumps(metadata, ensure_ascii=False), encoding="utf-8")
    return {"reference_genomes": len(references), "uploaded_genomes": len(uploads)}


def bound_resource_flags(flags, budget):
    """Bound both per-task threads and parallel tasks to a total thread budget."""
    values, remaining = {"-n": budget, "-M": 1, "-j": 1}, []
    i = 0
    while i < len(flags):
        if flags[i] in values:
            if i + 1 >= len(flags) or not flags[i+1].isdigit() or int(flags[i+1]) < 1:
                raise ValueError("Thread/task counts must be positive integers")
            values[flags[i]] = min(int(flags[i+1]), budget)
            i += 2
        else:
            remaining.append(flags[i]); i += 1
    values["-j"] = min(values["-j"], max(1, budget // max(values["-n"], values["-M"])))
    return remaining + [item for key in ("-n", "-M", "-j") for item in (key, str(values[key]))]
