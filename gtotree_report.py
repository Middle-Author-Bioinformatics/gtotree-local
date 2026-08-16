#!/usr/bin/env python3
"""
gtotree_report.py

Generate a self-contained, interactive HTML dashboard from GToTree outputs.
Styled after the Middle Author Bioinformatics DESeq report.

Panels
------
  1. Tree        - interactive phylogram/cladogram/radial viewer: pan/zoom,
                   reroot (click a tip), collapse/expand clades (click a node),
                   tip highlight/search, support-value toggle, export SVG/Newick.
  2. Alignment   - concatenated SCG alignment viewer: residue coloring,
                   position ruler, windowed scrolling, order-by-tree.
  3. SCG heatmap - per-genome single-copy-gene hit matrix (Plotly heatmap).
  4. Genomes     - sortable/filterable per-genome QC summary table.
  5. Taxonomy    - sunburst of NCBI/GTDB lineages (when -t/-D was used).

The output HTML is fully self-contained: CSS + dashboard JS are inlined, and
Plotly is bundled inline (falls back to a CDN <script> if plotly isn't
importable at generation time). Opens in any modern browser offline.

Usage
-----
  python gtotree_report.py -o report.html \
      --title "my_genomes_tree" \
      --tree GToTree_output/GToTree_output.tre \
      --alignment GToTree_output/Aligned_SCGs.faa \
      --summary GToTree_output/Genomes_summary_info.tsv \
      --hits GToTree_output/SCG_hit_counts.tsv

All inputs are optional; panels with no data are hidden automatically.

Dependencies: none required (pure stdlib). If `plotly` is installed the Plotly
JS bundle is inlined for a truly offline file; otherwise a CDN tag is used.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ASSET_DIR = Path(__file__).resolve().parent / "report_assets"


# ---------------------------------------------------------------------------
# parsers
# ---------------------------------------------------------------------------
def read_text(path: Path | None) -> str:
    if path and path.exists():
        return path.read_text(encoding="utf-8", errors="replace")
    return ""


def parse_fasta(path: Path | None):
    """Return {'ids': [...], 'seqs': [...]} for an aligned FASTA, or None."""
    if not path or not path.exists():
        return None
    ids, seqs, cur = [], [], []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur)); cur = []
                ids.append(line[1:].split()[0])
            elif line:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    if not ids:
        return None
    return {"ids": ids, "seqs": seqs}


def sniff_delim(sample: str) -> str:
    return "\t" if sample.count("\t") >= sample.count(",") else ","


def parse_table(path: Path | None):
    """Return {'columns': [...], 'rows': [[...]]} for a TSV/CSV, or None."""
    if not path or not path.exists():
        return None
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        return None
    delim = sniff_delim(text.splitlines()[0])
    reader = csv.reader(text.splitlines(), delimiter=delim)
    rows = [r for r in reader if r]
    if not rows:
        return None
    return {"columns": rows[0], "rows": rows[1:]}


def parse_hit_counts(path: Path | None):
    """SCG_hit_counts.tsv: first column = genome, remaining = per-gene counts.

    Returns {'genomes': [...], 'genes': [...], 'matrix': [[int]]} or None.
    """
    tbl = parse_table(path)
    if not tbl or len(tbl["columns"]) < 2:
        return None
    genes = tbl["columns"][1:]
    genomes, matrix = [], []
    for r in tbl["rows"]:
        if not r:
            continue
        genomes.append(r[0])
        vals = []
        for v in r[1:len(genes) + 1]:
            try:
                vals.append(int(float(v)))
            except (ValueError, TypeError):
                vals.append(0)
        # pad short rows
        while len(vals) < len(genes):
            vals.append(0)
        matrix.append(vals)
    if not genomes:
        return None
    return {"genomes": genomes, "genes": genes, "matrix": matrix}


# ---------------------------------------------------------------------------
# plotly bundle
# ---------------------------------------------------------------------------
def plotly_script_tag() -> str:
    try:
        import plotly  # noqa
        from plotly.offline import get_plotlyjs
        return "<script>" + get_plotlyjs() + "</script>"
    except Exception:
        # Fallback: CDN (report still fully works online; only the SCG heatmap
        # tab needs Plotly — tree/alignment/table/sunburst are dependency-free).
        return '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>'


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------
HTML_SHELL = """<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GToTree report &mdash; {title}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Space+Grotesk:wght@500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<style>
{css}
</style>
</head>
<body>
<header class="topbar">
  <div class="container topbar-inner">
    <div class="brand">
      <div class="brand-logo">MAB</div>
      <div>
        <div class="brand-title">GToTree Web Studio</div>
        <div class="brand-sub">Phylogenomics report</div>
      </div>
    </div>
    <nav class="tabs">
      <a class="nav-link active" data-tab="tree">Tree</a>
      <a class="nav-link" data-tab="alignment">Alignment</a>
      <a class="nav-link" data-tab="scg">SCG heatmap</a>
      <a class="nav-link" data-tab="genomes">Genomes</a>
      <a class="nav-link" data-tab="taxonomy">Taxonomy</a>
    </nav>
    <button id="themeBtn" class="icon-btn" title="Toggle theme">&#9790;</button>
  </div>
</header>

<div class="container">
  <section class="hero">
    <span class="pill">&#9679; GToTree &middot; phylogenomic tree</span>
    <h1>{title}</h1>
    <p class="lead">Interactive phylogenomic report. Explore and manipulate the tree, inspect the concatenated single-copy-gene alignment, review per-genome QC, and browse taxonomy &mdash; all in one offline-capable file.</p>
    <div class="hero-meta">
      <div><div class="m-k">Pipeline</div><div class="m-v">GToTree</div></div>
      <div><div class="m-k">SCG set</div><div class="m-v">{gene_set}</div></div>
      <div><div class="m-k">Generated</div><div class="m-v">{generated}</div></div>
    </div>
  </section>

  <div class="kpi-grid">
    <div class="kpi"><div class="k-label">&#127793; Tips in tree</div><div class="k-val" id="kpiTips">&mdash;</div><div class="k-sub">genomes placed</div></div>
    <div class="kpi"><div class="k-label">&#9906; Internal nodes</div><div class="k-val" id="kpiInternal">&mdash;</div><div class="k-sub">splits inferred</div></div>
    <div class="kpi"><div class="k-label">&#129516; Alignment</div><div class="k-val" id="kpiAln">&mdash;</div><div class="k-sub">concatenated SCGs</div></div>
    <div class="kpi"><div class="k-label">&#9881; SCG set</div><div class="k-val" id="kpiGeneSet" style="font-size:1.1rem">&mdash;</div><div class="k-sub">-H target</div></div>
  </div>

  <!-- TREE -->
  <div class="tabpane active" id="tab-tree">
    <div class="card card-pad">
      <div class="section-head"><h2>Interactive tree</h2><span class="meta">GToTree_output.tre</span></div>
      <p class="section-desc">Click a tip to re-root there. Click an internal node to collapse/expand its clade. Search highlights matching tips.</p>
      <div class="toolbar">
        <div class="btn-group">
          <button class="btn active" id="btnPhylo">Phylogram</button>
          <button class="btn" id="btnClado">Cladogram</button>
          <button class="btn" id="btnRadial">Radial</button>
        </div>
        <button class="btn" id="btnSupport">Support values</button>
        <button class="btn" id="btnExpand">Expand all</button>
        <button class="btn" id="btnResetRoot">Reset root</button>
        <div class="btn-group"><button class="btn" id="btnZoomOut">&minus;</button><button class="btn" id="btnZoomIn">+</button></div>
        <div class="spacer"></div>
        <div class="field">&#128269; <input id="treeSearch" placeholder="highlight tip..."></div>
        <button class="btn" id="btnNewick">&#8681; Newick</button>
        <button class="btn" id="btnTreeSVG">&#8681; SVG</button>
      </div>
      <div class="canvas" id="treeCanvas"></div>
    </div>
  </div>

  <!-- ALIGNMENT -->
  <div class="tabpane" id="tab-alignment">
    <div class="card card-pad">
      <div class="section-head"><h2>Alignment</h2><span class="meta" id="alnMeta"></span></div>
      <p class="section-desc">Concatenated single-copy-gene alignment. Residues are colored by biochemical group; rows follow the current tree order.</p>
      <div class="toolbar">
        <button class="btn active" id="alnColor">Color residues</button>
        <button class="btn active" id="alnSort">Order by tree</button>
        <div class="spacer"></div>
        <button class="btn" id="alnPrev">&#8592; prev</button>
        <input id="alnSlider" type="range" min="0" max="0" value="0" style="width:220px">
        <button class="btn" id="alnNext">next &#8594;</button>
      </div>
      <div class="canvas" id="alnCanvas"></div>
    </div>
  </div>

  <!-- SCG HEATMAP -->
  <div class="tabpane" id="tab-scg">
    <div class="card card-pad">
      <div class="section-head"><h2>SCG presence heatmap</h2><span class="meta">SCG_hit_counts.tsv</span></div>
      <p class="section-desc">Number of hits to each target single-copy gene per genome. Sparse rows can indicate incomplete genomes; multiple hits can indicate redundancy.</p>
      <div id="scgCanvas"></div>
    </div>
  </div>

  <!-- GENOMES -->
  <div class="tabpane" id="tab-genomes">
    <div class="card card-pad">
      <div class="section-head"><h2>Genome QC summary</h2><span class="meta" id="genomesMeta"></span></div>
      <p class="section-desc">Per-genome summary including estimated completeness/redundancy, hit counts, and labels. Click a header to sort.</p>
      <div class="toolbar"><div class="field">&#128269; <input id="tblFilter" placeholder="filter genomes..."></div></div>
      <div id="genomesCanvas"></div>
    </div>
  </div>

  <!-- TAXONOMY -->
  <div class="tabpane" id="tab-taxonomy">
    <div class="card card-pad">
      <div class="section-head"><h2>Taxonomy sunburst</h2><span class="meta">from lineage columns</span></div>
      <p class="section-desc">Hierarchical composition of the placed genomes across taxonomic ranks. Hover a wedge for its rank and genome count.</p>
      <div id="taxCanvas"></div>
    </div>
  </div>

  <footer>
    Generated by <strong>Middle Author Bioinformatics</strong> &middot; GToTree Web Studio. Please cite GToTree: Lee, M. D. (2019), <em>Bioinformatics</em> 35(20):4162&ndash;4164. This file runs entirely offline.
  </footer>
</div>

<script>window.GT_DATA = {data_json};</script>
{plotly}
<script>
{js}
</script>
</body>
</html>
"""


def build_html(title, gene_set, tree_nwk, alignment, summary, hits) -> str:
    css = read_text(ASSET_DIR / "dashboard.css")
    js = read_text(ASSET_DIR / "dashboard.js")
    data = {
        "title": title,
        "gene_set": gene_set or "",
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "newick": tree_nwk or "",
        "alignment": alignment,
        "summary": summary,
        "hits": hits,
    }
    data_json = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    # guard against </script> inside data breaking the inline script
    data_json = data_json.replace("</", "<\\/")
    return HTML_SHELL.format(
        title=html.escape(title),
        gene_set=html.escape(gene_set or "\u2014"),
        generated=data["generated"],
        css=css,
        js=js,
        plotly=plotly_script_tag(),
        data_json=data_json,
    )


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Build a self-contained interactive GToTree HTML dashboard.")
    p.add_argument("-o", "--output", type=Path, required=True, help="Output HTML path.")
    p.add_argument("--title", default="GToTree run", help="Report/run title.")
    p.add_argument("--gene-set", default="", help="SCG set (-H target) for display.")
    p.add_argument("--tree", type=Path, help="Newick tree file (.tre).")
    p.add_argument("--alignment", type=Path, help="Aligned SCGs FASTA (.faa).")
    p.add_argument("--summary", type=Path, help="Genomes_summary_info.tsv.")
    p.add_argument("--hits", type=Path, help="SCG_hit_counts.tsv.")
    args = p.parse_args(argv)

    tree_nwk = read_text(args.tree).strip() if args.tree else ""
    alignment = parse_fasta(args.alignment)
    summary = parse_table(args.summary)
    hits = parse_hit_counts(args.hits)

    gene_set = args.gene_set
    # try to recover the SCG set from the summary if not passed
    if not gene_set and summary:
        pass  # summary doesn't carry the SCG set; leave as-is

    if not (tree_nwk or alignment or summary or hits):
        print("Error: no usable GToTree outputs provided.", file=sys.stderr)
        return 2

    html_out = build_html(args.title, gene_set, tree_nwk, alignment, summary, hits)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html_out, encoding="utf-8")
    n_tips = alignment["ids"] if alignment else []
    print(f"Wrote {args.output}  (tree={'yes' if tree_nwk else 'no'}, "
          f"alignment={len(n_tips)} seqs, "
          f"summary={'yes' if summary else 'no'}, hits={'yes' if hits else 'no'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
