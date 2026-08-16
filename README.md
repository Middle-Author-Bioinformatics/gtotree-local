# GToTree backend (worker + report generator)

Two pieces that run on your GToTree host (the "tower"), modeled on the FeGenie
worker pattern:

## `gtotree_worker.py`
Polls the input S3 bucket for jobs (`<slug>/form-data.txt`), runs GToTree, and
publishes frontend-ready results to the results bucket. Same design as the
FeGenie worker: `--once` cron mode, lock file, `status.json` state machine
(running/complete/failed), structured logging, and a `raw-results.tar.gz`.

It also parses the web app's `form-data.txt`, re-validates the advanced
`GtotreeFlags` against an allow-list (defense against command injection), builds
the `-a/-A/-g/-f` inputs, runs GToTree, normalizes outputs to stable names, and
invokes `gtotree_report.py` to bake the interactive dashboard.

```bash
python3 gtotree_worker.py --once \
  --input-bucket midauthorbio-gtotree-input \
  --results-bucket midauthorbio-gtotree-results \
  --work-root /data/gtotree-worker \
  --command-prefix "conda run -n gtotree"
```

Cron (every 5 min):
```
*/5 * * * * /usr/bin/python3 /opt/gtotree/gtotree_worker.py --once --clean --continue-on-error >> /var/log/gtotree-worker.log 2>&1
```

Published to `s3://<results-bucket>/<slug>/`:
```
output/GToTree_output.tre          final Newick tree
output/Aligned_SCGs.faa            concatenated SCG alignment
output/Genomes_summary_info.tsv    per-genome QC + taxonomy
output/SCG_hit_counts.tsv          per-genome SCG hit matrix
output/gtotree-runlog.txt          run log  (+ Partitions.txt, citations.txt, iToL-colors.txt if present)
gtotree-report.html                self-contained interactive dashboard
result.json                        viewer manifest
status.json                        job state
raw-results.tar.gz                 full GToTree output dir
run.log
```

## `gtotree_report.py`
Bakes a **self-contained, interactive HTML dashboard** from the GToTree outputs,
styled like the MAB DESeq report. Panels:

- **Tree** — phylogram / cladogram / radial; click a tip to re-root, click a node
  to collapse/expand, tip highlight/search, support-value toggle, export SVG/Newick.
- **Alignment** — residue-colored concatenated SCG alignment, position ruler,
  windowed scrolling, order-by-tree.
- **SCG heatmap** — per-genome single-copy-gene hit matrix (Plotly).
- **Genomes** — sortable/filterable per-genome QC table.
- **Taxonomy** — sunburst of NCBI/GTDB lineages (present when run with `-t`/`-D`).

The tree, alignment, table, and sunburst are dependency-free JS (work offline).
Plotly is inlined when the `plotly` package is installed; otherwise a CDN tag is
used (only the heatmap tab then needs internet).

```bash
python3 gtotree_report.py -o report.html --title my_tree --gene-set Bacteria \
  --tree   GToTree_output/GToTree_output.tre \
  --alignment GToTree_output/Aligned_SCGs.faa \
  --summary   GToTree_output/Genomes_summary_info.tsv \
  --hits      GToTree_output/SCG_hit_counts.tsv
```
Keep `report_assets/` (dashboard.css + dashboard.js) beside `gtotree_report.py`.

## Install
```bash
pip install -r requirements.txt   # boto3 required; plotly optional
```

## The React app shows the same dashboard two ways
- **In-app**: the results page has a **Dashboard** tab that embeds
  `gtotree-report.html` (plus Tree / Run log / Files tabs).
- **Offline**: download `gtotree-report.html` from the Files tab and open it in
  any browser.
