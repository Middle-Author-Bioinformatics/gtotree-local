# Assembly Grasp → GToTree: TheBelly rollout

This is the first integration milestone: reference selection, uploaded genomes, metadata, reporting, and the lab-dashboard inquiry button. It preserves existing result-code/public-result access. Shared Cognito accounts, private report authorization, revocable shares, and administrator inquiry storage are a separate cloud milestone. Do not advertise private/account-owned reports on the basis of this change.

## Existing server setup

The observed deployment is:

- Checkout: `/home/ark/MAB/bin/gtotree-local`
- Worker interpreter: `/home/ark/miniconda3/bin/python3`
- Analysis/report environment: `/home/ark/miniconda3/bin/conda run -n gtotree`
- Cron launcher: every five minutes, protected by `/usr/bin/flock -n /tmp/gtotree-scraper.lock`
- Cron output: `/home/ark/MAB/bin/cronlogs/gtotree_cron.log`

The cron command omits `--once`, so the first successful invocation stays alive and polls every five minutes. Subsequent cron invocations exit while the outer lock is held. This remains supported; **no new Amplify application and no public web server on TheBelly are required**.

## Update after the coordinated backend PR is merged

Wait until the current worker has no active analysis/report subprocess before updating and restarting it. These commands are rollout instructions, not actions already performed on the server.

1. Review current activity:

```bash
pgrep -af 'gtotree_(processing|worker)\.py'
tail -n 60 /home/ark/MAB/bin/cronlogs/gtotree_cron.log
```

Use `pstree -ap WORKER_PID` with the **current Python worker PID** from the output. If it shows a running analysis, wait for that job to finish. Do not terminate GToTree or delete lock files to force an update.

2. Update the existing checkout and report environment:

```bash
cd /home/ark/MAB/bin/gtotree-local
git status --short
git pull --ff-only
/home/ark/miniconda3/bin/conda run --no-capture-output -n gtotree \
  python -m pip install -r /home/ark/MAB/bin/gtotree-local/requirements.txt
/home/ark/miniconda3/bin/python3 \
  /home/ark/MAB/bin/gtotree-local/gtotree_worker.py --help
```

Resolve any local modifications before pulling; do not discard them. The base interpreter already imports boto3 in the running worker. The requirements install above also puts Plotly in the environment used to generate reports, so the heatmap library can be embedded offline. Keep `job_validation.py`, `gtotree_report.py`, and the complete `report_assets/` directory together with the worker.

3. When the old worker is idle, stop **only its current Python worker PID**, then let the existing cron entry launch the updated worker on its next five-minute tick:

```bash
kill -TERM WORKER_PID
```

Replace `WORKER_PID`; do not copy a PID from an earlier conversation. The old worker predates graceful shutdown, so this first restart must happen while it is idle. The updated worker handles SIGTERM by finishing its current job and exiting before starting another. Future reloads can use that behavior, but allow time for the current job to finish.

The new internal lock is an advisory OS lock. Its file remains on disk intentionally and a stale PID inside it does not block startup. Do not delete it. The outer cron `flock` command can stay exactly as it is.

4. Verify the new worker and logs:

```bash
pgrep -af 'gtotree_worker\.py'
tail -n 60 /home/ark/MAB/bin/cronlogs/gtotree_cron.log
```

There should be one Python worker, plus its shell/flock launcher. A persistent process does not reload Python modules just because `git pull` ran; the restart is required.

## New defaults and operating behavior

| Setting | Default |
|---|---|
| Unique NCBI references | 100 |
| Uploaded genomes (one genome per file) | 100 |
| Combined genomes | 200 |
| Total S3 input bytes, including metadata | 100 MiB |
| Per-analysis thread budget | 16 |
| Analysis runtime | 21,600 seconds / 6 hours |
| Concurrent analyses per worker | 1 |
| Dashboard contact | `ark@midauthorbio.com` |

Lower the reference limit if needed (its maximum is 100). Override the other limits with the documented `--max-*`, `--job-timeout`, and `--contact-email` options. CPU-related user flags are bounded to the server budget, including the product of per-task threads and parallel tasks. GToTree and subprocesses are terminated together on timeout. A single failed job does not stop the queue by default; `--stop-on-error` is available for debugging. `--clean` removes successful scratch directories and preserves failed jobs. Otherwise both are retained, as before.

Worker checks protect analysis execution. Before a public account-based launch, enforce corresponding upload/quota policies in the cloud upload service as well; these worker checks happen after objects have reached S3.

## Targeted retry

Retrying must use the same outer lock as the cron launcher. First pause/stop the persistent worker when appropriate; a targeted command cannot run concurrently with it. The new worker finishes its current analysis before exiting on SIGTERM.

```bash
/usr/bin/flock -n /tmp/gtotree-scraper.lock \
  /home/ark/miniconda3/bin/python3 \
  /home/ark/MAB/bin/gtotree-local/gtotree_worker.py \
  --once --job RESULT_CODE --retry-failed \
  --command-prefix "/home/ark/miniconda3/bin/conda run -n gtotree"
```

Replace `RESULT_CODE`. If the persistent worker is still running, the lock prevents this command from starting. `--job` implies a single scan. `--retry-failed` accepts failed or interrupted/running states; completed jobs stay skipped. To deliberately rebuild a completed report/analysis, use `--job RESULT_CODE --force`. Untargeted `--force`/`--retry-failed` is rejected, so a maintenance command cannot accidentally rerun the entire bucket. Retrying regenerates that job's local working directory; copy any prior diagnostic files you want to retain first.

## Coordinated frontend deployment

Deploy backend first, then GToTree, then Assembly Grasp. The existing Lambda upload contract is unchanged; selected metadata is an additional uploaded file. Existing jobs without metadata still run.

Set these build-time variables on the **existing** Amplify apps and rebuild them:

| App | Variable | Value |
|---|---|---|
| Assembly Grasp | `VITE_GTOTREE_APP_URL` | Existing GToTree deployed HTTPS URL |
| GToTree | `VITE_ASSEMBLY_GRASP_ORIGINS` | Exact Assembly Grasp HTTPS origin; comma-separated if both custom and Amplify domains are used |
| GToTree | `VITE_DASHBOARD_CONTACT_EMAIL` | `ark@midauthorbio.com` or your chosen inquiry address |

An origin contains scheme + hostname + optional port, with no path (for example `https://grasp.example.org`). No AWS secret or password belongs in these values. There is no wildcard origin trust. The automatic handoff verifies origin, source window, and a random nonce; genome metadata is transferred with `postMessage`, not placed in a query string or fragment. A downloaded selection JSON can always be imported manually, including when a popup is blocked or opener isolation prevents the automatic flow.

New jobs retain the existing `<result-code>/` layout and add `output/selected-genomes.json`. The report includes a searchable reference-metadata tab. It labels the metadata as submitted references, not a claim that every selected genome passed QC or was placed in the tree. The viewer reads `status.json` so a tree uploaded before the dashboard does not prematurely mark the new job complete. Legacy jobs without status files still work.

The web inquiry form prepares an email with the result code and the user's message. It does not send automatically or create an admin record. The downloaded report also contains the inquiry email link. A future cloud endpoint can replace this transport after account and inquiry-storage configuration is selected.

## Validation

```bash
python -m unittest discover -s tests -v
```

Tests use a simulated S3 client and simulated GToTree outputs; the real report generator runs. They cover mixed public/uploaded input, metadata preservation, limits, path confinement, resource budgets, locking, failed-job retention/retry, and report-before-completion publication. No live AWS requests or real phylogenomic jobs are required by the tests.
