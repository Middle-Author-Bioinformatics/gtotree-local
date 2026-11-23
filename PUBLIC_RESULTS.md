# Scientific outputs and private diagnostics

The results bucket serves curated scientific downloads and the interactive
report. The worker does not publish run logs, raw output archives or exception
details. Its original outputs and diagnostics stay in its local work directory,
subject to the existing `--clean` option. The published copies replace absolute
POSIX, Windows and file-URI paths with filenames, including labels embedded in
the report's data. Citation URLs and biological values are preserved.

The frontend companion change removes Run log, limits Files to the approved
scientific outputs, uses generic errors and preserves the selected view across
refreshes. Updating the frontend alone does not remove old S3 objects.

## Update TheBelly and audit existing results

After merging this backend PR, pull the existing checkout:

```bash
cd /home/ark/MAB/bin/gtotree-local
git pull --ff-only
/home/ark/miniconda3/bin/conda run --no-capture-output -n gtotree python cleanup_public_results.py
```

The last command is read-only. It reports the number of scientific/status files
to sanitize and diagnostic/other files to remove for each result code. Use
`--job wNMFZMt46n` to target just the example report, or repeat `--job` for several.
The default covers recognized result folders in `midauthorbio-gtotree-results`.

## Apply once, with private backups

Stop the persistent worker gracefully, then apply after its current job has
finished. The apply command takes the **same outer lock as the cron launcher**;
it refuses to run while that lock is held by the worker. While cleanup holds
the lock, cron cannot launch a replacement worker.

```bash
for worker_pid in $(pgrep -f '^/home/ark/miniconda3/bin/python3 /home/ark/MAB/bin/gtotree-local/gtotree_worker[.]py'); do
    kill -TERM "$worker_pid"
done

/home/ark/miniconda3/bin/conda run --no-capture-output -n gtotree python cleanup_public_results.py --apply
```

If the command reports that the worker is still running, let its current job
finish before retrying. Do not delete lock files or force-kill an analysis.
Cron restarts the updated worker at a subsequent five-minute tick after cleanup
releases the lock. Frontend deployment remains the existing automatic Amplify
deployment after merging its PR; no Amplify CLI commands are needed.

For each result, cleanup first prepares replacements and backs up **all affected
original objects** locally. It then publishes sanitized data/report/metadata and
removes diagnostic objects, including logs and the old raw archive. No analyses
are rerun. Backup files and their `index.json` have owner-only permissions under
`~/.mab-private/gtotree-results/<timestamp>/<result-code>/`. The index maps each
original S3 key to its backup filename. `--backup-dir` can select another private
location. If backup fails, that result's remote objects are left untouched. Check
the command's errors before considering cleanup complete.

This keeps the existing result-code sharing model; scientific results remain
publicly retrievable. It does not add account-based authorization or erase copies
of files someone already downloaded.

## Refresh report appearance only

Report HTML is saved when it is generated. Pulling the generator updates new
reports; it does not change HTML already stored in S3. To remove the old nested
logo, run heading and summary cards from existing reports, stop the persistent
worker as above, then run:

```bash
/home/ark/miniconda3/bin/conda run --no-capture-output -n gtotree python cleanup_public_results.py --reports-only --apply
```

This uses private backups and the same cron lock. It rebuilds only report HTML
from the scientific data embedded in each report; it does not rerun analyses or
change trees, alignments, other downloads, result codes or statuses. Use `--job`
to refresh a single result. Cron resumes the worker at its next five-minute tick.

