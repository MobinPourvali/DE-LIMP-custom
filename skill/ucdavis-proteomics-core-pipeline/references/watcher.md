# Run watcher — monitor searches and auto-correct errors

A search is long (often hours) and runs detached (a background process locally, a
SLURM job on HIVE). **Always watch every search — this is mandatory, not optional.**
The orchestrator must never start a search and walk away, and must **auto-correct
stalls and failures on its own, without waiting for the user** (the user has standing
authorization for this — see golden rule #1's compute confirmation is for *starting*
the search, not for recovering a run that's already approved and in flight). Surface
what you did afterwards; don't ask permission to unstick a job.

`watch_run.sh` does one poll + diagnosis (state, **stall**, error class, fix); loop it
until the run finishes, and on `failed` **or** `stalled` apply the fix and resubmit.

## Two failure modes it catches
1. **Hard failure** — `sacct` reports FAILED / TIMEOUT / OUT_OF_MEMORY / NODE_FAIL /
   CANCELLED. Matched against known error signatures → a fix.
2. **Stall** — the job stays **RUNNING** but its log stops advancing (a hung file, e.g. a
   pathological DIA-NN run whose RT-window/search phase never returns). `sacct`/`squeue`
   will happily show RUNNING forever. `watch_run.sh --log <log> [--stall-min N]` flags
   this when the log's mtime is older than N minutes (default 15) while RUNNING, emitting
   `stalled: true`. **This is the NA41-class failure the first nail cohort hit.**

## Loop pattern
```
# start the search (background locally, or sbatch on HIVE), capture the log path / job id
while true; do
  status=$(bash scripts/watch_run.sh --slurm <jobid> --log <log> --hive)   # or --log <log> locally
  done=$(echo "$status" | python3 -c 'import sys,json;print(json.load(sys.stdin)["done"])')
  [ "$done" = "True" ] && break
  sleep 60        # for long SLURM jobs, schedule a wake-up instead of busy-waiting
done
# if failed: read error_class + fix, apply it, resubmit, watch again.
```
The agent itself is the watcher — `watch_run.sh` is its eyes. Surface progress and
any auto-fix you applied to the user.

## Error classes → fixes (what `watch_run.sh` detects)
| error_class | signal | fix |
|---|---|---|
| `out_of_memory` | oom-kill, `std::bad_alloc`, OUT_OF_MEMORY | raise sbatch `--mem` (e.g. 64G→128G); DIA-NN: fewer threads; resubmit |
| `timeout` | TIMEOUT / "DUE TO TIME LIMIT" | raise `--time`, or split the run; resubmit |
| `diann_no_dotnet` | `dotnet: not found` | wrong DIA-NN container (no .NET → `.raw` silently skipped). Use the HIVE **native** build `build_<v>/diann-<v>/diann-linux` or a .NET image |
| `empty_results` | 0 proteins / no fragment ions | FASTA/organism mismatch, mass-accuracy, or wrong acquisition type — check inputs |
| `gpu` | CUDA / no kernel image | AlphaDIA needs a GPU: submit to a GPU node (`--gres=gpu:1`) or reduce batch size |
| `sage_no_mzml` | msconvert not found / needs mzML | convert `.d`/`.raw` → mzML first (Linux/HIVE), then re-run Sage |
| `disk` | Disk quota / No space left | free space or point `--out` elsewhere; resubmit |
| `missing_input` | No such file / fasta/raw not found | a path is wrong — re-check Windows→WSL/HIVE path translation; resubmit |
| `stalled` | job **RUNNING** but log frozen > `--stall-min` (default 15 min) | a hung file. `scancel` the task/job, **retry it once on a fresh node**; if it stalls again, **drop that file** and continue (see playbook), note it in Data Quality Notes |
| `dependency_failed` | job PENDING with reason **`DependencyNeverSatisfied`** (an upstream step failed) | the chain is dead but **sits PENDING forever** (never leaves the queue). `sacct -j <arrayjob>` to find the failed step, fix it, resubmit downstream steps reusing completed outputs |
| `unknown_failure` | job FAILED with no known signature | read the full log; diagnose; fix; resubmit (and add the new signature here) |

## SLURM specifics (HIVE)
- State comes from `sacct -j <jobid> -o State` (falls back to `squeue`). Terminal
  states: COMPLETED (done), FAILED / TIMEOUT / OUT_OF_MEMORY / NODE_FAIL / CANCELLED
  (done + failed).
- The job log is the sbatch `--output` file (the skill's `emit_sbatch` names it
  `<out>/<job>_<jobid>.log`).
- When resubmitting after a fix, edit the sbatch (`--mem`/`--time`/`--gres`) and note
  the change in `commands.log` so the reproducibility bundle records what happened.

## Auto-recovery playbook (apply autonomously; log every action to `commands.log`)
- **`stalled` (hung file):** `scancel <arrayjob>_<taskid>` (or the job). Retry that one file
  once on a fresh node. If it stalls again, **drop it** — it's pathological (often a file
  the facility itself re-ran). In the 5-step parallel chain, step 4 already auto-skips a
  file with no step-2 `.quant`; just remove it from step 3/step 5's `--f` list and resume.
- **Broken `afterok` after you killed/dropped an array task:** the downstream steps go
  `DependencyNeverSatisfied`. Cancel them and **resubmit steps 3→4→5 fresh, reusing the
  completed `.quant`** (don't re-run step 1/library-prediction — reuse `step1.predicted.speclib`).
- **`out_of_memory`:** raise `--mem` and resubmit only the failed step (reuse prior `.quant`).
- **Missing `--temp` dir** (`cannot find the temp folder`): the temp dirs must pre-exist —
  `mkdir -p <out>/quant_step2 <out>/quant_step4` and resubmit (the generator now does this).
- **Whole-cohort restart is rarely needed** — resume from the earliest incomplete step,
  reusing everything already computed.

## Don't
- Don't declare a run "done" on a non-COMPLETED state, and don't proceed to DE until
  `report.parquet` exists.
- Don't silently retry forever — after 2 failed auto-fix attempts of the same class,
  stop and tell the user what's wrong. (Dropping 1 pathological file out of many and
  continuing is *success*, not a failed retry.)
- Don't wait for the user to notice a stall/failure — you are the monitor; recover, then
  report what you did.
- **Don't use a passive "wait until the job leaves the queue" monitor for a chain.** If an
  upstream step fails, the downstream job sits **PENDING forever** (`DependencyNeverSatisfied`)
  and such a monitor hangs indefinitely and never alerts — a whole overnight chain can die
  silently. The monitor must poll **STATE** (incl. `squeue -o %r` reason) and exit/act on
  `failed` / `dependency_failed` / `stalled`, not just on queue-exit.
- **A job's SLURM State=COMPLETED does NOT mean the tool succeeded.** If the sbatch's last
  line is `echo ...` (exit 0), FragPipe/DIA-NN can crash yet the job shows COMPLETED. Verify
  the expected OUTPUT file exists (`report.parquet`, `combined_peptide.tsv`), not just the state.
