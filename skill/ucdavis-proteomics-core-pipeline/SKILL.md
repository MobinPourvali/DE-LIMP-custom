---
name: ucdavis-proteomics-core-pipeline
description: >
  Run an end-to-end proteomics search + differential expression analysis from raw
  mass-spec data. Use this whenever the user wants to "analyze my proteomics data",
  "search these raw files", "run my DIA/DDA data", "find differentially expressed
  proteins", "process this timsTOF/Astral/Orbitrap run", or points at a folder of
  .raw / .d / .mzML files and asks what's in it. Detects acquisition + instrument,
  fetches a Brett-validated workflow from the DE-LIMP repo, downloads the pinned
  search engine, runs DIA-NN (DIA) or Sage (DDA), then limpa/limma DE — with full
  provenance back to the validated workflow. Also use it to "write the LC-MS methods
  section" / "generate a publication-ready methods section with the instrument grant
  acknowledgment" from facility raw data (UC Davis Proteomics Core).
---

# Proteomics Pipeline

Take a user from **raw MS files → differentially expressed proteins** using a
**validated** workflow, never ad-hoc parameters. The validated parameters live in
the public `bsphinney/DE-LIMP` repo under `workflows/`; this skill is the
orchestrator that reads them. **You hold orchestration logic; the repo holds the
science.** When they disagree, the repo wins.

All scripts are in `scripts/` next to this file. Reference detail is in
`references/` — read the relevant one before the step it covers; keep this file as
the spine.

## Golden rules (do not violate)

1. **Confirm before committing compute.** A search is multi-hour. Always show the
   auto-picked workflow and the organism/design and get an explicit "go" before
   running the engine.
2. **Never fabricate parameters.** If a value isn't in the workflow bundle or given
   by the user, say so — don't invent an FDR, organism, or instrument. (DE-LIMP
   architectural rule #2.)
3. **Never run computationally intensive work on a cluster login/head node — EVER.**
   On any cluster (HIVE/SLURM, or any other scheduler), **every** heavy step — the
   search, and any large DE / figure / conversion (e.g. msconvert) — must run as a
   **scheduled job** (`run_search.py --sbatch` → `sbatch`), not inline on the login
   node. Login nodes are shared; running compute there gets the user flagged/killed.
   The **only** things allowed on the login node are tiny orchestration commands —
   submitting jobs, `squeue`/`sacct` polling (`watch_run.sh`), and small file moves.
   If you're unsure whether a step is heavy, submit it as a job. (No SLURM but a big
   dataset locally? Warn the user it may be slow rather than hammering a shared host.)
4. **Organism is a hard filter** (it defines the FASTA). Instrument is only a
   tiebreaker. Acquisition is auto-detected and confirmed.
5. **Every run must be completely reproducible.** This is not optional. As you go,
   append every command you run (verbatim, with all arguments) to a `commands.log`.
   At the end you MUST produce a reproducibility bundle (step 9) that captures the
   pinned registry commit, exact tool + package versions, all parameters, input
   and output checksums, and a runnable `reproduce.sh`. Pin the registry to the
   commit SHA returned by `fetch_workflows.py` — never describe a result without
   the bundle that lets someone re-derive it. (DE-LIMP architectural rules #1, #4.)
6. **Assume nothing about the user's environment — this skill runs for anyone.** Most
   users are **not** UC Davis Core and **not** on HIVE: they're on macOS, Windows
   (**WSL2** recommended, but **native Windows works too** — the engines ship Windows
   builds), or generic Linux. HIVE / `/quobyte/proteomics-grp` is a Core **fast-path,
   never a requirement**. For **every** program you have the user run, give a **public
   source** anyone can obtain (download URL / `pip`/install command) and a path that
   works on *their* OS — never hardcode an internal `/quobyte/...` path without its
   public fallback. → `references/search-engines.md` ("Public program sources") +
   `references/install.md`.

## Audience: assume nothing is installed

The user may be a biologist on a fresh laptop with no R, no Python packages, no
Docker. **Do not ask them to install things by hand.** Run `setup.sh` — it
installs everything that can be installed without admin rights into one
self-contained conda env. Only one thing ever needs the user's hands (Docker
Desktop, and only for DIA-NN on macOS), and `setup.sh`/`build_diann_docker.sh`
print the exact step. Never dump a generic "please install R/limma/..." message;
the scripts handle it.

## Flow

### 0a. Access & environment — ASK THIS FIRST
**Claude Code runs locally on the user's machine.** Heavy steps either run locally or
are driven on **HIVE over SSH with the user's private key** — they do NOT run Claude
Code on HIVE. Ask two questions, then verify and pick a mode:

1. **"Do you have access to UC Davis HIVE (an account + an SSH private key)?"**
   If yes, **ask where their private key is** (e.g. `~/.ssh/id_ed25519`).
2. **"Are you a member of the UC Davis Proteomics Core?"**

Verify (don't just trust the answers):
```
bash scripts/check_access.sh <hive_user> <private_key_path>
```
Read `recommended_mode` + `facility_software_available`, then:

- **HIVE = yes → `hive_remote`:** drive HIVE over SSH from the local Claude Code
  (`export HIVE_USER=… HIVE_KEY=…`; use `scripts/hive_exec.sh '<cmd>'`). The search
  runs as a **SLURM job** (`run_search.py --sbatch` → `hive_exec.sh 'sbatch job.sh'`),
  never the login node. HIVE gives **compute**; the Core software is separate (next).
- **Core member = yes (with HIVE) → reuse the installed software** in
  `/quobyte/proteomics-grp/`: `acquire_tools.sh` finds the DIA-NN `.sif`,
  `fetch_fasta.py --hive` reuses pre-staged FASTAs. No rebuilding.
- **HIVE = yes but NOT Core → rebuild the toolchain in your own HIVE home** (you can't
  read the group dir). Follow `references/access.md` → "Rebuild on HIVE" (run
  `setup.sh` + `acquire_tools.sh` + `fetch_fasta.py` on HIVE via `hive_exec.sh`). **List
  these steps for the user — do not guess.**
- **Core = yes but NO HIVE →** the Core software is on HIVE, unreachable without HIVE
  access → fall back to **local**; suggest requesting a HIVE account.
- **Both = no → `local`:** `setup.sh` installs the toolchain on the user's machine
  (no admin); public engines (DIA-NN Academia, Sage). Self-contained, works fully.

Tell the user which mode you're using and why. → full runbook: `references/access.md`.

### 0. One-time setup (install everything that's missing)
```
bash scripts/setup.sh
source ~/.proteomics-pipeline/activate.sh        # puts R, python, sage on PATH
```
`setup.sh` installs micromamba (no admin), then a conda env with R + limpa +
limma + arrow + Sage + Python/pyarrow, and writes `~/.proteomics-pipeline/setup.json`.
Sourcing `activate.sh` makes every later step use those interpreters — **source it
in the same shell before running anything below** (or prefix later commands).

Read `setup.json` and **gate on `ready_for`**:
- `ready_for.de` false → DE can't run; re-run `setup.sh` and report any `notes`.
- `ready_for.dia` false → on macOS this means Docker is missing. Run
  `bash scripts/build_diann_docker.sh` and relay its instructions (install Docker
  Desktop, open it once), then continue. Don't silently fall back to a DDA engine.
- `ready_for.dda` false → Sage/R not ready, or (macOS) msconvert is unavailable for
  `.d`/`.raw` → mzML; tell the user and, if their data is DIA, route to DIA-NN.

This step is idempotent — on a machine that's already set up it just verifies and
returns in seconds. → detail: `references/install.md`.

### 0b. Detect the environment
```
bash scripts/detect_env.sh > /tmp/env.json
```
Read `platform_class` (mac|hpc|linux), `container_runtime`, `uc_davis_hive`. This
decides how tools are acquired and whether to submit via SLURM.
→ detail: `references/environment.md`.

### 0c. Is this a RESUME? — CHECK THIS BEFORE STARTING ANYTHING
A cluster search runs for hours and the user **will** close their laptop. The SLURM
jobs survive that; the conversation does not. So before doing any work, check whether
an earlier session already submitted something:
```
python3 scripts/checkpoint.py find --base <where results live>   # e.g. ~ or the project dir
python3 scripts/checkpoint.py status --session <session dir>     # live sacct re-query
```
`status` returns each stage's real state, whether the expected output file exists, and
`next_commands` — the exact commands still outstanding. If a run is found:

- **Still RUNNING** → say so, give the job ids, and offer to wait (`watch_run.sh`) rather
  than resubmitting. **Never resubmit a search that is already running** — it wastes hours
  of cluster time and the user's budget.
- **COMPLETED and the output exists** → skip the search entirely and continue from
  `next_commands` (usually DE onward).
- **FAILED or stalled** → apply the `watch_run.sh` playbook (step 7b) and resubmit only the
  incomplete steps; finished `.quant` files and the predicted library are reused.

`submit.sh` writes `RECOVERY.md` + `.recovery.json` into the session directory
automatically, so this works even if the previous session ended mid-sentence. Tell the
user, in plain words, that they can close their terminal and that the jobs will keep
going — they do not need `tmux` or `screen` for this.

### 1. Locate the raw files
Ask the user for a directory or file list if not given. Recognized: `.d` (Bruker),
`.raw` (Thermo), `.mzML[.gz]`, `.wiff` (convert first). Glob to a concrete list.

### 1b. Check for a prior analysis of this dataset
```
python3 scripts/session.py find-prior --raw /path/to/*.d
```
If a match is returned, this is a **re-analysis** of an existing dataset — note the
prior session dir; you'll pass it to `session.py init --reanalysis-of` (step 3b) so
the run nests under `<prior>/reanalysis/` and gets a `DIFFERENCES.md`, and you'll
run the Comparator at the end (step 12). If no match, it's a fresh analysis.

### 2. Detect acquisition + instrument
```
python3 scripts/detect_acquisition.py /path/to/*.d /path/to/*.raw
```
Returns per-file `acquisition` (DIA/DDA/unknown) + `confidence`, plus an overall
`instrument`. **If `needs_confirmation` is true, ask the user** before continuing
— mixed/unknown/low-confidence must not silently pick an engine.
→ detail: `references/search-engines.md`.

### 3. Ask organism + experimental design (auto-map conditions)
- **Organism** cannot be detected — ask, and resolve to a UniProt **taxid**
  (human = 9606). This is required and authoritative.
- **Conditions:** ask the user to either *tell you* the conditions in plain words
  ("the first three are control, the rest treated") **or** *upload a file* (any
  CSV/TSV with a sample column and a group column, however named). Don't make them
  hand-fill a template.

First get the real run names the conditions must map to:
```
python3 scripts/collect_conditions.py --list-runs --from-dir /path/to/raw --glob '*.d'
```
Then map their conditions onto those runs — the script does the fuzzy filename
matching so you don't guess:
```
# (a) user uploaded a file:
python3 scripts/collect_conditions.py --map conditions.csv \
    --from-dir /path/to/raw --glob '*.d' --from-file <their_file>
# (b) user described it in words: turn it into intent, then ground it against the runs:
python3 scripts/collect_conditions.py --map conditions.csv \
    --from-dir /path/to/raw --glob '*.d' \
    --mapping-json '{"groups": {"control": ["..."], "treated": ["..."]}}'
```
Read the returned `ambiguities` and **confirm every one with the user** before
proceeding — `unassigned_runs` (a raw file no condition matched), `conflicting_runs`
(a file matched to two groups), `unmatched_identifiers` (the user named something
with no matching file), `multi_match_identifiers` (one label hit several files —
usually fine, e.g. a replicate prefix), and `singleton_groups` (<2 replicates → no
within-group variance). Do **not** start a search while any run is unassigned or
conflicting. Finalize the CSV, then validate it:
```
python3 scripts/collect_conditions.py --validate conditions.csv --against report.parquet
```
(Fallback: if the user has nothing yet, `--emit-template` writes a blank
File.Name,Group sheet for them to fill.) → detail: `references/conditions.md`.

### 3b. Create the analysis session (organize all files)
**Ask the user where the results should go** — two choices:
- **In the folder with their raw data** (default, recommended — keeps results next
  to the data): pass `--raw` and omit `--base`.
- **A central location** (e.g. their Documents folder, or a path they give): pass
  `--base <that path>` (results land under `<path>/sessions/`).

Then scaffold the session and route **everything** into it:
```
# default — results live with the raw data:
python3 scripts/session.py init --name "<short study name>" --raw /path/to/*.d \
    [--reanalysis-of <prior session dir from step 1b>]
# or central, if the user chose one:
python3 scripts/session.py init --name "<short study name>" --raw /path/to/*.d \
    --base ~/Documents/DataAnalysis
```
The output's `placement` tells you which was used. **Use the printed `paths` map for
every later step** — put
`conditions.csv` and the FASTA in `paths.input_dir`, search output in
`paths.search_out`, DE results in `paths.de_dir` (= `output/tables`), the
reproducibility bundle in `paths.repro_dir`, the report in
`paths.analysis_report`, and append commands to `paths.commands_log`. Raw files
are recorded (in `input/raw_files.txt`), never copied. With `--reanalysis-of`, the
session nests under `<prior>/reanalysis/<date>_<name>/`.
→ detail: `references/outputs.md`.

### 4. Match a validated workflow (then CONFIRM)
```
python3 scripts/fetch_workflows.py match \
    --acquisition DIA --organism-taxid 9606 --instrument "Orbitrap Astral"
```
Hard-filters on acquisition+taxid, scores instrument, returns `selected` +
`candidates` + `needs_menu`. **Present `selected` to the user** — name, engine +
pinned version, FASTA, DE method, and the `validated` provenance — and get
confirmation. If `needs_menu` is true (no match / tie / no instrument info),
present `candidates` as a menu instead of auto-proceeding. If zero candidates,
tell the user no validated workflow exists for this acquisition+organism and stop
(offer to add one to the registry — see `workflows/README.md`).

**Record the `registry.commit` SHA from the match output** — it pins the exact
validated-params version for reproducibility. On confirm, pull the bundle's params
at that pinned commit:
```
python3 scripts/fetch_workflows.py pull --id <id> --ref <registry.commit> --dest ./wf
```
This writes the engine params file + `workflow.manifest.json` (engine, version,
fasta spec, de spec, and the pinned `registry` commit). Pulling at the commit (not
`main`) guarantees a future re-run gets byte-identical parameters.

### 5. Acquire the pinned engine
Honor the bundle's exact version — not "latest":
```
PIN_ENGINE=diann PIN_VERSION=2.6.0 bash scripts/acquire_tools.sh <platform_class>
```
Reads/writes `~/.proteomics-pipeline/tools/tools.json`. On HIVE it reuses the
existing `.sif`; on mac it uses Docker for DIA-NN. **Read `tools.json` `notes`** —
license gates (FragPipe) and missing-runtime warnings surface there.

### 6. Resolve the FASTA
```
python3 scripts/fetch_fasta.py --proteome UP000005640 --add-contaminants \
    --out ./search.fasta [--hive]
```
Pass `--hive` when `uc_davis_hive` is true to reuse pre-staged proteomes +
contaminants under `/quobyte/proteomics-grp/MRS/` instead of downloading.
→ detail: `references/environment.md` ("FASTA").

### 6b. Estimate search parameters from the data type
Search parameters are **derived from what the data is**, not hand-maintained.
If the bundle has a `params_file` (a validated SOP config), use it verbatim.
Otherwise (the default — `estimate_params: true`), generate them:
```
python3 scripts/estimate_params.py --engine <diann|sage> \
    --acquisition <DIA|DDA> --instrument "<detected instrument>" \
    --var-mods "<bundle var_mods>" --overrides '<bundle param_overrides as JSON>' \
    --out ./wf/params.<cfg|json>
```
The estimator keys mass tolerances on the instrument class from DIA-NN's
known-good table (Astral → MS1 4/MS2 10 ppm; timsTOF → 15 ppm; unidentified →
automatic calibration), sets DIA/DDA window mode from acquisition, and uses
standard trypsin/LFQ defaults for the rest. **It prints a `rationale` tagging
every value's provenance** — surface this to the user (and it flows into the
methods text), so a derived default is never mistaken for a confirmed setting.
Use the resulting file as `--params` below. → detail: `references/parameters.md`.

### 7. Run the search
```
python3 scripts/run_search.py --tools ~/.proteomics-pipeline/tools/tools.json \
    --bundle ./wf/workflow.manifest.json --params ./wf/params.<cfg|json> \
    --fasta ./search.fasta --out ./search_out --files /path/to/*.d --threads 16
```
- DIA → DIA-NN; DDA → Sage; FragPipe only if the bundle names it or the user asks.
- **FragPipe + diaTracer (timsTOF dia-PASEF, spectrum-centric DIA).** A second DIA route,
  used when the bundle names it or the user asks for it. diaTracer traces 3-D features and
  writes pseudo-MS/MS spectra, MSFragger searches those to build a library from the data
  itself, then **DIA-NN quantifies against it**. Worth reaching for as a cross-check on the
  DIA-NN route, or when the search needs semi-tryptic / nonspecific / open-modification
  behaviour that library-free DIA-NN handles poorly.
  - **⚠ Licensing — get this right, it is easy to state wrongly.** The gate is
    **MSFragger / IonQuant / diaTracer** (one shared academic licence; commercial users
    license via fragmatics.com). The DIA-NN that FragPipe bundles is explicitly *"available
    for both academic and commercial use within FragPipe"*, so DIA-NN is **not** the
    restriction here. A commercial user is not locked out — they must license the three
    gated tools. **Never tell them the route is closed to them.**
  - **⚠ The clause that matters for a core facility:** the academic licence forbids *"any
    work product, report or deliverable by LICENSEE, for a commercial or for-profit
    organization"* — which reads on fee-for-service work for an industry client. If asked
    whether their work qualifies, send them to info@fragmatics.com and their licensing
    office. Do **not** rule on it yourself. → `references/search-engines.md`.
  - Output is DIA-NN's own, in `<workdir>/dia-quant-output/`. With FragPipe's bundled
    DIA-NN 1.8.2 Beta 8 that is **`report.tsv` — there is no `report.parquet`** unless
    someone configured DIA-NN 2.x. `adapt_fragpipe` handles either.
  - **The pseudo-MS/MS mzML is written next to the `.d` file, not into the workdir**, so
    the raw-data directory must be **writable** — this fails on a read-only mount. And a
    manifest data type that isn't one of `DIA`/`GPF-DIA`/`DIA-Quant`/`DIA-Lib`/`DDA+`/`DDA`
    is **silently treated as DDA**, quietly turning a DIA run into a DDA one. Use `DIA`.
  - FragPipe 24.0 bundles **DIA-NN 1.8.2_beta_8**, older than the 2.x the `diann_*`
    workflows pin. Expect different numbers from that alone; say so rather than
    attributing every difference to the spectrum-centric approach.
  - Budget roughly **20 min per file for the diaTracer conversion alone** (paper: 34 files,
    32 cores, ~19.4 min each), before MSFragger and DIA-NN. The mzML are **reusable**, so a
    re-search with different parameters skips it. Needs Java 11+, MSFragger 4.4+,
    IonQuant 1.11.18+, diaTracer 2.2.1+. The three gated jars are **not** in the FragPipe
    zip and **cannot be scripted** — a human accepts the licence once and mints a key, then
    point `FRAGPIPE_TOOLS_FOLDER` at them. → `references/search-engines.md`.
- **⚠ DIA-NN licensing (ask before a DIA run):** DIA-NN's free "Academia" build is
  licensed for **academic / non-profit use only**. If the user is in a **commercial /
  non-academic** setting, do NOT use DIA-NN — **offer AlphaDIA** (`--engine alphadia`,
  Apache-2.0, commercial-OK; the open-source DIA alternative) or Sage (MIT, also
  commercial-OK). A quick question: *"Is this academic/non-profit, or commercial?"* —
  commercial → AlphaDIA. (AlphaDIA is deep-learning based; a GPU is strongly
  recommended — best on a HIVE GPU node or a GPU workstation.)
- **⚠ AlphaDIA limitation — timsTOF whole-proteome directDIA:** AlphaDIA does **not**
  work on **timsTOF** `.d` files for **whole-proteome library-free (directDIA)**
  searches. So for timsTOF + whole-proteome + library-free, AlphaDIA is NOT an option.
  In that case: academic → DIA-NN; commercial → use a **spectral-library** approach
  (predicted/empirical library) rather than directDIA, or Sage — and tell the user the
  constraint. (Non-timsTOF instruments, or library-based timsTOF runs, are fine for AlphaDIA.)
- **DIA-NN search settings:** always use the **estimated cfg** (step 6b) — it encodes
  DIA-NN's official recommended per-instrument mass tolerances (timsTOF → MS1/MS2 15
  ppm; Orbitrap Astral → 4/10 ppm; Orbitrap by resolution 240k→4, 120k→7, 60k→10,
  30k→15). Don't override these unless the user gives a validated SOP value.
- **DIA-NN parallel — AUTOMATIC above 5 files.** `run_search.py` routes to the **5-step
  parallel chain** by itself whenever the run is DIA-NN with **more than 5 files** on a
  machine that has SLURM and a `--cfg` with **fixed** mass accuracy. You don't decide
  this; it prints `parallel routing: YES|no -- <reason>` and records the reason in
  `search_provenance.json`. This matches facility practice in DE-LIMP: per-file passes
  run as a SLURM array instead of one long single-node job.
  - It falls back to a single search, with the reason printed, when there's **no SLURM**
    (the chain is job arrays — nothing to fall back to) or when **mass accuracy is on
    auto**. Mass accuracy is the one you can fix: steps 3/5 reuse the `.quant` files, so
    auto-calibration would differ between passes. Re-run `estimate_params.py` with the
    **real instrument** (step 6b) to pin it and parallel turns itself on.
  - Override either way with `--no-parallel` (force one job) or `--parallel-threshold N`.
  - It **generates** the chain but does not submit it. Submit `<out>/submit.sh` (over
    `hive_exec.sh` on HIVE), then watch the **step-5** job — that's the one that writes
    `report.parquet`. → `references/diann_parallel.md`.
- **On `hpc`:** add `--sbatch job.sh`, then `sbatch job.sh` (over `hive_exec.sh` for a
  remote HIVE run). Re-run with `--adapt-only` afterward for Sage/FragPipe/AlphaDIA to
  build `report.parquet`.
- Output is normalized to the **DE contract**: a DIA-NN-shaped `report.parquet`.
→ detail: `references/search-engines.md`.

### 7b. Watch the run — MANDATORY, auto-correct errors
**The moment a search starts, watch it — never leave a multi-hour search
unmonitored.** Poll with `watch_run.sh` in a loop until it finishes:
```
# local: python3 scripts/run_search.py ... run_in_background, then loop:
bash scripts/watch_run.sh --log <search log> --out <search out dir> --poll <n>
# SLURM on HIVE:
bash scripts/watch_run.sh --slurm <jobid> --log <job log> --out <out> --poll <n> --hive
```
**Always pass `--out` and increment `--poll` each time.** `--out` is what turns "still
running" into real progress: the chain drops a `.quant` per finished file, so the watcher
reports the actual step and file count rather than guessing. `--poll` rotates the note.

**Report progress to the user on every poll — a search is hours long and silence reads as
a hang.** Give them one short line from `progress.summary` plus `doing`, and pass along
`note` (with `note_source`) so there's something to read while waiting. For example:

> Step 2 of 5 — searching each file against the predicted library. 34 of 66 files done
> (52%). *While you wait: decoy sequences exist to be wrong on purpose. How often a search
> prefers a decoy is what calibrates the false discovery rate. (Elias & Gygi 2007)*

Keep it to a couple of lines; don't dump the JSON. Never present the note as a finding
about **their** data — it is general background, and it must not be mistaken for a result.
It returns `{state, done, failed, stalled, error_class, fix}`. While `done` is false,
keep polling (sleep ~60s between polls; for long SLURM jobs use a scheduled wake-up).
Pass `--log` so it can also catch **stalls** — a job stuck in `RUNNING` with its log
frozen (a hung file; `sacct` will show RUNNING forever). **If `failed` OR `stalled` is
true, apply the `fix` and resubmit AUTONOMOUSLY — do not wait for the user.** Starting
the search needed confirmation (golden rule #1); recovering a run that's already in
flight does not. Common auto-fixes (→ `references/watcher.md` playbook):
- **stalled/hung file:** `scancel` that task, retry once on a fresh node; if it stalls
  again, **drop that one file** and continue (step 4 of the parallel chain auto-skips a
  file with no `.quant`); note the dropped file in **Data Quality Notes**.
- **failed step in the 5-step chain:** resubmit from the earliest incomplete step,
  **reusing completed `.quant` and `step1.predicted.speclib`** — never restart the whole
  cohort; broken `afterok` after a killed task → resubmit steps 3→4→5 fresh.
- OOM → raise `--mem`; timeout → raise `--time`; missing temp dir → `mkdir -p` first.

**Every search — single-shot or parallel array — must run under this watch loop.** Stop
and tell the user only after 2 failed auto-fixes of the same class (dropping 1 pathological
file out of many is success, not a failure). Report what you recovered. Only proceed to DE
once `COMPLETED` and `report.parquet` exists. → detail: `references/watcher.md`.

### 8. Differential expression
```
Rscript scripts/run_de.R --input ./search_out/report.parquet \
    --metadata conditions.csv --method <dpc|maxlfq> --outdir ./de_results
```
Use the `de.method` from the bundle (`dpc` for DIA-NN/limpa, `maxlfq` for
Sage/FragPipe). Writes `DE_<method>_<contrast>.csv` + `Expression_Matrix.csv` +
`methods.txt` + `sessionInfo.txt` + `de_provenance.json` (exact R package versions).
→ detail: `references/de-analysis.md`.

### 8b. Generate figures
```
Rscript scripts/make_figures.R --de-dir ./de_results --conditions ./conditions.csv \
    --outdir ./figures --adjp 0.05 --logfc 1
```
Produces publication-quality volcano (per contrast), PCA, heatmap of top proteins,
p-value distributions, and a per-sample protein-count QC plot, plus `figures.json`
(captions). These get embedded in the report.

**Significance is `adj.P.Val` alone — there is no fold-change cutoff anywhere in this
pipeline.** `--logfc` only positions a labelled reference line on the volcano. Never
re-impose it when you count, rank, or describe significant proteins, and never justify a
result by "it's ≥2-fold". The moderated t-test tested H0: log2FC = 0, so that is the only
claim the FDR covers; filtering the list on the observed fold change afterwards would
report a stronger claim than the error rate supports, because an observed fold change is
a point estimate with no uncertainty attached. (Testing a fold-change threshold properly
needs `limma::treat()`'s interval null — McCarthy & Smyth 2009 — not a post-hoc cut.) You
may report how many significant proteins also exceed the reference line, as a description
of effect size.

### 8c. Audit the results for common mistakes (surface every issue)
```
python3 scripts/audit_results.py --out AUDIT.md --conditions ./conditions.csv \
    --de-dir ./de_results --acquisition-json /tmp/acq.json --adjp 0.05 --logfc 1
```
Checks for the classic new-user pitfalls: too few/no replicates, imbalanced or
**confounded** design, **mixed acquisition or mixed instruments** in one analysis,
suspiciously low ID depth, very high missingness, contaminant dominance, and DE
results that are too-empty or implausibly-large (batch/normalization artefacts).
**Surface every `WARN` to the user, and STOP on any `FAIL`** (e.g. a group with no
replicate, or a batch confounded with the biology) until they resolve it — don't
let a new user over-interpret a broken design. The findings also go into the
report's "Audit & caveats" section. → detail: `references/audit.md`.

Then run the **biological sample-quality** check (contamination that mimics biology —
hemolysis, tissue cross-contamination, skin) — and **read `references/anomaly-checks.md`
and walk its decision tree** (the judgment branches `audit_results.py` can't automate:
intensity-distribution outliers, replicate/condition ID overlap, mass/RT/CCS shifts,
p-value-distribution shape, largest-FC-in-low-abundance, log2FC>5 artefacts, PCA sample
swaps, GSEA background; plus XL-MS / phospho / non-model branches):
```
python3 scripts/sample_quality.py --matrix ./de_results/Expression_Matrix.csv \
    --conditions ./conditions.csv --report ./search_out/report.parquet --out SAMPLE_QUALITY.md
```
**If a contamination panel is confounded with a group, STOP** — DE may be contamination,
not biology, and protein-level filtering will not fix it (→ `references/anomaly-checks.md`).
**If the sample IS keratin** (nail/hair/wool/skin/feather), pass **`--keratin-sample`** to
`sample_quality.py` **and** `audit_results.py` — keratin is the analyte there, not a
contaminant, and must not be flagged.
Flag every anomaly in the 5-part format (*what · where · why · likely cause · fix*); the
report **always** gets a **Data Quality Notes** section, even if "nothing anomalous observed."

### 8d. Conditional expert review (spawn only on triggers)
Most analyses don't need a panel. If **any** trigger fires — a statistical claim goes to a
collaborator; any `log2FC > 5` or `p < 1e-10`; a new organism/instrument/reagent; cross-tool
discordance; or the user asks — spawn **three reviewers in parallel** (one message, three
`Agent` calls; each gets the report + data paths + context): **proteomics**, **biology**
("could contamination explain this?"), **statistics**. Consolidate their findings (dedupe,
sort by severity) into a `## Expert Review Notes` section and surface every **critical**
issue to the user before finalizing. Skip for format conversion / FASTA prep / exploratory
looks with no formal report. → detail: `references/anomaly-checks.md`.

### 9. Analyze the data (you write the report)
Generate the analysis brief (it lists the figures to embed), then **do the analysis
yourself** — you are the consultant the brief addresses:
```
python3 scripts/analysis_prompt.py --out ANALYSIS_PROMPT.md \
  --de-dir ./de_results --report ./search_out/report.parquet \
  --conditions ./conditions.csv --figures-dir ./figures [--qc ./QC_Metrics.csv] \
  --engine <engine> --acquisition <DIA|DDA> --instrument "<name>" \
  --workflow-manifest ./wf/workflow.manifest.json
```
Then **read `ANALYSIS_PROMPT.md` and every data file + figure it lists, and write a
complete `AI_Analysis_Report.md`** with ALL its OUTPUT sections (Overview, QC, Key
Findings Per Comparison, Cross-Comparison Biomarkers, High-Confidence Biomarkers,
Pathway/GSEA if present, Biological Interpretation, How This Analysis Works, Methods
& Reproducibility) **plus an "Audit & caveats" section from `AUDIT.md`**, a **Data
Quality Notes** section (the 5-part anomalies from steps 8c/8d + `SAMPLE_QUALITY.md`,
even if "nothing anomalous observed"), and — if a review panel ran (8d) — an **Expert
Review Notes** section. **Embed
every figure** (`![caption](figures/<file>.png)`) with an expert interpretation of
what it shows for THIS data. Compute significant proteins, up/down splits,
cross-comparison overlaps, and lowest-CV proteins directly from the CSVs — cite
specific proteins, never fabricate. Make it thorough and expert, like the DE-LIMP
AI export. The brief takes its pipeline description from `de_provenance.json`, so
the report stays correct for whichever engine/method ran.

Then **also save the report as a Word document** (both formats are required):
```
python3 scripts/to_docx.py --in <session>/output/AI_Analysis_Report.md \
    --out <session>/output/AI_Analysis_Report.docx
```
→ detail: `references/analysis.md`.

### 9d. Publication-ready Methods section + acknowledgment
Generate a drop-in LC-MS/MS Methods section straight from the facility raw data,
with the correct UC Davis Proteomics Core instrument-grant acknowledgment:
```
python3 scripts/make_methods.py --raw /path/to/*.d \
    --out <session>/output/methods.md --de-dir <session>/output/tables
python3 scripts/to_docx.py --in <session>/output/methods.md \
    --out <session>/output/methods.docx
```
It extracts the acquisition parameters from the raw metadata (Bruker `.d`
`analysis.tdf`; Thermo by facility filename prefix), fills the rest from facility
defaults **tagged `[facility default — confirm]`** (the LC column defaults to a
PepSep C18 10 cm × 150 µm, 1.5 µm — override with `--lc-column`), builds a
parameter table showing the source of each value, and appends the instrument's
grant acknowledgment (Fusion Lumos → S10OD021801; Exploris 480 → S10OD026918-01A1;
timsTOF → Dr. Neil Hunter / HHMI). **Verify the draft against the params table and
polish the prose; keep the acknowledgment exact** (confirm wording at the source
URL). This can also be run standalone — just point `--raw` at facility data, no
search/DE needed. → detail: `references/methods.md`.

### 10. Reproducibility bundle (mandatory)
Assemble the bundle that makes the whole analysis reproducible:
```
python3 scripts/provenance.py --outdir ./reproducibility \
  --workflow-manifest ./wf/workflow.manifest.json \
  --setup-json ~/.proteomics-pipeline/setup.json \
  --tools-json ~/.proteomics-pipeline/tools/tools.json \
  --params ./wf/<params_file> --conditions ./conditions.csv \
  --fasta ./search.fasta --fasta-info '<json from fetch_fasta>' \
  --raw /path/to/*.d --report ./search_out/report.parquet --de-dir ./de_results \
  --engine <engine> --de-method <dpc|maxlfq> --contrasts "<...>" \
  --q-cutoff 0.01 --logfc 1.0 --adjp 0.05 \
  --acquisition <DIA|DDA> --organism-taxid <taxid> --instrument "<name>" \
  --commands ./commands.log --timestamp "$(date -u +%FT%TZ)"
```
This writes `reproducibility/` with `run_manifest.json`, `reproduce.sh`,
`REPRODUCE.md`, the conda lock + pip freeze + R sessionInfo + tool versions, a
`skill.txt` recording **which skill produced this and how it was installed**, copies
of the params and conditions, and sha256 checksums of inputs/outputs. This step is
**mandatory and must always run** — code + versions are not optional. `provenance.py`
auto-discovers the env and Rscript even if `--setup-json` is absent, so versions are
always captured; still, **check the returned `skipped` count** and, if the conda
lock / R sessionInfo / checksums were skipped, fix the cause and re-run.

### 11. Output-files report
Catalog everything the run produced so the user knows what each file is:
```
python3 scripts/make_report.py --out OUTPUT_FILES.md \
  --search-out ./search_out --de-dir ./de_results --repro ./reproducibility \
  --extra ./conditions.csv ./search.fasta ./wf ./figures ./AUDIT.md \
          ./AI_Analysis_Report.md ./AI_Analysis_Report.docx
```
`OUTPUT_FILES.md` lists every file (figures, audit, search/DE outputs, the bundle)
with its size and a plain-language description, grouped by purpose, and flags
anything unrecognized (never silently omitted).

### 11a. Comparing two searches of the same files (tool bake-off)
When the same raw files have been searched more than once — DIA-NN library-free vs
FragPipe/diaTracer, two engine versions, two parameter sets — compare the **searches**,
not just the DE tables:
```
python3 scripts/compare_searches.py --out <session>/output/search_comparison --q 0.01 \
  --search "DIA-NN:<sessA>/output/search/report.parquet" \
  --search "diaTracer:<sessB>/output/search/report.parquet"
```
Writes `SEARCH_COMPARISON.md` + CSVs: protein groups and per-run depth, how many proteins
survive in **every** run, median CV across runs, the shared/unique split with Jaccard, and
the per-run correlation of log2 intensity on shared proteins. Unique-to-one-search proteins
are listed in full.

**Do not report "engine X found more proteins" as the verdict.** Depth, completeness and CV
have to be read together — a search can identify more while quantifying them less
reproducibly or in fewer runs. Two caveats to state whenever you present this: the CV is
computed across all runs with no group labels, so it only measures precision when those runs
are replicates; and the comparison assumes protein-group identifiers are comparable, which
fails if the searches used different FASTAs or inference rules. Then run 11b on the DE output
to see whether any of it changes which proteins come out significant.

### 11b. If this is a re-analysis: compare to the original
When step 1b found a prior analysis of the same dataset, compare the two with the
Comparator (a faithful port of DE-LIMP's Run Comparator):
```
Rscript scripts/compare_analyses.R --out <session>/output/comparison \
  --adjp 0.05 --logfc 0 \
  --analysis "Original:<prior>/output/tables" \
  --analysis "Reanalysis:<this session>/output/tables"
```
It writes `COMPARISON.md` + CSVs: protein-universe overlap, the 3×3 Up/Down/NS
concordance per shared contrast, direction concordance on co-significant proteins,
and logFC correlation. (`session.py finalize` already wrote `DIFFERENCES.md` —
what *settings* changed; the Comparator shows how the *results* changed.)

### 12. Finalize the session + report to the user
```
python3 scripts/session.py finalize --dir <session> --zip
```
This tidies `output/` (tables→`tables/`, figures→`figures/`), writes the session
`README.md` (and `DIFFERENCES.md` for a re-analysis), and zips the session for easy
sharing. Then summarize: workflow id + name, engine + **pinned version**,
**registry commit SHA**, FASTA source, DE method, per-contrast significant counts,
and the link to the validated workflow at that commit
(`https://github.com/bsphinney/DE-LIMP/tree/<commit>/<path>`). Point them at the
**session folder** and its `README.md`, then `AI_Analysis_Report.md` (the
interpretation), `OUTPUT_FILES.md` (what every file is), `tables/methods.txt`
verbatim (the Methods paragraph — don't paraphrase), `reproducibility/REPRODUCE.md`
(the re-run recipe), and — for a re-analysis — `DIFFERENCES.md` + the
`comparison/COMPARISON.md`.

## When something is missing
- Anything in the env (R, limpa, Sage, pyarrow) → re-run `setup.sh`; relay its
  `notes`. Never tell the biologist to "install R/limma/..." by hand.
- macOS + DIA data + no Docker → run `build_diann_docker.sh` and relay its exact
  steps (install Docker Desktop, open once). Don't silently switch to a DDA engine.
- macOS + DDA + Bruker/Thermo → msconvert is Linux-only; see `references/install.md`
  ("macOS + Sage"). Prefer DIA-NN if the data is DIA; else convert to mzML first.
- No validated workflow → stop, explain, offer to add one (`workflows/README.md`).
- FragPipe license (MSFragger/IonQuant) → surface the `tools.json` note and the fix.
- Acquisition/instrument unknown → ask the user; never guess into a multi-hour run.
