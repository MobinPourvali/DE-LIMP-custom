# Search engines reference

Acquisition detection, engine routing, and normalizing each engine's output to the
DE-input contract.

## Acquisition detection (`detect_acquisition.py`)

| format | how DIA/DDA is decided | instrument |
|---|---|---|
| Bruker `.d` | `analysis.tdf` SQLite: `DiaFrameMsMsInfo`/`DiaFrameMsMsWindowGroups` or `Frames.MsMsType==9` → DIA; `PasefFrameMsMsInfo`/`MsMsType==8` → DDA | `GlobalMetadata.InstrumentName` (e.g. "timsTOF Pro") |
| `.mzML[.gz]` | stream MS2 isolation windows: median width ≥3 Da over few centers → DIA; ≤2 Da, many centers → DDA | none (mzML rarely carries model reliably) |
| Thermo `.raw` | needs ThermoRawFileParser → mzML, else `unknown` | ThermoRawFileParser `metadata` model string, if present |
| `.wiff` | convert to mzML first | none |

Confidence is `high`/`medium`/`low`. **Anything not `high`, mixed, or with
disagreeing instruments sets `needs_confirmation` — ask the user.** Instrument may
be null; that's fine (matcher falls back to score 0 and the user confirms).

## Parameters

Search parameters are derived from the data type (instrument + acquisition) by
`estimate_params.py`, not hand-maintained per workflow. → `references/parameters.md`.

## Routing

Default by acquisition: **DIA → DIA-NN, DDA → Sage.** `--engine` overrides. The
bundle's `engine.name` is authoritative when present. FragPipe is opt-in only (the
bundle names it or the user asks) because MSFragger/IonQuant are license-gated.

## Licensing — pick a commercial-OK engine for non-academic users
| Engine | License | Commercial use | Notes |
|---|---|---|---|
| **DIA-NN** (Academia build) | academic / non-profit only | ❌ **no** | the free build the skill downloads is academic-only |
| **AlphaDIA** | Apache-2.0 | ✅ yes | open-source DIA alternative; deep-learning, **GPU recommended** |
| **Sage** | MIT | ✅ yes | fast, CPU; DDA + wide-window DIA |
| **FragPipe** (MSFragger/IonQuant) | license-gated | depends | requires the user's own license |

**For commercial / non-academic DIA users, route to AlphaDIA** (`--engine alphadia`)
instead of DIA-NN — same job, no license problem. Ask "academic or commercial?" before
a DIA run if it isn't already clear.

## Public program sources (anyone, any OS — cite these for non-Core / external users)

This skill is used **outside** UC Davis and the Core. Never assume the internal
`/quobyte/proteomics-grp` copy of a tool exists — **always point the user to the public
source** and a way to run it on *their* OS (macOS, Windows via WSL2 **or** native
Windows, or Linux). `setup.sh`/`acquire_tools.sh` fetch the free ones automatically.

| Program | Public source | Platforms | Install / note |
|---|---|---|---|
| **DIA-NN** (Academia) | https://github.com/vdemichev/DiaNN/releases | Windows, Linux | Windows GUI+CLI, or Linux binary; academic / non-profit only. `acquire_tools.sh` fetches the Linux build. DDA since 2.3 (`--dda`). |
| **Sage** | https://github.com/lazear/sage/releases | Win, macOS, Linux | MIT; single static binary. |
| **AlphaDIA** | https://github.com/MannLabs/alphadia | Win, macOS, Linux | Apache-2.0 (commercial-OK); `pip install alphadia`, GPU recommended. |
| **FragPipe** (MSFragger/IonQuant) | https://github.com/Nesvilab/FragPipe/releases | Win, macOS, Linux | Java GUI; MSFragger/IonQuant need the user's own (free-academic) license. |
| **ThermoRawFileParser** | https://github.com/compomics/ThermoRawFileParser/releases | Win, macOS, Linux | `.raw`→mzML; cross-platform .NET. |
| **ProteoWizard / msconvert** | https://proteowizard.sourceforge.io/ | Windows (native); Linux/macOS via Docker | vendor→mzML; Linux via the `chambm/pwiz-...` Docker image. |
| **.NET 8 runtime** | https://dotnet.microsoft.com/download/dotnet/8.0 · installer script https://dot.net/v1/dotnet-install.sh | Win, macOS, Linux | needed only so the **Linux** DIA-NN binary can read `.raw`; `ensure_dotnet8.sh` installs 8.0.latest. **Native Windows doesn't need it.** |
| **UniProt proteomes (FASTA)** | https://www.uniprot.org/proteomes/ | web / REST | e.g. `UP000005640` (human); `fetch_fasta.py` streams it. |
| **Contaminants (cRAP)** | https://www.thegpm.org/crap/ | web | `fetch_fasta.py --add-contaminants` appends a maintained set. |

**Rule:** whenever you tell a user to run a program, give the public link **+** the exact
command **+** an OS-appropriate path — not a HIVE/`/quobyte` path they can't reach.

## Per-engine invocation & output adapter (`run_search.py`)

### DIA-NN (native contract — no adapter)
```
<cmd> --cfg <bundle diann.cfg> --f <file> [--f <file> ...] \
      --fasta <fasta> --out report.parquet --threads N [--dda]
```
`report.parquet` is already the DE contract — `run_de.R` reads it directly.

**DDA data (DIA-NN 2.6+).** DIA-NN 2.3+ can search DDA / DDA-PASEF with **`--dda`**
(`run_search.py` adds it automatically when the bundle's acquisition is `DDA`). Notes:
- `--dda` **must** be used with DDA and **must not** be used with DIA data.
- QuantUMS is auto-disabled on DDA; PTM-localisation probabilities are unreliable on DDA.
- For DDA **quant**, DIA-NN recommends extra MS1 filtering on `Ms1.Global.Q.Value`
  (< 0.0001–0.01) and `Ms1.Global.Quality` (> 0.5–0.9), optionally `Ms1.Q.Value` /
  `Averagine`. Standard DIA/DDA q-value filters still apply.
- It's officially "beta", but performs strongly — on UC Davis nail (Exploris, 67 runs)
  it matched/beat the delivered FragPipe/Scaffold result at ~equal keratin coverage.
- Default routing still sends DDA → Sage (validated); use `--engine diann` (or a
  DDA+DIA-NN workflow) to search DDA with DIA-NN.

**Reading Thermo `.raw` directly (the .NET 8 gotcha — Linux only).** On **native
Windows** DIA-NN reads `.raw` out of the box (skip this whole section). On **Linux**
(incl. WSL2 and HIVE) the 2.6 *native* binary reads `.raw` via a bundled .NET
(ThermoFisher.CommonCore) component and needs a **.NET 8 runtime ≥ 8.0.17** on PATH. Without it, DIA-NN prints
`ERROR: cannot read .raw files ... .NET Runtime 8: 8.0.17 or later` and processes **0
files**. It does a HARD ≥ 8.0.17 check: an older 8.0.x (e.g. a cluster's `8.0.4`
module) is **rejected**, and a .NET 9 runtime is **not** used (rollForward is
LatestMinor — it won't cross the major version). Fix, done for you by
**`scripts/ensure_dotnet8.sh`** (idempotent; run it on a login node — it needs
internet, compute nodes usually don't):
```
export DOTNET_ROOT="$(bash scripts/ensure_dotnet8.sh | tail -1)"   # installs 8.0.latest if missing
export PATH="$DOTNET_ROOT:$PATH"
```
`run_search.py` calls this automatically whenever DIA-NN inputs include `.raw` (inline
run, and baked into the emitted sbatch as a preamble). When it's right, DIA-NN logs
`.NET runtime found, Thermo .raw support enabled`. Alternative with **no** .NET: feed
`.mzML` (DIA-NN reads those natively) — convert `.raw`→mzML with ThermoRawFileParser /
msconvert if you don't already have them.

### AlphaDIA (commercial-OK DIA; adapter required)
Apache-2.0 — the open-source alternative to DIA-NN for non-academic users. Library-free:
```
<cmd> -o <out> -f <raw> [-f ...] --fasta <fasta> [-c <config.yaml>]
```
- Install: `pip install alphadia` (into the conda env); `acquire_tools.sh` does this
  when AlphaDIA is pinned/requested. **GPU strongly recommended** (CPU works but is slow)
  — best on a HIVE GPU node. Reads `.raw`/`.d`/mzML; no msconvert needed.
- **Adapter:** AlphaDIA writes `pg.matrix.parquet` (protein-group × run matrix; also
  `precursors.parquet` with `raw.name`/`pg.name`/`pg.intensity`). `adapt_alphadia` melts
  the matrix → DIA-NN-shaped `report.parquet` (Run, Protein.Group, PG.MaxLFQ; Q-values
  zeroed since AlphaDIA already FDR-filtered). Like the Sage adapter, confirm it on real
  data the first time.
- **⚠ Limitation — timsTOF whole-proteome directDIA:** AlphaDIA does **not** work on
  **timsTOF** `.d` for **whole-proteome library-free (directDIA)** searches. For that
  case it is not a valid substitute for DIA-NN. Options: academic → DIA-NN; commercial →
  run AlphaDIA **library-based** (a predicted/empirical spectral library, `--library`)
  instead of directDIA, or use Sage. Non-timsTOF instruments and library-based timsTOF
  runs are unaffected.

### Sage (mzML-first; adapter required)
1. Convert `.d`/`.raw` → mzML with `msconvert` if needed (fails loudly if msconvert
   is absent and inputs aren't mzML).
2. `<cmd> <bundle sage_config.json> -f <fasta> -o <out> --parquet
   --disable-telemetry-i-dont-want-to-improve-sage <mzml...>`
3. **Adapter:** map `lfq.parquet` (protein, filename, intensity) →
   DIA-NN-shaped `report.parquet` with `Run, Protein.Group, PG.MaxLFQ` and zeroed
   Q-value columns (Sage already FDR-filtered). Q-values default to 0 so the
   MaxLFQ DE path keeps every row.
   - v0.14.x has no native protein grouping → may need `sage_protein_groups.py`
     post-hoc (DE-LIMP keeps it on HIVE). v0.15+ has IDPicker grouping.

### FragPipe (opt-in; adapter required)
- Build an `.fp-manifest` from the files + data type (DIA/DDA).
- `<cmd> --headless --workflow <.workflow> --manifest <m> --workdir <out>
  [--config-tools-folder $FRAGPIPE_TOOLS_FOLDER]`
- **Adapter:** `combined_protein.tsv` per-sample `MaxLFQ Intensity` columns →
  DIA-NN-shaped `report.parquet`.
- Needs Java 9+; MSFragger/IonQuant must already be licensed/present.
- **Reading Thermo `.raw`:** MSFragger reads `.raw` directly via
  `tools/ext/thermo/BatmassIoThermoServer.exe` — a **Windows .NET-Framework** exe that on
  **Linux needs `mono` on PATH** (the .NET 8 runtime does NOT run it — Framework ≠ Core).
  Native Windows/macOS: works out of the box. On a mono-less Linux cluster you get
  *"Could not find Batmass-IO Thermo binary"* + `Scans = 0` → either install mono
  (conda-forge `mono`, then `export PATH=<monoenv>/bin:$PATH` in the job) or pre-centroid
  to mzML. (On HIVE, mono 6.12 is installed at `/quobyte/proteomics-grp/conda_envs/mono`.)
- **Other DDA gotchas:** IonQuant MBR needs experiment groups in the manifest — set
  `ionquant.mbr=0` if you only want peptide/protein IDs; **profile-mode mzML** fail
  `CheckCentroid` (centroid during conversion, or read `.raw` so MSFragger vendor-centroids);
  the LFQ-MBR template ships **without** a `database.db-path` line (add one); and the
  headless run's exit code can read 0 even on a crash — **verify the output file exists**.
- **Which output file to verify (don't false-alarm):** the file depends on the run
  type. A **single-experiment / peptide-ID run** (`ionquant.mbr=0`, one experiment)
  writes **`psm.tsv` + `peptide.tsv` + `protein.tsv` + `ion.tsv`** — there is **no**
  `combined_peptide.tsv`/`combined_protein.tsv` (those only appear with **multi-experiment
  IonQuant LFQ**). A watcher that checks for `combined_*.tsv` will wrongly report
  "no output" on a perfectly good ID run. Success check = `peptide.tsv` (and `psm.tsv`)
  present and non-empty; use `combined_protein.tsv` only when the manifest defines ≥2
  experiment groups and you asked for the cross-run LFQ matrix.

## The adapters are the test surface

DIA-NN's path is fully wired and native. The **Sage and FragPipe → DE-contract
adapters are the part that needs real-data validation** (flagged in each
workflow's `VALIDATION.md`). Verify the protein×run matrix shape and that
intensities are in the expected (linear, pre-log) scale before trusting DE output.
