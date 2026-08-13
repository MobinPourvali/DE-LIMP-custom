# Three-engine DIA comparison — Poplar 6×3, Orbitrap

**Date:** 2026-08-04 · **Dataset:** Drakakaki lab poplar, 18 Thermo Orbitrap runs,
30 min gradient, 3 conditions × 2 biological × 3 technical (`Ex030322_MinminPlopi_fDia_30m_*`).

All three engines were measured on the **same 18 raw files**, verified by basename.
(FRAN records these with a synthetic `.d` extension although they are Thermo `.raw` —
match on basename, never on extension.)

## Results

| | Spectronaut (FRAN, 2022) | **DIA-NN 2.6** (2026) | Radiant + Fulcrum, **MBR on** | Radiant, MBR off (superseded) |
|---|---|---|---|---|
| Protein groups, total | 3,545 | **4,693** | 3,178 | 3,201 |
| Precursors, total | n/a¹ | 29,718 | 21,343 | 21,370 |
| Median PG per run | ~3,470 | **4,149** | 3,106 | 2,578 |
| PG quantified in **all 18** runs | 3,180 | **3,288** | 2,777 | 1,675 |
| — as % of that engine's total | 89.7% | 70.1% | **87.4%** | 52.3% |
| Matrix completeness | **98.1%** | 87.5% | 97.5% | 80.2% |

¹ FRAN stores Spectronaut rows without `precursor_id_diann`, so a distinct-precursor
count is not comparable on that side. Per-run precursor rows were 22,268–24,891.

## What MBR changed, and what it did not

Turning MBR on (SLURM 20105472, same 18 per-file searches reused for pass 1, so **only**
MBR differs):

| | MBR off | MBR on | change |
|---|---|---|---|
| Total PG | 3,201 | 3,178 | −0.7% (flat) |
| PG in all 18 runs | 1,675 | 2,777 | **+66%** |
| Matrix completeness | 80.2% | 97.5% | **+17.3 pts** |
| Median PG per run | 2,578 | 3,106 | +20% |

This is exactly the signature of match-between-runs: it does **not** find new proteins,
it finds already-identified ones in more runs. Radiant's data matrix is now as complete
as Spectronaut's (97.5% vs 98.1%) and more complete than DIA-NN's.

**The sparsity criticism in the first version of this document was an artefact of our
configuration, not a property of Radiant.** It has been withdrawn.

**The depth gap is real and survives.** DIA-NN finds 4,693 protein groups to Radiant's
3,178 (+48%) and 4,149 per run to Radiant's 3,106 (+34%). That is the substantive
finding, and MBR does not touch it.

On the percentage column: it is denominator-sensitive. DIA-NN's 70.1% is the lowest
precisely because its denominator is the largest; in absolute terms DIA-NN still has the
most proteins present in all 18 runs (3,288, vs Spectronaut 3,180 and Radiant 2,777).
Do not read that column as a quality ranking.

## Known defect in the DIA-NN column (not corrected — deliberate)

This DIA-NN run used an **inconsistent scan window**. DIA-NN optimises the radius per
file when `--window` is absent, and the cfg's `--window 0` was rejected outright
(`WARNING: scan window radius should be a positive integer`), so it fell back to auto:
**seventeen files got radius 7, one got 8**. Steps 3/5 then reused `.quant` files across
that split, which DIA-NN itself warns against:

> WARNING: combining reuse of .quant files with automatic optimisation of mass
> accuracies or scan window will lead to results that are different from those of the
> original analysis that produced the .quant files and is strongly not recommended

**Decision (2026-08-13): not re-run.** One file of eighteen differed by one radius unit,
and the DIA-NN-vs-Radiant gap is +48% — far too wide for that to flip the ranking. So the
*ordering* in this table stands.

What does NOT stand is the exact figure: **4,693 is not a reproducible number**, and it
should not be quoted as a DIA-NN benchmark, put in a methods section, or compared against
a future DIA-NN run. The chain now pins the window itself (step 1b), so any new run is
reproducible and would supersede this column.

## Confounds — do not quote these numbers as a clean engine benchmark

1. **The libraries are not the same.** Spectronaut used its own project library; DIA-NN
   and Radiant used a DIA-NN *predicted* library built from the poplar FASTA. Radiant's
   gap in particular is partly a library effect, not purely engine performance.
2. **2022 vs 2026 software generations.** The original Spectronaut version and settings
   are not recorded in FRAN.
3. FRAN's Spectronaut rows are already 1%-filtered on ingest (`q ≤ 0.01` count equals the
   total row count on every file). DIA-NN and Radiant were filtered at 1% here
   (`Q.Value`, `Lib.Q.Value`, `Lib.PG.Q.Value`). Comparable, but not identically derived.

A clean engine comparison would need all three run against the *same* library.

## DE-LIMP readiness

All three produce precursor-level output that limpa/DPC can consume:

| Engine | File | Route |
|---|---|---|
| DIA-NN | `report.parquet` | native `readDIANN()` |
| FragPipe | `report.tsv` | native |
| Radiant | `delimp_report.parquet` | via `scripts/radiant_to_delimp.py` |

limpa support is a property of the **file**, not the engine. The one file that does *not*
work is an adapted `report.parquet` that has lost `Precursor.Id` / `Precursor.Normalised`.

## Provenance

- DIA-NN: `/quobyte/proteomics-grp/brett/poplar_test/diann/out_q/report.parquet`, SLURM 20000423 (5-step parallel chain), finalpass ~21 min/file × 18.
- Radiant: `/quobyte/proteomics-grp/brett/poplar_test/radiant/out/delimp_report.parquet`.
- Spectronaut: FRAN `delimp_searches` id `62710070-705b-5edd-9742-f2ce342b3a4e`
  (`20220308_120552_Minmin_Poplar_March_2022_3conditions`), read-only query.
