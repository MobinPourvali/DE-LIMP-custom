# Three-engine DIA comparison — Poplar 6×3, Orbitrap

**Date:** 2026-08-04 · **Dataset:** Drakakaki lab poplar, 18 Thermo Orbitrap runs,
30 min gradient, 3 conditions × 2 biological × 3 technical (`Ex030322_MinminPlopi_fDia_30m_*`).

All three engines were measured on the **same 18 raw files**, verified by basename.
(FRAN records these with a synthetic `.d` extension although they are Thermo `.raw` —
match on basename, never on extension.)

## Results

| | Spectronaut (FRAN, 2022) | **DIA-NN 2.6** (2026) | Radiant + Fulcrum (2026) |
|---|---|---|---|
| Protein groups, total | 3,545 | **4,693** | 3,201 |
| Precursors, total | n/a¹ | 29,718 | 21,370 |
| Median PG per run | ~3,470 | **4,149** | 2,578 |
| PG quantified in **all 18** runs | 3,180 | **3,288** | 1,675 |
| — as % of that engine's total | **89.7%** | 70.1% | 52.3% |
| Matrix completeness | **98.1%** | 87.5% | 80.2% |

¹ FRAN stores Spectronaut rows without `precursor_id_diann`, so a distinct-precursor
count is not comparable on that side. Per-run precursor rows were 22,268–24,891.

## Reading this table

**DIA-NN wins on depth**: +32% protein groups over Spectronaut and +47% over Radiant,
and the highest median per run.

**Do not read the completeness percentages as a quality ranking.** That column is
denominator-sensitive — a deeper search necessarily reaches further into low-abundance
proteins that are not seen in every run, which lowers the percentage while *increasing*
the absolute number of fully-observed proteins. In absolute terms DIA-NN also has the
most proteins present in all 18 runs (3,288 vs Spectronaut's 3,180). Spectronaut's 98.1%
reflects a shallower, more conservative list, not a better matrix.

**Radiant is last on every measure**, and its sparsity is real rather than an artefact of
depth: it is both the shallowest (3,201 PG) *and* the least complete (52.3% in all runs).

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
