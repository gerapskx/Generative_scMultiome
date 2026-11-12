# AGENTS.md

R + Python workspace for the GSE193688 substantia nigra multiome project: QC →
annotation → scDesign3 simulation → ML fidelity → counterfactual reversal.

Two layers coexist:

- **Legacy scripts** at the root (`*.R`) — interactive RStudio drafts, run
  top-to-bottom, sharing state through the global environment.
- **Pipeline** (`R/`, `src/`, `workflow/`, `config/`) — the Snakemake build.
  This is what to extend.

## Layout

| Path | Role |
|---|---|
| `seuratR.R` | Legacy. Builds the sample sheet (`metadata`) and the ATAC recount path. |
| `06112026_Multiome.R` | Legacy. Current multiome + scDesign3 run. **Depends on `metadata` and `features` from `seuratR.R` in the same session.** |
| `scDesign3.R` | Legacy. GSE298008 download + caret/pROC/ranger classifier prototype (5-fold CV, ROC). |
| `scMultiome_Design.R` | Legacy, abandoned draft. Do not extend. |
| `R/00_ingest.R` | Stage 0. QC, doublets, consensus peaks, exact ATAC recount. |
| `src/assemble.py` | Stage 1. MTX → `data/processed/{rna,atac}.h5ad`. |
| `src/make_sample_sheet.py` | Rebuilds `config/samples.csv` from disk; preserves hand-edited `use`. |
| `config/config.yaml` | Analysis params. |
| `config/samples.csv` | Single source of truth for sample selection (`use` column). |
| `Snakefile` | Orchestration (repo root, so `conda:`/`script:` paths resolve). |

## Commands

```bash
conda activate multiome-py
python src/make_sample_sheet.py       # add --all to enable all 34 samples
snakemake -n -j 2                     # dry run
snakemake -j 8 --use-conda            # stage 0 (R) + stage 1 (python)
```

`snakemake` resolves `conda:`/`script:` paths relative to the Snakefile, which is
why it sits at the repo root rather than in `workflow/`.

## Environment

- **R is not installed on the Linux side** (no `R`/`Rscript`, no conda R env).
  The legacy scripts only ever ran in Windows RStudio with this folder as the
  working directory. `envs/r_env.yaml` is provided but unbuilt — verify before
  promising an R run.
- Scripts rely on the CWD containing the data (`list.files(pattern="\\.h5$")`,
  relative `GSE193688_RAW/`), so the working directory is load-bearing.
- `library(GOstats)` in both legacy scripts is never used and is deprecated in
  recent Bioconductor — the most likely cause of a `library()` failure. Delete it.
- Legacy code sets `plan(multisession, workers=22)` and also asks scDesign3 for
  `n_cores=12` and BiocParallel — that oversubscribes. Pick one backend.

## Data

- `GSE193688_RAW/` holds **all 34 samples** (28 GB): h5, fragments, peaks BED,
  prebuilt `.tbi` indexes. Do not re-index.
- `/mnt/c/Users/gerar/Downloads/GSE193688_RAW/` has only **2** of the 34. Not the
  data root.
- Group comes only from the donor prefix (`seuratR.R:34-44`): `^Y`→Young,
  `^C`→Aged, `^P`→PD. `YC4345` therefore counts as **Young**. No metadata file
  exists; this rule is the only source of truth. Cohort: 9 Aged, 15 PD, 10 Young.
- `GSM######_<donor>_filtered_feature_bc_matrix.h5`; donor = strip `^GSM[0-9]+_`.
- `GSE298008/` is a **failed download** — its `.mtx.gz` files are HTML error
  pages. The real data is in `GSE298008_11022026/`.
- `scDesignSCATAC.RData` is 2.5 GB; do not `load()` it casually.

## Signac / Seurat gotchas

- After `GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)` always set
  `seqlevelsStyle(annotations) <- "UCSC"` and `genome(annotations) <- "hg38"`, or
  the 10x `chr1:100-200` peak names never match.
- Peak-name separator differs by source: h5 → `StringToGRanges(sep = c(":", "-"))`,
  `FeatureMatrix()` → `StringToGRanges(sep = c("-", "-"))`. Mixing silently
  yields zero ranges.
- `merge()` on ChromatinAssays reduces peaks and leaves **structural zeros** for
  sample-specific peaks (125,610 consensus vs 54,342 one-donor-only for the 2
  pilot donors). `R/00_ingest.R` avoids it with `reduce()` + `FeatureMatrix()`;
  the legacy script still has the problem.
- `reduce()` default `min.gapwidth = 1` merges only touching intervals (gap 0 in
  0-based BED coordinates), not 1 bp gaps. Any reimplementation must use gap 0.
- Seurat v5: `JoinLayers()` before `GetAssayData(layer="counts")` on a merged object.
- `CreateFragmentObject(cells=)` needs a **named** vector mapping raw 10x
  barcodes → prefixed cell names.

## scDesign3 gotchas

- RNA and ATAC are fit **separately**; no joint copula. Cells are paired only
  through shared covariates.
- Returns `new_count` / `new_covariate` with **no dimnames** — assign names yourself.
- `select_variable_features()` is not a scDesign3 function; it is a local helper
  at `06112026_Multiome.R:315` using raw-count variance.
- ATAC calls need `parallelization="bpmapply"`; the RNA call does not.

## Known-broken lines in legacy scripts (do not copy)

- `06112026_Multiome.R:118` `sample(seq_along(sample.list), 5)` — errors with <5 samples.
- `:174-176` `group.by="celltype"` before `celltype` exists (`:225`).
- `:184` `features` undefined (lives in `seuratR.R:270`).
- `:213` `new_names` never applied; `"6"="OPcs"` typo.
- `:677` `MergedSample` → `MergedSamples`.
- `seuratR.R:300` `DimPlot(..., lab` truncated.

## Conventions

- Cell types are assigned manually from clusters (`RenameIdents`); `UN` = unassigned.
- Marker panel must be midbrain-specific: DA neurons `TH/SLC6A3/SLC18A2/NR4A2/SOX6`,
  GABAergic `GAD1/GAD2`. `PDGFRB` is pericyte, not endothelial; `TBR1` is cortical.
- `repo/` is an empty git repo (no commits). The Flask app at
  `/home/gerardoiga/2026/MPI-IE/Simon/Repo` is a mock UI with hardcoded fake runs —
  not the pipeline.
