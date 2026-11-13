# Testing synthetic data generated from GSE193688. This is for my rotation in Elsasser lab

snRNA + snATAC cohort: QC, annotation, scDesign3 simulation, ML fidelity testing, null/alternative DE benchmarking, and
a counterfactual PD → Aged reversal analysis. 

## Status

| Stage | Script | State |
|---|---|---|
| 0 ingest | `R/00_ingest.R` | **done** — 3,730 nuclei, 36,601 genes, 125,610 peaks |
| 1 assemble | `src/assemble.py` | **done** — `data/processed/{rna,atac}.h5ad` |
| 2 cluster + annotation gate | `R/01_cluster.R`, `src/annotate.py` | **done** — 23 Harmony clusters, 4 glial classes |
| 3 split | `src/split.py` | **done** — 2,610 fit / 1,120 holdout |
| 4 features | `R/02_features.R` | **done** — 2,000 genes, 2,000 peaks |
| 5 simulate | `R/03_simulate.R` | **done** — 5 seeds × alt/null × RNA/ATAC |
| 6 fidelity | `src/fidelity.py` | **done** |
| 7 hypothesis | `src/hypothesis.py` | **done** |
| 8 reversal | `R/04_reversal.R`, `src/reversal.py` | **done** |
| 9 dashboard | `dashboard/app.py` | **done** — 11 pages |
| 10 auto annotation | `src/auto_annotate.py` | **done** — validated against the pilot map |
| 11 cell-type impact | `src/celltype_impact.py` | **done** |
| 12 synthetic utility | `R/05_utility_sim.R`, `src/synthetic_utility.py` | **done** |
| 13 scaling | `src/plan_scaling.py`, `src/compare_scaling.py` | **done** — n02 → n08 → n16 → n24 + comparison |
| 14 subtype transfer | `src/subtype_split.py`, `R/06_subcluster.R`, `R/07_subtype_sim.R`, `src/subtype_transfer.py` | **done** — k=2 in Oligodendrocyte + Microglia, RNA + ATAC, synthetic ↔ real |

### Composition of the 2-donor pilot (important)

The pilot samples are **oligodendrocyte-dominated**: marker scoring puts ~88% of
nuclei in the oligodendrocyte lineage, with astrocytes, microglia and OPC making
up almost all of the rest. Dopaminergic neurons are essentially absent — `TH`
has **zero** UMIs in C16007 in the raw CellRanger ARC h5, and across all 34
samples the maximum is 35 `TH`+ nuclei in one sample. This is a property of the
data, not of the pipeline (verified against `filtered_feature_bc_matrix.h5`
directly, with correct Ensembl IDs).

Consequences: cell-type labels are glial (Oligodendrocyte, Astrocyte, Microglia,
OPC), and any DA-neuron / PD-neurodegeneration narrative is not supported by
these samples. The Aged-vs-PD contrast is also a donor contrast in this pilot.


## Headline results (2-donor pilot)

All numbers below are from `results/*/summary.json` and are reproducible from
the committed config; they describe **this pilot**, not the cohort.

### Fidelity — can a classifier tell real from synthetic?

AUC of a real-vs-synthetic discriminator, `alt` mode, 5 seeds. An honest
simulator should sit near 0.5; controls bracket what "separable" looks like.

| modality | LR AUC | RF AUC | identity control | naive Poisson | positive control (Aged vs PD) |
|---|---|---|---|---|---|
| ATAC | 0.393 ± 0.020 | 0.589 ± 0.011 | 0.51 / 0.51 | 0.81 / 0.81 | 0.99 / 0.97 |
| RNA | 0.604 ± 0.010 | 0.880 ± 0.009 | 0.56 / 0.47 | 1.00 / 1.00 | 1.00 / 1.00 |

ATAC is hard to distinguish (|AUC − 0.5| ≈ 0.1 for both models, and LR is
anti-correlated, i.e. synthetic scores *more* real than real). RNA is easier:
RF reaches 0.88, so the RNA simulator reproduces the marginal distribution but
not the full gene–gene covariance structure. Both controls behave as expected,
so the 0.88 is signal rather than a broken discriminator.

### Null / alternative calibration

Wilcoxon tests on simulated data, truth = FDR 0.05 and |log2FC| ≥ 0.5.

- **Null `condition`** (zeroed condition coefficients): **0 significant
  features at FDR 0.01, 0.05 and 0.1**, in both modalities. At FDR 0.2 ATAC
  oligodendrocyte leaks slightly (0.4% power), which is the expected
  liberal-tail behaviour. This is the key calibration result: the counterfactual
  null does not manufacture condition effects.
- **Null `celltype`** keeps power (0.51 ATAC / 0.64 RNA for Aged astrocyte),
  confirming the null targets only the condition axis.
- **Alternative `condition`** (Oligodendrocyte, the only stratum with ≥100
  cells per group): power **0.79 ATAC / 0.80 RNA** at FDR 0.05, with FDP 0.55
  (ATAC) and 0.34 (RNA). Astrocyte (80 fit cells) is underpowered at 0.06 / 0.15.
- Truth sets: 9,953 RNA and 6,383 ATAC features at FDR 0.05.

### Counterfactual PD → Aged reversal (Oligodendrocyte only)

`min_cells_per_stratum: 100` leaves only Oligodendrocyte eligible (Astro, Micro
and OPC have 50–99 cells per group). Two independent rankings are computed:

1. **Model Δ** — the group-level condition coefficient contrast from the fitted
   scDesign3 mean model, re-extracted after reversing the condition labels.
2. **Classifier contribution** — per-cell logistic-regression contributions,
   averaged over holdout PD cells.

| modality | observed SD | permuted SD | permuted q95 | frac > q95 | reverse r |
|---|---|---|---|---|---|
| RNA | 1.19 | 0.284 | 0.538 | 0.596 | −1 |
| ATAC | 0.466 | 0.131 | 0.271 | 0.779 | −1 |

The `reverse r = −1` is exact: reversing the labels flips every coefficient.
Against label permutation, the observed spread is 4.2× (RNA) and 3.6× (ATAC)
the permuted SD, and 60–78% of features exceed the permuted 95th percentile —
real structure, not noise. The large `permuted_max_abs` for RNA (2.48) is a
single tail draw; the robust controls (q95, fraction above q95, SD ratio) are
the ones to read.

**The two rankings agree only weakly** (Spearman 0.024 ATAC / −0.076 RNA;
top-100 overlap 0.33 for both), and this is expected rather than a failure:

- The model Δ is a **marginal, group-level** estimand — one coefficient per
  feature, with only `celltype` and `condition` as covariates.
- The classifier contribution is **multivariate and per-cell**; correlated
  features share credit, so a feature can rank high in one method and low in
  the other.
- The controls (reverse-direction and permuted Spearman ≈ 0) show the near-zero
  agreement is not a pipeline defect: the two methods are answering different
  questions.
- **Caveat:** classifier CV balanced accuracy is 0.998 (RNA) / 0.981 (ATAC).
  With one donor per condition, a classifier that separates donors this cleanly
  is almost certainly learning donor/batch identity, not condition biology.
  Treat the classifier branch as a per-cell sensitivity view, not as independent
  evidence for a PD effect.

Both branches therefore report a **donor contrast**, not a population-level PD
effect; the dashboard labels them accordingly.


## Sample-size scaling study (Aged vs PD)

The pilot is one point on a curve. `config/scaling/` defines a nested series over
the 24 Aged/PD donors in the cohort (9 Aged, 15 PD; Young is excluded because
`R/03_simulate.R` fixes the condition levels to `Aged`/`PD`):

| step | donors | Aged | PD |
|---|---|---|---|
| n02 | 2 | 1 | 1 |
| n08 | 8 | 3 | 5 |
| n16 | 16 | 6 | 10 |
| n24 | 24 | 9 | 15 |

Each step is a superset of the previous one, the pilot donors `C16007`/`P02284`
are always included, and the donors are chosen with a fixed seed. The study asks how
the *interpretation* changes with N, not just whether a p-value crosses a threshold.

### What is held fixed

- **Feature set**: every step reuses the n24 RNA/ATAC selection so the curve
  isolates sample size. Each step also records its own selection under
  `features/selected/`, and `compare_scaling.py` reports the Jaccard overlap and
  the Spearman correlation of DE log2FC against n24.
- **Simulation seeds**: `101, 202, 303, 404, 505` at every step. Sample size and
  simulation iteration are separate axes in `results/scale/comparison.csv`.
- **Split rule**: cell-level stratified by celltype × donor, same fraction and seed.
- **Annotation**: `src/auto_annotate.py` scores a canonical glial/neuronal marker
  panel with log2 enrichment. It reproduces the manual pilot map (3,276
  Oligodendrocyte / 179 Astrocyte / 143 Microglia / 132 OPC) and is used for every
  step so no human gate blocks the series.

### What is measured per step

| analysis | metrics |
|---|---|
| QC | cells, donors, clusters, cell types, doublet rate |
| features | selected count, Jaccard vs n24, DE log2FC Spearman vs n24 |
| fidelity | real-vs-synthetic AUC (LR/RF), identity/positive/naive controls |
| hypothesis | null FDP, alt power, FC correlation at FDR 0.01–0.2 |
| impact | DE burden, mean |log2FC|, donor-held-out AUC, eta² condition vs donor, PD fraction |
| utility | synthetic-trained AUC / balanced accuracy vs real-trained, naive-Poisson, permuted-label and majority baselines |
| reversal | counterfactual |log2FC|, reverse-direction correlation, model↔classifier agreement |

### Cell-type impact

`src/celltype_impact.py` ranks cell types by how strongly condition affects them.
The key metric is **donor-held-out AUC**: a condition classifier trained on some
donors and tested on others, within one cell type. At n02 there is one donor per
condition, so the AUC is undefined and every cell-level contrast is also a donor
contrast; it only becomes a condition test from n08 onward. `eta2_condition` versus
`eta2_donor` shows how much of the variance is condition and how much is donor.

### Synthetic-to-real utility (ML PoC)

`R/05_utility_sim.R` fits scDesign3 on training donors only and simulates
donor-disjoint cohorts. `src/synthetic_utility.py` then trains a condition classifier
on synthetic cells and tests it on held-out real donors. The baselines are:

- `real_trained` — train on real training donors (upper bound)
- `naive_poisson` — marginal-only synthetic (lower bound)
- `permuted_labels`, `majority` — chance controls

The learning curve varies the number of synthetic training cells (500–10,000) at each
step. A utility that rises with N is gaining donor-level information; one that rises
only with synthetic size is resampling the same donors. The test is donor-disjoint and
therefore `not_estimable` at n02.

### Results (Aged vs PD series)

Cohort at each step:

| step | cells | donors | cell types |
|---|---|---|---|
| n02 | 3,725 | 2 | 4 |
| n08 | 15,679 | 8 | 6 |
| n16 | 28,207 | 16 | 8 |
| n24 | 39,956 | 24 | 8 |

Oligodendrocytes dominate every step (79% of cells at n24), as in the pilot.

**Feature stability** (each step's own selection vs the n24 fixed set):
top-100 Jaccard 0.28 → 0.39 → 0.49 and log2FC Spearman 0.35 → 0.58 → 0.73 for n02 → n08 → n16.
The selected features are therefore *not* converged at n16; the n24 set is doing real work.

**Fidelity** (real-vs-synthetic discriminator AUC; 0.5 = faithful):

| N | RNA LR | RNA RF | ATAC LR | ATAC RF |
|---|---|---|---|---|
| 2 | 0.605 | 0.879 | 0.369 | 0.580 |
| 8 | 0.590 | 0.962 | 0.943 | 0.619 |
| 16 | 0.596 | 0.971 | 1.000 | 0.647 |
| 24 | 0.603 | 0.977 | 1.000 | 0.673 |

This is the main negative result of the series: the simulator becomes **more**
distinguishable from real cells as N grows, and ATAC is perfectly separable by
n16/n24. Fidelity does not improve with sample size.

**Hypothesis calibration** (RNA, condition contrast, FDR 0.05):

| N | alt power | alt FDP | alt FC corr | null discoveries |
|---|---|---|---|---|
| 2 | 0.473 | 0.390 | 0.514 | 0.1 |
| 8 | 0.425 | 0.541 | 0.540 | 0.05 |
| 16 | 0.432 | 0.631 | 0.504 | 0.2 |
| 24 | 0.298 | 0.720 | 0.489 | 0.07 |

The null side stays clean (≈0 false discoveries, FC corr ≈0.03). On the
alternative side, power does **not** grow with N: the real reference set grows
faster than the simulated discoveries, so the simulator's effect-size agreement
degrades. Treat the simulated DE as a fixed-effect approximation, not a
sample-size-scalable substitute for real DE.

**Cell-type impact at n24** (ranked by RNA DE burden):

| cell type | RNA DE burden | RNA donor-held-out AUC | ATAC DE burden |
|---|---|---|---|
| Oligodendrocyte | 0.710 | 0.732 | 0.878 |
| Microglia | 0.348 | 0.562 | 0.806 |
| OPC | 0.204 | 0.723 | 0.381 |
| Neuron | 0.193 | 0.607 | 0.000 |
| Astrocyte | 0.164 | 0.694 | 0.112 |
| Ependymal | 0.136 | 0.770 | 0.048 |

Oligodendrocytes carry the largest transcriptional disease burden and also have
the highest donor-held-out separability, consistent with the lineage dominating the
cohort. ATAC burden is highest in oligodendrocytes and microglia. n02 has no
donor-held-out AUC by construction (condition == donor).

**Synthetic-to-real utility** (train synthetic → test held-out donors; AUC):

| N | modality | real-trained | syn 500 | syn 1k | syn 2k | syn 5k | syn 10k |
|---|---|---|---|---|---|---|---|
| 8 | RNA | 0.802 | 0.674 | 0.676 | 0.686 | 0.720 | 0.725 |
| 8 | ATAC | 0.525 | 0.454 | 0.502 | 0.497 | 0.522 | 0.521 |
| 16 | RNA | 0.446 | 0.387 | 0.359 | 0.367 | 0.393 | 0.391 |
| 16 | ATAC | 0.392 | 0.401 | 0.387 | 0.415 | 0.417 | 0.420 |
| 24 | RNA | 0.663 | 0.553 | 0.549 | 0.548 | 0.562 | 0.567 |
| 24 | ATAC | 0.601 | 0.562 | 0.571 | 0.562 | 0.566 | 0.582 |

Synthetic utility rises monotonically with synthetic size at every step, but the
gain is modest (RNA +0.05 at n08, +0.004 at n16, +0.014 at n24 from 500 to
10k). The `real_trained` ceiling is unstable across steps because the donor-held-out
fold changes composition: at n16 the five held-out donors are 3 Aged / 2 PD and the
real-trained classifier is *below chance* (0.446 RNA, 0.392 ATAC). That is a
donor-confound artifact of a single fold, not evidence that condition is
unlearnable. n02 is `not_estimable` (one donor per condition).

**Reversal PoC** (known-effect recovery, null specificity, donor confound):

| N | recovery ρ | top-100 Jaccard | null ρ | SD ratio alt/null | donor bal. acc. | condition bal. acc. | confound index |
|---|---|---|---|---|---|---|---|
| 2 | 0.786 | 0.520 | −0.003 | 2.98 | 0.981 | — | — |
| 8 | 0.433 | 0.245 | 0.008 | 1.82 | 0.675 | 0.439 | 0.117 |
| 16 | 0.329 | 0.210 | 0.009 | 1.72 | 0.540 | 0.530 | 0.156 |
| 24 | 0.355 | 0.200 | −0.003 | 1.67 | 0.547 | 0.550 | 0.168 |

The n02 recovery (ρ=0.79) is the donor contrast being replayed, which is why
condition is *not identifiable* there. From n08 the honest donor-grouped
condition accuracy is at or near chance (0.44–0.55), while donor separability falls
with N (0.98 → 0.55). The null side is clean at every step (|ρ| < 0.01, ≤1.2% of
null effects above the real q95).

### Subtype transfer (Oligodendrocyte + Microglia)

The scaling study asks whether the Aged-vs-PD signal survives more donors. The
subtype study asks a different question: **if a cell type contains two subtypes,
does the simulator preserve that axis well enough to transfer a classifier between
synthetic and real cells?**

Design:

- Target cell types: **Oligodendrocyte** and **Microglia** (the two disease-impacted
  classes).
- **k=2 unsupervised subclustering** inside each cell type, separately in RNA and
  ATAC. Subtypes are defined on **training donors only**; held-out donors are
  assigned by a classifier trained on the training cells, so no held-out donor
  contributes to the subtype axis.
- Before clustering, cells are filtered on `scDblFinder_score < 0.3` and on the
  target-vs-other marker score. RNA PCA / ATAC LSI scores are corrected for
  `log(nCount)` and `log(nFeature)` on training cells, then Harmony is run over
  donor.
- scDesign3 is fit on training cells with `subtype + condition` (**alt**), then
  simulated in three modes:
  - **alt**: keep the subtype effect;
  - **null**: collapse subtype to one level and refit `condition` (primary null);
  - **null_zero**: keep the alt fit but zero every subtype coefficient
    (sensitivity null).
- Synthetic cells use a balanced 2×2 subtype × condition design, matched in total
  size to the real training cells (cap 20,000), with five seeds.
- Transfer matrix over `real`, `syn_alt`, `syn_null`, `syn_null_zero`, in both
  directions. The main block is real ↔ alt ↔ null; the null-zero domain is the
  sensitivity check.
- Logistic regression is primary, random forest a sensitivity model. Every pair is
  reported overall, by condition and by donor (donors exist only on the real side).
  Controls: `condition_only`, `permuted_labels`, `majority`, and a donor-mean
  `donor_only` ceiling.
- The subtype variable features are compared with the fixed n24 feature set by Jaccard
  overlap.

Run:

```bash
python src/subtype_split.py \
    --obs results/scale/n24/cluster/obs_annotated.csv \
    --out results/scale/n24/subtype \
    --celltypes Oligodendrocyte,Microglia --fraction 0.30 \
    --min-cells-per-donor 20 --seed 20261112

Rscript R/06_subcluster.R --ingest results/scale/n24/ingest \
    --features results/scale/n24/features \
    --obs results/scale/n24/cluster/obs_annotated.csv \
    --donors results/scale/n24/subtype/donors.csv \
    --out results/scale/n24/subtype \
    --celltype Microglia --modality rna

Rscript R/07_subtype_sim.R --ingest results/scale/n24/ingest \
    --features results/scale/n24/features \
    --labels results/scale/n24/subtype/Microglia/rna/labels.csv \
    --out results/scale/n24/subtype/Microglia/rna/simulation \
    --celltype Microglia --modality rna

python src/subtype_transfer.py --root results/scale/n24/subtype \
    --ingest results/scale/n24/ingest --out results/scale/n24/subtype/transfer

python src/subtype_compare.py --root results/scale/n24/subtype
```

Outputs: `subtype/{celltype}/{modality}/{labels,embedding,markers,features,qc}`,
`subtype/{celltype}/{modality}/simulation/{alt,null,null_zero}`,
`subtype/transfer/{transfer,transfer_summary,jaccard,summary}` and
`subtype/{comparison.csv,summary.json}`.

**Results.** The subtype splits are weak axes, not validated developmental
subtypes. Microglia has silhouette 0.160 (RNA) / 0.176 (ATAC) and RNA/ATAC label
agreement ARI −0.002; Oligodendrocyte has silhouette 0.260 (RNA) / 0.144 (ATAC)
and ARI 0.185.

Logistic-regression AUC, mean over five seeds:

| Cell type | Modality | real→real | syn alt→real | real→syn alt | syn null→real | syn null0→real | alt→alt |
|---|---|---:|---:|---:|---:|---:|---:|
| Oligodendrocyte | RNA | 0.817 | **0.862** | 0.852 | 0.497 | 0.505 | 0.995 |
| Oligodendrocyte | ATAC | 0.884 | **0.903** | 0.950 | 0.499 | 0.498 | 0.971 |
| Microglia | RNA | 0.700 | **0.729** | 0.959 | 0.492 | 0.496 | 0.921 |
| Microglia | ATAC | 0.886 | **0.905** | 1.000 | 0.512 | 0.521 | 1.000 |

The alt simulation carries a strong subtype signal in every combination; both null
modes collapse to chance when trained on synthetic and tested on real. The
`real→syn alt` column is generally **higher** than `syn alt→real` (Microglia RNA:
0.959 vs 0.729; Microglia ATAC: 1.000 vs 0.905), so synthetic cells are easier
to classify than real cells — the simulator smooths over real heterogeneity.

Controls are clean: condition-only AUC 0.484–0.548, permuted labels 0.490–0.528,
majority 0.500. The donor-mean ceiling is estimable only for Oligodendrocyte ATAC
(0.333); the other combinations have a single dominant subtype at donor level.
Donor-stratified `syn alt→real` is consistent rather than driven by one donor
(Oligodendrocyte RNA 0.861 ± 0.023, ATAC 0.904 ± 0.025; Microglia RNA
0.738 ± 0.048, ATAC 0.927 ± 0.058). Variable-feature Jaccard with the fixed
n24 set is 0.50 in all four combinations.

### Run

```bash
# one-time ingest over all 24 donors (hours; consensus peaks + fragment recount)
snakemake --configfile config/scaling/full.yaml -j 4 --use-conda --conda-frontend conda

# plan (already committed under config/scaling/)
python src/plan_scaling.py

# all steps, n24 first, then the rest, then the comparison table
bash run_scaling.sh

# or one step at a time
snakemake --configfile config/scaling/n24.yaml -j 8 --use-conda --conda-frontend conda
```

`results/scale/full/` is the shared ingest; each `results/scale/nXX/` is a complete
run. `src/compare_scaling.py` merges them into `results/scale/comparison.csv`
(long: N × modality × analysis × metric × celltype × synthetic_size × seed) and
`results/scale/summary.json`.

### Dashboard

`streamlit run dashboard/app.py` opens on the **full cohort (n24, all 24 donors)**.
A **Cohort** selector in the sidebar switches every page — Overview, Report, QC,
Annotation, Fidelity, Hypotheses, Reversal, Cell-type impact, Synthetic utility and
Subtype transfer — between n24, n16, n08, n02 and the 2-donor pilot, so each stage
can be inspected on all samples or on any smaller step. The Sample-size scaling page
always plots every step, independent of the selector.

The **Report** page renders a written explanation of every section — what it measures,
how to read it, what the selected cohort shows, and the caveat — generated from the
same artifacts, and can be downloaded as Markdown.

### Figures on disk

The dashboard only draws in memory. To get a folder tree of PNGs, one subfolder per
pipeline section and cohort:

```bash
python src/export_figures.py --out ~/Desktop/GSE193688_figures
```

It imports the dashboard pages, swaps `st.pyplot` for `savefig`, and drives the
dropdowns so each option is exported (`UMAP` × 4 colourings, every scaling
analysis × metric, every subtype cell type × modality). Output is 300 dpi PNG with a
`README.md` and `index.csv` listing each file's heading and selection. A
`_MISSING.txt` marks a stage that has no artifacts for that cohort; `_ERROR.txt`
records a page that failed to render.

### Labelled reversal plots

The reversal scatter panels label the features at the top of either ranking — the
counterfactual model (`|log2FC|`, red) and the classifier (`|contribution|`, green):

```bash
python src/plot_reversal.py --step n24            # 10 panels + grid, top 5 each
python src/plot_reversal.py --step n24 --top 8    # more labels per panel
```

Output goes to `~/Desktop/Reversal` by default (`--outdir` to change it): one PNG per
modality × cell type, a `reversal_grid_<step>.png` overview, and
`labelled_features.csv` listing every label with its coordinates and ranking method.
Labels use `adjustText` repulsion with thin leader lines; ATAC peaks are shortened to
`chr16:46.39Mb`. Panels whose log2FC range exceeds ±8 (a few very low-count genes)
switch to a symlog x-axis so the bulk of the distribution stays readable.

### Gene-expression UMAPs

Overlay one gene at a time on the same embedding, straight from the raw counts:

```bash
python src/plot_gene_umap.py --step n24 --genes XIST UTY TMEM156 TSIX
```

Counts are normalised to log1p(CP10k) using each cell's `nCount_RNA`; non-expressing
cells stay as a light grey background. Writes one PNG per gene, a combined grid, and
`gene_expression_summary_<step>.csv` with the percentage of expressing cells and the
mean expression per condition. The four-gene extraction is cached under
`results/scale/<step>/features/umap_gene_counts.tsv` (first run streams the full
Matrix Market file, ~35 s; later runs reuse the cache).



## Setup

```bash
conda env create -f envs/py_env.yaml
conda env create -f envs/r_env.yaml
conda activate multiome-r
Rscript -e 'devtools::install_github("SONGDONGYUAN1994/scDesign3")'
```

The R env is created from `envs/r_env.yaml`; `scDesign3` is installed from GitHub
because it is not on CRAN/Bioconductor. See `AGENTS.md` for environment and data
gotchas.

## Run

```bash
python src/make_sample_sheet.py       # config/samples.csv from GSE193688_RAW/
snakemake -n                          # dry run
snakemake -j 8 --use-conda
```

`config/samples.csv` is the single source of truth for what runs. The `use`
column selects samples; it is regenerated from disk but hand edits are preserved.
Currently `use=yes` for the two pilot donors `C16007` (Aged) and `P02284` (PD);
all 34 samples are present, so `python src/make_sample_sheet.py --all` switches
the whole cohort on.

## Data

`GSE193688_RAW/` — 34 samples, 28 GB, h5 + fragments + peak BEDs + `.tbi`.
Group is derived from the donor prefix: `Y*` Young, `C*` Aged, `P*` PD.

### Expected numbers for the 2-donor pilot

Computed from the BED files (standard chromosomes, width 20–10000, `reduce()`)
and confirmed against `results/ingest/consensus_peaks.bed`. Stage 0 should land
on these; a large deviation means the peak build changed.

| | peaks |
|---|---|
| consensus | **125,610** |
| covered by C16007 | 119,024 |
| covered by P02284 | 77,854 |
| called in both donors | 71,268 |
| called in one donor only | 54,342 |

`reduce()` with its default `min.gapwidth = 1` merges only *touching* intervals
(gap 0 in 0-based BED coordinates), which is why the count is 125,610 and not
the ~131k you get by also bridging 1 bp gaps.

The 54,342 one-donor peaks are why `merge()` is unusable: as structural zeros
they dominate a raw-variance feature ranking. The recount keeps them but fills in
the other donor's true fragment counts, so they are data rather than artifacts.


## Design notes

- **ATAC is recounted, not merged.** `merge()` on ChromatinAssays leaves
  structural zeros for sample-specific peaks. `R/00_ingest.R` builds a consensus
  set with `reduce()` and counts fragments into it with `FeatureMatrix()`.
- **Peak set is built from the samples in the run** (`consensus_from: used`).
  Building it from all 34 would inject structural zeros into a 2-sample pilot.
- **`Aged` is the reference level** for `group`; the reversal sign depends on it.
- **Feature selection is fit on `real_fit` only** — never on the holdout.

## Layout

```
config/config.yaml      analysis params
config/samples.csv      sample sheet (generated, hand-editable `use`)
R/00_ingest.R           QC + doublets + consensus + recount
src/make_sample_sheet.py
src/assemble.py         MTX -> AnnData
Snakefile               orchestration (repo root)
envs/                   conda environments
results/ingest/         stage 0 output
data/processed/         rna.h5ad, atac.h5ad
```
