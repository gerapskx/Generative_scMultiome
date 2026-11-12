# Subtype transfer — summary artifacts

These are the small, reviewable outputs of the n24 subtype-transfer experiment.
The full result tree (`results/scale/n24/subtype/`, ~2.4 GB of count matrices and
fits) is intentionally not committed.

## Contents

| Path | What it is |
|---|---|
| `summary.json` | Headline means, direction asymmetry, row counts |
| `comparison.csv` | Long table: QC, concordance, Jaccard, separability and transfer |
| `transfer/transfer.csv` | Raw transfer rows (cell type × modality × model × domain × seed × metric × stratum) |
| `transfer/transfer_summary.csv` | Mean/SD per grouping |
| `transfer/jaccard.csv` | Variable-feature overlap with the fixed n24 set |
| `*/qc.json` | Subclustering QC: silhouette, depth/doublet/condition controls, subtype sizes |
| `*/features.csv` | Variable features used for each cell type/modality |

## Headline transfer AUC (logistic regression, 5 seeds)

| Cell type | Modality | real→real | syn alt→real | real→syn alt | syn null→real | syn null0→real |
|---|---|---:|---:|---:|---:|---:|
| Oligodendrocyte | RNA | 0.817 | 0.862 | 0.852 | 0.497 | 0.505 |
| Oligodendrocyte | ATAC | 0.884 | 0.903 | 0.950 | 0.499 | 0.498 |
| Microglia | RNA | 0.700 | 0.729 | 0.959 | 0.492 | 0.496 |
| Microglia | ATAC | 0.886 | 0.905 | 1.000 | 0.512 | 0.521 |

## Caveats

- The k=2 splits are weak unsupervised axes (silhouette 0.14–0.26), not validated
  developmental subtypes.
- Microglia RNA/ATAC label agreement is ARI −0.002; Oligodendrocyte is 0.185.
- `donor_only` is estimable only for Oligodendrocyte ATAC (0.333); the other
  combinations are dominated by one subtype at donor level.
- Synthetic cells are generally easier to classify than real cells; the largest
  asymmetry is Microglia RNA (real→syn 0.959 vs syn→real 0.729).

Reproduce with `snakemake --configfile config/scaling/n24.yaml` or the commands in
the main README.
