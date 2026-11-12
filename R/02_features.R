#!/usr/bin/env Rscript
# Stage 5 - feature selection (fit cells only)
#
# RNA : SCTransform v2 residuals, ranked by residual variance, with a detection
#       floor (fraction of fit cells with a non-zero count).
# ATAC: TF-IDF on binarised counts, ranked by variance, with a detection floor
#       required in EVERY group so a peak private to one condition cannot be
#       selected just because it is all-zero elsewhere.
#
# Both lists are derived from the fit cells only. The holdout and the synthetic
# data are later subset to these same names, so no label leaks into selection.
#
# Run from the project root:
#   Rscript R/02_features.R
#
# Outputs (results/features/):
#   rna_features.tsv     selected genes, ranked
#   atac_features.tsv    selected peaks, ranked
#   rna_feature_stats.csv / atac_feature_stats.csv
#   feature_params.json

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

defaults <- list(
  ingest = "results/ingest",
  split = "results/split/split.csv",
  obs = "results/cluster/obs_annotated.csv",
  out = "results/features",
  n_rna = 2000,
  n_atac = 2000,
  min_fraction_rna = 0.05,
  min_fraction_atac = 0.05,
  seed = 20261110
)

parse_args <- function(defaults) {
  argv <- commandArgs(trailingOnly = TRUE)
  if (length(argv) %% 2 != 0) stop("arguments must be --key value pairs")
  out <- defaults
  if (length(argv)) {
    for (i in seq(1, length(argv), by = 2)) {
      key <- gsub("-", "_", sub("^--", "", argv[i]))
      if (!startsWith(argv[i], "--")) stop("expected --key, got: ", argv[i])
      if (is.null(out[[key]])) stop("unknown argument: --", key)
      out[[key]] <- argv[i + 1]
    }
  }
  out
}

args <- parse_args(defaults)
for (k in c("n_rna", "n_atac", "min_fraction_rna", "min_fraction_atac", "seed")) args[[k]] <- as.numeric(args[[k]])
dir.create(args$out, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

read_mtx <- function(stem) {
  gz <- paste0(stem, ".gz")
  path <- if (file.exists(gz)) gz else stem
  if (!file.exists(path)) stop("missing matrix: ", path)
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  Matrix::readMM(con)
}

# ------------------------------------------------------------------- inputs --
obs <- read.csv(args$obs, stringsAsFactors = FALSE, check.names = FALSE)
split <- read.csv(args$split, stringsAsFactors = FALSE, check.names = FALSE)
fit_cells <- split$cell[split$split == "fit"]
if (length(fit_cells) == 0) stop("no fit cells in ", args$split)

rna <- read_mtx(file.path(args$ingest, "rna_counts.mtx"))
atac <- read_mtx(file.path(args$ingest, "atac_counts.mtx"))
genes <- readLines(file.path(args$ingest, "rna_features.tsv"))
peaks <- readLines(file.path(args$ingest, "atac_features.tsv"))
rownames(rna) <- genes; rownames(atac) <- peaks
colnames(rna) <- obs$cell; colnames(atac) <- obs$cell

fit_cells <- intersect(fit_cells, obs$cell)
if (length(fit_cells) == 0) stop("fit cells do not match obs")
groups <- setNames(obs$group, obs$cell)[fit_cells]
msg(length(fit_cells), " fit cells; groups: ",
    paste(names(table(groups)), table(groups), sep = "=", collapse = " "))

# ---------------------------------------------------------------------- RNA --
msg("RNA: SCTransform residuals on ", nrow(rna), " genes")
set.seed(args$seed)
obj <- CreateSeuratObject(counts = rna[, fit_cells, drop = FALSE])
obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE)

fa <- SCTResults(obj, slot = "feature.attributes")
if (!all(c("residual_variance", "detection_rate") %in% colnames(fa))) {
  stop("unexpected SCTransform feature.attributes: ",
       paste(colnames(fa), collapse = ", "))
}
gene_stats <- data.frame(
  gene = rownames(fa),
  residual_variance = fa$residual_variance,
  detection_rate = fa$detection_rate,
  stringsAsFactors = FALSE
)
gene_stats <- gene_stats[
  is.finite(gene_stats$residual_variance) &
    gene_stats$residual_variance > 0 &
    gene_stats$detection_rate >= args$min_fraction_rna,
]
gene_stats <- gene_stats[order(gene_stats$residual_variance, decreasing = TRUE), ]
rna_features <- head(gene_stats$gene, args$n_rna)
msg("  ", nrow(gene_stats), " genes pass the floor; selected ", length(rna_features))

# --------------------------------------------------------------------- ATAC --
msg("ATAC: TF-IDF variance on ", nrow(atac), " peaks")
atac_fit <- atac[, fit_cells, drop = FALSE]
n_cells <- ncol(atac_fit)
tf <- atac_fit %*% Diagonal(x = 1 / pmax(Matrix::colSums(atac_fit), 1))
idf <- log(1 + n_cells / pmax(Matrix::rowSums(atac_fit > 0), 1))
tfidf <- Diagonal(x = idf) %*% tf

mu <- Matrix::rowMeans(tfidf)
mu2 <- Matrix::rowMeans(tfidf^2)
variance <- mu2 - mu^2
variance[!is.finite(variance)] <- 0

# detection floor in every group present among the fit cells
group_levels <- sort(unique(groups))
frac <- sapply(group_levels, function(g) {
  cells <- fit_cells[groups == g]
  Matrix::rowSums(atac_fit[, cells, drop = FALSE] > 0) / length(cells)
})
if (is.null(dim(frac))) frac <- matrix(frac, ncol = 1)
min_frac <- apply(frac, 1, min)

peak_stats <- data.frame(
  peak = peaks,
  tfidf_variance = variance,
  min_group_detection = min_frac,
  stringsAsFactors = FALSE
)
peak_stats <- peak_stats[
  peak_stats$tfidf_variance > 0 & peak_stats$min_group_detection >= args$min_fraction_atac,
]
peak_stats <- peak_stats[order(peak_stats$tfidf_variance, decreasing = TRUE), ]
atac_features <- head(peak_stats$peak, args$n_atac)
msg("  ", nrow(peak_stats), " peaks pass the floor; selected ", length(atac_features))

# ------------------------------------------------------------------- export --
writeLines(rna_features, file.path(args$out, "rna_features.tsv"))
writeLines(atac_features, file.path(args$out, "atac_features.tsv"))
utils::write.csv(gene_stats, file.path(args$out, "rna_feature_stats.csv"), row.names = FALSE)
utils::write.csv(peak_stats, file.path(args$out, "atac_feature_stats.csv"), row.names = FALSE)
writeLines(
  c(
    "{",
    sprintf('  "n_rna": %d,', length(rna_features)),
    sprintf('  "n_atac": %d,', length(atac_features)),
    sprintf('  "min_fraction_rna": %s,', format(args$min_fraction_rna)),
    sprintf('  "min_fraction_atac": %s,', format(args$min_fraction_atac)),
    sprintf('  "fit_cells": %d,', length(fit_cells)),
    sprintf('  "groups": [%s]', paste(sprintf('"%s"', group_levels), collapse = ", ")),
    "}"
  ),
  file.path(args$out, "feature_params.json")
)
msg("wrote ", length(rna_features), " genes and ", length(atac_features), " peaks to ", args$out)
