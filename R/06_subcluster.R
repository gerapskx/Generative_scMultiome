#!/usr/bin/env Rscript
# Stage 14 - k=2 subclustering for the subtype transfer experiment
#
# This is deliberately separate from the main cell-type clustering. It takes one
# target cell type (Oligodendrocyte or Microglia), one modality (RNA or ATAC), and
# the donor-disjoint split from `src/subtype_split.py`.
#
# Design
#   * cells are filtered on the ingest QC table before any subtype is defined:
#       - scDblFinder score below `--max-doublet-score` (default 0.3);
#       - for Microglia, microglial marker score > oligodendrocyte marker score;
#       - for Oligodendrocyte, oligodendrocyte marker score > microglial score.
#     This removes the doublet / ambient-RNA axis that otherwise dominates a
#     two-cluster solution in these two cell types.
#   * normalization / PCA / LSI / Harmony are fit on TRAINING donors only;
#   * k-means k=2 is fit on TRAINING donors only;
#   * held-out donor cells are assigned to a subtype by a classifier trained on
#     the training cells (multinomial logistic regression, nearest-centroid
#     fallback if the classifier collapses to one class);
#   * no held-out donor is used to define the subtype axis.
#   * RNA and ATAC embeddings are each corrected for sequencing depth by
#     regressing log library size and log detected features out of the PCA/LSI
#     scores on training cells, then applying that correction to all cells.
#     This is necessary because the raw two-cluster axis is otherwise a
#     depth / doublet gradient rather than a biological subtype axis.
#   * Harmony dimensions are winsorised at the training 0.1 / 99.9 percentiles
#     before k-means. This stops one outlying cell from becoming its own
#     cluster (which otherwise happens in Microglia ATAC).
#
# RNA:  log1p CPM -> top variable genes within the cell type -> PCA -> depth
#       correction -> Harmony by donor
# ATAC: TF-IDF -> top variable peaks within the cell type -> LSI (drop LSI 1)
#       -> depth correction -> Harmony by donor
#
# Outputs (out/<celltype>/<modality>/):
#   labels.csv          cell, donor, group, split, subtype, confidence, source
#   embedding.csv       cell, dim1..dimN, subtype, split
#   markers.csv         feature, mean_s1, mean_s2, log2fc, pct_s1, pct_s2, score
#   features.csv        feature, rank, variance (the variable features used)
#   filtered_cells.csv  cell, donor, group, split, subtype
#   qc.json             subtype sizes, donor/condition tables, silhouette,
#                        depth/doublet/condition/donor controls, filters
#
# Run from the project root:
#   Rscript R/06_subcluster.R --celltype Microglia --modality rna

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(harmony)
  library(jsonlite)
  library(MASS)
  library(nnet)
  library(cluster)
})

defaults <- list(
  ingest = "results/ingest",
  features = "results/features",
  obs = "results/cluster/obs_annotated.csv",
  donors = "results/subtype/donors.csv",
  qc = "",
  out = "results/subtype",
  celltype = "Microglia",
  modality = "rna",
  k = 2,
  n_features = 1000,
  dims = 30,
  seed = 20261112,
  min_cells_per_subtype = 50,
  max_silhouette_cells = 5000,
  max_doublet_score = 0.3,
  min_cells_per_donor = 10,
  marker_filter = "auto",
  depth_adjust = "true"
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

truthy <- function(x) {
  x <- tolower(trimws(as.character(x)))
  !is.na(x) & x %in% c("yes", "true", "t", "1", "y")
}

args <- parse_args(defaults)
for (k in c("k", "n_features", "dims", "seed", "min_cells_per_subtype",
             "max_silhouette_cells", "max_doublet_score",
             "min_cells_per_donor")) {
  args[[k]] <- as.numeric(args[[k]])
}
if (!args$modality %in% c("rna", "atac")) {
  stop("--modality must be rna or atac")
}
if (args$qc == "") args$qc <- file.path(args$ingest, "obs.csv")

dir.create(args$out, showWarnings = FALSE, recursive = TRUE)
out_dir <- file.path(args$out, args$celltype, args$modality)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

read_mtx <- function(stem) {
  gz <- paste0(stem, ".gz")
  path <- if (file.exists(gz)) gz else stem
  if (!file.exists(path)) stop("missing matrix: ", path)
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  Matrix::readMM(con)
}

load_selected <- function(ingest, features_dir, modality) {
  sel_gz <- file.path(ingest, paste0(modality, "_counts_sel.mtx.gz"))
  sel_tsv <- file.path(ingest, paste0(modality, "_features_sel.tsv"))
  if (file.exists(sel_gz) && file.exists(sel_tsv)) {
    mat <- read_mtx(file.path(ingest, paste0(modality, "_counts_sel.mtx")))
    feats <- readLines(sel_tsv)
    if (nrow(mat) != length(feats)) {
      stop(modality, ": selected matrix has ", nrow(mat),
           " rows but ", length(feats), " feature names")
    }
    want <- readLines(file.path(features_dir, paste0(modality, "_features.tsv")))
    if (!identical(feats, want)) {
      stop(modality, ": pre-subset features do not match ", features_dir)
    }
    rownames(mat) <- feats
    return(list(mat = mat, features = feats, source = "pre-subset"))
  }
  mat <- read_mtx(file.path(ingest, paste0(modality, "_counts.mtx")))
  all_feats <- readLines(file.path(ingest, paste0(modality, "_features.tsv")))
  rownames(mat) <- all_feats
  want <- readLines(file.path(features_dir, paste0(modality, "_features.tsv")))
  keep <- want[want %in% all_feats]
  if (length(keep) < 10) stop(modality, ": too few selected features present")
  list(mat = mat[keep, , drop = FALSE], features = keep, source = "full matrix")
}

auc <- function(y, p) {
  y <- as.integer(y)
  ok <- is.finite(p)
  y <- y[ok]
  p <- p[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(p)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

centroid_predict <- function(xtr, ytr, xte, levels) {
  cent <- do.call(rbind, lapply(levels, function(l) {
    colMeans(xtr[ytr == l, , drop = FALSE])
  }))
  rownames(cent) <- levels
  d2 <- outer(rowSums(xte^2), rowSums(cent^2), "+") - 2 * xte %*% t(cent)
  pred <- levels[apply(d2, 1, which.min)]
  conf <- 1 / (1 + sqrt(apply(d2, 1, min)))
  list(pred = pred, confidence = conf)
}

marker_sets <- list(
  Microglia = list(
    target = c("P2RY12", "CSF1R"),
    other = c("PLP1", "MBP", "MOBP", "CNP", "OLIG1")
  ),
  Oligodendrocyte = list(
    target = c("PLP1", "MBP", "MOBP", "CNP", "OLIG1"),
    other = c("P2RY12", "CSF1R")
  )
)

# ------------------------------------------------------------------- inputs --
obs <- read.csv(args$obs, stringsAsFactors = FALSE, check.names = FALSE)
donors <- read.csv(args$donors, stringsAsFactors = FALSE, check.names = FALSE)
qc <- read.csv(args$qc, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("celltype", "donor", "group", "n_cells", "split") %in% colnames(donors))) {
  stop("donors file must have celltype, donor, group, n_cells, split")
}
donors <- donors[donors$celltype == args$celltype, , drop = FALSE]
if (nrow(donors) == 0) stop("no donors for ", args$celltype)
if (!all(c("cell", "scDblFinder_score") %in% colnames(qc))) {
  stop("qc table must have cell and scDblFinder_score")
}

loaded <- load_selected(args$ingest, args$features, args$modality)
mat <- loaded$mat
if (ncol(mat) != nrow(obs)) stop("matrix columns != obs rows")
colnames(mat) <- obs$cell
if (anyDuplicated(obs$cell)) stop("duplicated cell ids in obs")
qc <- qc[match(obs$cell, qc$cell), ]
if (any(is.na(qc$scDblFinder_score))) {
  stop("some cells are missing from the QC table")
}

# ----------------------------------------------------------- cell filtering --
base <- obs$celltype == args$celltype
keep <- base
n_before <- sum(base)

dbl_ok <- rep(TRUE, nrow(obs))
if (args$max_doublet_score > 0) {
  dbl_ok <- is.na(qc$scDblFinder_score) |
    qc$scDblFinder_score < args$max_doublet_score
}
keep <- keep & dbl_ok

marker_ok <- rep(TRUE, nrow(obs))
marker_name <- "none"
if (!identical(tolower(args$marker_filter), "none")) {
  ms <- marker_sets[[args$celltype]]
  if (!is.null(ms)) {
    have_t <- ms$target[ms$target %in% rownames(mat)]
    have_o <- ms$other[ms$other %in% rownames(mat)]
    if (length(have_t) >= 1 && length(have_o) >= 1) {
      lib <- pmax(qc$nCount_RNA, 1)
      score <- function(genes) {
        log1p(Matrix::colSums(mat[genes, , drop = FALSE]) / lib * 1e4)
      }
      target_score <- score(have_t)
      other_score <- score(have_o)
      marker_ok <- target_score > other_score
      marker_name <- paste0("target>other:", paste(have_t, collapse = ","))
    }
  }
}
keep <- keep & marker_ok

ct_cells <- obs$cell[keep]
donor_ct <- obs$donor[keep]
if (length(ct_cells) < 40) {
  stop("too few cells after filtering: ", length(ct_cells))
}

donor_counts <- table(donor_ct)
ok_donors <- names(donor_counts)[donor_counts >= args$min_cells_per_donor]
dropped_donors <- setdiff(unique(donors$donor), ok_donors)
donors <- donors[donors$donor %in% ok_donors, , drop = FALSE]
keep <- keep & obs$donor %in% ok_donors
ct_cells <- obs$cell[keep]
donor_ct <- obs$donor[keep]
if (length(ct_cells) < 40) {
  stop("too few cells after dropping small donors: ", length(ct_cells))
}

train_donors <- donors$donor[donors$split == "train"]
test_donors <- donors$donor[donors$split == "test"]
if (length(train_donors) < 2 || length(test_donors) < 1) {
  stop("need at least two training donors and one test donor")
}

train_cells <- obs$cell[keep & obs$donor %in% train_donors]
test_cells <- obs$cell[keep & obs$donor %in% test_donors]
if (length(train_cells) < 40 || length(test_cells) < 20) {
  stop("too few cells: ", length(train_cells), " train / ",
       length(test_cells), " test")
}
msg(args$celltype, "/", toupper(args$modality), ": ",
    length(train_cells), " train cells, ", length(test_cells), " test cells; ",
    length(train_donors), " train donors, ", length(test_donors), " test donors")
if (length(dropped_donors)) {
  msg("dropped donors with < ", args$min_cells_per_donor, " filtered cells: ",
      paste(dropped_donors, collapse = ", "))
}

mat <- mat[, ct_cells, drop = FALSE]
lib <- pmax(qc$nCount_RNA[keep], 1)
nf <- pmax(qc$nFeature_RNA[keep], 1)
lib_atac <- pmax(qc$nCount_ATAC[keep], 1)
nf_atac <- pmax(qc$nFeature_ATAC[keep], 1)

train_idx <- match(train_cells, ct_cells)
test_idx <- match(test_cells, ct_cells)
train_pos <- seq_along(train_cells)
test_pos <- length(train_cells) + seq_along(test_cells)

# ---------------------------------------------------------------------- RNA --
if (args$modality == "rna") {
  obj <- CreateSeuratObject(counts = mat, project = args$celltype,
                            min.cells = 0, min.features = 0)
  obj <- NormalizeData(obj, verbose = FALSE)
  norm <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data"))
  v <- apply(norm[, train_idx, drop = FALSE], 1, var)
  v[!is.finite(v)] <- 0
  keep_features <- rownames(norm)[head(order(v, decreasing = TRUE),
                                       min(args$n_features, sum(v > 0)))]
  if (length(keep_features) < 10) stop("too few variable genes")
  norm_use <- norm[keep_features, , drop = FALSE]
  feature_variance <- v[keep_features]
  msg("RNA: ", length(keep_features), " variable genes")

  pca <- prcomp(t(norm_use[, train_idx, drop = FALSE]),
                 center = TRUE, scale. = FALSE, rank. = args$dims)
  emb_train <- pca$x[, seq_len(args$dims), drop = FALSE]
  emb_test <- scale(t(norm_use[, test_idx, drop = FALSE]),
                    center = pca$center, scale = FALSE) %*%
    pca$rotation[, seq_len(args$dims), drop = FALSE]
  emb_all <- rbind(emb_train, emb_test)
  rownames(emb_all) <- c(train_cells, test_cells)
  loadings <- pca$rotation[, seq_len(args$dims), drop = FALSE]
  norm_marker <- norm
  norm_marker_features <- rownames(norm)
  depth_use <- log(lib)
  depth_features <- log(nf)
  depth_name <- "log nCount_RNA"
} else {
  # --------------------------------------------------------------------- ATAC --
  lib_tf <- Matrix::colSums(mat)
  lib_tf[lib_tf == 0] <- 1
  tf <- mat %*% Matrix::Diagonal(x = 1 / lib_tf)
  rownames(tf) <- rownames(mat)
  colnames(tf) <- colnames(mat)
  idf <- log1p(ncol(mat) / pmax(Matrix::rowSums(mat > 0), 1))
  tfidf <- Matrix::Diagonal(x = idf) %*% tf
  rownames(tfidf) <- rownames(mat)
  colnames(tfidf) <- colnames(mat)

  tr <- tfidf[, train_idx, drop = FALSE]
  mu <- Matrix::rowMeans(tr)
  mu2 <- Matrix::rowMeans(tr^2)
  v <- mu2 - mu^2
  v[!is.finite(v)] <- 0
  keep_features <- rownames(tfidf)[head(order(v, decreasing = TRUE),
                                        min(args$n_features, sum(v > 0)))]
  if (length(keep_features) < 10) stop("too few variable peaks")
  tfidf_use <- tfidf[keep_features, , drop = FALSE]
  feature_variance <- v[keep_features]
  x <- t(as.matrix(tfidf_use))
  rownames(x) <- colnames(mat)
  msg("ATAC: ", length(keep_features), " variable peaks")

  # LSI 1 is depth-correlated; compute dims + 1 and drop it.
  pca <- prcomp(x[train_idx, , drop = FALSE], center = TRUE, scale. = FALSE,
                 rank. = args$dims + 1)
  emb_train <- pca$x[, 2:(args$dims + 1), drop = FALSE]
  emb_test <- scale(x[test_idx, , drop = FALSE], center = pca$center,
                    scale = FALSE) %*%
    pca$rotation[, 2:(args$dims + 1), drop = FALSE]
  emb_all <- rbind(emb_train, emb_test)
  rownames(emb_all) <- c(train_cells, test_cells)
  loadings <- pca$rotation[, 2:(args$dims + 1), drop = FALSE]
  norm_marker <- tfidf_use
  norm_marker_features <- rownames(tfidf_use)
  depth_use <- log(lib_atac)
  depth_features <- log(nf_atac)
  depth_name <- "log nCount_ATAC"
}

# ------------------------------------------------------ depth correction ----
if (truthy(args$depth_adjust)) {
  X <- cbind(1, depth_use[train_idx], depth_features[train_idx])
  B <- solve(crossprod(X), crossprod(X, emb_all[train_pos, , drop = FALSE]))
  Xall <- cbind(1, depth_use, depth_features)
  emb_all <- emb_all - Xall %*% B
  emb_train <- emb_all[train_pos, , drop = FALSE]
  msg("regressed depth out of the embedding")
}

# ------------------------------------------------------------------- harmony --
obj_h <- CreateSeuratObject(
  counts = matrix(0, 1, length(ct_cells), dimnames = list("x", ct_cells)),
  project = args$celltype, min.cells = 0, min.features = 0
)
obj_h$donor <- obs$donor[match(colnames(obj_h), obs$cell)]
obj_h[["pca"]] <- CreateDimReducObject(
  embeddings = emb_all, loadings = loadings,
  assay = DefaultAssay(obj_h), key = "PC_"
)
msg("harmonising over donor")
set.seed(args$seed)
obj_h <- RunHarmony(obj_h, group.by.vars = "donor", reduction.use = "pca",
                     dims.use = seq_len(ncol(emb_all)), project.dim = FALSE,
                     reduction.save = "harmony", verbose = FALSE)
emb_h <- Embeddings(obj_h, "harmony")
emb_h <- emb_h[c(train_cells, test_cells), , drop = FALSE]
if (ncol(emb_h) < 2) stop("harmony embedding has too few dimensions")

# Winsorise each Harmony dimension at the training 0.1 / 99.9 percentiles.
# A single outlying cell otherwise becomes its own k-means cluster and makes the
# k=2 split meaningless; clipping is fit on training cells only.
winsor_ref <- emb_h[train_pos, , drop = FALSE]
for (j in seq_len(ncol(emb_h))) {
  q <- stats::quantile(winsor_ref[, j], c(0.001, 0.999), names = FALSE,
                       na.rm = TRUE)
  emb_h[, j] <- pmin(pmax(emb_h[, j], q[1]), q[2])
}

# ------------------------------------------------------------------ k-means --
set.seed(args$seed)
dims_use <- seq_len(min(20, ncol(emb_h)))
km <- kmeans(emb_h[train_pos, dims_use, drop = FALSE], centers = args$k,
             nstart = 50, iter.max = 100)
lab_train_raw <- km$cluster
sizes <- table(lab_train_raw)
if (length(sizes) < 2) stop("k-means found fewer than two subtypes")
ord <- order(as.integer(sizes), decreasing = TRUE)
map <- setNames(paste0("s", seq_along(ord)), names(sizes)[ord])
subtype_train <- unname(map[as.character(lab_train_raw)])

# Orient the labels by the target cell type's marker score so s2 means the
# marker-high state in every modality and seed. This is a label convention only;
# the subtype itself is still defined by the unsupervised k-means partition.
orient_score <- NULL
orientation <- "size"
marker_ms <- marker_sets[[args$celltype]]
if (args$modality == "rna" && !is.null(marker_ms)) {
  have <- marker_ms$target[marker_ms$target %in% norm_marker_features]
  if (length(have) >= 1) {
    idx <- match(have, norm_marker_features)
    orient_score <- log1p(Matrix::colSums(norm_marker[idx, , drop = FALSE]) /
                           lib * 1e4)
    orientation <- "marker"
  }
}
if (!is.null(orient_score)) {
  train_score <- orient_score[train_idx]
  if (mean(train_score[subtype_train == "s2"]) <
      mean(train_score[subtype_train == "s1"])) {
    subtype_train <- ifelse(subtype_train == "s1", "s2", "s1")
  }
}
if (min(table(subtype_train)) < args$min_cells_per_subtype) {
  stop("subtype smaller than ", args$min_cells_per_subtype, " cells: ",
       paste(names(table(subtype_train)), table(subtype_train), sep = "=",
             collapse = ", "))
}

# --------------------------------------------------- held-out label assignment --
df_train <- data.frame(y = factor(subtype_train),
                       emb_h[train_pos, dims_use, drop = FALSE],
                       check.names = FALSE)
fit <- try(nnet::multinom(y ~ ., data = df_train, trace = FALSE,
                           decay = 1e-4, MaxNWts = 100000), silent = TRUE)
assignment <- "multinom"
if (inherits(fit, "try-error")) {
  assignment <- "centroid_fallback"
  cent <- centroid_predict(emb_h[train_pos, dims_use, drop = FALSE],
                           subtype_train,
                           emb_h[test_pos, dims_use, drop = FALSE],
                           levels(factor(subtype_train)))
  subtype_test <- cent$pred
  confidence_test <- cent$confidence
} else {
  pr <- predict(fit, newdata = data.frame(
    emb_h[test_pos, dims_use, drop = FALSE], check.names = FALSE),
    type = "probs")
  # nnet returns a vector of P(second level) for binary outcomes and a matrix
  # only for 3+ levels. Normalise both shapes to an n x k probability matrix.
  sub_levels <- levels(factor(subtype_train))
  if (is.null(dim(pr))) pr <- cbind(1 - as.numeric(pr), as.numeric(pr))
  if (ncol(pr) != length(sub_levels)) {
    stop("multinom returned ", ncol(pr), " probability columns for ",
         length(sub_levels), " subtypes")
  }
  colnames(pr) <- sub_levels
  subtype_test <- colnames(pr)[apply(pr, 1, which.max)]
  confidence_test <- apply(pr, 1, max)
}
if (length(unique(subtype_test)) < 2) {
  assignment <- paste0(assignment, "+centroid_fallback")
  cent <- centroid_predict(emb_h[train_pos, dims_use, drop = FALSE],
                           subtype_train,
                           emb_h[test_pos, dims_use, drop = FALSE],
                           levels(factor(subtype_train)))
  subtype_test <- cent$pred
  confidence_test <- cent$confidence
}
if (length(unique(subtype_test)) < 2) {
  stop("held-out cells collapse to one subtype after both assignments")
}
subtype_all <- c(subtype_train, subtype_test)
confidence_all <- c(rep(1, length(train_cells)), confidence_test)
source_all <- c(rep("cluster", length(train_cells)),
                 rep(assignment, length(test_cells)))

labels <- data.frame(
  cell = c(train_cells, test_cells),
  donor = obs$donor[match(c(train_cells, test_cells), obs$cell)],
  group = obs$group[match(c(train_cells, test_cells), obs$cell)],
  split = c(rep("train", length(train_cells)), rep("test", length(test_cells))),
  subtype = subtype_all,
  confidence = confidence_all,
  source = source_all,
  stringsAsFactors = FALSE
)
utils::write.csv(labels, file.path(out_dir, "labels.csv"), row.names = FALSE)

emb <- data.frame(cell = rownames(emb_h), emb_h, check.names = FALSE)
emb$subtype <- subtype_all[match(emb$cell, labels$cell)]
emb$split <- labels$split[match(emb$cell, labels$cell)]
utils::write.csv(emb, file.path(out_dir, "embedding.csv"), row.names = FALSE)

features <- data.frame(
  feature = keep_features,
  rank = seq_along(keep_features),
  variance = as.numeric(feature_variance),
  stringsAsFactors = FALSE
)
utils::write.csv(features, file.path(out_dir, "features.csv"), row.names = FALSE)

filtered <- data.frame(
  cell = ct_cells,
  donor = donor_ct,
  group = obs$group[match(ct_cells, obs$cell)],
  split = ifelse(donor_ct %in% train_donors, "train", "test"),
  subtype = subtype_all[match(ct_cells, labels$cell)],
  stringsAsFactors = FALSE
)
utils::write.csv(filtered, file.path(out_dir, "filtered_cells.csv"), row.names = FALSE)

# -------------------------------------------------------------------- markers --
train_cols <- train_idx
s1_cols <- train_cols[subtype_train == "s1"]
s2_cols <- train_cols[subtype_train == "s2"]
marker_mat <- norm_marker
m1 <- Matrix::rowMeans(marker_mat[, s1_cols, drop = FALSE])
m2 <- Matrix::rowMeans(marker_mat[, s2_cols, drop = FALSE])
p1 <- Matrix::rowMeans(marker_mat[, s1_cols, drop = FALSE] > 0)
p2 <- Matrix::rowMeans(marker_mat[, s2_cols, drop = FALSE] > 0)
lfc <- log2((m1 + 1e-6) / (m2 + 1e-6))
score <- abs(lfc) * sqrt(pmax(p1 * (1 - p2), 0))
markers <- data.frame(
  feature = rownames(marker_mat), mean_s1 = m1, mean_s2 = m2,
  log2fc = lfc, pct_s1 = p1, pct_s2 = p2, score = score,
  stringsAsFactors = FALSE
)
markers <- markers[order(markers$score, decreasing = TRUE), ]
markers <- head(markers, 200)
utils::write.csv(markers, file.path(out_dir, "markers.csv"), row.names = FALSE)

# ------------------------------------------------------------------------ QC --
set.seed(args$seed)
n_sil <- min(args$max_silhouette_cells, length(train_pos))
sil_idx <- sample(train_pos, n_sil, replace = FALSE)
sil <- cluster::silhouette(as.integer(factor(subtype_train))[sil_idx],
                            dist(emb_h[sil_idx, dims_use, drop = FALSE]))
silhouette_train <- mean(sil[, 3])

donor_train <- labels$donor[labels$split == "train"]
donor_test <- labels$donor[labels$split == "test"]
group_train <- labels$group[labels$split == "train"]
group_test <- labels$group[labels$split == "test"]

donor_tab <- as.data.frame.matrix(table(donor_train, subtype_train))
donor_tab$donor <- rownames(donor_tab)
donor_tab <- donor_tab[, c("donor", setdiff(colnames(donor_tab), "donor")),
                     drop = FALSE]
cond_tab <- as.data.frame.matrix(table(group_train, subtype_train))
cond_tab$group <- rownames(cond_tab)
cond_tab <- cond_tab[, c("group", setdiff(colnames(cond_tab), "group")),
                     drop = FALSE]

# donor-only control: leave-one-training-donor-out LDA on donor mean embeddings
donor_subtype <- colnames(donor_tab)[
  apply(donor_tab[, -1, drop = FALSE], 1, which.max)
]
names(donor_subtype) <- donor_tab$donor
donor_mean <- rowsum(emb_h[train_pos, seq_len(min(10, ncol(emb_h))),
                      drop = FALSE], group = donor_train)
donor_mean <- donor_mean / as.numeric(table(donor_train)[rownames(donor_mean)])
donor_auc <- NA_real_
if (length(unique(donor_subtype)) >= 2 && nrow(donor_mean) >= 4) {
  post <- rep(NA_real_, nrow(donor_mean))
  for (i in seq_len(nrow(donor_mean))) {
    fit_d <- try(MASS::lda(factor(donor_subtype[-i]) ~ .,
                           data = data.frame(donor_mean[-i, , drop = FALSE],
                                             check.names = FALSE)),
                  silent = TRUE)
    if (inherits(fit_d, "try-error")) next
    pr <- try(predict(fit_d, newdata = data.frame(
      donor_mean[i, , drop = FALSE], check.names = FALSE))$posterior[, 2],
      silent = TRUE)
    if (!inherits(pr, "try-error")) post[i] <- pr
  }
  donor_auc <- auc(donor_subtype == "s2", post)
}

# condition-only control: does the subtype predict condition on held-out cells?
condition_auc <- NA_real_
group_levels <- sort(unique(obs$group))
y_cond_train <- as.integer(factor(group_train, levels = group_levels)) - 1
x_sub_train <- as.integer(factor(subtype_train, levels = c("s1", "s2"))) - 1
x_sub_test <- as.integer(factor(subtype_test, levels = c("s1", "s2"))) - 1
fit_c <- try(glm(y_cond_train ~ x_sub_train, family = binomial()),
              silent = TRUE)
if (!inherits(fit_c, "try-error")) {
  pr <- predict(fit_c, newdata = data.frame(x_sub_train = x_sub_test),
                 type = "response")
  condition_auc <- auc(as.integer(factor(group_test, levels = group_levels)) - 1,
                        pr)
}

# technical controls: can the subtype be recovered from depth or doublet score?
depth_auc_train <- auc(subtype_train == "s2", depth_use[train_idx])
depth_auc_test <- auc(subtype_test == "s2", depth_use[test_idx])
doublet_auc_train <- auc(subtype_train == "s2",
                          qc$scDblFinder_score[train_idx])
doublet_auc_test <- auc(subtype_test == "s2",
                         qc$scDblFinder_score[test_idx])

marker_score <- function(genes) {
  have <- genes[genes %in% norm_marker_features]
  if (length(have) == 0) return(NULL)
  idx <- match(have, norm_marker_features)
  log1p(Matrix::colSums(norm_marker[idx, , drop = FALSE]) / lib * 1e4)
}
marker_score_train <- NULL
marker_score_test <- NULL
if (args$modality == "rna" && !is.null(marker_ms)) {
  target <- marker_score(marker_ms$target)
  other <- marker_score(marker_ms$other)
  if (!is.null(target)) {
    marker_score_train <- list(
      target = tapply(target[train_idx], subtype_train, mean),
      other = if (is.null(other)) NULL else tapply(other[train_idx],
                                                  subtype_train, mean)
    )
    marker_score_test <- list(
      target = tapply(target[test_idx], subtype_test, mean),
      other = if (is.null(other)) NULL else tapply(other[test_idx],
                                                  subtype_test, mean)
    )
  }
}

qc_out <- list(
  celltype = args$celltype,
  modality = args$modality,
  k = args$k,
  seed = args$seed,
  n_features = length(keep_features),
  dims = ncol(emb_h),
  depth_adjust = truthy(args$depth_adjust),
  winsor_quantiles = c(0.001, 0.999),
  n_cells_before_filter = n_before,
  n_cells_after_filter = length(ct_cells),
  n_dropped_doublet = n_before - sum(base & dbl_ok),
  n_dropped_marker = sum(base & dbl_ok) - sum(base & dbl_ok & marker_ok),
  n_dropped_donor = sum(base & dbl_ok & marker_ok) - length(ct_cells),
  max_doublet_score = args$max_doublet_score,
  marker_filter = marker_name,
  min_cells_per_donor = args$min_cells_per_donor,
  dropped_donors = as.list(dropped_donors),
  orientation = orientation,
  n_train_cells = length(train_cells),
  n_test_cells = length(test_cells),
  n_train_donors = length(train_donors),
  n_test_donors = length(test_donors),
  subtype_sizes_train = as.list(table(subtype_train)),
  subtype_sizes_test = as.list(table(subtype_test)),
  silhouette_train = silhouette_train,
  assignment = assignment,
  donor_table_train = split(donor_tab, seq_len(nrow(donor_tab))),
  condition_table_train = split(cond_tab, seq_len(nrow(cond_tab))),
  donor_auc = donor_auc,
  condition_auc = condition_auc,
  depth_auc_train = depth_auc_train,
  depth_auc_test = depth_auc_test,
  depth_covariate = depth_name,
  doublet_auc_train = doublet_auc_train,
  doublet_auc_test = doublet_auc_test,
  marker_score_train = marker_score_train,
  marker_score_test = marker_score_test,
  n_train_donors_mixed = sum(apply(donor_tab[, -1, drop = FALSE], 1,
                                   function(x) sum(x > 0) > 1)),
  n_test_donors_with_both = sum(apply(table(donor_test, subtype_test), 1,
                                       function(x) sum(x > 0) > 1)),
  confidence_mean_test = mean(confidence_test),
  confidence_q05_test = unname(quantile(confidence_test, 0.05)),
  confidence_q95_test = unname(quantile(confidence_test, 0.95))
)
jsonlite::write_json(qc_out, file.path(out_dir, "qc.json"), auto_unbox = TRUE,
                      pretty = TRUE, digits = NA)

msg("wrote ", out_dir, "/labels.csv, embedding.csv, markers.csv, features.csv, ",
    "filtered_cells.csv, qc.json")
msg("subtype sizes (train): ", paste(names(table(subtype_train)),
    table(subtype_train), sep = "=", collapse = ", "))
msg("subtype sizes (test): ", paste(names(table(subtype_test)),
    table(subtype_test), sep = "=", collapse = ", "))
