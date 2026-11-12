#!/usr/bin/env Rscript
# Stage 10a - donor-disjoint simulation for the synthetic-utility proof of concept
#
# The main simulation fits scDesign3 on the fit cells of every donor, so a
# synthetic-trained classifier can recognise donor identity. The utility test needs a
# donor-disjoint design: fit the simulator on training donors only, then train a
# classifier on synthetic cells and test it on held-out real donors.
#
# For each fold (one by default) this script:
#   1. holds out a fraction of donors in every condition,
#   2. fits scDesign3 on the remaining donors,
#   3. simulates `max(sizes)` cells with a condition x celltype design,
#   4. writes one count matrix per simulation seed.
#
# `src/synthetic_utility.py` subsamples these matrices to the requested sizes and
# trains the classifiers.
#
# Outputs (out/):
#   folds.csv                  fold, status, train/holdout donors, cell counts
#   manifest.csv               fold, modality, seed, path, n_features, n_cells
#   fold<f>/obs.csv           celltype/condition of the simulated cells
#   fold<f>/seed<s>_<modality>_counts.mtx.gz
#   fold<f>/fit_<modality>.rds

suppressPackageStartupMessages({
  library(scDesign3)
  library(SingleCellExperiment)
  library(Matrix)
  library(BiocParallel)
})

defaults <- list(
  ingest = "results/ingest",
  features = "results/features",
  obs = "results/cluster/obs_annotated.csv",
  out = "results/utility/sim",
  seed = 20261110,
  sizes = "500,1000,2000,5000,10000",
  sim_seeds = "101,202,303",
  holdout_fraction = 0.34,
  folds = 1,
  n_cores = 8,
  family_use = "nb",
  mu_formula = "celltype + condition + celltype:condition",
  sigma_formula = "1",
  corr_formula = "1"
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
for (k in c("seed", "holdout_fraction", "folds", "n_cores")) args[[k]] <- as.numeric(args[[k]])
sizes <- as.integer(strsplit(args$sizes, ",")[[1]])
sim_seeds <- as.integer(strsplit(args$sim_seeds, ",")[[1]])
max_size <- max(sizes)

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

write_mtx <- function(mat, path) {
  Matrix::writeMM(mat, path)
  R.utils::gzip(path, overwrite = TRUE, remove = TRUE)
  paste0(path, ".gz")
}

# Load only the selected features (see src/subset_features.py). Reading the full
# ATAC matrix for both modalities at once was the main OOM trigger at n24.
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
      stop(modality, ": pre-subset features do not match ", features_dir,
           " (rerun src/subset_features.py --force)")
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

# ------------------------------------------------------------------- inputs --
obs <- read.csv(args$obs, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("cell", "donor", "group", "celltype")
missing <- setdiff(required, colnames(obs))
if (length(missing)) stop("obs is missing: ", paste(missing, collapse = ", "))
obs$celltype <- as.character(obs$celltype)
obs$condition <- factor(obs$group, levels = c("Aged", "PD"))

group_of <- tapply(obs$group, obs$donor, function(x) x[1])
donors <- sort(unique(obs$donor))
donors_by_group <- split(donors, group_of[donors])
if (length(donors_by_group) < 2) {
  stop("utility needs at least two conditions")
}
min_donors <- min(vapply(donors_by_group, length, integer(1)))
if (min_donors < 2) {
  folds <- data.frame(
    fold = integer(0), status = character(0),
    train_donors = character(0), holdout_donors = character(0),
    n_train_cells = integer(0), n_holdout_cells = integer(0), n_celltypes = integer(0),
    stringsAsFactors = FALSE
  )
  utils::write.csv(folds, file.path(args$out, "folds.csv"), row.names = FALSE)
  utils::write.csv(
    data.frame(fold = integer(0), modality = character(0), seed = integer(0),
               path = character(0), n_features = integer(0),
               n_cells = integer(0), stringsAsFactors = FALSE),
    file.path(args$out, "manifest.csv"), row.names = FALSE
  )
  msg("not estimable: at least one condition has fewer than 2 donors")
  quit(save = "no", status = 0)
}

# --------------------------------------------------------------------- fits --
counts <- list()
features_all <- list()
for (modality in c("rna", "atac")) {
  loaded <- load_selected(args$ingest, args$features, modality)
  colnames(loaded$mat) <- obs$cell
  counts[[modality]] <- loaded$mat
  features_all[[modality]] <- loaded$features
  msg(modality, ": ", length(loaded$features), " features (", loaded$source, ")")
}

manifest <- list()
fold_rows <- list()

for (f in seq_len(args$folds)) {
  set.seed(args$seed + f)
  holdout <- unlist(lapply(donors_by_group, function(d) {
    k <- max(1L, floor(length(d) * args$holdout_fraction))
    k <- min(k, length(d) - 1L)
    sample(d, k)
  }), use.names = FALSE)
  train <- setdiff(donors, holdout)

  train_cells <- obs$cell[obs$donor %in% train]
  holdout_cells <- obs$cell[obs$donor %in% holdout]
  if (length(train_cells) < 20 || length(holdout_cells) < 10) {
    msg("fold ", f, ": too few cells, skipping")
    next
  }

  meta <- obs[match(train_cells, obs$cell), ]
  rownames(meta) <- meta$cell
  meta$celltype <- factor(meta$celltype)
  meta$condition <- factor(meta$group, levels = c("Aged", "PD"))

  tab <- table(meta$celltype, meta$condition)
  keep_types <- rownames(tab)[rowSums(tab > 0) == ncol(tab) & rowSums(tab) >= 10]
  if (length(keep_types) < 1) {
    msg("fold ", f, ": no cell type present in both conditions, skipping")
    next
  }
  meta <- meta[as.character(meta$celltype) %in% keep_types, ]
  meta$celltype <- droplevels(meta$celltype)
  meta$condition <- droplevels(meta$condition)
  train_cells <- rownames(meta)

  # synthetic design: celltype from the training cells, condition balanced 50/50
  set.seed(args$seed + 1000 + f)
  syn_celltype <- sample(as.character(meta$celltype), max_size, replace = TRUE)
  syn_condition <- sample(c("Aged", "PD"), max_size, replace = TRUE)
  new_covariate <- data.frame(
    celltype = factor(syn_celltype, levels = levels(meta$celltype)),
    condition = factor(syn_condition, levels = levels(meta$condition)),
    corr_group = 1,
    row.names = sprintf("syn%06d", seq_len(max_size))
  )

  fold_dir <- file.path(args$out, paste0("fold", f))
  dir.create(fold_dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(
    data.frame(cell = rownames(new_covariate),
               celltype = as.character(new_covariate$celltype),
               condition = as.character(new_covariate$condition),
               stringsAsFactors = FALSE),
    file.path(fold_dir, "obs.csv"), row.names = FALSE
  )

  for (modality in c("rna", "atac")) {
    msg("fold ", f, " ", toupper(modality), ": ", length(train_cells),
        " training cells")
    mat <- counts[[modality]][, train_cells, drop = FALSE]
    mat <- as(mat, "CsparseMatrix")
    feats <- rownames(mat)

    sce <- SingleCellExperiment(
      assays = list(counts = mat),
      colData = DataFrame(
        celltype = meta[train_cells, "celltype"],
        condition = meta[train_cells, "condition"],
        row.names = train_cells
      )
    )
    input_data <- construct_data(
      sce = sce,
      assay_use = "counts",
      celltype = "celltype",
      pseudotime = NULL,
      spatial = NULL,
      other_covariates = "condition",
      corr_by = args$corr_formula,
      ncell = ncol(sce)
    )

    fit_path <- file.path(fold_dir, paste0("fit_", modality, ".rds"))
    if (file.exists(fit_path)) {
      msg("reusing ", fit_path)
      fitted <- readRDS(fit_path)
    } else {
      msg("fitting marginal models")
      marginal <- fit_marginal(
        data = input_data,
        predictor = "gene",
        mu_formula = args$mu_formula,
        sigma_formula = args$sigma_formula,
        family_use = args$family_use,
        n_cores = args$n_cores,
        parallelization = "bpmapply",
        BPPARAM = MulticoreParam(workers = args$n_cores)
      )
      msg("fitting copula")
      copula <- fit_copula(
        sce = sce,
        assay_use = "counts",
        input_data = input_data$dat,
        marginal_list = marginal,
        family_use = args$family_use,
        copula = "gaussian",
        n_cores = args$n_cores,
        parallelization = "bpmapply",
        BPPARAM = MulticoreParam(workers = args$n_cores)
      )
      fitted <- list(input_data = input_data, marginal = marginal,
                     copula = copula, features = feats)
      saveRDS(fitted, fit_path)
    }

    seed_paths <- file.path(fold_dir, sprintf("seed%s_%s_counts.mtx", sim_seeds, modality))
    pending <- sim_seeds[!file.exists(paste0(seed_paths, ".gz"))]
    if (length(pending) == 0) {
      msg("  all ", length(sim_seeds), " seeds already simulated; reusing")
    } else {
      msg("  ", length(pending), "/", length(sim_seeds), " seeds to simulate")
    }

    para <- NULL
    if (length(pending) > 0) {
      para <- extract_para(
        sce = sce,
        assay_use = "counts",
        marginal_list = fitted$marginal,
        n_cores = args$n_cores,
        family_use = args$family_use,
        new_covariate = new_covariate,
        parallelization = "bpmapply",
        BPPARAM = MulticoreParam(workers = args$n_cores),
        data = input_data$dat
      )
    }

    for (i in seq_along(sim_seeds)) {
      s <- sim_seeds[[i]]
      path <- seed_paths[[i]]
      if (file.exists(paste0(path, ".gz"))) {
        msg("  seed ", s, " (cached)")
        manifest[[length(manifest) + 1]] <- data.frame(
          fold = f, modality = modality, seed = s, path = paste0(path, ".gz"),
          n_features = length(feats), n_cells = max_size,
          stringsAsFactors = FALSE
        )
        next
      }

      msg("  seed ", s)
      set.seed(s)
      sim <- simu_new(
        sce = sce,
        assay_use = "counts",
        mean_mat = para$mean_mat,
        sigma_mat = para$sigma_mat,
        zero_mat = para$zero_mat,
        copula_list = fitted$copula$copula_list,
        n_cores = args$n_cores,
        family_use = args$family_use,
        input_data = input_data$dat,
        new_covariate = new_covariate,
        important_feature = fitted$copula$important_feature,
        parallelization = "bpmapply",
        BPPARAM = MulticoreParam(workers = args$n_cores),
        filtered_gene = input_data$filtered_gene
      )
      sim <- as(sim, "CsparseMatrix")
      if (nrow(sim) != length(feats) || ncol(sim) != max_size) {
        stop("unexpected simulated shape: ", paste(dim(sim), collapse = " x "))
      }
      rownames(sim) <- feats
      colnames(sim) <- rownames(new_covariate)
      path <- write_mtx(sim, path)
      manifest[[length(manifest) + 1]] <- data.frame(
        fold = f, modality = modality, seed = s, path = path,
        n_features = nrow(sim), n_cells = ncol(sim),
        stringsAsFactors = FALSE
      )
    }
  }

  fold_rows[[length(fold_rows) + 1]] <- data.frame(
    fold = f, status = "ok",
    train_donors = paste(sort(train), collapse = ";"),
    holdout_donors = paste(sort(holdout), collapse = ";"),
    n_train_cells = length(train_cells),
    n_holdout_cells = length(holdout_cells),
    n_celltypes = length(keep_types),
    stringsAsFactors = FALSE
  )
}

if (length(fold_rows) == 0) {
  folds <- data.frame(
    fold = integer(0), status = character(0),
    train_donors = character(0), holdout_donors = character(0),
    n_train_cells = integer(0), n_holdout_cells = integer(0), n_celltypes = integer(0),
    stringsAsFactors = FALSE
  )
} else {
  folds <- do.call(rbind, fold_rows)
}
utils::write.csv(folds, file.path(args$out, "folds.csv"), row.names = FALSE)

if (length(manifest) == 0) {
  manifest <- data.frame(
    fold = integer(0), modality = character(0), seed = integer(0),
    path = character(0), n_features = integer(0), n_cells = integer(0),
    stringsAsFactors = FALSE
  )
} else {
  manifest <- do.call(rbind, manifest)
}
utils::write.csv(manifest, file.path(args$out, "manifest.csv"), row.names = FALSE)
msg("done: ", nrow(folds), " fold(s), ", nrow(manifest), " simulation(s)")
