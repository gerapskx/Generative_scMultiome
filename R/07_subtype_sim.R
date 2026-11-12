#!/usr/bin/env Rscript
# Stage 15 - scDesign3 simulation for the subtype transfer experiment
#
# Fits a subtype model on the training donors of one target cell type and one
# modality, then simulates balanced 2 x 2 subtype x condition data in three modes:
#
#   alt        : fit `subtype + condition`, keep the subtype effect
#   null       : collapse subtype to one level, refit `condition` (primary null)
#   null_zero  : use the alt fit but zero every `subtype` coefficient
#                (sensitivity null; keeps the alt copula and marginal shape)
#
# The synthetic cells carry a balanced subtype x condition design whose total size
# matches the real training cells (capped at `--max-cells`). The same subtype
# labels are written to obs.csv in every mode, so the transfer stage can evaluate
# a null simulation against the real subtype labels.
#
# Run from the project root:
#   Rscript R/07_subtype_sim.R --celltype Microglia --modality rna

suppressPackageStartupMessages({
  library(scDesign3)
  library(SingleCellExperiment)
  library(Matrix)
  library(BiocParallel)
})

defaults <- list(
  ingest = "results/ingest",
  features = "results/features",
  labels = "results/subtype/Microglia/rna/labels.csv",
  out = "results/subtype/Microglia/rna/simulation",
  celltype = "Microglia",
  modality = "rna",
  family_use = "nb",
  mu_formula_alt = "subtype + condition",
  sigma_formula = "1",
  corr_formula = "1",
  n_cores = 4,
  seeds = "101,202,303,404,505",
  modes = "alt,null,null_zero",
  max_cells = 20000,
  min_cells_per_stratum = 20
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
for (k in c("n_cores", "max_cells", "min_cells_per_stratum")) {
  args[[k]] <- as.numeric(args[[k]])
}
seeds <- as.integer(strsplit(args$seeds, ",")[[1]])
modes <- strsplit(args$modes, ",")[[1]]
if (!args$modality %in% c("rna", "atac")) stop("--modality must be rna or atac")
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

# ------------------------------------------------------------------- inputs --
labels <- read.csv(args$labels, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("cell", "donor", "group", "split", "subtype") %in% colnames(labels))) {
  stop("labels file must have cell, donor, group, split, subtype")
}
fit_rows <- labels$split == "train"
if (sum(fit_rows) < 40) stop("too few training cells: ", sum(fit_rows))

obs <- read.csv(file.path(args$ingest, "obs.csv"), stringsAsFactors = FALSE,
                  check.names = FALSE)
if (!all(c("cell", "donor", "group") %in% colnames(obs))) {
  stop("ingest obs.csv must have cell, donor, group")
}
loaded <- load_selected(args$ingest, args$features, args$modality)
counts <- loaded$mat
if (ncol(counts) != nrow(obs)) stop("matrix columns != ingest obs rows")
colnames(counts) <- obs$cell
counts <- counts[, labels$cell, drop = FALSE]
msg(length(loaded$features), " selected features (", loaded$source, ")")

meta <- labels
meta$subtype <- factor(meta$subtype, levels = c("s1", "s2"))
meta$condition <- factor(meta$group, levels = c("Aged", "PD"))
if (any(is.na(meta$subtype))) stop("subtype must be s1/s2")
if (any(is.na(meta$condition))) stop("group must be Aged/PD")

fit_cells <- meta$cell[fit_rows]
fit_mat <- counts[, fit_cells, drop = FALSE]
fit_mat <- as(fit_mat, "CsparseMatrix")
rm(counts)
invisible(gc())

sce <- SingleCellExperiment(
  assays = list(counts = fit_mat),
  colData = DataFrame(
    subtype = meta$subtype[fit_rows],
    condition = meta$condition[fit_rows],
    row.names = fit_cells
  )
)

input_alt <- construct_data(
  sce = sce,
  assay_use = "counts",
  celltype = "subtype",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = "condition",
  corr_by = args$corr_formula,
  ncell = ncol(sce)
)

# --------------------------------------------------------------- balanced design --
n_total <- min(nrow(labels[fit_rows, , drop = FALSE]), args$max_cells)
strata <- expand.grid(
  subtype = factor(c("s1", "s2"), levels = c("s1", "s2")),
  condition = factor(c("Aged", "PD"), levels = c("Aged", "PD")),
  stringsAsFactors = FALSE
)
n_strata <- nrow(strata)
base <- floor(n_total / n_strata)
rem <- n_total - base * n_strata
per <- rep(base, n_strata)
if (rem > 0) per[seq_len(rem)] <- per[seq_len(rem)] + 1
if (any(per < args$min_cells_per_stratum)) {
  stop("balanced design has fewer than ", args$min_cells_per_stratum,
       " cells per stratum: ", paste(per, collapse = ", "))
}
eval_subtype <- rep(as.character(strata$subtype), per)
eval_condition <- rep(as.character(strata$condition), per)
sim_cells <- sprintf("sim_%s_%d", args$celltype, seq_len(n_total))

new_covariate_alt <- data.frame(
  subtype = factor(eval_subtype, levels = c("s1", "s2")),
  condition = factor(eval_condition, levels = c("Aged", "PD")),
  corr_group = 1,
  row.names = sim_cells
)
new_covariate_null <- data.frame(
  subtype = factor(rep("all", n_total), levels = "all"),
  condition = factor(eval_condition, levels = c("Aged", "PD")),
  corr_group = 1,
  row.names = sim_cells
)
msg(n_total, " synthetic cells per seed; ", paste(per, collapse = "/"),
    " per subtype x condition stratum")

# ---------------------------------------------------------------------- fits --
fit_alt_path <- file.path(args$out, "fit_alt.rds")
if (file.exists(fit_alt_path)) {
  msg("reusing ", fit_alt_path)
  fit_alt <- readRDS(fit_alt_path)
} else {
  msg("fitting alt marginal models")
  marginal_alt <- fit_marginal(
    data = input_alt,
    predictor = "gene",
    mu_formula = args$mu_formula_alt,
    sigma_formula = args$sigma_formula,
    family_use = args$family_use,
    n_cores = args$n_cores,
    parallelization = "bpmapply",
    BPPARAM = MulticoreParam(workers = args$n_cores)
  )
  msg("fitting alt copula")
  copula_alt <- fit_copula(
    sce = sce,
    assay_use = "counts",
    input_data = input_alt$dat,
    marginal_list = marginal_alt,
    family_use = args$family_use,
    copula = "gaussian",
    n_cores = args$n_cores,
    parallelization = "bpmapply",
    BPPARAM = MulticoreParam(workers = args$n_cores)
  )
  fit_alt <- list(input_data = input_alt, marginal = marginal_alt,
                   copula = copula_alt, features = loaded$features)
  saveRDS(fit_alt, fit_alt_path)
}

fit_null_path <- file.path(args$out, "fit_null.rds")
if (file.exists(fit_null_path)) {
  msg("reusing ", fit_null_path)
  fit_null <- readRDS(fit_null_path)
} else {
  msg("fitting null (subtype collapsed) marginal models")
  sce_null <- SingleCellExperiment(
    assays = list(counts = fit_mat),
    colData = DataFrame(
      subtype = factor(rep("all", ncol(fit_mat)), levels = "all"),
      condition = meta$condition[fit_rows],
      row.names = fit_cells
    )
  )
  input_null <- construct_data(
    sce = sce_null,
    assay_use = "counts",
    celltype = "subtype",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = "condition",
    corr_by = args$corr_formula,
    ncell = ncol(sce_null)
  )
  marginal_null <- fit_marginal(
    data = input_null,
    predictor = "gene",
    mu_formula = "condition",
    sigma_formula = args$sigma_formula,
    family_use = args$family_use,
    n_cores = args$n_cores,
    parallelization = "bpmapply",
    BPPARAM = MulticoreParam(workers = args$n_cores)
  )
  msg("fitting null copula")
  copula_null <- fit_copula(
    sce = sce_null,
    assay_use = "counts",
    input_data = input_null$dat,
    marginal_list = marginal_null,
    family_use = args$family_use,
    copula = "gaussian",
    n_cores = args$n_cores,
    parallelization = "bpmapply",
    BPPARAM = MulticoreParam(workers = args$n_cores)
  )
  fit_null <- list(input_data = input_null, marginal = marginal_null,
                   copula = copula_null, features = loaded$features)
  saveRDS(fit_null, fit_null_path)
}

# Null sensitivity: keep the alt copula and marginal family, but zero the
# subtype coefficients in the alt marginal fits.
zero_subtype <- function(marginal) {
  lapply(marginal, function(g) {
    if (is.null(g) || is.null(g$fit)) return(g)
    fit <- g$fit
    if (inherits(fit, "gamlss")) {
      cf <- fit$mu.coefficients
      if (!is.null(cf)) {
        hit <- grepl("subtype", names(cf))
        if (any(hit)) {
          cf[hit] <- 0
          fit$mu.coefficients <- cf
        }
      }
    } else if (inherits(fit, "glm")) {
      cf <- fit$coefficients
      hit <- grepl("subtype", names(cf))
      if (any(hit)) {
        cf[hit] <- 0
        fit$coefficients <- cf
      }
    }
    g$fit <- fit
    g
  })
}

# -------------------------------------------------------------------- simulate --
manifest <- list()
for (mode in modes) {
  if (!mode %in% c("alt", "null", "null_zero")) {
    stop("unknown mode: ", mode)
  }
  msg("mode ", mode)
  mode_dir <- file.path(args$out, mode)
  dir.create(mode_dir, showWarnings = FALSE, recursive = TRUE)

  if (mode == "null") {
    fitted <- fit_null
    new_covariate <- new_covariate_null
    eval_labels <- data.frame(
      cell = sim_cells, subtype = eval_subtype, condition = eval_condition,
      subtype_model = "all", stringsAsFactors = FALSE
    )
  } else {
    fitted <- fit_alt
    new_covariate <- new_covariate_alt
    eval_labels <- data.frame(
      cell = sim_cells, subtype = eval_subtype, condition = eval_condition,
      subtype_model = as.character(new_covariate_alt$subtype),
      stringsAsFactors = FALSE
    )
  }

  marginal_use <- fitted$marginal
  if (mode == "null_zero") marginal_use <- zero_subtype(marginal_use)

  para <- extract_para(
    sce = sce,
    assay_use = "counts",
    marginal_list = marginal_use,
    n_cores = args$n_cores,
    family_use = args$family_use,
    new_covariate = new_covariate,
    parallelization = "bpmapply",
    BPPARAM = MulticoreParam(workers = args$n_cores),
    data = fitted$input_data$dat
  )

  utils::write.csv(eval_labels, file.path(mode_dir, "obs.csv"), row.names = FALSE)

  for (i in seq_along(seeds)) {
    seed <- seeds[[i]]
    path <- file.path(mode_dir, sprintf("seed%s_counts.mtx", seed))
    if (file.exists(paste0(path, ".gz"))) {
      msg("  seed ", seed, " (cached)")
      manifest[[length(manifest) + 1]] <- data.frame(
        modality = args$modality, celltype = args$celltype, mode = mode,
        seed = seed, path = paste0(path, ".gz"),
        n_features = length(loaded$features), n_cells = n_total,
        n_train_cells = sum(fit_rows), stringsAsFactors = FALSE
      )
      next
    }
    msg("  seed ", seed)
    set.seed(seed)
    sim <- simu_new(
      sce = sce,
      assay_use = "counts",
      mean_mat = para$mean_mat,
      sigma_mat = para$sigma_mat,
      zero_mat = para$zero_mat,
      copula_list = fitted$copula$copula_list,
      n_cores = args$n_cores,
      family_use = args$family_use,
      input_data = fitted$input_data$dat,
      new_covariate = new_covariate,
      important_feature = fitted$copula$important_feature,
      parallelization = "bpmapply",
      BPPARAM = MulticoreParam(workers = args$n_cores),
      filtered_gene = fitted$input_data$filtered_gene
    )
    sim <- as(sim, "CsparseMatrix")
    if (nrow(sim) != length(loaded$features) || ncol(sim) != n_total) {
      stop("unexpected simulated shape: ", paste(dim(sim), collapse = " x "))
    }
    rownames(sim) <- loaded$features
    colnames(sim) <- sim_cells
    path <- write_mtx(sim, path)
    manifest[[length(manifest) + 1]] <- data.frame(
      modality = args$modality, celltype = args$celltype, mode = mode,
      seed = seed, path = path, n_features = nrow(sim), n_cells = ncol(sim),
      n_train_cells = sum(fit_rows), stringsAsFactors = FALSE
    )
  }
}

utils::write.csv(do.call(rbind, manifest), file.path(args$out, "manifest.csv"),
                 row.names = FALSE)
msg("wrote ", args$out, "/manifest.csv and ", length(modes), " mode directories")
