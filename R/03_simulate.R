#!/usr/bin/env Rscript
# Stage 6 - scDesign3 simulation
#
# Fits one scDesign3 model per modality on the fit cells and reuses it for every
# mode/seed. RNA and ATAC are fitted independently: scDesign3 has no joint
# cross-modality copula, so cells are paired only through shared covariates.
#
# Modes
#   alt   : the fitted condition effect is kept (the observed world)
#   null  : every `condition` coefficient in the marginal models is zeroed, so
#           the simulation has no condition effect. Used by the hypothesis module.
#
# The simulated cells carry the covariates of the holdout cells
# (simulation.match_composition: true), so the fidelity discriminator compares
# real and synthetic cells with the same cell-type/condition composition.
#
# Run from the project root:
#   Rscript R/03_simulate.R
#
# Outputs (results/simulation/):
#   fit_<modality>.rds                     reusable fitted objects
#   <modality>/<mode>/obs.csv              covariates of the simulated cells
#   <modality>/<mode>/seed<seed>_counts.mtx.gz
#   manifest.csv

suppressPackageStartupMessages({
  library(scDesign3)
  library(SingleCellExperiment)
  library(Matrix)
  library(BiocParallel)
})

defaults <- list(
  ingest = "results/ingest",
  features = "results/features",
  split = "results/split/split.csv",
  obs = "results/cluster/obs_annotated.csv",
  out = "results/simulation",
  family_use = "nb",
  mu_formula = "celltype + condition + celltype:condition",
  sigma_formula = "1",
  corr_formula = "1",
  n_cores = 4,
  seeds = "101,202,303,404,505",
  modes = "alt,null",
  match_composition = "true"
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
args$n_cores <- as.numeric(args$n_cores)
seeds <- as.integer(strsplit(args$seeds, ",")[[1]])
modes <- strsplit(args$modes, ",")[[1]]
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

# Load only the selected features. The pre-subset copy (written by
# src/subset_features.py) avoids reading the full ATAC matrix, which is large
# enough to trigger the OOM killer. Falls back to the full matrix if absent.
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
split <- read.csv(args$split, stringsAsFactors = FALSE, check.names = FALSE)
fit_cells <- split$cell[split$split == "fit"]
holdout_cells <- split$cell[split$split == "holdout"]
if (length(fit_cells) == 0) stop("no fit cells")

meta <- obs[match(c(fit_cells, holdout_cells), obs$cell), ]
rownames(meta) <- meta$cell
meta$celltype <- factor(meta$celltype)
meta$condition <- factor(meta$group, levels = c("Aged", "PD"))

# every celltype must be present in every condition, otherwise the interaction
# term is unidentifiable and gamlss will fail
tab <- table(meta$celltype, meta$condition)
keep_types <- rownames(tab)[rowSums(tab > 0) == ncol(tab)]
if (length(keep_types) == 0) stop("no cell type is present in both conditions")
dropped <- setdiff(rownames(tab), keep_types)
if (length(dropped)) {
  msg("dropping cell types absent from a condition (unidentifiable interaction): ",
      paste(dropped, collapse = ", "))
}
fit_cells <- fit_cells[as.character(meta[fit_cells, "celltype"]) %in% keep_types]
holdout_cells <- holdout_cells[as.character(meta[holdout_cells, "celltype"]) %in% keep_types]
if (length(fit_cells) < 10) stop("too few fit cells after dropping unidentifiable cell types")

target_cells <- if (truthy(args$match_composition)) holdout_cells else fit_cells
msg(length(fit_cells), " fit cells; ", length(target_cells), " simulated cells; ",
    length(keep_types), " cell types")

new_covariate <- data.frame(
  celltype = factor(as.character(meta[target_cells, "celltype"]), levels = levels(meta$celltype)),
  condition = factor(as.character(meta[target_cells, "condition"]), levels = levels(meta$condition)),
  corr_group = 1,
  row.names = target_cells
)

# ---------------------------------------------------------------- simulation --
manifest <- list()
for (modality in c("rna", "atac")) {
  msg("=== ", toupper(modality), " ===")
  loaded <- load_selected(args$ingest, args$features, modality)
  counts <- loaded$mat
  keep_features <- loaded$features
  colnames(counts) <- obs$cell

  if (length(keep_features) < 10) stop(modality, ": too few selected features present")
  msg(length(keep_features), " selected features (", loaded$source, ")")

  fit_mat <- counts[, fit_cells, drop = FALSE]
  fit_mat <- as(fit_mat, "CsparseMatrix")
  rm(counts)
  invisible(gc())

  sce <- SingleCellExperiment(
    assays = list(counts = fit_mat),
    colData = DataFrame(
      celltype = meta[fit_cells, "celltype"],
      condition = meta[fit_cells, "condition"],
      row.names = fit_cells
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

  fit_path <- file.path(args$out, paste0("fit_", modality, ".rds"))
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
    fitted <- list(input_data = input_data, marginal = marginal, copula = copula,
                   features = keep_features)
    saveRDS(fitted, fit_path)
  }

  if (!identical(fitted$features, keep_features)) {
    stop(modality, ": cached fit uses a different feature set; delete ", fit_path)
  }

  for (mode in modes) {
    msg("mode ", mode)
    mode_dir <- file.path(args$out, modality, mode)
    dir.create(mode_dir, showWarnings = FALSE, recursive = TRUE)

    seed_paths <- file.path(mode_dir, sprintf("seed%s_counts.mtx", seeds))
    pending <- seeds[!file.exists(paste0(seed_paths, ".gz"))]
    if (length(pending) == 0) {
      msg("  all ", length(seeds), " seeds already simulated; reusing")
    } else {
      msg("  ", length(pending), "/", length(seeds), " seeds to simulate")
    }

    para <- NULL
    if (length(pending) > 0) {
      marginal_use <- fitted$marginal
      if (mode == "null") {
        marginal_use <- lapply(marginal_use, function(g) {
          if (is.null(g) || is.null(g$fit)) return(g)
          fit <- g$fit
          # scDesign3 uses a plain glm/glm.nb for simple formulas and gamlss for
          # flexible ones; the condition coefficients live in different slots.
          if (inherits(fit, "gamlss")) {
            cf <- fit$mu.coefficients
            if (!is.null(cf)) {
              hit <- grepl("condition", names(cf))
              if (any(hit)) {
                cf[hit] <- 0
                fit$mu.coefficients <- cf
              }
            }
          } else if (inherits(fit, "glm")) {
            cf <- fit$coefficients
            hit <- grepl("condition", names(cf))
            if (any(hit)) {
              cf[hit] <- 0
              fit$coefficients <- cf
            }
          }
          g$fit <- fit
          g
        })
      } else if (mode != "alt") {
        stop("unknown mode: ", mode)
      }

      para <- extract_para(
        sce = sce,
        assay_use = "counts",
        marginal_list = marginal_use,
        n_cores = args$n_cores,
        family_use = args$family_use,
        new_covariate = new_covariate,
        parallelization = "bpmapply",
        BPPARAM = MulticoreParam(workers = args$n_cores),
        data = input_data$dat
      )
    }

    utils::write.csv(
      data.frame(cell = target_cells, celltype = new_covariate$celltype,
                 condition = new_covariate$condition, stringsAsFactors = FALSE),
      file.path(mode_dir, "obs.csv"), row.names = FALSE
    )

    for (i in seq_along(seeds)) {
      seed <- seeds[[i]]
      path <- seed_paths[[i]]
      if (file.exists(paste0(path, ".gz"))) {
        msg("  seed ", seed, " (cached)")
        manifest[[length(manifest) + 1]] <- data.frame(
          modality = modality, mode = mode, seed = seed, path = paste0(path, ".gz"),
          n_features = length(keep_features), n_cells = length(target_cells),
          stringsAsFactors = FALSE
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
        input_data = input_data$dat,
        new_covariate = new_covariate,
        important_feature = fitted$copula$important_feature,
        parallelization = "bpmapply",
        BPPARAM = MulticoreParam(workers = args$n_cores),
        filtered_gene = input_data$filtered_gene
      )
      sim <- as(sim, "CsparseMatrix")
      if (nrow(sim) != length(keep_features) || ncol(sim) != length(target_cells)) {
        stop("unexpected simulated shape: ", paste(dim(sim), collapse = " x "))
      }
      rownames(sim) <- keep_features
      colnames(sim) <- sprintf("sim%s_%s", seed, seq_len(ncol(sim)))
      path <- write_mtx(sim, path)
      manifest[[length(manifest) + 1]] <- data.frame(
        modality = modality, mode = mode, seed = seed, path = path,
        n_features = nrow(sim), n_cells = ncol(sim),
        stringsAsFactors = FALSE
      )
    }
  }
}

utils::write.csv(do.call(rbind, manifest), file.path(args$out, "manifest.csv"),
                 row.names = FALSE)
msg("done")
