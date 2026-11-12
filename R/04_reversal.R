#!/usr/bin/env Rscript
# Stage 9a - counterfactual reversal, deterministic model branch
#
# Question: if a PD nucleus were moved to the Aged state, which features would
# move, and in which direction? This is a counterfactual, not an observed
# contrast, so it is answered two independent ways:
#
#   (a) this script: take the scDesign3 marginal fit, evaluate the expected
#       count under the observed condition and under the flipped condition with
#       extract_para(), and rank features by the counterfactual log2 fold change;
#   (b) src/reversal.py: train an Aged-vs-PD classifier and rank features by how
#       much they push a PD cell toward the Aged class.
#
# The two methods share the selected features but not their assumptions: (a) is a
# per-feature marginal model, (b) is a multivariate discriminative model.
#
# Because the fitted mu_formula is `celltype + condition + celltype:condition`,
# the marginal model has no cell-level covariates. The counterfactual mean is
# therefore a group-level quantity (one value per cell type and condition), which
# is exactly what extract_para returns. The per-cell variation is what (b) adds.
#
# Controls
#   reverse_direction   Aged -> PD must mirror PD -> Aged (log2FC sign flips)
#   permuted_condition  condition labels shuffled within cell type, models refit;
#                       the counterfactual effect must collapse toward zero
#
# Only cell types with >= --min-cells-per-stratum cells in BOTH conditions are
# reported. In the 2-donor pilot that is Oligodendrocyte only.
#
# Run from the project root:
#   Rscript R/04_reversal.R
#
# Outputs (results/reversal/):
#   model_delta.csv     feature-level counterfactual effect per modality
#   model_summary.json  target/reference, strata, controls

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
  sim = "results/simulation",
  out = "results/reversal",
  target = "PD",
  reference = "Aged",
  min_cells_per_stratum = 100,
  min_target_cells = 20,
  family_use = "nb",
  mu_formula = "celltype + condition + celltype:condition",
  sigma_formula = "1",
  n_cores = 4,
  seed = 20261110,
  permute = "true"
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
args$min_cells_per_stratum <- as.numeric(args$min_cells_per_stratum)
args$min_target_cells <- as.numeric(args$min_target_cells)
args$n_cores <- as.numeric(args$n_cores)
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
obs$condition <- factor(obs$group, levels = c(args$reference, args$target))
obs$celltype <- factor(obs$celltype)
split_cells <- split$cell[split$split == "holdout"]

tab <- table(obs$celltype, obs$condition)
eligible <- rownames(tab)[apply(tab, 1, min) >= args$min_cells_per_stratum]
if (length(eligible) == 0) {
  stop("no cell type has >= ", args$min_cells_per_stratum,
       " cells in both conditions")
}
msg("eligible cell types (>= ", args$min_cells_per_stratum, "/condition): ",
    paste(eligible, collapse = ", "))

meta <- obs
rownames(meta) <- meta$cell

# ------------------------------------------------------------------- model ----
rows <- list()
summary_list <- list()

for (modality in c("rna", "atac")) {
  msg("=== ", toupper(modality), " ===")
  fit_path <- file.path(args$sim, paste0("fit_", modality, ".rds"))
  if (!file.exists(fit_path)) stop("missing fitted model: ", fit_path,
                                   " (run R/03_simulate.R first)")
  fitted <- readRDS(fit_path)
  features <- fitted$features

  # rebuild the fit sce exactly as R/03_simulate.R did, because extract_para
  # needs the same feature x cell matrix and the same colData
  counts <- read_mtx(file.path(args$ingest, paste0(modality, "_counts.mtx")))
  features_all <- readLines(file.path(args$ingest, paste0(modality, "_features.tsv")))
  rownames(counts) <- features_all
  colnames(counts) <- obs$cell
  fit_cells <- split$cell[split$split == "fit"]
  sce <- SingleCellExperiment(
    assays = list(counts = as(counts[features, fit_cells, drop = FALSE], "CsparseMatrix")),
    colData = DataFrame(
      celltype = meta[fit_cells, "celltype"],
      condition = meta[fit_cells, "condition"],
      row.names = fit_cells
    )
  )

  BPPARAM <- MulticoreParam(workers = args$n_cores)

  # a two-row covariate grid per cell type: observed condition and flipped
  # condition. Rownames are arbitrary; extract_para only uses the values.
  for (ct in eligible) {
    grid <- data.frame(
      celltype = factor(c(ct, ct), levels = levels(obs$celltype)),
      condition = factor(c(args$target, args$reference),
                         levels = c(args$reference, args$target)),
      row.names = c(paste0(ct, "_", args$target), paste0(ct, "_", args$reference))
    )

    para <- extract_para(
      sce = sce, assay_use = "counts", marginal_list = fitted$marginal,
      n_cores = args$n_cores, family_use = args$family_use,
      new_covariate = grid, parallelization = "bpmapply", BPPARAM = BPPARAM,
      data = fitted$input_data$dat
    )
    mu_target <- as.numeric(para$mean_mat[1, ])   # observed PD
    mu_ref <- as.numeric(para$mean_mat[2, ])      # counterfactual Aged

    eps <- 1e-8
    log2fc <- log2((mu_ref + eps) / (mu_target + eps))
    delta_count <- mu_ref - mu_target

    # ---- control: reverse direction (Aged -> PD) ------------------------
    # the model is symmetric, so this must be the exact mirror image
    log2fc_reverse <- log2((mu_target + eps) / (mu_ref + eps))

    # ---- control: permuted condition ------------------------------------
    log2fc_perm <- rep(NA_real_, length(features))
    if (truthy(args$permute)) {
      msg("  ", ct, ": refitting with condition permuted within cell type")
      set.seed(args$seed)
      perm_input <- fitted$input_data
      perm_input$dat$condition <- perm_input$dat$condition[
        ave(seq_len(nrow(perm_input$dat)), perm_input$dat$celltype,
            FUN = function(i) sample(i))
      ]
      perm_marginal <- fit_marginal(
        data = perm_input, predictor = "gene",
        mu_formula = args$mu_formula, sigma_formula = args$sigma_formula,
        family_use = args$family_use, n_cores = args$n_cores,
        parallelization = "bpmapply", BPPARAM = BPPARAM
      )
      para_perm <- extract_para(
        sce = sce, assay_use = "counts", marginal_list = perm_marginal,
        n_cores = args$n_cores, family_use = args$family_use,
        new_covariate = grid, parallelization = "bpmapply", BPPARAM = BPPARAM,
        data = perm_input$dat
      )
      perm_target <- as.numeric(para_perm$mean_mat[1, ])
      perm_ref <- as.numeric(para_perm$mean_mat[2, ])
      log2fc_perm <- log2((perm_ref + eps) / (perm_target + eps))
    }

    n_target <- sum(obs$celltype == ct & obs$condition == args$target)
    n_ref <- sum(obs$celltype == ct & obs$condition == args$reference)
    n_holdout <- sum(split_cells %in% obs$cell[obs$celltype == ct &
                                                obs$condition == args$target])

    rows[[length(rows) + 1]] <- data.frame(
      modality = modality, celltype = ct, feature = features,
      n_target = n_target, n_reference = n_ref, n_target_holdout = n_holdout,
      mu_target = mu_target, mu_reference = mu_ref,
      delta_count = delta_count, log2fc = log2fc,
      log2fc_reverse = log2fc_reverse, log2fc_permuted = log2fc_perm,
      stringsAsFactors = FALSE
    )
    summary_list[[length(summary_list) + 1]] <- list(
      modality = modality, celltype = ct, n_target = n_target,
      n_reference = n_ref, n_target_holdout = n_holdout,
      reverse_correlation = if (sd(log2fc) > 0 && sd(log2fc_reverse) > 0)
        cor(log2fc, log2fc_reverse) else NA_real_,
      observed_sd = sd(log2fc),
      permuted_sd = if (all(is.na(log2fc_perm))) NA_real_
        else sd(log2fc_perm, na.rm = TRUE),
      permuted_q95 = if (all(is.na(log2fc_perm))) NA_real_
        else as.numeric(quantile(abs(log2fc_perm), 0.95, na.rm = TRUE)),
      frac_observed_above_perm_q95 = if (all(is.na(log2fc_perm))) NA_real_
        else mean(abs(log2fc) > quantile(abs(log2fc_perm), 0.95, na.rm = TRUE)),
      permuted_max_abs = if (all(is.na(log2fc_perm))) NA_real_
        else max(abs(log2fc_perm), na.rm = TRUE),
      permuted_correlation = if (all(is.na(log2fc_perm))) NA_real_
        else cor(log2fc, log2fc_perm, use = "complete.obs")
    )
  }
}

model <- do.call(rbind, rows)
model$direction <- ifelse(model$log2fc > 0, "up_in_reference", "down_in_reference")
model$rank <- ave(-abs(model$log2fc), model$modality, model$celltype,
                  FUN = function(x) rank(x, ties.method = "first"))
utils::write.csv(model, file.path(args$out, "model_delta.csv"), row.names = FALSE)

summary <- list(
  target = args$target,
  reference = args$reference,
  min_cells_per_stratum = args$min_cells_per_stratum,
  eligible_celltypes = eligible,
  permuted = truthy(args$permute),
  permute_seed = args$seed,
  strata = summary_list
)
jsonlite::write_json(summary, file.path(args$out, "model_summary.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = NA)
msg("wrote ", nrow(model), " feature rows to ", args$out)
for (s in summary_list) {
  msg("  ", s$modality, " ", s$celltype, ": reverse r = ",
      round(s$reverse_correlation, 3), ", permuted max|log2FC| = ",
      round(s$permuted_max_abs, 3))
}
