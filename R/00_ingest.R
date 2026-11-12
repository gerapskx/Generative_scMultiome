#!/usr/bin/env Rscript
# Stage 0 - ingest
#
# Reads the 10x Multiome h5 files listed in config/samples.csv, applies RNA QC
# and per-sample doublet detection, builds a consensus peak set by reduce(), and
# recounts ATAC fragments on that consensus with FeatureMatrix() so the exported
# matrix has true zeros instead of the structural zeros produced by merge().
#
# Run from the project root:
#   Rscript R/00_ingest.R
#
# Outputs (results/ingest/):
#   rna_counts.mtx[.gz]   genes x cells, raw integer counts
#   atac_counts.mtx[.gz]  consensus peaks x cells, raw integer counts
#   rna_features.tsv      gene names, one per line, matching mtx rows
#   atac_features.tsv     peak names (chr1-100-200), matching mtx rows
#   obs.csv               per-cell metadata, matching mtx columns, in order
#   qc_per_sample.csv     per-sample cell accounting
#   consensus_peaks.bed   the peak set used
#   sessionInfo.txt

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(Rsamtools)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(EnsDb.Hsapiens.v86)
})

# ---------------------------------------------------------------- arguments --
defaults <- list(
  samples = "config/samples.csv",
  out = "results/ingest",
  min_features_rna = 200,
  max_features_rna = 6000,
  min_counts_rna = 500,
  max_percent_mt = 5,
  min_counts_atac = 1000,
  peak_min_width = 20,
  peak_max_width = 10000,
  consensus_from = "used",
  doublets = "true",
  gzip = "true",
  validate_fragments = "false"
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
for (k in c("min_features_rna", "max_features_rna", "min_counts_rna",
            "max_percent_mt", "min_counts_atac", "peak_min_width",
            "peak_max_width")) {
  args[[k]] <- as.numeric(args[[k]])
}

dir.create(args$out, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

# ------------------------------------------------------------------ samples --
sheet <- read.csv(args$samples, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("sample", "gsm", "donor", "group", "h5", "fragments", "peaks", "use")
missing <- setdiff(required, colnames(sheet))
if (length(missing)) stop("samples sheet missing columns: ", paste(missing, collapse = ", "))

is_used <- truthy(tolower(as.character(sheet$use)))
sheet$use <- ifelse(is_used, "yes", "no")
used <- sheet[is_used, , drop = FALSE]
if (nrow(used) == 0) stop("no samples with use=yes in ", args$samples)

groups <- sort(unique(used$group))
if (length(groups) < 2) warning("only one group in this run: ", paste(groups, collapse = ", "))
msg(nrow(used), " sample(s) selected: ", paste(used$donor, collapse = ", "))

for (col in c("h5", "fragments", "peaks")) {
  absent <- used[[col]][!file.exists(used[[col]])]
  if (length(absent)) stop("missing ", col, ": ", paste(absent, collapse = ", "))
}
if (anyDuplicated(used$sample)) stop("duplicated sample ids in ", args$samples)

# -------------------------------------------------------------- annotations --
msg("loading EnsDb annotations")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"   # 10x peaks are chr1:100-200
genome(annotations) <- "hg38"

# ---------------------------------------------------------- consensus peaks --
peak_source <- if (identical(args$consensus_from, "all")) sheet else used
msg("building consensus peaks from ", nrow(peak_source), " sample(s)")
peak_list <- lapply(peak_source$peaks, function(f) {
  gr <- rtracklayer::import(f)
  GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
})
consensus <- reduce(do.call(c, peak_list))
consensus <- consensus[
  width(consensus) >= args$peak_min_width &
    width(consensus) <= args$peak_max_width
]
msg("consensus peaks: ", length(consensus))
if (length(consensus) == 0) stop("empty consensus peak set")
export(consensus, file.path(args$out, "consensus_peaks.bed"), format = "BED")

# ------------------------------------------------------------------- ingest --
rna_list <- list()
atac_list <- list()
obs_list <- list()
qc_rows <- list()

# Per-sample checkpoints make the recount resumable. The fragment recount is the
# expensive part, so losing a 3-hour run to a restart is worth avoiding.
ckpt_dir <- file.path(args$out, "checkpoints")
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)
expected_peaks <- paste0(as.character(seqnames(consensus)), "-",
                         start(consensus), "-", end(consensus))

for (i in seq_len(nrow(used))) {
  s <- used[i, ]
  ckpt <- file.path(ckpt_dir, paste0(s$sample, ".rds"))
  if (file.exists(ckpt)) {
    msg("reusing checkpoint ", s$sample)
    saved <- readRDS(ckpt)
    if (!identical(saved$peak_names, expected_peaks)) {
      stop("checkpoint ", ckpt,
           " was built on a different peak set; delete ", ckpt_dir)
    }
    rna_list[[s$sample]] <- saved$rna
    atac_list[[s$sample]] <- saved$atac
    obs_list[[s$sample]] <- saved$obs
    qc_rows[[s$sample]] <- saved$qc
    next
  }
  msg("processing ", s$sample)

  counts <- Read10X_h5(s$h5)
  if (!all(c("Gene Expression", "Peaks") %in% names(counts))) {
    stop(s$sample, ": h5 is not a multiome matrix")
  }

  # --- RNA ---
  obj <- CreateSeuratObject(
    counts = counts[["Gene Expression"]],
    project = s$sample,
    min.cells = 0,
    min.features = 0
  )
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- RenameCells(obj, add.cell.id = s$sample)

  # --- doublets (per sample, on raw counts) ---
  dbl_class <- rep("singlet", ncol(obj))
  dbl_score <- rep(NA_real_, ncol(obj))
  if (truthy(args$doublets)) {
    set.seed(1)
    dbl <- tryCatch(
      scDblFinder(SingleCellExperiment(
        list(counts = GetAssayData(obj, layer = "counts"))
      )),
      error = function(e) {
        warning(s$sample, ": scDblFinder failed (", conditionMessage(e),
                "); all cells treated as singlets", call. = FALSE)
        NULL
      }
    )
    if (!is.null(dbl)) {
      dbl_class <- as.character(dbl$scDblFinder.class)
      dbl_score <- dbl$scDblFinder.score
    }
  }

  # --- ATAC recount on the consensus set ---
  frag <- s$fragments
  if (!file.exists(paste0(frag, ".tbi"))) {
    msg("  indexing ", basename(frag))
    Rsamtools::indexTabix(frag, format = "bed")
  }
  cell_map <- setNames(
    sub(paste0("^", s$sample, "_"), "", colnames(obj)),
    colnames(obj)
  )
  fragments <- CreateFragmentObject(
    path = frag,
    cells = cell_map,
    validate.fragments = truthy(args$validate_fragments)
  )
  atac <- FeatureMatrix(
    fragments = fragments,
    features = consensus,
    cells = colnames(obj),
    process_n = 2000
  )
  atac <- atac[, colnames(obj), drop = FALSE]   # FeatureMatrix may reorder

  # --- QC filter ---
  keep <- obj$nFeature_RNA >= args$min_features_rna &
    obj$nFeature_RNA <= args$max_features_rna &
    obj$nCount_RNA >= args$min_counts_rna &
    obj$percent.mt <= args$max_percent_mt &
    dbl_class == "singlet" &
    Matrix::colSums(atac) >= args$min_counts_atac
  keep[is.na(keep)] <- FALSE

  qc_rows[[s$sample]] <- data.frame(
    sample = s$sample,
    donor = s$donor,
    group = s$group,
    cells_in = ncol(obj),
    doublets = sum(dbl_class == "doublet", na.rm = TRUE),
    median_features_rna = stats::median(obj$nFeature_RNA),
    median_percent_mt = stats::median(obj$percent.mt),
    median_counts_atac = stats::median(Matrix::colSums(atac)),
    cells_out = sum(keep),
    stringsAsFactors = FALSE
  )
  msg("  ", ncol(obj), " cells -> ", sum(keep), " after QC")

  rna_list[[s$sample]] <- GetAssayData(obj, layer = "counts")[, keep, drop = FALSE]
  atac_list[[s$sample]] <- atac[, keep, drop = FALSE]
  obs_list[[s$sample]] <- data.frame(
    cell = colnames(obj)[keep],
    sample = s$sample,
    donor = s$donor,
    group = s$group,
    nCount_RNA = obj$nCount_RNA[keep],
    nFeature_RNA = obj$nFeature_RNA[keep],
    percent_mt = obj$percent.mt[keep],
    nCount_ATAC = Matrix::colSums(atac)[keep],
    nFeature_ATAC = Matrix::colSums(atac > 0)[keep],
    scDblFinder_score = dbl_score[keep],
    scDblFinder_class = dbl_class[keep],
    stringsAsFactors = FALSE
  )

  saveRDS(
    list(rna = rna_list[[s$sample]], atac = atac_list[[s$sample]],
         obs = obs_list[[s$sample]], qc = qc_rows[[s$sample]],
         peak_names = rownames(atac_list[[s$sample]])),
    ckpt, compress = FALSE
  )
}

gene_sets <- lapply(rna_list, rownames)
if (length(gene_sets) > 1 &&
    !all(vapply(gene_sets[-1], identical, logical(1), gene_sets[[1]]))) {
  stop("RNA gene sets differ between samples - they must share a reference")
}

rna <- do.call(cbind, unname(rna_list))
atac <- do.call(cbind, unname(atac_list))
obs <- do.call(rbind, unname(obs_list))

if (!identical(colnames(rna), colnames(atac))) stop("RNA and ATAC cell order diverged")
if (!identical(colnames(rna), obs$cell)) stop("obs.csv order does not match counts")

# ------------------------------------------------------------------- export --
write_mtx <- function(mat, path) {
  Matrix::writeMM(mat, path)
  if (truthy(args$gzip) && requireNamespace("R.utils", quietly = TRUE)) {
    R.utils::gzip(path, overwrite = TRUE, remove = TRUE)
    paste0(path, ".gz")
  } else {
    path
  }
}

rna_path <- write_mtx(rna, file.path(args$out, "rna_counts.mtx"))
atac_path <- write_mtx(atac, file.path(args$out, "atac_counts.mtx"))

writeLines(rownames(rna), file.path(args$out, "rna_features.tsv"))
writeLines(rownames(atac), file.path(args$out, "atac_features.tsv"))
utils::write.csv(obs, file.path(args$out, "obs.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, unname(qc_rows)),
                 file.path(args$out, "qc_per_sample.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(args$out, "sessionInfo.txt"))

msg("wrote ", rna_path, " (", nrow(rna), " genes x ", ncol(rna), " cells)")
msg("wrote ", atac_path, " (", nrow(atac), " peaks x ", ncol(atac), " cells)")
msg("group counts: ", paste(names(table(obs$group)), table(obs$group), sep = "=", collapse = " "))
