#!/usr/bin/env Rscript
# Stage 2 - cluster + annotate (gate)
#
# Builds one Seurat object holding both modalities, then clusters RNA (PCA), ATAC
# (LSI) and the joint weighted-neighbour graph. Writes the cluster assignment plus
# the marker table a human needs to fill in config/celltype_map.csv.
#
# With one donor per condition, donor and condition are perfectly confounded and an
# un-integrated embedding splits cells by donor. `--integrate harmony` (the default)
# harmonises each modality over donor so clusters are cell-type driven; the
# embedding is only used for clustering, all statistics use raw counts.
# Use `--integrate none` to inspect the raw, donor-confounded clusters.
#
# Clustering is unsupervised and runs on ALL cells (fit + holdout); the split and
# every feature-selection step happen afterwards, so no label information leaks.
#
# Run from the project root:
#   Rscript R/01_cluster.R
#
# Outputs (results/cluster/):
#   clusters.csv         cell, rna_cluster, atac_cluster, cluster, umap_1, umap_2
#   cluster_markers.csv  top RNA markers per joint cluster
#   cluster_sizes.csv    cells per joint cluster, split by donor/group
#   sessionInfo.txt
# and a template config/celltype_map.csv (only created if absent).

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(Matrix)
})

defaults <- list(
  ingest = "results/ingest",
  out = "results/cluster",
  map = "config/celltype_map.csv",
  resolution = 0.8,
  dims_rna = 30,
  dims_atac = 30,
  n_hvg = 3000,
  integrate = "harmony",
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
for (k in c("resolution", "dims_rna", "dims_atac", "n_hvg", "seed")) {
  args[[k]] <- as.numeric(args[[k]])
}
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
msg("reading stage-0 export from ", args$ingest)
obs <- read.csv(file.path(args$ingest, "obs.csv"), stringsAsFactors = FALSE,
                check.names = FALSE)
rownames(obs) <- obs$cell
rna <- read_mtx(file.path(args$ingest, "rna_counts.mtx"))
atac <- read_mtx(file.path(args$ingest, "atac_counts.mtx"))
genes <- readLines(file.path(args$ingest, "rna_features.tsv"))
peaks <- readLines(file.path(args$ingest, "atac_features.tsv"))

if (nrow(rna) != length(genes)) stop("rna matrix / feature list mismatch")
if (nrow(atac) != length(peaks)) stop("atac matrix / feature list mismatch")
if (ncol(rna) != nrow(obs) || ncol(atac) != nrow(obs)) stop("cell count mismatch")
rownames(rna) <- genes
colnames(rna) <- obs$cell
rownames(atac) <- peaks
colnames(atac) <- obs$cell
msg(nrow(rna), " genes x ", ncol(rna), " cells; ", nrow(atac), " peaks")

# --------------------------------------------------------------------- object --
set.seed(args$seed)
obj <- CreateSeuratObject(counts = rna, assay = "RNA", project = "gse193688",
                          min.cells = 0, min.features = 0)
# rownames from FeatureMatrix are chr1-100-200
gr <- Signac::StringToGRanges(peaks, sep = c("-", "-"))
obj[["ATAC"]] <- CreateChromatinAssay(counts = atac, ranges = gr, genome = "hg38")
obj$donor <- obs$donor
obj$group <- factor(obs$group, levels = c("Aged", "PD", "Young"))
obj$sample <- obs$sample

# ---------------------------------------------------------------------- RNA --
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = args$n_hvg, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = args$dims_rna, verbose = FALSE)

# --------------------------------------------------------------------- ATAC --
DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj, verbose = FALSE)
obj <- FindTopFeatures(obj, min.cutoff = "q0", verbose = FALSE)
obj <- RunSVD(obj, n = args$dims_atac, verbose = FALSE)
# LSI 1 correlates with sequencing depth; drop it from every downstream use
atac_dims <- 2:args$dims_atac

# ----------------------------------------------------------- per-modality ----
msg("clustering RNA and ATAC separately")
if (args$integrate == "harmony") {
  # With one donor per condition, donor and condition are perfectly confounded,
  # so an un-integrated embedding splits cells by donor. Harmony by donor gives
  # cell-type clusters that both donors share. The embedding is used ONLY for
  # clustering and annotation; every downstream statistic uses raw counts, so the
  # condition effect is not removed from the analysis.
  msg("harmonising RNA PCA and ATAC LSI over donor")
  obj <- harmony::RunHarmony(obj, group.by.vars = "donor", reduction.use = "pca",
                             dims.use = 1:args$dims_rna, project.dim = FALSE,
                             reduction.save = "harmony_rna", verbose = FALSE)
  obj <- harmony::RunHarmony(obj, group.by.vars = "donor", reduction.use = "lsi",
                             dims.use = atac_dims, project.dim = FALSE,
                             reduction.save = "harmony_atac", verbose = FALSE)
  rna_reduction <- "harmony_rna"
  atac_reduction <- "harmony_atac"
  # harmony returns only the requested dims, renumbered from 1
  rna_dims_use <- seq_len(args$dims_rna)
  atac_dims_use <- seq_along(atac_dims)
} else if (args$integrate == "none") {
  rna_reduction <- "pca"
  atac_reduction <- "lsi"
  rna_dims_use <- seq_len(args$dims_rna)
  atac_dims_use <- atac_dims
} else {
  stop("unknown --integrate value: ", args$integrate, " (use harmony or none)")
}

obj <- FindNeighbors(obj, reduction = rna_reduction, dims = rna_dims_use,
                     graph.name = "rna_nn", verbose = FALSE)
obj <- FindClusters(obj, graph.name = "rna_nn", resolution = args$resolution,
                    verbose = FALSE)
obj$rna_cluster <- as.character(Idents(obj))

obj <- FindNeighbors(obj, reduction = atac_reduction, dims = atac_dims_use,
                     graph.name = "atac_nn", verbose = FALSE)
obj <- FindClusters(obj, graph.name = "atac_nn", resolution = args$resolution,
                    verbose = FALSE)
obj$atac_cluster <- as.character(Idents(obj))

# -------------------------------------------------------------------- joint --
msg("building the weighted-neighbour graph")
obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list(rna_reduction, atac_reduction),
  dims.list = list(rna_dims_use, atac_dims_use),
  verbose = FALSE
)
obj <- FindClusters(obj, graph.name = "wsnn", resolution = args$resolution,
                    verbose = FALSE)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
               reduction.key = "wnnUMAP_", verbose = FALSE)

joint_cluster <- as.character(Idents(obj))
# stable, sortable labels: c01, c02, ... (numeric order of the raw factor)
lv <- levels(Idents(obj))
relabel <- setNames(sprintf("c%02d", seq_along(lv)), lv)
joint_cluster <- unname(relabel[joint_cluster])

umap <- Embeddings(obj, "wnn.umap")
clusters <- data.frame(
  cell = colnames(obj),
  rna_cluster = obj$rna_cluster,
  atac_cluster = obj$atac_cluster,
  cluster = joint_cluster,
  umap_1 = umap[, 1],
  umap_2 = umap[, 2],
  donor = obs[colnames(obj), "donor"],
  group = obs[colnames(obj), "group"],
  stringsAsFactors = FALSE
)
utils::write.csv(clusters, file.path(args$out, "clusters.csv"), row.names = FALSE)

sizes <- as.data.frame.matrix(table(clusters$cluster, clusters$donor))
sizes$total <- rowSums(sizes)
sizes <- cbind(cluster = rownames(sizes), sizes)
utils::write.csv(sizes, file.path(args$out, "cluster_sizes.csv"), row.names = FALSE)

# ------------------------------------------------------------------ markers --
msg("finding cluster markers")
Idents(obj) <- factor(joint_cluster, levels = sort(unique(joint_cluster)))
DefaultAssay(obj) <- "RNA"
markers <- FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, logfc.threshold = 0.25,
                          min.pct = 0.1, verbose = FALSE)
markers <- markers[order(markers$cluster, markers$avg_log2FC, decreasing = TRUE), ]
utils::write.csv(markers, file.path(args$out, "cluster_markers.csv"), row.names = FALSE)

# ------------------------------------------------------------- annotation ----
# The gate: a human fills `celltype` in the map. Until every cluster has a label
# the downstream stages must refuse to run.
if (!file.exists(args$map)) {
  template <- data.frame(
    cluster = sort(unique(clusters$cluster)),
    celltype = "",
    note = "",
    stringsAsFactors = FALSE
  )
  utils::write.csv(template, args$map, row.names = FALSE)
  msg("wrote annotation template ", args$map, " - fill in `celltype` before continuing")
} else {
  msg("keeping existing ", args$map)
}

capture.output(sessionInfo(), file = file.path(args$out, "sessionInfo.txt"))
msg("wrote ", file.path(args$out, "clusters.csv"), " (", nrow(clusters), " cells, ",
    length(unique(joint_cluster)), " clusters)")
print(table(clusters$cluster))
