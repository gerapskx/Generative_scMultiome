library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)
library(GOstats)
library(AnnotationDbi)
library(GenomeInfoDb)
library(SingleCellExperiment)
library(scDesign3)
library(foreach)
library(doParallel)
library(future.apply)
library(harmony)

cores <- detectCores()
cl <- makeCluster(cores-2)
plan(multisession, workers=22)

h5_files <- list.files(
  pattern = "\\.h5$",
  full.names = TRUE
)

fragment_files <- list.files(
  pattern = "fragments\\.tsv\\.gz$",
  full.names = TRUE
)

for (f in fragment_files) {
  
  if (!file.exists(paste0(f, ".tbi"))) {
    
    Rsamtools::indexTabix(
      f,
      format = "bed"
    )
    
  }
  
}

annotations <- GetGRangesFromEnsDb(
  ensdb = EnsDb.Hsapiens.v86
)

seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"



sample.list <- list()

for (i in seq_along(h5_files)) {
  
  sample <- sub(
    "_filtered_feature_bc_matrix\\.h5$",
    "",
    basename(h5_files[i])
  )
  
  frag_file <- fragment_files[
    grepl(sample, basename(fragment_files), fixed = TRUE)
  ]
  
  if (length(frag_file) != 1) {
    stop(paste(
      "Could not uniquely identify fragment file for",
      sample,
      "\nMatching files:",
      paste(basename(frag_file), collapse = ", ")
    ))
  }
  
  message("Processing ", sample)

  inputdata.10x <- Read10X_h5(h5_files[i])
  
  rna_counts <- inputdata.10x[["Gene Expression"]]
  atac_counts <- inputdata.10x[["Peaks"]]
  
  samples <- CreateSeuratObject(
    counts = rna_counts,
    project = sample
  )
    samples[["percent.mt"]] <- PercentageFeatureSet(
    samples,
    pattern = "^MT-"
  )
  
  grange.counts <- StringToGRanges(
    rownames(atac_counts),
    sep = c(":", "-")
  )
  
  grange.use <- seqnames(grange.counts) %in%
    standardChromosomes(grange.counts)
  
  atac_counts <- atac_counts[as.vector(grange.use), ]
  

  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = frag_file,
    min.cells = 10,
    annotation = annotations
  )
  
  samples[["ATAC"]] <- chrom_assay
  samples$sample <- sample
  sample.list[[sample]] <- samples
}


set.seed(123) 

sample.list.subset <- sample.list[
  sample(seq_along(sample.list), 5)
]


MergedSamples <- merge(
  x = sample.list.subset[[1]],
  y = sample.list.subset[-1],
  add.cell.ids = names(sample.list.subset),
  project = "Multiome"
)


sample_id <- sub(
  "_filtered_feature_bc_matrix\\.h5$",
  "",
  basename(h5_files)
)

donor_id <- sub(
  "^GSM[0-9]+_",
  "",
  sample_id
)

GSM <- sub(
  "_.*$",
  "",
  sample_id
)

# Cohort assignment from donor naming convention
group <- ifelse(
  grepl("^Y", donor_id), "Young",
  ifelse(
    grepl("^C", donor_id), "Aged",
    ifelse(
      grepl("^P", donor_id), "PD",
      NA
    )
  )
)

metadata <- data.frame(
  GSM = GSM,
  sample = sample_id,
  donor = donor_id,
  group = group,
  h5 = h5_files
)

metadata$group <- factor(
  metadata$group,
  levels = c("Young", "Aged", "PD")
)

print(metadata)


MergedSamples$group



# RNA dimensions
dim(MergedSamples[["RNA"]])

# ATAC dimensions
dim(MergedSamples[["ATAC"]])


VlnPlot(MergedSamples, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()


cell_metadata <- metadata[
  match(MergedSamples$sample, metadata$sample),
]


Merged <- AddMetaData(
  object = MergedSamples,
  metadata = cell_metadata[, c("GSM", "donor", "group", "h5")]
)


  
###########################################################################

DefaultAssay(MergedSamples) <- "RNA"

MergedSamples <- SCTransform(
  MergedSamples,
  verbose = FALSE
)

MergedSamples <- RunPCA(
  MergedSamples,
  npcs = 50,
  verbose = FALSE
)

MergedSamples <- harmony::RunHarmony(
  object = MergedSamples,
  group.by.vars = "sample",
  reduction.use = "pca",
  dims.use = 1:50,
  reduction.save = "harmony.rna",
  verbose = FALSE
)

MergedSamples <- RunUMAP(
  MergedSamples,
  reduction = "harmony.rna",
  dims = 1:50,
  reduction.name = "umap.rna",
  reduction.key = "rnaUMAP_"
)

####################################################################################

DefaultAssay(MergedSamples) <- "ATAC"
MergedSamples
MergedSamples <- RunTFIDF(MergedSamples)

MergedSamples <- FindTopFeatures(MergedSamples, min.cutoff = 'q0')

MergedSamples <- RunSVD(MergedSamples)


MergedSamples <- harmony::RunHarmony(
  object = MergedSamples,
  group.by.vars = "sample",
  reduction.use = "lsi",
  dims.use = 2:50,
  reduction.save = "harmony.atac",
  project.dim = FALSE,
  verbose = FALSE
)


MergedSamples <- RunUMAP(
  MergedSamples,
  reduction = "harmony.atac",
  dims = 1:49,
  reduction.name = "umap.atac",
  reduction.key = "atacUMAP_"
)


MergedSamples

MergedSamples <- FindMultiModalNeighbors(MergedSamples, reduction.list = list("harmony.rna", "harmony.atac"), dims.list = list(1:50, 1:49))
MergedSamples <- RunUMAP(MergedSamples, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
MergedSamples <- FindClusters(MergedSamples, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = FALSE)


p1 <- DimPlot(MergedSamples, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(MergedSamples, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(MergedSamples, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 +p3 + theme(plot.title = element_text(hjust = 0.5))

DefaultAssay(MergedSamples) <- "RNA"

MergedSamples <- FindClusters(
  MergedSamples,graph.name = "wsnn",
  resolution = 0.1)


cell_markers <- list(
  
  Neurons = c("STMN2", "TBR1"),
  
  Oligodendrocytes = c("MBP", "MOG"),
  
  Astrocytes = c("GFAP", "AQP4"),
  
  OPCs = c("OLIG1", "OLIG2"),
  
  Microglia = c("GPR34", "TREM2", "C1QC"),
  
  Endothelial = c("FLT1", "VWF", "PDGFRB"),
  
  T_cells = c("CD3E", "CD8A")
  
)


features <- unique(unlist(cell_markers))


features <- intersect(
  features,
  rownames(MergedSamples)
)


DotPlot(MergedSamples,
  features = features, assay = "RNA",
  group.by = "seurat_clusters"
) +
  RotatedAxis()

MergedSamples <- FindClusters(
  MergedSamples,graph.name = "wsnn",
  resolution = 0.1)

Graphs(MergedSamples)

DimPlot(
  MergedSamples,
  reduction = "wnn.umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

new_names <- c(
  "0" = "Oligodendrocytes0",
  "1" = "Oligodendrocytes1",
  "2" = "Oligodendrocytes2",
  "3" = "Microglia",
  "4" = "Astrocytes",
  "5" = "OPCs1",
  "6" = "OPCs2",
  "7" = "Neurons"
)


MergedSamples$celltype <- Idents(MergedSamples)

MergedSamples <- RenameIdents(MergedSamples, new_names)
                               
p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))


CoveragePlot(MergedSamples, region = 'GFAP', features = 'GFAP', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)


saveRDS(
  MergedSamples,
  file = "MergedSamples.rds"
)
#####################################

MergedSamples <- JoinLayers(
  MergedSamples,
  assay = "RNA"
)

# Extract raw RNA counts
rna_counts <- GetAssayData(
  MergedSamples,
  assay = "RNA",
  layer = "counts"
)

# Extract raw ATAC counts
atac_counts <- GetAssayData(
  MergedSamples,
  assay = "ATAC",
  layer = "counts"
)

# Confirm dimensions
dim(rna_counts)
dim(atac_counts)

setequal(
  colnames(rna_counts),
  colnames(atac_counts)
)

atac_counts <- atac_counts[
  ,
  colnames(atac_counts),
  drop = FALSE
]

# Confirm that cells are identical and in the same order
stopifnot(
  identical(colnames(rna_counts), colnames(rna_counts))
)

cell_metadata <- MergedSamples@meta.data
cell_metadata
# Align metadata to count matrices
cell_metadata <- cell_metadata[
  colnames(rna_counts),
  ,
  drop = FALSE
]

# Check alignment
stopifnot(
  identical(rownames(cell_metadata), colnames(rna_counts))
)

# Keep relevant variables
cell_metadata <- cell_metadata[, c(
  "celltype",
  "sample",
  "donor",
  "group"
), drop = FALSE]

# Remove cells with missing cell-type annotation
valid_cells <- !is.na(cell_metadata$celltype)

cell_metadata <- cell_metadata[valid_cells, , drop = FALSE]

rna_counts <- rna_counts[, valid_cells, drop = FALSE]
atac_counts <- atac_counts[, valid_cells, drop = FALSE]

# Cell type as a factor
cell_metadata$celltype <- factor(cell_metadata$celltype)

table(cell_metadata$celltype)
table(cell_metadata$group)





select_variable_features <- function(
    counts,
    n_features = 100,
    min_fraction = 0.05
) {
  
  counts <- as(counts, "dgCMatrix")
  
  # Minimum number of cells expressing a feature
  min_cells <- ceiling(min_fraction * ncol(counts))
  
  # Remove features detected in too few cells
  detected_cells <- Matrix::rowSums(counts > 0)
  
  keep <- detected_cells >= min_cells
  
  counts <- counts[keep, , drop = FALSE]
  
  # Calculate variance without converting the entire
  # sparse matrix into a dense matrix
  feature_means <- Matrix::rowMeans(counts)
  
  feature_variance <- Matrix::rowMeans(counts^2) -
    feature_means^2
  
  # Remove zero-variance and invalid features
  valid <- is.finite(feature_variance) &
    feature_variance > 0
  
  counts <- counts[valid, , drop = FALSE]
  feature_variance <- feature_variance[valid]
  
  # Select top variable features
  top_features <- order(
    feature_variance,
    decreasing = TRUE
  )
  
  top_features <- head(top_features, n_features)
  
  counts <- counts[top_features, , drop = FALSE]
  
  return(counts)
}



# RNA
rna_counts_filtered <- select_variable_features(
  rna_counts,
  n_features = 100,
  min_fraction = 0.05
)

# ATAC
atac_counts_filtered <- select_variable_features(
  atac_counts,
  n_features = 100,
  min_fraction = 0.05
)

# Dimensions
dim(rna_counts_filtered)
dim(atac_counts_filtered)





#########################################3

# RNA SCE
sce_rna <- SingleCellExperiment(
  assays = list(
    counts = rna_counts_filtered
  ),
  colData = cell_metadata
)

# ATAC SCE
sce_atac <- SingleCellExperiment(
  assays = list(
    counts = atac_counts_filtered
  ),
  colData = cell_metadata
)

# Verify
sce_rna
sce_atac

stopifnot(
  identical(colnames(sce_rna), colnames(sce_atac))
)

stopifnot(
  identical(
    as.character(sce_rna$celltype),
    as.character(sce_atac$celltype)
  )
)



set.seed(123)

library(BiocParallel)

bp_param <- SnowParam(
  workers = 20,
  progressbar = TRUE
)

set.seed(123)

sim_rna <- scdesign3(
  sce = sce_rna,
  assay_use = "counts",
  celltype = "celltype",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  mu_formula = "celltype",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 20,
  usebam = FALSE,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = TRUE,
  parallelization = "pbmapply",
  BPPARAM = bp_param
)



# Check dimensions
dim(sim_rna$new_count)


logcounts(sce_rna) <- log1p(counts(sce_rna))
simu_RNA <- SingleCellExperiment(list(counts = sim_rna$new_count), colData = sim_rna$new_covariate)
logcounts(simu_RNA) <- log1p(counts(simu_RNA))

set.seed(123)


compare_figure <- plot_reduceddim(ref_sce = sce_rna, 
                                  sce_list = list(simu_RNA), 
                                  name_vec = c("Reference", "scDesign3"),
                                  assay_use = "logcounts", 
                                  if_plot = TRUE, 
                                  color_by = "celltype", 
                                  n_pc = 20)

plot(compare_figure$p_umap)


set.seed(456)

library(BiocParallel)

sim_atac <- scdesign3(
  sce = sce_atac,
  assay_use = "counts",
  celltype = "celltype",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  mu_formula = "celltype",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 12,
  usebam = FALSE,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = TRUE, 
  parallelization = "pbmapply",
  BPPARAM = bp_param
)

# Check dimensions
dim(sim_atac$new_count)


logcounts(sce_atac) <- log1p(counts(sce_atac))
simu_ATAC <- SingleCellExperiment(list(counts = sim_atac$new_count), colData = sim_atac$new_covariate)
logcounts(simu_ATAC) <- log1p(counts(simu_ATAC))

compare_figureATAC <- plot_reduceddim(ref_sce = sce_atac, 
                                  sce_list = list(simu_ATAC), 
                                  name_vec = c("Reference", "scDesign3"),
                                  assay_use = "logcounts", 
                                  if_plot = TRUE, 
                                  color_by = "celltype", 
                                  n_pc = 20)

p1 <- plot(compare_figure$p_umap)
p2 <- plot(compare_figureATAC$p_umap)

p1+p2


######################################################################

# Extract synthetic count matrices
sim_rna_counts <- sim_rna$new_count
sim_atac_counts <- sim_atac$new_count

# Generate synthetic cell names
synthetic_cells <- paste0(
  "Synthetic_Cell_",
  seq_len(ncol(sim_rna_counts))
)

# Verify that both simulations generated the same number of cells
stopifnot(
  ncol(sim_rna_counts) == ncol(sim_atac_counts)
)

stopifnot(
  nrow(sim_rna_counts) == nrow(rna_counts_filtered)
)

stopifnot(
  nrow(sim_atac_counts) == nrow(atac_counts_filtered)
)

# Assign feature and cell names
dimnames(sim_rna_counts) <- list(
  rownames(rna_counts_filtered),
  synthetic_cells
)

dimnames(sim_atac_counts) <- list(
  rownames(atac_counts_filtered),
  synthetic_cells
)

# Construct synthetic metadata
sim_metadata <- cell_metadata

rownames(sim_metadata) <- synthetic_cells

# Confirm metadata alignment
stopifnot(
  identical(rownames(sim_metadata), colnames(sim_rna_counts))
)













# Create synthetic RNA Seurat object
SyntheticMultiome <- CreateSeuratObject(
  counts = sim_rna_counts,
  assay = "RNA",
  meta.data = sim_metadata,
  project = "Synthetic_Multiome"
)

# Create synthetic ATAC chromatin assay
synthetic_chrom_assay <- CreateChromatinAssay(
  counts = sim_atac_counts,
  sep = c(":", "-"),
  genome = "hg38",
  min.cells = 0,
  min.features = 0
)

# Add ATAC assay
SyntheticMultiome[["ATAC"]] <- synthetic_chrom_assay

# Inspect synthetic object
SyntheticMultiome

# Check dimensions
dim(SyntheticMultiome[["RNA"]])
dim(SyntheticMultiome[["ATAC"]])

# Check metadata
head(SyntheticMultiome@meta.data)

# Check cell-type composition
table(SyntheticMultiome$celltype)
table(MergedSamples$celltype)

DefaultAssay(SyntheticMultiome) <- "RNA"

SyntheticMultiome <- NormalizeData(
  SyntheticMultiome
)

SyntheticMultiome <- FindVariableFeatures(
  SyntheticMultiome
)

SyntheticMultiome <- ScaleData(
  SyntheticMultiome
)

SyntheticMultiome <- RunPCA(
  SyntheticMultiome,
  npcs = 50
)

SyntheticMultiome <- RunUMAP(
  SyntheticMultiome,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap.rna",
  reduction.key = "rnaUMAP_"
)

DefaultAssay(SyntheticMultiome) <- "ATAC"

SyntheticMultiome <- RunTFIDF(
  SyntheticMultiome
)

SyntheticMultiome <- FindTopFeatures(
  SyntheticMultiome,
  min.cutoff = "q0"
)

SyntheticMultiome <- RunSVD(
  SyntheticMultiome,
  n = 50
)

SyntheticMultiome <- RunUMAP(
  SyntheticMultiome,
  reduction = "lsi",
  dims = 2:30,
  reduction.name = "umap.atac",
  reduction.key = "atacUMAP_"
)

SyntheticMultiome <- FindMultiModalNeighbors(
  SyntheticMultiome,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:30)
)

SyntheticMultiome <- RunUMAP(
  SyntheticMultiome,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

SyntheticMultiome <- FindClusters(
  SyntheticMultiome,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.2
)


VlnPlot(
  MergedSamples,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "celltype"
)

VlnPlot(
  SyntheticMultiome,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "celltype"
)


UMAP1 <- DimPlot(MergedSamples, 
                 reduction = "wnn.umap", 
                 group.by = "celltype")
UMAP2 <- DimPlot(SyntheticMultiome, 
                 reduction = "wnn.umap", 
                 group.by = "celltype")
UMAP1 + UMAP2


MergedSamples

CoveragePlot(SyntheticMultiome, region = 'STMN2', features = 'STMN2', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

