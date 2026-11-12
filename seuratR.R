
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)
library(GOstats)
library(org.Hs.eg.db)
library(AnnotationDbi)


h5_files <- list.files(
  pattern = "\\.h5",
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

write.csv(
  metadata, "sample_metadata.csv",
  row.names = FALSE
)


rna.list <- list()

for (i in seq_along(h5_files)) {
  
  sample_name <- metadata$sample[i]
    counts <- Read10X_h5(h5_files[i])
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      stop("Gene Expression matrix not found in ", sample_name)
    }
    
  }
    obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 0,
    min.features = 0
  )
  
  obj <- RenameCells(
    obj,
    add.cell.id = sample_name
  )
  
  obj$sample <- sample_name
  obj$donor <- metadata$donor[i]
  obj$group <- metadata$group[i]
  
  rna.list[[sample_name]] <- obj
  
}

length(rna.list)


for (i in seq_along(rna.list)) {
  
  obj <- rna.list[[i]]
  
  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )
  
  rna.list[[i]] <- obj
  
}


for (i in seq_along(rna.list)) {
  
  print(
    VlnPlot(
      rna.list[[i]],
      features = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt"
      ),
      ncol = 3,
      pt.size = 0
    )
  )
  
}

for (i in seq_along(rna.list)) {
  

}



for (i in seq_along(rna.list)) {
  
  obj <- rna.list[[i]]
  
  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  
  # Identify the top 2,000 highly variable genes
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    loess.span = 0.3,
    clip.max = "auto"
  )
  
  # Save the processed object
  rna.list[[i]] <- obj
  
}


for (i in seq_along(rna.list)) {
  
  rna.list[[i]] <- FindVariableFeatures(
    rna.list[[i]],
    selection.method = "vst",
    nfeatures = 2000,
    loess.span = 0.3,
    clip.max = "auto",
    verbose = FALSE
  )
  
}

anchors <- FindIntegrationAnchors(
  object.list = rna.list,
  anchor.features = 2000,
  reduction = "cca",
  dims = 1:20
)

rna.integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:20
)










DefaultAssay(rna.integrated) <- "integrated"

rna.integrated <- ScaleData(rna.integrated)

rna.integrated <- RunPCA(
  rna.integrated,
  npcs = 30
)

rna.integrated <- RunUMAP(
  rna.integrated,
  dims = 1:20
)

DimPlot(
  rna.integrated,
  reduction = "umap",
  group.by = "seurat_clusters"
)


rna.integrated <- FindNeighbors(
  rna.integrated,
  reduction = "pca",
  dims = 1:30
)

rna.integrated <- FindClusters(
  rna.integrated,graph.name = "integrated_nn",
  resolution = 0.5)

table(Idents(rna.integrated))


FeaturePlot(
  rna.integrated,
  features = c(
    "STMN2",
    "MBP",
    "OLIG1",
    "OLIG2",
    "TREM2",
    "GPR34"
  ),
  ncol = 3
)



cell_markers <- list(
  
  Neurons = c("STMN2", "TBR1"),
  
  Oligodendrocytes = c("MBP", "MOG"),
  
  Astrocytes = c("GFAP", "AQP4"),
  
  OPCs = c("OLIG1", "OLIG2"),
  
  Microglia = c("GPR34", "TREM2", "C1QC"),
  
  Endothelial = c("FLT1", "VWF", "PDGFRB"),
  
  T_cells = c("CD3E", "CD8A")
  
)

DefaultAssay(rna.integrated) <- "RNA"

features <- unique(unlist(cell_markers))

features <- intersect(
  features,
  rownames(rna.integrated)
)

DotPlot(
  rna.integrated,
  features = features,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

new_names <- c(
  "0" = "Oligodendrocytes",
  "1" = "Oligodendrocytes1",
  "2" = "Oligodendrocytes2",
  "3" = "Oligodendrocytes3",  
  "4" = "Astrocytes",
  "5" = "UN",
  "6" = "Microglia",
  "7" = "OPCs"
)

rna.integrated
rna.integrated <- RenameIdents(rna.integrated, new_names)

DimPlot(
  rna.integrated,
  reduction = "umap", lab
)

rna.integrated$celltype <- Idents(rna.integrated)

DimPlot(
  rna.integrated,
  reduction = "umap",
  #group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 0.5,
  cols = c(
    "Oligodendrocytes"  = "#3B75AF",
    "Oligodendrocytes1" = "#74ADD1",
    "Oligodendrocytes2" = "#ABD9E9",
    "Oligodendrocytes3" = "#4575B4",
    "Astrocytes"        = "#E69F00",
    "UN"                = "#999999",
    "Microglia"         = "#009E73",
    "OPCs"              = "#CC79A7"
  )
) +
  ggtitle("Cell Type Distribution") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "right"
  )



#############################################################################


bed_files <- list.files(
  pattern = "\\.bed.gz$",
)

peak_list <- lapply(
  bed_files,
  rtracklayer::import
)


combined.peaks <- reduce(
  do.call(c, peak_list)
)

combined.peaks <- combined.peaks[
  width(combined.peaks) >= 20 &
    width(combined.peaks) <= 10000
]

length(combined.peaks)

summary(width(combined.peaks))

annotations <- GetGRangesFromEnsDb(
  ensdb = EnsDb.Hsapiens.v86
)

seqlevelsStyle(annotations) <- "UCSC"

genome(annotations) <- "hg38"


fragment_files <- list.files(
  pattern = "\\.tsv.gz$",
)

for (f in fragment_files) {
  
  if (!file.exists(paste0(f, ".tbi"))) {
    
    Rsamtools::indexTabix(
      f,
      format = "bed"
    )
    
  }
  
}


atac.list <- list()

for (sample_name in names(rna.list)) {
  
  cat("Processing:", sample_name, "\n")
  
  # Find corresponding sample information
  i <- match(sample_name, metadata$sample)
  
  if (is.na(i)) next
  
  frag_path <- file.path(
    paste0(sample_name, "_atac_fragments.tsv.gz")
  )
  
  if (!file.exists(frag_path)) {
    warning("Fragment file missing: ", sample_name)
    next
  }
  
  # Get cells retained by RNA QC
  rna_obj <- rna.list[[sample_name]]
  
  new_cells <- colnames(rna_obj)
  
  # Recover original 10x barcodes
  raw_barcodes <- sub(
    paste0("^", sample_name, "_"),
    "",
    new_cells
  )
  
  # Map fragment-file barcodes to Seurat cell names
  cell_map <- raw_barcodes
  names(cell_map) <- new_cells
  
  # Create fragment object
  fragments <- CreateFragmentObject(
    path = frag_path,
    cells = cell_map,
    validate.fragments = TRUE
  )
  
  # Count fragments in common peaks
  peak_counts <- FeatureMatrix(
    fragments = fragments,
    features = combined.peaks,
    cells = new_cells,
    process_n = 2000
  )
  
  # Convert genomic coordinates to GRanges
  peak_ranges <- StringToGRanges(
    rownames(peak_counts),
    sep = c("-", "-")
  )
  
  # Create ChromatinAssay
  chrom_assay <- CreateChromatinAssay(
    counts = peak_counts,
    ranges = peak_ranges,
    fragments = fragments,
    genome = "hg38",
    annotation = annotations
  )
  
  # Create ATAC Seurat object
  atac_obj <- CreateSeuratObject(
    counts = peak_counts,
    assay = "ATAC",
    project = sample_name
  )
  
  atac_obj[["ATAC"]] <- chrom_assay
  
  # Transfer metadata from RNA object
  atac_obj <- AddMetaData(
    atac_obj,
    metadata = rna_obj@meta.data[
      colnames(atac_obj),
      ,
      drop = FALSE
    ]
  )
  
  atac.list[[sample_name]] <- atac_obj
  
}

length(atac.list)

atac.list[[1]]


for (i in seq_along(atac.list)) {
  
  obj <- atac.list[[i]]
  
  # Nucleosome signal
  obj <- NucleosomeSignal(obj)
  
  # TSS enrichment
  obj <- TSSEnrichment(
    obj,
    fast = FALSE
  )
  
  # Blacklist fraction
  obj$blacklist_fraction <- FractionCountsInRegion(
    obj,
    assay = "ATAC",
    regions = Signac::blacklist_hg38
  )
  
  atac.list[[i]] <- obj
  
}


atac.qc <- VlnPlot(
  atac.list[[2]],
  features = c(
    "nCount_ATAC",
    "nFeature_ATAC",
    "TSS.enrichment",
    "nucleosome_signal",
    "blacklist_fraction"
  ),
  ncol = 3,
  pt.size = 0
)

atac.qc


atac.integrated <- merge(
  x = atac.list[[1]],
  y = atac.list[-1]
)

atac.integrated

table(atac.integrated$group)

DefaultAssay(atac.integrated) <- "ATAC"

atac.integrated <- RunTFIDF(
  atac.integrated
)

atac.integrated <- FindTopFeatures(
  atac.integrated,
  min.cutoff = "q0"
)

atac.integrated <- RunSVD(
  atac.integrated,
  n = 30
)

DepthCor(atac.integrated)

atac.integrated <- harmony::RunHarmony(
  object = atac.integrated,
  group.by.vars = "sample",
  reduction.use = "lsi",
  dims.use = 2:30,
  project.dim = FALSE
)

atac.integrated <- FindNeighbors(
  object = atac.integrated,
  reduction = "harmony",
  dims = 2:30
)


atac.integrated <- RunUMAP(
  object = atac.integrated,
  reduction = "harmony",
  dims = 1:20,
  n.components = 2
)

DimPlot(
  atac.integrated,
  reduction = "umap",
  group.by = "group"
)

DimPlot(
  c(atac.integrated),
  reduction = "umap",
  group.by = "sample"
)

atac.integrated <- FindNeighbors(
  atac.integrated,
  reduction = "harmony",
  dims = 1:20
)

atac.integrated <- FindClusters(
  atac.integrated,
  resolution = 0.1
)

DimPlot(
  atac.integrated,
  reduction = "umap",
  group.by = "seurat_clusters"
)
########################################################################################################

p1 <- DimPlot(rna.integrated, label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(atac.integrated, group.by = "sample", label = FALSE) + NoLegend() + ggtitle("ATAC")
p1 + p2


gene.activities <- GeneActivity(atac.integrated, features = VariableFeatures(rna.integrated))

# add gene activities as a new assay
atac.integrated[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(atac.integrated) <- "ACTIVITY"
atac.integrated <- NormalizeData(atac.integrated)
atac.integrated <- ScaleData(atac.integrated, features = rownames(atac.integrated))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = rna.integrated, 
                                        query = atac.integrated, 
                                        features = VariableFeatures(object = rna.integrated),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna.integrated$celltype,
                                     weight.reduction = atac.integrated[["lsi"]], dims = 2:30)

atac.integrated <- AddMetaData(atac.integrated, metadata = celltype.predictions)


p1 <- DimPlot(atac.integrated, 
              group.by = "predicted.id", 
              label = TRUE) 

p2 <- DimPlot(rna.integrated, 
              group.by = "celltype", 
              label = TRUE) 
p1+p2




# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(rna.integrated)



rna.integrated[["RNA"]] <- JoinLayers(
  rna.integrated[["RNA"]]
)


refdata <- GetAssayData(rna.integrated, assay = "RNA", layer = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac.integrated[["lsi"]],
                           dims = 2:30)

atac.integrated[["RNA"]] <- imputation

coembed <- merge(x = rna.integrated, y = atac.integrated)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

DimPlot(coembed, group.by = c("orig.ident", "predicted.id"))


########################################################################################################


obj
