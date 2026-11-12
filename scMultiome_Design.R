###This script is to analyze scMULTIOME DATA in a normal way, then do simulations#########


library(Seurat)
library(Signac)
library(GenomicRanges)
library(Rsamtools)
library(EnsDb.Hsapiens.v86)

files <- list.files(
  pattern = "\\.h5$",
  full.names = TRUE
)


counts_list <- lapply(files, Read10X_h5)

seurat_list <- list()

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

seqlevelsStyle(annotations) <- 'UCSC'


for (i in seq_along(counts_list)) {
  rna_counts <- counts_list[[i]]$`Gene Expression`
  atac_counts <- counts_list[[i]]$Peaks
  
  # Initialize object with RNA
  obj <- CreateSeuratObject(
    counts = rna_counts, 
    assay = "RNA", 
    project = paste0("Sample_", i)
  )
  
  # Create and add the Chromatin Assay (ATAC)
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    annotation = annotations
  )
  obj[["ATAC"]] <- chrom_assay
  
  # RNA QC Metric: Mitochondrial percentage
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  seurat_list[[i]] <- obj
}

seurat_list


if (length(seurat_list) > 1) {
  seurat_obj <- merge(
    x = seurat_list[[1]], 
    y = seurat_list[-1],
      add.cell.ids = names(seurat_list),
      project = "PD-Samples"
    )
} else {
  seurat_obj <- seurat_list[[1]]
}




rm(seurat_list, counts_list) # Clean up memory






seurat_list <- lapply(seq_along(files), function(i) {
  
  counts <- Read10X_h5(files[i])
  
  # RNA assay
  obj <- CreateSeuratObject(
    counts = counts[["Gene Expression"]],
    min.cells = 3,
    min.features = 200,
    project = basename(files[i])
  )
  
  # ATAC assay
  peaks <- CreateChromatinAssay(
    counts = counts[["Peaks"]],
    sep = c(":", "-"),
    genome = "hg38"
  )
  
  obj[["ATAC"]] <- peaks
  
  obj
})

names(seurat_list) <- basename(files)



names(seurat_list) <- tools::file_path_sans_ext(basename(files))


seurat_obj <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "MiProyectoSingleCell"
)


View(seurat_obj@meta.data)


View(seurat_list$GSM5818649_C16007_filtered_feature_bc_matrix@neighbors)


#######