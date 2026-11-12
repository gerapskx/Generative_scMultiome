################################################################################

if (!require("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("SONGDONGYUAN1994/scDesign3")
install.packages('Seurat')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("Signac", quietly = TRUE)) {
  install.packages("Signac")
}

BiocManager::install("MultiAssayExperiment")
BiocManager::install("SingleCellExperiment")
BiocManager::install("SingleCellMultiModal")

library("scDesign3")
library("SingleCellExperiment")
library("MultiAssayExperiment")
library("SingleCellMultiModal")



#################################################################################


data_dir <- "GSE298008_11022026"

base_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE298nnn/GSE298008/suppl/"
)

files <- c(
  "GSE298008_ATAC_counts_matrix.mtx.gz",
  "GSE298008_ATAC_regions.tsv.gz",
  "GSE298008_RNA_counts_matrix.mtx.gz",
  "GSE298008_SCT_counts_matrix.mtx.gz",
  "GSE298008_barcodes.tsv.gz",
  "GSE298008_cell_metadata.csv.gz",
  "GSE298008_genes.tsv.gz"
)

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

for (f in files) {
  
  dest <- file.path(data_dir, f)
  
  if (file.exists(dest) && file.info(dest)$size > 0) {
    message("Already downloaded: ", f)
    next
  }
  
  message("Downloading: ", f)
  
  download.file(
    url = paste0(base_url, utils::URLencode(f, reserved = TRUE)),
    destfile = dest,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )
  
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    stop("Download failed: ", f)
  }
  
  message(
    "Completed: ", f,
    " (", round(file.info(dest)$size / 1024^2, 2), " MB)"
  )
}

message("Download process completed.")


#################################################################3
  


# ============================================================
# GSE298008: Create a Seurat object from RNA counts
# ============================================================

library(Seurat)


# 1. Read input files -----------------------------------------

rna_file <- file.path(
  data_dir, "GSE298008_RNA_counts_matrix.mtx.gz"
)

barcode_file <- file.path(
  data_dir, "GSE298008_barcodes.tsv.gz"
)

metadata_file <- file.path(
  data_dir, "GSE298008_cell_metadata.csv.gz"
)

genes_file <- file.path(
  data_dir, "GSE298008_genes.tsv.gz"
)

required_files <- c(
  rna_file, barcode_file, metadata_file, genes_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing files: ",
    paste(basename(missing_files), collapse = ", ")
  )
}

# 2. Load barcodes, genes and metadata ------------------------

barcodes <- read.delim(
  gzfile(barcode_file),
  header = FALSE,
  stringsAsFactors = FALSE
)[[1]]

genes <- read.delim(
  gzfile(genes_file),
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata <- read.csv(
  gzfile(metadata_file),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Inspect file structure
print(dim(genes))
print(head(genes))
print(names(metadata))
print(head(metadata))

stopifnot(!anyDuplicated(barcodes))

# 3. Read the sparse RNA count matrix -------------------------

con <- gzfile(rna_file, open = "rt")
rna_counts <- tryCatch(
  Matrix::readMM(con),
  finally = close(con)
)

rna_counts <- as(rna_counts, "CsparseMatrix")

message(
  "RNA matrix dimensions: ",
  nrow(rna_counts), " features x ",
  ncol(rna_counts), " cells"
)

# 4. Check feature and barcode dimensions ---------------------

stopifnot(
  nrow(rna_counts) == nrow(genes),
  ncol(rna_counts) == length(barcodes)
)

# Use the first column as feature IDs.
# If the file contains gene symbols in another column,
# inspect head(genes) and change this selection if appropriate.
gene_ids <- as.character(genes[[1]])

if (anyNA(gene_ids) || any(gene_ids == "")) {
  stop("Gene identifiers contain missing or empty values.")
}

# Ensure unique feature names for Seurat
gene_ids <- make.unique(gene_ids)

rownames(rna_counts) <- gene_ids
colnames(rna_counts) <- as.character(barcodes)

# 5. Identify the barcode column in metadata ------------------

# Find the metadata column with exact matches to the barcodes.
overlap <- vapply(
  metadata,
  function(x) {
    sum(as.character(x) %in% barcodes, na.rm = TRUE)
  },
  numeric(1)
)

print(sort(overlap, decreasing = TRUE))

barcode_col <- names(which.max(overlap))

if (length(barcode_col) == 0 ||
    max(overlap) != length(barcodes)) {
  stop(
    "Could not identify a metadata column matching every barcode. ",
    "Inspect names(metadata), head(metadata), and head(barcodes)."
  )
}

message("Barcode column identified: ", barcode_col)

# Reorder metadata to match the RNA matrix exactly
metadata_barcodes <- as.character(metadata[[barcode_col]])

if (anyDuplicated(metadata_barcodes)) {
  stop("The identified metadata barcode column has duplicates.")
}

metadata <- metadata[
  match(barcodes, metadata_barcodes),
  ,
  drop = FALSE
]

stopifnot(
  nrow(metadata) == length(barcodes),
  identical(
    as.character(metadata[[barcode_col]]),
    as.character(barcodes)
  )
)

rownames(metadata) <- as.character(barcodes)

# Keep the barcode column in metadata only if useful.
# It is usually redundant because cell names are already barcodes.
metadata[[barcode_col]] <- NULL

# 6. Create the Seurat object -------------------------------

seurat_obj <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  project = "GSE298008",
  meta.data = metadata,
  min.cells = 0,
  min.features = 0
)

# 7. Basic quality control -----------------------------------

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)

print(seurat_obj)
print(Assays(seurat_obj))
print(dim(seurat_obj))
head(seurat_obj[[]])

# 8. Save the object -----------------------------------------

saveRDS(
  seurat_obj,
  file = file.path(data_dir, "GSE298008_RNA_seurat.rds"),
  compress = FALSE
)

message("Seurat object saved successfully.")


#########################################################

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

seurat_obj
DimPlot(seurat_obj, group.by = "seurat_clusters")

library(Seurat)
library(SingleCellExperiment)

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 1000
)

# Extract the selected genes
top_genes <- VariableFeatures(seurat_obj)

# Subset the Seurat object
seurat_obj_short <- subset(
  seurat_obj,
  features = top_genes
)

# Convert to SingleCellExperiment
sce_filtered <- as.SingleCellExperiment(seurat_obj_short)

dim(sce_filtered)



#########################################################


mae <- scMultiome(DataType = "mouse_gastrulation",
  mod = "*",
  dry.run = FALSE,
  format = "MTX"
)

mae
experiments(mae)
mae

rna <- experiments(mae)[["rna"]]
rna
# 1. Check the dataset and its metadata
class(rna)


colnames(as.data.frame(SummarizedExperiment::colData(rna)))

rna <- rna[1:200, ]
rna

######################################################


counts <- assay(sce_filtered, "counts")

gene_stats <- data.frame(
  gene = rownames(counts),
  mean = Matrix::rowMeans(counts),
  detected_cells = Matrix::rowSums(counts > 0),
  variance = apply(as.matrix(counts), 1, var)
)

gene_stats <- gene_stats[
  order(gene_stats$variance, decreasing = TRUE),
]

head(gene_stats, 20)

#########################################################

example_simu <- scdesign3(
  sce = sce_filtered,
  assay_use = "counts",
  celltype = "seurat_clusters",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
 mu_formula = "seurat_clusters",    
  sigma_formula = "1",      
  family_use = "nb",
  n_cores = 12,
  corr_formula = "1",
  copula = "gaussian"
)


dim(example_simu$new_count)

logcounts(sce_filtered) <- log1p(counts(sce_filtered))
simu_sce <- SingleCellExperiment(list(counts = example_simu$new_count), colData = example_simu$new_covariate)
logcounts(simu_sce) <- log1p(counts(simu_sce))

set.seed(123)


compare_figure <- plot_reduceddim(ref_sce = sce_filtered, 
                                  sce_list = list(simu_sce), 
                                  name_vec = c("Reference", "scDesign3"),
                                  assay_use = "logcounts", 
                                  if_plot = TRUE, 
                                  color_by = "seurat_clusters", 
                                  n_pc = 20)

plot(compare_figure$p_umap)


library(SingleCellExperiment)
library(Matrix)

# Match genes shared by both datasets
common_genes <- intersect(
  rownames(sce_filtered),
  rownames(simu_sce)
)

ref_counts <- assay(sce_filtered, "counts")[common_genes, , drop = FALSE]
sim_counts <- assay(simu_sce, "counts")[common_genes, , drop = FALSE]

ref_log <- assay(sce_filtered, "logcounts")[common_genes, , drop = FALSE]
sim_log <- assay(simu_sce, "logcounts")[common_genes, , drop = FALSE]

# Calculate metrics for each gene
gene_metrics <- data.frame(
  gene = common_genes,
  ref_mean = Matrix::rowMeans(ref_log),
  sim_mean = Matrix::rowMeans(sim_log),
  ref_var = apply(ref_log, 1, var),
  sim_var = apply(sim_log, 1, var),
  ref_zero = Matrix::rowMeans(ref_counts == 0),
  sim_zero = Matrix::rowMeans(sim_counts == 0)
)

# Relative agreement of gene means
gene_metrics$mean_difference <-
  gene_metrics$sim_mean - gene_metrics$ref_mean

gene_metrics$abs_mean_difference <-
  abs(gene_metrics$mean_difference)

# Difference in zero proportions
gene_metrics$zero_difference <-
  gene_metrics$sim_zero - gene_metrics$ref_zero

# Correlation of gene-level summaries
cor(
  gene_metrics$ref_mean,
  gene_metrics$sim_mean,
  method = "spearman"
)

cor(
  gene_metrics$ref_var,
  gene_metrics$sim_var,
  method = "spearman"
)

# Summarize errors across genes
summary(gene_metrics$abs_mean_difference)
summary(abs(gene_metrics$zero_difference))

par(mfrow = c(1, 3))

plot(
  gene_metrics$ref_mean,
  gene_metrics$sim_mean,
  xlab = "Reference mean log-expression",
  ylab = "Simulated mean log-expression",
  main = "Gene means",
  pch = 16, cex = 0.5
)
abline(0, 1, col = "red")

plot(
  gene_metrics$ref_var,
  gene_metrics$sim_var,
  xlab = "Reference variance",
  ylab = "Simulated variance",
  main = "Gene variances",
  pch = 16, cex = 0.5
)
abline(0, 1, col = "red")

plot(
  gene_metrics$ref_zero,
  gene_metrics$sim_zero,
  xlab = "Reference zero fraction",
  ylab = "Simulated zero fraction",
  main = "Zero proportions",
  pch = 16, cex = 0.5
)
abline(0, 1, col = "red")

par(mfrow = c(1, 1))

ks_results <- lapply(common_genes, function(g) {
  x <- as.numeric(ref_log[g, ])
  y <- as.numeric(sim_log[g, ])
  
  test <- ks.test(x, y, exact = FALSE)
  
  data.frame(
    gene = g,
    ks_statistic = unname(test$statistic),
    p_value = test$p.value
  )
})

ks_results <- do.call(rbind, ks_results)

# Correct for testing many genes
ks_results$FDR <- p.adjust(
  ks_results$p_value,
  method = "BH"
)

head(
  ks_results[order(-ks_results$ks_statistic), ],
  20
)

#######################################################

library(caret)

set.seed(123)

# Use a manageable set of shared genes
# Ideally, choose highly variable genes using the reference data only
n_genes <- min(100, length(common_genes))

selected_genes <- common_genes[
  order(
    apply(ref_log, 1, var),
    decreasing = TRUE
  )[seq_len(n_genes)]
]

# Cells x genes
ref_x <- t(as.matrix(ref_log[selected_genes, , drop = FALSE]))
sim_x <- t(as.matrix(sim_log[selected_genes, , drop = FALSE]))

# Balance the number of cells
n <- min(nrow(ref_x), nrow(sim_x))

ref_x <- ref_x[sample(nrow(ref_x), n), , drop = FALSE]
sim_x <- sim_x[sample(nrow(sim_x), n), , drop = FALSE]

X <- rbind(ref_x, sim_x)

y <- factor(
  c(rep("Reference", n), rep("Simulated", n)),
  levels = c("Reference", "Simulated")
)

# Scale features
X <- scale(X)
X[!is.finite(X)] <- 0

library(caret)
# Cross-validated classification
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

classifier <- train(
  x = X,
  y = y,
  method = "glm",
  family = binomial(),
  metric = "ROC",
  trControl = ctrl
)

classifier

############################################

library(caret)
library(pROC)
library(ranger)

set.seed(123)

# X: cells x genes numeric matrix or data frame
# y: class labels for each cell

y <- factor(y, levels = c("Reference", "Simulated"))

stopifnot(nrow(X) == length(y))
stopifnot(all(table(y) > 0))

# Stratified 5-fold cross-validation
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Logistic regression
fit_glm <- train(
  x = X,
  y = y,
  method = "glm",
  family = binomial(),
  metric = "ROC",
  trControl = ctrl
)

# Random forest
fit_rf <- train(
  x = X,
  y = y,
  method = "ranger",
  metric = "ROC",
  trControl = ctrl,
  num.trees = 500,
  importance = "impurity"
)

fit_glm
fit_rf

results <- resamples(list(
  LogisticRegression = fit_glm,
  RandomForest = fit_rf
))

summary(results)
bwplot(results, metric = "ROC")



###11_02-03_2026#####################################################################


library(Matrix)

atac_matrix <- readMM("GSE298008_11022026/GSE298008_ATAC_counts_matrix.mtx.gz")

regions <- read.table("GSE298008_11022026/GSE298008_ATAC_regions.tsv.gz", header = FALSE, stringsAsFactors = FALSE)

barcodes <- read.table("GSE298008_11022026/GSE298008_barcodes.tsv.gz", header = FALSE, stringsAsFactors = FALSE)


rownames(atac_matrix) <- regions[, 1]
colnames(atac_matrix) <- barcodes[, 1]

atac_matrix

cell_metadata <- read.delim(
  "GSE298008_11022026/GSE298008_cell_metadata.csv.gz",
  row.names = 1, sep = "",
  stringsAsFactors = FALSE
)
cell_metadata <- cell_metadata[colnames(atac_matrix), , drop = FALSE]

cell_metadata <- cell_metadata[
  colnames(atac_matrix),
  ,
  drop = FALSE
]


all(colnames(atac_matrix) %in% rownames(cell_metadata))



sce_atac <- SingleCellExperiment(
  assays = list(counts = atac_matrix),
  colData = cell_metadata,
  rowData = data.frame(
    region = rownames(atac_matrix)
  )
)

print(sce_atac)

BiocManager::install("DelayedMatrixStats")

atac_counts <- assay(sce_atac, "counts")
atac_counts <- as(atac_counts, "dgCMatrix")

# Calculate variance
peak_variance <- rowVars(atac_counts)

# Identify valid peaks
valid_peaks <- which(
  is.finite(peak_variance) &
    peak_variance > 0 &
    Matrix::rowSums(atac_counts) > 0
)

# Select top 2000 variable peaks
top_peaks <- valid_peaks[
  order(peak_variance[valid_peaks], decreasing = TRUE)
]

top_peaks <- head(top_peaks, 2000)

# Generate reduced SCE
sce_atac_subset <- sce_atac[top_peaks, ]



sce_atac_subset

# Calculate the threshold (e.g., peak must be in at least 5% of cells)
min_cells <- 0.05 * ncol(sce_atac_subset)

# Identify which peaks meet this threshold
keep_peaks <- rowSums(counts(sce_atac_subset) > 0) > min_cells

# Subset the SCE object
sce_atac_subset_filtered <- sce_atac_subset[keep_peaks, ]
sce_atac_subset_filtered
sum(is.na(sce_atac_subset_filtered$celltype_ordered))
set.seed(123)
sce_atac_subset_filtered

sim_result <- scdesign3(
  sce = sce_atac_subset_filtered,
  assay_use = "counts",
  celltype = "celltype_ordered", # Replace with the actual column name in your metadata
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  mu_formula = "celltype_ordered", 
  sigma_formula = "1",
  family_use = "nb", # ATAC-seq is highly sparse; 'poisson' or 'nb' are standard
  n_cores = 12,
  usebam = FALSE,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = TRUE
)


dim(sim_result$new_count)

logcounts(sce_atac_subset_filtered) <- log1p(counts(sce_atac_subset_filtered))
simu_sce <- SingleCellExperiment(list(counts = sim_result$new_count), colData = sim_result$new_covariate)
logcounts(sim_result) <- log1p(counts(sim_result))

sim_counts <- as.matrix(sim_result$new_count)

# Extract simulated metadata
sim_metadata <- as.data.frame(sim_result$new_covariate)

# Check dimensions
dim(sim_counts)
dim(sim_metadata)

# Create simulated SingleCellExperiment
simu_sce <- SingleCellExperiment(
  assays = list(
    counts = sim_counts,
    logcounts = log1p(sim_counts)
  ),
  colData = S4Vectors::DataFrame(sim_metadata)
)

# Check output
class(simu_sce)
dim(simu_sce)
assayNames(simu_sce)

class(sim_result)
class(sim_result$new_count)
class(sim_result$new_covariate)


set.seed(123)
sim_result


sim_counts <- sim_result$new_count
sim_metadata <- S4Vectors::DataFrame(sim_result$new_covariate)

# Check dimensions
stopifnot(ncol(sim_counts) == nrow(sim_metadata))

# Create SingleCellExperiment
simu_sce_experiment <- SingleCellExperiment(
  assays = list(
    counts = sim_counts,
    logcounts = log1p(sim_counts)
  ),
  colData = sim_metadata
)

# Verify
simu_sce
assayNames(simu_sce)
dim(counts(simu_sce))
dim(logcounts(simu_sce))



compare_figure <- plot_reduceddim(ref_sce = sce_atac_subset_filtered, 
                                  sce_list = list(simu_sce_experiment), 
                                  name_vec = c("Reference", "scDesign3"),
                                  assay_use = "logcounts", 
                                  if_plot = TRUE, 
                                  color_by = "celltype_ordered", 
                                  n_pc = 20)

plot(compare_figure$p_umap)



# Match genes shared by both datasets
common_genesATAC <- intersect(
  rownames(sce_atac_subset_filtered),
  rownames(simu_sce_experiment)
)

ref_counts <- assay(sce_atac_subset_filtered, "counts")[common_genesATAC, , drop = FALSE]
sim_counts <- assay(simu_sce_experiment, "counts")[common_genesATAC, , drop = FALSE]

ref_log <- assay(sce_atac_subset_filtered, "logcounts")[common_genesATAC, , drop = FALSE]
sim_log <- assay(simu_sce_experiment, "logcounts")[common_genesATAC, , drop = FALSE]

# Calculate metrics for each gene
gene_metricsATAC <- data.frame(
  gene = common_genesATAC,
  ref_mean = Matrix::rowMeans(ref_log),
  sim_mean = Matrix::rowMeans(sim_log),
  ref_var = apply(ref_log, 1, var),
  sim_var = apply(sim_log, 1, var),
  ref_zero = Matrix::rowMeans(ref_counts == 0),
  sim_zero = Matrix::rowMeans(sim_counts == 0)
)

# Relative agreement of gene means
gene_metricsATAC$mean_difference <-
  gene_metricsATAC$sim_mean - gene_metricsATAC$ref_mean

gene_metricsATAC$abs_mean_difference <-
  abs(gene_metricsATAC$mean_difference)

# Difference in zero proportions
gene_metricsATAC$zero_difference <-
  gene_metricsATAC$sim_zero - gene_metricsATAC$ref_zero

# Correlation of gene-level summaries
cor(
  gene_metricsATAC$ref_mean,
  gene_metricsATAC$sim_mean,
  method = "spearman"
)

cor(
  gene_metricsATAC$ref_var,
  gene_metricsATAC$sim_var,
  method = "spearman"
)

# Summarize errors across genes
summary(gene_metricsATAC$abs_mean_difference)
summary(abs(gene_metricsATAC$zero_difference))

par(mfrow = c(1, 3))

plot(
  gene_metricsATAC$ref_mean,
  gene_metricsATAC$sim_mean,
  xlab = "Reference mean log-expression",
  ylab = "Simulated mean log-expression",
  main = "Gene means",
  pch = 16, cex = 0.5
)
abline(0, 1, col = "red")

plot(
  gene_metricsATAC$ref_var,
  gene_metricsATAC$sim_var,
  xlab = "Reference variance",
  ylab = "Simulated variance",
  main = "Gene variances",
  pch = 16, cex = 0.5
)
abline(0, 1, col = "red")

plot(
  gene_metricsATAC$ref_zero,
  gene_metricsATAC$sim_zero,
  xlab = "Reference zero fraction",
  ylab = "Simulated zero fraction",
  main = "Zero proportions",
  pch = 16, cex = 0.5
)
abline(0, 1, col = "red")

par(mfrow = c(1, 1))



ks_results <- lapply(gene_metricsATAC, function(g) {
  x <- as.numeric(ref_log[g, ])
  y <- as.numeric(sim_log[g, ])
  
  test <- ks.test(x, y, exact = FALSE)
  
  data.frame(
    gene = g,
    ks_statistic = unname(test$statistic),
    p_value = test$p.value
  )
})

ks_results <- do.call(rbind, ks_results)

# Correct for testing many genes
ks_results$FDR <- p.adjust(
  ks_results$p_value,
  method = "BH"
)

head(
  ks_results[order(-ks_results$ks_statistic), ],
  20
)

#######################################################

library(caret)

set.seed(123)

# Use a manageable set of shared genes
# Ideally, choose highly variable genes using the reference data only
n_genes <- min(100, length(common_genesATAC))

selected_genes <- common_genesATAC[
  order(
    apply(ref_log, 1, var),
    decreasing = TRUE
  )[seq_len(n_genes)]
]

# Cells x genes
ref_x <- t(as.matrix(ref_log[selected_genes, , drop = FALSE]))
sim_x <- t(as.matrix(sim_log[selected_genes, , drop = FALSE]))

# Balance the number of cells
n <- min(nrow(ref_x), nrow(sim_x))

ref_x <- ref_x[sample(nrow(ref_x), n), , drop = FALSE]
sim_x <- sim_x[sample(nrow(sim_x), n), , drop = FALSE]

X <- rbind(ref_x, sim_x)

y <- factor(
  c(rep("Reference", n), rep("Simulated", n)),
  levels = c("Reference", "Simulated")
)

# Scale features
X <- scale(X)
X[!is.finite(X)] <- 0

library(caret)
# Cross-validated classification
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

classifier <- train(
  x = X,
  y = y,
  method = "glm",
  family = binomial(),
  metric = "ROC",
  trControl = ctrl
)

classifier

############################################

library(caret)
library(pROC)
library(ranger)

set.seed(123)

# X: cells x genes numeric matrix or data frame
# y: class labels for each cell

y <- factor(y, levels = c("Reference", "Simulated"))

stopifnot(nrow(X) == length(y))
stopifnot(all(table(y) > 0))

# Stratified 5-fold cross-validation
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Logistic regression
fit_glm <- train(
  x = X,
  y = y,
  method = "glm",
  family = binomial(),
  metric = "ROC",
  trControl = ctrl
)

# Random forest
fit_rf <- train(
  x = X,
  y = y,
  method = "ranger",
  metric = "ROC",
  trControl = ctrl,
  num.trees = 500,
  importance = "impurity"
)

fit_glm
fit_rf

results <- resamples(list(
  LogisticRegression = fit_glm,
  RandomForest = fit_rf
))

summary(results)
bwplot(results, metric = "ROC")





