# Set this to the working directory that containing data files
# Dataset used is dowloaded from GEO (GSE107185_up-tf-stem.counts.tsv.gz, with accession number GSE107185)
setwd("/Users/Wanying/Desktop/IGP_2020_spring/MM2-1_Quantitative Systems Biology in Single Cells/QSBSC_Class_2020/Users/Wanying/Final_project/Remake_plot/")

# Please notice that there are some big changes in Seurat V3 compare to Seurat V2
# Seurat 3.1.4 is used in this assignment
library(Seurat)

# ---------------- Copy and paste some functions from swne package due to installation issues --------------
# Cannot install swne package, so pasted a few functions form the package here
ReadData <- function(matrix.dir, make.sparse = T) {
  if (dir.exists(matrix.dir)) {
    if (!grepl("\\/$", matrix.dir)) { matrix.dir <- paste(matrix.dir, "/", sep = "") }
    
    barcode.loc <- paste0(matrix.dir, "barcodes.tsv.gz")
    gene.loc <- paste0(matrix.dir, "features.tsv.gz")
    matrix.loc <- paste0(matrix.dir, "matrix.mtx.gz")
    
    if (!file.exists(barcode.loc)) {
      barcode.loc <- paste0(matrix.dir, "barcodes.tsv")
      if (!file.exists(barcode.loc)) stop("Barcode file missing")
    }
    if (!file.exists(gene.loc)) {
      genes.loc <- paste0(matrix.dir, "features.tsv")
      if (!file.exists(genes.loc)) {
        genes.loc <- paste0(matrix.dir, "genes.tsv")
        if (!file.exists(genes.loc)) stop("Gene name file missing")
      }
    }
    if (!file.exists(matrix.loc)) {
      matrix.loc <- paste0(matrix.dir, "matrix.mtx")
      if (!file.exists(matrix.loc)) stop("Expression matrix file missing")
    }
    
    counts <- Matrix::readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if (all(grepl(pattern = "\\-1$", cell.names))) {
      cell.names <- as.vector(as.character(sapply(cell.names, ExtractField, field = 1, delim = "-")))
    }
    rownames(counts) <- make.unique(names = as.character(sapply(gene.names, ExtractField, field = 2, delim = "\\t")))
    colnames(counts) <- cell.names
    
  } else {
    counts <- as.matrix(read.table(matrix.dir, sep = "\t", header = T, row.names = 1))
  }
  colnames(counts) <- make.names(colnames(counts))
  
  if (make.sparse) {
    counts <- as(counts, "dgCMatrix")
  }
  return(counts)
}


FilterData <- function(x, min.samples.frac, trim, min.nonzero.features = 500, max.sample.sum = 50000) {
  if (is.data.frame(x) | is.matrix(x)) {
    x <- as(as.matrix(x), "dgCMatrix")
  }
  
  min.cells <- round(ncol(x)*min.samples.frac)
  x <- x[ , Matrix::colSums(x) < max.sample.sum]
  x <- x[ , Matrix::colSums(x > 0) > min.nonzero.features]
  x <- x[Matrix::rowSums(x > 0) > min.cells, ]
  # if (trim > 0) {
  #   x <- t(winsorize_matrix(t(x), trim = trim))
  # }
  return(x)
}

winsorize_matrix <- function(mat, trim) {
  if(trim < 1) {
    trim <- trim*ncol(mat)
  }
  
  mat <- t(mat)
  inplaceWinsorizeSparseCols(mat, trim)
  return(t(mat))
}
# -------------------- End of pasted functions -------------------------

# Code is modified from github (hESC_run_clustering.R)
matrix.file <- "up-tf-stem.counts.tsv" ## Input counts matrix
regress.batch <- TRUE

min.cells.frac <- 0.01 ## Fraction of cells expressing gene in order to include
min.genes.exp <- 300 ## Minimum number of expressed genes for a cell

regress.model <- "negbinom" ## Model for regressing out unwanted effects
min.expr <- 0.1 ## Minimum average expression for variable genes
min.var.z <- 0.5 ## Minimum Z-score for variable genes

pcs.use <- 30 ## Number of PCs to use for clustering
cluster.res <- 1.8 ## Cluster resolution
min.conn <- 1e-4 ## Minimum SNN connectivity between clusters to calculate accuracy
acc.cutoff <- 0.95 ## Minimum classification accuracy to keep clusters separate

## Load the dataset, filter and trim counts matrix
counts <- ReadData(matrix.file)
counts <- FilterData(counts, min.cells.frac,trim, min.genes.exp)
nUMI <- Matrix::colSums(counts)

## Check for batch effects
batch <- factor(sapply(colnames(counts), function(x) strsplit(x, split = "\\.")[[1]][[2]]))
pd <- data.frame(nUMI, batch)

## Mitochondrial fraction
mito.genes <- grep("^MT-", rownames(counts), value = T)
pd$percent.mt <- Matrix::colSums(counts[mito.genes, ])/Matrix::colSums(counts)
counts <- counts[!rownames(counts) %in% mito.genes,]

## Run clustering pipeline via Seurat
se.obj <- CreateSeuratObject(counts = counts, project = "10x", meta.data = pd, min.cells = 3, min.features = 200)

## Find overdispersed genes for PCA
se.obj <- FindVariableFeatures(se.obj)

## Regress out confounders
# Scale the data
se.obj <- ScaleData(se.obj, vars.to.regress = c("nUMI", "percent.mt"),
                    features = se.obj@assays$RNA@var.features, model.use=regress.model)

## Run PCA and clustering
se.obj <- RunPCA(se.obj, npcs=40, features = se.obj@assays$RNA@var.features)
se.obj <- FindNeighbors(se.obj, dims = 1:30, k.param = 30, prune.SNN = 1/20)
se.obj <- FindClusters(se.obj, resolution = 0.38)

# se.obj <- RunUMAP(se.obj, dims = 1:10)
# se.obj <- RunTSNE(se.obj, features = se.obj@assays$RNA@var.features)
se.obj <- RunTSNE(se.obj, dims.use = 1:30) # Default is set to 30 by author
saveRDS(se.obj, file = "up-tf-stem.counts.rds") # Save the object for future use

# DimPlot(se.obj, reduction = "umap")
jpeg(filename = 'output.jpg', width=600, height=500, res=100, quality=400)
DimPlot(se.obj, reduction = "tsne")
dev.off()
