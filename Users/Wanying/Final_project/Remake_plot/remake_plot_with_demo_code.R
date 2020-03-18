# This code is modified from demo code of instructor
# Load .Robj object directly for plotting instead of filtering data from scratch

# Please set this to the directory containing the .Robj file
setwd("/Users/Wanying/Desktop/IGP_2020_spring/MM2-1_Quantitative Systems Biology in Single Cells/QSBSC_Class_2020/Users/Wanying/Final_project/Remake_plot/")

# Please notice that there are some big changes in Seurat V3 compare to Seurat V2
# Seurat 3.1.4 is used in this assignment
library(Seurat)

## Load clustering
r <- readRDS("up-tf-stem_clustering.Robj")

# update v2 Seurat object to v3
r <- UpdateSeuratObject(r)

## tSNE plots (Figure 1c-e)
plot.seed <- 32907098
tsne.emb <- Embeddings(object = r, reduction = "tsne")

jpeg(filename = 'output_with_demo_code.jpg', width=600, height=500, res=100, quality=400)
DimPlot(r, reduction = "tsne")
dev.off()
