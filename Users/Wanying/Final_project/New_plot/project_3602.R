library(Seurat)
library(dplyr)

# data_dir = '/home/libraseq/wennie/3602/RNA_analysis/outs/raw_feature_bc_matrix/'
data_dir = '/home/libraseq/wennie/3602/RNA_analysis/outs/filtered_feature_bc_matrix/'
output_dir = paste("/home/libraseq/wennie/3602_SeuratOutput/RNA_analysis_with_filtered_data/Some_removal/", Sys.Date(), sep='')

libraseq.data <- Read10X(data.dir = data_dir)
# Only include features that are detected in >=3 cells
libraseq <- CreateSeuratObject(counts = libraseq.data, project = "3602", min.cells = 3)

libraseq[["percent.mt"]] = PercentageFeatureSet(libraseq, pattern='^MT-')
# Remove samples with more than 5% mitochondrial counts
libraseq <- subset(libraseq, subset = percent.mt < 5)

libraseq <- NormalizeData(libraseq, normalization.method = "LogNormalize", scale.factor = 10000)
libraseq <- FindVariableFeatures(libraseq, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(libraseq)
libraseq <- ScaleData(libraseq, features = all.genes)

libraseq <- RunPCA(libraseq, features = VariableFeatures(object = libraseq))

# Cluster the cells
# run with the first 10 dimention of PCs
libraseq <- FindNeighbors(libraseq, dims = 1:10)
libraseq <- FindClusters(libraseq, resolution = 0.5)

# Run non-linear dimensional reduction (UMAP or tSNE)
# libraseq <- RunTSNE(libraseq, dims = 1:10) # can also run on certain features
# DimPlot(libraseq, reduction = "tsne")
libraseq <- RunUMAP(libraseq, dims = 1:10)
jpeg(filename = paste(output_dir,'some_removal_cluster_UMAP.jpg'), width=900, height=480, res=150)
DimPlot(libraseq, reduction = "umap")
dev.off()

saveRDS(libraseq, file = paste(output_dir, "some_removal_3602_object.rds"))

# use write.csv() funciton to save each slot of the seurat object, such as assay$RNA$counts
# Can further process exported file in python

output_dir =  "/home/libraseq/wennie/3602_SeuratOutput/RNA_analysis_with_filtered_data/Some_removal/01 Export_csv_files/2019-12-11"

file_to_save = libraseq@assays$RNA@counts
write.csv(file_to_save, file = paste(output_dir,'metaData_RNA_counts.csv'))

file_to_save=libraseq@assays$RNA@counts
write.csv(file_to_save, file=paste(output_dir,'RNA_counts.csv'))

file_to_save=libraseq@reductions$pca@cell.embeddings
write.csv(file_to_save, file=paste(output_dir,'cell_embeddings_pca.csv'))

file_to_save=libraseq@reductions$umap@cell.embeddings
write.csv(file_to_save, file=paste(output_dir,'reduction_umap.csv'))

file_to_save=libraseq@reductions$tsne@cell.embeddings
write.csv(file_to_save, file=paste(output_dir,'reduction_tsne.csv'))



# Finding differentially expressed features (cluster biomarkers)
# Find all markers of cluster 0
# Features can be sorted by avg_logFC
cluster0.markers <- FindMarkers(libraseq, ident.1 = 0, min.pct = 0.25)

# Find variable markers for each cluster
libraseq.markers <- FindAllMarkers(libraseq, min.pct = 0.25, logfc.threshold = 0.25)
# Write to a .csv file for later use

jpeg(filename = paste(output_dir,'markerExpression_allMarkers_Cluster1.jpg'), width=4000, height=6000, res=150)
FeaturePlot(libraseq, features=cluster0.markers, reduction='tsne')

# Further separate top features, get top 147 features by adj_p_val (>0.05)
# Then separate positive and negative changes by avg_logFC
top147_by_p_val_adj = ordered_markers[1:147,]
neg_only = top147_by_p_val_adj[order(top147_by_p_val_adj$avg_logFC),][1:50,]
pos_only = top147_by_p_val_adj[order(top147_by_p_val_adj$avg_logFC, decreasing = TRUE),][1:50,]

write.csv(top147_by_p_val_adj, "2019-12-12 cluster0_top147_by_p_val_adj.csv")
write.csv(neg_only, "2019-12-12 cluster0_top50_neg_by_p_val_adj.csv")
write.csv(pos_only, "2019-12-12 cluster0_top50_pos_by_p_val_adj.csv")

# Again plot those top markers
jpeg(filename = paste(output_dir,'markerExpression_top50_neg_cluster0.jpg'), width=4000, height=6000, res=150)
FeaturePlot(libraseq, features=neg_markers, reduction='tsne')
dev.off()

