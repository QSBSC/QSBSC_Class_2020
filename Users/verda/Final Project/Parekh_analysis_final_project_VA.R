#########################PART 1#########################
Here, I have reproduced Figure 2B (Expression level changes of Key EMT Genes) from the Parekh paper. 
This code was pulled from https://github.com/yanwu2014/SEUSS-Analysis and specifically copied from `hESC_EMT_analysis.R`.
#########################PART 1#########################

## Analyze SNAI2/KLF4 effects on genes related to Epithelial-Mesenchymal transition in hPSC medium

## Barplot of key EMT markers -- Figure 2B 
library(ggplot2)
library(swne)

## Name of the sample to be analyzed
file.handle <- "/Users/vagan/Downloads/Parekh_paper_data/up-tf-stem"

## Load regression results 
coefs.df <- read.table(paste(file.handle, "regression.pvals.tsv", sep = "_"), sep = "\t", 
                       header = T, stringsAsFactors = F)
coefs <- UnflattenDataframe(coefs.df, "cf")

emt.genes.plot <- c("CDH1", "EPCAM", "LAMC1", "SPP1", "THY1", "VIM", "TPM2")

coefs.plot <- c(coefs[emt.genes.plot, "KLF4"], coefs[emt.genes.plot, "SNAI2"])
tfs.factor <- c(rep("KLF4", length(emt.genes.plot)), rep("SNAI2", length(emt.genes.plot)))
genes.factor <- rep(emt.genes.plot, 2)

barplot.df <- data.frame(coefs = coefs.plot, TF = tfs.factor, Gene = genes.factor)

pdf("up-tf-stem_EMT_barplot.pdf", width = 5, height = 4.5)
ggplot(data = barplot.df, aes(x = Gene, y = coefs)) + geom_bar(aes(fill = TF), position = "dodge", stat = "identity") + 
  theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                          axis.text.x = element_text(angle = 90, hjust = 1, size = 20), 
                          legend.text = element_text(size = 16))
dev.off()

------------------------------------------------------------------------------------------------------------------------
#########################PART 2#########################
Here, I have regenerated the tSNE projection for hESCs across media conditions as a PCA (linear dimensionality reduction)
and UMAP (nonlinear dimensionality reduction). You will see these plots at the bottom of the script in Step 7.
This code was pulled from https://github.com/yanwu2014/SEUSS-Analysis and specifically copied from `hESC_summarize_results.R`.
#########################PART 2#########################

#### NOTE
The following packages must be installed prior to running (many needed only to load `swne` package. Uncomment to execute.
                                                           
```{r, results=FALSE}
getPackageIfNeeded <- function(pkg) {
  if (!require(pkg, character.only=TRUE, quietly =TRUE))
    install.packages(pkgs=pkg, dependencies=TRUE)
}
getBiocPackageIfNeeded <- function(pkg, vers = NULL) {
  if (!require(pkg, character.only=TRUE, quietly =TRUE))
    BiocManager::install(pkgs=pkg, version=vers)
}

bcpkgs <- list(list(pkg="BiocManager",vers="3.10"),
               list(pkg="TxDb.Hsapiens.UCSC.hg38.knownGene",vers="3.10"),
               list(pkg="org.Hs.eg.db"),
               list(pkg="GO.db"),
               list(pkg="qvalue", vers="3.10")
)
invisible(lapply(bcpkgs, function(x) 
  do.call(getBiocPackageIfNeeded, args=x)))
pkgs	<-	c("devtools","uwot","umap","Seurat")
invisible(sapply(pkgs,getPackageIfNeeded))
if (!require("NNLM")) 
  devtools::install_version("NNLM", version = "0.4.3")
if (!require("swne")) 
  devtools::install_github("yanwu2014/swne")
if (!require("perturbLM")) 
  devtools::install_github("yanwu2014/perturbLM", upgrade = "never")
```

#### Step 1
Load all required libraries.
```{r}
#### Create t-SNE plots, differential expression heatmaps, and cluster enrichment heatmaps
library(perturbLM)
library(swne)
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(grid)
library(gtable)
```
#### Step 2
Define file locations and load data. Example below includes path information for my computer. Must be modified to work appropriately on your computer.

```{r}
## Name of the sample to be analyzed
file.handle <- "/Users/vagan/Downloads/Parekh_paper_data/up-tf-stem"
#file.handle <- "up-tf-stem"
# file.handle <- "up-tf-endo"
# file.handle <- "up-tf-multi"
# file.handle <- "up-tf-klf"
# file.handle <- "up-tf-myc"
# file.handle <- "up-tf-neuron-nohygro"

## Load regression results
coefs.df <- read.table(paste(file.handle, "regression.pvals.tsv", sep = "_"), sep = "\t", 
                       header = T, stringsAsFactors = F)
top.coefs.df <- subset(coefs.df, FDR < 0.05 & abs(cf) > 0.025)
print(sort(table(top.coefs.df$Group), decreasing = T))

## Load clustering
r <- readRDS(paste(file.handle, "clustering.Robj", sep = "_"))


clusters <- r@ident; names(clusters) <- r@cell.names;
levels(clusters) <- paste("C", levels(clusters), sep = "")
clusters.list <- UnflattenGroups(clusters)

# update v2 Seurat object to v3
r <- UpdateSeuratObject(r)

## Load genotypes
# genotypes.list <- ReadGenotypes("up-tf-klf-myc_pheno_dict.csv")
genotypes.list <- ReadGenotypes(paste(file.handle, "pheno_dict.csv", sep = "_"))
genotypes.list <- lapply(genotypes.list, function(x) x[x %in% names(clusters)])
genotypes.sizes <- sapply(genotypes.list, length)
```

#### Step 3
Calculate enrichment scores and summarize results.
```{r}
##  Calculate genotype enrichment
genotypes.cl.tbl <- GenotypeClusterCounts(genotypes.list, clusters.list)
genotypes.cl.pvals <- GenotypeClusterPvals(genotypes.cl.tbl)
genotypes.min.p <- apply(genotypes.cl.pvals, 1, function(x) min(abs(x)) * length(x))
metap::sumlog(genotypes.min.p)
## Summarize analysis results
genotype.counts <- sort(table(top.coefs.df$Group), decreasing = T)
genotype.counts <- genotype.counts[names(genotype.counts) %in% rownames(genotypes.cl.pvals)]
genotype.avg_cf <- tapply(top.coefs.df$cf, top.coefs.df$Group, function(x) mean(abs(x)))[names(genotype.counts)]
genotype.summary <- data.frame(n_diff_genes = as.integer(genotype.counts), avg_cf = genotype.avg_cf)
genotype.summary$min_cluster_pval <- apply(genotypes.cl.pvals[rownames(genotype.summary), ], 1, function(x) min(abs(x)))
genotype.summary$min_cluster <- apply(genotypes.cl.pvals[rownames(genotype.summary), ], 1, function(x) names(x[which.min(abs(x))]))
genotype.summary$n_cells <- genotypes.sizes[rownames(genotype.summary)]
```

Uncomment code below if you want to save the output data to a file.
```{r}
# write.table(genotype.summary, file = paste(file.handle, "summary.tsv", sep = "_"), sep = "\t")
```

#### Step 4
Continue processing to identify significant genotypes.
```{r}
## Determine significant genotypes
min.cl.pval <- 1e-12
min.diff.genes <- 40
ctrl.cl <- genotype.summary["mCherry", "min_cluster"]
sig.genotypes <- rownames(subset(genotype.summary, (abs(min_cluster_pval) < min.cl.pval & min_cluster != ctrl.cl) 
                                 | n_diff_genes > min.diff.genes))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry-int-ctrl")), decreasing = T); print(sig.genotypes);
```

#### Step 5
Generate barplot graph.
```{r}
## Differential expression barplot
n.diff.genes <- genotype.summary[sig.genotypes, "n_diff_genes"]
names(n.diff.genes) <- sig.genotypes
gg.bar <- ggBarplot(n.diff.genes, fill.color = "purple")

# pdf(paste(file.handle, "diff_genes_barplot.pdf", sep = "_"), width = 3.5, height = 1.5)
# print(gg.bar)
# dev.off()

## Cluster enrichment heatmap (Figure 1c-e)
max.lp <- 50

genotypes.cl.lp <- apply(genotypes.cl.pvals[sig.genotypes,], 1:2, function(x) ifelse(x > 0, -log10(x), log10(-1 * x)))
genotypes.cl.lp[genotypes.cl.lp > max.lp] <- max.lp
genotypes.cl.lp[genotypes.cl.lp < -1 * max.lp] <- -1 * max.lp

# pdf(paste(file.handle, "cl_enrich_heatmap.pdf", sep = "_"), width = 4.5, height = 5)
# ggHeat(genotypes.cl.lp, clustering = "none", x.lab.size = 14, y.lab.size = 14)
# dev.off()

gg.heat <- ggHeat(t(genotypes.cl.lp), clustering = "none", x.lab.size = 14, y.lab.size = 14)
gg.bar <- gg.bar + theme(axis.text.x = element_blank())
gg.heat <- gg.heat + theme(legend.position = "none")

gg.1 <- ggplotGrob(gg.heat)
gg.2 <- ggplotGrob(gg.bar)

gg.aligned <- rbind(gg.2, gg.1, size = "first")
gg.aligned$widths <- grid::unit.pmax(gg.1$widths, gg.2$widths)

#pdf(paste(file.handle, "stacked_diff_genes_cl_enrich_nolabels.pdf", sep = "_"), width = 5.5, height = 6)
grid::grid.newpage()
grid::grid.draw(gg.aligned)
#dev.off()
```

#### Step 6
Generate tSNE plots
```{r}
## tSNE plots (Figure 1c-e)
plot.seed <- 32907098
tsne.emb <- Embeddings(object = r, reduction = "tsne")

#pdf(paste(file.handle, "tsne_plot.pdf", sep = "_"), width = 4, height = 4)
PlotDims(tsne.emb, sample.groups = clusters, pt.size = 0.75, alpha.plot = 0.5, label.size = 6, do.label = T,
         show.legend = F, seed = plot.seed)
#dev.off()

#pdf(paste(file.handle, "tsne_plot_nolabel.pdf", sep = "_"), width = 4, height = 4)
PlotDims(tsne.emb, sample.groups = clusters, pt.size = 0.75, alpha.plot = 0.5, label.size = 0, do.label = F,
         show.legend = F, seed = plot.seed)
#dev.off()

#pdf(paste(file.handle, "tsne_plot_batch.pdf", sep = "_"), width = 4, height = 4)
batch <- factor(r@meta.data$batch)
# names(batch) <- r@cell.names
names(batch) <- colnames(r)

batch <- batch[rownames(tsne.emb)]

PlotDims(tsne.emb, sample.groups = batch, pt.size = 1, alpha.plot = 0.5, label.size = 8, do.label = T,
         show.legend = F, seed = plot.seed)
#dev.off()


## tSNE plot overlay
# plot.seed <- 32907098
# tsne.emb <- Embeddings(object = r, reduction = "tsne")

TF <- "NEUROD1"

tf.groups <- as.character(clusters); names(tf.groups) <- names(clusters);
tf.groups[names(tf.groups) %in% genotypes.list[[TF]]] <- TF
tf.groups[!names(tf.groups) %in% genotypes.list[[TF]]] <- ""
tf.groups <- factor(tf.groups)

#pdf(paste0(file.handle, "_tsne_plot_", TF, ".pdf"), width = 3, height = 3)
PlotDims(tsne.emb, sample.groups = tf.groups, pt.size = 0.35, alpha.plot = 0.4, label.size = 0, do.label = T,
         show.legend = F, seed = plot.seed) + scale_color_manual(values = c("grey", "red")) + ggtitle(TF)
#dev.off()
```

#### Step 7 

## Using UMAP visualization in the Seurat package 
r <- RunUMAP(r, reduction = "pca", dims = 1:10)

## This is swne on the UMAP embeddings, which captures both the local and global dataset structure 
## Embeddings -- Dimensional Reduction Cell Embeddings Accessor Function
umap.emb <- Embeddings(r, reduction = "umap")
PlotDims(umap.emb, sample.groups = clusters, pt.size = 0.75, alpha.plot = 0.5, label.size = 6, do.label = T, show.legend = F, seed = plot.seed) 

## These are the 3 types of dimensionality reduction we learned: UMAP, t-SNE, and PCA 
## Both UMAP and t-SNE identify high dimensional struture (i.e., nonlinear dimensionality reduction), visualizing cell clusters, although distances are generally not interpretable
## PCA identifies structure through linear transformation and shows clusters of samples based on their similarity); distances are interpretable (i.e., linear dimensionality reduction ) 
DimPlot(r, reduction = "umap")
DimPlot(r, reduction = "tsne")
DimPlot(r, reduction = "pca")
