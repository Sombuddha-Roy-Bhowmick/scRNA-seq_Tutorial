library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(reshape)
library(ggplot2)

#Load the NSCLC counts matrix
counts_matrix_filename = "filtered_gene_bc_matrices/GRCh38"
counts <- Read10X(data.dir = counts_matrix_filename)

#Filtering low-quality cells
counts_per_cell <- Matrix::colSums(counts)
counts_per_gene <- Matrix::rowSums(counts)
genes_per_cell <- Matrix::colSums(counts>0) # count a gene only if it has non-zero reads mapped
cells_per_gene <- Matrix::rowSums(counts>0) # only count cells where the gene is expressed
png(filename = "counts_per_cell_histogram.png") 
hist(log10(counts_per_cell + 1), main = 'counts per cell', col = 'blue') 
dev.off()

png(filename ="genes_per_cell_histogram.png")
hist(log10(genes_per_cell+1), main='genes per cell', col='blue')
dev.off()

png(filename ="counts_Vs_genes_per_cell.png")
plot(counts_per_cell, genes_per_cell, log='xy', col='blue') 
dev.off()

png(filename ="counts_per_gene_histogram.png")
hist(log10(counts_per_gene+1),main='counts per gene',col='blue')
dev.off()

png(filename ="genes_per_cell_ordered.png")
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)') # Plotting cells ranked by their number of detected genes
dev.off()

#Cell Filtering
MIN_GENES_PER_CELL <- 300
MAX_GENES_PER_CELL <- 5000
png(filename ="genes_per_cell_ordered_with_threshold.png")
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
abline(h=MIN_GENES_PER_CELL, col='blue')
abline(h=MAX_GENES_PER_CELL, col='yellow') 
dev.off()

#Examining percent mitochondrial read content
mito_genes <- grep("^mt-", rownames(counts) , ignore.case=T, value=T) #Defining mitochondrial genes
mito_gene_read_counts = Matrix::colSums(counts[mito_genes,]) 
pct_mito = mito_gene_read_counts / counts_per_cell * 100 #computing percentage mitochondrial reads
png(filename ="cells_sorted_by_percentage_mitochondrial_counts.png")
plot(sort(pct_mito), xlab = "cells sorted by percentage mitochondrial counts", ylab = "percentage mitochondrial counts")
dev.off()
MAX_PCT_MITO <- 10 #Setting the maximum allowed percent mitochondrial reads
png(filename ="maximum_mitochondrial_reads_allowed.png")
plot(sort(pct_mito))
abline(h=MAX_PCT_MITO, col='red')
dev.off()

#Beginning Analysis with Seurat
seurat <- CreateSeuratObject(counts, min.cells = 3, min.features = 350, project = "10X_NSCLC")
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-") # adding  % mit content to the metadata
png(filename ="Violin_Plot.png")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1) + NoLegend()
dev.off()
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
png(filename ="feature_scatter_plot.png")
CombinePlots(plots = list(plot1, plot2))
dev.off()


png(filename ="violin_plot_nFeature_RNA.png")
VlnPlot(object = seurat, features = c("nFeature_RNA"), group.by = c('orig.ident'), pt.size = 0.1) +NoLegend()
dev.off()


seurat <- subset(seurat, subset = nFeature_RNA > 350 & nFeature_RNA < 5000 & percent.mt < 10) #filtering

#Normalization
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#Detection of variable genes across the single cells
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(seurat), 20) # identify the 20 most highly variable genes
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
png(filename ="Top20_genes.png")
CombinePlots(plots = list(plot1, plot2))
dev.off()

#Scaling the data
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

#Linear dimensionality reduction
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
png(filename ="PCA_Dimensions.png")
VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
dev.off()
png(filename ="PC1_PC2_plot.png")
DimPlot(seurat, reduction = "pca") + NoLegend()
dev.off()
png(filename ="PCA_Heatmap.png")
DimHeatmap(seurat, dims = 1:3, cells = 500, balanced = TRUE)
dev.off()
png(filename ="PCA_feature_plot.png")
FeaturePlot(seurat,  reduction='pca', features=c("nFeature_RNA",  "percent.mt"))
dev.off()

save(seurat, file = "nsclc.Rda") #saving the seurat object

#Determining the ‘dimensionality’ of the dataset
png(filename ="Elbow_plot.png")
ElbowPlot(seurat)
dev.off()

#Visualization by tSNE and UMAP
seurat <- RunTSNE(seurat,dims = 1:20)
png(filename ="tSNE_plot.png")
dev.off()

seurat <- RunUMAP(seurat, dims = 1:10)
png(filename ="UMAP_plot.png")
DimPlot(seurat, reduction = "umap")  + NoLegend()
dev.off()

#Clustering the cells
seurat <- FindNeighbors(object = seurat, dims = 1:20)
seurat <- FindClusters(object = seurat, resolution = 0.5)

png(filename ="tSNE1_tSNE2_plot.png")
TSNEPlot(object = seurat, label=T)+ NoLegend()
dev.off()

#Comparing tSNE and UMAP visualizations
p1 <- TSNEPlot(object = seurat, label=T) + NoLegend()
p2 <- UMAPPlot(object = seurat, label=T) + NoLegend()
png(filename ="Combining_tSNE_UMAP_plots.png")
CombinePlots(plots = list(p1, p2))
dev.off()

#Mapping gene expression onto cell clusters and assigning identity to clusters

genes1 <- c("CD3D", "CD3E", "CD4","CD8A") # T-cell markers
genes2 <- c( "FOXP3", "IL7R")  # T-reg markers
genes3 <- c("GZMA", "GZMB", "PRF1", "NKG7", "GNLY") # NK cell markers
genes4 <-  c("CD19", "MS4A1") # B-cell markers
genes5 <- c("CD14", "LYZ", "FCGR3A") # Monocyte/Macrophage markers

png(filename ="Tcell_Markers_umap.png")
FeaturePlot(object = seurat, features = genes1, cols = c("grey", "blue"), reduction = "umap")
dev.off()
png(filename ="Tcell_Markers_tsne.png")
FeaturePlot(object = seurat, features = genes1, cols = c("grey", "blue"), reduction = "tsne")
dev.off()
png(filename ="Treg_Markers_umap.png")
FeaturePlot(object = seurat, features = genes2, cols = c("grey", "blue"), reduction = "umap")
dev.off()
png(filename ="Treg_Markers_tsne.png")
FeaturePlot(object = seurat, features = genes2, cols = c("grey", "blue"), reduction = "tsne")
dev.off()
png(filename ="NKcell_Markers_umap.png")
FeaturePlot(object = seurat, features = genes3, cols = c("grey", "blue"), reduction = "umap")
dev.off()
png(filename ="NKcell_Markers_tsne.png")
FeaturePlot(object = seurat, features = genes3, cols = c("grey", "blue"), reduction = "tsne")
dev.off()
png(filename ="Bcell_Markers_umap.png")
FeaturePlot(object = seurat, features = genes4, cols = c("grey", "blue"), reduction = "umap")
dev.off()
png(filename ="NKcell_Markers_tsne.png")
FeaturePlot(object = seurat, features = genes4, cols = c("grey", "blue"), reduction = "tsne")
dev.off()
png(filename ="Monocyte_Macrophage_Markers_umap.png")
FeaturePlot(object = seurat, features = genes5, cols = c("grey", "blue"), reduction = "umap")
dev.off()
png(filename ="Monocyte_Macrophage_Markers_tsne.png")
FeaturePlot(object = seurat, features = genes5, cols = c("grey", "blue"), reduction = "tsne")
dev.off()

#Labeling cell subsets
new.cluster.ids <- c("B cells", "CD4+ T cells", "2", "CD8+ T cells", "Monocytes/Macrophages", "Monocytes/Macrophages", "NK cells", "Monocytes/Macrophages", "Treg cells" , "9", "10", "11", "12", "13", "14") 
names(new.cluster.ids) <- levels(seurat)
seurat <- RenameIdents(seurat, new.cluster.ids)
png(filename ="Labelled_cell_subset_tsne.png")
DimPlot(seurat, reduction = "tsne", label = TRUE, pt.size = 0.1) + NoLegend()
dev.off()

#Using differentially expressed genes to help identify cell subsets
MIN_LOGFOLD_CHANGE = 2 # set to minimum required average log fold change in gene expression.
MIN_PCT_CELLS_EXPR_GENE = .25  # minimum percent of cells that must express gene in either cluster.

all.markers = FindAllMarkers(seurat, min.pct = MIN_PCT_CELLS_EXPR_GENE, logfc.threshold = MIN_LOGFOLD_CHANGE, only.pos = TRUE) # Finding all markers

all.markers.sortedByPval = all.markers[order(all.markers$p_val),] # sorting all markers by p-value

top10 <- all.markers.sortedByPval %>%  group_by(cluster)  %>% do(head(., n=20)) 
png(filename ="Top20_marker_genes_each_cluster_heatmap.png")
DoHeatmap(object = seurat, features = top10$gene)+ NoLegend() + theme(axis.text.y = element_text(color = "black", size = 5))
dev.off()
save(seurat, all.markers, file = "~/nsclc.Rda")

png(filename ="Top5_markers_in_tsne_plots.png")
FeaturePlot(seurat, features = all.markers.sortedByPval$gene[1:4], reduction = "tsne")
dev.off()

# Identification of genes uniquely differentially expressed in each cluster
genes_uniquely_DE = all.markers.sortedByPval %>% dplyr::filter(avg_log2FC >= MIN_LOGFOLD_CHANGE) %>% group_by(gene) %>%  summarize(n=n()) %>%  dplyr::filter(n==1) # differentially expressed genes that are only differentially expressed in one cluster of cells
genes_uniquely_DE.markers.sortedByPval = all.markers.sortedByPval[all.markers.sortedByPval$gene %in% genes_uniquely_DE$gene & all.markers.sortedByPval$avg_log2FC >= MIN_LOGFOLD_CHANGE,] # Sort DE genes by p value
top_marker_each = genes_uniquely_DE.markers.sortedByPval %>% dplyr::group_by(cluster) %>% do(head(., n=5))  # top 5 markers for each cluster
png(filename ="Heatmap_Top5_markers_each_cluster.png")
DoHeatmap(object = seurat, features = top_marker_each$gene)+ NoLegend() + theme(axis.text.y = element_text(color = "black", size = 5))
dev.off()

#Visualisation of marker genes : violin plots
for (i in 1:5) {
  plot <- VlnPlot(seurat, features = top_marker_each$gene[i], pt.size = 0.1) + NoLegend()
  print(plot)
  ggsave(filename = paste0("violin_plot_marker_gene_", i, ".png"), plot = plot, width = 6, height = 4)
}  # first five marker genes

#Visualisation of marker genes : dot plots
png(filename ="Dot_Plot_Marker_Genes_With_Labels.png")
DotPlot(seurat, features=unique(top_marker_each$gene),  dot.scale = 6)  + RotatedAxis() +theme(axis.text.x = element_text(color = "black", size = 7))
dev.off()

png(filename ="Dot_Plot_Marker_Genes_Without_Labels.png")
DotPlot(seurat, features=unique(top_marker_each$gene),  dot.scale = 6)  + RotatedAxis() +theme(axis.text.x = element_text(color = "black", size = 7)) + NoLegend()
dev.off()







