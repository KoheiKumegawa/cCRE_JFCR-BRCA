#----------------------------------------------------------------------------
# 03_epiclustering.R
#----------------------------------------------------------------------------
library(ArchR)
library(ggplot2)
library(ggrepel)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
addArchRThreads(threads = 22)
addArchRGenome("hg19")

pt_ls <- read.csv("ref/patient_list_20231201.csv")
arc <- readRDS("rds/02_arc_ep.rds")

#----- clustering -----#
arc <- addIterativeLSI(arc, useMatrix = "TileMatrix", name = "IterativeLSI", force = T)
arc <- addUMAP(arc, reducedDims = "IterativeLSI", force = T)
arc <- addClusters(arc, reducedDims = "IterativeLSI", maxClusters = 40, resolution = 1.2, force = T)

# % of sample x cluster
sample_colors <- readRDS("rds/02_sample_colors.rds")
sampletype_colors <- readRDS("rds/02_sampletype_colors.rds")
subtype_colors <- c("ER+/HER2-" = "blue", "ER-/PGR+/HER-" = "green4", "ER+/HER2+" = "yellow2", "HER2+" = "orange", "TN" = "red", "NA" = "gray")


#UMAP visualization
p1 <- plotEmbedding(arc, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", plotAs = "points", size = 0.5)
p2 <- plotEmbedding(arc, colorBy = "cellColData", name = "Sample", embedding = "UMAP", plotAs = "points", size = 0.5, labelMeans = F, pal = sample_colors)
p3 <- plotEmbedding(arc, colorBy = "cellColData", name = "SampleType", embedding = "UMAP", plotAs = "points", size = 0.5, labelMeans = F, pal = sampletype_colors)
p3_2 <- plotEmbedding(arc, colorBy = "cellColData", name = "Receptor", embedding = "UMAP", plotAs = "points", size = 0.5, labelMeans = F, pal = subtype_colors)
plotPDF(p1, p2, p3, p3_2, name = "03_UMAP_EP.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

cluster_colors <- p1$plot_env$pal
names(cluster_colors) <- str_split(names(cluster_colors), pattern = "-", simplify = T)[,2]
saveRDS(cluster_colors, "rds/03_epicluster_colors.rds")

cM <- confusionMatrix(arc$Sample, arc$Clusters)[pt_ls$ID, paste0("C", c(1:35))] %>% as.matrix()
write.csv(cM, "output/Tables/03_ClusterSampleMatrix_EP.csv")

p5 <- ggplot(reshape2::melt(cM), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") + ArchR::theme_ArchR() + scale_y_continuous(labels = percent) +
  scale_fill_manual(values = cluster_colors) + labs(x = "Cluster", y = "%Sample") +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))
p6 <- ggplot(reshape2::melt(cM), aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") + ArchR::theme_ArchR() + scale_y_continuous(labels = percent) +
  scale_fill_manual(values = sample_colors) + labs(x = "Sample", y = "%Cluster") +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))

pdf("output/Plots/03_BarPlot_ClusterSample_EP.pdf", height = 4, width = 8)
p5
p6
dev.off()

#----- gene activity -----#
markers <- read.csv("ref/markers.csv", header = F)[1,] %>% as.character
arc <- addImputeWeights(arc, reducedDims = "IterativeLSI")

# UMAP overlay
p7 <- plotEmbedding(arc,
                    colorBy = "GeneScoreMatrix", 
                    name = markers, 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = 0.5, rastr = T)
plotPDF(p7, name = "03_UMAP_GSMarkers_EP.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

# Box/Violin plots
p8 <- plotGroups(
  ArchRProj = arc, 
  groupBy = "Clusters", 
  colorBy = "GeneScoreMatrix", 
  name = markers,
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
pdf("output/Plots/03_VlnBoxplot_GSMarkers_EP.pdf", width = 5, height = 2.5)
p8
dev.off()

saveRDS(arc, "rds/03_arc_ep.rds")

#----- cluster sample matrix heatmap -----#
cM <- read.csv("output/Tables/03_ClusterSampleMatrix_EP.csv", row.names = 1)
cM <- t(t(cM) / colSums(cM))*100
colnames(cM) <- gsub("C", "EP", colnames(cM))
col_fun1 <- colorRamp2(c(0,50,100), c("white","#8c96c6","#4d004b"))
ht1 <- Heatmap(cM, name = "%Sample in each cluster", cluster_columns = F, cluster_rows = F, col = col_fun1, border = "black",
               rect_gp = gpar(col = "gray", lwd = 0.5))
p9 <- draw(ht1)
pdf("output/Plots/03_Heatmap_confusion_EP.pdf", width = 8, height = 6)
p9
dev.off()

#----- gene activity heatmap -----#
GeneScoreClusters <- getMarkerFeatures(arc, 
                                       useMatrix = "GeneScoreMatrix", 
                                       groupBy = "Clusters",
                                       bias = c("TSSEnrichment", "log10(nFrags)"),
                                       testMethod = "wilcoxon")

mtx <- assays(GeneScoreClusters)$Mean
rownames(mtx) <- rowData(GeneScoreClusters)$name
mtx <- mtx[c("EPCAM","KRT5","KRT17","CDH1","CDH2","CD44","MYC","EGFR","TP53","CCND1","FOXA1","GRHL2","AR","VIM","ACTA2","ESR1","PGR","BCL2","SCUBE2","MKI67","AURKA","BIRC5","CCNB1","ERBB2","GRB7"),]
mtx <- t(scale(t(mtx)))


col_fun2 <- colorRamp2(c(-2,-1,0,1,2), c(viridis(5, option = "A")))
fh <- function(x) hclust(dist(x), method = "ward.D2")
ht2 <- Heatmap(t(mtx), name = "Mean GeneScore z-score", cluster_columns = fh, cluster_rows = F, col = col_fun2, border = "black",
               rect_gp = gpar(col = "gray", lwd = 0.5))
p10 <- draw(ht2)
pdf("output/Plots/03_Heatmap_markerGS_EP.pdf", width = 8, height = 6)
p10
dev.off()
