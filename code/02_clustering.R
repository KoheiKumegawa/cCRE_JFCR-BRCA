#----------------------------------------------------------------------------
# 02_clustering.R
#----------------------------------------------------------------------------
library(ArchR)
library(ggplot2)
library(ggrepel)
library(scales)
addArchRThreads(threads = 22)
addArchRGenome("hg19")

pt_ls <- read.csv("ref/patient_list_20231201.csv")
arc <- readRDS("rds/01_arc.rds")

#----- clustering -----#
arc <- addIterativeLSI(arc, useMatrix = "TileMatrix", name = "IterativeLSI", force = T)
arc <- addUMAP(arc, reducedDims = "IterativeLSI", force = T)
arc <- addClusters(arc, reducedDims = "IterativeLSI", maxClusters = 40, resolution = 1.2, force = T)

#UMAP visualization
p1 <- plotEmbedding(arc, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", plotAs = "points", size = 0.5)
p2 <- plotEmbedding(arc, colorBy = "cellColData", name = "Sample", embedding = "UMAP", plotAs = "points", size = 0.5, labelMeans = F)
p3 <- plotEmbedding(arc, colorBy = "cellColData", name = "SampleType", embedding = "UMAP", plotAs = "points", size = 0.5, labelMeans = F)
plotPDF(p1, p2, p3, name = "02_UMAP.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

# % of sample x cluster
cluster_colors <- p1$plot_env$pal
names(cluster_colors) <- str_split(names(cluster_colors), pattern = "-", simplify = T)[,2]
sample_colors <- p2$plot_env$pal
names(sample_colors) <- str_split(names(sample_colors), pattern = "-", simplify = T)[,2]
sampletype_colors <- p3$plot_env$pal
names(sampletype_colors) <- str_split(names(sampletype_colors), pattern = "-", simplify = T)[,2]

saveRDS(sample_colors, "rds/02_sample_colors.rds")
saveRDS(cluster_colors, "rds/02_cluster_colors.rds")
saveRDS(sampletype_colors, "rds/02_sampletype_colors.rds")

cM <- confusionMatrix(arc$Sample, arc$Clusters)[pt_ls$ID, paste0("C", c(1:39))] %>% as.matrix()
write.csv(cM, "output/Tables/02_ClusterSampleMatrix.csv")

p5 <- ggplot(reshape2::melt(cM), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") + ArchR::theme_ArchR() + scale_y_continuous(labels = percent) +
  scale_fill_manual(values = cluster_colors) + labs(x = "Cluster", y = "%Sample") +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))
p6 <- ggplot(reshape2::melt(cM), aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") + ArchR::theme_ArchR() + scale_y_continuous(labels = percent) +
  scale_fill_manual(values = sample_colors) + labs(x = "Sample", y = "%Cluster") +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))

pdf("output/Plots/02_BarPlot_ClusterSample.pdf", height = 4, width = 8)
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
plotPDF(p7, name = "02_UMAP_GSMarkers.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

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
pdf("output/Plots/02_VlnBoxplot_GSMarkers.pdf", width = 5, height = 2.5)
p8
dev.off()

#cluster-specific activated genes
GeneScoreClusters <- getMarkerFeatures(arc, 
                                       useMatrix = "GeneScoreMatrix", 
                                       groupBy = "Clusters",
                                       bias = c("TSSEnrichment", "log10(nFrags)"),
                                       testMethod = "wilcoxon")
MarkerClusters <- getMarkers(seMarker = GeneScoreClusters, cutOff = "FDR < 0.01 & Log2FC >= 1")
MarkerClustersTop100 <- lapply(names(MarkerClusters), function(i) MarkerClusters[[i]]$name[c(1:100)]) %>% do.call(cbind, .)
colnames(MarkerClustersTop100) <- names(MarkerClusters)

#export
write.csv(MarkerClustersTop100, "output/Tables/02_MarkerClustersTop100Genes.csv")
saveRDS(GeneScoreClusters, "rds/02_GeneScoreClusters.rds")

#----- Epithelial score -----#
krt_type2 <- paste0("KRT", c(1,2,3,4,5,"6A","6B","6C",7,8,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86))
krt_type1 <- paste0("KRT", c(9,10,12,13,14,15,16,17,18,19,20,23,24,25,26,27,28,31,32,"33A","33B",34,35,36,37,38,39,40))

idx <- which(rowData(GeneScoreClusters)$name %in% c(krt_type1, krt_type2, "EPCAM", "CDH1"))
episcore <- assays(GeneScoreClusters)$Mean[idx,] %>% colSums
episcore <- sort(episcore)
episcore <- data.frame(cluster = names(episcore), score = episcore)

p9 <- ggplot(episcore, aes(x = reorder(cluster, score), y = score)) + geom_point() + theme_ArchR() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "Cluster", y = "Epithelial score") +
  geom_hline(yintercept = 10, lty = "dashed")
pdf("output/Plots/02_DotEpiScore.pdf", height = 4, width = 6)
p9
dev.off()

#----- Cluster annotation -----#
cellTypeAnno <- list(T_NK_ILC = paste0("C", c(5,6,8,9,10)),
                     B = paste0("C", c(11,7)),
                     Plasma = paste0("C", c(12)),
                     Myeloid = paste0("C", c(18,19)),
                     Endo = paste0("C", c(15)),
                     Fib = paste0("C", c(13,14)),
                     Epithelial = paste0("C", c(1:4,16:17,20:39)))
cellType <- arc$Clusters
for(i in names(cellTypeAnno)){
  idx <- cellTypeAnno[[i]]
  cellType[which(cellType %in% idx)] <- i
}
arc$cellType <- cellType

col1 <- c(Endo = "#D51F26", B = "#272E6A", Epithelial = "#208A42", Myeloid = "#89288F", T_NK_ILC =  "#F47D2B", Fib = "#FEE500", Plasma = "#8A9FD1")
p9 <- plotEmbedding(arc, colorBy = "cellColData", name = "cellType", embedding = "UMAP", plotAs = "points", size = 0.5, labelMeans = F, pal = col1)
col2 <- c("ER+/HER2-" = "blue", "ER-/PGR+/HER-" = "green4", "ER+/HER2+" = "yellow2", "HER2+" = "orange", "TN" = "red", "NA" = "gray")
p10 <- plotEmbedding(arc, colorBy = "cellColData", name = "Receptor", embedding = "UMAP", plotAs = "points", size = 0.5, labelMeans = F, pal = col2)
plotPDF(p9,p10, name = "02_UMAP_Type.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

arc_ep <- arc[arc$cellType == "Epithelial"]
saveRDS(arc_ep, "rds/02_arc_ep.rds")

#----- genomeTrack -----#
markers <- c(markers, "IGLL5")
markers <- markers[-27]

gr <- ArchR::geneAnnoHg19$genes[ArchR::geneAnnoHg19$genes$symbol %in% markers]
gr <- extendGR(gr, upstream = 5000, downstream = 5000)
p11 <- plotBrowserTrack(arc, groupBy = "Clusters", region = gr, 
                        plotSummary = c("bulkTrack", "geneTrack"),
                        useGroups = paste0("C", c(1:4,16:17,20:39,15,13:14,5:6,8:10,7,11,12,18:19)))
# Epi -> Endo -> Fib -> T -> B -> Plasma -> Myeloid
plotPDF(p11, name = "02_GenomeTrackMarkers.pdf", ArchRProj = arc, addDOC = FALSE, width = 6, height = 8)
