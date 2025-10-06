#----------------------------------------------------------------------------
# 02_clustering_all.R
#----------------------------------------------------------------------------
library(ArchR)
library(scales)
library(ggsignif)
library(viridis)
library(ComplexHeatmap)
library(circlize)
addArchRThreads(threads = 24)
addArchRGenome("hg19")

arc <- readRDS("rds/01_arc.rds")
pt_ls <- read.csv("ref/patient_list_v2.csv")
SampleTypeColors <- c("Primary" = "#666666", "LocalRecurrence" = "#1B9E77", "LymphNode" = "#E7298A", "Pleural" = "#7570B3", "Ascitis" = "#E6AB02")

#----- Clustering visualization (same as 01_qc) -----#
p1 <- plotEmbedding(arc, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", plotAs = "points", size = 0.5)
p2 <- plotEmbedding(arc, colorBy = "cellColData", name = "Sample", embedding = "UMAP", plotAs = "points", size = 0.5, labelMeans = F)
p3 <- plotEmbedding(arc, colorBy = "cellColData", name = "SampleType", embedding = "UMAP", plotAs = "points", pal = SampleTypeColors, size = 0.5, labelMeans = F)
plotPDF(p1,p2,p3, name = "02_UMAP.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

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

cM <- confusionMatrix(arc$Sample, arc$Clusters)[pt_ls$ID, paste0("C", c(1:33))] %>% as.matrix()
write.csv(cM, "output/Tables/02_ClusterSampleMatrix.csv")

p4 <- ggplot(reshape2::melt(cM), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") + ArchR::theme_ArchR() + scale_y_continuous(labels = percent) +
  scale_fill_manual(values = cluster_colors) + labs(x = "Cluster", y = "%Sample") +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))
p5 <- ggplot(reshape2::melt(cM), aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") + ArchR::theme_ArchR() + scale_y_continuous(labels = percent) +
  scale_fill_manual(values = sample_colors) + labs(x = "Sample", y = "%Cluster") +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))

cM <- reshape2::melt(cM)
cM$Var2 <- factor(cM$Var2, levels = c("C1", "C2", "C3", "C4", "C14", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27", "C28", "C29", "C30", "C31", "C32", "C33", "C5", "C6", "C10", "C11", "C9", "C7", "C8", "C12", "C13", "C15", "C16"))
p5.2 <- ggplot(cM, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill", width = 1) + ArchR::theme_ArchR() + scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(values = sample_colors) + labs(x = "Sample", y = "%Cluster") +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))

pdf("output/Plots/02_Barplot_ClusterSample.pdf", height = 4, width = 8)
p4
p5
p5.2
dev.off()

#----- Cell-type marker gene activity visualization -----#
arc <- addImputeWeights(arc, reducedDims = "IterativeLSI")
markers <- c("EPCAM", "KRT7", "KRT8", "KRT15", "KRT16", "KRT18", "KRT19", "KRT5", "KRT14", "KRT17", 
             "CDH1", "CDH2", "CD44", "MKI67", "MYC", "EGFR", "TP53", "CCNB1", "CCND1", "ESR1", "FOXA1", 
             "GRHL2", "AR", "ERBB2", "VIM", "PTPRC", "CD19", "MS4A1", "PAX5", "IGLL5",
             "SDC1", "CD38", "CD3D", "BCL11B", "CD4", "CD8A", "CD14", "ITGAM", "ITGAX", "CD68", 
             "PECAM1", "KDR", "ENG", "FLT4", "THY1", "ACTA2", "CXCL12", "FAP", "NCAM1", "KLRD1", 
             "NCR1", "FOXP3", "PDCD1", "CTLA4", "TBX21", "EOMES")

### UMAP overlay ###
p6 <- plotEmbedding(arc, 
                    colorBy = "GeneScoreMatrix", 
                    name = markers, 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", 
                    size = 0.5, rastr = T)
plotPDF(p6, name = "02_UMAPOL_MarkerGA.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

### Cluster-specific activated genes ###
GAClusters <- getMarkerFeatures(arc, 
                                useMatrix = "GeneScoreMatrix", 
                                groupBy = "Clusters",
                                bias = c("TSSEnrichment", "log10(nFrags)"),
                                testMethod = "wilcoxon")
saveRDS(GAClusters, "rds/02_GAClusters.rds")
MarkerClusters <- getMarkers(seMarker = GAClusters, cutOff = "FDR < 0.01 & Log2FC >= 1") #C22: No genes

MarkerClustersDF <- lapply(paste0("C", c(1:21,23:33)), function(x){
  out <- MarkerClusters[[x]]
  out$Clusters <- x
  return(out)
})
MarkerClustersDF <- do.call(rbind, MarkerClustersDF)
write.csv(MarkerClustersDF, "output/Tables/02_MarkerClusterGA.csv")
MarkerClusters_SYMBOL <- lapply(MarkerClusters, function(x) x$name)
lapply(names(MarkerClusters_SYMBOL), function(x) write.table(data.frame(gene = MarkerClusters_SYMBOL[[x]]), paste0("output/Tables/02_MarkerClusterSYMBOL/02_MarkerClusterSYMBOL_", x, ".txt"), col.names = F, row.names = F, quote = F))


### gene enrichment analysis ###
source("code/enrichment_func.R")
gmt_files <- paste0("ref/CellTypeDB/", list.files("ref/CellTypeDB/", pattern = ".gmt"))
DIR <- "output/clusterProfiler_Celltype"

for (i in 1:length(gmt_files)){
  skip_to_next <- FALSE
  gmt_file <- gmt_files[i]
  tryCatch(EnrichmentAnalysis(gmt_file, MarkerClusters_SYMBOL, DIR),
           error = function(e) {skip_to_next <<- TRUE})
  if(skip_to_next) { next } 
}

saveRDS(arc, "rds/02_arc.rds")

#----- genomeTrack -----#
markers3 <- c("EPCAM", "CDH1", "ESR1", "FOXA1", "ERBB2", "VIM", "KDR", "THY1", "CD3D", "PAX5", "IGLL5", "ITGAX")

gr <- ArchR::geneAnnoHg19$genes[ArchR::geneAnnoHg19$genes$symbol %in% markers3]
names(gr) <- gr$symbol
gr <- gr[markers3]
gr <- extendGR(gr, upstream = 5000, downstream = 5000)

final_anno <- factor(final_anno, names(celltype_colors))
idx3 <- sort(final_anno) %>% names

p12 <- plotBrowserTrack(arc, groupBy = "Clusters", region = gr, 
                        plotSummary = c("bulkTrack", "geneTrack"),
                        useGroups = idx3)
plotPDF(p12, name = "02_GenomeTrack_Marker.pdf", ArchRProj = arc, addDOC = FALSE, width = 6, height = 8)

#----- Sample x cell type -----#
rownames(pt_ls) <- pt_ls$ID

tmp <- as.matrix(confusionMatrix(arc$Sample, arc$CellTypeFinal))
tmp <- tmp / rowSums(tmp)
tmp <- tmp[order(tmp[,"Epithelial"], decreasing = F),] %>% as.data.frame()
tmp$SampleType <- pt_ls[rownames(tmp),]$Type
tmp$Receptor <- pt_ls[rownames(tmp),]$receptor_status
tmp$ID <- rownames(tmp)

idy1 <- pt_ls$ID[which(pt_ls$Type %in% c("Primary", "LocalRecurrence") & pt_ls$receptor_status %in% c("ER+/HER2-", "ER-/PGR+/HER-"))]
idy2 <- pt_ls$ID[which(pt_ls$Type %in% c("Primary", "LocalRecurrence") & pt_ls$receptor_status %in% c("HER2+", "ER+/HER2+"))]
idy3 <- pt_ls$ID[which(pt_ls$Type %in% c("Primary", "LocalRecurrence") & pt_ls$receptor_status %in% "TN")]
idy4 <- pt_ls$ID[pt_ls$Type %ni% c("Primary", "LocalRecurrence")]

idy5 <- c(rownames(tmp)[rownames(tmp) %in% idy1],
          rownames(tmp)[rownames(tmp) %in% idy2],
          rownames(tmp)[rownames(tmp) %in% idy3],
          rownames(tmp)[rownames(tmp) %in% idy4])

cM4 <- reshape2::melt(tmp)
cM4$ID <- factor(cM4$ID, levels = idy5)
cM4$variable <- factor(cM4$variable, levels = c("Epithelial", "T", "B", "Plasma", "Macrophage","Fibroblast", "Endothelial"))

p13 <- ggplot(cM4, aes(x = ID, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") + ArchR::theme_ArchR() + 
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  labs(x = "Cluster", y = "%Sample") + scale_fill_manual(values = celltype_colors) +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))

cM4$SampleType[cM4$SampleType %in% c("Pleural", "Ascitis")] <- "Pleural/Ascitis"
cM4 <- cM4[cM4$SampleType %ni% "LocalRecurrence",]
cM4$SampleType <- factor(cM4$SampleType, levels = c("Primary", "LymphNode", "Pleural/Ascitis"))

p14 <- ggplot(cM4, aes(x = SampleType, y = value*100, fill = variable)) +
  geom_boxplot(outlier.shape = 3, alpha = 0.8) + ArchR::theme_ArchR() + 
  labs(x = "Cluster", y = "%Sample") + scale_fill_manual(values = celltype_colors) +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))

cM5 <- cM4[which(cM4$SampleType == "Primary" & cM4$Receptor %in% c("ER+/HER2-", "HER2+", "TN")),]
p15 <- ggplot(cM5, aes(x = Receptor, y = value*100, fill = variable)) +
  geom_boxplot(outlier.shape = 3, alpha = 0.8) + ArchR::theme_ArchR() + 
  labs(x = "Cluster", y = "%Sample") + scale_fill_manual(values = celltype_colors) +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))

pdf("output/Plots/02_BarBoxplot_CelltypeComposition_v2.pdf", height = 4, width = 6)
p13
p14
p15
dev.off()
