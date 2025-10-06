#----------------------------------------------------------------------------
# 04_clustering_epi.R
#----------------------------------------------------------------------------
library(ArchR)
library(scales)
library(ggsignif)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(ggrastr)
library(karyoploteR)
library(ggrepel)
addArchRThreads(threads = 24)
addArchRGenome("hg19")

arc <- readRDS("rds/02_arc.rds")
pt_ls <- read.csv("ref/patient_list_v2.csv")
SampleTypeColors <- c("Primary" = "#666666", "LocalRecurrence" = "#1B9E77", "LymphNode" = "#E7298A", "Pleural" = "#7570B3", "Ascitis" = "#E6AB02")

#----- Clustering epithelial cells -----#
arc_epi <- arc[arc$CellTypeFinal == "Epithelial"]

arc_epi <- addIterativeLSI(arc_epi, useMatrix = "TileMatrix", name = "IterativeLSI", force = T)
arc_epi <- addUMAP(arc_epi, reducedDims = "IterativeLSI", force = T)
arc_epi <- addClusters(arc_epi, reducedDims = "IterativeLSI", maxClusters = 40, resolution = 1.2, force = T)

sample_colors <- readRDS("rds/02_sample_colors.rds")

p1 <- plotEmbedding(arc_epi, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", plotAs = "points", size = 0.5, rastr = T)
p2 <- plotEmbedding(arc_epi, colorBy = "cellColData", name = "Sample", embedding = "UMAP", plotAs = "points", size = 0.5, pal = sample_colors, labelMeans = F, rastr = T)
p3 <- plotEmbedding(arc_epi, colorBy = "cellColData", name = "SampleType", embedding = "UMAP", plotAs = "points", pal = SampleTypeColors, size = 0.5, labelMeans = F, rastr = T)
plotPDF(p1,p2,p3, name = "04_UMAP_epi.pdf", ArchRProj = arc_epi, addDOC = FALSE, width = 5, height = 5)

arc_epi$Clusters <- gsub("C", "Ep", arc_epi$Clusters)

cM <- confusionMatrix(arc_epi$Sample, arc_epi$Clusters)[pt_ls$ID, paste0("Ep", c(1:35))] %>% as.matrix()
write.csv(cM, "output/Tables/04_ClusterSampleMatrix_Ep.csv")

cluster_colors <- p1$plot_env$pal
names(cluster_colors) <- str_split(names(cluster_colors), pattern = "-", simplify = T)[,2]
names(cluster_colors) <- gsub("C", "Ep", names(cluster_colors))

p4 <- ggplot(reshape2::melt(cM), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") + ArchR::theme_ArchR() + scale_y_continuous(labels = percent) +
  scale_fill_manual(values = cluster_colors) + labs(x = "Sample", y = "%Cluster") +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))
p5 <- ggplot(reshape2::melt(cM), aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") + ArchR::theme_ArchR() + scale_y_continuous(labels = percent) +
  scale_fill_manual(values = sample_colors) + labs(x = "Cluster", y = "%Sample") +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))

pdf("output/Plots/04_Barplot_ClusterSample_Ep.pdf", height = 4, width = 8)
p4
p5
dev.off()

saveRDS(arc_epi, "rds/04_arc_epi.rds")

#----- Normal vs Tumor cells (epiAneufinder results) -----#
# arc_epi <- readRDS("rds/04_arc_epi.rds")
epiAneufinder_se_ls <- lapply(pt_ls$ID, function(i){
  dt <- data.table::fread(paste0("/mnt/host_mnt/Volumes/Shared/Kume/Analysis/S23-3_scATAC-BRCA/output/AllSample1Mb/",
                                 i, "/epiAneufinder_results/results_table.tsv")) %>% data.frame
  gr <- GRanges(seqnames = dt$seq, IRanges(start=dt$start, end=dt$end))
  dt <- dt[,c(5:ncol(dt))]
  colnames(dt) <- gsub("cell.", "", colnames(dt))
  colnames(dt) <- paste0(i, "#", colnames(dt))
  
  dt <- dt[, intersect(colnames(dt), arc$cellNames)]
  se <- SummarizedExperiment(assays = list(CNV = dt), rowRanges = gr)
  
  return(se)
})
names(epiAneufinder_se_ls) <- pt_ls$ID

### Non-specific amplified regions determined by non-epithelial cells ###
idx1 <- arc[arc$CellTypeFinal != "Epithelial"]$cellNames

# % of amplified cells per each window
prop_NonEpi_ampCell <- lapply(epiAneufinder_se_ls[-2], function(x){
  mtx <- assay(x)[,which(colnames(x) %in% idx1)]
  out <- apply(mtx, 1, function(i) length(which(i == 2)))
  out <- out / ncol(mtx)
  return(out)
})

# amp windows (>25% of cells with amplified "2")
gr_NonEpi_ampWindow <- lapply(names(epiAneufinder_se_ls)[-2], function(i){
  out <- rowRanges(epiAneufinder_se_ls[[i]])[which(prop_NonEpi_ampCell[[i]] > 0.25)]
  if(length(out) > 0){
    out$Sample <- i
  }
  return(out)
})
gr_NonEpi_ampWindow <- unlist(GRangesList(gr_NonEpi_ampWindow))
gr_NonEpi_ampWindow <- sort(gr_NonEpi_ampWindow)
names(gr_NonEpi_ampWindow) <- NULL

# redundunt amp window (at least 2 samples)
gr_NonEpi_ampWindow_redun <- gr_NonEpi_ampWindow[which(countOverlaps(gr_NonEpi_ampWindow,gr_NonEpi_ampWindow) >= 2)] %>% unique

# subsetByOverlaps(geneAnnoHg19$genes, gr_NonEpi_ampWindow_redun)$symbol %>% as.character() %>% sort

gr1 <- subsetByOverlaps(gr_NonEpi_ampWindow, gr_NonEpi_ampWindow_redun, invert = T)

pdf("output/Plots/04_Karyotype_NormalAmpWindows.pdf", width = 6, height = 6)
kp <- plotKaryotype(genome = "hg19", chromosomes = paste0("chr", c(1:22, "X")))
kpPlotRegions(kp, data=gr1)
kpPlotRegions(kp, data=gr_NonEpi_ampWindow_redun, col = "red")
dev.off()

#epithelial clusters
idx2 <- arc[arc$CellTypeFinal == "Epithelial"]$cellNames
prop_Epi_ampWindow <- lapply(epiAneufinder_se_ls, function(x){
  mtx <- assay(x)[,which(colnames(x) %in% idx2)]
  out <- apply(mtx, 1, function(i) length(which(i == 2)))
  out <- out / ncol(mtx)
  return(out)
})
gr_Epi_ampWindow <- lapply(names(epiAneufinder_se_ls), function(i){
  out <- rowRanges(epiAneufinder_se_ls[[i]])[which(prop_Epi_ampWindow[[i]] > 0.25)]
  if(length(out) > 0){
    out$Sample <- i
  }
  return(out)
})

# amplified windows in 25% of epithelial cells per each sample
gr_Epi_ampWindow <- unlist(GRangesList(gr_Epi_ampWindow)) #2002 regions
# gr_Epi_ampWindow_redun <- gr_Epi_ampWindow[which(countOverlaps(gr_Epi_ampWindow,gr_Epi_ampWindow) >= 5)] %>% unique #359 region
gr_Epi_ampWindow <- subsetByOverlaps(gr_Epi_ampWindow, gr_NonEpi_ampWindow_redun, invert = T) #1834 region

pdf("output/Plots/04_Karyotype_TumorAmpWindows.pdf", width = 6, height = 8)
kp <- plotKaryotype(genome = "hg19", chromosomes = paste0("chr", c(1:22, "X")))
kpPlotRegions(kp, data=gr_Epi_ampWindow)
dev.off()

gr_Epi_ampWindowER <- gr_Epi_ampWindow[gr_Epi_ampWindow$Sample %in% pt_ls$ID[which(pt_ls$receptor_status == "ER+/HER2-")]]
gr_Epi_ampWindowER_redun <- gr_Epi_ampWindowER[which(countOverlaps(gr_Epi_ampWindowER,gr_Epi_ampWindowER) >= 7)] %>% unique

gr_Epi_ampWindowTN <- gr_Epi_ampWindow[gr_Epi_ampWindow$Sample %in% pt_ls$ID[which(pt_ls$receptor_status == "TN")]]
gr_Epi_ampWindowTN_redun <- gr_Epi_ampWindowTN[which(countOverlaps(gr_Epi_ampWindowTN,gr_Epi_ampWindowTN) >= 2)] %>% unique

pdf("output/Plots/04_Karyotype_TumorAmpWindowsER.pdf", width = 6, height = 8)
kp <- plotKaryotype(genome = "hg19", chromosomes = paste0("chr", c(1:22, "X")))
kpPlotRegions(kp, data=gr_Epi_ampWindowER)
kpPlotRegions(kp, data=gr_Epi_ampWindowER_redun, col = "red")
dev.off()

pdf("output/Plots/04_Karyotype_TumorAmpWindowsTN.pdf", width = 6, height = 8)
kp <- plotKaryotype(genome = "hg19", chromosomes = paste0("chr", c(1:22, "X")))
kpPlotRegions(kp, data=gr_Epi_ampWindowTN)
kpPlotRegions(kp, data=gr_Epi_ampWindowTN_redun, col = "red")
dev.off()

saveRDS(gr_Epi_ampWindow, "rds/04_gr_Epi_ampWindow.rds")

options(scipen = 10000)
d <- data.frame(seqnames = seqnames(gr_Epi_ampWindowER_redun), start = start(gr_Epi_ampWindowER_redun)-1, end = end(gr_Epi_ampWindowER_redun))
write.table(d, "output/output_bed/04_AmpRegions_ER_redun.bed", row.names = F, col.names = F, quote = F, sep = "\t")
d <- data.frame(seqnames = seqnames(gr_Epi_ampWindowTN_redun), start = start(gr_Epi_ampWindowTN_redun)-1, end = end(gr_Epi_ampWindowTN_redun))
write.table(d, "output/output_bed/04_AmpRegions_TN_redun.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# calculating CNV score: amplified windows / all tumor amplified windows
prop_Epi_ampCell <- lapply(epiAneufinder_se_ls, function(x){
  qh1 <- queryHits(findOverlaps(rowRanges(x), gr_Epi_ampWindow))
  mtx <- assay(x)[qh1,which(colnames(x) %in% idx2)]
  out <- apply(mtx, 2, function(i) length(which(i == 2)))
  out <- out / nrow(mtx)
  return(out)
})
names(prop_Epi_ampCell) <- NULL
prop_Epi_ampCell <- unlist(prop_Epi_ampCell)
plot(sort(prop_Epi_ampCell))

# visualize CNV score
arc_epi$CNVscore <- prop_Epi_ampCell[arc_epi$cellNames]
arc_epi$CNVscore[is.na(arc_epi$CNVscore)] <- 0

df1 <- arc_epi@embeddings$UMAP$df
colnames(df1) <- c("UMAP1", "UMAP2")
df1$CNVscore <- arc_epi$CNVscore
df1 <- df1[order(df1$CNVscore),]

p6 <- ggplot(df1, aes(x = UMAP1, y = UMAP2, color = CNVscore)) + geom_point_rast(size = 0.5) + scale_color_viridis_c(option = "G") + theme_ArchR()
pdf("output/Plots/04_UMAPOL_CNVscore.pdf", width = 5, height = 5.5)
p6
dev.off()

df2 <- data.frame(Clusters = arc_epi$Clusters, CNV = arc_epi$CNVscore)
df2$Clusters <- factor(df2$Clusters, levels = paste0("Ep", 1:35))

p7 <- ggplot(df2, aes(x = Clusters, y = CNV)) + geom_boxplot(outlier.shape = 3) + theme_ArchR() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "Clusters", y = "CNV score")
idx3 <- sapply(paste0("Ep", 1:35), function(x) median(df2$CNV[df2$Clusters==x])) %>% sort %>% names
p7.2 <- ggplot(df2, aes(x = factor(Clusters, levels = idx3),y = CNV)) + geom_boxplot(outlier.alpha = .25) + theme_ArchR() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "Clusters", y = "CNV score")
pdf("output/Plots/04_Boxplot_CNVscore_Clusters.pdf", width = 8, height = 4)
p7
p7.2
dev.off()

cM2 <- t(t(cM) / colSums(cM))*100
cM2_value <- apply(cM2, 2, max)
MeanCNVscore <- lapply(paste0("Ep",c(1:35)), function(i) mean(df2$CNV[df2$Clusters == i])) %>% unlist
df3 <- data.frame(MaxPercSample = cM2_value, MeanCNVscore = MeanCNVscore, Clusters = names(cM2_value))

p8 <- ggplot(df3, aes(x=MeanCNVscore,y=MaxPercSample,label=Clusters)) + geom_point(shape = 4) + theme_ArchR() +
  ylim(0,100) + geom_text_repel() + labs(x = "Mean CNV score", y = "%Cell of the most frequent sample")
pdf("output/Plots/04_Scatter_CNVscore_MostFreqSample.pdf", width = 6, height = 6)
p8
dev.off()
write.csv(df3, "output/Tables/04_Scatter_CNVscore_MostFreqSample.csv")

df3 <- read.csv("output/Tables/04_Scatter_CNVscore_MostFreqSample.csv")


#----- gene activity heatmap -----#
GeneScoreClusters <- getMarkerFeatures(arc_epi, 
                                       useMatrix = "GeneScoreMatrix", 
                                       groupBy = "Clusters",
                                       bias = c("TSSEnrichment", "log10(nFrags)"),
                                       testMethod = "wilcoxon")

mtx <- assays(GeneScoreClusters)$Mean
rownames(mtx) <- rowData(GeneScoreClusters)$name
mtx <- mtx[c("KRT18","KRT19","KRT5","KRT14","CDH1","CDH2","MYC","EGFR","CCND1","FOXA1","AR","VIM","ESR1","PGR","ERBB2","GRB7"),]
mtx <- t(scale(t(mtx)))

col_fun2 <- colorRamp2(c(-2,-1,0,1,2), c(viridis(5, option = "A")))
fh <- function(x) hclust(dist(x), method = "ward.D2")
ht1 <- Heatmap(t(mtx), name = "Mean GeneScore z-score", cluster_columns = fh, cluster_rows = fh, col = col_fun2, border = "black",
               rect_gp = gpar(col = "gray", lwd = 0.5))
p9 <- draw(ht1)
pdf("output/Plots/04_Heatmap_markerGS_Epi.pdf", width = 6, height = 7)
p9
dev.off()

#----- gene activity differential -----#
arc_epi <- readRDS("rds/04_arc_epi.rds")
ClusterSubtype <- list(HER2 = c(3,4), ER = c(27,23,22,29,34,35,32,30,12,16,28,21,33,15,20,14,25), TN = c(5,6,1,2))
ClusterSubtype <- lapply(ClusterSubtype, function(i) paste0("Ep",i))

arc_epi$ClusterSubtype <- "Others"
tmp1 <- lapply(names(ClusterSubtype), function(x){
  tmp <- ClusterSubtype[[x]]
  out <- which(arc_epi$Clusters %in% tmp)
  return(out)
})
names(tmp1) <- names(ClusterSubtype)
arc_epi$ClusterSubtype[tmp1$HER2] <- "HER2"
arc_epi$ClusterSubtype[tmp1$ER] <- "ER"
arc_epi$ClusterSubtype[tmp1$TN] <- "TN"

GeneScoreClusters2 <- getMarkerFeatures(arc_epi, 
                                       useMatrix = "GeneScoreMatrix", 
                                       groupBy = "ClusterSubtype",
                                       bias = c("TSSEnrichment", "log10(nFrags)"),
                                       testMethod = "wilcoxon")
MarkersClusters2 <- getMarkers(GeneScoreClusters2, cutOff = "FDR < 0.05 & Log2FC > 1")

MarkerClustersDF <- lapply(c("ER", "HER2", "TN"), function(x){
  out <- MarkersClusters2[[x]]
  out$Clusters <- x
  return(out)
})
MarkerClustersDF <- do.call(rbind, MarkerClustersDF)
write.csv(MarkerClustersDF, "output/Tables/04_MarkerSubtypeGA.csv")

p10 <- plotMarkers(GeneScoreClusters2, name = "ER", cutOff = "FDR < 0.05 & Log2FC > 1")
p11 <- plotMarkers(GeneScoreClusters2, name = "HER2", cutOff = "FDR < 0.05 & Log2FC > 1")
p12 <- plotMarkers(GeneScoreClusters2, name = "TN", cutOff = "FDR < 0.05 & Log2FC > 1")
plotPDF(p10,p11,p12, name = "04_Volcano_SubtypeGA.pdf", ArchRProj = arc_epi, addDOC = FALSE, width = 5, height = 5)

