#----------------------------------------------------------------------------
# 04_CNA.R
#----------------------------------------------------------------------------
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ComplexHeatmap)
library(circlize)
source("code/func_scATAC_CNV.R")
addArchRThreads(threads = 22)
addArchRGenome("hg19")

arc <- readRDS("rds/02_arc.rds")
arc_ep <- readRDS("rds/03_arc_ep.rds")

#----- CNV estimation -----#
fragments <- getFragmentsFromProject(ArchRProj = arc) #per sample
blacklist <- import.bed("ref/hg19-blacklist.v2.bed")

#analysis
windows <- makeWindows(genome = BSgenome.Hsapiens.UCSC.hg19, blacklist = blacklist)
cna_samples <- mclapply(fragments, function(x){
  out <- scCNA(windows, x, neighbors = 100, LFC = 1.5, FDR = 0.1, force = TRUE, remove = c("chrM","chrY"))
  return(out)
}, mc.cores = 12)

saveRDS(cna_samples, "rds/04_cna_samples.rds")

#----- extract cna info -----#
#matrix
mtx_cna <- lapply(cna_samples, function(x){
  out <- t(assays(x)$CNA) %>% `rownames<-`(., rownames(x@colData))
  return(out)
})
mtx_cna <- do.call(rbind, mtx_cna)
mtx_cna <- mtx_cna[arc$cellNames,]

window_info <- rowRanges(cna_samples[[1]])
window_info$name2 <- paste0("wd_", c(1:1074))
names(window_info) <- window_info$name2
colnames(mtx_cna) <- window_info$name2
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

mtx_cna_EP <- mtx_cna[which(arc$cellType == "Epithelial"),]
mtx_cna_nEP <- mtx_cna[which(arc$cellType != "Epithelial"),]

#annotation
cellAnno <- data.frame(Cell = arc$cellNames, Sample = arc$Sample, Cluster = arc$Clusters)

col_fun1 = colorRamp2(c(0,1), c("#E6E7E8", "#8816A7"))

sample_colors <- readRDS("rds/02_sample_colors.rds")
cluster_colors <- readRDS("rds/02_cluster_colors.rds")
chr <- as.factor(seqnames(cna_samples[[1]]@rowRanges))
ha1 <- HeatmapAnnotation(chr = chr, col = list(chr = c(ArchRPalettes$stallion2, ArchRPalettes$calm)[c(1:23)] %>% `names<-`(., levels(chr))))
ra1 <- rowAnnotation(Cluster = cellAnno$Cluster, Sample = cellAnno$Sample,
                     col = list(Cluster = cluster_colors,
                                Sample = sample_colors))
samples <- names(sample_colors)

#----- Non-EP CNA -----#
cluster_rowsplit <- arc$Clusters[which(arc$cellType != "Epithelial")]
ha2 <- HeatmapAnnotation(frequency = anno_points(colSums(mtx_cna_nEP)/nrow(mtx_cna_nEP)))
ht1 <- Heatmap(mtx_cna_nEP, name = "CNA estimation", 
               cluster_rows = F, cluster_columns = F, 
               show_row_dend = F, show_row_names = F, show_column_names = F,
               col = col_fun1, 
               row_split = factor(cluster_rowsplit, levels = paste0("C", 1:39)),
               top_annotation = c(ha1,ha2), 
               left_annotation = ra1[which(arc$cellType != "Epithelial")], 
               use_raster = T)
p1 <- draw(ht1)
pdf("output/Plots/04_Heatmap_CNV_NonEP.pdf", width = 8, height = 8)
p1
dev.off()

#----- EP CNA all windows -----#
cluster_rowsplit <- arc$Clusters[which(arc$cellType == "Epithelial")]
ha3 <- HeatmapAnnotation(frequency = anno_points(colSums(mtx_cna_EP)/nrow(mtx_cna_EP)))
ht2 <- Heatmap(mtx_cna_EP, name = "CNA estimation", 
               cluster_rows = F, cluster_columns = F, 
               show_row_dend = F, show_row_names = F, show_column_names = F,
               col = col_fun1, 
               row_split = factor(cluster_rowsplit, levels = paste0("C", 1:39)),
               top_annotation = c(ha1,ha3), 
               left_annotation = ra1[which(arc$cellType == "Epithelial")], 
               use_raster = T)
p2 <- draw(ht2)
pdf("output/Plots/04_Heatmap_CNV_EP_allwindow.pdf", width = 8, height = 8)
p2
dev.off()

#----- non EP CNA frequent windows -----#
tmp1 <- colSums(mtx_cna_nEP)/nrow(mtx_cna_nEP)
tmp2 <- colSums(mtx_cna_EP)/nrow(mtx_cna_EP)
plot(tmp1)
which(tmp1 > 0.1)
which(tmp2 > 0.1)
 
#----- EP CNA remove frequent windows in non EP -----#
cluster_rowsplit <- arc$Clusters[which(arc$cellType == "Epithelial")]
mtx_cna_EP2 <- mtx_cna_EP[,which(tmp1 <= 0.1)]
ncol(mtx_cna_EP2)
ha4 <- HeatmapAnnotation(frequency = anno_points(colSums(mtx_cna_EP2)/nrow(mtx_cna_EP2)))
ht3 <- Heatmap(mtx_cna_EP2, name = "CNA estimation", 
               cluster_rows = F, cluster_columns = F, 
               show_row_dend = F, show_row_names = F, show_column_names = F,
               col = col_fun1, 
               row_split = factor(cluster_rowsplit, levels = paste0("C", 1:39)),
               top_annotation = c(ha1[which(tmp1 <= 0.1)],ha4), 
               left_annotation = ra1[which(arc$cellType == "Epithelial")], 
               use_raster = T)
p3 <- draw(ht3)
pdf("output/Plots/04_Heatmap_CNV_EP_selectedWindows.pdf", width = 8, height = 8)
p3
dev.off()

#epithelial clusters
epicluster_colors <- readRDS("rds/03_epicluster_colors.rds")
cellAnno2 <- data.frame(Cell = arc_ep$cellNames, Sample = arc_ep$Sample, Cluster = arc_ep$Clusters)
ra2 <- rowAnnotation(Cluster = cellAnno2$Cluster, Sample = cellAnno2$Sample,
                     col = list(Cluster = epicluster_colors,
                                Sample = sample_colors))
mtx_cna_EP2 <- mtx_cna_EP2[cellAnno2$Cell,]

cluster_rowsplit <- arc_ep$Clusters
ha5 <- HeatmapAnnotation(frequency = anno_points(colSums(mtx_cna_EP2)/nrow(mtx_cna_EP2)))
ht4 <- Heatmap(mtx_cna_EP2, name = "CNA estimation", 
               cluster_rows = F, cluster_columns = F, 
               show_row_dend = F, show_row_names = F, show_column_names = F,
               col = col_fun1, 
               row_split = factor(cluster_rowsplit, levels = paste0("C", 1:35)),
               top_annotation = c(ha1[which(tmp1 <= 0.1)],ha5), 
               left_annotation = ra2, 
               use_raster = T)
p4 <- draw(ht4)
pdf("output/Plots/04_Heatmap_CNV_EP_selectedWindows_EPClust.pdf", width = 8, height = 10)
p4
dev.off()

#----- EP CNA frequent amplified region -----#
#calculate % of CNA in each cluster
percentCNA_cluster <- lapply(paste0("C", c(1:35)), function(i){
  idx <- cellAnno2$Cell[which(cellAnno2$Cluster == i)]
  mtx <- mtx_cna_EP2[idx,]
  out <- colSums(mtx)/nrow(mtx)
  return(out)
})
percentCNA_cluster <- do.call(rbind, percentCNA_cluster)
percentCNA_cluster_max <- apply(percentCNA_cluster,2,max)
percentCNA_cluster_mean <- apply(percentCNA_cluster,2,mean)

df <- data.frame(Max = percentCNA_cluster_max, Mean = percentCNA_cluster_mean)
p5 <- ggplot(df, aes(x = Mean, y = Max)) + geom_point() + theme_ArchR() + 
  labs(x = "Mean", y = "Max") +
  geom_hline(yintercept = 0.9) + geom_vline(xintercept = 0.1)
pdf("output/Plots/04_ScatterPlot_CNV_EP_MAX_Mean.pdf", width = 5, height = 5)
p5
dev.off()

pAR_specific <- window_info[rownames(df)[which(df$Max > 0.9 & df$Mean <= 0.1)]]
pAR_common <- window_info[rownames(df)[which(df$Mean > 0.1)]]

export.bed(pAR_specific, "output/output_bed/04_pAR_specific.bed")
export.bed(pAR_common, "output/output_bed/04_pAR_common.bed")

#common
idx <- lapply(seq_along(pAR_common), function(i) findOverlaps(pAR_common[i], pAR_common) %>% subjectHits)
idx <- idx[unlist(lapply(idx, length)) > 1] %>% unique()
idx <- idx[c(1,2,4)]
pAR_common_merge <- lapply(idx, function(i) GenomicRanges::reduce(pAR_common[i])) 
pAR_common_merge <- unlist(GRangesList(pAR_common_merge))

#pAR_specific
idx <- lapply(seq_along(pAR_specific), function(i) findOverlaps(pAR_specific[i], pAR_specific) %>% subjectHits)
idx <- idx[unlist(lapply(idx, length)) > 1] %>% unique()
pAR_specific_merge <- lapply(idx, function(i) GenomicRanges::reduce(pAR_specific[i])) 
pAR_specific_merge <- unlist(GRangesList(pAR_specific_merge))

export.bed(pAR_specific_merge, "output/output_bed/04_pAR_specific_merge.bed")
export.bed(pAR_common_merge, "output/output_bed/04_pAR_ccommon_merge.bed")

#----- Karyotype plot -----#
pAR_specific_merge <- import.bed("output/output_bed/04_pAR_specific_merge.bed")
pAR_common_merge <- import.bed("output/output_bed/04_pAR_ccommon_merge.bed")

library(karyoploteR)

kp <- plotKaryotype(genome = "hg19", chromosomes=c("chr1","chr8"))
kpPlotRegions(kp, data=pAR_common_merge)

pdf("output/Plots/04_CNA_karyotypePlot.pdf",width = 6, height = 8)
kp <- plotKaryotype(genome = "hg19", chromosomes=c("chr1","chr5","chr8","chr12","chr17","chr19"))
kpPlotRegions(kp, data=pAR_common_merge, col="darkblue", layer.margin = 0.01, border=NA, r0=0, r1=0.5)
kpPlotRegions(kp, data=pAR_specific_merge, col="red", layer.margin = 0.05, border=NA, r0=0.6, r1=1)
dev.off()

#----- ERBB2 -----#
mtx_cna_log <- lapply(cna_samples, function(x){
  out <- t(assays(x)$log2FC) %>% `rownames<-`(., rownames(x@colData))
  return(out)
})
mtx_cna_log <- do.call(rbind, mtx_cna_log)
mtx_cna_log <- mtx_cna_log[arc_ep$cellNames,]

colnames(mtx_cna_log) <- window_info$name2
idy <- findOverlaps(window_info, geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == "ERBB2")]) %>% queryHits()
erbb2_log2FC <- data.frame(apply(mtx_cna_log[,idy],1,mean))

erbb2_log2FC$Clusters <- gsub("C", "EP", arc_ep$Clusters)
colnames(erbb2_log2FC) <- c("log2FC", "Clusters")
p9 <-  ggplot(erbb2_log2FC, aes(x = factor(Clusters, levels = paste0("EP", c(1:35))), y = log2FC, fill = Clusters)) + geom_boxplot() + ylim(-3.5,3.5) +
  theme_ArchR() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = epicluster_colors %>% `names<-`(.,paste0("EP", c(1:35))))
pdf("output/Plots/04_Boxplot_ERBB2_window_log2FC.pdf", width = 8, height = 6)
p9
dev.off()
