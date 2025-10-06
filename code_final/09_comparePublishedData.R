#------------------------------------------------------------------------------
# 09_comparePublishedData.R
#------------------------------------------------------------------------------
library(ArchR)
library(rtracklayer)
library(ggplot2)
library(scales)
library(ComplexHeatmap)
library(circlize)
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}

#----- liftover hg19 to hg38 -----#
arc <- readRDS("rds/06_arc.rds")
gr_valid <- readRDS("rds/05_gr_valid.rds")

gr <- arc@peakSet
gr <- subsetByOverlaps(gr, gr_valid)
names(gr) <- c(paste0("cCRE_", c(1:length(gr))))

d <- data.frame(seqnames = seqnames(gr), start = start(gr)-1, end = end(gr), name = names(gr))
write.table(d, "output/output_bed/09_peakSet_hg19.bed", row.names = F, col.names = F, quote = F, sep = "\t")

ClusterPeaksList <- readRDS("rds/06_ClusterPeaksList.rds")

ClusterPeaksList <- lapply(names(ClusterPeaksList), function(x){
  gr_out <- ClusterPeaksList[[x]]
  gr_out$name <- paste0(x, "_DA", c(1:length(gr_out)))
  return(gr_out)
})
gr_DA <- unlist(GRangesList(ClusterPeaksList))
d <- data.frame(seqnames = seqnames(gr_DA), start = start(gr_DA)-1, end = end(gr_DA), name = gr_DA$name)
write.table(d, "output/output_bed/09_DACRE_hg19.bed", row.names = F, col.names = F, quote = F, sep = "\t")

#----- hg38 liftovered regions -----#
gr2 <- import.bed("output/output_bed/09_peakSet_hg38.bed")
gr_DA2 <- import.bed("output/output_bed/09_DACRE_hg38.bed")
gr_DA2$cluster <- stringr::str_split(gr_DA2$name, "_", simplify = T)[,1]

gr_DA2_ls <- lapply(unique(gr_DA2$cluster), function(x) gr_DA2[gr_DA2$cluster==x])
names(gr_DA2_ls) <- unique(gr_DA2$cluster)

#----- Peak overlaps, all peaks -----#
arc_T <- readRDS("../A2501_EpiEncode/Terekhanova_v1/rds/02_arc.rds")
gr_T <- arc_T@peakSet
arc_R <- readRDS("../A2501_EpiEncode/regner_v1//rds/01_arc.rds")
gr_R <- arc_R@peakSet

#proportion of overlaps
fo1 <- findOverlaps(gr2, gr_T)
length(unique(queryHits(fo1)))
fo2 <- findOverlaps(gr2, gr_R)
length(unique(queryHits(fo2)))

df1 <- data.frame(data = factor(c("Terekhanova", "Terekhanova", "Regner", "Regner"), levels = c("Terekhanova", "Regner")),
                  overlap = rep(c("overlap", "Nonoverlap"),2),
                  peaks = c(length(unique(queryHits(fo1))), length(gr2) - length(unique(queryHits(fo1))),
                            length(unique(queryHits(fo2))), length(gr2) - length(unique(queryHits(fo2)))))
# data    overlap  peaks
# 1 Terekhanova    overlap 189188
# 2 Terekhanova Nonoverlap  35348
# 3      Regner    overlap 175766
# 4      Regner Nonoverlap  48770
p1 <- ggplot(df1, aes(x=data,y=peaks,fill=overlap)) + geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent, expand = c(0,0)) + theme_cb() + scale_fill_manual(values = c("overlap"="orange","Nonoverlap"="gray"))

pdf("output/Plots/09_Barplot_overlapPeaks_acrossDataSet.pdf", width = 4, height = 4)
p1
dev.off()

#output for Table
gr2_out <- gr2
gr2_out$Terekhanova <- "NonOverlap"
gr2_out$Terekhanova[queryHits(fo1)] <- "Overlap"
gr2_out$Regner <- "NonOverlap"
gr2_out$Regner[queryHits(fo2)] <- "Overlap"

d <- data.frame(seqnames = seqnames(gr2_out), start = start(gr2_out)-1, end = end(gr2_out), ID = gr2_out$name,
                Terekhanova_Peaks = gr2_out$Terekhanova, Regner_Peaks = gr2_out$Regner)
write.csv(d, "output/Tables/ST13_cCRE_hg38.csv")

d <- data.frame(seqnames = seqnames(gr_DA2), start = start(gr_DA2)-1, end = end(gr_DA2), cluster = gr_DA2$cluster)
write.csv(d, "output/Tables/ST14_DAcCRE_hg38.csv")

#----- Peak overlaps, DA regions -----#
#Terekhanova
ClusterPeaksList_T <- readRDS("../A2501_EpiEncode/Terekhanova_v1/rds/02_ClusterPeaks.rds")
ClusterPeaksList_T <- getMarkers(ClusterPeaksList_T, cutOff = "FDR < 0.25 & Log2FC > 0", returnGR = T)
names(ClusterPeaksList_T)

PeakOvlpJI_T <- mclapply(names(gr_DA2_ls), function(x){
  out <- lapply(names(ClusterPeaksList_T), function(y){
    intersect_num <- findOverlaps(gr_DA2_ls[[x]], ClusterPeaksList_T[[y]]) %>% length()
    union_num <- length(gr_DA2_ls[[x]]) + length(ClusterPeaksList_T[[y]]) - intersect_num
    return(intersect_num/union_num)
  })
  out <- unlist(out)
  return(out)
}, mc.cores = 24)
PeakOvlpJI_T <- do.call(rbind, PeakOvlpJI_T)
colnames(PeakOvlpJI_T) <- names(ClusterPeaksList_T)
rownames(PeakOvlpJI_T) <- names(gr_DA2_ls)
PeakOvlpJI_T

fh = function(x) hclust(dist(x), method = "ward.D2")
col_fun1 <- colorRamp2(c(0,0.1,0.3), c("white", "orange", "red"))

row_ha1 = rowAnnotation(No = anno_barplot(lapply(gr_DA2_ls, length) %>% unlist))
column_ha1 = HeatmapAnnotation(No = anno_barplot(lapply(ClusterPeaksList_T, length) %>% unlist))

ht1 <- Heatmap(PeakOvlpJI_T, name = "Jaccard similarity", 
               row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6), 
               top_annotation = column_ha1, right_annotation = row_ha1,
               column_title = "Terekhanova et al",
               cluster_rows = F, cluster_columns = fh, col = col_fun1, border = "black")
p2 <- draw(ht1)

#Regner
ClusterPeaksList_R <- readRDS("../A2501_EpiEncode/regner_v1/rds/01_ClusterPeaks.rds")
ClusterPeaksList_R <- getMarkers(ClusterPeaksList_R, cutOff = "FDR < 0.25 & Log2FC > 0", returnGR = T)
names(ClusterPeaksList_R)

PeakOvlpJI_R <- mclapply(names(gr_DA2_ls), function(x){
  out <- lapply(names(ClusterPeaksList_R), function(y){
    intersect_num <- findOverlaps(gr_DA2_ls[[x]], ClusterPeaksList_R[[y]]) %>% length()
    union_num <- length(gr_DA2_ls[[x]]) + length(ClusterPeaksList_R[[y]]) - intersect_num
    return(intersect_num/union_num)
  })
  out <- unlist(out)
  return(out)
}, mc.cores = 24)
PeakOvlpJI_R <- do.call(rbind, PeakOvlpJI_R)
colnames(PeakOvlpJI_R) <- names(ClusterPeaksList_R)
rownames(PeakOvlpJI_R) <- names(gr_DA2_ls)
PeakOvlpJI_R
max(PeakOvlpJI_R)

column_ha2 = HeatmapAnnotation(No = anno_barplot(lapply(ClusterPeaksList_R, length) %>% unlist))

ht2 <- Heatmap(PeakOvlpJI_R, name = "Jaccard similarity", 
               row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6), 
               top_annotation = column_ha2, right_annotation = row_ha1,
               column_title = "Regner et al",
               cluster_rows = F, cluster_columns = fh, col = col_fun1, border = "black")
p3 <- draw(ht2)
pdf("output/Plots/09_Heatmap_JI_DAcCREs_acrossDataSet.pdf", width = 6, height = 6)
p2
p3
dev.off()

#----- Compare TCGA data -----#
se <- readRDS("/mnt/host_mnt/Volumes/NEXTSSD1/Analysis/A2_CA_BC/v3/rds/TCGA_BRCA_ATAC_UQ_se.rds")
TCGA_ATACpeak <- rowRanges(se)

#ours vs TCGA
gr2_TCGA <- subsetByOverlaps(gr2, TCGA_ATACpeak)
gr2_noTCGA <- subsetByOverlaps(gr2, TCGA_ATACpeak, invert = T)

tmp1 <- subsetByOverlaps(gr_DA2,gr2_TCGA)
tmp2 <- subsetByOverlaps(gr_DA2,gr2_noTCGA)
df2.1 <- data.frame(cluster=names(table(tmp1$cluster)),number=as.numeric(table(tmp1$cluster)),overlap="Overlap w/ TCGA peaks")
df2.2 <- data.frame(cluster=names(table(tmp2$cluster)),number=as.numeric(table(tmp2$cluster)),overlap="NonOverlap")
df2.3 <- data.frame(cluster=names(table(gr_DA2$cluster)),number=as.numeric(table(gr_DA2$cluster)),overlap="Total DA-cCREs")
df2 <- rbind(df2.1,df2.2)
df2 <- rbind(df2, df2.3)
df2$cluster <- factor(df2$cluster, levels = c(paste0("C",c(5:13,15:16)), paste0("Ep",c(1:35))))

p3.2 <- ggplot(df2, aes(x=cluster,y=number,fill=overlap)) + geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(expand = c(0,0)) + theme_cb() + labs(x="Cluster",y="# cCREs") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("Overlap w/ TCGA peaks"="orange","NonOverlap"="gray", "Total DA-cCREs"="black"))

pdf("output/Plots/09_Barplot_TCGAoverlaps_cCRE.pdf", width = 8, height = 3)
p3.2
dev.off()

#----- Compare pEPN -----#
pEPN_final <- readRDS("rds/06_pEPN_final.rds")
gr <- arc@peakSet

gr_pEPN <- lapply(names(pEPN_final), function(x){
  idx <- pEPN_final[[x]]
  out <- gr[idx]
  mcols(out) <- NULL
  out$gene <- gsub("Promoter_","",x)
  return(out)
})
gr_pEPN <- unlist(GRangesList(gr_pEPN))
d <- data.frame(seqnames = seqnames(gr_pEPN), start = start(gr_pEPN)-1, end = end(gr_pEPN), gene = gr_pEPN$gene)
write.table(d, "output/output_bed/09_pEPN_hg19.bed", row.names = F, col.names = F, quote = F, sep = "\t")

gr_pEPN <- import.bed("output/output_bed/09_pEPN_hg38.bed")
length(unique(gr_pEPN$name))
# 89633 peaks, 9793
d <- data.frame(seqnames = seqnames(gr_pEPN), start = start(gr_pEPN)-1, end = end(gr_pEPN), gene = gr_pEPN$name)
write.csv(d, "output/Tables/ST15_pEnhancer_hg38.csv")

#Terekhanova
arc_T
p2g_500k <- getPeak2GeneLinks(arc_T, corCutOff = 0.4, FDRCutOff = 0.01, resolution = 500, returnLoops = F) #28423 peak-gene associations

p2g_500k_gr <- mclapply(unique(p2g_500k$idxRNA), function(i){
  idx1 <- p2g_500k$idxATAC[which(p2g_500k$idxRNA==i)]
  out <- metadata(p2g_500k)$peakSet[idx1]
  out$gene <- metadata(p2g_500k)$geneSet[i]$name
  return(out)
}, mc.cores = 16)
p2g_500k_gr <- unlist(GRangesList(p2g_500k_gr))
length(unique(p2g_500k_gr$gene))
# 276400 peaks, 12201 genes

genes_overlap1 <- intersect(unique(gr_pEPN$name), unique(p2g_500k_gr$gene)) #5916
JI_pEPN_p2g <- mclapply(genes_overlap1, function(x){
  gr_tmp1 <- gr_pEPN[which(gr_pEPN$name==x)]
  gr_tmp2 <- p2g_500k_gr[which(p2g_500k_gr$gene==x)]
  
  intersect_num <- findOverlaps(gr_tmp1, gr_tmp2) %>% length()
  union_num <- length(gr_tmp1) + length(gr_tmp2) - intersect_num
  return(intersect_num/union_num)
}, mc.cores = 24)
JI_pEPN_p2g <- unlist(JI_pEPN_p2g)
names(JI_pEPN_p2g) <- genes_overlap1
JI_pEPN_p2g <- sort(JI_pEPN_p2g, decreasing = T)

Count_pEPN_p2g <- mclapply(genes_overlap1, function(x){
  gr_tmp1 <- gr_pEPN[which(gr_pEPN$name==x)]
  gr_tmp2 <- p2g_500k_gr[which(p2g_500k_gr$gene==x)]
  intersect_num <- findOverlaps(gr_tmp1, gr_tmp2) %>% length()
  return(intersect_num)
}, mc.cores = 24)
Count_pEPN_p2g <- unlist(Count_pEPN_p2g)
names(Count_pEPN_p2g) <- genes_overlap1
Count_pEPN_p2g <- sort(Count_pEPN_p2g, decreasing = T)

saveRDS(Count_pEPN_p2g, "rds/09_Count_pEPN_p2g.rds")

### Regnar Basal
basal_p2g <- read.csv("ref/mmc5.csv")
head(basal_p2g)

str_split(str_split(basal_p2g$peakName, ":", simplify = T)[,2], "-", simplify = T)[,1]
str_split(str_split(basal_p2g$peakName, ":", simplify = T)[,2], "-", simplify = T)[,2]

gr_basal_p2g <- GRanges(seqnames = str_split(basal_p2g$peakName, ":", simplify = T)[,1],
                        IRanges(start=as.numeric(str_split(str_split(basal_p2g$peakName, ":", simplify = T)[,2], "-", simplify = T)[,1]),
                                end=as.numeric(str_split(str_split(basal_p2g$peakName, ":", simplify = T)[,2], "-", simplify = T)[,2])),
                        cancer = basal_p2g$significant_effect_size_in_cancer_condition,
                        normal = basal_p2g$significant_effect_size_in_normal_condition,
                        gene = basal_p2g$geneName)
length(which(gr_basal_p2g$cancer==T & gr_basal_p2g$normal == T))
length(which(gr_basal_p2g$cancer==T & gr_basal_p2g$normal == F))
length(which(gr_basal_p2g$cancer==F & gr_basal_p2g$normal == T))

genes_overlap2 <- intersect(unique(gr_pEPN$name), unique(gr_basal_p2g$gene)) #7659
Count_pEPN_RegnerBasal <- mclapply(genes_overlap2, function(x){
  gr_tmp1 <- gr_pEPN[which(gr_pEPN$name==x)]
  gr_tmp2 <- gr_basal_p2g[which(gr_basal_p2g$gene==x)]
  intersect_num <- findOverlaps(gr_tmp1, gr_tmp2) %>% length()
  return(intersect_num)
}, mc.cores = 24)
Count_pEPN_RegnerBasal <- unlist(Count_pEPN_RegnerBasal)
names(Count_pEPN_RegnerBasal) <- genes_overlap2
Count_pEPN_RegnerBasal <- sort(Count_pEPN_RegnerBasal, decreasing = T)
plot(Count_pEPN_RegnerBasal)
saveRDS(Count_pEPN_RegnerBasal, "rds/09_Count_pEPN_RegnerBasal.rds")

### Regnar Luminal
luminal_p2g <- read.csv("ref/mmc8.csv")
head(luminal_p2g)

gr_luminal_p2g <- GRanges(seqnames = str_split(luminal_p2g$peakName, ":", simplify = T)[,1],
                        IRanges(start=as.numeric(str_split(str_split(luminal_p2g$peakName, ":", simplify = T)[,2], "-", simplify = T)[,1]),
                                end=as.numeric(str_split(str_split(luminal_p2g$peakName, ":", simplify = T)[,2], "-", simplify = T)[,2])),
                        cancer = luminal_p2g$significant_effect_size_in_cancer_condition,
                        normal = luminal_p2g$significant_effect_size_in_normal_condition,
                        gene = luminal_p2g$geneName)
length(which(gr_luminal_p2g$cancer==T & gr_luminal_p2g$normal == T))
length(which(gr_luminal_p2g$cancer==T & gr_luminal_p2g$normal == F))
length(which(gr_luminal_p2g$cancer==F & gr_luminal_p2g$normal == T))

genes_overlap3 <- intersect(unique(gr_pEPN$name), unique(gr_luminal_p2g$gene)) #7659
Count_pEPN_Regnerluminal <- mclapply(genes_overlap3, function(x){
  gr_tmp1 <- gr_pEPN[which(gr_pEPN$name==x)]
  gr_tmp2 <- gr_luminal_p2g[which(gr_luminal_p2g$gene==x)]
  intersect_num <- findOverlaps(gr_tmp1, gr_tmp2) %>% length()
  return(intersect_num)
}, mc.cores = 8)
Count_pEPN_Regnerluminal <- unlist(Count_pEPN_Regnerluminal)
names(Count_pEPN_Regnerluminal) <- genes_overlap3
Count_pEPN_Regnerluminal <- sort(Count_pEPN_Regnerluminal, decreasing = T)
plot(Count_pEPN_Regnerluminal)

saveRDS(Count_pEPN_Regnerluminal, "rds/09_Count_pEPN_Regnerluminal.rds")

Count_pEPN_p2g <- readRDS("rds/09_Count_pEPN_p2g.rds")
Count_pEPN_RegnerBasal <- readRDS("rds/09_Count_pEPN_RegnerBasal.rds")
Count_pEPN_Regnerluminal <- readRDS("rds/09_Count_pEPN_Regnerluminal.rds")

df4 <- c("Terekhanova" = length(which(Count_pEPN_p2g > 0)),
  "Regner_Basal"  = length(which(Count_pEPN_RegnerBasal > 0)),
  "Regner_Luminal" = length(which(Count_pEPN_Regnerluminal > 0)))
p4 <- lapply(names(df4), function(x){
  df_tmp <- data.frame(group=c("overlap", "nonoverlap"), value=c(df4[x], 9793-df4[x]))
  out <-ggplot(df_tmp, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1, color="white") + ggtitle(x) + 
    coord_polar("y", start=0) + theme_cb() +scale_fill_manual(values = c("overlap"="orange","nonoverlap"="gray"))
  return(out)
  })
pdf("output/Plots/09_Pie_NumOvlpGenes_P2G.pdf", width = 6, height = 6)
p4
dev.off()

df5.1 <- data.frame(gene = names(Count_pEPN_p2g), overlap_count = as.numeric(Count_pEPN_p2g))
df5.2 <- data.frame(gene = names(Count_pEPN_RegnerBasal), overlap_count = as.numeric(Count_pEPN_RegnerBasal))
df5.3 <- data.frame(gene = names(Count_pEPN_Regnerluminal), overlap_count = as.numeric(Count_pEPN_Regnerluminal))
write.csv(df5.1, "output/Tables/09_OverlapEnhancerCounts_Terekhanova.csv")
write.csv(df5.2, "output/Tables/09_OverlapEnhancerCounts_RegnerBasal.csv")
write.csv(df5.3, "output/Tables/09_OverlapEnhancerCounts_RegnerLuminal.csv")
