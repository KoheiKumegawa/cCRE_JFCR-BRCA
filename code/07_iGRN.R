#----------------------------------------------------------------------------
# 07_iGRN.R
#----------------------------------------------------------------------------
library(ArchR)
library(rtracklayer)
library(ggsignif)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(stringr)
library(ggrastr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
addArchRGenome("hg19")

arc <- readRDS("rds/06r_arc.rds")

#----- co-accessibility +-500kb -----#
arc <- addCoAccessibility(arc, reducedDims = "IterativeLSI", maxDist = 500000)
saveRDS(arc, "rds/07_arc_coA.rds")

coA_all <- getCoAccessibility(arc, corCutOff = 0, resolution = 500, returnLoops = FALSE)
coA_all_corr <- data.frame(corr = coA_all$correlation, FDR = coA_all$FDR) #20828478 -> 10414239 co-accessibility

# coA_all_corr <- coA_all_corr[order(coA_all_corr$FDR),]
# coA_all_corr$corr[max(which(coA_all_corr$FDR < 0.01))]

p1 <- ggplot(coA_all_corr, aes(x=corr)) + 
  geom_density(color = "lightgray", fill="red",alpha=.2) + theme_ArchR() +
  geom_vline(xintercept = 0.4, lty = "dotted")
pdf("output/Plots/07_density_coA_all.pdf", width = 5, height = 5)
p1
dev.off()

coA_signif <- getCoAccessibility(arc, corCutOff = 0.4, resolution = 500, returnLoops = FALSE) #1396326 -> 698163 co-accessibility
coA_signif_gr <- getCoAccessibility(arc, corCutOff = 0.4, resolution = 500, returnLoops = T)
coA_signif_gr <- coA_signif_gr$CoAccessibility

#profile coA
gr <- arc@peakSet
names(gr) <- c(paste0("cCRE_", c(1:229238)))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak(gr, tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
mcols(gr) <- as.GRanges(peakAnno)

km2 <- readRDS("rds/06-R2_km2.rds")
gr$module <- NA
gr[names(km2$cluster)]$module <- paste0("rm", km2$cluster)

start_gr <- resize(coA_signif_gr, width = 1, fix = "start")
end_gr <- resize(coA_signif_gr, width = 1, fix = "end")

genomic_anno <- gr$annotation
genomic_anno <- str_split(genomic_anno, pattern = "\\(", simplify = T)[,1] %>% gsub(" ", "", .)

df <- data.frame(peak1 = genomic_anno[findOverlaps(gr, start_gr) %>% queryHits()],
                 peak2 = genomic_anno[findOverlaps(gr, end_gr) %>% queryHits()])
df <- table(df) %>% melt()
df$value %>% sum
write.csv(df, "output/Tables/07_pEG_regionProlile.csv")
df2 <- read.csv("output/Tables/07_pEG_regionProlile3.csv", header = F)

df2 <- df2[order(df2$V1, decreasing = T),] %>% as.data.frame %>% `row.names<-`(.,NULL)
pie_colors <- c(ArchRPalettes$stallion, ArchRPalettes$calm)[c(1:28)] %>% `names<-`(., df2$V2)
p2 <- ggplot(df2, aes(x = "x", y = V1, fill = reorder(V2, V1))) +
  geom_bar(stat = "identity", position = "stack", color = "white") + 
  coord_polar(theta = "y") + theme_ArchR() + scale_fill_manual(values = pie_colors)
pdf("output/Plots/07_class_coA_all.pdf", width = 8, height = 5)
p2
dev.off()

df3 <- df2
df3$V1 <- (df3$V1 / sum(df3$V1))*100
df3$V1[df3$V1 < 1]

df4 <- df2[c(1:13),]
df4$V1[13] <- sum(df2$V1[df3$V1 < 1])
df4$V2[13] <- "Others"
pie_colors2 <- c(ArchRPalettes$stallion, ArchRPalettes$calm)[c(1:13)] %>% `names<-`(., df4$V2)
p3 <- ggplot(df4, aes(x = "x", y = V1, fill = factor(V2, levels = df4$V2))) +
  geom_bar(stat = "identity", position = "stack", color = "white") + 
  coord_polar(theta = "y") + theme_ArchR() + scale_fill_manual(values = pie_colors2)
pdf("output/Plots/07_class_coA_all_v2.pdf", width = 8, height = 5)
p3
dev.off()

#----- only promoter -----#
coA_signif_promoter <- coA_signif[which(genomic_anno[coA_signif$queryHits] == "Promoter"),]
coA_signif_promoter2peak <- mclapply(unique(coA_signif_promoter$queryHits), function(i) coA_signif_promoter$subjectHits[coA_signif_promoter$queryHits == i], mc.cores = 12)
names(coA_signif_promoter2peak) <- paste0("cCRE_", unique(coA_signif_promoter$queryHits))

saveRDS(coA_signif_promoter2peak, "rds/07_coA_signif_promoter2peak.rds")
# coA_signif_promoter2peak <- readRDS("rds/07_coA_signif_promoter2peak.rds")

coA_signif_gene2peak <- mclapply(unique(gr$SYMBOL), function(x){
  idx <- names(gr)[which(gr$SYMBOL == x)]
  tmp <- which(names(coA_signif_promoter2peak) %in% idx)
  out <- unlist(coA_signif_promoter2peak[tmp]) %>% unique %>% sort
  return(out)
}, mc.cores = 10)
names(coA_signif_gene2peak) <- unique(gr$SYMBOL)
coA_signif_gene2peak <- coA_signif_gene2peak[!is.na(names(coA_signif_gene2peak))]

# coA_signif_gene2peak <- readRDS("rds/07_coA_signif_gene2peak.rds")

Module_order <- c("T","T-C6","B","Plasma","Fibroblast","Endothelial","Myeloid","ER",
                  "EP1","EP3","EP4","EP5","EP6","EP8","EP9","EP11","EP16","EP18","EP19","EP21","EP22","EP23","EP24-EP25","EP26","EP29","EP30")  
gr_modules <- lapply(Module_order, function(i) import.bed(paste0("output/output_bed/06-R2_modules/", i, "_module.bed")))
names(gr_modules) <- Module_order

module_tgt_genes <- mclapply(Module_order, function(i){
  mod_gr <- import.bed(paste0("output/output_bed/06-R2_modules/", i, "_module.bed"))
  id <- names(gr)[findOverlaps(gr, mod_gr) %>% queryHits(.)]
  res <- lapply(coA_signif_gene2peak, function(n) intersect(paste0("cCRE_", n), id))
  res <- lapply(res, length) %>% unlist() %>% sort(., decreasing = T)
  out <- data.frame(gene = names(res), peakno = as.numeric(res))
  return(out)
}, mc.cores = 10)
names(module_tgt_genes) <- Module_order

saveRDS(coA_signif_promoter, "rds/07_coA_signif_promoter.rds")
saveRDS(coA_signif_gene2peak, "rds/07_coA_signif_gene2peak.rds")
saveRDS(module_tgt_genes, "rds/07_module_tgt_genes.rds")
lapply(names(module_tgt_genes), function(i) write.csv(module_tgt_genes[[i]], paste0("output/Tables/module_tgt_genes/", i, "_module.csv")))

#----- iGRN and target genes of each module -----#
coA_gene_df <- lapply(coA_signif_gene2peak, length) %>% unlist %>% sort(.,decreasing = T)
coA_gene_df <- data.frame(gene = names(coA_gene_df), peak = as.numeric(coA_gene_df))
which(coA_gene_df$peak >= 2) %>% length # 14361 inferred gene regulatory network, iGRN

iGRN_gene_module_df <- lapply(module_tgt_genes, function(x) which(x$peakno >= 5) %>% length) %>% unlist
iGRN_gene_module_df <- data.frame(module = names(iGRN_gene_module_df), tgtGeneNo = as.numeric(iGRN_gene_module_df))

p4 <- ggplot(iGRN_gene_module_df, aes(x=factor(module,levels = Module_order), y=tgtGeneNo)) + geom_bar(stat = "identity") + 
      theme_ArchR() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("output/Plots/07_Barplot_moduleTgtGenesNo.pdf", width = 5, height = 2.5)
p4
dev.off()

for(i in seq_along(module_tgt_genes)){
  module_tgt_genes[[i]]$label <- ""
  module_tgt_genes[[i]]$label[c(1:10)] <- module_tgt_genes[[i]]$gene[c(1:10)]
}
p5 <- lapply(names(module_tgt_genes), function(i){
  df <- module_tgt_genes[[i]]
  df$rank <- c(1:nrow(df))
  out <- ggplot(df, aes(x = rank, y = peakno, label = label)) + 
    geom_point_rast() + theme_ArchR() + geom_text_repel(max.overlaps = 1000) +
    labs(x = "Ranked genes", y = "Number of correlated cCREs") + ggtitle(i)
  return(out)
})
pdf("output/Plots/07_RankPlot_Module_iGRN_overlaps.pdf", width = 6, height = 6)
p5
dev.off()

#----- Example, FOXA1 -----#
gr_interest <- geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == "FOXA1")]
gr_interest <- resize(gr_interest, fix = "start", width = 1)
gr_interest <- extendGR(gr_interest, upstream = 500000, downstream = 500000)

idx <- names(gr)[which(gr$SYMBOL == "FOXA1" & genomic_anno == "Promoter")]
loop_gr <- coA_signif_gr[c(findOverlaps(start_gr, gr[idx]) %>% queryHits(),
                           findOverlaps(end_gr, gr[idx]) %>% queryHits()) %>% unique()]

FOXA1_iGRN <- gr[paste0("cCRE_", coA_signif_gene2peak$FOXA1)]

epicluster_colors <- readRDS("rds/03_epicluster_colors.rds")
cluster_colors <- readRDS("rds/02_cluster_colors.rds")
names(epicluster_colors) <- gsub("C", "EP", names(epicluster_colors))
colors_j <- c(cluster_colors[paste0("C", c(5,6,7,8,9,10,11,12,13,14,15,18,19))], epicluster_colors)
p6 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "Clusters2",
  pal = colors_j,
  region = gr_interest,
  loops = loop_gr,
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"))
p7 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "Clusters2",
  pal = colors_j,
  region = gr_interest,
  loops = loop_gr,
  features = FOXA1_iGRN,
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"))
p7.2 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "Clusters2",
  pal = colors_j,
  region = gr_interest,
  loops = loop_gr,
features = FOXA1_iGRN[FOXA1_iGRN$module == "rm10"],
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"))
plotPDF(p6,p7,p7.2, name = "07_GenomeTrack_FOXA1_allClusters.pdf", ArchRProj = arc, addDOC = FALSE, width = 8, height = 10)

lapply(c("rm10", "rm12", "rm14", "rm18","rm2","rm23","rm6"), function(i) export.bed(FOXA1_iGRN[FOXA1_iGRN$module == i], paste0("output/output_bed/07_FOXA1_iGRN/07_FOXA1_iGRN_", i, ".bed")))
export.bed(FOXA1_iGRN, "output/output_bed/07_FOXA1_iGRN/07_FOXA1_iGRN_all.bed")

#----- ER+ tumors heterogeneity -----#
er_response_genes <- read.csv("ref/MSigDB_ER_response_genes_v2.csv", header = F)[1,] %>% as.character() %>% unique %>% sort

er_gep <- lapply(c("ER", "EP11", "EP16", "EP22", "EP29", "EP30"), function(i){
  tmp <- module_tgt_genes[[i]]
  tmp <- na.omit(tmp)
  rownames(tmp) <- tmp$gene
  out <- tmp[er_response_genes,"peakno"]
  return(out)
})
er_gep <- do.call(rbind, er_gep)
rownames(er_gep) <- c("ER", "EP11", "EP16", "EP22", "EP29", "EP30")
colnames(er_gep) <- er_response_genes

er_gep <- er_gep[, order(er_gep[1,], decreasing = T)]
write.csv(er_gep, "output/Tables/07_Heatmap_ERResponseGenes_NocCRE.csv")
er_gep2 <- er_gep[c(2:6),]

col_fun1 <- colorRamp2(c(0,5,10,15,20), colors = c("white", rev(viridis(n = 4, option = "G"))))
ha1 <- HeatmapAnnotation(ER = anno_points(er_gep[1,]))

ht1 <- Heatmap(er_gep2, cluster_columns = F, cluster_rows = F, col = col_fun1, top_annotation = ha1,
               name = "No.cCRE", column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 6), border = "black")
p8 <- draw(ht1)
pdf("output/Plots/07_Heatmap_ERResponseGenes_NocCRE.pdf", width = 10, height = 2)
p8
dev.off()

#all gene tracks
lapply(er_response_genes, function(i){
  gr_interest <- geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == i)]
  gr_interest <- resize(gr_interest, fix = "start", width = 1)
  gr_interest <- extendGR(gr_interest, upstream = 500000, downstream = 500000)
  
  idx <- names(gr)[which(gr$SYMBOL == i & genomic_anno == "Promoter")]
  loop_gr <- coA_signif_gr[c(findOverlaps(start_gr, gr[idx]) %>% queryHits(), findOverlaps(end_gr, gr[idx]) %>% queryHits()) %>% unique()]
  out <- plotBrowserTrack(
    ArchRProj = arc, 
    groupBy = "Clusters2", useGroups = c("EP11", "EP16", "EP22", "EP29", "EP30"),
    pal = epicluster_colors,
    region = gr_interest,
    loops = loop_gr,
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"))
  plotPDF(out, name = paste0("07_GenomeTrack_ERResponseGenes_", i, ".pdf"), ArchRProj = arc, addDOC = FALSE, width = 8, height = 10)
  return(NULL)
})

#TFF1
gr_interest <- geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == "TFF1")]
gr_interest <- resize(gr_interest, fix = "start", width = 1)
gr_interest <- extendGR(gr_interest, upstream = 500000, downstream = 500000)

idx <- names(gr)[which(gr$SYMBOL == "TFF1" & genomic_anno == "Promoter")]
loop_gr <- coA_signif_gr[c(findOverlaps(start_gr, gr[idx]) %>% queryHits(),
                           findOverlaps(end_gr, gr[idx]) %>% queryHits()) %>% unique()]

TFF1_iGRN <- gr[paste0("cCRE_", coA_signif_gene2peak$TFF1)]
p9 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "Clusters2", useGroups = c("EP11", "EP16", "EP22", "EP29", "EP30"),
  pal = epicluster_colors,
  region = gr_interest,
  loops = loop_gr,
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"))
p10 <- lapply(c("ER","EP11", "EP16", "EP22", "EP29", "EP30"), function(i)
  plotBrowserTrack(
    ArchRProj = arc, 
    groupBy = "Clusters2", useGroups = c("EP11", "EP16", "EP22", "EP29", "EP30"),
    pal = epicluster_colors,
    region = gr_interest,
    features = gr_modules[[i]],
    loops = loop_gr, title = paste0(i, " cCRE module"),
    plotSummary = c("bulkTrack", "featureTrack")))
p10.2 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "Clusters2", useGroups = c("EP11", "EP16", "EP22", "EP29", "EP30"),
  pal = epicluster_colors,
  region = gr_interest,
  loops = loop_gr,
  features = TFF1_iGRN, title = "TFF1 iGRN (47 cCREs)",
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"))
plotPDF(p9,p10[[1]],p10[[2]],p10[[3]],p10[[4]],p10[[5]],p10[[6]],p10.2, name = "07_GenomeTrack_TFF1_v2.pdf", ArchRProj = arc, addDOC = FALSE, width = 8, height = 8)

#----- CNA interaction -----#
iGRN <- coA_gene_df[coA_gene_df$peak >=2,]

# TCGA Nature 2012
amplified_genes <- read.csv("ref/TCGA_nature2012.csv", header = F)[1,] %>% as.character()

#amplified regions
pAR_common_merge <- import.bed("output/output_bed/04_pAR_ccommon_merge.bed")
pAR_specific_merge <- import.bed("output/output_bed/04_pAR_specific_merge.bed")

gene_gr <- geneAnnoHg19$genes
names(gene_gr) <- gene_gr$symbol
gene_gr <- gene_gr[intersect(iGRN$gene, names(gene_gr))]

idx <- findOverlaps(gene_gr, c(pAR_common_merge, pAR_specific_merge)) %>% queryHits()
idx <- as.character(gene_gr$symbol[idx])

iGRN$TCGAamp <- "No"
iGRN$TCGAamp[iGRN$gene %in% amplified_genes] <- "Yes"
iGRN$scATACamp <- "No"
iGRN$scATACamp[iGRN$gene %in% idx] <- "Yes"
iGRN$amp <- "No"
iGRN$amp[iGRN$TCGAamp == "Yes" | iGRN$scATACamp == "Yes"] <- "Yes"

iGRN$rank <- c(1:nrow(iGRN))
iGRN$label <- ""
iGRN$label[c(1:10)] <- iGRN$gene[c(1:10)]

iGRN$label[which(iGRN$gene %in% amplified_genes)] <- iGRN$gene[which(iGRN$gene %in% amplified_genes)]

p13 <- ggplot(iGRN, aes(x = rank, y = peak, label = label)) + 
  geom_point_rast() + theme_ArchR() + geom_text_repel(max.overlaps = 1000) +
  labs(x = "Ranked iGRNs", y = "Number of cCREs in each iGRN")
p14 <- ggplot(iGRN, aes(x = amp, y = peak, fill = amp)) + 
  geom_violin(alpha = 0.4) + geom_boxplot(outlier.size = 0, alpha = 0.4) + theme_ArchR() + 
  geom_signif(comparisons = list(c("No", "Yes")), test = "wilcox.test") + 
  scale_fill_manual(values = c("No" = "#E6E7E8", "Yes" = "#8816A7")) +
  labs(x = "Amplification", y = "Number of cCREs in each iGRN")
write.csv(coA_gene_df, "output/Tables/07_coA_gene_peakNo_ampOverlaps.csv")
pdf("output/Plots/07_RankPlot_correlated_cCRE_perGenes.pdf", width = 6, height = 6)
p13
dev.off()
pdf("output/Plots/07_BoxPlot_correlated_cCRE_Amp_v2.pdf", width = 3, height = 5)
p14
dev.off()

#----- ERBB2 -----#
gr_interest <- geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == "ERBB2")]
gr_interest <- resize(gr_interest, fix = "start", width = 1)
gr_interest <- extendGR(gr_interest, upstream = 500000, downstream = 500000)

idx <- names(gr)[which(gr$SYMBOL == "ERBB2" & genomic_anno == "Promoter")]
loop_gr <- coA_signif_gr[c(findOverlaps(start_gr, gr[idx]) %>% queryHits(),
                           findOverlaps(end_gr, gr[idx]) %>% queryHits()) %>% unique()]

names(epicluster_colors) <- gsub("C", "EP", names(epicluster_colors))
p15 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "Clusters2", useGroups = c("EP5", "EP6"),
  pal = epicluster_colors,
  region = gr_interest,
  loops = loop_gr,
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"))
p16 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "Clusters2", useGroups = c("EP5", "EP6"),
  pal = epicluster_colors,
  region = gr_interest,
  features = gr_modules$EP5,
  loops = loop_gr, title = "EP5 cCRE module",
  plotSummary = c("bulkTrack", "featureTrack"))
p17 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "Clusters2", useGroups = c("EP5", "EP6"),
  pal = epicluster_colors,
  region = gr_interest,
  features = gr_modules$EP6,
  loops = loop_gr, title = "EP6 cCRE module",
  plotSummary = c("bulkTrack", "featureTrack"))
plotPDF(p15,p16,p17, name = "07_GenomeTrack_ERBB2.pdf", ArchRProj = arc, addDOC = FALSE, width = 8, height = 8)

#----- CCND1 -----#
gr_interest <- geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == "CCND1")]
gr_interest <- resize(gr_interest, fix = "start", width = 1)
gr_interest <- extendGR(gr_interest, upstream = 500000, downstream = 500000)

idx <- names(gr)[which(gr$SYMBOL == "CCND1" & genomic_anno == "Promoter")]
loop_gr <- coA_signif_gr[c(findOverlaps(start_gr, gr[idx]) %>% queryHits(),
                           findOverlaps(end_gr, gr[idx]) %>% queryHits()) %>% unique()]

CCND1_iGRN <- gr[paste0("cCRE_", coA_signif_gene2peak$CCND1)]
p18 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "Clusters2",
  region = gr_interest,
  loops = loop_gr, pal = colors_j,
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"))
plotPDF(p18, name = "07_GenomeTrack_CCND1.pdf", ArchRProj = arc, addDOC = FALSE, width = 8, height = 8)
p19 <- plotBrowserTrack(
    ArchRProj = arc, 
    groupBy = "Clusters2",
    region = gr_interest,
    features = CCND1_iGRN,
    loops = loop_gr,
    plotSummary = c("bulkTrack", "featureTrack"))
plotPDF(p19, name = "07_GenomeTrack_CCND1_iGRN.pdf", ArchRProj = arc, addDOC = FALSE, width = 8, height = 8)

#----- TFF1 details -----#
MarkerPeakClusters2 <- getMarkerFeatures(
  ArchRProj = arc, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
## Caution! Line up chromosomal order!! (chr1 -> chr10 -> chr11... to chr1 -> chr2 -> chr3)
rowranges1 <- GRanges(seqnames = rowData(MarkerPeakClusters2)$seqnames, 
                      IRanges(start = rowData(MarkerPeakClusters2)$start, end = rowData(MarkerPeakClusters2)$end))
idx <- findOverlaps(gr, rowranges1) %>% subjectHits()
MarkerPeakClusters2_ordered <- MarkerPeakClusters2[idx, ]
rownames(MarkerPeakClusters2_ordered) <- names(gr)
rowData(MarkerPeakClusters2_ordered)$seqnames


mtx <- assays(MarkerPeakClusters2_ordered)$Mean
mtx <- mtx[paste0("cCRE_", coA_signif_gene2peak$TFF1), paste0("EP",c(11,16,22,29,30))]
z_mtx <- t(scale(t(mtx)))
z_mtx <- na.omit(z_mtx)

rownames(z_mtx)
module_info <- gr[rownames(z_mtx)]$module
module_info[module_info %ni% paste0("rm", c(10,25,6,26,23,12))] <- "Others"
module_colors <- c(c(ArchRPalettes$bear)[c(1:6)], "darkgray") %>% `names<-`(., c(paste0("rm", c(10,25,6,26,23,12)), "Others"))

ra1 <- rowAnnotation(Module = module_info, col = list(Module = module_colors))

col_fun1 <- colorRamp2(c(-2,-1,0,1,2), colors = viridis(5, option = "D"))
ht1 <- Heatmap(z_mtx, cluster_rows = F, cluster_columns = F, col = col_fun1, name = "Accessibility z-score",
               show_row_dend = F, show_column_dend = F, show_row_names = F, 
               column_names_gp = gpar(fontsize = 8),
               left_annotation = ra1,
               use_raster = T, raster_by_magick = F)
p20 <- draw(ht1)
plotPDF(p20, name = "07_Heatmap_TFF1_iGRNs.pdf", ArchRProj = arc, addDOC = FALSE, width = 8, height = 8)

#----- iGRN info -----#
iGRN_cCRE <- coA_signif_gene2peak[iGRN$gene]
iGRN_cCRE <- lapply(iGRN_cCRE, function(x){
  out <- paste0("cCRE_", x)
  out <- paste(out, collapse = ";")
  return(out)
})
iGRN_cCRE <- do.call(rbind, iGRN_cCRE)
write.csv(iGRN_cCRE, "output/Tables/07_iGRN_cCRE_names.csv")
