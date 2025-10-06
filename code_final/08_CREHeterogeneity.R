#------------------------------------------------------------------------------
# 08_CREHeterogeneity.R
#------------------------------------------------------------------------------
library(ArchR)
library(ggplot2)
library(ggrastr)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(viridis)

arc <- readRDS("rds/06_arc.rds")
gr <- arc@peakSet
names(gr) <- c(paste0("cCRE_", c(1:length(gr))))

pEPN_final <- readRDS("rds/06_pEPN_final.rds")
IntegDACRE_pEPN <- read.csv("output/Tables/06_IntegDACRE_pEPN_nameSortedDF.csv", row.names = 1)
IntegDACRE_pEPN_gr <- readRDS("rds/06_IntegDACRE_pEPN_gr.rds")

cluster_info <- data.frame(
  Cluster = c("Ep9", "Ep8", "Ep26", "Ep2", "Ep1", "Ep6", "Ep5", "Ep11", "Ep10", "Ep7", "Ep19", "Ep25", "Ep24", "Ep13", "Ep14", "Ep20", "Ep17", "Ep15", "Ep33", "Ep21", "Ep28", "Ep18", "Ep16", "Ep12", "Ep31", "Ep30", "Ep32", "Ep35", "Ep34", "Ep29", "Ep22", "Ep23", "Ep27", "Ep4", "Ep3"),
  Sample = c("P200A", "P209PE", "P51", "P44", "P41", "P190", "P190", "Mixed", "Mixed", "P20LN", "Mixed", "P181", "P93", "P175PE", "P189", "P185", "Mixed", "P191", "P49", "P179", "P33", "P211LN", "P178", "P180", "Mixed", "P39", "P40", "P50", "P50", "P35", "P165", "P210LR", "P64", "P166", "P122"),
  Subtype = c("Metastatic", "Metastatic", "ER–/PGR+/HER–", "TN", "TN", "TN", "TN", "N/A", "N/A", "Metastatic", "N/A", "ER+/HER2–", "ER+/HER2+", "Metastatic", "ER+/HER2–", "ER+/HER2–", "N/A", "ER+/HER2–", "ER+/HER2–", "ER+/HER2–", "ER+/HER2–", "Metastatic", "ER+/HER2–", "ER+/HER2–", "N/A", "ER+/HER2–", "ER+/HER2–", "ER+/HER2–", "ER+/HER2–", "ER+/HER2–", "ER+/HER2–", "ER+/HER2–", "ER+/HER2–", "HER2+", "HER2+"),
  Malignancy = c("Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Non-malignant", "Undetermined", "Malignant", "Undetermined", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Non-malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Non-malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant", "Malignant")
)
rownames(cluster_info) <- cluster_info$Cluster

#----- Intra-tumor heterogeneity Ep5 vs Ep6 -----#
# promoter ≥ 5 pEnhancer overlapped with DA-cCREs
gene1 <- rownames(IntegDACRE_pEPN)[IntegDACRE_pEPN$Ep5 >= 5]
gene2 <- rownames(IntegDACRE_pEPN)[IntegDACRE_pEPN$Ep6 >= 5]
intersect(gene1, gene2)
VennDiagram::venn.diagram(list(Ep5 = gene1, Ep6 = gene2), filename="output/Plots/08_Venn_Ep5vsEp6.svg", scaled=T, imagetype = "svg")

promoters1 <- paste0("Promoter_", intersect(gene1, gene2))
df2 <- lapply(promoters1, function(x){
  gr1 <- subsetByOverlaps(IntegDACRE_pEPN_gr$Ep5, gr[pEPN_final[[x]]])
  gr2 <- subsetByOverlaps(IntegDACRE_pEPN_gr$Ep6, gr[pEPN_final[[x]]])
  idx1 <- findOverlaps(gr1, gr2) %>% queryHits() %>% length
  idx2 <- length(gr1) - idx1
  idx3 <- length(gr2) - idx1
  out <- c("Overlap" = idx1, "Ep5_only" = idx2, "Ep6_only" = idx3)
})
df2 <- do.call(rbind, df2)
rownames(df2) <- intersect(gene1, gene2)
write.csv(df2, "output/Tables/08_pEnhancerNo_Ep5vsEp6_v2.csv")

df2 <- reshape2::melt(df2)
p2 <- ggplot(df2, aes(x = reorder(Var1, -value), y = value, fill = Var2)) + geom_bar(stat = "identity") + theme_ArchR() +
  scale_fill_manual(values = c("Overlap" = "darkgray", "Ep5_only" = "darkblue", "Ep6_only" = "darkred")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("output/Plots/08_Barplot_pEnhancer_Ep5vsEp6_v2.pdf", width = 5, height = 4)
p2
dev.off()

#example, CRELD2
dt1 <- data.frame(seqnames = "chr22", 
                  start = start(geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == "CRELD2")]),
                  end = resize(gr[pEPN_final$Promoter_CRELD2], width = 1, fix = "center") %>% start)
idx1 <- which(dt1$start > dt1$end)
start1 <- dt1[idx1,]$end
end1 <- dt1[idx1,]$start
dt1[idx1,]$start <- start1
dt1[idx1,]$end <- end1

gr1 <- subsetByOverlaps(IntegDACRE_pEPN_gr$Ep5, gr[pEPN_final$Promoter_CRELD2])
gr2 <- subsetByOverlaps(IntegDACRE_pEPN_gr$Ep6, gr[pEPN_final$Promoter_CRELD2])
gr3 <- subsetByOverlaps(gr1, gr2)
gr1 <- subsetByOverlaps(gr1, gr3, invert = T)
gr2 <- subsetByOverlaps(gr2, gr3, invert = T)

gr_loop1 <- GRanges(seqnames = dt1$seqnames, IRanges(start = dt1$start, end = dt1$end)) # 44 enhancers
p3.1 <- plotBrowserTrack(
  ArchRProj = arc,
  groupBy = "Clusters2",
  useGroups = c("Ep5", "Ep6"),
  region = GRanges(seqnames = "chr22", ranges = IRanges(start = 49820001, end = 50760000)),
  features = gr3,
  loops = gr_loop1,
  pal = c("darkblue", "darkred"))
grid::grid.newpage()
grid::grid.draw(p3.1)

p3.2 <- plotBrowserTrack(
  ArchRProj = arc,
  groupBy = "Clusters2",
  useGroups = c("Ep5", "Ep6"),
  region = GRanges(seqnames = "chr22", ranges = IRanges(start = 49820001, end = 50760000)),
  features = gr1,
  loops = gr_loop1,
  pal = c("darkblue", "darkred"))
grid::grid.newpage()
grid::grid.draw(p3.2)

p3.3 <- plotBrowserTrack(
  ArchRProj = arc,
  groupBy = "Clusters2",
  useGroups = c("Ep5", "Ep6"),
  region = GRanges(seqnames = "chr22", ranges = IRanges(start = 49820001, end = 50760000)),
  features = gr2,
  loops = gr_loop1,
  pal = c("darkblue", "darkred"))
grid::grid.newpage()
grid::grid.draw(p3.3)

pdf("output/Plots/08_GenomeTrack_CRELD2.pdf", width = 10, height = 4)
grid::grid.newpage()
grid::grid.draw(p3.1)
grid::grid.newpage()
grid::grid.draw(p3.2)
grid::grid.newpage()
grid::grid.draw(p3.3)
dev.off()

#----- Intra-tumor heterogeneity Ep34 vs Ep35 -----#
# promoter ≥ 5 pEnhancer overlapped with DA-cCREs
gene3 <- rownames(IntegDACRE_pEPN)[IntegDACRE_pEPN$Ep34 >= 5]
gene4 <- rownames(IntegDACRE_pEPN)[IntegDACRE_pEPN$Ep35 >= 5]
intersect(gene3, gene4)
VennDiagram::venn.diagram(list(Ep34 = gene3, Ep35 = gene4), filename="output/Plots/08_Venn_Ep34vsEp35.tiff", scaled=T, imagetype = "tiff")

promoters2 <- paste0("Promoter_", intersect(gene3, gene4))
df4 <- lapply(promoters2, function(x){
  gr1 <- subsetByOverlaps(IntegDACRE_pEPN_gr$Ep34, gr[pEPN_final[[x]]])
  gr2 <- subsetByOverlaps(IntegDACRE_pEPN_gr$Ep35, gr[pEPN_final[[x]]])
  idx1 <- findOverlaps(gr1, gr2) %>% queryHits() %>% length
  idx2 <- length(gr1) - idx1
  idx3 <- length(gr2) - idx1
  out <- c("Overlap" = idx1, "Ep34_only" = idx2, "Ep35_only" = idx3)
})
df4 <- do.call(rbind, df4)
rownames(df4) <- intersect(gene3, gene4)
write.csv(df4, "output/Tables/08_pEnhancerNo_Ep34vsEp35_v2.csv")

df4 <- reshape2::melt(df4)
p5 <- ggplot(df4, aes(x = reorder(Var1, -value), y = value, fill = Var2)) + geom_bar(stat = "identity") + theme_ArchR() +
  scale_fill_manual(values = c("Overlap" = "darkgray", "Ep34_only" = "darkblue", "Ep35_only" = "darkred")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("output/Plots/08_Barplot_pEnhancer_Ep34vsEp35_v2.pdf", width = 12, height = 4)
p5
dev.off()

#example, PTK6
dt1 <- data.frame(seqnames = "chr20",
                  start = end(geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == "PTK6")]),
                  end = resize(gr[pEPN_final$Promoter_PTK6], width = 1, fix = "center") %>% start)
idx1 <- which(dt1$start > dt1$end)
start1 <- dt1[idx1,]$end
end1 <- dt1[idx1,]$start
dt1[idx1,]$start <- start1
dt1[idx1,]$end <- end1

gr1 <- subsetByOverlaps(IntegDACRE_pEPN_gr$Ep34, gr[pEPN_final$Promoter_PTK6])
gr2 <- subsetByOverlaps(IntegDACRE_pEPN_gr$Ep35, gr[pEPN_final$Promoter_PTK6])
gr3 <- subsetByOverlaps(gr1, gr2)
gr1 <- subsetByOverlaps(gr1, gr3, invert = T)
gr2 <- subsetByOverlaps(gr2, gr3, invert = T)

gr_loop2 <- GRanges(seqnames = dt1$seqnames, IRanges(start = dt1$start, end = dt1$end)) # 48 enhancers
p6.1 <- plotBrowserTrack(
  ArchRProj = arc,
  groupBy = "Clusters2",
  useGroups = c("Ep34", "Ep35"),
  region = GRanges(seqnames = "chr20", ranges = IRanges(start = 61680001, end = 62640000)),
  features = gr3,
  loops = gr_loop2,
  pal = c("darkblue", "darkred"))
grid::grid.newpage()
grid::grid.draw(p6.1)

p6.2 <- plotBrowserTrack(
  ArchRProj = arc,
  groupBy = "Clusters2",
  useGroups = c("Ep34", "Ep35"),
  region = GRanges(seqnames = "chr20", ranges = IRanges(start = 61680001, end = 62640000)),
  features = gr1,
  loops = gr_loop2,
  pal = c("darkblue", "darkred"))
grid::grid.newpage()
grid::grid.draw(p6.2)

p6.3 <- plotBrowserTrack(
  ArchRProj = arc,
  groupBy = "Clusters2",
  useGroups = c("Ep34", "Ep35"),
  region = GRanges(seqnames = "chr20", ranges = IRanges(start = 61680001, end = 62640000)),
  features = gr2,
  loops = gr_loop2,
  pal = c("darkblue", "darkred"))
grid::grid.newpage()
grid::grid.draw(p6.3)

pdf("output/Plots/08_GenomeTrack_PTK6.pdf", width = 10, height = 4)
grid::grid.newpage()
grid::grid.draw(p6.1)
grid::grid.newpage()
grid::grid.draw(p6.2)
grid::grid.newpage()
grid::grid.draw(p6.3)
dev.off()

gr_446prom <- gr[pEPN_final[paste0("Promoter_", intersect(gene3, gene4))] %>% unlist %>% unique] #5576 cCREs
gr1 <- subsetByOverlaps(IntegDACRE_pEPN_gr$Ep34, gr_446prom)
gr2 <- subsetByOverlaps(IntegDACRE_pEPN_gr$Ep35, gr_446prom)
gr3 <- subsetByOverlaps(gr1, gr2)
gr1 <- subsetByOverlaps(gr1, gr3, invert = T)
gr2 <- subsetByOverlaps(gr2, gr3, invert = T)
gr_out <- list(both = gr3, Ep34 = gr1, Ep35 = gr2)

lapply(names(gr_out), function(x){
  gr <- gr_out[[x]] %>% sort
  dt <- data.frame(seqnames = seqnames(gr), start = start(gr)-1, end = end(gr))
  write.table(dt, paste0("output/output_bed/08_Ovlp_promoters_pEnhancers_", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
})
write.table(data.frame(gene = intersect(gene3, gene4)), "output/Tables/08_promoters_ovlp_pEnhancer5_Ep34vsEp35.txt", quote = F, col.names = F, row.names = F)

df <- lapply(gr_out, length) %>% unlist
df <- reshape2::melt(df)
df$name <- rownames(df)
p6.4 <- ggplot(df, aes(x = "x", y = value, fill = name)) +
  geom_bar(stat = "identity", position = "stack", color = "white") + 
  coord_polar(theta = "y") + theme_ArchR() +
  scale_fill_manual(values = c("both" = "darkgray", "Ep34" = "darkblue", "Ep35" = "darkred"))
pdf("output/Plots/08_Pie_Ep34vsEp35_pEnhancer.pdf", width = 5, height = 5)
p6.4
dev.off()

df1 <- fread("output/homer_motif/08_Ovlp_promoters_pEnhancers_Ep34vsEp35/Ep34/knownResults.txt")[,c(1,4)]
df2 <- fread("output/homer_motif/08_Ovlp_promoters_pEnhancers_Ep34vsEp35/Ep35/knownResults.txt")[,c(1,4)]
df1 <- df1[order(df1$`Motif Name`),]
df2 <- df2[order(df2$`Motif Name`),]
df <- cbind(df1, df2)
colnames(df) <- c("Motif", "Ep34", "Motif2", "Ep35")
df <- data.frame(Motif = str_split(df$Motif, "\\(", simplify = T)[,1], Motif2 = df$Motif, Ep34 = -(df$Ep34 / 2.303), Ep35 = -(df$Ep35 / 2.303))

p6.5 <- ggplot(df, aes(x = Ep34, y = Ep35, label = Motif)) + geom_point() + geom_text_repel(size=3) + theme_ArchR() +
  labs(x = "Motif enrichment [-log10(P-value)] pEnhancer only overlapped with Ep34 DA-cCREs", 
       y = "Motif enrichment [-log10(P-value)] pEnhancer only overlapped with Ep35 DA-cCREs") +
  ggtitle("Homer TF motifs")

pdf("output/Plots/08_ScatterPlot_MotifEnrichment_pEnhancer_Ep34vsEp35.pdf", width = 6, height = 6)
p6.5
dev.off()
write.csv(df, "output/Tables/08_ScatterPlot_MotifEnrichment_pEnhancer_Ep34vsEp35.csv")

arc <- addMotifAnnotations(arc, motifSet = "homer", annoName = "homer")
arc <- addBgdPeaks(arc) %>% addDeviationsMatrix(., peakAnnotation = "homer")
homerMarkerEp34Ep35 <- getMarkerFeatures(arc, groupBy = "Clusters2", useGroups = "Ep35", bgdGroups = "Ep34", useMatrix = "homerMatrix", useSeqnames = "z")

df <- data.frame(Motif = rowData(homerMarkerEp34Ep35)$name, 
                 FDR = as.numeric(assays(homerMarkerEp34Ep35)$FDR[,1]) %>% -log10(.),
                 MeanDiff = as.numeric(assays(homerMarkerEp34Ep35)$MeanDiff[,1]))
p6.6 <- ggplot(df, aes(x = MeanDiff, y = FDR, label = Motif)) + geom_point() + geom_text_repel(size=3) + theme_ArchR() +
  labs(x = "MeanDiff [Ep35 - Ep34]", y = "-log10(FDR)") + geom_vline(xintercept = 0, lty = "dotted") + ggtitle("Homer TF motifs")
pdf("output/Plots/08_VolcanoPlot_MotifScore_Ep34vsEp35.pdf", width = 6, height = 6)
p6.6
dev.off()
write.csv(df, "output/Tables/08_VolcanoPlot_MotifScore_Ep34vsEp35.csv")

#----- Inter-tumor heterogeneity ER+/HER2- -----#
ER_clusters <- cluster_info$Cluster[cluster_info$Subtype == "ER+/HER2–"]
IntegDACRE_pEPN_ER <- IntegDACRE_pEPN[,ER_clusters] %>% as.matrix() #17 clusters
#df5 <- data.frame(promoter = rownames(IntegDACRE_pEPN_ER), variance = rowVars(IntegDACRE_pEPN_ER), median = rowMedians(IntegDACRE_pEPN_ER))
#write.csv(df5, "output/Tables/08_ER_pEPN_profile.csv")

IntegDACRE_pEPN_ER2 <- IntegDACRE_pEPN_ER >= 5
IntegDACRE_pEPN_ER3 <- lapply(c(1:nrow(IntegDACRE_pEPN_ER2)), function(i) which(IntegDACRE_pEPN_ER2[i,] == TRUE) %>% length()) %>% unlist
names(IntegDACRE_pEPN_ER3) <- rownames(IntegDACRE_pEPN_ER2)
gene5 <- which(IntegDACRE_pEPN_ER3 >= 6) %>% names
write.table(data.frame(gene = gene5), "output/Tables/08_promoters_ovlp_pEnhancer5_ERclusters.txt", quote = F, col.names = F, row.names = F)

IntegDACRE_pEPN_gr <- readRDS("rds/06_IntegDACRE_pEPN_gr.rds")

df6 <- mclapply(names(pEPN_final), function(y){
  out <- lapply(ER_clusters, function(x) countOverlaps(gr[pEPN_final[[y]]], IntegDACRE_pEPN_gr[[x]]))
  out <- do.call(rbind, out)
  rownames(out) <- ER_clusters
  colnames(out) <- paste0("pEnhancer#", c(1:ncol(out)))
  return(out)
}, mc.cores = 24)
names(df6) <- names(pEPN_final)
saveRDS(df6, "rds/08_pEnhancer_DAcCRE_ERclusters_TableList.rds")
# df6 <- readRDS("rds/08_pEnhancer_DAcCRE_ERclusters_TableList.rds")

df6$Promoter_ZNF703
col_fun1 <- colorRamp2(c(0,1), colors = c("lightgray", "darkorchid4"))
ht1 <- Heatmap(df6$Promoter_ZNF703, cluster_rows = F, cluster_columns = F, col = col_fun1, rect_gp = gpar(col = "black", lwd = 0.5), column_title = "ZNF703")
p7 <- draw(ht1)
pdf("output/Plots/08_Heatmap_pEnhancer_ERclusters_ZNF703.pdf", width = 8, height = 4)
p7
dev.off()

ht2 <- Heatmap(df6$Promoter_TFF1, cluster_rows = F, cluster_columns = F, col = col_fun1, 
               rect_gp = gpar(col = "black", lwd = 0.5), column_title = "TFF1")
p8 <- draw(ht2)
pdf("output/Plots/08_Heatmap_pEnhancer_ERclusters_TFF1.pdf", width = 8, height = 4)
p8
dev.off()

arc <- addModuleScore(arc, useMatrix = "GeneScoreMatrix", name = "Module1", features = list(ERgene = intersect(getFeatures(arc), gene5)))
df7 <- data.frame(Module = arc$Module1.ERgene, Cluster = arc$Clusters2)
df7 <- df7[which(df7$Cluster %in% c(ER_clusters, paste0("Ep", c(3,4,26,24,5,6,1,2)))),]
df7$Cluster <- factor(df7$Cluster, levels = c(ER_clusters, paste0("Ep", c(3,4,26,24,5,6,1,2))))

p9 <- ggplot(df7, aes(x = Cluster, y = Module)) + geom_violin(fill = "lightgray") + geom_boxplot(outlier.size = 0, alpha = 0.5) +
  theme_ArchR() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("output/Plots/08_Vlnplot_ERclusterGenes_GSModule.pdf", width = 8, height = 4)
p9
dev.off()

#example, ZNF703
dt1 <- data.frame(seqnames = "chr8",
                  start = start(geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == "ZNF703")]),
                  end = resize(gr[pEPN_final$Promoter_ZNF703], width = 1, fix = "center") %>% start)
idx1 <- which(dt1$start > dt1$end)
start1 <- dt1[idx1,]$end
end1 <- dt1[idx1,]$start
dt1[idx1,]$start <- start1
dt1[idx1,]$end <- end1

gr_loop3 <- GRanges(seqnames = dt1$seqnames, IRanges(start = dt1$start, end = dt1$end)) # 49 enhancers
p10 <- plotBrowserTrack(
  ArchRProj = arc,
  groupBy = "Clusters2",
  useGroups = cluster_info$Cluster[cluster_info$Subtype == "ER+/HER2–"],
  region = GRanges(seqnames = "chr8", ranges = IRanges(start = 37150001, end = 38030000)),
  loops = gr_loop3)
grid::grid.newpage()
grid::grid.draw(p10)

pdf("output/Plots/08_GenomeTrack_ZNF703.pdf", width = 10, height = 6)
grid::grid.newpage()
grid::grid.draw(p10)
dev.off()

#----- Summarize motif analysis -----#
dirlist <- cluster_info$Cluster[cluster_info$Subtype == "ER+/HER2–"]
dirlist <- dirlist[which(dirlist != "Ep22")]
motif_DF <- lapply(dirlist, function(x){
  out <- data.table::fread(paste0("output/homer_motif/08_ER_clusters/", x, "/knownResults.txt"), header = F, skip = 1)[, c(1,3,7)] %>% data.frame %>% `colnames<-`(., c("Motif", "P", "TgtFreq"))
  out$mlog10P <- -log10(as.numeric(out$P))
  out$mlog10P[which(out$mlog10P == Inf)] <- as.numeric(gsub("1e-", "", out$P[which(out$mlog10P == Inf)]))
  out$TgtFreq <- gsub("%", "", out$TgtFreq) %>% as.numeric()
  out <- out[order(out$Motif), ]
  return(out)
}) %>% `names<-`(., dirlist)

motif_p <- lapply(motif_DF, function(x) x$mlog10P)
motif_p <- do.call(cbind, motif_p) %>% `rownames<-`(., motif_DF[[1]]$Motif)
motif_p_scale <- apply(motif_p, 2, function(x) scales::rescale(x, to = c(0,100)))
idy <- lapply(c(1:ncol(motif_p_scale)), function(i) order(motif_p_scale[,i], decreasing = T)[c(1:10)]) %>% unlist %>% unique

n1 <- str_split(rownames(motif_p_scale), "/", simplify = T)[,1]
n2 <- str_split(rownames(motif_p_scale), "/", simplify = T)[,2]
n2 <- str_split(n2, "-", simplify = T)[,1]
rownames(motif_p_scale) <- paste0(n1, "/", n2)

col_fun2 <- colorRamp2(c(0,10,20,30,40,50,75,100), colors = c("white", rev(viridis(n = 7, option = "G"))))
fh <- function(x) hclust(dist(x), method = "ward.D2")
ht3 <- Heatmap(t(motif_p_scale[idy,]), cluster_rows = fh, cluster_columns =fh, col = col_fun2, name = "Norm.Enrichment -log10(P-value) [0-MAX]",
               show_row_dend = T, show_column_dend = T, show_row_names = T, 
               row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
               use_raster = T, raster_by_magick = F)
p12 <- draw(ht3)
pdf("output/Plots/08_Heatmap_homerMotif_NormPvalue.pdf", height = 6, width = 8)
p12
dev.off()
