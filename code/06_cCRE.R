#----------------------------------------------------------------------------
# 06_cCRE.R
#----------------------------------------------------------------------------
library(ArchR)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(ggrepel)
library(scales)
library(ggrastr)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
addArchRThreads(threads = 20)
addArchRGenome("hg19")

pt_ls <- readRDS("rds/01_pt_ls_2.rds")
arc <- readRDS("rds/02_arc.rds")
arc_ep <- readRDS("rds/03_arc_ep.rds")

#----- cluster re-assignment ----#
df1 <- data.frame(row.names = arc$cellNames, Cluster = arc$Clusters, Cluster2 = arc$Clusters)
df2 <- data.frame(row.names = arc_ep$cellNames, EPCluster = arc_ep$Clusters)
df2$EPCluster <- gsub("C", "EP", df2$EPCluster)
df1[rownames(df2), "Cluster2"] <- df2$EPCluster

arc$Clusters2 <- df1$Cluster2
p1 <- plotEmbedding(arc, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP", plotAs = "points", size = 0.5)
plotPDF(p1, name = "06r_UMAP_Cluster2.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#----- peak call -----#
pathToMacs2 <- findMacs2()
arc <- addGroupCoverages(arc, groupBy = "Clusters2", force = T)
arc <- addReproduciblePeakSet(arc, groupBy = "Clusters2", pathToMacs2 = pathToMacs2, method = "q", cutOff = 0.1, force = T)
arc <- addPeakMatrix(arc, force = T)

saveRDS(arc, "rds/06r_arc.rds")

#----- cCRE genomic positions -----#
gr <- arc@peakSet
names(gr) <- c(paste0("cCRE_", c(1:229238)))
 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak(gr, tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
df <- peakAnno@annoStat
 
pie_colors <- ArchRPalettes$stallion %>% `names<-`(., df$Feature)
p1 <- ggplot(df, aes(x = "x", y = Frequency, fill = Feature)) +
   geom_bar(stat = "identity", position = "stack", color = "white") + 
   coord_polar(theta = "y") + theme_ArchR() + scale_fill_manual(values = pie_colors)
pdf("output/Plots/06r_cCRE_all_annotation.pdf", width = 5, height = 5)
p1
plotAnnoPie(peakAnno)
dev.off()

#----- modules -----#
# arc <- readRDS("rds/06r_arc.rds")
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
rownames(mtx) <- names(gr)

mtx <- mtx[which(rowSums(mtx) != 0),]
mtx <- t(scale(t(mtx)))

km2 <- kmeans(mtx, centers = 26, iter.max = 20)
km_clusters <- km2$cluster
saveRDS(km2, "rds/06-R2_km2.rds")

col_fun1 <- colorRamp2(c(-2,-1,0,1,2,3,4), colors = viridis(7, option = "D"))
ht1 <- Heatmap(mtx, cluster_rows = F, cluster_columns = F, col = col_fun1, name = "Accessibility z-score",
               show_row_dend = F, show_column_dend = F, show_row_names = F, column_names_gp = gpar(fontsize = 6),
               row_split = km_clusters, use_raster = T, raster_by_magick = F)
p8 <- draw(ht1)
pdf("output/Plots/06-R2_Heatmap_kmeans25.pdf", width = 6, height = 8)
p8
dev.off()

df <- data.frame(module = c(1:26), first_anno = lapply(c(1:26), function(i) colSums(mtx[which(km_clusters == i), ]) %>% sort(., decreasing = T) %>% .[1] %>% names) %>% unlist)
df$final_anno <- df$first_anno
df$final_anno[c(5,11,15,17,19,20,21)] <- c("Plasma", "Endothelial", "B", "Fibroblast", "Myeloid", "T-C6", "T")
df$final_anno[9] <- "EP24-EP25"
df$final_anno[10] <- "ER"
# colSums(mtx[which(km_clusters == 10), ]) %>% sort(., decreasing = T)

gr_modules <- lapply(c(1:26), function(i) gr[names(which(km_clusters == i))])
names(gr_modules) <- df$final_anno

lapply(names(gr_modules), function(i){
  g <- gr_modules[[i]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/06-R2_modules/", i, "_module.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})
saveRDS(gr_modules, "rds/06-R2_gr_modules.rds")

Module_order <- c("T","T-C6","B","Plasma","Fibroblast","Endothelial","Myeloid","ER",
                  "EP1","EP3","EP4","EP5","EP6","EP8","EP9","EP11","EP16","EP18","EP19","EP21","EP22","EP23","EP24-EP25","EP26","EP29","EP30")  

km_clusters_rename <- km_clusters
for(i in c(1:26)){km_clusters_rename[km_clusters == i] <- df$final_anno[i]}

ht2 <- Heatmap(mtx, cluster_rows = F, cluster_columns = F, col = col_fun1, name = "Accessibility z-score",
               show_row_dend = F, show_column_dend = F, show_row_names = F, column_names_gp = gpar(fontsize = 6),
               row_split = factor(km_clusters_rename, levels = Module_order), 
               row_title_gp = gpar(fontsize = 6), row_title_rot = 0,
               use_raster = T, raster_by_magick = F)
p9 <- draw(ht2)
pdf("output/Plots/06-R2_Heatmap_kmeans25_ordered.pdf", width = 6, height = 8)
p9
dev.off()

colnames(mtx)
df_ha <- data.frame(clusters = colnames(mtx),
                    mixtype = c(rep(NA,13), rep("less-mixed", 7), rep("moderately-mixed", 2), rep("less-mixed", 2), "moderately-mixed", "highly-mixed",
                                rep("less-mixed", 4), rep(c("moderately-mixed", "less-mixed"), 3), "moderately-mixed", "highly-mixed",
                                rep("less-mixed", 5), "highly-mixed", rep("less-mixed", 4)),
           samples = c(rep(NA,13),
                       "P190","P190","P41","P44","P122","P166","P189","P209PE","P200A","P210",
                       "P180","P178",NA,"P211LN","P185","P165","P179",NA,"P175PE","P191",
                       "P93","P181","P20LN",NA,NA,"P51","P64","P33","P35","P39",NA,"P40","P49","P50","P50"))
rownames(pt_ls) <- pt_ls$ID
df_ha$Type <- pt_ls[df_ha$samples,]$Type
df_ha$receptor <- pt_ls[df_ha$samples,]$receptor_status

sample_colors <- readRDS("rds/02_sample_colors.rds")
sampletype_colors <- readRDS("rds/02_sampletype_colors.rds")
subtype_colors <- c("ER+/HER2-" = "blue", "ER-/PGR+/HER-" = "green4", "ER+/HER2+" = "yellow2", "HER2+" = "orange", "TN" = "red", "NA" = "gray")

ha1 <- HeatmapAnnotation(Sample = df_ha$samples, SampleType = df_ha$Type, Receptor = df_ha$receptor, Mix = df_ha$mixtype,
                         col = list(Sample = sample_colors, SampleType = sampletype_colors, Receptor = subtype_colors, Mix = ArchRPalettes$kelly[c(1:3)] %>% `names<-`(., c("moderately-mixed", "less-mixed", "highly-mixed"))))
ht3 <- Heatmap(mtx, cluster_rows = F, cluster_columns = F, col = col_fun1, name = "Accessibility z-score",
               show_row_dend = F, show_column_dend = F, show_row_names = F, column_names_gp = gpar(fontsize = 6),
               row_split = factor(km_clusters_rename, levels = Module_order), 
               row_title_gp = gpar(fontsize = 6), row_title_rot = 0, top_annotation = ha1,
               use_raster = T, raster_by_magick = F)
p9.2 <- draw(ht3)
pdf("output/Plots/06-R2_Heatmap_kmeans25_ordered_ha1.pdf", width = 7, height = 8)
p9.2
dev.off()

ht4 <- Heatmap(mtx, cluster_rows = F, cluster_columns = F, col = col_fun1, name = "Accessibility z-score",
               show_row_dend = F, show_column_dend = F, show_row_names = F, column_names_gp = gpar(fontsize = 6),
               row_split = factor(km_clusters_rename, levels = Module_order), 
               row_title = df[Module_order,]$module,
               row_title_gp = gpar(fontsize = 6), row_title_rot = 0, top_annotation = ha1,
               use_raster = T, raster_by_magick = F)
p9.3 <- draw(ht4)
pdf("output/Plots/06-R2_Heatmap_kmeans25_ordered_ha1_kmeansno.pdf", width = 7, height = 8)
p9.3
dev.off()

pred_sample <- c(rep(NA,13),
  "P190","P190","P41","P44","P122","P166","P189","P209PE","P200A","P210",
  "P180","P178",NA,"P211LN","P185","P165","P179",NA,"P175PE","P191",
  "P93","P181","P20LN",NA,NA,"P51","P64","P33","P35","P39",NA,"P40","P49","P50","P50")
rownames(pt_ls) <- pt_ls$ID
pred_samplesub <- pt_ls[pred_sample,"receptor_status"]
names(pred_samplesub) <- colnames(mtx)

mtx_rm10 <- mtx[names(which(km2$cluster == "10")),]
mtx_rm10 <- reshape2::melt(mtx_rm10)
mtx_rm10$Subtype <- pred_samplesub[mtx_rm10$Var2]

p9.4 <- ggplot(mtx_rm10, aes(x = Var2, y = value, fill = Subtype)) + geom_violin(alpha = 0.5) + geom_boxplot(outlier.size = 0, alpha = 0) +
  scale_fill_manual(values = subtype_colors) + theme_ArchR() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("output/Plots/06-R2_Boxplot_rm10_zscore.pdf", width = 10, height = 6)
p9.4
dev.off()

#----- module annotation -----#
ModuleAnno <- mclapply(gr_modules, function(x) annotatePeak(x, tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db"), mc.cores = 1)

Anno_df <- lapply(names(gr_modules), function(i){
  tmp <- ModuleAnno[[i]]@annoStat$Frequency * length(ModuleAnno[[i]]@anno) * 0.01
  out <- data.frame(Feature = ModuleAnno[[i]]@annoStat$Feature, PeakNo = tmp, Sample = i)
  return(out)
})
Anno_df <- do.call(rbind, Anno_df)
p10 <- ggplot(Anno_df, aes(x = factor(Sample, levels = as.character(Module_order)), y = PeakNo, fill = Feature)) + geom_bar(stat = "identity") +
  theme_ArchR() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = pie_colors)
pdf("output/Plots/06-R2_Barplot_ModuleNumAnno.pdf", width = 5, height = 5)
p10
dev.off()

p10.2 <- lapply(Module_order, function(i) ggplot(Anno_df[which(Anno_df$Sample == i),], aes(x = "x", y = PeakNo, fill = Feature)) +
  geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y") + theme_ArchR() + scale_fill_manual(values = pie_colors) + ggtitle(i))
pdf("output/Plots/06-R2_Pie_ModuleNumAnno.pdf", width = 5, height = 5)
p10.2
dev.off()

chr_df <- lapply(gr_modules, function(x) table(seqnames(x))) %>% do.call(cbind, .)
chr_df <- chr_df[,Module_order]
# colSums(chr_df)
chr_df <- reshape2::melt(chr_df)
chr_colors <- c(ArchRPalettes$stallion2, ArchRPalettes$calm)[c(1:23)] %>% `names<-`(., paste0("chr", c(1:22,"X")))
p10.3 <- lapply(Module_order, function(i) ggplot(chr_df[which(chr_df$Var2 == i),], aes(x = "x", y = value, fill = Var1)) +
                  geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y") + theme_ArchR() + scale_fill_manual(values = chr_colors) + ggtitle(i))
pdf("output/Plots/06-R2_Pie_ModuleNumChr.pdf", width = 5, height = 5)
p10.3
dev.off()

p10.4 <- ggplot(chr_df, aes(x = factor(Var2, levels = as.character(Module_order)), y = value, fill = Var1)) + geom_bar(stat = "identity") +
  theme_ArchR() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = chr_colors)
pdf("output/Plots/06-R2_Barplot_ModuleNumChr.pdf", width = 5, height = 5)
p10.4
dev.off()

#----- Summarize motif analysis -----#
dirlist <- list.dirs("output/homer_motif/06-R2_modules/", recursive = F, full.names = F)
motif_DF <- lapply(dirlist, function(x){
  out <- data.table::fread(paste0("output/homer_motif/06-R2_modules/", x, "/knownResults.txt"), header = F, skip = 1)[, c(1,3,7)] %>% data.frame %>% `colnames<-`(., c("Motif", "P", "TgtFreq"))
  out$mlog10P <- -log10(as.numeric(out$P))
  out$mlog10P[which(out$mlog10P == Inf)] <- as.numeric(gsub("1e-", "", out$P[which(out$mlog10P == Inf)]))
  out$TgtFreq <- gsub("%", "", out$TgtFreq) %>% as.numeric()
  out <- out[order(out$Motif), ]
  return(out)
}) %>% `names<-`(., dirlist)

motif_p <- lapply(motif_DF, function(x) x$mlog10P)
motif_p <- do.call(cbind, motif_p) %>% `rownames<-`(., motif_DF[[1]]$Motif)
motif_p <- motif_p[,as.character(Module_order)]

motif_p
motif_p_scale <- apply(motif_p, 2, function(x) scales::rescale(x, to = c(0,100)))
idy <- lapply(c(1:ncol(motif_p_scale)), function(i) order(motif_p_scale[,i], decreasing = T)[c(1:10)]) %>% unlist %>% unique

n1 <- str_split(rownames(motif_p_scale), "/", simplify = T)[,1]
n2 <- str_split(rownames(motif_p_scale), "/", simplify = T)[,2]
n2 <- str_split(n2, "-", simplify = T)[,1]
rownames(motif_p_scale) <- paste0(n1, "/", n2)

col_fun2 <- colorRamp2(c(0,10,20,30,40,50,75,100), colors = c("white", rev(viridis(n = 7, option = "G"))))
fh <- function(x) hclust(dist(x), method = "ward.D2")
ht3 <- Heatmap(motif_p_scale[idy,], cluster_rows = fh, cluster_columns =F, col = col_fun2, name = "Norm.Enrichment -log10(P-value) [0-MAX]",
               show_row_dend = T, show_column_dend = T, show_row_names = T, 
               row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
               use_raster = T, raster_by_magick = F)
p11 <- draw(ht3)
pdf("output/Plots/06-R2_HomerMotif_Heatmap_normPvalue_v1.pdf", height = 7, width = 8)
p11
dev.off()

idy <- lapply(c(1:ncol(motif_p_scale)), function(i) order(motif_p_scale[,i], decreasing = T)[c(1:10)]) %>% unlist %>% unique
rownames(motif_p_scale)

motif_of_interest <- c("ERE(NR),IR3/MCF7",
                       "FOXA1(Forkhead)/MCF7",
                       "GATA3(Zf)/iTreg",
                       "ARE(NR)/LNCAP",
                       "CTCF(Zf)/CD4+",
                       "ETS1(ETS)/Jurkat",
                       "ETS:RUNX(ETS,Runt)/Jurkat",
                       "RUNX(Runt)/HPC7",
                       "IRF8(IRF)/BMDM",
                       "EWS:FLI1-fusion(ETS)/SK_N_MC",
                       "Tcf21(bHLH)/ArterySmoothMuscle",
                       "TEAD(TEA)/Fibroblast",
                       "AP-1(bZIP)/ThioMac",
                       "AP-2alpha(AP2)/Hela",
                       "Sox3(HMG)/NPC",
                       "KLF3(Zf)/MEF",
                       "NF1(CTF)/LNCAP")

ht4 <- Heatmap(motif_p_scale[motif_of_interest,], cluster_rows = F, cluster_columns = F, col = col_fun2, 
               name = "Norm.Enrichment -log10(P-value) [0-MAX]",
               show_row_dend = T, show_column_dend = T, show_row_names = T, 
               row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
               use_raster = T, raster_by_magick = F)
p11.2 <- draw(ht4)
pdf("output/Plots/06-R2_HomerMotif_Heatmap_normPvalue_v2.pdf", height = 4, width = 9)
p11.2
dev.off()
