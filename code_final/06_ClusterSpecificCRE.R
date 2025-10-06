#------------------------------------------------------------------------------
# 06_ClusterSpecificCRE.R
#------------------------------------------------------------------------------
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(karyoploteR)
library(ChIPseeker)
library(ggrastr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(ggrepel)
addArchRThreads(threads = 24)
addArchRGenome("hg19")

arc <- readRDS("rds/05_arc.rds")
gr_valid <- readRDS("rds/05_gr_valid.rds")
pt_ls <- read.csv("ref/patient_list_v2.csv")

#----- Cluster-specific CREs -----#
ClusterPeaks <- getMarkerFeatures(
  ArchRProj = arc, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(ClusterPeaks, "rds/06_ClusterPeaks.rds")

ClusterPeaksList <- getMarkers(ClusterPeaks, cutOff = "FDR < 0.25 & Log2FC > 0", returnGR = T)
ClusterPeaksList <- lapply(ClusterPeaksList, function(x) subsetByOverlaps(x, gr_valid))
ClusterPeaksList <- ClusterPeaksList[lapply(ClusterPeaksList, length) %>% unlist != 0]
# Ep13 Ep19 : without specific CREs

#----- Cluster-specific CREs -----#
amp_gr <- readRDS("rds/04_gr_Epi_ampWindow.rds")
ClusterSampleCorrespond <- c("Ep9" = "P200A", 
                             "Ep8" = "P209PE", 
                             "Ep26" = "P51", 
                             "Ep2" = "P44", 
                             "Ep1" = "P41", 
                             "Ep6" = "P190", 
                             "Ep5" = "P190", 
                             "Ep11" = "Mixed", 
                             "Ep10" = "Mixed", 
                             "Ep7" = "P20LN", 
                             # "Ep19" = "Mixed", 
                             "Ep25" = "P181", 
                             "Ep24" = "P93", 
                             # "Ep13" = "P175PE", 
                             "Ep14" = "P189", 
                             "Ep20" = "P185", 
                             "Ep17" = "Mixed", 
                             "Ep15" = "P191", 
                             "Ep33" = "P49", 
                             "Ep21" = "P179", 
                             "Ep28" = "P33", 
                             "Ep18" = "P211LN", 
                             "Ep16" = "P178", 
                             "Ep12" = "P180", 
                             "Ep31" = "Mixed", 
                             "Ep30" = "P39", 
                             "Ep32" = "P40", 
                             "Ep35" = "P50", 
                             "Ep34" = "P50", 
                             "Ep29" = "P35", 
                             "Ep22" = "P165", 
                             "Ep23" = "P210LR", 
                             "Ep27" = "P64", 
                             "Ep4" = "P166", 
                             "Ep3" = "P122")

ClusterPeaksList_AmpOverlap <- lapply(names(ClusterSampleCorrespond), function(x){
  gr <- ClusterPeaksList[[x]]
  n <- ClusterSampleCorrespond[x]
  out <- subsetByOverlaps(gr, amp_gr[amp_gr$Sample == n])
  return(out)
})
names(ClusterPeaksList_AmpOverlap) <- names(ClusterSampleCorrespond)
ClusterPeaksList_AmpOverlap <- ClusterPeaksList_AmpOverlap[unlist(lapply(ClusterPeaksList_AmpOverlap, length)) !=0]
ClusterPeaksList_AmpOverlap

ClusterPeaksList <- lapply(ClusterPeaksList, function(x){
  x$Amp <- "N"
  return(x)
})

for(i in names(ClusterPeaksList_AmpOverlap)){
  fo1 <- findOverlaps(ClusterPeaksList[[i]], ClusterPeaksList_AmpOverlap[[i]])
  ClusterPeaksList[[i]]$Amp[queryHits(fo1)] <- "Y"
}

NumCRETb <- lapply(names(ClusterPeaksList), function(x){
  out <- table(ClusterPeaksList[[x]]$Amp) %>% as.data.frame()
  out$Var1 <- as.character(out$Var1)
  out$Cluster <- x
  if(length(which(out$Var1 == "Y")) == 0){
    out <- rbind(out, c("Y", 0, x))
  }
  return(out)
}) %>% do.call(rbind, .)
df_tmp <- data.frame(Var1 = c("N", "Y", "N", "Y"), Freq = c(0,0,0,0), Cluster = c("Ep13","Ep13","Ep19","Ep19"))
NumCRETb <- rbind(NumCRETb, df_tmp)

NumCRETb$Freq <- as.numeric(NumCRETb$Freq)
NumCRETb$Cluster <- factor(NumCRETb$Cluster, levels = colnames(ClusterPeaks))
NumCRETb$Var1 <- factor(NumCRETb$Var1, levels = c("Y", "N"))

p1 <- ggplot(NumCRETb, aes(x = Cluster, y = Freq, fill = Var1)) +　geom_bar(stat = "identity") + theme_ArchR() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("Y" = "red", "N" = "darkgray")) +
  labs(y = "# of cluster-specific cCREs (FDR < 25%)")
pdf("output/Plots/06_Barplot_NumClusterSpecificCRE.pdf", width = 8, height = 4)
p1
dev.off()
write.csv(NumCRETb, "output/Tables/06_Barplot_NumClusterSpecificCRE.csv")

gr1 <- ClusterPeaksList$Ep3[which(ClusterPeaksList$Ep3$Amp == "N")]
gr2 <- ClusterPeaksList$Ep3[which(ClusterPeaksList$Ep3$Amp == "Y")]
gr3 <- geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == "ERBB2")]
pdf("output/Plots/06_Karyotype_Ep3_cCRE_Amp.pdf", width = 6, height = 6)
kp <- plotKaryotype(genome = "hg19", chromosomes = paste0("chr", c(1:22, "X")))
kpPlotRegions(kp, data=gr1)
kpPlotRegions(kp, data=gr2, col = "red")
kpPlotMarkers(kp, chr = seqnames(gr3), x = start(gr3), labels = gr3$symbol)
dev.off()

gr1 <- ClusterPeaksList$Ep28[which(ClusterPeaksList$Ep28$Amp == "N")]
gr2 <- ClusterPeaksList$Ep28[which(ClusterPeaksList$Ep28$Amp == "Y")]
pdf("output/Plots/06_Karyotype_Ep28_cCRE_Amp.pdf", width = 6, height = 6)
kp <- plotKaryotype(genome = "hg19", chromosomes = paste0("chr", c(1:22, "X")))
kpPlotRegions(kp, data=gr1)
kpPlotRegions(kp, data=gr2, col = "red")
dev.off()

gr1 <- ClusterPeaksList$Ep27[which(ClusterPeaksList$Ep27$Amp == "N")]
gr2 <- ClusterPeaksList$Ep27[which(ClusterPeaksList$Ep27$Amp == "Y")]
pdf("output/Plots/06_Karyotype_Ep27_cCRE_Amp.pdf", width = 6, height = 6)
kp <- plotKaryotype(genome = "hg19", chromosomes = paste0("chr", c(1:22, "X")))
kpPlotRegions(kp, data=gr1)
kpPlotRegions(kp, data=gr2, col = "red")
dev.off()

gr_tmp <- unlist(GRangesList(ClusterPeaksList))
pdf("output/Plots/06_Karyotype_DAcCRE_CNAoverlap.pdf", width = 6, height = 6)
kp <- plotKaryotype(genome = "hg19", chromosomes = paste0("chr", c(1:22, "X")))
kpPlotRegions(kp, data=amp_gr2, col = "indianred1")
kpPlotRegions(kp, data=gr_tmp)
dev.off()

saveRDS(ClusterPeaksList, "rds/06_ClusterPeaksList.rds")

#----- Enriched motifs for each cCRE -----#
arc <- addMotifAnnotations(ArchRProj = arc, motifSet = "cisbp", annoName = "cisbp")
arc <- addMotifAnnotations(ArchRProj = arc, motifSet = "homer", annoName = "homer")

enrichMotifs <- peakAnnoEnrichment(
  seMarker = ClusterPeaks,
  ArchRProj = arc,
  pe = "cisbp",
  cutOff = "FDR < 0.25 & Log2FC > 0"
)
enrichMotifsHomer <- peakAnnoEnrichment(
  seMarker = ClusterPeaks,
  ArchRProj = arc,
  pe = "homer",
  cutOff = "FDR < 0.25 & Log2FC > 0"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 3, clusterCols = F, transpose = TRUE)
heatmapEM2 <- plotEnrichHeatmap(enrichMotifsHomer, n = 3, clusterCols = F, transpose = TRUE)

p2 <- ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
p3 <- ComplexHeatmap::draw(heatmapEM2, heatmap_legend_side = "bot", annotation_legend_side = "bot")

pdf("output/Plots/06_Heatmap_MotifEnrichment_ClusterSpecific.pdf", width = 6, height = 8)
p2
p3
dev.off()

tmp1 <- list.files("output/PeakCalls/", pattern = "gr.rds")
peaks_clusters <- lapply(tmp1, function(x) readRDS(paste0("output/PeakCalls/", x)))
names(peaks_clusters) <- gsub("-reproduciblePeaks.gr.rds", "", tmp1)
lapply(names(peaks_clusters), function(x){
  g <- peaks_clusters[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/06_PeakClusters/06_PeakClusters_", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
})

tmp2 <- list.files("output/homer_motif/06_PeakClusters/")[-1]
motifDF <- lapply(tmp2, function(x){
  out <- fread(paste0("output/homer_motif/06_PeakClusters/", x, "/knownResults.txt")) %>% data.frame()
  out[,1] <- paste0(stringr::str_split(out[,1], pattern = "/", simplify = T)[,1], "#", 
                    stringr::str_split(stringr::str_split(out[,1], pattern = "\\/", simplify = T)[,2], pattern = "-", simplify = T)[,1])
  out <- out[out$Motif.Name %in% c("ERE(NR),IR3#MCF7","FOXA1(Forkhead)#MCF7"),]
  out <- out[,c(1,3)]
  out$Cluster <- x
  return(out)
})
motifDF <- do.call(rbind,motifDF)
motifDF$P.value <- -log10(as.numeric(motifDF$P.value))
motifDF$P.value[is.infinite(motifDF$P.value)] <- 350

motifDF$Cluster <- factor(motifDF$Cluster, levels = c(paste0("C",c(5:13,15:16)), paste0("Ep", c(1:35))))

p3.2 <- ggplot(motifDF[motifDF$Motif.Name=="FOXA1(Forkhead)#MCF7",], aes(x=Cluster, y=P.value, fill=Motif.Name)) + geom_bar(stat = "identity", position = "dodge") +
  theme_ArchR() + scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y="Motif enrichment [-log10(P-value)]", y="") +
  scale_fill_manual(values = c("FOXA1(Forkhead)#MCF7"="darkred"))
p3.3 <- ggplot(motifDF[motifDF$Motif.Name=="ERE(NR),IR3#MCF7",], aes(x=Cluster, y=P.value, fill=Motif.Name)) + geom_bar(stat = "identity", position = "dodge") +
  theme_ArchR() + scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y="Motif enrichment [-log10(P-value)]", y="") +
  scale_fill_manual(values = c("ERE(NR),IR3#MCF7"="darkblue"))
pdf("output/Plots/06_Barplot_ERFOXA1_MotifEnrichment_allClusterPeak.pdf", width = 8, height = 4)
p3.2
p3.3
dev.off()

#----- co-accessiblity, 500kb -----#
arc <- addCoAccessibility(arc, reducedDims = "IterativeLSI", maxDist = 500000)
saveRDS(arc, "rds/06_arc.rds")
# arc <- readRDS("rds/06_arc.rds")

coA_all <- getCoAccessibility(arc, corCutOff = 0, resolution = 500, returnLoops = FALSE) # 20570714 > 10285357 co-accessible pairs
fo1 <- findOverlaps(arc@peakSet, gr_valid)
coA_all <- coA_all[which(coA_all$queryHits %in% queryHits(fo1) & coA_all$subjectHits %in% queryHits(fo1)),] 
# 20536978 > 10268489 co-accessible pairs, 16868 pairs removed (sparse peaks)

coA_all_corr <- data.frame(corr = coA_all$correlation, FDR = coA_all$FDR)

# coA_all_corr <- coA_all_corr[order(coA_all_corr$FDR),]
# coA_all_corr$corr[max(which(coA_all_corr$FDR < 0.01))]

p4 <- ggplot(coA_all_corr, aes(x=corr)) + geom_density(color = "lightgray", fill="lightgray",alpha=.2) + theme_ArchR() + geom_vline(xintercept = 0.4, lty = "dotted")
pdf("output/Plots/06_DensityPlot_coA.pdf", width = 5, height = 5)
p4
dev.off()

coA_sf <- getCoAccessibility(arc, corCutOff = 0.4, resolution = 500, returnLoops = FALSE) # 1364538 > 682269 co-accessible pairs
coA_sf <- coA_sf[which(coA_sf$queryHits %in% queryHits(fo1) & coA_sf$subjectHits %in% queryHits(fo1)),] # 1361398 > 680699 co-accessible pairs

#not filter sparse peaks (coA_sf is filtered)
gr <- arc@peakSet
names(gr) <- c(paste0("cCRE_", c(1:length(gr))))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak(gr, tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
mcols(gr) <- as.GRanges(peakAnno)

genomic_anno <- gr$annotation
genomic_anno <- str_split(genomic_anno, pattern = "\\(", simplify = T)[,1] %>% gsub(" ", "", .)
table(genomic_anno)
# 3'UTR            5'UTR DistalIntergenic       Downstream             Exon           Intron         Promoter 
# 5279             1214            70297             2734            10970            68661            66132 
gr$annotation2 <- genomic_anno
gr$annotation3 <- paste0(gr$annotation2, "_", gr$SYMBOL)

#promoter cCREs
promoter_name <- unique(gr$annotation3[which(gr$annotation2 == "Promoter")]) #18513 promoters
promoter_name <- intersect(promoter_name, gr[unique(coA_sf$queryHits)]$annotation3) #15601 promoters
promoter_name <- grep("Promoter_NA$", promoter_name, value = T, invert = T) #15600 promoters

#----- Defining pEPN, putative enhancer-promoter networks -----#
## Raw pEPNs ##
pEPN_raw <- mclapply(promoter_name, function(x){
  idx1 <- which(gr$annotation3 == x)
  out <- sort(coA_sf[which(coA_sf$queryHits %in% idx1),]$subjectHits)
  out <- unique(out)
  return(out)
}, mc.cores = 24)
names(pEPN_raw) <- promoter_name
saveRDS(pEPN_raw, "rds/06_pEPN_raw.rds")

pEPN_raw_anno <- mclapply(pEPN_raw, function(x){
  out <- gr$annotation2[x]
  return(out)
}, mc.cores = 24)
pEPN_raw_anno <- unlist(pEPN_raw_anno)
table(pEPN_raw_anno)
# 3'UTR            5'UTR DistalIntergenic       Downstream             Exon           Intron         Promoter 
# 6618             1453            53279             3584            12257            59403           124235 

## 1. Removing promoter cCREs ##
pEPN_NP <- mclapply(pEPN_raw, function(x){ #non promoter
  anno <- gr$annotation2[x]
  out <- x[which(anno != "Promoter")]
  return(out)
}, mc.cores = 24)
pEPN_NP <- pEPN_NP[which(lapply(pEPN_NP, length) %>% unlist != 0)] # 12676 promoters
saveRDS(pEPN_NP, "rds/06_pEPN_NP.rds")

## 2. Removing cCREs overlapping with epiAneufinder CN gain region ##
amp_gr_idx <- findOverlaps(gr, amp_gr) %>% queryHits() %>% unique %>% sort #50594
pEPN_CNAOL_DT <- mclapply(pEPN_NP, function(x){
  res1 <- which(x %ni% amp_gr_idx) %>% length
  res2 <- length(x) - res1
  out <- c("CNA" = res2, "nonCNA" = res1)
  return(out)
}, mc.cores = 24)
pEPN_CNAOL_DT <- do.call(rbind, pEPN_CNAOL_DT) %>% as.data.frame
colSums(pEPN_CNAOL_DT)
# CNA nonCNA 
# 34208 102386
write.csv(pEPN_CNAOL_DT, "output/Tables/06_pEPN_CNAoverlaps_Numbers.csv")

pEPN_NCNA <- mclapply(pEPN_NP, function(x) x[x %ni% amp_gr_idx], mc.cores = 24)
pEPN_NCNA <- pEPN_NCNA[which(lapply(pEPN_NCNA, length) %>% unlist != 0)] # 9977 promoters with at least one enhancer
saveRDS(pEPN_NCNA, "rds/06_pEPN_NCNA.rds")

#karyotype plot: CN gain region
pEPN_CNAOL <- mclapply(pEPN_NP, function(x) x[x %in% amp_gr_idx], mc.cores = 24)
pEPN_CNAOL <- pEPN_CNAOL[which(lapply(pEPN_CNAOL, length) %>% unlist != 0)] #3626 promoters
gr_CNAoverlaps <- gr[as.numeric(unlist(pEPN_CNAOL))] #34208
amp_gr_uq <- unique(amp_gr)
gr_pEPN_NP <- gr[unlist(pEPN_NP)]
pdf("output/Plots/06_Karyotype_CorrelatedCRE_CNA.pdf", width = 6, height = 6)
kp <- plotKaryotype(genome = "hg19", chromosomes = paste0("chr", c(1:22, "X")))
kpPlotRegions(kp, data=amp_gr_uq, col = "indianred1")
kpPlotRegions(kp, data=gr_pEPN_NP)
kp <- plotKaryotype(genome = "hg19", chromosomes = paste0("chr", c(1:22, "X")))
kpPlotRegions(kp, data=amp_gr_uq, col = "indianred1")
kpPlotRegions(kp, data=gr_CNAoverlaps)
dev.off()

#rank plot: # of overlaps CN gain region
TCGACNA_genes <- fread("ref/TCGA_nature2012.csv", header = F)[1,] %>% as.character()
df1 <- data.frame(gene = gsub("Promoter_", "", rownames(pEPN_CNAOL_DT)), CNA_overlaps = pEPN_CNAOL_DT$CNA, CNA_NonOverlaps = pEPN_CNAOL_DT$nonCNA)
df1$TCGA_AmpGene <- "No"
df1$TCGA_AmpGene[which(df1$gene %in% TCGACNA_genes)] <- "Yes"
df1 <- df1[order(df1$CNA_overlaps, decreasing = T),]
df1$rank <- c(1:nrow(df1))
df1 <- df1[order(df1$TCGA_AmpGene),]

p5 <- ggplot(df1, aes(x = rank, y = CNA_overlaps, color = TCGA_AmpGene)) + geom_point_rast(size = 1) + theme_ArchR() +
  labs(x = "Rank-sorted promoters", y = "# of correlated cCREs") + scale_color_manual(values = c("Yes" = "red", "No" = "gray"))
pdf("output/Plots/06_Rankplot_correlatedCREs_overlapping_epiAneufinderCNA.pdf", width = 5, height = 5)
p5
dev.off()
write.csv(df1, "output/Tables/06_pEPN_CNAoverlaps_Numbers.csv")

## 3. Removing promoters with extremely high number of correlated cCREs ##
ExtHigh_genes <- which(lapply(pEPN_NCNA, length) %>% unlist %>% sort(.,decreasing = T) > 50) %>% names() %>% gsub("Promoter_", "", .)
pEPN_final <- pEPN_NCNA[names(pEPN_NCNA) %ni% paste0("Promoter_", ExtHigh_genes)]
saveRDS(pEPN_final, "rds/06_pEPN_final.rds")

#visualize
tmp <- lapply(pEPN_NCNA, length) %>% unlist %>% sort(.,decreasing = T)
df2 <- data.frame(Promoter = names(tmp), No_cCREs = as.numeric(tmp))
df2$rank <- c(1:nrow(df2))
p6 <- ggplot(df2, aes(x = rank, y = No_cCREs)) + geom_point_rast(size = 1) + theme_ArchR() +
  labs(x = "Rank-sorted promoters", y = "# of correlated cCREs") + ggtitle("After removing promoter and overlap cCREs with CN gain regions")

tmp <- lapply(pEPN_final, length) %>% unlist %>% sort(.,decreasing = T)
df3 <- data.frame(Promoter = names(tmp), No_cCREs = as.numeric(tmp))
df3$rank <- c(1:nrow(df3))
p7 <- ggplot(df3, aes(x = rank, y = No_cCREs)) + geom_point_rast(size = 1) + theme_ArchR() +
  labs(x = "Rank-sorted promoters", y = "# of correlated cCREs") + ggtitle("After removing promoters with > 50 correlated cCRE")

df4 <- data.frame(No_cCREs = c(df2$No_cCREs[c(1:5000)], df3$No_cCREs[c(1:5000)]), 
                  rank = rep(c(1:5000),2),
                  Class = c(rep("Before",5000), rep("After",5000)))
p8 <- ggplot(df4, aes(x = rank, y = No_cCREs, color = Class)) + geom_point_rast(size = 1) + theme_ArchR() +
  labs(x = "Rank-sorted promoters", y = "# of correlated cCREs") +
  scale_color_manual(values = c("Before" = "darkgray", "After" = "black")) + ggtitle("Merging before and after filteration")

pdf("output/Plots/06_Rankplot_pEPN_final.pdf", width = 6, height = 6)
p6
p7
p8
dev.off()
write.csv(df2, "output/Tables/06_Rankplot_pEPN_final.csv")

#----- final pEPN properties -----#
# total inferred enhancers: 260829
# promoter 124235
# CNA 34208 
# Extremely high correlated cCREs 12735 (184 promoters)
# final 89651 pEnhancers-9793 promoters

df5 <- data.frame(Class = factor(c("Promoter","CNgain","ExtremeCorr","Final"), 
                                 levels = c("Promoter","CNgain","ExtremeCorr","Final")), 
                  No = c(124235,34208,12735,89651))
pie_colors <- ArchRPalettes$stallion %>% `names<-`(., df1$Class)
p9 <- ggplot(df5, aes(x = "x", y = No, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", color = "white") + 
  coord_polar(theta = "y") + theme_ArchR() + scale_fill_manual(values = pie_colors)
pdf("output/Plots/06_Pie_promoter_cCRE_class.pdf", width = 5, height = 5)
p9
dev.off()

df5$Perc <- df5$No / sum(df5$No)*100
# Class     No      Perc
# 1    Promoter 124235 47.630823
# 2      CNgain  34208 13.115106
# 3 ExtremeCorr  12735  4.882509
# 4       Final  89651 34.371561

#example, EPCAM: chr2:47,161,000-47,880,000
dt1 <- data.frame(seqnames = "chr2", 
                  start = start(geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == "EPCAM")]),
                  end = resize(gr[pEPN_final$Promoter_EPCAM], width = 1, fix = "center") %>% start)
idx1 <- which(dt1$start > dt1$end)
start1 <- dt1[idx1,]$end
end1 <- dt1[idx1,]$start
dt1[idx1,]$start <- start1
dt1[idx1,]$end <- end1

gr_loop1 <- GRanges(seqnames = dt1$seqnames, IRanges(start = dt1$start, end = dt1$end)) # 17 enhancers
p10 <- plotBrowserTrack(
  ArchRProj = arc,
  groupBy = "Clusters2",
  region = GRanges(seqnames = "chr2", ranges = IRanges(start = 47161001, end = 47880000)),
  features = gr[pEPN_final$Promoter_EPCAM],
  loops = gr_loop1)
grid::grid.newpage()
grid::grid.draw(p10)
pdf("output/Plots/06_GenomeTrack_EPCAM.pdf", width = 10, height = 8)
grid::grid.draw(p10)
dev.off()

#Distance from TSS
promoter_geneName <- gsub("Promoter_", "",names(pEPN_final))
promoter_gene_gr <- geneAnnoHg19$genes[geneAnnoHg19$genes$symbol %in% promoter_geneName]

distTSS <- mclapply(promoter_gene_gr$symbol, function(x){
  tss1 <- resize(promoter_gene_gr[which(promoter_gene_gr$symbol == x)], width = 1, fix = "start") %>% start %>% as.numeric()
  enhancer1 <- start(resize(gr[pEPN_final[[paste0("Promoter_", x)]]], width = 1, fix = "center")) %>% as.numeric()
  return(enhancer1 - tss1)
}, mc.cores = 24)
distTSS <- unlist(distTSS)
plot(sort(distTSS))
distTSSDF <- data.frame(Dist = distTSS)
summary(distTSS)

p11 <- ggplot(distTSSDF, aes(x = Dist)) + geom_histogram(bins = 100) + theme_ArchR() + xlim(-1000000,1000000) +
  labs(x = "Distance from TSS", y = "Frequency")
p12 <- ggplot(distTSSDF, aes(x = Dist)) + geom_line(stat="density", adjust = .25) + theme_ArchR() + xlim(-500000,500000) +
  labs(x = "Distance from TSS", y = "Density")

#No.of pEnhancers Mapped per promoters
freqEP <- lapply(pEPN_final, length) %>% unlist
freqEPDF <- data.frame(Number = freqEP)
summary(freqEPDF)
df6 <- table(freqEPDF) %>% as.data.frame()
p13 <- ggplot(df6, aes(x = freqEPDF, y = Freq)) + geom_bar(stat = "identity", width = 1) + theme_ArchR() + 
  labs(x = "# of pEnhancers per promoter", y = "Frequency")

#No.of promoters Mapped per pEnhancers
freqPE <- unlist(pEPN_final)
freqPE <- table(freqPE) %>% as.data.frame()
df7 <- table(freqPE$Freq) %>% as.data.frame()
p14 <- ggplot(df7, aes(x = Var1, y = Freq)) + geom_bar(stat = "identity", width = 1) + theme_ArchR() + 
  labs(x = "# of promoters per pEnhancer", y = "Frequency")

pdf("output/Plots/06_Proteries_pEPN_final.pdf", width = 6, height = 6)
p11
p12
p13
p14
dev.off()

#example, CCND1: chr11:68243000-70705000
dt2 <- data.frame(seqnames = "chr11", 
                  start = start(geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol == "CCND1")]),
                  end = resize(gr[pEPN_NCNA$Promoter_CCND1], width = 1, fix = "center") %>% start)
idx1 <- which(dt2$start > dt2$end)
start1 <- dt2[idx1,]$end
end1 <- dt2[idx1,]$start
dt2[idx1,]$start <- start1
dt2[idx1,]$end <- end1

gr_loop2 <- GRanges(seqnames = dt2$seqnames, IRanges(start = dt2$start, end = dt2$end)) # 153 linked cCREs
p15 <- plotBrowserTrack(
  ArchRProj = arc,
  groupBy = "Clusters2",
  region = GRanges(seqnames = "chr11", ranges = IRanges(start = 68243001, end = 70705000)),
  features = gr[pEPN_NCNA$Promoter_CCND1],
  loops = gr_loop2)
grid::grid.newpage()
grid::grid.draw(p15)
pdf("output/Plots/06_GenomeTrack_CCND1.pdf", width = 10, height = 8)
grid::grid.draw(p15)
dev.off()

gr_pEPN_final <- gr[unlist(pEPN_final)]
pdf("output/Plots/06_Karyotype_pEPN_final.pdf", width = 6, height = 6)
kp <- plotKaryotype(genome = "hg19", chromosomes = paste0("chr", c(1:22, "X")))
kpPlotRegions(kp, data=gr_pEPN_final)
dev.off()

#----- Integration of Cluster DA-cCREs and pEPN -----#
IntegDACRE_pEPN <- lapply(ClusterPeaksList, function(x){
  out <- mclapply(pEPN_final, function(y) findOverlaps(gr[y], x) %>% queryHits() %>% length, mc.cores = 24)
  out <- unlist(out)
  out <- sort(out, decreasing = T)
  return(out)
})
saveRDS(IntegDACRE_pEPN, "rds/06_IntegDACRE_pEPN.rds")

IntegDACRE_pEPN_nameSorted <- lapply(IntegDACRE_pEPN, function(x) x[sort(names(x))])
IntegDACRE_pEPN_nameSortedDF <- do.call(rbind, IntegDACRE_pEPN_nameSorted)
df_out <- t(IntegDACRE_pEPN_nameSortedDF)
rownames(df_out) <- gsub("Promoter_", "", rownames(df_out))
write.csv(df_out, "output/Tables/06_IntegDACRE_pEPN_nameSortedDF.csv")
idx2 <- apply(IntegDACRE_pEPN_nameSortedDF, 1, max)
df8 <- IntegDACRE_pEPN_nameSortedDF / idx2

col_fun1 <- colorRamp2(c(0,0.1,0.2,0.3,0.4), viridis(5,option = "D"))
fh <- function(x) hclust(dist(x), method = "ward.D2")

ht1 <- Heatmap(df8, name = "Prop.overlaps",col = col_fun1, cluster_rows = fh, cluster_columns = fh,
               show_column_names = F, row_names_gp = gpar(fontsize = 6))
p16 <- draw(ht1)
pdf("output/Plots/06_Heatmap_Overlaps_DA-cCRE_pEnhancers_Prop_v1.pdf", width = 10, height = 5)
p16
dev.off()

column_dend1 <- column_dend(p16)
idx3 <- list(column_dend1[[1]][[1]], # pan-cancer
             column_dend1[[1]][[2]], # ER+
             column_dend1[[2]][[1]], # T cell
             column_dend1[[2]][[2]][[1]][[1]][[1]], # B cell
             column_dend1[[2]][[2]][[1]][[1]][[2]], # Macrophage
             column_dend1[[2]][[2]][[1]][[2]][[2]], # Fib/Endo
             column_dend1[[2]][[2]][[2]][[2]]) # Fib/Endo
promoter_DAcCRE_genelist <- list(
  "1" = colnames(df8)[unlist(idx3[[1]])] %>% gsub("Promoter_", "", .) %>% sort(),
  "2" = colnames(df8)[unlist(idx3[[2]])] %>% gsub("Promoter_", "", .) %>% sort(),
  "3" = colnames(df8)[unlist(idx3[[3]])] %>% gsub("Promoter_", "", .) %>% sort(),
  "4" = colnames(df8)[unlist(idx3[[4]])] %>% gsub("Promoter_", "", .) %>% sort(),
  "5" = colnames(df8)[unlist(idx3[[5]])] %>% gsub("Promoter_", "", .) %>% sort(),
  "6" = colnames(df8)[unlist(idx3[[6]])] %>% gsub("Promoter_", "", .) %>% sort(),
  "7" = colnames(df8)[unlist(idx3[[7]])] %>% gsub("Promoter_", "", .) %>% sort() 
)
lapply(names(promoter_DAcCRE_genelist), function(x){
  df <- data.frame(gene = promoter_DAcCRE_genelist[[x]])
  write.table(df, paste0("output/Tables/06_Heatmap_Overlaps_DA-cCRE_pEnhancers_Prop_v1_gene_", x, ".txt"), row.names = F, col.names = F, quote = F)
  return(NULL)
})

# Rankplot all clusters
p17 <- lapply(names(IntegDACRE_pEPN), function(x){
  df_tmp <- data.frame(gene = gsub("Promoter_", "", names(IntegDACRE_pEPN[[x]])), pEnhancer = as.numeric(IntegDACRE_pEPN[[x]]))
  df_tmp$rank <- c(1:nrow(df_tmp))
  out <- ggplot(df_tmp, aes(x=rank,y=pEnhancer)) + geom_point_rast(size = 1) + theme_ArchR() + 
    ggtitle(x) + labs(x = "Rank-sorted promoters", y = "# of linked pEnhancers")
  return(out)
})
pdf("output/Plots/06_Rankplot_Overlaps_DAcCRE_pEnhancers.pdf", width = 5, height = 5)
p17
dev.off()


#----- epiAneufinder Amplified regions output -----#
amp_gr
options(scipen=100)
d3 <- data.frame(seqnames = seqnames(amp_gr), start = start(amp_gr)-1, end = end(amp_gr))
write.table(d3, "output/output_bed/06_epiAneufinder_CNgainRegions.bed", row.names = F, col.names = F, quote = F, sep = "\t")

#----- for ROSE analysis -----#
fragments <- getFragmentsFromProject(arc)
fragments <- unlist(fragments)
clusters <- unique(arc$Clusters2)

lapply(clusters, function(x){
  idx_cell <- arc$cellNames[which(arc$Clusters2 == x)]
  idx_frag <- which(fragments$RG %in% idx_cell)
  gr1 <- fragments[idx_frag]
  inserts <- c(
    GRanges(seqnames = seqnames(gr1), ranges = IRanges(start(gr1), start(gr1))),
    GRanges(seqnames = seqnames(gr1), ranges = IRanges(end(gr1), end(gr1)))
  )

  #save bed file (2bp)
  dt <- data.frame(seqnames = seqnames(inserts), start = start(inserts)-1, end = end(inserts))
  dt <- dt[order(dt$seqnames, dt$start, dt$end), ]
  dt <- dt[dt$end > dt$start,]
  dt$ncol <- paste0("insert_", c(1:length(inserts)))
  write.table(dt, paste0("output/inserts_bam/", x, "_inserts.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
})

options(scipen = 10000000)
peaknames <- c(paste0("C", c(5:13,15:16)), paste0("Ep", c(1:18,20:35)))
gr_perCluster <- lapply(peaknames, function(x) readRDS(paste0("output/PeakCalls/", x, "-reproduciblePeaks.gr.rds")))
names(gr_perCluster) <- peaknames
mclapply(peaknames, function(x){
  gr1 <- gr_perCluster[[x]]
  dt <- data.frame(seqnames = seqnames(gr1), start = start(gr1)-1, end = end(gr1))
  dt$ncol <- paste0("peak_", c(1:length(gr1)))
  write.table(dt, paste0("output/inserts_bam/", x, "_peaks.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
})

#----- Comparison between SE and pEPN -----#
gr_SE <- lapply(peaknames, function(x){
  print(x)
  out <- import.bed(paste0("/mnt/host_mnt/Users/nextganken/tools/ROSE/results_20241004/",x,"/", x, "_peaks_Gateway_SuperEnhancers.bed"))
  return(out)
}) 
names(gr_SE) <- peaknames

#check overlaps with epiAneufinder-defined CN gain regions
gr_SE_AmpOverlap <- lapply(names(ClusterSampleCorrespond), function(x){
  gr <- gr_SE[[x]]
  n <- ClusterSampleCorrespond[x]
  out <- subsetByOverlaps(gr, amp_gr[amp_gr$Sample == n])
  return(out)
})
names(gr_SE_AmpOverlap) <- names(ClusterSampleCorrespond)
gr_SE_AmpOverlap <- gr_SE_AmpOverlap[unlist(lapply(gr_SE_AmpOverlap, length)) !=0]
gr_SE_AmpOverlap

gr_SE <- lapply(gr_SE, function(x){
  x$Amp <- "N"
  return(x)
})

for(i in names(gr_SE_AmpOverlap)){
  fo1 <- findOverlaps(gr_SE[[i]], gr_SE_AmpOverlap[[i]])
  gr_SE[[i]]$Amp[queryHits(fo1)] <- "Y"
}

NumSETb <- lapply(names(gr_SE), function(x){
  out <- table(gr_SE[[x]]$Amp) %>% as.data.frame()
  out$Var1 <- as.character(out$Var1)
  out$Cluster <- x
  if(length(which(out$Var1 == "Y")) == 0){
    out <- rbind(out, c("Y", 0, x))
  }
  return(out)
}) %>% do.call(rbind, .)
NumSETb$Freq <- as.numeric(NumSETb$Freq)
NumSETb$Var1 <- factor(NumSETb$Var1, levels = c("Y", "N"))
NumSETb$Cluster <- factor(NumSETb$Cluster, levels = peaknames)

p18 <- ggplot(NumSETb, aes(x = Cluster, y = Freq, fill = Var1)) +　geom_bar(stat = "identity") + theme_ArchR() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("Y" = "red", "N" = "darkgray")) +
  labs(y = "# of Super-Enhancers called by ROSE using ATAC signals per cluster")
pdf("output/Plots/06_Barplot_NumClusterSuperEnhancers.pdf", width = 8, height = 4)
p18
dev.off()
write.csv(NumSETb, "output/Tables/06_Barplot_NumClusterSuperEnhancers.csv")

#Promoter with ≥5 pEnhancers and SEs
idx4 <- lapply(IntegDACRE_pEPN, function(x) which(x >= 5) %>% length) %>% unlist
idx4 <- names(which(idx4 > 100))
gene_gr <- geneAnnoHg19$genes
strand(gene_gr) <- "*"

SE_genes <- lapply(peaknames, function(x){
  gr_SE1 <- extendGR(gr_SE[[x]], upstream = 500000, downstream = 500000)
  fo1 <- findOverlaps(gene_gr, gr_SE1)
  out <- gene_gr$symbol[queryHits(fo1)] %>% as.character() %>% sort %>% unique
  return(out)
})
names(SE_genes) <- peaknames

SE_genes_Jaccard <- mclapply(peaknames, function(x){
  out <- lapply(peaknames, function(y){
    intersect_num <- intersect(SE_genes[[x]], SE_genes[[y]]) %>% length()
    union_num <- length(SE_genes[[x]]) + length(SE_genes[[y]]) - intersect_num
    return(intersect_num/union_num)
  })
  out <- unlist(out)
  return(out)
}, mc.cores = 24)
SE_genes_Jaccard <- do.call(rbind, SE_genes_Jaccard)

geneNo <- lapply(SE_genes, length) %>% unlist %>% as.numeric()
colnames(SE_genes_Jaccard) <- peaknames
rownames(SE_genes_Jaccard) <- paste0(peaknames, " (", geneNo, ")")

fh = function(x) hclust(dist(x), method = "ward.D2")
col_fun1 <- colorRamp2(c(0.2,0.6,1), c("white", "orange", "red"))
ht2 <- Heatmap(SE_genes_Jaccard, name = "Jaccard similarity", 
               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), 
               cluster_rows = fh, cluster_columns = fh, col = col_fun1, border = "black")
p19 <- draw(ht2)

pdf("output/Plots/06_Heatmap_Jaccard_geneNearbySE.pdf", width = 9, height = 8)
p19
dev.off()

df9 <- data.frame(Cluster = names(SE_genes), geneNo = as.numeric(geneNo), Class = "SE")
df10 <- lapply(IntegDACRE_pEPN, function(x) length(which(x >= 5))) %>% unlist
df10 <- rbind(data.frame(Cluster = names(df10), geneNo = as.numeric(df10), Class = "DAcCRE_pEnhancer"), c("Ep13", 0, "DAcCRE_pEnhancer"))
df11 <- rbind(df9, df10)
df11$Cluster <- factor(df11$Cluster, levels = peaknames)
df11$geneNo <- as.numeric(df11$geneNo)

p20 <- ggplot(df11, aes(x = Cluster, y = geneNo, fill = Class)) + geom_bar(stat = "identity", position = "dodge") + theme_ArchR() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("DAcCRE_pEnhancer" = "#D51F26", "SE" = "#272E6A"))
pdf("output/Plots/06_Barplot_NumberGeneNearbySE_OverlapsDACREpEPN.pdf", width = 8, height =4)
p20
dev.off()

#----- DA-cCRE overlapped with pEPN -----#
ClusterPeaksList <- readRDS("rds/06_ClusterPeaksList.rds")
pEPN_final <- readRDS("rds/06_pEPN_final.rds")
IntegDACRE_pEPN <- readRDS("rds/06_IntegDACRE_pEPN.rds")

IntegDACRE_pEPN_gr <- lapply(ClusterPeaksList, function(x){
  idx <- mclapply(pEPN_final, function(y) findOverlaps(gr[y], x) %>% subjectHits(), mc.cores = 24)
  idx <- unlist(idx) %>% unique %>% sort %>% as.numeric
  out <- x[idx]
  return(out)
})
saveRDS(IntegDACRE_pEPN_gr, "rds/06_IntegDACRE_pEPN_gr.rds")

df12 <- data.frame(NonOverlap = lapply(ClusterPeaksList, length) %>% unlist - lapply(IntegDACRE_pEPN_gr, length) %>% unlist,
                   Overlap = lapply(IntegDACRE_pEPN_gr, length) %>% unlist)
df12$Cluster <- factor(rownames(df12), levels = rownames(df12))
df12 <- reshape2::melt(df12)
df12$variable <- factor(df12$variable, levels = c("Overlap", "NonOverlap"))

p21 <- ggplot(df12, aes(x = Cluster, y = value, fill = variable)) +　geom_bar(stat = "identity") + theme_ArchR() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("NonOverlap" = "darkgray", "Overlap" = "orange")) +
  labs(y = "# of DA-cCREs")
pdf("output/Plots/06_Barplot_NumClusterSpecificCRE_Ovlp_pEnhancer.pdf", width = 8, height = 4)
p21
dev.off()

lapply(names(IntegDACRE_pEPN_gr), function(x){
  gr_out <- IntegDACRE_pEPN_gr[[x]]
  d <- data.frame(seqnames = seqnames(gr_out), start = start(gr_out)-1, end = end(gr_out))
  write.table(d, paste0("output/output_bed/ovlp_pEnhancer_DAcCRE/06_Overlap_DAcCRE_pEnhancer_",x,".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
})
