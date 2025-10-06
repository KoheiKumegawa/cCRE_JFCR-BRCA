#------------------------------------------------------------------------------
# 05_PeakCall.R
#------------------------------------------------------------------------------
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggrastr)
library(ChIPseeker)
library(ComplexHeatmap)
library(circlize)
library(viridis)
addArchRThreads(threads = 24)
addArchRGenome("hg19")
source("code/bulkATAC_functions.R")

arc <- readRDS("rds/02_arc.rds")
arc_epi <- readRDS("rds/04_arc_epi.rds")
pt_ls <- read.csv("ref/patient_list_v2.csv")

#----- Pseudo-bulk peak calling -----#
fragments <- parallel::mclapply(pt_ls$ID, function(x){
  fr <- fread(paste0("data/", x, "_scatac.fragments.tsv.gz")) %>% data.frame(.)
  out <- GRanges(seqnames = fr[,1], IRanges(start = fr[,2]+1, end = fr[,3]))
  return(out)
}, mc.cores = 10)
names(fragments) <- pt_ls$ID

#fragments to inserts
parallel::mclapply(names(fragments), function(x){
  FragmentGRToInsert(gr = fragments[[x]], 
                     name = x, 
                     outputDirBED = "output/inserts_bed/")
}, mc.cores = 10)

#Peak call using inserts
insertFiles <- list.files("output/inserts_bed/") %>% `names<-`(., gsub("_inserts.bed", "", .))
parallel::mclapply(names(insertFiles), function(x){
  RunMacs2(inputBedPATH = paste0("output/inserts_bed/", insertFiles[x]),
           inputBedName = x,
           genome = "hg19",
           outputDir = "output/sample_peaks/",
           shift = -75,
           extsize = 150,
           method = "q",
           cutoff = 0.01)
}, mc.cores = 10)

#Make sample peak set
BSgenome   <- BSgenome.Hsapiens.UCSC.hg19
chromSizes <- GRanges(names(seqlengths(BSgenome)), IRanges(1, seqlengths(BSgenome))) %>% 
  GenomeInfoDb::keepStandardChromosomes(., pruning.mode = "coarse")
blacklist  <- rtracklayer::import.bed("ref/hg19-blacklist.v2.bed")
gr_ls <- GenomicRanges::GRangesList(lapply(list.files("output/sample_peaks_q001///", pattern = "summits.bed", full.names = T), function(x) import.bed(x)))
names(gr_ls) <- gsub("_summits.bed", "", list.files("output/sample_peaks_q001//", pattern = "summits.bed"))

gr_ls_proc <- parallel::mclapply(gr_ls, function(x){
  gr <- resize(x, width = 501, fix = "center") %>%
    subsetByOverlaps(., chromSizes, type = "within") %>%
    subsetByOverlaps(., blacklist, invert=TRUE) %>%
    MakeSamplePeakSet(., by = "score")
  mcols(gr)$scorePerMillion <- mcols(gr)$score / (sum(mcols(gr)$score) / 1000000)
  return(gr)
}, mc.cores = 10) %>% GenomicRanges::GRangesList(.)

#cancer type specific peakset
gr_cumulative <- MakeSamplePeakSet(unlist(gr_ls_proc), by = "scorePerMillion")
mcols(gr_cumulative)$sampleOverlap <- countOverlaps(gr_cumulative, gr_ls_proc)
reproduciblePeaks <- gr_cumulative[which(mcols(gr_cumulative)$sampleOverlap >= 2 & seqnames(gr_cumulative) %ni% "chrY")]

#----- ArchR Cluster-based peak calling -----#
# adding epithelial cluster info
newClusters <- arc$Clusters
names(newClusters) <- arc$cellNames
newClusters[arc_epi$cellNames] <- arc_epi$Clusters

arc$Clusters2 <- newClusters[arc$cellNames]
p1 <- plotEmbedding(arc, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP", plotAs = "points", size = 0.5)
plotPDF(p1, name = "05_UMAP_newClusters.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#peak call 
pathToMacs2 <- findMacs2()
arc <- addGroupCoverages(arc, groupBy = "Clusters2", force = T)
arc <- addReproduciblePeakSet(arc, groupBy = "Clusters2", pathToMacs2 = pathToMacs2, method = "q", cutOff = 0.1, force = T)
arc <- addPeakMatrix(arc, force = T)

pmat01 <- getMatrixFromProject(arc, useMatrix = "PeakMatrix", binarize = F)
pmat01_binarize <- getMatrixFromProject(arc, useMatrix = "PeakMatrix", binarize = T)

df1 <- data.frame(NFrag = log10(rowSums(assay(pmat01))), NCell = log10(rowSums(assay(pmat01_binarize))))
df1$Eval <- "Valid"
df1$Eval[which(df1$NFrag < 1 & df1$NCell < 1)] <- "Sperse"
table(df1$Eval)
# Sperse  Valid 
# 702 224585 

p2 <- ggplot(df1, aes(x = NFrag, y = NCell, color = Eval)) + geom_point_rast(shape = 1) + theme_ArchR() + 
  scale_color_manual(values = c("Valid" = "gray", "Sperse" = "blue")) +
  labs(x = "# of inserts per each peak [log10]", y = "# of cells with inserts in each peak [log10]")
pdf("output/Plots/05_Scatterplot_peakSperseness.pdf", width = 5, height = 5)
p2
dev.off()

gr_valid <- rowRanges(pmat01)[which(df1$Eval == "Valid")]

saveRDS(arc, "rds/05_arc.rds")
saveRDS(gr_valid, "rds/05_gr_valid.rds")

#----- cCRE genomic positions -----#
gr <- arc@peakSet
gr <- subsetByOverlaps(gr, gr_valid)
names(gr) <- c(paste0("cCRE_", c(1:length(gr))))

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak(gr, tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
df <- peakAnno@annoStat

pie_colors <- ArchRPalettes$stallion %>% `names<-`(., df$Feature)
p3 <- ggplot(df, aes(x = "x", y = Frequency, fill = Feature)) +
  geom_bar(stat = "identity", position = "stack", color = "white") + 
  coord_polar(theta = "y") + theme_ArchR() + scale_fill_manual(values = pie_colors)
pdf("output/Plots/05_cCRE_GenomicPosition.pdf", width = 5, height = 5)
p3
plotAnnoPie(peakAnno)
dev.off()

#----- Overlaps between cluster peaks -----#
peaknames <- list.files("output/PeakCalls/", pattern = "gr.rds") %>% gsub("-reproduciblePeaks.gr.rds", "", .)
peaknames <- c(paste0("C", c(5:13,15:16)), paste0("Ep", c(1:18,20:35)))

gr_perCluster <- lapply(peaknames, function(x) readRDS(paste0("output/PeakCalls/", x, "-reproduciblePeaks.gr.rds")))
names(gr_perCluster) <- peaknames
peakNo <- lapply(gr_perCluster, length) %>% unlist

perClusterOverlapsJaccard <- mclapply(peaknames, function(x){
  out <- lapply(peaknames, function(y){
    intersect_num <- findOverlaps(gr_perCluster[[x]], gr_perCluster[[y]]) %>% length()
    union_num <- length(gr_perCluster[[x]]) + length(gr_perCluster[[y]]) - intersect_num
    return(intersect_num/union_num)
  })
  out <- unlist(out)
  return(out)
}, mc.cores = 24)
perClusterOverlapsJaccard <- do.call(rbind, perClusterOverlapsJaccard)
colnames(perClusterOverlapsJaccard) <- peaknames
rownames(perClusterOverlapsJaccard) <- paste0(peaknames, " (", peakNo, ")")
perClusterOverlapsJaccard

fh = function(x) hclust(dist(x), method = "ward.D2")
col_fun1 <- colorRamp2(c(0.2,0.6,1), c("white", "orange", "red"))
ht1 <- Heatmap(perClusterOverlapsJaccard, name = "Jaccard similarity", 
               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), 
               cluster_rows = fh, cluster_columns = fh, col = col_fun1, border = "black")
p4 <- draw(ht1)
pdf("output/Plots/05_Heatmap_Jaccard_PeakperClusters.pdf", width = 9, height = 8)
p4
dev.off()

#----- Cluster peak annotation -----#
ClusterPeakAnno <- mclapply(gr_perCluster, function(x) annotatePeak(x, tssRegion=c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db"), mc.cores = 12)

Anno_df <- lapply(names(gr_perCluster), function(i){
  tmp <- ClusterPeakAnno[[i]]@annoStat$Frequency * length(ClusterPeakAnno[[i]]@anno) * 0.01
  out <- data.frame(Feature = ClusterPeakAnno[[i]]@annoStat$Feature, PeakNo = tmp, Sample = i)
  return(out)
})
Anno_df <- do.call(rbind, Anno_df)
p5 <- ggplot(Anno_df, aes(x = factor(Sample, levels = peaknames), y = PeakNo, fill = Feature)) + geom_bar(stat = "identity") +
  theme_ArchR() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = pie_colors) +
  labs(x = "Cluster", y = "# of cCREs")
pdf("output/Plots/05_Barplot_PeakperClusters_Anno.pdf", width = 7, height = 4)
p5
dev.off()

chr_df <- lapply(gr_perCluster, function(x) table(seqnames(x))) %>% do.call(cbind, .)
chr_df <- chr_df[,peaknames]
chr_df <- reshape2::melt(chr_df)
chr_colors <- c(ArchRPalettes$stallion2, ArchRPalettes$calm)[c(1:23)] %>% `names<-`(., paste0("chr", c(1:22,"X")))

p6 <- ggplot(chr_df, aes(x = factor(Var2, levels = peaknames), y = value, fill = Var1)) + geom_bar(stat = "identity") +
  theme_ArchR() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = chr_colors)
pdf("output/Plots/05_Barplot_PeakperClusters_Chr.pdf",  width = 7, height = 4)
p6
dev.off()

#----- Peak overlaps -----#
TCGA_ATAC_peaks_hg19 <- import.bed("ref/06_TCGApeaks_hg19.bed")
complist_gr <- list(scATAC = gr, pseudoBulk = reproduciblePeaks, TCGA = TCGA_ATAC_peaks_hg19)

ChIPseeker::vennplot(complist_gr)

pdf("output/Plots/05_Vennplot_sc_pseudoBulk_TCGA.pdf", width = 5, height = 5)
ChIPseeker::vennplot(complist_gr)
dev.off()

lapply(names(complist_gr), function(x){
  gr1 <- complist_gr[[x]]
  dt <- data.frame(seqnames = seqnames(gr1), start = start(gr1)-1, end = end(gr1))
  dt$ncol <- paste0("peak_", c(1:length(gr1)))
  write.table(dt, paste0("output/output_bed/05_PeaksForCompare_", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
})
