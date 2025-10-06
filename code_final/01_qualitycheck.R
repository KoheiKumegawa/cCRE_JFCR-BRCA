#----------------------------------------------------------------------------
# 01_qualitycheck.R
#----------------------------------------------------------------------------
library(ArchR)
library(ggsignif)
library(viridis)
library(ComplexHeatmap)
library(circlize)
addArchRThreads(threads = 24)
addArchRGenome("hg19")

arc <- readRDS("rds/00_arc.rds")
SampleTypeColors <- RColorBrewer::brewer.pal(5, "Dark2")

#----- Clustering (No batch correction) -----#
arc <- addIterativeLSI(arc, useMatrix = "TileMatrix", name = "IterativeLSI", force = T)
arc <- addUMAP(arc, reducedDims = "IterativeLSI", force = T)
arc <- addClusters(arc, reducedDims = "IterativeLSI", maxClusters = 40, force = T)

p1 <- plotEmbedding(arc, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", plotAs = "points", size = 0.5)
p2 <- plotEmbedding(arc, colorBy = "cellColData", name = "Sample", embedding = "UMAP", plotAs = "points", size = 0.5, labelMeans = F)
p3 <- plotEmbedding(arc, colorBy = "cellColData", name = "SampleType", embedding = "UMAP", plotAs = "points", pal = SampleTypeColors, size = 0.5, labelMeans = F)
p4 <- plotEmbedding(arc, colorBy = "cellColData", name = "Method", embedding = "UMAP", plotAs = "points", size = 0.5, pal = "", labelMeans = F)
plotPDF(p1,p2,p3,p4, name = "01_UMAP_NoCorrect.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#----- Clustering (Harmony) -----#
arc <- addHarmony(
  ArchRProj = arc,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Method"
)
arc <- addUMAP(arc, reducedDims = "Harmony", name = "UMAP_Harmony", force = T)
arc <- addClusters(arc, reducedDims = "Harmony", name = "Clusters_Harmony", maxClusters = 40, force = T)

p5 <- plotEmbedding(arc, colorBy = "cellColData", name = "Clusters_Harmony", embedding = "UMAP_Harmony", plotAs = "points", size = 0.5)
p6 <- plotEmbedding(arc, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Harmony", plotAs = "points", size = 0.5, labelMeans = F)
p7 <- plotEmbedding(arc, colorBy = "cellColData", name = "SampleType", embedding = "UMAP_Harmony", plotAs = "points", pal = SampleTypeColors, size = 0.5, labelMeans = F)
p8 <- plotEmbedding(arc, colorBy = "cellColData", name = "Method", embedding = "UMAP_Harmony", plotAs = "points", size = 0.5, labelMeans = F)
plotPDF(p5,p6,p7,p8, name = "01_UMAP_Harmony.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)


saveRDS(arc, "rds/01_arc.rds")

#----- No Batch Correction VS Harmony-----#
df1 <- table(arc$Clusters, arc$Method) %>% as.data.frame()
df1$Var1 <- factor(df1$Var1, paste0("C", c(1:33)))
p9 <- ggplot(df1, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat = "identity", position = "fill") + theme_ArchR() +
  scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("DSC" = "#D51F26", "DSCI" = "#272E6A")) + 
  theme(axis.text = element_text(size= 6)) + labs(x = "Clusters", y = "Proportion") + ggtitle("No batch effect control")

df2 <- table(arc$Clusters_Harmony, arc$Method) %>% as.data.frame()
df2$Var1 <- factor(df2$Var1, paste0("C", c(1:35)))
p10 <- ggplot(df2, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat = "identity", position = "fill") + theme_ArchR() +
  scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c("DSC" = "#D51F26", "DSCI" = "#272E6A")) + 
  theme(axis.text = element_text(size= 6)) + labs(x = "Clusters", y = "Proportion") + ggtitle("Harmony batch effect control")
pdf("output/Plots/01_Barplot_CompareHarmony.pdf",width = 8, height = 4)
p9
p10
dev.off()

CS_NBC <- table(arc$Clusters, arc$Method) %>% as.data.frame.matrix()
CS_NBC <- CS_NBC / rowSums(CS_NBC)
CS_NBC <- apply(CS_NBC, 1, max)
CS_NBC <- 1 - CS_NBC

CS_HBC <- table(arc$Clusters_Harmony, arc$Method) %>% as.data.frame.matrix()
CS_HBC <- CS_HBC / rowSums(CS_HBC)
CS_HBC <- apply(CS_HBC, 1, max)
CS_HBC <- 1 - CS_HBC

df3 <- data.frame(CS = c(CS_NBC, CS_HBC), BC = factor(c(rep("NoCorrection", length(CS_NBC)), rep("Harmony", length(CS_HBC))), levels = c("NoCorrection", "Harmony")))
p11 <- ggplot(df3, aes(x = BC, y = CS)) + geom_jitter(height = 0, width = 0.1, shape = 3) + geom_boxplot(outlier.shape = NA, alpha = 0) +
  theme_ArchR() + labs(x="", y = "Complexity Score [1-max(Prop.)]") + geom_signif(comparisons = list(c("NoCorrection", "Harmony")), test = "wilcox.test")
pdf("output/Plots/01_Boxplot_ClusterComplexityScore.pdf",width = 3, height = 5)
p11
dev.off()

#----- QC metrics -----#
p12 <- plotGroups(arc, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment",plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p13 <- plotGroups(arc, groupBy = "Clusters", colorBy = "cellColData", name = "TSSEnrichment",plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p14 <- plotGroups(arc, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)",plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p15 <- plotGroups(arc, groupBy = "Clusters", colorBy = "cellColData", name = "log10(nFrags)",plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p16 <- plotGroups(arc, groupBy = "Sample", colorBy = "cellColData", name = "PromoterRatio",plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p17 <- plotGroups(arc, groupBy = "Clusters", colorBy = "cellColData", name = "PromoterRatio",plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

p18 <- plotEmbedding(arc, colorBy = "cellColData", name = "TSSEnrichment", pal = viridis(256), embedding = "UMAP", plotAs = "points", size = 0.5)
p19 <- plotEmbedding(arc, colorBy = "cellColData", name = "log10(nFrags)", pal = viridis(256), embedding = "UMAP", plotAs = "points", size = 0.5)
p20 <- plotEmbedding(arc, colorBy = "cellColData", name = "PromoterRatio", pal = viridis(256), embedding = "UMAP", plotAs = "points", size = 0.5)

p21 <- plotFragmentSizes(arc)
p22 <- plotTSSEnrichment(arc)
 
pdf("output/Plots/01_Boxplot_QCMetrics.pdf",width = 8, height = 4)
p12
p13
p14
p15
p16
p17
dev.off()
pdf("output/Plots/01_UMAP_QCMetrics.pdf",width = 8, height = 4)
p18
p19
p20
dev.off()
plotPDF(p18,p19,p20, name = "01_UMAP_QCMetrics.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)
plotPDF(p21,p22, name = "01_Density_QCMetrics.pdf", ArchRProj = arc, addDOC = FALSE, width = 8, height = 8)

#----- Nucleosome-free fragments output -----#
fragments <- getFragmentsFromProject(arc)
saveRDS(fragments, "rds/01_fragments.rds")
# fragments <- readRDS("rds/01_fragments.rds")

for(i in names(fragments)){
  outFile <- paste0("output/output_bed/nucfreefrags/", i, ".nucfreefrags.tsv.gz")
  print(outFile)
  frag <- fragments[[i]]
  idx <- which(width(frag) < 150)
  frag <- frag[idx]
  dt <- data.frame(seqnames = seqnames(frag), start = start(frag), end = end(frag), RG = stringr::str_split(frag$RG, "#", simplify = T)[,2])
  
  rm(frag)
  gc(reset = T)
  gc(reset = T)
  
  dt <- unique(dt[order(dt$seqnames, dt$start, dt$end), ])
  dt <- dt[dt$end > dt$start,] #this can sometimes be a weird aligner error!
  data.table::fwrite(dt, stringr::str_split(outFile, pattern="\\.gz", simplify=TRUE)[,1], sep = "\t", col.names = FALSE)
  Rsamtools::bgzip(stringr::str_split(outFile, pattern="\\.gz", simplify=TRUE)[,1], outFile)
  
  rm(dt)
  gc(reset = T)
  gc(reset = T)
}

#----- Nucleosome-free fragments clustering -----#
setwd("anal_nucfree/")

#assign sample tsv files
pt_ls <- read.csv("../ref/patient_list_v2.csv")
sampleName <- paste0(pt_ls$ID, ".nucfreefrags.tsv.gz")
names(sampleName) <- pt_ls$ID
outFile <- as.character(sampleName)

#make ArrowFiles
ArrowFiles <- character(length(sampleName))
ArrowFiles <- createArrowFiles(inputFiles = paste0("../output/output_bed/nucfreefrags/", outFile),
                               sampleNames = names(sampleName),
                               minTSS = 0, 
                               minFrags = 1, 
                               addTileMat = TRUE, addGeneScoreMat = F, 
                               force = TRUE)

#pre-filtered ArchR project
arc_nucfree <- ArchRProject(ArrowFiles, outputDirectory = "output", copyArrows = F) #34403 cells

arc <- readRDS("../rds/01_arc.rds")
arc_nucfree <- arc_nucfree[arc$cellNames] # the same cells 22773 cells

arc_nucfree <- addIterativeLSI(arc_nucfree, useMatrix = "TileMatrix", name = "IterativeLSI", force = T)
arc_nucfree <- addUMAP(arc_nucfree, reducedDims = "IterativeLSI", force = T)
arc_nucfree <- addClusters(arc_nucfree, reducedDims = "IterativeLSI", maxClusters = 40, force = T)

p21 <- plotEmbedding(arc_nucfree, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", plotAs = "points", size = 0.5)
p22 <- plotEmbedding(arc_nucfree, colorBy = "cellColData", name = "Sample", embedding = "UMAP", plotAs = "points", size = 0.5, labelMeans = F)
pdf("../output/Plots/01_UMAP_nucleosomefreefragments.pdf", width = 8, height = 8)
p21
p22
dev.off()

#cluster-overlap analysis
cM <- data.frame(All = arc$Clusters, NucFree = arc_nucfree[arc$cellNames]$Clusters)
cM <- table(cM)[paste0("C",c(1:33)), paste0("C", c(1:31))] %>% as.data.frame.matrix()

mtx <- cM / rowSums(cM) * 100
apply(mtx, 1, function(i) order(i,decreasing = T)[1])
idx1 <- c(16,15,1,2,4,3,5,7,6,11,10,13,8,9,29,17,20,23,26,12,31,30,18,19,14,28,27,22,21,25,24)
sort(idx1)

col_fun1 <- colorRamp2(c(0,50,100), c("white", "yellow", "red"))
ht1 <- Heatmap(mtx[,idx1], name = "%Overlap (Nucleosome-free clust in All frag. clust)", cluster_rows = F, cluster_columns = F, col = col_fun1, 
               row_title = "All fragments", column_title = "Nucleosome-free fragments only", border_gp = gpar(col = "black", lty = 1))
p23 <- draw(ht1)

pdf("../output/Plots/01_Heatmap_Overlap_Allfrag_vs_NucFreefrag_Clusters.pdf", width = 9, height = 6)
p23
dev.off()

saveRDS(arc_nucfree, "../rds/01_arc_nucfree.rds")
