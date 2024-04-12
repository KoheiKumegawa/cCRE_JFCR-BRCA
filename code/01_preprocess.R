#----------------------------------------------------------------------------
# 01_preprocess.R
#----------------------------------------------------------------------------
library(ArchR)
addArchRThreads(threads = 22)
addArchRGenome("hg19")

#----- make ArchRProj and remove doublets -----#
#assign sample tsv files
pt_ls <- read.csv("ref/patient_list_20231201.csv")
sampleName <- paste0(pt_ls$prefix, ".fragments.tsv.gz")
names(sampleName) <- pt_ls$ID
outFile <- as.character(sampleName)

#make ArrowFiles
ArrowFiles <- character(length(sampleName))
ArrowFiles <- createArrowFiles(inputFiles = paste0("data/", outFile),
                               sampleNames = names(sampleName),
                               minTSS = 4, 
                               minFrags = 1000, 
                               addTileMat = TRUE, addGeneScoreMat = TRUE, 
                               force = TRUE)

#doublet
doubScores <- addDoubletScores(ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)

#pre-filtered ArchR project
pre_arc <- ArchRProject(ArrowFiles, outputDirectory = "output", copyArrows = T) #total:34745 cells
pre_arc <- filterDoublets(pre_arc)

#add sample info
pre_arc$SampleType <- NA
pre_arc$Receptor <- NA

for(i in c(1:nrow(pt_ls))){
  idx <- which(pre_arc$Sample == pt_ls$ID[i])
  pre_arc$SampleType[idx] <- pt_ls$Type[i]
  pre_arc$Receptor[idx] <- pt_ls$receptor_status[i]
}

#----- quality filter per sample -----#
QualityCells <- lapply(rev(seq(1500, 3000, 100)), function(i){
  out <- pre_arc$Sample[pre_arc$TSSEnrichment > 4 & pre_arc$nFrags > i] %>% table
  return(out)
}) %>% do.call(cbind, .) %>% `colnames<-`(., rev(seq(1500, 3000, 100)))

p1 <- lapply(c(1:nrow(QualityCells)), function(i){
  df <- data.frame(cellnum = QualityCells[i,], fragment = colnames(QualityCells))
  out <- ggplot(df, aes(x = fragment, y = cellnum)) + geom_point() + ggtitle(rownames(QualityCells)[i])
  return(out)
})
pdf("output/Plots/01_Dotplot_QC_cellnumbers.pdf", width = 6, height = 4)
p1
dev.off()

#1500 frags : =<100 cells
sample_1500frags <- which(QualityCells[,16] <= 100) %>% as.numeric()

#others: decreacing rate 80% 
QualityCellsFilter <- lapply(c(1:nrow(QualityCells)), function(i){
  tmp <- QualityCells[i,] / QualityCells[i,16]
  idx <- which(tmp > 0.8)
  out <- colnames(QualityCells)[min(idx)]
  return(out)
}) %>% unlist()

QualityCellsFilter[sample_1500frags] <- 1500
QualityCellsFilter <- as.numeric(QualityCellsFilter)
names(QualityCellsFilter) <- rownames(QualityCells)

FilteredCells <- lapply(names(QualityCellsFilter), function(i){
  filter <- QualityCellsFilter[i]
  out <- pre_arc$cellNames[pre_arc$TSSEnrichment > 4 & pre_arc$nFrags > filter & pre_arc$Sample == i]
  return(out)
}) %>% unlist()

arc <- pre_arc[pre_arc$cellNames %in% FilteredCells]
sort(table(arc$Sample))

#-----------------
# 38 samples
# 22775 cells
#-----------------

saveRDS(pre_arc, "rds/01_pre_arc.rds")
saveRDS(arc, "rds/01_arc.rds")
