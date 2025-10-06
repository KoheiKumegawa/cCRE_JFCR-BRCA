#----------------------------------------------------------------------------
# PalData_Projection.R
#----------------------------------------------------------------------------
library(Seurat)
library(ArchR)
library(ggplot2)
library(scales)
library(ggrastr)
library(parallel)
library(dplyr)
library(viridisLite)

seu_ER <- readRDS("rds/SeuratObject_ERTotal.rds")
seu_HR <- readRDS("rds/SeuratObject_HER2.rds")
seu_TN <- readRDS("rds/SeuratObject_TNBC.rds")

DimPlot(seu_ER, cells.highlight = colnames(seu_ER)[which(seu_ER$seurat_clusters %in% c(0,4,6,5))])
DimPlot(seu_HR, cells.highlight = colnames(seu_HR)[which(seu_HR$seurat_clusters %in% c(0,3,7))])
DimPlot(seu_TN, cells.highlight = colnames(seu_TN)[which(seu_TN$seurat_clusters %in% c(0,2))])

df1 <- data.frame(cellnames = c(colnames(seu_ER)[which(seu_ER$seurat_clusters %in% c(0,4,6,5))],
                                colnames(seu_HR)[which(seu_HR$seurat_clusters %in% c(0,3,7))],
                                colnames(seu_TN)[which(seu_TN$seurat_clusters %in% c(0,2))]),
                  celltypes = "Epithelial")
                  
  
seu_ERSub <- readRDS("rds/SeuratObject_ERTotalSub.rds")
seu_HRSub <- readRDS("rds/SeuratObject_HER2Sub.rds")
seu_TNSub <- readRDS("rds/SeuratObject_TNBCSub.rds")

DimPlot(seu_ERSub, cols = ArchR::ArchRPalettes$stallion)
seu_ERSub_clusters <- c("0"="T","1"="TAMs","2"="CAFs","3"="Pericytes","4"=NA,"5"="Endothelial","6"="TAMs","7"="B","8"="Myeloid","9"="CAFs","10"="Plasma","11"=NA,"12"=NA)
DimPlot(seu_HRSub, cols = ArchR::ArchRPalettes$stallion)
seu_HRSub_clusters <- c("0"="TAMs","1"="T","2"="Plasma","3"="CAFs","4"="Endothelial","5"="B","6"="T","7"="Pericytes","8"=NA,"9"=NA,"10"=NA)
DimPlot(seu_TNSub, cols = ArchR::ArchRPalettes$stallion)
seu_TNSub_clusters <- c("0"="T","1"="TAMs","2"="Plasma","3"="CAFs","4"="T","5"="B","6"="DCs","7"="Endothelial","8"="Pericytes","9"="Myeloid")

df2 <- data.frame(cellnames = c(colnames(seu_ERSub),colnames(seu_HRSub),colnames(seu_TNSub)),
                  celltypes = c(seu_ERSub_clusters[seu_ERSub@meta.data[,"seurat_clusters"]],
                                seu_HRSub_clusters[seu_HRSub@meta.data[,"seurat_clusters"]],
                                seu_TNSub_clusters[seu_TNSub@meta.data[,"seurat_clusters"]]))


df_celltype <- rbind(df1, df2)
rownames(df_celltype) <- df_celltype$cellnames

seu <- readRDS("rds/R_pre_seu.rds")
tmp <- gsub("-","_",colnames(seu))
tmp <- gsub("_1","-1",tmp)
seu <- RenameCells(seu, new.names = tmp)

seu2 <- seu[,df_celltype$cellnames]

#----- Clustering -----#
seu2 <- NormalizeData(seu2, normalization.method = "LogNormalize", scale.factor = 10000)
seu2 <- FindVariableFeatures(seu2, selection.method = "vst", nfeatures = 2000)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
ElbowPlot(seu2, ndims = 50)

seu2 <- FindNeighbors(seu2, dims = 1:50)
seu2 <- RunUMAP(seu2, dims = 1:50)
seu2 <- FindClusters(seu2, resolution = 0.2)


seu2$annotations <- factor(df_celltype$celltypes, levels = c("Epithelial", "T", "B", "Plasma", "TAMs","CAFs", "Endothelial","Myeloid","Pericytes","DCs"))
celltype_colors <- c("#FFB300","#803E75","#817066","#00538A","#8FC31F","#C10020","#A6BDD7","#7DD06F","#844081","#688EC1") %>% 
  `names<-`(.,c("Epithelial", "T", "B", "Plasma", "TAMs","CAFs", "Endothelial","Myeloid","Pericytes","DCs"))

p1 <- DimPlot(seu2, group.by = "annotations", cols = celltype_colors)

#----- Integration of scRNAseq data -----#
arc <- readRDS("rds/01_arc.rds")
arc <- addGeneIntegrationMatrix(
  ArchRProj = arc, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seu2,
  addToArrow = FALSE,
  groupRNA = "annotations",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
arc$annotations <- factor(arc$predictedGroup_Un, levels = names(celltype_colors))

df3 <- arc@embeddings$UMAP$df
df3$annotations <- arc$annotations
colnames(df3) <- c("UMAP1","UMAP2","annotations")
p2 <- ggplot(df3, aes(x=UMAP1,y=UMAP2,color=annotations)) + geom_point_rast(size=0.5) + theme_ArchR() +
  scale_color_manual(values = celltype_colors)

cM <- as.matrix(confusionMatrix(arc$Clusters, arc$predictedGroup_Un))
cM <- reshape2::melt(cM)
cM$Var1 <- factor(cM$Var1, levels=paste0("C", 1:33))

p3 <- ggplot(cM, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") + ArchR::theme_ArchR() + 
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  labs(x = "Cluster", y = "%Sample") + scale_fill_manual(values = celltype_colors) +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))

plotPDF(p1, name = "UMAP_PalData_AuthorsAnno.pdf", ArchRProj = arc, addDOC = FALSE, width = 5.5, height = 5)
plotPDF(p2, name = "UMAP_Proj_PalData_AuthorsAnno.pdf", ArchRProj = arc, addDOC = FALSE, width = 5.5, height = 5)
pdf("output/Plots/Barplot_Proj_PalData_AuthorsAnno.pdf", height = 4, width = 8)
p3
dev.off()

krt_type2 <- paste0("KRT", c(1,2,3,4,5,"6A","6B","6C",7,8,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86))
krt_type1 <- paste0("KRT", c(9,10,12,13,14,15,16,17,18,19,20,23,24,25,26,27,28,31,32,"33A","33B",34,35,36,37,38,39,40))
epi_markers <- c(krt_type1, krt_type2, "EPCAM", "CDH1")
fibroblast_markers <-  c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FBLN1", "FN1", "ACTA2", "TAGLN", "PDGFRB", "FAP")

epifib_markers <- list(epi=epi_markers, fib=fibroblast_markers)
arc <- addModuleScore(arc, useMatrix = "GeneScoreMatrix", features = epifib_markers, name = "score")

plotGroups(arc, groupBy = "Clusters", name = "score.epi", plotAs = "violin")
plotGroups(arc, groupBy = "Clusters", name = "score.fib", plotAs = "violin")

df4 <- data.frame(episcore=sapply(paste0("C",c(1:33)), function(x) mean(arc$score.epi[arc$Clusters == x])),
                  fibscore=sapply(paste0("C",c(1:33)), function(x) mean(arc$score.fib[arc$Clusters == x])),
                  cluster=factor(paste0("C",c(1:33)), levels = paste0("C",c(1:33))))

df4_melt <- reshape2::melt(df4)
p4 <- ggplot(df4, aes(x=reorder(cluster,episcore),y=episcore)) + geom_point() + theme_ArchR() + theme(axis.text.x = element_text(angle = 90,hjust = 1))
p5 <- ggplot(df4, aes(x=reorder(cluster,fibscore),y=fibscore)) + geom_point() + theme_ArchR() + theme(axis.text.x = element_text(angle = 90,hjust = 1))
pdf("output/Plots/Dotplot_MeanEpiFibScore.pdf", height = 3, width = 8)
p4
p5
dev.off()

saveRDS(seu2, "rds/proj_seu2.rds")
saveRDS(arc, "rds/proj_arc.rds")
