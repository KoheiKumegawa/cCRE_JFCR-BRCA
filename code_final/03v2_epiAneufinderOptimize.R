#----------------------------------------------------------------------------
# 03v2_epiAneufinderOptimize.R
#----------------------------------------------------------------------------
library(ArchR)
library(scales)
library(ggsignif)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(ggrastr)
addArchRThreads(threads = 24)
addArchRGenome("hg19")

arc <- readRDS("rds/05_arc.rds")
pt_ls <- read.csv("ref/patient_list_v2.csv")
SampleTypeColors <- c("Primary" = "#666666", "LocalRecurrence" = "#1B9E77", "LymphNode" = "#E7298A", "Pleural" = "#7570B3", "Ascitis" = "#E6AB02")
options(scipen=100)

#----- CNVkit bin-level visualization -----#
# P210
d <- data.table::fread("/mnt/host_mnt/Volumes/Shared2/NGS/Analysis_1/20240315_kume1/pb210/pb210_cnvkit_output/pb210.cs.rmdup.cnr.gz")
d <- d[which(d$chromosome %in% paste0("chr", c(1:22, "X"))),]
d$chromosome <- factor(d$chromosome, levels = paste0("chr", c(1:22, "X")))
d <- unique(d[order(d$chromosome, d$start, d$end), ])

P210_binlevel_gr <- GRanges(seqnames = d$chromosome, IRanges(start = d$start+1, end = d$end), log2CR = d$log2)

df1 <- data.frame(log2CR = P210_binlevel_gr$log2CR, position = seq_along(P210_binlevel_gr))
idx1 <- lapply(paste0("chr", c(1:22, "X")), function(i) which(as.character(seqnames(P210_binlevel_gr)) == i)[1]) %>% unlist

p1 <- ggplot(df1, aes(x=position, y=log2CR)) + geom_point_rast(color="gray", size = 0.5) + theme_ArchR() +
  geom_vline(xintercept = c(idx1, length(P210_binlevel_gr)), lty = "dashed", color = "lightgray") +
  geom_hline(yintercept = 0) + labs(x = "Genomic position", y = "log2(Copy ratio)") + ylim(-4, 4)
pdf("output/Plots/03v2_Dotplot_cnvkit_bin_P210.pdf", width = 12, height = 4)
p1
dev.off()

# P190
d <- data.table::fread("/mnt/host_mnt/Volumes/Shared2/NGS/Analysis_1/20240315_kume1/pb190/pb190_cnvkit_output/pb190.cs.rmdup.cnr.gz")
d <- d[which(d$chromosome %in% paste0("chr", c(1:22, "X"))),]
d$chromosome <- factor(d$chromosome, levels = paste0("chr", c(1:22, "X")))
d <- unique(d[order(d$chromosome, d$start, d$end), ])

P190_binlevel_gr <- GRanges(seqnames = d$chromosome, IRanges(start = d$start+1, end = d$end), log2CR = d$log2)

df2 <- data.frame(log2CR = P190_binlevel_gr$log2CR, position = seq_along(P190_binlevel_gr))
idx2 <- lapply(paste0("chr", c(1:22, "X")), function(i) which(as.character(seqnames(P190_binlevel_gr)) == i)[1]) %>% unlist

p2 <- ggplot(df2, aes(x=position, y=log2CR)) + geom_point_rast(color="gray", size = 0.5) + theme_ArchR() +
  geom_vline(xintercept = c(idx2, length(P190_binlevel_gr)), lty = "dashed", color = "lightgray") +
  geom_hline(yintercept = 0) + labs(x = "Genomic position", y = "log2(Copy ratio)") + ylim(-4, 4)
pdf("output/Plots/03v2_Dotplot_cnvkit_bin_P190.pdf", width = 12, height = 4)
p2
dev.off()

P190_chr8_binlevel_gr <- P190_binlevel_gr[which(seqnames(P190_binlevel_gr)=="chr8")]
P190_chr8_binlevel_gr$MYC <- "N"
P190_chr8_binlevel_gr$MYC[queryHits(findOverlaps(P190_chr8_binlevel_gr, geneAnnoHg19$genes[which(geneAnnoHg19$genes$symbol=="MYC")]))] <- "Y"

df2.2 <- data.frame(log2CR = P190_chr8_binlevel_gr$log2CR, position = seq_along(P190_chr8_binlevel_gr), MYC =P190_chr8_binlevel_gr$MYC)
df2.2 <- df2.2[order(df2.2$MYC),]
p2.2 <- ggplot(df2.2, aes(x=position, y=log2CR, color=MYC)) + geom_point_rast(size = 0.5) + theme_ArchR() +
  scale_color_manual(values = c("N"="gray","Y"="red")) + 
  geom_hline(yintercept = 0) + labs(x = "Genomic position", y = "log2(Copy ratio)") + ylim(-4, 4)
pdf("output/Plots/03v2_Dotplot_cnvkit_bin_P190_chr8_MYChighlight.pdf", width = 8, height = 4)
p2.2
dev.off()


#----- epiAneufinder summary -----#
### P210
cnv_se_p210 <- lapply(paste0(c(1:10), "Mb"), function(i){
  dt <- data.table::fread(paste0("/mnt/host_mnt/Volumes/Shared/Kume/Analysis/S23-3_scATAC-BRCA/output/P210LR_wdsize2/",
                                 i, "/epiAneufinder_results/results_table.tsv")) %>% data.frame
  gr <- GRanges(seqnames = dt$seq, IRanges(start=dt$start, end=dt$end))
  dt <- dt[,c(5:ncol(dt))]
  colnames(dt) <- gsub("cell.", "", colnames(dt))
  colnames(dt) <- paste0("P210LR#", colnames(dt))
  
  dt <- dt[, intersect(colnames(dt), arc$cellNames)]
  se <- SummarizedExperiment(assays = list(CNV = dt), rowRanges = gr)
  
  return(se)
})

#500kb
dt <- data.table::fread(paste0("/mnt/host_mnt/Volumes/Shared/Kume/Analysis/S23-3_scATAC-BRCA/output/epiAneufinder_results/P210LR/epiAneufinder_results/results_table.tsv")) %>% data.frame
gr <- GRanges(seqnames = dt$seq, IRanges(start=dt$start, end=dt$end))
dt <- dt[,c(5:ncol(dt))]
colnames(dt) <- gsub("cell.", "", colnames(dt))
colnames(dt) <- paste0("P210LR#", colnames(dt))
dt <- dt[, intersect(colnames(dt), arc$cellNames)]
se <- SummarizedExperiment(assays = list(CNV = dt), rowRanges = gr)
cnv_se_p210 <- c(list(se), cnv_se_p210)

#100kb
dt <- data.table::fread(paste0("/mnt/host_mnt/Volumes/Shared/Kume/Analysis/S23-3_scATAC-BRCA/output/epiAneufinder_results_default/P210LR/epiAneufinder_results/results_table.tsv")) %>% data.frame
gr <- GRanges(seqnames = dt$seq, IRanges(start=dt$start, end=dt$end))
dt <- dt[,c(5:ncol(dt))]
colnames(dt) <- gsub("cell.", "", colnames(dt))
colnames(dt) <- paste0("P210LR#", colnames(dt))
dt <- dt[, intersect(colnames(dt), arc$cellNames)]
se <- SummarizedExperiment(assays = list(CNV = dt), rowRanges = gr)
cnv_se_p210 <- c(list(se), cnv_se_p210)

names(cnv_se_p210) <- c("100kb","500kb",paste0(c(1:10), "Mb"))

cnv_se_p210 <- lapply(cnv_se_p210, function(x){
  out <- x[,arc$cellNames[which(arc$Sample == "P210LR")]]
  return(out)
})

idy1 <- arc[colnames(cnv_se_p210[[1]])]$Clusters2
idy1[idy1 != "Ep23"] <- "Others"
idy1[idy1 == "Ep23"] <- "Ep23"
idy1 <- factor(idy1, levels = c("Ep23", "Others"))

p3 <- lapply(names(cnv_se_p210), function(i){
  mtx <- assay(cnv_se_p210[[i]]) %>% t()
  idx3 <- rowRanges(cnv_se_p210[[i]]) %>% seqnames()
  col_fun2 <- colorRamp2(c(0,1,2), c("blue", "lightgray", "red"))
  ht2 <- Heatmap(mtx, name = paste0(i, "(", ncol(mtx),")"), col = col_fun2, 
                 column_title_gp = gpar(fontsize = 6), row_title_gp = gpar(fontsize = 6),
                 cluster_columns = F, cluster_rows = F, 
                 show_row_names = F, column_split = idx3, row_split = idy1, use_raster = T)
  out <- draw(ht2)
  return(out)
})

pdf("output/Plots/03v2_Heatmap_epiAneufinder_multipleWSize_P210.pdf", width = 12, height = 4)
p3
dev.off()

### P190
cnv_se_p190 <- lapply(paste0(c(1:10), "Mb"), function(i){
  dt <- data.table::fread(paste0("/mnt/host_mnt/Volumes/Shared2/Kume/Analysis/S23-3_scATAC-BRCA/output/P190_wdsize/",
                                 i, "/epiAneufinder_results/results_table.tsv")) %>% data.frame
  gr <- GRanges(seqnames = dt$seq, IRanges(start=dt$start, end=dt$end))
  dt <- dt[,c(5:ncol(dt))]
  colnames(dt) <- gsub("cell.", "", colnames(dt))
  colnames(dt) <- paste0("P190#", colnames(dt))
  
  dt <- dt[, intersect(colnames(dt), arc$cellNames)]
  se <- SummarizedExperiment(assays = list(CNV = dt), rowRanges = gr)
  
  return(se)
})

#500kb
dt <- data.table::fread(paste0("/mnt/host_mnt/Volumes/Shared/Kume/Analysis/S23-3_scATAC-BRCA/output/epiAneufinder_results/P190/epiAneufinder_results/results_table.tsv")) %>% data.frame
gr <- GRanges(seqnames = dt$seq, IRanges(start=dt$start, end=dt$end))
dt <- dt[,c(5:ncol(dt))]
colnames(dt) <- gsub("cell.", "", colnames(dt))
colnames(dt) <- paste0("P190#", colnames(dt))
dt <- dt[, intersect(colnames(dt), arc$cellNames)]
se <- SummarizedExperiment(assays = list(CNV = dt), rowRanges = gr)
cnv_se_p190 <- c(list(se), cnv_se_p190)

#100kb
dt <- data.table::fread(paste0("/mnt/host_mnt/Volumes/Shared/Kume/Analysis/S23-3_scATAC-BRCA/output/epiAneufinder_results_default/P190/epiAneufinder_results/results_table.tsv")) %>% data.frame
gr <- GRanges(seqnames = dt$seq, IRanges(start=dt$start, end=dt$end))
dt <- dt[,c(5:ncol(dt))]
colnames(dt) <- gsub("cell.", "", colnames(dt))
colnames(dt) <- paste0("P190#", colnames(dt))
dt <- dt[, intersect(colnames(dt), arc$cellNames)]
se <- SummarizedExperiment(assays = list(CNV = dt), rowRanges = gr)
cnv_se_p190 <- c(list(se), cnv_se_p190)

names(cnv_se_p190) <- c("100kb","500kb",paste0(c(1:10), "Mb"))

cnv_se_p190 <- lapply(cnv_se_p190, function(x){
  out <- x[,arc$cellNames[which(arc$Sample == "P190")]]
  return(out)
})

idy2 <- arc[colnames(cnv_se_p190[[1]])]$Clusters2
idy2[idy2 %ni% c("Ep5","Ep6")] <- "Others"
idy2 <- factor(idy2, levels = c("Ep5","Ep6","Others"))

p4 <- lapply(names(cnv_se_p190), function(i){
  mtx <- assay(cnv_se_p190[[i]]) %>% t()
  idx3 <- rowRanges(cnv_se_p190[[i]]) %>% seqnames()
  col_fun2 <- colorRamp2(c(0,1,2), c("blue", "lightgray", "red"))
  ht2 <- Heatmap(mtx, name = paste0(i, "(", ncol(mtx),")"), col = col_fun2, 
                 column_title_gp = gpar(fontsize = 6), row_title_gp = gpar(fontsize = 6),
                 cluster_columns = F, cluster_rows = F, 
                 show_row_names = F, column_split = idx3, row_split = idy2, use_raster = T)
  out <- draw(ht2)
  return(out)
})
pdf("output/Plots/03v2_Heatmap_epiAneufinder_multipleWSize_P190.pdf", width = 12, height = 4)
p4
dev.off()

#----- CNVkit VS epiAneufinder, window size -----#
cnv_call_p210 <- data.table::fread("/mnt/host_mnt/Volumes/Shared2/NGS/Analysis_1/20240315_kume1/pb210/pb210_cnvkit_output/pb210.cs.rmdup_CNV_CALLS.bed")
cnv_call_p210 <- cnv_call_p210[which(cnv_call_p210$V5 > 8)]
cnv_call_p210_gr <- GRanges(seqnames = cnv_call_p210$V1, IRanges(start = cnv_call_p210$V2+1, end = cnv_call_p210$V3))

df5 <- lapply(names(cnv_se_p210), function(x){
  print(x)
  obj <- cnv_se_p210[[x]]
  gr <- rowRanges(obj)

  HighCN_region <- queryHits(findOverlaps(gr, cnv_call_p210_gr))
  
  mtx <- assay(obj[, which(idy1 == "Ep23")]) %>% t
  met_x <- apply(mtx, 1, function(y) sum(y == 2))
  
  if (length(HighCN_region) == 1) {
    met_y <- mtx[, HighCN_region]
    if (is.null(dim(met_y))) {  # met_y is a vector
      met_y <- as.matrix(met_y)
    }
    met_y <- apply(met_y, 1, function(y) sum(y == 2))
  } else if (length(HighCN_region) > 1) {
    met_y <- apply(mtx[, HighCN_region, drop = FALSE], 1, function(y) sum(y == 2))
  } else {
    return(NA)
  }
  out <- met_y / met_x
  return(out)
})
df5 <- do.call(cbind,df5)
df5_mean <- colMeans(df5, na.rm = T)
names(df5_mean) <- names(cnv_se_p210)

df5.2 <- lapply(cnv_se_p210, function(x){
  gr <- rowRanges(x)
  HighCN_region <- findOverlaps(gr, cnv_call_p210_gr) %>% queryHits
  mtx <- assay(x[,which(idy1 == "Others")]) %>% t
  met_x <- ncol(mtx)
  met_y <- apply(mtx, 1, function(x) which(x == 1) %>% length)
  out <- met_y / met_x
  return(out)
})
df5.2 <- do.call(cbind,df5.2)
df5.2_mean <- colMeans(df5.2, na.rm = T)
plot(df5.2_mean)

df5.3 <- data.frame(accuracy = df5_mean, norm = df5.2_mean, window = factor(names(df5_mean), levels = names(df5_mean)))
df5.3 <- reshape2::melt(df5.3)
p5 <- ggplot(df5.3[which(df5.3$variable == "accuracy"),], aes(x = window, y = value, group = 1)) + 
  geom_point(color = "red") + geom_line(color = "red") + theme_ArchR() + 
  labs(x = "Window Size", y = "Accuracy index")
p5.2 <- ggplot(df5.3[which(df5.3$variable == "norm"),], aes(x = window, y = value, group = 1)) + 
  geom_point(color = "blue") + geom_line(color = "blue") + theme_ArchR() +
  labs(x = "Window Size", y = "Normal window ratio")

pdf("output/Plots/03v2_Boxplot_epiAneufinder_optimize_pb210.pdf", width = 6, height = 3)
p5
p5.2
dev.off()

#----- P190 -----#
arc_epi <- readRDS("rds/04_arc_epi.rds")
arc_epi_v2 <- arc_epi[arc_epi$Clusters %in% c("Ep5","Ep6")]
arc_epi_v2 <- addImputeWeights(arc_epi_v2)
p5.3 <- plotGroups(arc_epi_v2, groupBy = "Clusters", name = "MYC", colorBy = "GeneScoreMatrix", plotAs = "violin", addBoxplot=T)
pdf("output/Plots/03v2_Boxplot_MYC_GeneScore_Ep5Ep6.pdf", width = 3, height = 3)
p5.3
dev.off()

cnv_call_p190 <- data.table::fread("/mnt/host_mnt/Volumes/Shared2/NGS/Analysis_1/20240315_kume1/pb190/pb190_cnvkit_output/pb190.cs.rmdup_CNV_CALLS.bed")
cnv_call_p190 <- cnv_call_p190[which(cnv_call_p190$V5 > 8)]
cnv_call_p190_gr <- GRanges(seqnames = cnv_call_p190$V1, IRanges(start = cnv_call_p190$V2+1, end = cnv_call_p190$V3), CN=cnv_call_p190$V5)

df6 <- lapply(names(cnv_se_p190), function(x){
  print(x)
  obj <- cnv_se_p190[[x]]
  gr <- rowRanges(obj)
  
  HighCN_region <- queryHits(findOverlaps(gr, cnv_call_p190_gr))
  
  mtx <- assay(obj[, which(idy2 != "Others")]) %>% t
  met_x <- apply(mtx, 1, function(y) sum(y == 2))
  
  if (length(HighCN_region) == 1) {
    met_y <- mtx[, HighCN_region]
    if (is.null(dim(met_y))) {  # met_y is a vector
      met_y <- as.matrix(met_y)
    }
    met_y <- apply(met_y, 1, function(y) sum(y == 2))
  } else if (length(HighCN_region) > 1) {
    met_y <- apply(mtx[, HighCN_region, drop = FALSE], 1, function(y) sum(y == 2))
  } else {
    return(NA)
  }
  
  out <- met_y / met_x
  return(out)
})
df6 <- do.call(cbind,df6)
df6_mean <- colMeans(df6, na.rm = T)
plot(df6_mean)

df6.2 <- lapply(cnv_se_p190, function(x){
  gr <- rowRanges(x)
  HighCN_region <- findOverlaps(gr, cnv_call_p190_gr) %>% queryHits
  mtx <- assay(x[,which(idy2 == "Others")]) %>% t
  met_x <- ncol(mtx)
  met_y <- apply(mtx, 1, function(x) which(x == 1) %>% length)
  out <- met_y / met_x
  return(out)
})
df6.2 <- do.call(cbind,df6.2)
df6.2_mean <- colMeans(df6.2, na.rm = T)
plot(df6.2_mean)

#----- epiAneufinder vs cnvkit -----#
d <- data.table::fread("output/output_bed/03_cnvkit_HighCN_P210LR_hg19Liftover.bed") %>% data.frame
d <- d[which(d$V1 %in% paste0("chr", c(1:22, "X"))),]
d$V1 <- factor(d$V1, levels = paste0("chr", c(1:22, "X")))
d <- unique(d[order(d$V1, d$V2, d$V3), ])
cnvkit_highCN_gr <- GRanges(seqnames = d$V1, IRanges(start = d$V2+1, end = d$V3), log2CR = d$V4)

df3 <- lapply(cnv_se_P210LR, function(x){
  gr <- rowRanges(x)
  HighCN_region <- findOverlaps(gr, cnvkit_highCN_gr) %>% queryHits
  mtx <- assay(x[,which(idy1 == "Tumor_cell_C20")]) %>% t
  met_x <- apply(mtx, 1, function(x) which(x == 2) %>% length)
  met_y <- apply(mtx[,HighCN_region], 1, function(x) which(x == 2) %>% length)
  out <- met_y / met_x
  return(out)
})
df3 <- do.call(cbind,df3)
df3_mean <- colMeans(df3, na.rm = T)
df3 <- reshape2::melt(df3)

p4 <- ggplot(df3, aes(x = Var2, y = value)) + geom_boxplot(outlier.shape = 3) + theme_ArchR() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Window size", y = "Signal-noise ratio (DNA-seq gain window / All epiAneufinder gain window)")

df4 <- lapply(cnv_se_P210LR, function(x){
  gr <- rowRanges(x)
  HighCN_region <- findOverlaps(gr, cnvkit_highCN_gr) %>% queryHits
  mtx <- assay(x[,which(idy1 == "Others")]) %>% t
  met_x <- ncol(mtx)
  met_y <- apply(mtx, 1, function(x) which(x == 1) %>% length)
  out <- met_y / met_x
  return(out)
})
df4 <- do.call(cbind,df4)
df4_mean <- colMeans(df4, na.rm = T)
df4 <- reshape2::melt(df4)

p5 <- ggplot(df4, aes(x = Var2, y = value)) + geom_boxplot(outlier.shape = 3) + theme_ArchR() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Window size", y = "Normal window ratio (epiAneufinder normal window in normal cells / all windows)")

df5 <- data.frame(SNratio = df3_mean, Norm = df4_mean, window = factor(names(df3_mean), levels = names(df3_mean)))
df5 <- reshape2::melt(df5)
p6 <- ggplot(df5[which(df5$variable == "SNratio"),], aes(x = window, y = value, group = 1)) + 
  geom_point(color = "red") + geom_line(color = "red") + theme_ArchR() + 
  labs(x = "Window Size", y = "Signal-noise Ratio")
p7 <- ggplot(df5[which(df5$variable == "Norm"),], aes(x = window, y = value, group = 1)) + 
  geom_point(color = "blue") + geom_line(color = "blue") + theme_ArchR() +
  labs(x = "Window Size", y = "Normal window ratio")

pdf("output/Plots/03_Boxplot_epiAneufinder_optimize.pdf", width = 8, height = 4)
p4
p5
p6
p7
dev.off()
