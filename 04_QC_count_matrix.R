library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(scater)
library(cowplot)

source("in_house_scripts/wrapper.R")

mx.sample <- read.delim2("count_matrix/sample_counts.tsv", header = T, row.names = 1) # Replace path by count matrix generated in Step 3

sce.sample <- SingleCellExperiment(assays = list(counts = as.matrix(mx.sample))

umi.qc <- sce.sample

umi.qc <- calculateQCMetrics(
    umi.qc
)

par(mfrow=c(1,2))
hist(
    umi.qc$total_counts,
    breaks = 150,xlab = c("Total counts"), main = ""
)
hist(
    umi.qc$total_features_by_counts,
    breaks = 150,xlab = c("Total features counts"), main = ""
)

## Thresholds are set manually according to the graph
hist(
    umi.qc$total_counts,
    breaks = 150, ylim = c(0,25),xlab = c("Total counts"), main = "Sample"
); abline(v = 200, col = "red"); abline(v = 2200, col = "red")

filter_by_total_counts <- (umi.qc$total_counts >= 200 & umi.qc$total_counts <= 2200)
table(filter_by_total_counts)

## Thresholds are set manually according to the graph
hist(
    umi.qc$total_features_by_counts,
    breaks = 150, ylim = c(0,15),xlab = c("Total features counts"), main = "Sample"
); abline(v = 100, col = "red"); abline(v = 1600, col = "red")

filter_by_total_features_by_counts <- (umi.qc$total_features_by_counts >= 100 & umi.qc$total_features_by_counts <= 1600)
table(filter_by_total_features_by_counts)

umi.qc$use <- (
    filter_by_total_features_by_counts &
    filter_by_total_counts
)


keep_feature <- nexprs(
  umi.qc[,colData(umi.qc)$use], 
  byrow = TRUE, 
  detection_limit = 1
) >= 1
rowData(umi.qc)$use <- keep_feature

table(keep_feature)


plotHighestExprs(umi.qc, exprs_values = "counts")

qc.seurat <- CreateSeuratObject(assay = "RNA", counts = as.matrix(home_mx[,intersect(colnames(home_mx), colnames(dge_mx))]), project = "YeastDropseq_Sample")

plot_grid(
  VlnPlot(qc.seurat, features = c("nFeature_RNA"), ncol = 1) + NoLegend(),
  VlnPlot(qc.seurat, features = c("nCount_RNA"), ncol = 1) + NoLegend(),
  FeatureScatter(qc.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend(), 
  ncol = 3
)

plot_grid(
  VlnPlot(qc.seurat[rowData(umi.qc)$use, colData(umi.qc)$use], features = c("nFeature_RNA"), ncol = 1) + NoLegend(),
  VlnPlot(qc.seurat[rowData(umi.qc)$use, colData(umi.qc)$use], features = c("nCount_RNA"), ncol = 1) + NoLegend(),
  FeatureScatter(qc.seurat[rowData(umi.qc)$use, colData(umi.qc)$use], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend(), 
  ncol = 3
)

qc.seurat <- NormalizeData(qc.seurat, normalization.method = "LogNormalize", scale.factor = 10000)

qc.seurat <- FindVariableFeatures(qc.seurat, selection.method = "vst", nfeatures = 100)

all.genes <- rownames(qc.seurat)
qc.seurat <- ScaleData(qc.seurat, features = all.genes)

qc.seurat <- RunPCA(qc.seurat, features = VariableFeatures(object = qc.seurat))

sce.qc <- as.SingleCellExperiment(qc.seurat)
sce.qc$use <- umi.qc$use
mat <- Seurat::GetAssayData(qc.seurat, assay = "RNA", slot = "scale.data")
pca <- qc.seurat[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  
attributes(reducedDim(sce.qc, "PCA") )$percentVar <- eigValues / total_variance


plotPCA(sce.qc, colour_by = "use", size_by = "nFeature_RNA")
DimHeatmap(qc.seurat, dims = 1:2, cells = ncol(qc.seurat), balanced = TRUE)

sce.qc <- sce.qc[rowData(umi.qc)$use, colData(umi.qc)$use]
sce.qc.sample <- calculateQCMetrics(sce.qc)


dim(sce.qc.sample)
dim(sce.qc)

## Repeated the code above for each sample and generate the sce.qc.sample object for the analysis below:

sce.qc.dmso <- calculateQCMetrics(sce.qc.dmso)
rowData(sce.qc.dmso)$cell_feature_count <- rowSums(as.matrix(counts(sce.qc.dmso) > 0 + 0))

sce.qc.guan <- calculateQCMetrics(sce.qc.guan)
rowData(sce.qc.guan)$cell_feature_count <- rowSums(as.matrix(counts(sce.qc.guan) > 0 + 0))

sce.qc.mpa <- calculateQCMetrics(sce.qc.mpa)
rowData(sce.qc.mpa)$cell_feature_count <- rowSums(as.matrix(counts(sce.qc.mpa) > 0 + 0))

sce.qc.gmpa <- calculateQCMetrics(sce.qc.gmpa)
rowData(sce.qc.gmpa)$cell_feature_count <- rowSums(as.matrix(counts(sce.qc.gmpa) > 0 + 0))

df.sce.feature <- data.frame(rbind(
  data.frame(count = log10(rowData(sce.qc.dmso)$cell_feature_count), condition = "DMSO"),
  data.frame(count = log10(rowData(sce.qc.guan)$cell_feature_count), condition = "Guanine"),
  data.frame(count = log10(rowData(sce.qc.mpa)$cell_feature_count), condition = "MPA"),
  data.frame(count = log10(rowData(sce.qc.gmpa)$cell_feature_count), condition = "Guanine + MPA")
))


g1 <- ggplot(
 df.sce.feature, 
  aes(x = count, colour = condition)
) + theme_bw() + theme(legend.position = "none")  + 
  geom_histogram(alpha=.2, binwidth=.05) + facet_grid(vars(condition))  + labs(y = "Frequency", x = "Feature counts")  #+ xlim(0, 60) + ylim(0, 800)


df.sce.feature <- data.frame(rbind(
  data.frame(count = rowData(sce.qc.dmso)$log10_total_counts, condition = "DMSO"),
  data.frame(count = rowData(sce.qc.guan)$log10_total_counts, condition = "Guanine"),
  data.frame(count = rowData(sce.qc.mpa)$log10_total_counts, condition = "MPA"),
  data.frame(count = rowData(sce.qc.gmpa)$log10_total_counts, condition = "Guanine + MPA")
))

g2 <- ggplot(
 df.sce.feature, 
  aes(x = count, colour = condition)
) + theme_bw() + theme(legend.position = "none")  +
  geom_histogram(alpha=.2, binwidth=.05) + facet_grid(vars(condition)) + labs(y = "Frequency", x = "Counts") 

plot_grid(
  g1, g2, ncol = 2, nrow = 1
)
                                   

                                   
# Calculate CV per sample
df.dmso <- calculateCV(counts(sce.qc.dmso))
df.guan <- calculateCV(counts(sce.qc.guan))
df.mpa <- calculateCV(counts(sce.qc.mpa))
df.gmpa <- calculateCV(counts(sce.qc.gmpa))

# Plot summary
df <- rbind(
  data.frame(df.dmso, sample = "DMSO"),
  data.frame(df.guan, sample = "Guanine"),
  data.frame(df.mpa, sample = "MPA"),
  data.frame(df.gmpa, sample = "Guanine+MPA"))
df$CV <- log10(df$CV)

ggplot(
    df,
    aes(x = MEAN, y = CV)
    ) + 
      geom_jitter() + 
      xlab("avg") + ylab("CV") + main + facet_wrap(~ sample,ncol = 2)

ggplot(
    df,
    aes(x = MEAN, y = CV)
    ) + 
      geom_jitter() + 
      scale_x_log10() +
      geom_smooth(se = FALSE, method = "lm",formula = y ~ splines::bs(x, 3)) + 
      xlab("log10(avg)") + ylab("log10(CV)") + main + facet_wrap(~ sample,ncol = 2)

# summary of CV
summary(df.dmso$CV)
summary(df.guan$CV)
summary(df.mpa$CV)
summary(df.gmpa$CV)
