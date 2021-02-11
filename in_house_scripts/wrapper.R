knitr::opts_chunk$set(fig.width=12,fig.retina = 6,  message = FALSE, warning = FALSE) 

require(SingleCellExperiment)
require(SC3)
require(Seurat)
require(ggplot2)
require(scales)
require(ggpubr)
require(cowplot)
require(RColorBrewer)
require(networkD3)
require(dplyr)
require(patchwork)
library(reshape2)
library(clusterProfiler)
library(org.Sc.sgd.db)


options(stringsAsFactors = FALSE)

# standardize theme
main <- theme_minimal() + theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=12), text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) 
dotsize = .7
linesize = .3
scale_y_log <- scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

                                                   
# in-house functions 
calculateQCMetrics <- function(object)  ## latest version of scater has defunc this
{
  colData(object) <- setNames(perCellQCMetrics(object), c("total_counts", "total_features_by_counts", "total"))
  rowData(object) <- data.frame(perFeatureQCMetrics(object), log10_total_counts = log10(rowSums(counts(object))))
  return(object)
}  
                                         
calculateCV <- function(mx)
{
  SD <- apply(mx, 1, sd)
  MEAN <- apply(mx, 1, mean)
  CV <- SD/MEAN
  
  return(data.frame(SD, MEAN, CV))
}

plotCVLogMeanCurve <- function(df)
{
  df$CV <- log10(df$CV)
  

    g <- ggplot(
    df,
    aes(x = MEAN, y = CV)
    ) + 
      geom_jitter() + 
      scale_x_log10() +
      geom_smooth(se = FALSE, method = "lm",formula = y ~ splines::bs(x, 3)) + 
      xlab("log10(avg)") + ylab("log10(CV)") + 
      theme_minimal() + theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=12), text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) 

    return(g)
}

plotCVMeanCurve <- function(df)
{
    g <- ggplot(
    df,
    aes(x = MEAN, y = CV)
    ) + 
      geom_jitter() + 
      geom_smooth(se = FALSE, method = "lm",formula = y ~ splines::bs(x, 3)) + 
      xlab("avg") + ylab("CV") + 
      theme_minimal() + theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=12), text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) 

    return(g)
}
                                                   
calculatePCAVariance <- function(seurat.object, assay = "RNA")
{
  mat <- Seurat::GetAssayData(seurat.object, assay = assay, slot = "scale.data")
  pca <- seurat.object[["pca"]]
  total_variance <- sum(matrixStats::rowVars(mat))
  eigValues = (pca@stdev)^2
  return (eigValues / total_variance)
}
          