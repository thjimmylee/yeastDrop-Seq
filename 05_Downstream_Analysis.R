source("in_house_scripts/wrapper.R")

# Retrieve raw count matrix from GEO Browser (Accession Number: GSE165686)

# creating seurat object
seurat.dmso <- readRDS("~/data/yeastdropseq/seurat/seurat_DMSO.rds")
seurat.guan <- readRDS("~/data/yeastdropseq/seurat/seurat_Guanine.rds")
seurat.mpa <- readRDS("~/data/yeastdropseq/seurat/seurat_MPA.rds")
seurat.gmpa <- readRDS("~/data/yeastdropseq/seurat/seurat_Guanine-MPA.rds")

seurat.dmso@meta.data$sample <- "D"
seurat.guan@meta.data$sample <- "G"
seurat.mpa@meta.data$sample <- "M"
seurat.gmpa@meta.data$sample <- "MG"

seurat.yds <- merge(seurat.dmso, seurat.guan, project = "yeastdropseq")
seurat.yds <- merge(seurat.yds, seurat.mpa, project = "yeastdropseq")
seurat.yds <- merge(seurat.yds, seurat.gmpa, project = "yeastdropseq")


# preprocessing seurat object
seurat.yds <- NormalizeData(seurat.yds, normalization.method = "LogNormalize", scale.factor = 10000 ,assay = "RNA")

seurat.yds <- FindVariableFeatures(seurat.yds, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seurat.yds)
seurat.yds <- ScaleData(seurat.yds, features = all.genes)

seurat.yds <- RunPCA(seurat.yds, features = VariableFeatures(object = seurat.yds))
seurat.yds <- RunUMAP(seurat.yds, reduction = "pca", dims = 1:40, n.neighbors = 4, min.dist = .05)


# Plotting reduced dimensions
plotly::plot_ly(data = data.frame(sample = seurat.yds@meta.data$sample, seurat.yds@reductions$pca@cell.embeddings), x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~sample,size = .3)

pcplot_before <- DimPlot(seurat.yds, reduction = "pca",group.by = "sample") + 
  xlab(paste0("PC1 (", percent(calculatePCAVariance(seurat.yds), accuracy = .1)[1], ")")) + 
  ylab(paste0("PC2 (", percent(calculatePCAVariance(seurat.yds), accuracy = .1)[2], ")"))

plot_grid(
  pcplot_before + main,
  DimPlot(seurat.yds, reduction = "umap",group.by = "sample") + main
)

### Running SC3 on all sample

###############################
### SC3 routine starts here ###
###############################

# create SCEset 
sce.yds <- as.SingleCellExperiment(seurat.yds, assay = "RNA")
rowData(sce.yds)$feature_symbol <- rownames(sce.yds)

sce.yds <- sc3_estimate_k(sce.yds)
str(metadata(sce.yds)$sc3)

counts(sce.yds) <- as.matrix(counts(sce.yds))
logcounts(sce.yds) <- as.matrix(logcounts(sce.yds))

sce.yds <- sc3(sce.yds, ks = 2:10, biology = TRUE)

attributes(reducedDim(sce.yds, "PCA") )$percentVar <- eigValues / total_variance

# SC3 analysis
k = 4
plot_grid(
  ggplot(
    data.frame(
      reducedDim(sce.yds, "PCA")[,c(1,2)],
      sample = sce.yds$sample,
      sc3_cluster = colData(sce.yds)[,paste0("sc3_", k, "_clusters")],
      sc3_outlier = colData(sce.yds)[,paste0("sc3_", k, "_log2_outlier_score")]
      ), aes( x = PC_1, y = PC_2, colour = sample)
  ) + geom_point() + main + 
    xlab(paste0("PC1 (", percent(attributes(reducedDim(sce.yds, "PCA"))$percentVar, .1)[1], ")")) +
    ylab(paste0("PC2 (", percent(attributes(reducedDim(sce.yds, "PCA"))$percentVar, .1)[2], ")")),
  
  
  ggplot(
    data.frame(
      reducedDim(sce.yds, "PCA")[,c(1,2)],
      sample = sce.yds$sample,
      sc3_cluster = colData(sce.yds)[,paste0("sc3_", k, "_clusters")],
      sc3_outlier = colData(sce.yds)[,paste0("sc3_", k, "_log2_outlier_score")]
      ), aes( x = PC_1, y = PC_2, colour = sc3_cluster, size = sc3_outlier)
  ) + geom_point() + main + 
    xlab(paste0("PC1 (", percent(attributes(reducedDim(sce.yds, "PCA"))$percentVar, .1)[1], ")")) +
    ylab(paste0("PC2 (", percent(attributes(reducedDim(sce.yds, "PCA"))$percentVar, .1)[2], ")"))

)

# Sankey plot
nodes = data.frame(name = c(c("D", "G", "M", "GM"),paste("sc3_4", as.character(1:4)) ) )

df <- data.frame(table(paste0(sce.yds$sample, "=", sce.yds$sc3_4_clusters)))

df$source <- as.character(sub("=.*", "", df$Var1))
df$source <- match(df$source, nodes$name[1:4]) - 1
df$target <- as.numeric(sub(".*=", "", df$Var1)) - 1 + 4

sankeyNetwork(Links = df, Nodes = nodes,
 Source = "source", Target = "target",
 Value = "Freq", NodeID = "name", fontSize= 12)

sc3_plot_consensus(
      sce.yds, k = k,
      show_pdata = c(
          "sample", 
          paste0("sc3_", k, "_clusters"), 
          paste0("sc3_", k, "_log2_outlier_score")
      )
  )
sc3_plot_silhouette(sce.yds, k = k)

sc3_plot_expression(
      sce.yds, k = k,
      show_pdata = c(
          "sample", 
          paste0("sc3_", k, "_clusters"), 
          paste0("sc3_", k, "_log2_outlier_score")
      )
  )
sc3_plot_cluster_stability(sce.yds, k = k)

# reporting cell number in each cluster
table(colData(sce.yds)[,paste0("sc3_", k, "_clusters")])

#############################
### SC3 routine ends here ###
#############################

### Running SC3 on each sample

# repeat pre-processing for each sample as following:
seurat.sample <- NormalizeData(seurat.sample, normalization.method = "LogNormalize", scale.factor = 10000 ,assay = "RNA")

seurat.sample <- FindVariableFeatures(seurat.sample, selection.method = "vst", nfeatures = 500)

all.genes <- rownames(seurat.sample)
seurat.seurat.sampledmso <- ScaleData(seurat.sample, features = all.genes)

seurat.sample <- RunPCA(seurat.sample, features = VariableFeatures(object = seurat.sample))

# with this sample specific seurat object, perform the SC3 routine as showed above

# plot heatmap that shows the marker genes in each cluster
sc3_plot_markers(
    sce.sample, k = k,
    show_pdata = c(
        paste0("sc3_", k, "_clusters"),
        paste0("sc3_", k, "_log2_outlier_score")
    ),auroc = .75 # fine tune this for a better result.
)

### Re-labeling and plotting the 8 clusters:
sce.dmso$final_clusters <- paste0("D.", sce.dmso$sc3_2_clusters)
sce.guan$final_clusters <- paste0("G.", sce.guan$sc3_2_clusters)
sce.mpa$final_clusters <- "M.1"
sce.gmpa$final_clusters <- paste0("MG.", sce.gmpa$sc3_3_clusters)
df.final_clusters <- rbind(colData(sce.dmso)[,"final_clusters", drop = F],
      colData(sce.guan)[,"final_clusters", drop = F],
      colData(sce.mpa)[,"final_clusters", drop = F],
      colData(sce.gmpa)[,"final_clusters", drop = F])
sce.yds$sc3_clusters <- df.final_clusters$final_clusters

plotly::plot_ly(data = data.frame(sample =sce.yds$sc3_clusters, reducedDim(sce.yds, "PCA")[,c(1:3)]), x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~sample,size = .3)

plot_grid(
  ggplot(
  data.frame(
    reducedDim(sce.yds, "PCA")[,c(1,2)],
    sc3_clusters = sce.yds$sc3_clusters
    ), aes( x = PC_1, y = PC_2, colour = sc3_clusters)
) + geom_point() + main + 
  xlab(paste0("PC1 (", percent(attributes(reducedDim(sce.yds, "PCA"))$percentVar, .1)[1], ")")) +
  ylab(paste0("PC2 (", percent(attributes(reducedDim(sce.yds, "PCA"))$percentVar, .1)[2], ")")),

ggplot(
  data.frame(
    reducedDim(sce.yds, "UMAP"),
    sc3_clusters = sce.yds$sc3_clusters
    ), aes( x = UMAP_1, y = UMAP_2, colour = sc3_clusters)
) + geom_point() + main + 
  xlab("UMAP1") +
  ylab("UMAP2")  
)

table(sce.yds$sc3_clusters)

k = 8
plot_grid(
  ggplot(
  data.frame(
    reducedDim(sce.yds, "PCA")[,c(1,2)],
    sc3_clusters = sce.yds$sc3_clusters
    ), aes( x = PC_1, y = PC_2, colour = sc3_clusters)
) + geom_point() + main + 
  xlab(paste0("PC1 (", percent(attributes(reducedDim(sce.yds, "PCA"))$percentVar, .1)[1], ")")) +
  ylab(paste0("PC2 (", percent(attributes(reducedDim(sce.yds, "PCA"))$percentVar, .1)[2], ")")),

  
  ggplot(
    data.frame(
      reducedDim(sce.yds, "PCA")[,c(1,2)],
      sample = sce.yds$sample,
      sc3_cluster = colData(sce.yds)[,paste0("sc3_", k, "_clusters")],
      sc3_outlier = colData(sce.yds)[,paste0("sc3_", k, "_log2_outlier_score")]
      ), aes( x = PC_1, y = PC_2, colour = sc3_cluster, size = sc3_outlier)
  ) + geom_point() + main + 
    xlab(paste0("PC1 (", percent(attributes(reducedDim(sce.yds, "PCA"))$percentVar, .1)[1], ")")) +
    ylab(paste0("PC2 (", percent(attributes(reducedDim(sce.yds, "PCA"))$percentVar, .1)[2], ")"))

)


plot_grid(
  ggplot(
  data.frame(
    reducedDim(sce.yds, "UMAP"),
    sc3_clusters = sce.yds$sc3_clusters
    ), aes( x = UMAP_1, y = UMAP_2, colour = sc3_clusters)
) + geom_point() + main + 
  xlab("UMAP1") +
  ylab("UMAP2"),
  
  
  ggplot(
  data.frame(
    reducedDim(sce.yds, "UMAP"),
    sample = sce.yds$sample,
    sc3_cluster = colData(sce.yds)[,paste0("sc3_", k, "_clusters")],
    sc3_outlier = colData(sce.yds)[,paste0("sc3_", k, "_log2_outlier_score")]
    ), aes( x = UMAP_1, y = UMAP_2, colour = sc3_cluster, size = sc3_outlier)
) + geom_point() + main + 
  xlab("UMAP1") +
  ylab("UMAP2")

)

nodes = data.frame(name = c(c("D.1", "D.2", "MG.1", "MG.2", "MG.3", "G.1", "G.2", "M.1"),paste("sc3_8", as.character(1:8)) ) )

df <- data.frame(table(paste0(sce.yds$sc3_clusters, "_", sce.yds$sc3_8_clusters)))

df$source <- as.character(sub("_.*", "", df$Var1))
df$source <- match(df$source, nodes$name[1:8]) - 1
df$target <- as.numeric(sub(".*_", "", df$Var1)) - 1 + 8

sankeyNetwork(Links = df, Nodes = nodes,
 Source = "source", Target = "target",
 Value = "Freq", NodeID = "name", fontSize= 12)

### DE gene analysis
sample.markers <- FindAllMarkers(seurat.yds, only.pos = TRUE, min.pct = .5, logfc.threshold = .5, test.use = "wilcox")
sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(seurat.yds, features = sample.markers$gene) + NoLegend()

FeaturePlot(seurat.yds, features = rbind(
  sample.markers$gene[head(which(sample.markers$cluster %in% "YeastDropseq_DMSO"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "YeastDropseq_Guanine"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "YeastDropseq_MPA"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "YeastDropseq_GM"), 3)]
  ),ncol = 3,reduction = "pca",) 

seurat.yds@active.ident <- setNames(as.factor(sce.yds$sc3_clusters), colnames(sce.yds))

sample.markers <- FindAllMarkers(seurat.yds, only.pos = TRUE, min.pct = .5, logfc.threshold = .5, test.use = "wilcox")
sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


DoHeatmap(seurat.yds, features = sample.markers$gene) + NoLegend()

FeaturePlot(seurat.yds, features = rbind(
  sample.markers$gene[head(which(sample.markers$cluster %in% "D.1"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "D.2"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "MG.1"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "MG.2"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "MG.3"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "G.1"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "G.2"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "M.1"), 3)]
  ),ncol = 3,reduction = "pca",)

seurat.yds@active.ident <- setNames(as.factor(sce.yds$sc3_4_clusters), colnames(sce.yds))

sample.markers <- FindAllMarkers(seurat.yds, only.pos = TRUE, min.pct = .5, logfc.threshold = .5, test.use = "wilcox")
sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(seurat.yds, features = sample.markers$gene) + NoLegend()

FeaturePlot(seurat.yds, features = rbind(
  sample.markers$gene[head(which(sample.markers$cluster %in% "1"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "2"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "3"), 3)],
  sample.markers$gene[head(which(sample.markers$cluster %in% "4"), 3)]
  ),ncol = 3,reduction = "pca",) 


### GO Term analysis of the 4 samples with `clusterProfiler`

object = readRDS("seurat_yeastdropseq_4samples.rds") ## Replace by Seurat Object created
df.gene.candidates <- read.delim2("Supplementary_DEgenes_4samples.tsv", header = T, sep = "\t")
df.gene.candidates
table(df.gene.candidates$sample)

minGeneAnnotation = 3
pvalueThres = .01
qvalueThres = .05

df.sample.go <- data.frame()

sample = setNames(c("DMSO", "Guanine", "MG", "MPA"), c("D", "G", "MG", "M"))
expression = setNames(c("negative", "positive"), c("DOWN", "UP"))


for( x in names(sample))
{
  for( y in names(expression) )
  {
    
    gene <- df.gene.candidates$gene[df.gene.candidates$sample %in% sample[x] & df.gene.candidates$expression %in% expression[y]]

    ego2 <- enrichGO(
      gene = paste0(gene, "_mRNA"),
      universe = paste0(rownames(object), "_mRNA"),
      OrgDb = org.Sc.sgd.db,
      keyType = "ENSEMBLTRANS",
      ont = "BP",
      pAdjustMethod = "BH",
      minGSSize = minGeneAnnotation,
      pvalueCutoff = pvalueThres,
      qvalueCutoff = qvalueThres
    )

    df.sample.go <- rbind(
      df.sample.go,
      data.frame(
        ego2@result[,-grep("Count", colnames(ego2@result))],
        sample = x,
        regulation = expression[y]
      )
    )
  }
}

df.sample.go[df.sample.go$p.adjust <= .05, ]


### GO Term analysis of the 8 SC3 clusters with `clusterProfiler`

object = readRDS("seurat_yeastdropseq_8clusters.rds") ## Replace by Seurat Object created
df.gene.candidates <- read.delim2("Supplementary_DEgenes_sc3_8clusters.tsv", header = T, sep = "\t")
df.gene.candidates
table(df.gene.candidates$sc3_cluster)

minGeneAnnotation = 3
pvalueThres = .01
qvalueThres = .05

df.sample.go <- data.frame()

sample = setNames(c("DMSO.1", "DMSO.2", "G.1", "G.2", "M", "MG.1", "MG.2", "MG.3"), c("DMSO.1", "DMSO.2", "G.1", "G.2", "M", "MG.1", "MG.2", "MG.3"))
expression = setNames(c("negative", "positive"), c("DOWN", "UP"))


for( x in names(sample))
{
  for( y in names(expression) )
  {
    
    gene <- df.gene.candidates$gene[df.gene.candidates$sc3_cluster %in% sample[x] & df.gene.candidates$expression %in% expression[y]]
    # print(gene)
    ego2 <- enrichGO(
      gene = paste0(gene, "_mRNA"),
      universe = paste0(rownames(object), "_mRNA"),
      OrgDb = org.Sc.sgd.db,
      keyType = "ENSEMBLTRANS",
      ont = "BP",
      pAdjustMethod = "BH",
      minGSSize = minGeneAnnotation,
      pvalueCutoff = pvalueThres,
      qvalueCutoff = qvalueThres
    )

    df.sample.go <- rbind(
      df.sample.go,
      data.frame(
        ego2@result[,-grep("Count", colnames(ego2@result))],
        sample = x,
        regulation = expression[y]
      )
    )
  }
}

df.sample.go[df.sample.go$p.adjust <= .05, ]
