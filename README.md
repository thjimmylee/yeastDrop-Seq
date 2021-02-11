## yeastDrop-Seq Computational Pipeline

This repository contains scripts for data pre-processing and downstream analysis in the [yeastDrop-Seq paper](#).


Dataset is available at GEO viewer (Accession Number: [GSE165686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165686))

### Required softwares

#### Data preprocessing
[UMI-tools](https://github.com/CGATOxford/UMI-tools)
[samtools](https://github.com/samtools/samtools)
[Fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
[STAR](https://github.com/alexdobin/STAR)
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
[subread](http://subread.sourceforge.net/)

#### Downstream analysis
[SC3](https://bioconductor.org/packages/release/bioc/html/SC3.html)
[SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
[Seurat](https://cran.r-project.org/web/packages/Seurat/index.html)
[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
[scales](https://cran.r-project.org/web/packages/scales/index.html)
[ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)
[cowplot](https://cran.r-project.org/web/packages/cowplot/index.html)
[RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/)
[networkD3](https://cran.r-project.org/web/packages/networkD3/index.html)
[dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
[patchwork](https://cran.r-project.org/web/packages/patchwork/)
[reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
[clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
[org.Sc.sgd.db](https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html)