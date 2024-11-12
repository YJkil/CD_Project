# Molecular characterization of subdomain specification of the cochlear duct based on Foxg1 and Gata3
### Yongjin Gil 1, †, Jiho Ryu 1, †, Hayoung Yang 1, †, Yechan Ma 2 , Ki-Hoan Nam 3 , Sung-Wuk Jang 2,* and Sungbo Shim 
---
#### This repository contains the code used for the analysis and visualization in the paper.

<img width="1384" alt="Foxg1_Summary" src="https://github.com/user-attachments/assets/b2c508f7-633d-4b84-b48e-824f4d49535b">

#### [Schematic summary.]

<img width="1384" alt="Foxg1_Pipeline" src="https://github.com/user-attachments/assets/74d3671c-2db8-420e-bd54-ae6e9629ef6c">

#### [Schematic diagram of the data analysis process.]
---

+ ### 1_Seurat_Pipeline.R
  + #### This script includes code for data preprocessing (QC, filtering, etc.), integration, dimensionality reduction, clustering, marker gene identification, and cell type assignment, following the Seurat pipeline.

+ ### 2_Subset_Clustering.R
  + #### This script includes code for subsetting and clustering E13.5 otic-associated cells.

+ ### 3_hdWGCNA.R
  + #### This script includes code for co-expression network analysis, module eigengenes and connectivity and visualization using hdWGCNA.
 
+ ### 4_Foxg1_Gata3_Network.R
  + #### This script includes code for correlation analysis of Foxg1 and Gata3, network analysis using modules from the previous hdWGCNA, and visualization.

+ ### 5_Gene_Ontology.R
  + #### This script includes code for Gene Ontology analysis (Biological Process) and visualization of genes related to modules, as well as genes highly correlated with Foxg1 and Gata3.

---
## R Packages Used

Below is a list of R packages used in this analysis:

### Core Packages
- **clusterProfiler**: v4.13.4
- **AnnotationHub**: v3.13.3
- **BiocFileCache**: v2.13.2
- **dbplyr**: v2.5.0
- **org.Mm.eg.db**: v3.20.0
- **AnnotationDbi**: v1.67.0
- **IRanges**: v2.38.1
- **S4Vectors**: v0.42.1
- **Biobase**: v2.64.0
- **BiocGenerics**: v0.51.3
- **GOSemSim**: v2.31.2
- **ggraph**: v2.2.1
- **fmsb**: v0.7.6
- **pheatmap**: v1.0.12
- **reshape2**: v1.4.4
- **readxl**: v1.4.3
- **UCell**: v2.10.1
- **corrplot**: v0.95
- **hdWGCNA**: v0.3.02
- **igraph**: v2.0.3
- **ggrepel**: v0.9.6
- **WGCNA**: v1.73
- **fastcluster**: v1.2.6
- **dynamicTreeCut**: v1.63-1
- **magrittr**: v2.0.3
- **patchwork**: v1.3.0
- **cowplot**: v1.1.3
- **future**: v1.34.0
- **glmGamPoi**: v1.18.0
- **Matrix**: v1.7-0
- **harmony**: v1.2.1
- **Rcpp**: v1.0.13
- **gridExtra**: v2.3
- **lubridate**: v1.9.3
- **forcats**: v1.0.0
- **stringr**: v1.5.1
- **purrr**: v1.0.2
- **readr**: v2.1.5
- **tidyr**: v1.3.1
- **tibble**: v3.2.1
- **ggplot2**: v3.5.1
- **tidyverse**: v2.0.0
- **dplyr**: v1.1.4
- **Seurat**: v5.1.0
- **SeuratObject**: v5.0.2
- **sp**: v2.1-4

### Other Packages
- **fs**: v1.6.4
- **matrixStats**: v1.4.1
- **spatstat.sparse**: v3.1-0
- **enrichplot**: v1.26.1
- **httr**: v1.4.7
- **RColorBrewer**: v1.1-3
- **doParallel**: v1.0.17
- **tools**: v4.4.1
- **sctransform**: v0.4.1
- **backports**: v1.5.0
- **utf8**: v1.2.4
- **R6**: v2.5.1
- **lazyeval**: v0.2.2
- **uwot**: v0.2.2
- **withr**: v3.0.2
- **preprocessCore**: v1.68.0
- **progressr**: v0.15.0
- **cli**: v3.6.3
- **spatstat.explore**: v3.3-2
- **fastDummies**: v1.7.4
- **labeling**: v0.4.3
- **spatstat.data**: v3.1-2
- **proxy**: v0.4-27
- **ggridges**: v0.5.6
- **pbapply**: v1.7-2
- **yulab.utils**: v0.1.7
- **gson**: v0.1.0
- **foreign**: v0.8-87
- **DOSE**: v4.0.0
- **R.utils**: v2.12.3
- **parallelly**: v1.38.0
- **rstudioapi**: v0.17.1
- **impute**: v1.80.0
- **RSQLite**: v2.3.7
- **FNN**: v1.1.4.1
- **gridGraphics**: v0.5-1
- **generics**: v0.1.3
- **ica**: v1.0-3
- **spatstat.random**: v3.3-2
- **GO.db**: v3.20.0
- **fansi**: v1.0.6
- **abind**: v1.4-8
- **R.methodsS3**: v1.8.2
- **lifecycle**: v1.0.4
- **yaml**: v2.3.10
- **SummarizedExperiment**: v1.36.0
- **qvalue**: v2.38.0
- **SparseArray**: v1.6.0
- **Rtsne**: v0.17
- **blob**: v1.2.4
- **promises**: v1.3.0
- **crayon**: v1.5.3
- **ggtangle**: v0.0.4
- **miniUI**: v0.1.1.1
- **lattice**: v0.22-6
- **KEGGREST**: v1.46.0
- **pillar**: v1.9.0
- **knitr**: v1.48
- **fgsea**: v1.31.6
- **GenomicRanges**: v1.58.0
- **future.apply**: v1.11.3
- **codetools**: v0.2-20
- **fastmatch**: v1.1-4
- **leiden**: v0.4.3.1
- **glue**: v1.7.0
- **ggfun**: v0.1.7
- **spatstat.univar**: v3.0-1
- **data.table**: v1.16.2
- **treeio**: v1.30.0
- **vctrs**: v0.6.5
- **png**: v0.1-8
- **spam**: v2.11-0
- **cellranger**: v1.1.0
- **gtable**: v0.3.6
- **cachem**: v1.1.0
- **xfun**: v0.49
- **S4Arrays**: v1.6.0
- **mime**: v0.12
- **tidygraph**: v1.3.1
- **survival**: v3.7-0
- **SingleCellExperiment**: v1.28.0
- **iterators**: v1.0.14
- **fitdistrplus**: v1.2-1
- **ROCR**: v1.0-11
- **nlme**: v3.1-165
- **ggtree**: v3.14.0
- **bit64**: v4.5.2
- **filelock**: v1.0.3
- **RcppAnnoy**: v0.0.22
- **GenomeInfoDb**: v1.42.0
- **irlba**: v2.3.5.1
- **KernSmooth**: v2.23-24
- **rpart**: v4.1.23
- **colorspace**: v2.1-1
- **DBI**: v1.2.3
- **Hmisc**: v5.2-0
- **nnet**: v7.3-19
- **tidyselect**: v1.2.1
- **curl**: v5.2.3
- **bit**: v4.5.0
- **compiler**: v4.4.1
- **htmlTable**: v2.4.3
- **BiocNeighbors**: v2.0.0
- **DelayedArray
