# **Molecular characterization of subdomain specification of the cochlear duct based on Foxg1 and Gata3**
### Yongjin Gil 1, †, Jiho Ryu 1, †, Hayoung Yang 1, †, Yechan Ma 2 , Ki-Hoan Nam 3 , Sung-Wuk Jang 2,* and Sungbo Shim 
---
#### This repository contains the code used for the analysis and visualization in the paper.

<img width="1384" alt="Foxg1_Summary2" src="https://github.com/user-attachments/assets/ed58bf65-d798-4d33-9265-7652a9c28818">

#### [Schematic summary.]

<img width="1384" alt="Foxg1_Pipeline2" src="https://github.com/user-attachments/assets/35bccedb-450c-4809-b323-df24623161f4">


#### [Schematic diagram of the data analysis process.]
---

+ ### **1_Seurat_Pipeline.R**
+ 
  + #### This script includes code for data preprocessing (QC, filtering, etc.), integration, dimensionality reduction, clustering, marker gene identification, and cell type assignment, following the Seurat pipeline.

+ ### **2_Subset_Clustering.R**
  + #### This script includes code for subsetting and clustering E13.5 otic-associated cells.

+ ### **3_hdWGCNA.R**
  + #### This script includes code for co-expression network analysis, module eigengenes and connectivity, and visualization using hdWGCNA.
 
+ ### **4_Foxg1_Gata3_Network.R**
  + #### This script includes code for correlation analysis of Foxg1 and Gata3, network analysis using modules from the previous hdWGCNA, and visualization.

+ ### **5_Gene_Ontology.R**
  + #### This script includes code for Gene Ontology analysis (Biological Process) and visualization of genes related to modules, as well as genes highly correlated with Foxg1 and Gata3.

+ ### **6_HL_Association.R**
  + #### This script includes code for enrichment analysis to identify overlapping genes between genes in each module and those in the Hearing Loss database (Clingen), pairwise correlation analysis among HL-related genes, and visualization.

+ ### **7_Visualization.R**
  + #### This script includes code for UMAP visualization by cell type, stage, and reference, dot plots displaying marker gene expression by cell type, and feature plots showing the expression of otic-associated markers.
 
+ ### **Hearing_Loss_Genes_Clingen.csv**
  + #### This file includes hearing loss-related genes from the Clingen database.
  + #### DiStefano, M.T.; Hughes, M.Y.; Patel, M.J.; Wilcox, E.H.; Oza, A.M. Expert interpretation of genes and variants in hereditary hearing loss. Medizinische Genetik 2020, 32, 109-115, doi:10.1515/medgen-2020-2018.
---
## Session Info

- **R version**: 4.4.1 (2024-06-14 ucrt)  
- **Platform**: Windows 11 x64 (build 22631)  
- **Time Zone**: Asia/Seoul

### Loaded Packages

- **Seurat**: 5.1.0
- **hdWGCNA**: 0.3.02
- **clusterProfiler**: 4.13.4
- **AnnotationHub**: 3.13.3
- **org.Mm.eg.db**: 3.20.0
- **ggraph**: 2.2.1
- **ggplot2**: 3.5.1
- **pheatmap**: 1.0.12
- **dplyr**: 1.1.4
- **tidyverse**: 2.0.0
- **igraph**: 2.0.3
- **Matrix**: 1.7-0

### Other Notable Libraries

- **AnnotationDbi**, **Biobase**, **BiocGenerics**
- **IRanges**, **S4Vectors**, **BiocFileCache**
- **magrittr**, **patchwork**, **cowplot**
- **UCell**, **corrplot**, **WGCNA**, **dynamicTreeCut**
- **lubridate**, **ggrepel**, **gridExtra**, **readxl**
- **fastcluster**, **future**, **glmGamPoi**
