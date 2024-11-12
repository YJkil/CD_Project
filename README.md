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
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=Korean_Korea.utf8  LC_CTYPE=Korean_Korea.utf8    LC_MONETARY=Korean_Korea.utf8
[4] LC_NUMERIC=C                  LC_TIME=Korean_Korea.utf8    

time zone: Asia/Seoul
tzcode source: internal

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] clusterProfiler_4.13.4 AnnotationHub_3.13.3   BiocFileCache_2.13.2   dbplyr_2.5.0          
 [5] org.Mm.eg.db_3.20.0    AnnotationDbi_1.67.0   IRanges_2.38.1         S4Vectors_0.42.1      
 [9] Biobase_2.64.0         BiocGenerics_0.51.3    GOSemSim_2.31.2        ggraph_2.2.1          
[13] fmsb_0.7.6             pheatmap_1.0.12        reshape2_1.4.4         readxl_1.4.3          
[17] UCell_2.10.1           corrplot_0.95          hdWGCNA_0.3.02         igraph_2.0.3          
[21] ggrepel_0.9.6          WGCNA_1.73             fastcluster_1.2.6      dynamicTreeCut_1.63-1 
[25] magrittr_2.0.3         patchwork_1.3.0        cowplot_1.1.3          future_1.34.0         
[29] glmGamPoi_1.18.0       Matrix_1.7-0           harmony_1.2.1          Rcpp_1.0.13           
[33] gridExtra_2.3          lubridate_1.9.3        forcats_1.0.0          stringr_1.5.1         
[37] purrr_1.0.2            readr_2.1.5            tidyr_1.3.1            tibble_3.2.1          
[41] ggplot2_3.5.1          tidyverse_2.0.0        dplyr_1.1.4            Seurat_5.1.0          
[45] SeuratObject_5.0.2     sp_2.1-4              

loaded via a namespace (and not attached):
  [1] fs_1.6.4                    matrixStats_1.4.1           spatstat.sparse_3.1-0       enrichplot_1.26.1          
  [5] httr_1.4.7                  RColorBrewer_1.1-3          doParallel_1.0.17           tools_4.4.1                
  [9] sctransform_0.4.1           backports_1.5.0             utf8_1.2.4                  R6_2.5.1                   
 [13] lazyeval_0.2.2              uwot_0.2.2                  withr_3.0.2                 preprocessCore_1.68.0      
 [17] progressr_0.15.0            cli_3.6.3                   spatstat.explore_3.3-2      fastDummies_1.7.4          
 [21] labeling_0.4.3              spatstat.data_3.1-2         proxy_0.4-27                ggridges_0.5.6             
 [25] pbapply_1.7-2               yulab.utils_0.1.7           gson_0.1.0                  foreign_0.8-87             
 [29] DOSE_4.0.0                  R.utils_2.12.3              parallelly_1.38.0           rstudioapi_0.17.1          
 [33] impute_1.80.0               RSQLite_2.3.7               FNN_1.1.4.1                 gridGraphics_0.5-1         
 [37] generics_0.1.3              ica_1.0-3                   spatstat.random_3.3-2       GO.db_3.20.0               
 [41] fansi_1.0.6                 abind_1.4-8                 R.methodsS3_1.8.2           lifecycle_1.0.4            
 [45] yaml_2.3.10                 SummarizedExperiment_1.36.0 qvalue_2.38.0               SparseArray_1.6.0          
 [49] Rtsne_0.17                  blob_1.2.4                  promises_1.3.0              crayon_1.5.3               
 [53] ggtangle_0.0.4              miniUI_0.1.1.1              lattice_0.22-6              KEGGREST_1.46.0            
 [57] pillar_1.9.0                knitr_1.48                  fgsea_1.31.6                GenomicRanges_1.58.0       
 [61] future.apply_1.11.3         codetools_0.2-20            fastmatch_1.1-4             leiden_0.4.3.1             
 [65] glue_1.7.0                  ggfun_0.1.7                 spatstat.univar_3.0-1       data.table_1.16.2          
 [69] treeio_1.30.0               vctrs_0.6.5                 png_0.1-8                   spam_2.11-0                
 [73] cellranger_1.1.0            gtable_0.3.6                cachem_1.1.0                xfun_0.49                  
 [77] S4Arrays_1.6.0              mime_0.12                   tidygraph_1.3.1             survival_3.7-0             
 [81] SingleCellExperiment_1.28.0 iterators_1.0.14            fitdistrplus_1.2-1          ROCR_1.0-11                
 [85] nlme_3.1-165                ggtree_3.14.0               bit64_4.5.2                 filelock_1.0.3             
 [89] RcppAnnoy_0.0.22            GenomeInfoDb_1.42.0         irlba_2.3.5.1               KernSmooth_2.23-24         
 [93] rpart_4.1.23                colorspace_2.1-1            DBI_1.2.3                   Hmisc_5.2-0                
 [97] nnet_7.3-19                 tidyselect_1.2.1            curl_5.2.3                  bit_4.5.0                  
[101] compiler_4.4.1              htmlTable_2.4.3             BiocNeighbors_2.0.0         DelayedArray_0.32.0        
[105] plotly_4.10.4               checkmate_2.3.2             scales_1.3.0                lmtest_0.9-40              
[109] rappdirs_0.3.3              digest_0.6.36               goftest_1.2-3               spatstat.utils_3.1-0       
[113] rmarkdown_2.29              XVector_0.44.0              RhpcBLASctl_0.23-42         htmltools_0.5.8.1          
[117] pkgconfig_2.0.3             base64enc_0.1-3             sparseMatrixStats_1.18.0    MatrixGenerics_1.18.0      
[121] fastmap_1.2.0               rlang_1.1.4                 htmlwidgets_1.6.4           UCSC.utils_1.2.0           
[125] shiny_1.9.1                 DelayedMatrixStats_1.28.0   farver_2.1.2                zoo_1.8-12                 
[129] jsonlite_1.8.9              BiocParallel_1.39.0         R.oo_1.27.0                 ggplotify_0.1.2            
[133] Formula_1.2-5               GenomeInfoDbData_1.2.13     dotCall64_1.2               munsell_0.5.1              
[137] ape_5.8                     viridis_0.6.5               reticulate_1.39.0           stringi_1.8.4              
[141] zlibbioc_1.50.0             MASS_7.3-61                 plyr_1.8.9                  parallel_4.4.1             
[145] listenv_0.9.1               deldir_2.0-4                Biostrings_2.73.2           graphlayouts_1.2.0         
[149] splines_4.4.1               tensor_1.5                  hms_1.1.3                   spatstat.geom_3.3-3        
[153] RcppHNSW_0.6.0              BiocVersion_3.20.0          evaluate_1.0.1              tester_0.2.0               
[157] BiocManager_1.30.25         tzdb_0.4.0                  foreach_1.5.2               tweenr_2.0.3               
[161] httpuv_1.6.15               RANN_2.6.2                  polyclip_1.10-7             scattermore_1.2            
[165] ggforce_0.4.2               xtable_1.8-4                tidytree_0.4.6              RSpectra_0.16-2            
[169] later_1.3.2                 viridisLite_0.4.2           aplot_0.2.3                 memoise_2.0.1              
[173] cluster_2.1.6               timechange_0.3.0            globals_0.16.3  
