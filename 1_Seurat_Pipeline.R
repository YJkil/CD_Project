## Import packages
library(Seurat)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(harmony)
library(Matrix)
library(tidyr)
library(glmGamPoi)
library(future)

##########################################################################################################
#####################################################################################################
## Cao et al., 2019 (The single-cell transcriptional landscape of mammalian organogenesis, Nature) ##
#####################################################################################################

## Change working directory
setwd("./Cao/sub_trajectory_summary")
getwd()

## Read CDS object (Auditory Epithelial Trajectory)
Cao_auditory_raw = readRDS("Auditory_epithelial_trajectory_cds.RDS")

########################################
## Change CDS obejct -> Seurat object ##
########################################
# Extract Expression matrix
exp_mat <- exprs(Cao_auditory_raw)

# Extract Cell metadata
cell_meta <- pData(Cao_auditory_raw)

# Extract Feature metadata
feature_meta <- fData(Cao_auditory_raw)

# Change Ensemble name to short gene name
rownames(exp_mat) <- feature_meta$gene_short_name
rownames(exp_mat) <- make.unique(rownames(exp_mat)) # remove duplicate rownames

# Create seurat object
Cao_seu <- CreateSeuratObject(counts = exp_mat, meta.data = cell_meta, project = "Cao")
# head(Cao_seu@assays$RNA$counts) # check counts matrix
# head(Cao_seu@meta.data) # check metadata

# Extract metadata
Cao_seu@meta.data = Cao_seu@meta.data[,c(1:3,13)]
Cao_seu@meta.data$ref = "Cao_2019" # Add ref to metadata

##########################################################################################################
##########################################################################################################################################
## Sun et al., 2022 (Single-cell transcriptomic landscapes of the otic neuronal lineage at multiple early embryonic ages, Cell Reports) ##
##########################################################################################################################################

## Change working directory
setwd("./dataset")

##########################
## Create Seurat object ##
##########################
# Function to read data and create Seurat object
create_seurat_object <- function(file_path) {
  data <- read.csv(gzfile(file_path), header = TRUE)
  
  genes <- data[, 1]
  rownames(data) <- genes
  data <- data[, -1]
  
  data_sparse <- as(as.matrix(data), "dgCMatrix")
  
  CreateSeuratObject(counts = data_sparse, assay = "RNA", min.cells = 0, min.features = 0, project = "Sun")
}

# File paths
file_paths <- c("GSM5401072_E9_5_raw_count.csv.gz", 
                "GSM5401073_E11_5_1_raw_count.csv.gz", 
                "GSM5401074_E11_5_2_raw_count.csv.gz", 
                "GSM5401075_E11_5_3_raw_count.csv.gz", 
                "GSM5401076_E13_5_1_raw_count.csv.gz", 
                "GSM5401077_E13_5_2_raw_count.csv.gz", 
                "GSM5401078_E13_5_3_raw_count.csv.gz", 
                "GSM5401079_E13_5_4_raw_count.csv.gz")

# Create Seurat objects
Sun_seus <- lapply(file_paths, create_seurat_object)

# Separate each Seurat object
Sun_seu_9_1 <- Sun_seus[[1]]
Sun_seu_11_1 <- Sun_seus[[2]]
Sun_seu_11_2 <- Sun_seus[[3]]
Sun_seu_11_3 <- Sun_seus[[4]]
Sun_seu_13_1 <- Sun_seus[[5]]
Sun_seu_13_2 <- Sun_seus[[6]]
Sun_seu_13_3 <- Sun_seus[[7]]
Sun_seu_13_4 <- Sun_seus[[8]]

## Add metadata
# Stage
Sun_seu_9_1@meta.data$day = "9.5"
Sun_seu_11_1@meta.data$day = "11.5"
Sun_seu_11_2@meta.data$day = "11.5"
Sun_seu_11_3@meta.data$day = "11.5"
Sun_seu_13_1@meta.data$day = "13.5"
Sun_seu_13_2@meta.data$day = "13.5"
Sun_seu_13_3@meta.data$day = "13.5"
Sun_seu_13_4@meta.data$day = "13.5"

# Reference
Sun_seu_9_1@meta.data$ref = "Sun_2022"
Sun_seu_11_1@meta.data$ref = "Sun_2022"
Sun_seu_11_2@meta.data$ref = "Sun_2022"
Sun_seu_11_3@meta.data$ref = "Sun_2022"
Sun_seu_13_1@meta.data$ref = "Sun_2022"
Sun_seu_13_2@meta.data$ref = "Sun_2022"
Sun_seu_13_3@meta.data$ref = "Sun_2022"
Sun_seu_13_4@meta.data$ref = "Sun_2022"

##########################################################################################################
#############
## Merging ##
#############
# Merge all Seurat object
CD_seu <- merge(x = Cao_seu, 
                y = c(Sun_seu_9_1, Sun_seu_11_1, Sun_seu_11_2, Sun_seu_11_3,
                      Sun_seu_13_1, Sun_seu_13_2, Sun_seu_13_3, Sun_seu_13_4),
                merge.data = TRUE)

#####################
## Quality Control ##
#####################
# Calculate Mitochondria percent
CD_seu[["mitoRatio"]] <- PercentageFeatureSet(CD_seu, pattern = "^mt-")
table(CD_seu$mitoRatio)

# Filtering
CD_seu_filtered <- subset(CD_seu, subset = nCount_RNA > 800 & nFeature_RNA > 200 & mitoRatio < 10)
CD_seu_filtered

####################
## Pre-processing ##
####################
# Normalization & Features selection & SCT normalization
options(future.globals.maxSize = 2 * 1024^3)  # 2 GiB
CD_seu_filtered <- CD_seu_filtered %>%
  SCTransform(vars.to.regress = c("mitoRatio"))
DefaultAssay(CD_seu_filtered) <- "SCT"

# Run PCA & UMAP
CD_seu_filtered <- RunPCA(CD_seu_filtered, assay = "SCT", npcs = 50)
CD_seu_filtered <- RunUMAP(CD_seu_filtered, dims = 1:30, reduction = 'pca')

# Dimplot(UMAP) before integration
# DimPlot(CD_seu_filtered, reduction = "umap", group.by = "ref")

#################
## Integration ##
#################
# Find variable features to integrate
inte_features <- SelectIntegrationFeatures5(CD_seu_filtered, assay = "SCT", nfeatures = 3000) 

# Set variable features to integrate
VariableFeatures(CD_seu_filtered) <- inte_features

# Run Harmony
CD_seu_harmony <- IntegrateLayers(
  object = CD_seu_filtered,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  assay = "SCT",
  features = inte_features,
  new.reduction = "harmony",
)

# Clustering
CD_seu_harmony <- FindNeighbors(CD_seu_harmony, reduction = "harmony", dims = 1:30)
CD_seu_harmony <- FindClusters(CD_seu_harmony, resolution = 0.03)


# Run UMAP
CD_seu_harmony <- RunUMAP(CD_seu_harmony, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
# CD_seu_harmony <- RunTSNE(CD_seu_harmony, reduction = "harmony", dims = 1:30, reduction.name = "tsne_harmony") # tSNE

# Dimplot(UMAP) after integration
# DimPlot(CD_seu_harmony, reduction = "umap_harmony", group.by = "seurat_clusters")
# DimPlot(CD_seu_harmony, reduction = "umap_harmony", group.by = "day")
# DimPlot(CD_seu_harmony, reduction = "umap_harmony", group.by = "ref")

# Rename metadata (day)
CD_seu_harmony@meta.data$day <- paste0("E", CD_seu_harmony@meta.data$day)
CD_seu_harmony@meta.data$day <- factor(CD_seu_harmony@meta.data$day, levels = c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5"))

# Find Clusters marker genes
CD_seu_harmony <- PrepSCTFindMarkers(CD_seu_harmony)
CD.markers <- FindAllMarkers(CD_seu_harmony, only.pos = TRUE, assay = "SCT")

# Extract Top 5 marker genes
CD.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

top5 <- top5 %>%
  distinct(gene, .keep_all = TRUE)

# Assign Cell types
new.cluster.ids <- c("C0 - Osteoblasts", 
                     "C1 - Radial glial progenitor cells", 
                     "C2 - Otic epithelial cells", 
                     "C3 - Neuronal cells", 
                     "C4 - Inner ear glial cells",
                     "C5 - Blood cells",
                     "C6 - Endothelial-haematopoietic cells", 
                     "C7 - Macrophages")
names(new.cluster.ids) <- levels(CD_seu_harmony)
CD_seu_harmony <- RenameIdents(CD_seu_harmony, new.cluster.ids)
CD_seu_harmony$Cell_type = Idents(CD_seu_harmony)

# Save to RDS file
# saveRDS(CD_seu_harmony, "./CD_seu_harmony.rds")