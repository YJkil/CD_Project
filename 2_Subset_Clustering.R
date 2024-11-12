####################################
## Subset 1 (Otic-associated cell)##
####################################
CD_seu_harmony_sub1 = subset(CD_seu_harmony, subset = SCT_snn_res.0.03 %in% c(2))
CD_seu_harmony_sub1@meta.data = CD_seu_harmony_sub1@meta.data[,c(1:8,10:11)]

# SCTransform
CD_seu_harmony_sub1 <- CD_seu_harmony_sub1 %>%
  SCTransform(vars.to.regress = c("mitoRatio"))
DefaultAssay(CD_seu_harmony_sub1) <- "SCT"

# Find variable features for integrate
integ_features <- SelectIntegrationFeatures5(CD_seu_harmony_sub1, assay = "SCT", nfeatures = 3000) 
VariableFeatures(CD_seu_harmony_sub1) <- integ_features

# Run PCA
CD_seu_harmony_sub1 <- RunPCA(CD_seu_harmony_sub1, assay = "SCT", npcs = 50)
# CD_seu_harmony_sub1 <- RunUMAP(object = CD_seu_harmony_sub1, dims = 1:30, reduction = 'pca')
# CD_seu_harmony_sub1 <- RunTSNE(CD_seu_harmony_sub1, reduction = "harmony", dims = 1:30, reduction.name = "tsne_harmony")

# Integration (Harmony)
CD_seu_harmony_sub1 <- IntegrateLayers(
  object = CD_seu_harmony_sub1,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  assay = "SCT",
  features = integ_features,
  new.reduction = "harmony"
)

# Clustering
CD_seu_harmony_sub1 <- FindNeighbors(CD_seu_harmony_sub1, reduction = "harmony", dims = 1:30)
CD_seu_harmony_sub1 <- FindClusters(CD_seu_harmony_sub1, resolution = 0.03)

# Run UMAP
CD_seu_harmony_sub1 <- RunUMAP(CD_seu_harmony_sub1, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")

# DimPlot(CD_seu_harmony_sub1, reduction = "umap_harmony", group.by = c("seurat_clusters"), pt.size = 1.5)

# Save to RDS file
# saveRDS(CD_seu_harmony_sub1, file = "CD_seu_harmony_sub1.rds")

##########################################
## Subset 2 (E13.5 Otic-associated cell)##
##########################################
CD_seu_harmony_sub2 = subset(CD_seu_harmony_sub1, subset = SCT_snn_res.0.03 %in% c(2, 3))
CD_seu_harmony_sub2 = subset(CD_seu_harmony_sub2, subset = day %in% c("E13.5"))

# SCTransform
CD_seu_harmony_sub2 <- CD_seu_harmony_sub2 %>%
  SCTransform(vars.to.regress = c("mitoRatio"))
DefaultAssay(CD_seu_harmony_sub2) <- "SCT"

# Find most variable features across samples to integrate
integ_features <- SelectIntegrationFeatures5(CD_seu_harmony_sub2, assay = "SCT", nfeatures = 3000) 
VariableFeatures(CD_seu_harmony_sub2) <- integ_features

# Run PCA
CD_seu_harmony_sub2 <- RunPCA(CD_seu_harmony_sub2, assay = "SCT", npcs = 50)
# CD_seu_harmony_sub2 <- RunUMAP(object = CD_seu_harmony_sub2, dims = 1:30, reduction = 'pca')

# Integration (Harmony)
CD_seu_harmony_sub2 <- IntegrateLayers(
  object = CD_seu_harmony_sub2,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  assay = "SCT",
  features = integ_features,
  new.reduction = "harmony"
)

# Clustering
CD_seu_harmony_sub2 <- FindNeighbors(CD_seu_harmony_sub2, reduction = "harmony", dims = 1:30)
CD_seu_harmony_sub2 <- FindClusters(CD_seu_harmony_sub2, resolution = 0.03)
# CD_seu_harmony_sub2 <- FindClusters(CD_seu_harmony_sub2, resolution = 0.1)
# CD_seu_harmony_sub2 <- FindClusters(CD_seu_harmony_sub2, resolution = 0.5)
# CD_seu_harmony_sub2 <- FindClusters(CD_seu_harmony_sub2, resolution = 1)

# Run UMAP
CD_seu_harmony_sub2 <- RunUMAP(CD_seu_harmony_sub2, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
# DimPlot(CD_seu_harmony_sub2, reduction = "umap_harmony", group.by = c("SCT_snn_res.0.03"), pt.size = 1.5)

# Save to RDS file
# saveRDS(CD_seu_harmony_sub2, "CD_seu_harmony_sub2.rds")
