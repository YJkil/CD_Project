# Import packages
library(cowplot)

#### Plot 1 - UMAP (Reference)
# pdf("1_UMAP_Ref_umap.pdf", width = 7, height = 7)
DimPlot(CD_seu_harmony, reduction = "umap_harmony", group.by = c("ref"), pt.size = 0.000001) + 
  coord_fixed(ratio = 1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30)
  ) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  ggtitle("Reference")
# dev.off()

##### Plot 2 - UMAP (Cell Type)
# pdf("2_UMAP_Cell_Type_umap.pdf", width = 7, height = 7)
DimPlot(CD_seu_harmony, reduction = "umap_harmony", group.by = c("Cell_type"), pt.size = 0.000001) + 
  coord_fixed(ratio = 1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10)
  ) +
  scale_color_manual(values = c(
    "#E41A1C",  # Red
    "#FF7F00",  # Orange
    "#377EB8",  # Blue
    "#4DAF4A",  # Green
    "#984EA3",  # Purple
    "#A65628",  # Brown
    "#F781BF",  # Pink
    "#66C2A5"  # Cyan
  )) +
  ggtitle("Cell Type")
# dev.off()

##### Plot 3 - UMAP (Stage)
# pdf("3_UMAP_Stage_umap.pdf", width = 7, height = 7)
DimPlot(CD_seu_harmony, reduction = "umap_harmony", group.by = c("day"), pt.size = 0.000001) + 
  coord_fixed(ratio = 1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30)
  ) +
  scale_color_manual(values = c("#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A", "#984EA3")) +
  ggtitle("Stage")
# dev.off()

##### Plot 4 - Dot Plot (Marker gene Expression)
CD_seu_harmony@meta.data$seurat_clusters <- paste0("C", CD_seu_harmony@meta.data$seurat_clusters)

# pdf("4_DotPlot_umap.pdf", width = 5, height = 7)
DotPlot(CD_seu_harmony, assay = "SCT",
        features = top5$gene,
        group.by = "seurat_clusters") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.2),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
#    aspect.ratio = 4
  ) +
  scale_color_gradient(low = "white", high = "black") + coord_flip()
# dev.off()

##### Plot 5 - Feature Plot
p1 = FeaturePlot(CD_seu_harmony, features = "Col1a2",
            reduction = "umap_harmony", pt.size = 0.000001) + 
  coord_fixed(ratio = 1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30)
  ) +
  scale_color_gradient(low = "lightgrey", high = "#E41A1C") + 
  NoLegend()


p2 = FeaturePlot(CD_seu_harmony, features = c("Sox2"),
            reduction = "umap_harmony", pt.size = 0.000001) + 
  coord_fixed(ratio = 1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30)
  ) +
  scale_color_gradient(low = "lightgrey", high = "#FF7F00") + 
  NoLegend()

p3 = FeaturePlot(CD_seu_harmony, features = c("Epcam"),
                 reduction = "umap_harmony", pt.size = 0.000001) + 
  coord_fixed(ratio = 1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30)
  ) +
  scale_color_gradient(low = "lightgrey", high = "#377EB8") + 
  NoLegend()

p4 = FeaturePlot(CD_seu_harmony, features = c("Dcx"),
                 reduction = "umap_harmony", pt.size = 0.000001) + 
  coord_fixed(ratio = 1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30)
  ) +
  scale_color_gradient(low = "lightgrey", high = "#4DAF4A") + 
  NoLegend()

p5 = FeaturePlot(CD_seu_harmony, features = c("Plp1"),
                 reduction = "umap_harmony", pt.size = 0.000001) + 
  coord_fixed(ratio = 1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30)
  ) +
  scale_color_gradient(low = "lightgrey", high = "#984EA3") + 
  NoLegend()

p6 = FeaturePlot(CD_seu_harmony, features = c("Hba-a1"),
                 reduction = "umap_harmony", pt.size = 0.000001) + 
  coord_fixed(ratio = 1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30)
  ) +
  scale_color_gradient(low = "lightgrey", high = "#A65628") + 
  NoLegend()

p7 = FeaturePlot(CD_seu_harmony, features = c("Ramp2"),
                 reduction = "umap_harmony", pt.size = 0.000001) + 
  coord_fixed(ratio = 1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30)
  ) +
  scale_color_gradient(low = "lightgrey", high = "#F781BF") + 
  NoLegend()

p8 = FeaturePlot(CD_seu_harmony, features = c("Fcer1g"),
                 reduction = "umap_harmony", pt.size = 0.000001) + 
  coord_fixed(ratio = 1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30)
  ) +
  scale_color_gradient(low = "lightgrey", high = "#66C2A5") + 
  NoLegend()

# pdf("5_Feature_marker_vertical_umap.pdf", width = 28, height = 14)
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, nrow = 2)
# dev.off()

##### Plot 6 - Feature Plot (blend Foxg1 & Gata3)
# pdf("15_Feature_UMAP_Foxg1_Gata3.pdf", width = 12, height = 3)
FeaturePlot(CD_seu_harmony_sub2, features = c("Foxg1", "Gata3"), blend = T, reduction = "umap_harmony", pt.size = 0.005, cols = c("lightgrey", "#E41A1C", "#4DAF4A"), blend.threshold = 1)
# dev.off()
