# Import packages
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(magrittr)
library(WGCNA)
library(hdWGCNA)
library(corrplot)
library(Seurat)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(UCell)
library(fmsb)

# Using the cowplot theme for ggplot
theme_set(theme_cowplot())

# Set random seed for reproducibility
set.seed(12345)

# Optionally enable multi-threading
enableWGCNAThreads(nThreads = 8)

####################################
## Set up Seurat object for WGCNA ##
####################################

# Load the data (E13.5 Otic-associated cells)
CD_seu_harmony_sub2 <- SetupForWGCNA(
  CD_seu_harmony_sub2,
  wgcna_name = "CD"
)

#########################
## Construct metacells ##
#########################

# Construct metacells  in each group
CD_seu_harmony_sub2 <- MetacellsByGroups(
  seurat_obj = CD_seu_harmony_sub2,
  reduction = 'umap_harmony',
  k = 25,
  max_shared = 10
)

# Normalize metacell expression matrix
CD_seu_harmony_sub2 <- NormalizeMetacells(CD_seu_harmony_sub2)

####################################
## Co-expression network analysis ##
####################################

# Set up the expression matrix
CD_seu_harmony_sub2 <- SetDatExpr(
  CD_seu_harmony_sub2,
  assay = 'SCT' # using SCT assay
)

# Select soft-power threshold
CD_seu_harmony_sub2 <- TestSoftPowers(
  CD_seu_harmony_sub2,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
plot_list <- PlotSoftPowers(CD_seu_harmony_sub2)
# pdf("Soft_Power.pdf", width = 10, height = 9)
wrap_plots(plot_list, ncol=2)
# dev.off()
power_table <- GetPowerTable(CD_seu_harmony_sub2)

# construct co-expression network:
CD_seu_harmony_sub2 <- ConstructNetwork(
  CD_seu_harmony_sub2,
  tom_name = 'CD',
  overwrite_tom = TRUE
)
# pdf("Module_dendrogram.pdf", width = 10, height = 2.5)
PlotDendrogram(CD_seu_harmony_sub2, main='CD hdWGCNA Dendrogram')
# dev.off()

########################################
## Module Eigengenes and Connectivity ##
########################################

# Compute harmonized module eigengenes
CD_seu_harmony_sub2 <- ScaleData(CD_seu_harmony_sub2, features=VariableFeatures(CD_seu_harmony_sub2))
CD_seu_harmony_sub2 <- ModuleEigengenes(
  CD_seu_harmony_sub2
  #  group.by.vars="batch"
)

# harmonized module eigengenes:
hMEs <- GetMEs(CD_seu_harmony_sub2)

# module eigengenes:
MEs <- GetMEs(CD_seu_harmony_sub2, harmonized=FALSE)

# Compute module connectivity (kME)
CD_seu_harmony_sub2 <- ModuleConnectivity(
  CD_seu_harmony_sub2
)

# Rename the modules
CD_seu_harmony_sub2 <- ResetModuleNames(
  CD_seu_harmony_sub2,
  new_name = "CD-M"
)

# change module colors
new_colors <- list('CD-M1' = "#FF7F00",
                   'CD-M2' = "#4DAF4A",
                   'CD-M3' = "#E41A1C",
                   'CD-M4' = '#33CCCC',
                   'CD-M5' = '#003333',
                   'CD-M6' = '#FF9999',
                   'CD-M7' = '#333366',
                   'CD-M8' = '#6666CC',
                   'CD-M9' = '#CC6699',
                   'CD-M10' = '#660066')
CD_seu_harmony_sub2 <- ResetModuleColors(CD_seu_harmony_sub2, new_colors)

# Plot genes ranked by kME for each module
p <- PlotKMEs(CD_seu_harmony_sub2, ncol=5)
#pdf(file = "Modules_top_genes.pdf", width = 12, height = 4.5)
p
#dev.off()

# Get the module assignment table
modules <- GetModules(CD_seu_harmony_sub2) %>% subset(module != 'grey')
# Get hub genes
hub_df <- GetHubGenes(CD_seu_harmony_sub2, n_hubs = 10)

# Compute hub gene signature scores with UCell method
CD_seu_harmony_sub2 <- ModuleExprScore(
  CD_seu_harmony_sub2,
  n_genes = 25,
  method='UCell'
)

###################
## Visualization ##
###################
# 1. Feature Plot
# Make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  CD_seu_harmony_sub2,
  features='hMEs', 
  order=TRUE,
  reduction = "umap_harmony",
  plot_ratio = 1.3,
  point_size = 0.7
)

#pdf(file = "Modules_feature_plot.pdf", width = 20, height = 15)
wrap_plots(plot_list, ncol=5)
#dev.off()

# 2. Individual module network plot
ModuleNetworkPlot(
  CD_seu_harmony_sub2,
  outdir = 'ModuleNetworks'
)

# 3. UMAP (co-expression networks)
pdf("Module_UMAP_OV.pdf", width = 7, height = 7)
ModuleUMAPPlot(
  CD_seu_harmony_sub2,
  edge.alpha=1,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=0,# how many hub genes to plot per module?
  keep_grey_edges=FALSE,
  return_graph = FALSE
)
dev.off()

# 4. Radar Chart (using fmsb packages)
# Get modules
modules <- GetModules(CD_seu_harmony_sub2)

# Create df for radar chart
radar_df <- modules[,c(1,4:14)]
radar_df[, -1] <- apply(radar_df[, -1], 2, function(x) pmax(x, 0))
radar_df = radar_df[,2:12]
radar_df <- radar_df[c("Foxg1", "Gata3", "Otx2", "Bmp4", "Sox2", "Fgf10"), ] # CD marker genes

# Rename cols
col_name = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "Grey", "M8", "M9", "M10")
colnames(radar_df) <- col_name

# Manipulate df
max_min <- data.frame(
  M1 = c(1, 0),
  M2 = c(1, 0),
  M3 = c(1, 0),
  M4 = c(1, 0),
  M5 = c(1, 0),
  M6 = c(1, 0),
  M7 = c(1, 0),
  Grey = c(1, 0),
  M8 = c(1, 0),
  M9 = c(1, 0),
  M10 = c(1, 0)
)
radar_df <- rbind(max_min, radar_df)

# Assign colors to modules
colors_border=c('#FF9999', "#4DAF4A", "#FF7F00", "#E41A1C", '#33CCCC', '#6666CC')

colors_in <- c(
  rgb(255, 153, 153, maxColorValue = 255, alpha = 0.3*255),  # '#FF9999'
  rgb(77, 175, 74, maxColorValue = 255, alpha = 0.3*255),    # "#4DAF4A"
  rgb(255, 127, 0, maxColorValue = 255, alpha = 0.3*255),    # "#FF7F00"
  rgb(228, 26, 28, maxColorValue = 255, alpha = 0.3*255),    # "#E41A1C"
  rgb(51, 204, 204, maxColorValue = 255, alpha = 0.3*255),   # '#33CCCC'
  rgb(102, 102, 204, maxColorValue = 255, alpha = 0.3*255)   # '#6666CC'
)

# pdf("module_score_foxg1_gata3_radar.pdf", width = 8, height = 6)
radarchart(radar_df, axistype=1, 
           pcol=colors_border, pfcol=colors_in, plwd=3 , plty=1,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0, 100, 25), cglwd=0.6,
           vlcex=1.4
) 
legend(x=1.3, y=1, legend = rownames(radar_df3[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.2, pt.cex=3)
# dev.off()
