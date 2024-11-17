# Import packages
library(readxl)
library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(Seurat)
library(pheatmap)

#######################################
## Disease genes enrichment analysis ##
#######################################

# Load Disease dataset (Clingen)
HL_db <- read_csv("Hearing_Loss_Genes_Clingen.csv")

# Change gene name (human -> mouse)
HL_db$Gene <- str_to_title(tolower(HL_db$Gene))

# Modify data.frame
HL_modules <- modules %>% filter(gene_name %in% HL_db$Gene) # for module
HL_modules$module <- factor(HL_modules$module,
                            levels = c("CD-M1", "CD-M2", "CD-M3", "CD-M4", "CD-M5", "CD-M6", "CD-M7", "CD-M8", "CD-M9", "CD-M10"))

# Sort df
HL_modules_sorted <- HL_modules[order(HL_modules$module), ]
HL_modules_sorted

# Select kME columns
kME_columns <- c("kME_CD-M1", "kME_CD-M2", "kME_CD-M3", "kME_CD-M4", "kME_CD-M5", 
                 "kME_CD-M6", "kME_CD-M7", "kME_CD-M8", "kME_CD-M9", "kME_CD-M10")

# Melt the df for ggplot
HL_df_melted <- melt(HL_modules_sorted, id.vars = c("gene_name", "module"), measure.vars = kME_columns)
HL_df_melted <-  na.omit(HL_df_melted)

# Function of Z-score Normalization
z_score <- function(x) {
  return((x - mean(x)) / sd(x))
}

# Z-score Normalization
HL_df_melted_z <- HL_df_melted %>%
  mutate(value_zscore = z_score(value))

# Modify df
HL_df_melted_z$variable <- gsub("kME_", "", HL_df_melted_z$variable)
HL_df_melted_z$module <- factor(HL_df_melted_z$module, levels = c("CD-M1", "CD-M2", "CD-M3", "CD-M4", "CD-M5", "CD-M6", "CD-M7", "CD-M8", "CD-M9", "CD-M10"))
HL_df_melted_z$variable <- factor(HL_df_melted_z$variable, levels = c("CD-M1", "CD-M2", "CD-M3", "CD-M4", "CD-M5", "CD-M6", "CD-M7", "CD-M8", "CD-M9", "CD-M10"))

# Visualzization (kME)
# pdf("HL_heatmap_kME.pdf", width = 4.8, height = 10)
ggplot(HL_df_melted_z, aes(x = variable, y = gene_name, fill = value_zscore)) + 
  geom_tile(color = "white", linewidth = 0.5) +  
  scale_fill_viridis_c(option = "viridis") +  
  labs(
    title = "Heatmap of kME Values",  
    x = "kME",  
    y = "Genes (Grouped by Module)",
    fill = "Z-Score"
  ) + 
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.4, color = "black"),  
    axis.text.y = element_text(size = 10, face = "bold", hjust = 1, vjust = 0.4, color = "black"),
    axis.title.x = element_text(margin = margin(t = 10)),  
    axis.title.y = element_text(margin = margin(r = 10)),  
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.5, "lines"),
    strip.text.y = element_text(angle = 0, face = "bold", size = 12),
    strip.text = element_blank(),
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12)
  ) +
  facet_grid(rows = vars(module), scales = "free_y", space = "free_y")
# dev.off()

# Create df for binary heatmap
module_disease_genes = unique(HL_df_melted_z$gene_name)
module_disease_df = subset(HL_db, subset = Gene %in% module_disease_genes)
module_disease_df = module_disease_df[,1:2]

# Rename column
disease_col = c("Gene", "Disease")
colnames(module_disease_df) = disease_col
module_disease_df$Disease = as.factor(module_disease_df$Disease)

genes <- unique(module_disease_df$Gene)
diseases <- levels(module_disease_df$Disease)
expand_disease_df <- expand.grid(Gene = genes, Disease = diseases)

# Add column to match gene and disease
expand_disease_df <- expand_disease_df %>%
  mutate(Match = ifelse(paste(Gene, Disease) %in% paste(module_disease_df$Gene, module_disease_df$Disease), 1, 0))

df_unique <- HL_df_melted_z %>%
  distinct(gene_name, module)

df_wide <- df_unique %>%
  group_by(gene_name) %>%
  summarize(variable = paste(module, collapse = ", "), .groups = 'drop')

# Create df for binary heatmap
expand_disease_df <- expand_disease_df %>%
  left_join(df_wide %>% dplyr::select(gene_name, variable), by = c("Gene" = "gene_name"))

# Visualzization (binary)
# pdf("Disease_term_genes_binary.pdf", width = 4.5, height = 10)
ggplot(expand_disease_df, aes(x = Disease, y = Gene, fill = factor(Match))) + 
  geom_tile(color = "white", linewidth = 0.5) + 
  scale_fill_manual(values = c("0" = "gray", "1" = "black"), name = "Match") +  
  theme_minimal(base_size = 14) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
    axis.text.x = element_text(size = 9, face = "bold", angle = 90, hjust = 1, vjust = 1, color = "black"),  
    axis.text.y = element_text(size = 9, face = "bold", hjust = 1, vjust = 0.5, color = "black"),
    axis.title.x = element_text(margin = margin(t = 10)), 
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid = element_blank(), 
    panel.spacing.y = unit(0.5, "lines"),
    strip.text.y = element_text(angle = 0, face = "bold", size = 12),
    strip.text = element_blank(),
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12),
  ) +
  facet_grid(rows = vars(variable), scales = "free_y", space = "free_y") + NoLegend()
# dev.off()

########################################
## Disease genes correlation analysis ##
########################################

# Df : genes & modules
anno_df <- HL_df_melted_z[1:35 ,1:2]

# Info of Gene lists & modules
genes_anno <- anno_df$gene_name
modules_anno <- anno_df$module

# Extract expression values from seurat object
expr_data <- FetchData(CD_seu_harmony_sub2, vars = genes_anno)

# Calculate pairwise correlation
cor_matrix <- cor(expr_data, method = "pearson")

# Assign module colors
module_colors <- c('CD-M1' = "#FF7F00",
                   'CD-M2' = "#4DAF4A",
                   'CD-M3' = "#E41A1C",
                   'CD-M4' = '#33CCCC',
                   'CD-M5' = '#003333',
                   'CD-M6' = '#FF9999',
                   'CD-M7' = '#333366',
                   'CD-M8' = '#6666CC',
                   'CD-M9' = '#CC6699')
annotation_colors <- list(module = module_colors)

# Add module info to annotation
annotation <- data.frame(module = modules_anno)
rownames(annotation) <- genes_anno
annotation

# Correlation matrix of Disease associated genes
p1 = pheatmap(cor_matrix,
              annotation_row = annotation,
              annotation_colors = annotation_colors,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              clustering_method = "complete",
              display_numbers = TRUE,         # 상관 계수를 히트맵에 표시
              number_format = "%.2f",         # 소수점 두 자리까지 표시
              main = "Pairwise Gene Expression Correlation")

# pdf("HL_genes_Correlation_matrix.pdf", width = 12, height = 11)
p1
# dev.off()