# Import packages
library(igraph)
library(ggraph)

##################################
## Correlation analysis (Foxg1) ##
##################################

# Subset all gene names
gene_list = rownames(CD_seu_harmony_sub2@assays$SCT$counts)

# Make empty data.frame
Corr_df_Foxg1 = data.frame(Gene = character(0), Corr = numeric(0), P_Value = numeric(0))

# Calculate Foxg1 expression and all other genes
for (gene in gene_list) {
  Foxg1_expression = CD_seu_harmony_sub2@assays$SCT$counts["Foxg1",]
  other_gene_expression = CD_seu_harmony_sub2@assays$SCT$counts[gene, ]
  correlation_result = cor.test(Foxg1_expression, other_gene_expression)
  Corr_df_Foxg1 = rbind(Corr_df_Foxg1, data.frame(Gene = gene, Correlation = correlation_result$estimate, P_Value = correlation_result$p.value))
}

# Calculate Bonferonni P value
Corr_df_Foxg1$Bon_P = p.adjust(Corr_df_Foxg1$P_Value, method = "bonferroni")

# Calculate -log10 P value
Corr_df_Foxg1$Negative_log10_P = -log10(Corr_df_Foxg1$P_Value)

# Subset Bon_p < 0.05 genes
Corr_df_Foxg1_Bon_P = subset(Corr_df_Foxg1, Bon_P < 0.05)
head(Corr_df_Foxg1_Bon_P, n = 50)

# Sorting
Corr_df_Foxg1_Bon_P_sorted = Corr_df_Foxg1_Bon_P[order(Corr_df_Foxg1_Bon_P$Bon_P), ]
head(Corr_df_Foxg1_Bon_P_sorted, n = 20)

###################################
########## Foxg1 Network ##########
###################################

# Subset df Corr >0.15
Corr_df_Foxg1_0.15 = subset(Corr_df_Foxg1, subset = Correlation > 0.15)

# Subset Foxg1 associated module & correlation > 0.15 genes 
modules = GetModules(CD_seu_harmony_sub2)
module_Foxg1 = subset(modules, subset = gene_name %in% Corr_df_Foxg1_0.15$Gene) # Overlapping genes
# module_Foxg1 = subset(module_Foxg1, subset = module %in% c("CD-M3", "CD-M6")) # Foxg1 associated module : M3 & M6

# Merge df
Corr_df_Foxg1_0.15 = merge(Corr_df_Foxg1_0.15, module_Foxg1, by.x = "Gene", by.y = "gene_name", all.x = TRUE)
Corr_df_Foxg1_0.15 = Corr_df_Foxg1_0.15[,c(1,3:7)]

# Sorting
Corr_df_Foxg1_0.15 = Corr_df_Foxg1_0.15[order(Corr_df_Foxg1_0.15$Correlation, decreasing = T), ]

# Subset Foxg1 associated modules & Top 20 Correlated genes
Corr_df_Foxg1_0.15 = subset(Corr_df_Foxg1_0.15, subset = module %in% c("CD-M3", "CD-M6"))
Corr_df_Foxg1_0.15_top21 = Corr_df_Foxg1_0.15[1:21,]
Corr_df_Foxg1_0.15_top21

# Create Foxg1 Networks
# Edges
edges = Corr_df_Foxg1_0.15_top21 %>%
  filter(Gene != "Foxg1") %>%
  mutate(from = "Foxg1", to = Gene) %>%
  select(from, to)

# Nodes
nodes = Corr_df_Foxg1_0.15_top21 %>%
  mutate(color = ifelse(module == "CD-M3", "#E41A1C", "#FF9999")) %>%
  select(name = Gene, color)

# Create igraph object
g = graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
g

# Network Plot
# pdf("Foxg1_network_box_padding.pdf", width = 5, height = 5)
ggraph(g, layout = 'fr') + 
  geom_edge_link(aes(edge_alpha = 0.2), show.legend = FALSE) + 
  geom_node_point(aes(x = x, y = y, color = color, fill = color, size = 5), 
                  show.legend = FALSE, shape = 21, color = "black", stroke = 1) +  
  geom_label_repel(aes(x = x, y = y, label = name), 
                   box.padding = unit(0.5, "lines"),
                   size = 5,  
                   fontface = "bold",  
                   label.size = 0.2,  
                   fill = "white",  
                   color = "black",  
                   label.padding = unit(0.15, "lines"), 
                   segment.color = "grey50") +  
  theme_void() +
  scale_color_identity() +
  scale_fill_identity() +  
  scale_size_continuous(range = c(5, 5)) +  
  labs(title = "Foxg1 Network")
# dev.off()

##################################
## Correlation analysis (Gata3) ##
##################################

# Subset all gene names
gene_list = rownames(CD_seu_harmony_sub2@assays$SCT$counts)

# Make empty data.frame
Corr_df_Gata3 = data.frame(Gene = character(0), Corr = numeric(0), P_Value = numeric(0))

# Calculate Gata3 expression and all other genes
for (gene in gene_list) {
  Gata3_expression = CD_seu_harmony_sub2@assays$SCT$counts["Gata3",]
  other_gene_expression = CD_seu_harmony_sub2@assays$SCT$counts[gene, ]
  correlation_result = cor.test(Gata3_expression, other_gene_expression)
  Corr_df_Gata3 = rbind(Corr_df_Gata3, data.frame(Gene = gene, Correlation = correlation_result$estimate, P_Value = correlation_result$p.value))
}

# Valculate Bonferonni P value
Corr_df_Gata3$Bon_P = p.adjust(Corr_df_Gata3$P_Value, method = "bonferroni")

# Calculate -log10 P value
Corr_df_Gata3$Negative_log10_P = -log10(Corr_df_Gata3$P_Value)

# Subset Bon_p < 0.05 genes
Corr_df_Gata3_Bon_P = subset(Corr_df_Gata3, Bon_P < 0.05)
head(Corr_df_Gata3_Bon_P, n = 50)

# Sorting
Corr_df_Gata3_Bon_P_sorted = Corr_df_Gata3_Bon_P[order(Corr_df_Gata3_Bon_P$Bon_P), ]
head(Corr_df_Gata3_Bon_P_sorted, n = 20)

###################################
########## Gata3 Network ##########
###################################

# Subset df Corr >0.15
Corr_df_Gata3_0.15 = subset(Corr_df_Gata3, subset = Correlation > 0.15)

# Subset Gata3 associated module & correlation > 0.15 genes 
modules = GetModules(CD_seu_harmony_sub2)
module_Gata3 = subset(modules, subset = gene_name %in% Corr_df_Gata3_0.15$Gene) # Overlapping genes
# module_Gata3 = subset(module_Gata3, subset = module %in% c("CD-M3", "CD-M6")) # Gata3 associated module : M3 & M6

# Merge df
Corr_df_Gata3_0.15 = merge(Corr_df_Gata3_0.15, module_Gata3, by.x = "Gene", by.y = "gene_name", all.x = TRUE)
Corr_df_Gata3_0.15 = Corr_df_Gata3_0.15[,c(1,3:7)]

# Sorting
Corr_df_Gata3_0.15 = Corr_df_Gata3_0.15[order(Corr_df_Gata3_0.15$Correlation, decreasing = T), ]
Corr_df_Gata3_0.15
# Subset Gata3 associated modules & Top 20 Correlated genes
Corr_df_Gata3_0.15 = subset(Corr_df_Gata3_0.15, subset = module %in% c("CD-M2"))
Corr_df_Gata3_0.15_top21 = Corr_df_Gata3_0.15[1:21,]
Corr_df_Gata3_0.15_top21

# Create Foxg1 Networks
# Edges
edges = Corr_df_Gata3_0.15_top21 %>%
  filter(Gene != "Gata3") %>%
  mutate(from = "Gata3", to = Gene) %>%
  select(from, to)

# Nodes
nodes = Corr_df_Gata3_0.15_top21 %>%
  mutate(color = ifelse(module == "CD-M2", "#4DAF4A", "white")) %>%
  select(name = Gene, color)

# Create igraph object
g = graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
g

# Network Plot
# pdf("Gata3_network_box_padding.pdf", width = 5, height = 5)
ggraph(g, layout = 'fr') + 
  geom_edge_link(aes(edge_alpha = 0.2), show.legend = FALSE) + 
  geom_node_point(aes(x = x, y = y, color = color, fill = color, size = 5), 
                  show.legend = FALSE, shape = 21, color = "black", stroke = 1) +  
  geom_label_repel(aes(x = x, y = y, label = name), 
                   box.padding = unit(0.5, "lines"),
                   size = 5,  
                   fontface = "bold",  
                   label.size = 0.2,  
                   fill = "white",  
                   color = "black",  
                   label.padding = unit(0.15, "lines"), 
                   segment.color = "grey50") +  
  theme_void() +
  scale_color_identity() +
  scale_fill_identity() +  
  scale_size_continuous(range = c(5, 5)) +  
  labs(title = "Gata3 Network")
# dev.off()
