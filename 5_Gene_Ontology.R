# Import library
library(GOSemSim)
library(org.Mm.eg.db) # Mus musculus
library(AnnotationHub)
library(clusterProfiler)

#################################
## GO analysis for the modules ##
#################################

# Subset each modules
modules_M1 = subset(modules, subset = module %in% c("CD-M1"))
modules_M2 = subset(modules, subset = module %in% c("CD-M2"))
modules_M3 = subset(modules, subset = module %in% c("CD-M3"))
modules_M4 = subset(modules, subset = module %in% c("CD-M4"))
modules_M5 = subset(modules, subset = module %in% c("CD-M5"))
modules_M6 = subset(modules, subset = module %in% c("CD-M6"))
modules_M7 = subset(modules, subset = module %in% c("CD-M7"))
modules_M8 = subset(modules, subset = module %in% c("CD-M8"))
modules_M9 = subset(modules, subset = module %in% c("CD-M9"))

# GO analysis (BP; Biological Process)
M1_BP = enrichGO(gene          = modules_M1$gene_name,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

M2_BP = enrichGO(gene          = modules_M2$gene_name,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

M3_BP = enrichGO(gene          = modules_M3$gene_name,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

M4_BP = enrichGO(gene          = modules_M4$gene_name,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

M5_BP = enrichGO(gene          = modules_M5$gene_name,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

M6_BP = enrichGO(gene          = modules_M6$gene_name,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

M7_BP = enrichGO(gene          = modules_M7$gene_name,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

M8_BP = enrichGO(gene          = modules_M8$gene_name,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

M9_BP = enrichGO(gene          = modules_M9$gene_name,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

# Visualization
# Dotplot
p1 = dotplot(M1_BP, font.size = 15, showCategory = 20) +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = "#FF7F00") +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
  )

p2 = dotplot(M2_BP, font.size = 15, showCategory = 20) +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = "#4DAF4A") +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
  )

p3 = dotplot(M3_BP, font.size = 15, showCategory = 20) +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = "#E41A1C") +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
  )

p4 = dotplot(M4_BP, font.size = 15, showCategory = 20) +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = '#33CCCC') +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
  )

p5 = dotplot(M5_BP, font.size = 15, showCategory = 20) +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = '#003333') +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
  )

p6 = dotplot(M6_BP, font.size = 15, showCategory = 20) +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = '#FF9999') +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
  )

p7 = dotplot(M7_BP, font.size = 15, showCategory = 20) +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = '#333366') +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
  )

p8 = dotplot(M8_BP, font.size = 15, showCategory = 20) +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = '#6666CC') +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
  )

p9 = dotplot(M9_BP, font.size = 15, showCategory = 20) +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = '#CC6699') +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
  )

# pdf("Dotplot_module.pdf", width = 6.5, height = 8)
p1
p2
p3
p4
p5
p6
p7
p8
p9
# dev.off()

################################################
## GO analysis for the Foxg1 & Gata3 Networks ##
################################################

# GO analysis (BP; Biological Process)
Foxg1_BP = enrichGO(gene          = Corr_df_Foxg1_0.15$Gene[1:50],
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
Gata3_BP = enrichGO(gene          = Corr_df_Gata3_0.15$Gene[1:50],
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)

# pdf("Dotplot_Foxg1_Gata3.pdf", width = 6, height = 8)
dotplot(Foxg1_BP, font.size = 15, showCategory = 20) +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = "#E41A1C") +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
  )

dotplot(Gata3_BP, font.size = 15, showCategory = 20) +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = "#4DAF4A") +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10, face = "bold"),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")
  )
# dev.off()
