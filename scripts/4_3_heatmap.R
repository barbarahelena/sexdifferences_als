## Heatmap metabolites
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(rio)
library(ggsci)

# Load data
met <- readRDS("data/metabolome/metabolomics.RDS")
meta <- readRDS("data/metadata.RDS")

# Compute Sex x Genotype interaction per metabolite (two-way ANOVA) on all metabolites
met_with_meta <- met %>%
  rownames_to_column("ID") %>%
  left_join(meta %>% select(ID, Sex, Intervention), by = "ID")

interaction_results <- map_dfr(colnames(met), function(m) {
  df <- met_with_meta %>% select(ID, Sex, Intervention, value = all_of(m))
  fit <- lm(value ~ Sex * Intervention, data = df)
  p_int <- anova(fit)["Sex:Intervention", "Pr(>F)"]
  tibble(metabolite = m, p_interaction = p_int)
}) %>%
  mutate(
    q_interaction = p.adjust(p_interaction, method = "BH"),
    sig_interaction = case_when(
      p_interaction < 0.001 ~ "***",
      p_interaction < 0.01  ~ "**",
      p_interaction < 0.05  ~ "*",
      TRUE ~ ""
    )
  )
write.csv(interaction_results, "results/metabolomics/ttests/interaction_lm_results.csv", row.names = FALSE)

# Filter: metabolites with significant Sex x Genotype interaction (nominal p < 0.05)
# Note: with n=5/group FDR correction yields no significant results
sig_metabolites <- interaction_results %>%
  filter(p_interaction < 0.05) %>%
  pull(metabolite)

# Prepare metabolite data matrix
met_matrix <- met %>%
  select(all_of(sig_metabolites)) %>%
  mutate(across(everything(.), ~ scale(.x)))%>%
  as.matrix() %>% t()

# Get sample metadata
sample_info <- meta %>%
  select(ID, Sex, Intervention) %>%
  arrange(Sex, Intervention)
met_matrix <- met_matrix[,sample_info$ID]
pval <- interaction_results[match(rownames(met_matrix), interaction_results$metabolite), ]

# Define color functions for heatmap
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c(pal_nejm()(6)[6], "white", pal_nejm()(3)[3])) # color mapping for blocks heatmap
pvalue_col_fun = colorRamp2(c(0, 1, 3), c("lightgrey", "white", pal_nejm()(5)[5])) # color mapping for -log10(pvalue)

(pl1 <- Heatmap(met_matrix, name = "Z-score", 
                col = col_fun, 
                rect_gp = gpar(col = "white", lwd = 2),
                column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                column_split = c(rep("WT (M)", 5), rep("TDP43 (M)", 5),
                                    rep("WT (F)", 5), rep("TDP43 (F)", 5)),,
                column_order = sample_info$ID,
                cluster_columns = FALSE,
                clustering_method_rows = "ward.D",
                left_annotation = rowAnnotation(
                                                interaction = anno_simple(-log10(pval$p_interaction),
                                                                     which = 'row',
                                                                     col = pvalue_col_fun,
                                                                     pch = pval$sig_interaction,
                                                                     pt_gp = gpar(fontsize = 7)
                                                                     ),
                                                annotation_name_side = "bottom",
                                                width = unit(0.8, "cm"),
                                                annotation_name_gp = grid::gpar(fontsize = 8)),
                # width = unit(6, "cm"),
                # height = unit(0.25*28, "cm"),
                row_names_side = "left",
                column_gap = unit(2, "mm"),
                row_dend_side = "right",
                show_column_names = FALSE,
                row_names_gp = grid::gpar(fontsize = 10))
                )
               
lgd_pvalue = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3), 
                    labels = c("1", "0.1", "0.01", "0.001"))
lgd_sig = Legend(pch = c("*","**","***"), type = "points",
                    labels = c("<0.05", "<0.01", "<0.001"),
                    legend_gp = gpar(fontsize = 8))
pl1_anno <- draw(pl1, annotation_legend_list = list(lgd_pvalue, lgd_sig))

pdf("results/metabolomics/heatmaps/heatmap.pdf", width = 8, height = 10)
pl1_anno
dev.off()
svg("results/metabolomics/heatmaps/heatmap.svg", width = 8, height = 10)
pl1_anno
dev.off()

pdf("results/metabolomics/heatmaps/heatmap_wide.pdf", width = 8, height = 8)
pl1_anno
dev.off()
svg("results/metabolomics/heatmaps/heatmap_wide.svg", width = 10, height = 10)
pl1_anno
dev.off()


# With colnames as check
(pl2 <- Heatmap(met_matrix, name = "Z-score", 
                col = col_fun, 
                rect_gp = gpar(col = "white", lwd = 2),
                column_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                column_split = c(rep("WT (M)", 5), rep("TDP43 (M)", 5),
                                    rep("WT (F)", 5), rep("TDP43 (F)", 5)),,
                column_order = sample_info$ID,
                cluster_columns = FALSE,
                clustering_method_rows = "average",
                left_annotation = rowAnnotation(
                                                interaction = anno_simple(-log10(pval$p_interaction),
                                                                     which = 'row',
                                                                     col = pvalue_col_fun,
                                                                     pch = pval$sig_interaction,
                                                                     pt_gp = gpar(fontsize = 7)
                                                                     ),
                                                annotation_name_side = "bottom",
                                                width = unit(0.8, "cm"),
                                                annotation_name_gp = grid::gpar(fontsize = 8)),
                # width = unit(6, "cm"), 
                # height = unit(0.25*28, "cm"), 
                row_names_side = "left", 
                column_gap = unit(2, "mm"),
                row_dend_side = "right",
                show_column_names = TRUE,
                row_names_gp = grid::gpar(fontsize = 7),
                column_names_gp = grid::gpar(fontsize = 7)
                )
)               
lgd_pvalue = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3), 
                    labels = c("1", "0.1", "0.01", "0.001"))
lgd_sig = Legend(pch = c("*","**","***"), type = "points",
                    labels = c("<0.05", "<0.01", "<0.001"),
                    legend_gp = gpar(fontsize = 8))
pl2_anno <- draw(pl2, annotation_legend_list = list(lgd_pvalue, lgd_sig))

pdf("results/metabolomics/heatmaps/heatmap_samplenames.pdf", width = 8, height = 10)
pl2_anno
dev.off()
svg("results/metabolomics/heatmaps/heatmap_samplenames.svg", width = 8, height = 10)
pl2_anno
dev.off()

## Ggarrange all pathway and metabolomics plots ##
heatmap_grob <- grid.grabExpr(draw(
      pl1_anno,
      padding = unit(c(2, 10, 10, 2), "mm"), # top, right, bottom, left padding
  )
)

(metabolomics <- ggarrange(pca1, pca2, as_ggplot(heatmap_grob), nrow = 1, labels = LETTERS[4:6], widths = c(1, 1.5, 3.0)))


# fig <- ggarrange(ggarrange(pca1, pca2, nrow = 1, labels = LETTERS[1:2], widths = c(1,1.5)),
#           ggarrange(plgeno, plwt, pltdp, nrow = 1, labels = LETTERS[3:5]),
#           as_ggplot(heatmap_grob), 
#           boxplots,
#           nrow = 4, heights = c(1, 1, 1.5, 1.5), labels = c("", "", LETTERS[6], ""))
# ggsave("results/metabolomics/assembled_figure.pdf", plot = fig, width = 14, height = 20)

# fig <- ggarrange(ggarrange(pca1, pca2, nrow = 1, labels = LETTERS[1:2], widths = c(1,1.5)),
#           ggarrange(plgeno, plwt, pltdp, nrow = 1, labels = LETTERS[3:5]),
#           ggarrange(boxplots, as_ggplot(heatmap_grob), labels = c("", LETTERS[18]), nrow = 1),
#           nrow = 3, heights = c(0.7, 0.7, 1.8))
# ggsave("results/metabolomics/assembled_figure_2.pdf", plot = fig, width = 18, height = 21)
