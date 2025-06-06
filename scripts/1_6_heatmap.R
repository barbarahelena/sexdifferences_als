## Heatmap metabolites
## Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(rio)
library(ggsci)

# Load data
met <- readRDS("data/metabolome/metabolomics.RDS")
meta <- readRDS("data/metadata.RDS")

res_ctrl <- rio::import("results/metabolomics/ttests/metabolites_welcht_ctrl_sex_sig.csv")
res_tdp43 <- rio::import("results/metabolomics/ttests/metabolites_welcht_tdp_sex_sig.csv")
sig_metabolites <- unique(c(res_ctrl$metabolite, res_tdp43$metabolite))

# Join the results for female mice only (to deduplicate tables)
resfull_ctrl <- rio::import("results/metabolomics/ttests/metabolites_welcht_ctrl_sex_diff.csv")
resfull_tdp43 <- rio::import("results/metabolomics/ttests/metabolites_welcht_tdp_sex_sig.csv")
combined_results <- resfull_ctrl %>% 
  filter(Sex == "Female" & metabolite %in% sig_metabolites) %>% 
  select(metabolite, p.value, q.value, Sex, sig) %>%
  full_join(
    resfull_tdp43 %>% 
      filter(Sex == "Female" & metabolite %in% sig_metabolites) %>% 
      select(metabolite, p.value, q.value, sig),
    by = "metabolite", 
    suffix = c("_ctrl", "_tdp43")
  )
write.csv(combined_results$metabolite, "results/metabolomics/ttests/sig_metabolites.csv", row.names = FALSE)

# Prepare metabolite data matrix
met_matrix <- met %>%
  select(all_of(sig_metabolites)) %>%
  mutate(across(everything(.), ~ scale(.x)))%>%
  as.matrix() %>% t()

# Create p-value annotation matrix
pval <- combined_results %>% select(metabolite, p.value_ctrl, p.value_tdp43, sig_ctrl, sig_tdp43)

# Get sample metadata
sample_info <- meta %>%
  select(ID, Sex, Intervention) %>%
  arrange(Sex, Intervention)
met_matrix <- met_matrix[,sample_info$ID]
pval <- pval[match(rownames(met_matrix), pval$metabolite),]

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
                left_annotation = rowAnnotation(pvalue_tdp43 = anno_simple(-log10(pval$p.value_tdp43), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pval$sig_tdp43,
                                                                     pt_size = unit(1, "snpc")*0.7
                                                                     ),
                                                pvalue_ctrl = anno_simple(-log10(pval$p.value_ctrl), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pval$sig_ctrl,
                                                                     pt_size = unit(1, "snpc")*0.7
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
lgd_sig = Legend(pch = c("*","**","***","****"), type = "points", 
                    labels = c("0.05", "0.01", "0.001", "<0.001"),
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
                left_annotation = rowAnnotation(pvalue_tdp43 = anno_simple(-log10(pval$p.value_tdp43), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pval$sig_tdp43,
                                                                     pt_size = unit(1, "snpc")*0.7
                                                                     ),
                                                pvalue_ctrl = anno_simple(-log10(pval$p.value_ctrl), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pval$sig_ctrl,
                                                                     pt_size = unit(1, "snpc")*0.7
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
lgd_sig = Legend(pch = c("*","**","***","****"), type = "points", 
                    labels = c("0.05", "0.01", "0.001", "<0.001"),
                    legend_gp = gpar(fontsize = 8))
pl2_anno <- draw(pl2, annotation_legend_list = list(lgd_pvalue, lgd_sig))

pdf("results/metabolomics/heatmaps/heatmap_samplenames.pdf", width = 8, height = 10)
pl2_anno
dev.off()
svg("results/metabolomics/heatmaps/heatmap_samplenames.svg", width = 8, height = 10)
pl2_anno
dev.off()

## Ggarrange all metabolomics plots ##
heatmap_grob <- grid.grabExpr(draw(
      pl1_anno,
      padding = unit(c(2, 10, 10, 2), "mm"), # top, right, bottom, left padding
  )
)
ggarrange(ggarrange(pca1, pca2, nrow = 1, labels = LETTERS[1:2], widths = c(1,1.5)),
          ggarrange(plgeno, plwt, pltdp, nrow = 1, labels = LETTERS[3:5]),
          as_ggplot(heatmap_grob), 
          boxplots,
          nrow = 4, heights = c(1, 1, 1.5, 1.5), labels = c("", "", LETTERS[6], ""))
ggsave("results/metabolomics/assembled_figure.pdf", width = 14, height = 20)

ggarrange(ggarrange(pca1, pca2, nrow = 1, labels = LETTERS[1:2], widths = c(1,1.5)),
          ggarrange(plgeno, plwt, pltdp, nrow = 1, labels = LETTERS[3:5]),
          ggarrange(boxplots, as_ggplot(heatmap_grob), labels = c("", LETTERS[18])),
          nrow = 3, heights = c(0.7, 0.7, 1.8))
ggsave("results/metabolomics/assembled_figure_2.pdf", width = 18, height = 21)
