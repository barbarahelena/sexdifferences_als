## Heatmap metabolites
## Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(rio)
library(ggsci)

## Data
df <- readRDS("data/metadata.RDS")
met <- readRDS("data/metabolomics.RDS")
res1 <- rio::import("results/metabolomics/ttests/metabolites_welcht_ctrl_sig.csv")
res2 <- rio::import("results/metabolomics/ttests/metabolites_welcht_sig.csv")
totmetlist <- c(res1$metabolite, res2$metabolite)
resfull1 <- rio::import("results/metabolomics/ttests/metabolites_welcht_ctrl_diff.csv")
resfull2 <- rio::import("results/metabolomics/ttests/metabolites_welcht_diff.csv")
totres <- full_join(resfull1 %>% filter(Sex == "Female" & metabolite %in% totmetlist) %>% select(1:4,7), 
                    resfull2 %>% filter(Sex == "Female" & metabolite %in% totmetlist) %>% select(1:4,7), 
                    by = "metabolite", suffix = c("_ctrl", "_tdp43"))
totres
met <- met %>% select(any_of(totmetlist))
metanno <- met %>% mutate(ID = rownames(.)) %>% left_join(df, by = "ID") %>% arrange(Sex, Intervention)

## Heatmap with all mice and metabolites with sex differences in ctrl/tdp43
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c(pal_nejm()(2)[2], "white", pal_nejm()(1))) # color mapping for blocks heatmap
pvalue_col_fun = colorRamp2(c(0, 1, 3), c("lightgrey", "white", pal_nejm()(5)[5])) # color mapping for -log10(pvalue)
matmet <- t(as.matrix(met)) # matrix of all mets (only works with matrices)
mat <- as.data.frame(matmet)
mat$metabolite <- rownames(mat)
pvalanno <- left_join(mat, totres, by = "metabolite")
matmet <- matmet[,metanno$ID]

(pl1 <- Heatmap(matmet, name = "Z-score", 
                col = col_fun, 
                rect_gp = gpar(col = "white", lwd = 2),
                column_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                column_split = c(rep("WT (M)", 5), rep("TDP43 (M)", 5),
                                    rep("WT (F)", 5), rep("TDP43 (F)", 5)),  
                column_order = metanno$ID,
                cluster_columns = FALSE,
                clustering_method_rows = "average",
                top_annotation = HeatmapAnnotation(Group = anno_block(gp = gpar(fill = c(pal_nejm()(4)))), 
                                                    height = unit(0.2, "cm")),
                left_annotation = rowAnnotation(pvalue_tdp43 = anno_simple(-log10(pvalanno$p.value_tdp43), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pvalanno$sig_tdp43,
                                                                     pt_size = unit(1, "snpc")*0.5
                                                                     ),
                                                pvalue_ctrl = anno_simple(-log10(pvalanno$p.value_ctrl), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pvalanno$sig_ctrl,
                                                                     pt_size = unit(1, "snpc")*0.5
                                                                     ),
                                                annotation_name_side = "bottom",
                                                width = unit(0.5, "cm"),
                                                annotation_name_gp = grid::gpar(fontsize = 8)),
                width = unit(6, "cm"), 
                height = unit(0.25*28, "cm"), 
                row_names_side = "left", 
                column_gap = unit(2, "mm"),
                row_dend_side = "right",
                #show_column_names = FALSE,
                row_names_gp = grid::gpar(fontsize = 7))
                )

lgd_pvalue = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3), 
                    labels = c("1", "0.1", "0.01", "0.001"))
lgd_sig = Legend(pch = c("*","**","***","****"), type = "points", labels = c("0.05", "0.01", "0.001", "<0.001"))
pl1_anno <- draw(pl1, annotation_legend_list = list(lgd_pvalue, lgd_sig))

pdf("results/heatmaps/heatmap_fourgroups.pdf", width = 8, height = 10)
pl1_anno
dev.off()
svg("results/heatmaps/heatmap_fourgroups.svg", width = 8, height = 10)
pl1_anno
dev.off()

# Different color scheme
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c(pal_nejm()(6)[6], "white", pal_nejm()(3)[3])) # color mapping for blocks heatmap
pvalue_col_fun = colorRamp2(c(0, 1, 3), c("lightgrey", "white", pal_nejm()(5)[5])) # color mapping for -log10(pvalue)
matmet <- matmet[,order(metanno$ID)]

(pl1 <- Heatmap(matmet, name = "Z-score", 
                col = col_fun, 
                rect_gp = gpar(col = "white", lwd = 2),
                column_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                column_split = c(rep("WT (M)", 5), rep("TDP43 (M)", 5),
                                    rep("WT (F)", 5), rep("TDP43 (F)", 5)),,
                column_order = metanno$ID,
                cluster_columns = FALSE,
                clustering_method_rows = "average",
                top_annotation = HeatmapAnnotation(
                        Group = anno_block(gp = gpar(fill = c(pal_nejm()(7)[c(4,7,4,7)]),
                                                        alpha = rep(0.9, 4))), 
                        Sex = anno_block(gp = gpar(fill = c(pal_nejm()(2)[c(2,2,1,1)]),
                                                        alpha = rep(0.9, 4))),
                        height = unit(0.4, "cm")
                        ),
                left_annotation = rowAnnotation(pvalue_tdp43 = anno_simple(-log10(pvalanno$p.value_tdp43), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pvalanno$sig_tdp43,
                                                                     pt_size = unit(1, "snpc")*0.7
                                                                     ),
                                                pvalue_ctrl = anno_simple(-log10(pvalanno$p.value_ctrl), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pvalanno$sig_ctrl,
                                                                     pt_size = unit(1, "snpc")*0.7
                                                                     ),
                                                annotation_name_side = "bottom",
                                                width = unit(0.8, "cm"),
                                                annotation_name_gp = grid::gpar(fontsize = 8)),
                width = unit(6, "cm"), 
                height = unit(0.25*28, "cm"), 
                row_names_side = "left", 
                column_gap = unit(2, "mm"),
                row_dend_side = "right",
                show_column_names = FALSE,
                row_names_gp = grid::gpar(fontsize = 7))
                )
                
lgd_pvalue = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3), 
                    labels = c("1", "0.1", "0.01", "0.001"))
lgd_sex = Legend(title = "Sex", labels = c("Female", "Male"), 
legend_gp = gpar(fill = c(rev(pal_nejm()(2)))))
lgd_group = Legend(title = "Mice", labels = c("Control", "TDP43"), 
legend_gp = gpar(fill = c(pal_nejm()(7)[c(4,7)])))
lgd_sig = Legend(pch = c("*","**","***","****"), type = "points", 
                    labels = c("0.05", "0.01", "0.001", "<0.001"), grid_width = unit(5, "mm"))
pl1_anno <- draw(pl1, annotation_legend_list = list(lgd_pvalue, lgd_group, lgd_sex))

pdf("results/heatmaps/heatmap_diffcolors.pdf", width = 8, height = 10)
pl1_anno
dev.off()
svg("results/heatmaps/heatmap_diffcolors.svg", width = 8, height = 10)
pl1_anno
dev.off()

(pl1 <- Heatmap(matmet, name = "Z-score", 
                col = col_fun, 
                rect_gp = gpar(col = "white", lwd = 2),
                column_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                column_split = c(rep("WT (M)", 5), rep("TDP43 (M)", 5),
                                    rep("WT (F)", 5), rep("TDP43 (F)", 5)),,
                column_order = metanno$ID,
                cluster_columns = FALSE,
                clustering_method_rows = "average",
                left_annotation = rowAnnotation(pvalue_tdp43 = anno_simple(-log10(pvalanno$p.value_tdp43), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pvalanno$sig_tdp43,
                                                                     pt_size = unit(1, "snpc")*0.7
                                                                     ),
                                                pvalue_ctrl = anno_simple(-log10(pvalanno$p.value_ctrl), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pvalanno$sig_ctrl,
                                                                     pt_size = unit(1, "snpc")*0.7
                                                                     ),
                                                annotation_name_side = "bottom",
                                                width = unit(0.8, "cm"),
                                                annotation_name_gp = grid::gpar(fontsize = 8)),
                width = unit(6, "cm"), 
                height = unit(0.25*28, "cm"), 
                row_names_side = "left", 
                column_gap = unit(2, "mm"),
                row_dend_side = "right",
                show_column_names = FALSE,
                row_names_gp = grid::gpar(fontsize = 7))
                )
               
lgd_pvalue = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3), 
                    labels = c("1", "0.1", "0.01", "0.001"))
lgd_sig = Legend(pch = c("*","**","***","****"), type = "points", 
                    labels = c("0.05", "0.01", "0.001", "<0.001"),
                    legend_gp = gpar(fontsize = 8))
pl1_anno <- draw(pl1, annotation_legend_list = list(lgd_pvalue, lgd_sig))

pdf("results/heatmaps/heatmap_noanno.pdf", width = 8, height = 10)
pl1_anno
dev.off()
svg("results/heatmaps/heatmap_noanno.svg", width = 8, height = 10)
pl1_anno
dev.off()

## Heatmap with means per group and metabolites with sex differences in ctrl/tdp43
totres2 <- full_join(resfull1, resfull2, 
                    by = c("metabolite", "Sex"), suffix = c("_ctrl", "_tdp43")) %>% 
                    pivot_wider(., names_from = "Sex", 
                        values_from = c("mean_ctrl", "p.value_ctrl", "mean_tdp43", "p.value_tdp43"))
rownames(totres2) <- totres2$metabolite
totres2 <- totres %>% select(contains("mean"), contains("p.value"))
(pl1 <- Heatmap(matmet, name = "Z-score", 
                col = col_fun, 
                rect_gp = gpar(col = "white", lwd = 2),
                column_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                column_split = c(rep("WT (M)", 5), rep("TDP43 (M)", 5),
                                    rep("WT (F)", 5), rep("TDP43 (F)", 5)),,
                column_order = metanno$ID,
                cluster_columns = FALSE,
                clustering_method_rows = "average",
                left_annotation = rowAnnotation(pvalue_tdp43 = anno_simple(-log10(pvalanno$p.value_tdp43), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pvalanno$sig_tdp43,
                                                                     pt_size = unit(1, "snpc")*0.7
                                                                     ),
                                                pvalue_ctrl = anno_simple(-log10(pvalanno$p.value_ctrl), 
                                                                     which = 'row',
                                                                     col = pvalue_col_fun, 
                                                                     pch = pvalanno$sig_ctrl,
                                                                     pt_size = unit(1, "snpc")*0.7
                                                                     ),
                                                annotation_name_side = "bottom",
                                                width = unit(0.8, "cm"),
                                                annotation_name_gp = grid::gpar(fontsize = 8)),
                width = unit(6, "cm"), 
                height = unit(0.25*28, "cm"), 
                row_names_side = "left", 
                column_gap = unit(2, "mm"),
                row_dend_side = "right",
                show_column_names = FALSE,
                row_names_gp = grid::gpar(fontsize = 7))
                )