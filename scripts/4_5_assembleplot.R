# Assemble plot - cross-dataset figure (mouse pathways/CAZymes + human metabolomics)
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

#### Libraries ####
library(tidyverse)
library(ggsci)
library(ggpubr)
library(ComplexHeatmap)
library(ggthemes)
library(Cairo)

#### Define global variables required by 3_3_cazymes.R ####
# Mouse experiment timepoints (weeks); adjust if a subset had insufficient n
timepoints        <- c(6, 8, 10, 12, 14, 16)
tp_sufficient_tdp <- c(6, 8, 10, 12, 14, 16)  # TDP43 timepoints with enough samples per sex

#### Source scripts ####
source("scripts/3_2_pathways_microbes.R")  # -> heatmap_top, lgd_sig
lgd_sig_path <- lgd_sig                    # rename before next source overwrites lgd_sig
source("scripts/3_3_cazymes.R")           # -> line_plist_sex (GH78=[[1]], GH106=[[2]])
source("scripts/3_4_cazyme_richness.R")   # -> res_tdp43$plot (richness line plot)
source("scripts/4_2_pca_metabolomics.R")  # -> pca1, pca2
source("scripts/4_3_heatmap.R")           # -> pl1 (met heatmap), lgd_pvalue, lgd_sig
lgd_sig_met    <- lgd_sig
lgd_pvalue_met <- lgd_pvalue

#### Convert ComplexHeatmaps to ggplot grobs ####
heatmap_gg_path <- as_ggplot(grid.grabExpr(
  draw(
    heatmap_bugs_pathways,
    annotation_legend_list = list(lgd_sig_path),
    padding = unit(c(15, 65, 10, 2), "mm")
  )
))

heatmap_gg_met <- as_ggplot(grid.grabExpr(
  draw(
    pl1,
    annotation_legend_list = list(lgd_pvalue_met, lgd_sig_met),
    padding = unit(c(5, 10, 10, 2), "mm")
  )
))

#### Extract CAZyme sex line plots ####
pl_gh78_sex  <- line_plist_sex[[1]]
pl_gh106_sex <- line_plist_sex[[2]]

#### Extract CAZyme Shannon diversity line plots ####
pl_shannon_tdp43 <- res_tdp43_shannon$plot
pl_shannon_wt    <- res_wt_shannon$plot

#### Assemble figure ####
# Left col:  A) full pathway-microbe heatmap (top)
#            H) metabolite heatmap (bottom)
# Right col: B) GH78 | C) GH106
#            D) CAZyme diversity TDP43 | E) CAZyme diversity WT
#            F) PCA total
#            G) PCA by sex

left_col <- ggarrange(
  heatmap_gg_path,
  heatmap_gg_met,
  nrow    = 2,
  labels  = c("A", "H"),
  heights = c(2.5, 0.8)
)

right_col <- ggarrange(
  ggarrange(pl_gh78_sex, pl_gh106_sex, ncol = 2, labels = LETTERS[2:3]),
  ggarrange(pl_shannon_tdp43, pl_shannon_wt, ncol = 2, labels = c("D", "E"),
            common.legend = TRUE, legend = "bottom"),
  pca1,
  pca2,
  nrow    = 4,
  labels  = c("", "", "F", "G"),
  heights = c(0.8, 0.8, 1, 1)
)

assembled <- ggarrange(
  left_col,
  right_col,
  ncol   = 2,
  widths = c(2.0, 1)
)

#### Save ####
ggsave("results/metabolomics/assembled_figure.pdf",
       plot   = assembled,
       width  = 18, height = 18)

cairo_pdf("results/metabolomics/assembled_figure_cairo.pdf",
          width = 18, height = 18, family = "Helvetica")
print(assembled)
dev.off()
