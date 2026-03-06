# Assemble plot human cohort 2
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

#### Libraries ####
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rio)
library(tableone)
library(gt)
library(ComplexHeatmap)
library(Cairo)
library(vegan)
library(MicrobiomeStat)
library(ape)
library(ggthemes)

#### Source scripts ####
# Each script is sourced to load its plot objects; all paths are relative to project root.
# Note: this will re-run all analyses and re-save intermediate files.
source("scripts/1_2_tableone.R")      # -> table_one_csv, tab_gt, meta
source("scripts/1_4_diversity.R")     # -> pl1 (ALS PCoA), pl2 (Control PCoA)
source("scripts/1_5_diffabundance.R") # -> bar_q_li_als_sex (list: $plot, $n)
source("scripts/1_6_pathwaycorr.R")   # -> heatmap_top, lgd_sig_path
source("scripts/1_8_cazymes.R")       # -> pl_caz1, pl_caz2, plist, families_present

#### Population table as ggplot ####
tab_df <- as.data.frame(table_one_csv) |>
  tibble::rownames_to_column("Variable")

# Replace NA / NaN cells with "-"
tab_df <- tab_df |>
  mutate(across(everything(), ~ ifelse(str_detect(., "(NA|NaN)"), "-", .)))

tab_gg <- ggtexttable(
  tab_df, rows = NULL,
  theme = ttheme(
    "mBlue",
    base_size   = 24,
    padding     = unit(c(3, 3), "mm"),
    colnames.style = colnames_style(fill = "#f0f0f0", face = "bold", size = 18, hjust = 0, x = 0.05),
    tbody.style    = tbody_style(fill = c("white", "#f9f9f9"), size = 18, hjust = 0, x = 0.05)
  )
) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))

#### Pathway heatmap as ggplot ####
heatmap_gg <- as_ggplot(grid.grabExpr(
  draw(
    heatmap_top,
    annotation_legend_list = list(lgd_sig_path),
    padding = unit(c(5, 30, 5, 10), "mm")
  )
))

#### GH78 and GH106 boxplots ####
gh78_idx  <- which(families_present == "GH78")
gh106_idx <- which(families_present == "GH106")
pl_gh78  <- plist[[gh78_idx]]
pl_gh106 <- plist[[gh106_idx]]

#### Assemble figure (portrait) ####
# Row 1: A) population table | B+C) Bray-Curtis PCoA ALS + Control
# Row 2: D+E) LInDA volcanos ALS + Control (stacked) | F) pathway heatmap
# Row 3: G+H) CAZy PCoA ALS + Control | I) GH78 boxplot | J) GH106 boxplot

# Microbiome PCoA: ALS and Control side by side
pcoa_mb <- ggarrange(
  pl1, pl2,
  ncol          = 2,
  labels        = LETTERS[2:3],
  common.legend = TRUE,
  legend        = "bottom"
)

# LInDA volcanos: ALS and Control stacked
linda_volc <- ggarrange(
  p_li_als_sex,
  p_li_ctrl_sex,
  nrow   = 2,
  heights = c(1.2, 1.0),
  labels = LETTERS[4:5]
)

# CAZy PCoA: ALS and Control side by side
pcoa_caz <- ggarrange(
  pl_caz1, pl_caz2,
  ncol          = 2,
  labels        = LETTERS[7:8],
  common.legend = TRUE,
  legend        = "bottom"
)

row1 <- ggarrange(
  ggarrange(tab_gg, NULL, nrow = 2, heights = c(1, 0.2)),
  pcoa_mb,
  ncol   = 2,
  labels = c("A", ""),
  widths = c(1.7, 2.0)
)

row2 <- ggarrange(
  linda_volc,
  heatmap_gg,
  ncol   = 2,
  labels = c("", "F"),
  widths = c(1.5, 2.25)
)

row3 <- ggarrange(
  pcoa_caz,
  pl_gh78,
  pl_gh106,
  ncol   = 3,
  labels = c("", "I", "J"),
  widths = c(2, 1, 1)
)

assembled <- ggarrange(
  row1, row2, row3,
  nrow    = 3,
  heights = c(1, 2, 1)
)

#### Save ####
ggsave("results/humancohort2/assembled_figure.pdf",
       plot   = assembled,
       width  = 18, height = 22)

cairo_pdf("results/humancohort2/assembled_figure_cairo.pdf",
          width = 18, height = 22, family = "Helvetica")
print(assembled)
dev.off()

