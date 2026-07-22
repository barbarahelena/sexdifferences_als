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
source("scripts/1_5_diffabundance.R") # -> p_li_als_sex, p_li_ctrl_sex (microbiome sex-diff volcanoes)
p_li_als_sex_mb  <- p_li_als_sex
p_li_ctrl_sex_mb <- p_li_ctrl_sex
source("scripts/1_6_pathwaycorr.R")   # -> heatmap_top, lgd_sig_path
source("scripts/1_8_cazymes.R")       # -> pl_caz1, pl_caz2, p_li_als_sex, p_li_ctrl_sex (CAZyme sex-diff volcanoes; overwrites the microbiome versions above)

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

#### Assemble figure (portrait) ####
# Row 1: A) population table | B+C) Bray-Curtis PCoA ALS + Control
# Row 2: D+E) LInDA volcanos (microbiome) ALS + Control (stacked) | F) pathway heatmap
# Row 3: G+H) CAZy PCoA ALS + Control | I+J) LInDA volcanos (CAZymes) ALS + Control

# Microbiome PCoA: ALS and Control side by side
pcoa_mb <- ggarrange(
  pl1, pl2,
  ncol          = 2,
  labels        = LETTERS[2:3],
  common.legend = TRUE,
  legend        = "bottom"
)

# LInDA volcanos (microbiome): ALS and Control stacked
linda_volc <- ggarrange(
  p_li_als_sex_mb,
  p_li_ctrl_sex_mb,
  nrow   = 2,
  heights = c(1.2, 1.0),
  labels = LETTERS[4:5]
)

# LInDA volcanos (CAZymes): ALS and Control side by side
caz_volc <- ggarrange(
  p_li_als_sex,
  p_li_ctrl_sex,
  ncol   = 2,
  labels = LETTERS[9:10]
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
  NULL,
  heatmap_gg,
  ncol   = 3,
  labels = c("", "", "F"),
  widths = c(1.5, 0.15, 2.25)
)

row3 <- ggarrange(
  pcoa_caz,
  caz_volc,
  ncol   = 2,
  labels = c("", ""),
  widths = c(2, 2)
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

