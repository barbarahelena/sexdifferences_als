# Test separate microbes – LinDA
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggrepel)
library(MicrobiomeStat)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.1), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1.0)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(),
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))

}

# colors
grp_colors <- c(
  "Higher in ALS"   = pal_nejm()(8)[3], "Higher in Control" = pal_nejm()(8)[6],
  "Not significant" = "grey70"
)
sex_colors <- c(
  "Higher in Male"   = pal_nejm()(2)[2], "Higher in Female" = pal_nejm()(2)[1],
  "Not significant"  = "grey70"
)

# plot helpers
make_wx_plotdf <- function(wx_res, pos_label, neg_label) {
  wx_res |>
    mutate(
      sig       = case_when(
        qval < 0.05 ~ "q < 0.05",
        pval < 0.05 ~ "p < 0.05",
        TRUE        ~ "ns"
      ),
      sig       = factor(sig, levels = c("q < 0.05", "p < 0.05", "ns")),
      direction = ifelse(diff > 0, pos_label, neg_label),
      color_grp = ifelse(qval < 0.05, direction, "Not significant"),
      feature   = fct_reorder(feature, diff)
    )
}

wx_volcano <- function(plot_df, title, xlab, col_vals) {
  ggplot(plot_df, aes(x = diff, y = -log10(pval), color = color_grp)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.7, size = 2) +
    geom_text_repel(aes(label = ifelse(qval < 0.05, as.character(feature), NA)),
                    size = 3, max.overlaps = 20) +
    scale_color_manual(values = col_vals, name = NULL) +
    labs(x = xlab, y = expression(-log[10](p)), title = title) +
    theme_Publication()
}

wx_barplot <- function(plot_df, title, xlab, col_vals, q_only = TRUE) {
  plot_df_sig <- if (q_only) filter(plot_df, sig == "q < 0.05") else filter(plot_df, sig != "ns")
  if (nrow(plot_df_sig) == 0) return(list(plot = NULL, n = 0))
  p <- ggplot(plot_df_sig, aes(x = diff, y = feature, fill = color_grp)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_col() +
    scale_fill_manual(values = col_vals, name = NULL) +
    labs(x = xlab, y = NULL, title = title) +
    theme_Publication()
  list(plot = p, n = nrow(plot_df_sig))
}

wx_lollipop <- function(plot_df, title, xlab, col_vals, q_only = FALSE) {
  plot_df_sig <- if (q_only) filter(plot_df, sig == "q < 0.05") else filter(plot_df, sig != "ns")
  if (nrow(plot_df_sig) == 0) return(list(plot = NULL, n = 0))
  p <- ggplot(plot_df_sig, aes(x = diff, y = feature, color = color_grp)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.5) +
    scale_color_manual(values = col_vals, name = NULL) +
    labs(x = xlab, y = NULL, title = title) +
    theme_Publication()
  list(plot = p, n = nrow(plot_df_sig))
}

# Get data
mb_full <- readRDS("data/human_cohort2/microbiome.RDS")
df <- readRDS("data/human_cohort2/metadata.RDS")
common_samples <- intersect(rownames(mb_full), df$ID)
mb_full <- mb_full[common_samples, ]
rownames(df) <- df$ID
df$Group <- relevel(factor(df$Group), ref = "Control")
df$Sex   <- relevel(factor(df$Sex),   ref = "Male")
metadata_subset <- df[common_samples, ]

# ── LinDA ─────────────────────────────────────────────────────────────────────
# LinDA (Linear model for Differential Abundance) from MicrobiomeStat.
# log2FoldChange direction follows factor reference levels set at data loading:
#   Group: ref = "Control" → positive log2FC = Higher in ALS
#   Sex:   ref = "Male"    → positive log2FC = Higher in Female

mb_linda <- t(as.matrix(mb_full))  # taxa × samples

# Helper: prevalence filter + LinDA linear model.
# term: the coefficient to extract (e.g. "GroupALS", "SexFemale").
run_linda <- function(counts, col_data, formula, term) {
  cts  <- counts[, rownames(col_data)]
  keep <- rowMeans(cts > 0) >= 0.1
  cts  <- cts[keep, ]
  res  <- MicrobiomeStat::linda(feature.dat      = cts,
                               meta.dat          = col_data,
                               formula           = formula,
                               feature.dat.type  = "proportion",
                               prev.filter       = 0.3,
                               mean.abund.filter = 0.001,
                               alpha             = 0.05,
                               verbose           = FALSE)
  res$output[[term]] |>
    rownames_to_column("feature") |>
    dplyr::select(feature, diff = log2FoldChange, pval = pvalue, qval = padj)
}

li_main     <- run_linda(mb_linda, metadata_subset,                                          "~ Group + Age", "GroupALS")
li_female   <- run_linda(mb_linda, metadata_subset[metadata_subset$Sex   == "Female",  ], "~ Group + Age", "GroupALS")
li_male     <- run_linda(mb_linda, metadata_subset[metadata_subset$Sex   == "Male",    ], "~ Group + Age", "GroupALS")
li_als_sex  <- run_linda(mb_linda, metadata_subset[metadata_subset$Group == "ALS",     ], "~ Sex + Age",   "SexFemale")
li_ctrl_sex <- run_linda(mb_linda, metadata_subset[metadata_subset$Group == "Control", ], "~ Sex + Age",   "SexFemale")

# Save combined results table
list(main = li_main, female = li_female, male = li_male,
     als_sex = li_als_sex, ctrl_sex = li_ctrl_sex) |>
  imap(function(df, nm) rename_with(df, \(x) paste0(x, "_", nm), -feature)) |>
  purrr::reduce(left_join, by = "feature") |>
  write.csv("results/humancohort2/linda_results.csv", row.names = FALSE)

# Plot data frames
pdf_li_main     <- make_wx_plotdf(li_main,     "Higher in ALS",    "Higher in Control")
pdf_li_female   <- make_wx_plotdf(li_female,   "Higher in ALS",    "Higher in Control")
pdf_li_male     <- make_wx_plotdf(li_male,     "Higher in ALS",    "Higher in Control")
pdf_li_als_sex  <- make_wx_plotdf(li_als_sex,  "Higher in Female", "Higher in Male")
pdf_li_ctrl_sex <- make_wx_plotdf(li_ctrl_sex, "Higher in Female", "Higher in Male")

# Volcano plots
(p_li_main     <- wx_volcano(pdf_li_main,     "ALS vs Control",               "log2FC (ALS / Control)",  grp_colors))
(p_li_female   <- wx_volcano(pdf_li_female,   "ALS vs Control in Females",    "log2FC (ALS / Control)",  grp_colors))
(p_li_male     <- wx_volcano(pdf_li_male,     "ALS vs Control in Males",      "log2FC (ALS / Control)",  grp_colors))
(p_li_als_sex  <- wx_volcano(pdf_li_als_sex,  "Sex differences in ALS",       "log2FC (Female / Male)",  sex_colors))
(p_li_ctrl_sex <- wx_volcano(pdf_li_ctrl_sex, "Sex differences in Controls", "log2FC (Female / Male)",  sex_colors))

ggsave(p_li_main,     filename = "results/humancohort2/linda_main_volcano.pdf",       width = 7, height = 6)
ggsave(p_li_female,   filename = "results/humancohort2/linda_female_volcano.pdf",     width = 7, height = 6)
ggsave(p_li_male,     filename = "results/humancohort2/linda_male_volcano.pdf",       width = 7, height = 6)
ggsave(p_li_als_sex,  filename = "results/humancohort2/linda_als_sex_volcano.pdf",    width = 7, height = 6)
ggsave(p_li_ctrl_sex, filename = "results/humancohort2/linda_ctrl_sex_volcano.pdf",   width = 7, height = 6)

# Lollipop plots – q < 0.05 only (skipped when no hits)
(lollipop_q_li_main     <- wx_lollipop(pdf_li_main,     "ALS vs Control (q < 0.05)",              "log2FC (ALS / Control)",  grp_colors, q_only = TRUE))
(lollipop_q_li_female   <- wx_lollipop(pdf_li_female,   "ALS vs Control in Females (q < 0.05)",   "log2FC (ALS / Control)",  grp_colors, q_only = TRUE))
(lollipop_q_li_male     <- wx_lollipop(pdf_li_male,     "ALS vs Control in Males (q < 0.05)",     "log2FC (ALS / Control)",  grp_colors, q_only = TRUE))
(lollipop_q_li_als_sex  <- wx_lollipop(pdf_li_als_sex,  "Sex differences in ALS (q < 0.05)",      "log2FC (Female / Male)",  sex_colors, q_only = TRUE))
(lollipop_q_li_ctrl_sex <- wx_lollipop(pdf_li_ctrl_sex, "Sex differences in Controls (q < 0.05)", "log2FC (Female / Male)",  sex_colors, q_only = TRUE))

if (lollipop_q_li_main$n     > 0) ggsave(lollipop_q_li_main$plot,     filename = "results/humancohort2/linda_main_lollipop_q05.pdf",       width = 7, height = 0.25 * lollipop_q_li_main$n     + 2)
if (lollipop_q_li_female$n   > 0) ggsave(lollipop_q_li_female$plot,   filename = "results/humancohort2/linda_female_lollipop_q05.pdf",     width = 7, height = 0.25 * lollipop_q_li_female$n   + 2)
if (lollipop_q_li_male$n     > 0) ggsave(lollipop_q_li_male$plot,     filename = "results/humancohort2/linda_male_lollipop_q05.pdf",       width = 7, height = 0.25 * lollipop_q_li_male$n     + 2)
if (lollipop_q_li_als_sex$n  > 0) ggsave(lollipop_q_li_als_sex$plot,  filename = "results/humancohort2/linda_als_sex_lollipop_q05.pdf",    width = 7, height = 0.25 * lollipop_q_li_als_sex$n  + 2)
if (lollipop_q_li_ctrl_sex$n > 0) ggsave(lollipop_q_li_ctrl_sex$plot, filename = "results/humancohort2/linda_ctrl_sex_lollipop_q05.pdf",   width = 7, height = 0.25 * lollipop_q_li_ctrl_sex$n + 2)

# Barplots – q < 0.05 only (skipped when no hits)
(bar_q_li_main     <- wx_barplot(pdf_li_main,     "ALS vs Control (q < 0.05)",              "log2FC (ALS / Control)",  grp_colors))
(bar_q_li_female   <- wx_barplot(pdf_li_female,   "ALS vs Control in Females (q < 0.05)",   "log2FC (ALS / Control)",  grp_colors))
(bar_q_li_male     <- wx_barplot(pdf_li_male,     "ALS vs Control in Males (q < 0.05)",     "log2FC (ALS / Control)",  grp_colors))
(bar_q_li_als_sex  <- wx_barplot(pdf_li_als_sex,  "Sex differences in ALS (q < 0.05)",      "log2FC (Female / Male)",  sex_colors))
(bar_q_li_ctrl_sex <- wx_barplot(pdf_li_ctrl_sex, "Sex differences in Controls (q < 0.05)", "log2FC (Female / Male)",  sex_colors))

if (bar_q_li_main$n     > 0) ggsave(bar_q_li_main$plot,     filename = "results/humancohort2/linda_main_barplot_q05.pdf",       width = 7, height = 0.25 * bar_q_li_main$n     + 2)
if (bar_q_li_female$n   > 0) ggsave(bar_q_li_female$plot,   filename = "results/humancohort2/linda_female_barplot_q05.pdf",     width = 7, height = 0.25 * bar_q_li_female$n   + 2)
if (bar_q_li_male$n     > 0) ggsave(bar_q_li_male$plot,     filename = "results/humancohort2/linda_male_barplot_q05.pdf",       width = 7, height = 0.25 * bar_q_li_male$n     + 2)
if (bar_q_li_als_sex$n  > 0) ggsave(bar_q_li_als_sex$plot,  filename = "results/humancohort2/linda_als_sex_barplot_q05.pdf",    width = 7, height = 0.25 * bar_q_li_als_sex$n  + 2)
if (bar_q_li_ctrl_sex$n > 0) ggsave(bar_q_li_ctrl_sex$plot, filename = "results/humancohort2/linda_ctrl_sex_barplot_q05.pdf",   width = 7, height = 0.25 * bar_q_li_ctrl_sex$n + 2)
