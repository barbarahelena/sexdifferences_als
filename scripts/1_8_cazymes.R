# Cayman analyses human cohort Calgary

# Library
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rio)
library(rstatix)
library(ggrepel)
library(vegan)
library(MicrobiomeStat)

# Plot theme
theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(family = 'Helvetica'),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.9)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text.x = element_text(size = rel(0.9)),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.7, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
} 

res_dir <- "results/humancohort2"

# Data
cayman <- readRDS("data/human_cohort2/cayman_families_merged.RDS")
head(cayman)[1:5,1:5]
rownames(cayman) <- cayman$family
cayman$family <- NULL
cayman <- as.data.frame(t(as.matrix(cayman)))
cayman[is.na(cayman)] <- 0
names(cayman)
dim(cayman)
head(cayman)[1:5,1:5]
cayman$sampleID <- rownames(cayman)
fam <- ncol(cayman)
clinical <- readRDS("data/human_cohort2/metadata.RDS")

# --- CAZyme richness by sex ---
richness <- readRDS("data/human_cohort2/cayman_sample_statistics.RDS")

richness_plot <- richness |>
  dplyr::rename(sampleID = sample) |>
  dplyr::left_join(clinical |> dplyr::select(ID, Sex, Group), by = c("sampleID" = "ID")) |>
  dplyr::filter(!is.na(Group))

# LM interaction test: Sex * Group
model_int_richness <- lm(richness ~ Sex * Group, data = richness_plot)
int_row_richness <- grep(":", rownames(summary(model_int_richness)$coefficients), value = TRUE)
p_int_richness <- if (length(int_row_richness) > 0) summary(model_int_richness)$coefficients[int_row_richness[1], 4] else NA
int_text_richness <- paste0("Sex × Group interaction: p = ", format.pval(p_int_richness, digits = 2))

(pl_richness <- ggplot(richness_plot, aes(x = Sex, y = richness)) +
    geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.5, alpha = 0.9) +
    geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
    geom_pwc(method = "wilcox.test", p.adjust.method = "fdr", label = "p.signif", hide.ns = TRUE) +
    scale_fill_manual(guide = "none", values = pal_nejm()(2)) +
    facet_wrap(~ Group) +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          plot.caption = element_text(size = 14)) +
    labs(y = "CAZyme richness", x = "", title = "Sex differences in CAZyme richness",
         caption = int_text_richness))

ggsave(file.path(res_dir, "cayman_richness_by_sex.pdf"), pl_richness, width = 5, height = 5)

# --- Beta diversity: CAZyme families ---
# Prepare numeric matrix (all families, unfiltered)
caz_mat <- cayman[, setdiff(names(cayman), "sampleID")]
caz_mat <- caz_mat[rownames(caz_mat) %in% clinical$ID, ]
caz_mat <- log10(caz_mat + 1)

# Beta diversity in ALS
dfals <- clinical |> filter(Group == "ALS")
caz_als <- caz_mat[which(rownames(caz_mat) %in% dfals$ID), ]
bray_als <- vegan::vegdist(caz_als, method = 'bray')
set.seed(14)
pcoord_als <- ape::pcoa(bray_als, correction = "cailliez")
str(pcoord_als$values)
expl_variance_als <- pcoord_als$values$Relative_eig * 100
head(expl_variance_als)
dbray_als <- pcoord_als$vectors[, c('Axis.1', 'Axis.2')]
dbray_als <- as.data.frame(dbray_als)
dbray_als$ID <- rownames(dbray_als)
dbray_als <- left_join(dbray_als, clinical, by = 'ID')
set.seed(14)
res_als <- adonis2(bray_als ~ Sex, data = dbray_als)
res_als <- as.data.frame(res_als)

(pl_caz1 <- ggplot(data = dbray_als, aes(Axis.1, Axis.2)) +
        geom_point(aes(color = Sex), size = 3, alpha = 0.7) +
        xlab(paste0('PCo1 (', round(expl_variance_als[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_als[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_nejm()(2)) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        theme_Publication() +
        labs(color = "", fill = "", title = "CAZymes: Patients with ALS") +
        stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
            alpha = 0.1, linewidth = 1.0) +
        theme(legend.position = "top"))

(pl_caz1 <- pl_caz1 + annotate("text", x = Inf, y = Inf,
                    label = paste0("Sex: R\u00b2 = ", round(res_als$R2[1], 3), ", p = ", round(res_als$`Pr(>F)`[1], 3)),
                    hjust = 1.1, vjust = 1.1, size = 5))

# Beta diversity in controls
dfctrl <- clinical |> filter(Group == "Control")
caz_ctrl <- caz_mat[which(rownames(caz_mat) %in% dfctrl$ID), ]
bray_ctrl <- vegan::vegdist(caz_ctrl, method = 'bray')
pcoord_ctrl <- ape::pcoa(bray_ctrl, correction = "cailliez")
str(pcoord_ctrl$values)
expl_variance_ctrl <- pcoord_ctrl$values$Relative_eig * 100
head(expl_variance_ctrl)
dbray_ctrl <- pcoord_ctrl$vectors[, c('Axis.1', 'Axis.2')]
dbray_ctrl <- as.data.frame(dbray_ctrl)
dbray_ctrl$ID <- rownames(dbray_ctrl)
dbray_ctrl <- left_join(dbray_ctrl, clinical, by = 'ID')
res_ctrl <- adonis2(bray_ctrl ~ Sex, data = dbray_ctrl)
res_ctrl <- as.data.frame(res_ctrl)

(pl_caz2 <- ggplot(data = dbray_ctrl, aes(Axis.1, Axis.2)) +
        geom_point(aes(color = Sex), size = 3, alpha = 0.7) +
        xlab(paste0('PCo1 (', round(expl_variance_ctrl[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_ctrl[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_nejm()(2)) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        theme_Publication() +
        labs(color = "", fill = "", title = "CAZymes: Controls") +
        stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
            alpha = 0.1, linewidth = 1.0) +
        theme(legend.position = "top"))

pl_caz2 <- pl_caz2 + annotate("text", x = Inf, y = Inf,
                    label = paste0("Sex: R\u00b2 = ", round(res_ctrl$R2[1], 3), ", p = ", round(res_ctrl$`Pr(>F)`[1], 3)),
                    hjust = 1.1, vjust = 1.1, size = 5)
pl_caz2
pca_combined <- ggarrange(pl_caz1, pl_caz2, nrow = 1, labels = LETTERS[1:2], common.legend = TRUE, legend = "bottom")
ggsave(file.path(res_dir, "cayman_PCoA_BrayCurtis.pdf"), plot = pca_combined, width = 12, height = 6)

# --- CAZy FAMILY ABUNDANCE DIFFERENCES: ALS vs Control (LinDA) ---
# LinDA (Linear model for Differential Abundance) from MicrobiomeStat.
# log2FoldChange direction follows factor reference levels set below:
#   Group: ref = "Control" → positive log2FC = Higher in ALS
#   Sex:   ref = "Male"    → positive log2FC = Higher in Female

grp_colors <- c(
  "Higher in ALS"     = pal_nejm()(8)[3],
  "Higher in Control" = pal_nejm()(8)[6],
  "Not significant"   = "grey70"
)
sex_colors <- c(
  "Higher in Male"    = pal_nejm()(8)[2],
  "Higher in Female"  = pal_nejm()(8)[1],
  "Not significant"   = "grey70"
)

make_wx_plotdf <- function(wx_res, pos_label, neg_label) {
  wx_res |>
    mutate(
      direction = ifelse(diff > 0, pos_label, neg_label),
      color_grp = ifelse(qval < 0.05, direction, "Not significant"),
      sig       = case_when(
        qval < 0.05 ~ "q < 0.05",
        pval < 0.05 ~ "p < 0.05",
        TRUE        ~ "ns"
      ),
      sig       = factor(sig, levels = c("q < 0.05", "p < 0.05", "ns")),
      feature   = fct_reorder(as.character(feature), diff)
    )
}

wx_volcano <- function(plot_df, title, xlab, col_vals) {
  ggplot(plot_df, aes(x = diff, y = -log10(pval), color = color_grp)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.7, size = 2) +
    geom_text_repel(data = subset(plot_df, qval < 0.05), aes(label = as.character(feature)),
                     size = 3, max.overlaps = 20) +
    scale_color_manual(values = col_vals, name = NULL) +
    labs(x = xlab, y = expression(-log[10]("p-value")), title = title) +
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

# Metadata aligned to the CAZyme samples (cayman rows are already keyed by sampleID)
metadata_subset <- clinical
rownames(metadata_subset) <- metadata_subset$ID
common_samples_linda <- intersect(cayman$sampleID, metadata_subset$ID)
metadata_subset <- metadata_subset[common_samples_linda, ]
metadata_subset$Group <- relevel(factor(metadata_subset$Group), ref = "Control")
metadata_subset$Sex   <- relevel(factor(metadata_subset$Sex),   ref = "Male")

# CAZyme families (rows) x samples (columns) matrix for LinDA
mb_linda <- t(as.matrix(cayman[common_samples_linda, setdiff(names(cayman), "sampleID")]))

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
                                mean.abund.filter = 0,
                                alpha             = 0.05,
                                p.adj.method      = "BH",
                                verbose           = TRUE)
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
  write.csv(file.path(res_dir, "cayman_linda_results.csv"), row.names = FALSE)

# Plot data frames
pdf_li_main     <- make_wx_plotdf(li_main,     "Higher in ALS",    "Higher in Control")
pdf_li_female   <- make_wx_plotdf(li_female,   "Higher in ALS",    "Higher in Control")
pdf_li_male     <- make_wx_plotdf(li_male,     "Higher in ALS",    "Higher in Control")
pdf_li_als_sex  <- make_wx_plotdf(li_als_sex,  "Higher in Female", "Higher in Male")
pdf_li_ctrl_sex <- make_wx_plotdf(li_ctrl_sex, "Higher in Female", "Higher in Male")

# Volcano plots
(p_li_main     <- wx_volcano(pdf_li_main,     "CAZymes: ALS vs Control",               "log2FC (ALS / Control)",  grp_colors))
(p_li_female   <- wx_volcano(pdf_li_female,   "CAZymes: ALS vs Control in Females",    "log2FC (ALS / Control)",  grp_colors))
(p_li_male     <- wx_volcano(pdf_li_male,     "CAZymes: ALS vs Control in Males",      "log2FC (ALS / Control)",  grp_colors))
(p_li_als_sex  <- wx_volcano(pdf_li_als_sex,  "CAZymes: Sex differences in ALS",       "log2FC (Female / Male)",  sex_colors))
(p_li_ctrl_sex <- wx_volcano(pdf_li_ctrl_sex, "CAZymes: Sex differences in Controls",  "log2FC (Female / Male)",  sex_colors))

# Harmonize the y-axis (-log10 p) between ALS and Control sex-diff volcanoes so the two panels are comparable
y_max_sex <- max(-log10(pdf_li_als_sex$pval), -log10(pdf_li_ctrl_sex$pval), na.rm = TRUE)
p_li_als_sex  <- p_li_als_sex  + coord_cartesian(ylim = c(0, y_max_sex * 1.05))
p_li_ctrl_sex <- p_li_ctrl_sex + coord_cartesian(ylim = c(0, y_max_sex * 1.05))

ggsave(p_li_main,     filename = file.path(res_dir, "cayman_linda_main_volcano.pdf"),     width = 7, height = 6)
ggsave(p_li_female,   filename = file.path(res_dir, "cayman_linda_female_volcano.pdf"),   width = 7, height = 6)
ggsave(p_li_male,     filename = file.path(res_dir, "cayman_linda_male_volcano.pdf"),     width = 7, height = 6)
ggsave(p_li_als_sex,  filename = file.path(res_dir, "cayman_linda_als_sex_volcano.pdf"),  width = 7, height = 6)
ggsave(p_li_ctrl_sex, filename = file.path(res_dir, "cayman_linda_ctrl_sex_volcano.pdf"), width = 7, height = 6)

# Lollipop plots – q < 0.05 only (skipped when no hits)
(lollipop_q_li_main     <- wx_lollipop(pdf_li_main,     "CAZymes: ALS vs Control (q < 0.05)",              "log2FC (ALS / Control)",  grp_colors, q_only = TRUE))
(lollipop_q_li_female   <- wx_lollipop(pdf_li_female,   "CAZymes: ALS vs Control in Females (q < 0.05)",   "log2FC (ALS / Control)",  grp_colors, q_only = TRUE))
(lollipop_q_li_male     <- wx_lollipop(pdf_li_male,     "CAZymes: ALS vs Control in Males (q < 0.05)",     "log2FC (ALS / Control)",  grp_colors, q_only = TRUE))
(lollipop_q_li_als_sex  <- wx_lollipop(pdf_li_als_sex,  "CAZymes: Sex differences in ALS (q < 0.05)",      "log2FC (Female / Male)",  sex_colors, q_only = TRUE))
(lollipop_q_li_ctrl_sex <- wx_lollipop(pdf_li_ctrl_sex, "CAZymes: Sex differences in Controls (q < 0.05)", "log2FC (Female / Male)",  sex_colors, q_only = TRUE))

if (lollipop_q_li_main$n     > 0) ggsave(lollipop_q_li_main$plot,     filename = file.path(res_dir, "cayman_linda_main_lollipop_q05.pdf"),     width = 7, height = 0.25 * lollipop_q_li_main$n     + 2)
if (lollipop_q_li_female$n   > 0) ggsave(lollipop_q_li_female$plot,   filename = file.path(res_dir, "cayman_linda_female_lollipop_q05.pdf"),   width = 7, height = 0.25 * lollipop_q_li_female$n   + 2)
if (lollipop_q_li_male$n     > 0) ggsave(lollipop_q_li_male$plot,     filename = file.path(res_dir, "cayman_linda_male_lollipop_q05.pdf"),     width = 7, height = 0.25 * lollipop_q_li_male$n     + 2)
if (lollipop_q_li_als_sex$n  > 0) ggsave(lollipop_q_li_als_sex$plot,  filename = file.path(res_dir, "cayman_linda_als_sex_lollipop_q05.pdf"),  width = 7, height = 0.25 * lollipop_q_li_als_sex$n  + 2)
if (lollipop_q_li_ctrl_sex$n > 0) ggsave(lollipop_q_li_ctrl_sex$plot, filename = file.path(res_dir, "cayman_linda_ctrl_sex_lollipop_q05.pdf"), width = 7, height = 0.25 * lollipop_q_li_ctrl_sex$n + 2)

# Barplots – q < 0.05 only (skipped when no hits)
(bar_q_li_main     <- wx_barplot(pdf_li_main,     "CAZymes: ALS vs Control (q < 0.05)",              "log2FC (ALS / Control)",  grp_colors))
(bar_q_li_female   <- wx_barplot(pdf_li_female,   "CAZymes: ALS vs Control in Females (q < 0.05)",   "log2FC (ALS / Control)",  grp_colors))
(bar_q_li_male     <- wx_barplot(pdf_li_male,     "CAZymes: ALS vs Control in Males (q < 0.05)",     "log2FC (ALS / Control)",  grp_colors))
(bar_q_li_als_sex  <- wx_barplot(pdf_li_als_sex,  "CAZymes: Sex differences in ALS (q < 0.05)",      "log2FC (Female / Male)",  sex_colors))
(bar_q_li_ctrl_sex <- wx_barplot(pdf_li_ctrl_sex, "CAZymes: Sex differences in Controls (q < 0.05)", "log2FC (Female / Male)",  sex_colors))

if (bar_q_li_main$n     > 0) ggsave(bar_q_li_main$plot,     filename = file.path(res_dir, "cayman_linda_main_barplot_q05.pdf"),     width = 7, height = 0.25 * bar_q_li_main$n     + 2)
if (bar_q_li_female$n   > 0) ggsave(bar_q_li_female$plot,   filename = file.path(res_dir, "cayman_linda_female_barplot_q05.pdf"),   width = 7, height = 0.25 * bar_q_li_female$n   + 2)
if (bar_q_li_male$n     > 0) ggsave(bar_q_li_male$plot,     filename = file.path(res_dir, "cayman_linda_male_barplot_q05.pdf"),     width = 7, height = 0.25 * bar_q_li_male$n     + 2)
if (bar_q_li_als_sex$n  > 0) ggsave(bar_q_li_als_sex$plot,  filename = file.path(res_dir, "cayman_linda_als_sex_barplot_q05.pdf"),  width = 7, height = 0.25 * bar_q_li_als_sex$n  + 2)
if (bar_q_li_ctrl_sex$n > 0) ggsave(bar_q_li_ctrl_sex$plot, filename = file.path(res_dir, "cayman_linda_ctrl_sex_barplot_q05.pdf"), width = 7, height = 0.25 * bar_q_li_ctrl_sex$n + 2)

# --- Boxplots for LinDA-significant CAZy families (q < 0.05) ---
df_tot <- cayman |>
  left_join(clinical, by = c("sampleID" = "ID")) |>
  filter(!is.na(Group))

# Helper: one boxplot + Wilcoxon p-value per significant family, arranged into a grid.
sig_family_boxplots <- function(features, data, group_var, group_colors) {
  if (length(features) == 0) return(list(plot = NULL, n = 0, nrow = 0, ncol = 0))
  plist <- lapply(features, function(gf) {
    data$cazy_val <- log10(data[[gf]] + 1)
    ggplot(data, aes(x = .data[[group_var]], y = cazy_val)) +
      geom_boxplot(aes(fill = .data[[group_var]]), outlier.shape = NA, width = 0.4) +
      geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
      stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
      scale_fill_manual(values = group_colors, guide = "none") +
      labs(title = gf, y = "log10(cpm + 1)", x = "") +
      theme_Publication() +
      theme(legend.position = "none")
  })
  n    <- length(plist)
  ncol <- min(4, n)
  nrow <- ceiling(n / ncol)
  p    <- ggarrange(plotlist = plist, nrow = nrow, ncol = ncol, labels = LETTERS[seq_len(n)])
  list(plot = p, n = n, nrow = nrow, ncol = ncol)
}

sig_main     <- pdf_li_main     |> filter(sig == "q < 0.05") |> pull(feature) |> as.character()
sig_female   <- pdf_li_female   |> filter(sig == "q < 0.05") |> pull(feature) |> as.character()
sig_male     <- pdf_li_male     |> filter(sig == "q < 0.05") |> pull(feature) |> as.character()
sig_als_sex  <- pdf_li_als_sex  |> filter(sig == "q < 0.05") |> pull(feature) |> as.character()
sig_ctrl_sex <- pdf_li_ctrl_sex |> filter(sig == "q < 0.05") |> pull(feature) |> as.character()

box_li_main     <- sig_family_boxplots(sig_main,     df_tot,                                     "Group", pal_nejm()(2))
box_li_female   <- sig_family_boxplots(sig_female,   df_tot |> filter(Sex == "Female"),           "Group", pal_nejm()(2))
box_li_male     <- sig_family_boxplots(sig_male,     df_tot |> filter(Sex == "Male"),             "Group", pal_nejm()(2))
box_li_als_sex  <- sig_family_boxplots(sig_als_sex,  df_tot |> filter(Group == "ALS"),            "Sex",   pal_nejm()(2))
box_li_ctrl_sex <- sig_family_boxplots(sig_ctrl_sex, df_tot |> filter(Group == "Control"),        "Sex",   pal_nejm()(2))

if (box_li_main$n     > 0) ggsave(box_li_main$plot,     filename = file.path(res_dir, "cayman_linda_main_sig_boxplots.pdf"),     width = 3.5 * box_li_main$ncol,     height = 3.5 * box_li_main$nrow)
if (box_li_female$n   > 0) ggsave(box_li_female$plot,   filename = file.path(res_dir, "cayman_linda_female_sig_boxplots.pdf"),   width = 3.5 * box_li_female$ncol,   height = 3.5 * box_li_female$nrow)
if (box_li_male$n     > 0) ggsave(box_li_male$plot,     filename = file.path(res_dir, "cayman_linda_male_sig_boxplots.pdf"),     width = 3.5 * box_li_male$ncol,     height = 3.5 * box_li_male$nrow)
if (box_li_als_sex$n  > 0) ggsave(box_li_als_sex$plot,  filename = file.path(res_dir, "cayman_linda_als_sex_sig_boxplots.pdf"),  width = 3.5 * box_li_als_sex$ncol,  height = 3.5 * box_li_als_sex$nrow)
if (box_li_ctrl_sex$n > 0) ggsave(box_li_ctrl_sex$plot, filename = file.path(res_dir, "cayman_linda_ctrl_sex_sig_boxplots.pdf"), width = 3.5 * box_li_ctrl_sex$ncol, height = 3.5 * box_li_ctrl_sex$nrow)

# --- Boxplots for CAZy families of interest ---
# Prepare: join cayman abundance with clinical metadata
gene_families <- setdiff(names(cayman), c("sampleID"))

# Filter CAZyme families: present in ≥10% of samples AND mean CPM > 1
prevalence <- colMeans(cayman[, gene_families] > 0)
mean_abundance <- colMeans(cayman[, gene_families])
gene_families <- gene_families[prevalence >= 0.30 & mean_abundance > 20]
message(length(gene_families), " CAZyme families retained after filtering (prevalence ≥10%, mean CPM >1)")

# Keep only families that survived the join as columns
gene_families <- gene_families[gene_families %in% names(df_tot)]

# --- Boxplots for specific CAZy families of interest ---
families_of_interest <- c("GH78", "GH106")
families_present <- families_of_interest[families_of_interest %in% names(df_tot)]

plist <- list()
for (i in seq_along(families_present)) {
  gf <- families_present[i]
  df_tot$cazy_val <- log10(df_tot[[gf]] + 1)

  # LM interaction test: Sex * Group
  model_int <- lm(cazy_val ~ Sex * Group, data = df_tot)
  int_row <- grep(":", rownames(summary(model_int)$coefficients), value = TRUE)
  p_int <- if (length(int_row) > 0) summary(model_int)$coefficients[int_row[1], 4] else NA
  int_text <- paste0("Sex \u00d7 Group interaction: p = ", format.pval(p_int, digits = 2))

  pl <- ggplot(df_tot, aes(x = Sex, y = cazy_val)) +
    geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
    labs(title = gf,
         caption = int_text,
         y = "log10(cpm + 1)",
         x = "") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
          plot.caption = element_text(size = 14),
          legend.position = "none") +
    facet_wrap(~ Group)

  plist[[i]] <- pl
}

(plots_cazy <- ggarrange(plotlist = list(plist[[1]], plist[[2]]), common.legend = TRUE, legend = "bottom",
                         labels = LETTERS[1:2],
                         nrow = 1, ncol = 2))
ggsave(file.path(res_dir, "cayman_families_of_interest_boxplots.pdf"), plots_cazy,
       width = 8, height = 5)
