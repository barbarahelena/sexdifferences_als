# CAZymes in ALS mouse experiments
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Load libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rio)
library(vegan)
library(ggrepel)
library(lme4)
library(lmerTest)

# Theme
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
                legend.key.size= unit(0.7, "cm"),
                legend.spacing  = unit(0, "cm"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.7), face = "italic")
        ))
}

# CLR (centred log-ratio) transformation: log(x / geometric_mean(x)) per sample
clr_transform <- function(mat, pseudocount = 0.5) {
  mat <- as.matrix(mat) + pseudocount
  log_mat <- log(mat)
  sweep(log_mat, 1, rowMeans(log_mat), "-")
}

res_dir <- "results/microbiome/cazymes"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

## Load and prepare data
cayman <- rio::import("data/microbiome/families_cpm_table.tsv")
meta <- readRDS("data/microbiome/meta_microbiome.RDS")
stat <- rio::import("data/microbiome/sample_statistics.tsv")

head(cayman)[1:5,1:5]
rownames(cayman) <- cayman$family
cayman$family <- NULL
cayman[is.na(cayman)] <- 0
cayman <- log10(cayman + 1)
cayman <- as.data.frame(t(as.matrix(cayman)))
names(cayman)
dim(cayman)
head(cayman)[1:5,1:5]
cayman$sampleID <- rownames(cayman)
fam <- ncol(cayman)

# Filter out low-count samples (< 1M reads, as in script 2_1)
lowcountids <- meta$ID[meta$LowCount == TRUE]
message("Removing ", length(lowcountids), " low-count samples: ", paste(lowcountids, collapse = ", "))

# Join with metadata, filter low counts and restrict to weeks 6-18
gene_families <- setdiff(names(cayman), "sampleID")
df_tot <- cayman |>
  left_join(meta, by = c("sampleID" = "ID")) |>
  left_join(stat |> select(sample, passed_reads), by = c("sampleID" = "sample")) |>
  filter(!is.na(Genotype)) |>
  filter(!sampleID %in% lowcountids) |>
  filter(Age_ints >= 6 & Age_ints <= 18)
message(nrow(df_tot), " samples after filtering")

# Set reference levels so LMM coefficients are named GenotypeTDP43 and SexMale
df_tot$Genotype <- relevel(factor(df_tot$Genotype), ref = "WT")
df_tot$Sex      <- relevel(factor(df_tot$Sex),      ref = "Female")

# Check sample counts per group
print(table(interaction(df_tot$Sex, df_tot$Genotype), df_tot$Age_ints))

# Define timepoints
timepoints <- sort(unique(df_tot$Age_ints))

# Timepoints with sufficient TDP43 samples (>=3 per sex)
tp_sufficient_tdp <- timepoints[sapply(timepoints, function(tp) {
  d <- df_tot |> filter(Genotype == "TDP43" & Age_ints == tp)
  nrow(d) > 0 && all(table(d$Sex) >= 3)
})]

# ============================================================================
# 1a. BRAY-CURTIS PCoA PLOTS PER TIMEPOINT - Genotype (TDP43 vs WT)
# ============================================================================

## --- PCoA on all mice, faceted by timepoint ---
caz_mat_all <- as.matrix(df_tot[, gene_families])
rownames(caz_mat_all) <- df_tot$sampleID

bray_all <- vegan::vegdist(caz_mat_all, method = "bray")
pcoord_all <- ape::pcoa(bray_all, correction = "cailliez")
expl_variance_all <- pcoord_all$values$Rel_corr_eig * 100

dbray_all <- as.data.frame(pcoord_all$vectors[, c("Axis.1", "Axis.2")])
dbray_all$sampleID <- rownames(dbray_all)
dbray_all <- left_join(dbray_all, df_tot |> select(sampleID, Sex, Genotype, Age_weeks, Age_ints), by = "sampleID")

# PERMANOVA per timepoint: Genotype effect
bray_per_tp_genotype <- function(tp) {
  set.seed(14)
  dfsel <- df_tot |> filter(Age_ints == tp)
  mat <- as.matrix(dfsel[, gene_families])
  rownames(mat) <- dfsel$sampleID
  bray <- vegan::vegdist(mat, method = "bray")
  res <- adonis2(bray ~ Genotype, data = dfsel)
  return(res$`Pr(>F)`[1])
}

res_geno <- data.frame(
  pvalue = sapply(timepoints, bray_per_tp_genotype),
  Age_ints = timepoints
) |> left_join(df_tot |> distinct(Age_ints, Age_weeks), by = "Age_ints")
res_geno$Age_weeks <- factor(res_geno$Age_weeks, levels = sort(unique(res_geno$Age_weeks)))

dbray_all$Age_weeks <- factor(dbray_all$Age_weeks, levels = levels(res_geno$Age_weeks))
dbray_all_filt <- dbray_all |> filter(Age_ints %in% timepoints)

(braycurt_geno <- ggplot(dbray_all_filt, aes(Axis.1, Axis.2)) +
    geom_point(aes(color = Genotype, shape = Sex), size = 2, alpha = 0.7) +
    xlab(paste0("PCo1 (", round(expl_variance_all[1], 1), "%)")) +
    ylab(paste0("PCo2 (", round(expl_variance_all[2], 1), "%)")) +
    scale_color_manual(values = pal_nejm()(6)[c(3, 6)]) +
    scale_fill_manual(values = pal_nejm()(6)[c(3, 6)], guide = "none") +
    theme_Publication() +
    labs(color = "", fill = "", title = "CAZyme PCoA Bray-Curtis: Genotype") +
    stat_ellipse(geom = "polygon", aes(color = Genotype, fill = Genotype), type = "norm",
                 alpha = 0.1, linewidth = 1.0) +
    theme(legend.position = "top") +
    geom_text(data = res_geno, aes(x = Inf, y = Inf, label = paste0("p = ", round(pvalue, 3))),
              hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
    facet_wrap(~ Age_weeks))
ggsave(file.path(res_dir, "PCoA_BrayCurtis_genotype.pdf"), braycurt_geno, width = 14, height = 8)

# ============================================================================
# Line plots for GH78 and GH106
# ============================================================================

# Filter CAZyme families: present in >=30% of samples AND mean CPM > 20
prevalence <- colMeans(df_tot[, gene_families] > 0)
gene_families <- gene_families[prevalence >= 0.30]
message(length(gene_families), " CAZyme families retained after filtering")

# Precompute CLR-transformed matrix once (samples x gene_families)
clr_mat <- clr_transform(df_tot[, gene_families])
rownames(clr_mat) <- df_tot$sampleID

## --- Line plots over time for GH78 and GH106: LMM with Age as factor + emmeans ---
# Age as a factor: no linearity assumption. Age_fac is coded with the earliest timepoint
# (week 6) as the reference level.
df_tot$Age_fac <- factor(df_tot$Age_ints)

# Difference-in-differences vs baseline (week 6): is the between-group difference at age X
# different from the between-group difference at baseline? Tested by refitting the model on
# just that pair of timepoints (baseline + age X) and reading off the group:age interaction
# coefficient directly - under this treatment coding, that coefficient *is* the DiD estimate.
# BH-adjusted across the non-baseline ages.
pairwise_contrast_vs_baseline <- function(data, response, group_var, baseline_age, ages, stratum = "All") {
  rows <- lapply(ages, function(age) {
    d2 <- data |> filter(Age_ints %in% c(baseline_age, age))
    d2$Age_fac <- relevel(droplevels(factor(d2$Age_ints)), ref = as.character(baseline_age))
    f <- as.formula(paste(response, "~", group_var, "* Age_fac + scale(passed_reads) + (1 | MouseID)"))
    m <- tryCatch(lmerTest::lmer(f, data = d2, REML = FALSE), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    coefs <- summary(m)$coefficients
    int_row <- grep(paste0("^", group_var, ".*:Age_fac"), rownames(coefs), value = TRUE)
    if (length(int_row) == 0) return(NULL)
    data.frame(Age_ints = age,
               estimate = coefs[int_row[1], "Estimate"],
               p.value  = coefs[int_row[1], "Pr(>|t|)"])
  })
  df <- bind_rows(rows)
  df$q.value <- p.adjust(df$p.value, method = "fdr")
  df$grouping_variable <- group_var
  df$Genotype_stratum  <- stratum
  df$sig <- case_when(
    df$q.value < 0.001 ~ "***",
    df$q.value < 0.01  ~ "**",
    df$q.value < 0.05  ~ "*",
    TRUE ~ "ns"
  )
  df
}

line_plist_geno <- list()
line_plist_sex  <- list()
lmm_anova_results    <- list()
lmm_contrast_results <- list()
gene_families2 <- c("GH78", "GH106")

for (i in seq_along(gene_families2)) {
  gf <- gene_families2[i]
  df_tot$cazy_val <- clr_mat[df_tot$sampleID, gf]
  df_lmm <- df_tot |> filter(Age_ints %in% tp_sufficient_tdp)

  # ---- Genotype model: Genotype * Age_fac + (1 | MouseID) ----
  # scale(passed_reads) adjusts for residual sequencing-depth differences between samples.
  means_geno <- df_lmm |>
    group_by(Genotype, Age_ints) |>
    summarise(mean = mean(cazy_val, na.rm = TRUE),
              se = sd(cazy_val, na.rm = TRUE) / sqrt(n()), .groups = "drop")

  lmm_geno <- tryCatch(
    lmerTest::lmer(cazy_val ~ Genotype * Age_fac + scale(passed_reads) + (1 | MouseID),
                   data = df_lmm, REML = FALSE),
    error = function(e) NULL
  )
  if (!is.null(lmm_geno)) {
    anova_geno <- anova(lmm_geno)
    contr_geno <- pairwise_contrast_vs_baseline(df_lmm, "cazy_val", "Genotype",
                                                 baseline_age = min(timepoints),
                                                 ages = setdiff(tp_sufficient_tdp, min(timepoints)))
    y_fixed_geno <- max(means_geno$mean + means_geno$se, na.rm = TRUE) + 0.15
    contr_geno$y_pos <- y_fixed_geno
    lmm_anova_results[[paste0(gf, "_genotype")]] <- data.frame(
      anova_geno, term = rownames(anova_geno), model = paste0(gf, "_genotype")
    )
    lmm_contrast_results[[paste0(gf, "_genotype")]] <- data.frame(
      contr_geno, model = paste0(gf, "_genotype")
    )
  } else {
    contr_geno <- NULL
  }

  pl_geno <- ggplot(means_geno, aes(x = Age_ints, y = mean, group = Genotype, color = Genotype)) +
    geom_line() +
    geom_point(size = 1) +
    geom_jitter(data = df_lmm, aes(x = Age_ints, y = cazy_val, color = Genotype),
                width = 0.2, alpha = 0.5, inherit.aes = FALSE) +
    scale_color_manual(values = pal_nejm()(8)[c(6, 3)]) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3) +
    labs(title = gf, x = "Age (weeks)", y = "CLR(CPM)") +
    scale_x_continuous(breaks = timepoints) +
    theme_Publication()
  if (!is.null(contr_geno)) {
    pl_geno <- pl_geno +
      geom_text(data = contr_geno |> filter(sig != "ns"), aes(x = Age_ints, y = y_pos, label = sig),
                inherit.aes = FALSE, size = 6, vjust = 0)
  }
  line_plist_geno[[length(line_plist_geno) + 1]] <- pl_geno

  # ---- Sex model within TDP43: Sex * Age_fac + (1 | MouseID) ----
  # scale(passed_reads) adjusts for residual sequencing-depth differences between samples.
  df_tdp <- df_lmm |> filter(Genotype == "TDP43")
  means_sex <- df_tdp |>
    group_by(Sex, Age_ints) |>
    summarise(mean = mean(cazy_val, na.rm = TRUE),
              se = sd(cazy_val, na.rm = TRUE) / sqrt(n()), .groups = "drop")

  lmm_sex <- tryCatch(
    lmerTest::lmer(cazy_val ~ Sex * Age_fac + scale(passed_reads) + (1 | MouseID),
                   data = df_tdp, REML = FALSE),
    error = function(e) NULL
  )
  if (!is.null(lmm_sex)) {
    anova_sex <- anova(lmm_sex)
    # Week 8 excluded from testing (only 2 TDP43 males at that age); still shown in the plot.
    contr_sex <- pairwise_contrast_vs_baseline(df_tdp, "cazy_val", "Sex",
                                                baseline_age = min(timepoints),
                                                ages = setdiff(tp_sufficient_tdp, c(min(timepoints), 8)),
                                                stratum = "TDP43")
    y_fixed_sex <- max(means_sex$mean + means_sex$se, na.rm = TRUE) + 0.02
    contr_sex$y_pos <- y_fixed_sex
    lmm_anova_results[[paste0(gf, "_sex_TDP43")]] <- data.frame(
      anova_sex, term = rownames(anova_sex), model = paste0(gf, "_sex_TDP43")
    )
    lmm_contrast_results[[paste0(gf, "_sex_TDP43")]] <- data.frame(
      contr_sex, model = paste0(gf, "_sex_TDP43")
    )
  } else {
    contr_sex <- NULL
  }

  pl_sex <- ggplot(means_sex, aes(x = Age_ints, y = mean, group = Sex, color = Sex)) +
    geom_line() +
    geom_point(size = 1) +
    geom_jitter(data = df_tdp, aes(x = Age_ints, y = cazy_val, color = Sex),
                width = 0.2, alpha = 0.5, inherit.aes = FALSE) +
    scale_color_nejm() +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3) +
    labs(title = gf, x = "Age (weeks)", y = "CLR(CPM)") +
    scale_x_continuous(breaks = timepoints) +
    theme_Publication()
  if (!is.null(contr_sex)) {
    pl_sex <- pl_sex +
      geom_text(data = contr_sex |> filter(sig != "ns"), aes(x = Age_ints, y = y_pos, label = sig),
                inherit.aes = FALSE, size = 6, vjust = 0)
  }
  line_plist_sex[[length(line_plist_sex) + 1]] <- pl_sex
}

# Genotype: two panels (A = GH78, B = GH106), shared legend
line_plots_geno <- ggarrange(plotlist = line_plist_geno, ncol = 2, nrow = 1,
                              labels = c("A", "B"), common.legend = TRUE, legend = "bottom")
line_plots_geno <- annotate_figure(line_plots_geno,
                                    top = text_grob("TDP43 vs WT", face = "bold", size = 14))
ggsave(file.path(res_dir, "GH78_GH106_lineplots_genotype.pdf"), line_plots_geno, width = 12, height = 6)

# Sex: two panels (A = GH78, B = GH106), shared legend
line_plots_sex <- ggarrange(plotlist = line_plist_sex, ncol = 2, nrow = 1,
                             labels = c("A", "B"), common.legend = TRUE, legend = "bottom")
line_plots_sex <- annotate_figure(line_plots_sex,
                                   top = text_grob("TDP43: males vs females", face = "bold", size = 14))
ggsave(file.path(res_dir, "GH78_GH106_lineplots_sex.pdf"), line_plots_sex, width = 12, height = 6)

# Save LMM results: F-tests and per-timepoint contrasts separately
write.csv2(bind_rows(lmm_anova_results),    file.path(res_dir, "lmm_GH78_GH106_anova.csv"),     row.names = FALSE)
write.csv2(bind_rows(lmm_contrast_results), file.path(res_dir, "lmm_GH78_GH106_contrasts.csv"), row.names = FALSE)

## --- LMM: GH78 ~ Genotype * Sex + Age_fac + (1 | MouseID) ---
df_tot$GH78_log <- clr_mat[df_tot$sampleID, "GH78"]
model_gh78 <- lmerTest::lmer(GH78_log ~ Genotype * Sex + Age_fac + (1 | MouseID), data = df_tot, REML = FALSE)
summary(model_gh78)
confint(model_gh78)

## --- LMM: GH106 ~ Genotype * Sex + Age_fac + (1 | MouseID) ---
df_tot$GH106_log <- clr_mat[df_tot$sampleID, "GH106"]
model_gh106 <- lmerTest::lmer(GH106_log ~ Genotype * Sex + Age_fac + (1 | MouseID), data = df_tot, REML = FALSE)
summary(model_gh106)
confint(model_gh106)

# ============================================================================
# 3a. VOLCANO PLOTS PER TIMEPOINT: Genotype (TDP43 vs WT)
# ============================================================================

volcano_geno_plist <- list()
volcano_geno_results <- list()

for (tp in timepoints) {
  df_tp <- df_tot |> filter(Age_ints == tp)
  tp_label <- unique(df_tp$Age_weeks)

  statres <- data.frame()
  for (gf in gene_families) {
    df_tp$mb <- clr_mat[df_tp$sampleID, gf]
    tryCatch({
      model <- lm(mb ~ Genotype, data = df_tp)
      res <- summary(model)
      if (nrow(res$coefficients) >= 2) {
        ci <- confint(model)
        statres <- rbind(statres, data.frame(
          gene_family = gf,
          pval = res$coefficients[2, 4],
          estimate = res$coefficients[2, 1],
          se = res$coefficients[2, 2],
          conflow = ci[2, 1],
          confhigh = ci[2, 2]
        ))
      }
    }, error = function(e) NULL)
  }

  statres <- statres |>
    arrange(pval) |>
    mutate(padj = p.adjust(pval, method = "fdr"),
           timepoint = tp_label)

  volcano_geno_results[[as.character(tp)]] <- statres

  statres_sig <- statres |>
    mutate(sig_level = case_when(
      padj < 0.05 & abs(estimate) > 0.5 ~ "Significant & Large Effect",
      padj < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"))

  pl <- ggplot(statres_sig, aes(x = estimate, y = -log10(pval), color = sig_level)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_text_repel(data = statres_sig |> filter(padj < 0.05),
                    aes(label = gene_family), size = 3, max.overlaps = 20,
                    box.padding = 0.5, point.padding = 0.3, color = "black") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Significant & Large Effect" = "red",
                                   "Significant" = "orange",
                                   "Not Significant" = "gray70")) +
    theme_Publication() +
    labs(x = "Estimate (TDP43 vs WT)",
         y = "-log10(p-value)",
         title = paste0("CAZymes: TDP43 vs WT - ", tp_label), color = "")

  volcano_geno_plist[[as.character(tp)]] <- pl
}

# Save individual volcano plots
for (tp in names(volcano_geno_plist)) {
  ggsave(file.path(res_dir, paste0("volcano_genotype_", tp, "wk.pdf")),
         volcano_geno_plist[[tp]], width = 7, height = 7)
}

# Combined volcano plot
n_volc_g <- length(volcano_geno_plist)
n_cols_g <- min(3, n_volc_g)
n_rows_g <- ceiling(n_volc_g / n_cols_g)
(volcano_geno_combined <- ggarrange(plotlist = volcano_geno_plist, common.legend = TRUE, legend = "bottom",
                                     labels = LETTERS[1:n_volc_g],
                                     nrow = n_rows_g, ncol = n_cols_g))
ggsave(file.path(res_dir, "volcano_genotype_all_timepoints.pdf"), volcano_geno_combined,
       width = 7 * n_cols_g, height = 7 * n_rows_g)

# Save volcano stats
volcano_geno_all <- bind_rows(volcano_geno_results)
write.csv2(volcano_geno_all, file.path(res_dir, "lm_cazyme_genotype_per_timepoint.csv"), row.names = FALSE)

# ============================================================================
# 3b. VOLCANO PLOTS PER TIMEPOINT: Sex differences within TDP43
# ============================================================================

volcano_plist <- list()
volcano_results <- list()

for (tp in tp_sufficient_tdp) {
  df_tp <- df_tot |> filter(Genotype == "TDP43" & Age_ints == tp)
  tp_label <- unique(df_tp$Age_weeks)

  statres <- data.frame()
  for (gf in gene_families) {
    df_tp$mb <- clr_mat[df_tp$sampleID, gf]
    tryCatch({
      model <- lm(mb ~ Sex, data = df_tp)
      res <- summary(model)
      if (nrow(res$coefficients) >= 2) {
        ci <- confint(model)
        statres <- rbind(statres, data.frame(
          gene_family = gf,
          pval = res$coefficients[2, 4],
          estimate = res$coefficients[2, 1],
          se = res$coefficients[2, 2],
          conflow = ci[2, 1],
          confhigh = ci[2, 2]
        ))
      }
    }, error = function(e) NULL)
  }

  statres <- statres |>
    arrange(pval) |>
    mutate(padj = p.adjust(pval, method = "fdr"),
           timepoint = tp_label)

  volcano_results[[as.character(tp)]] <- statres

  statres_sig <- statres |>
    mutate(sig_level = case_when(
      padj < 0.05 & abs(estimate) > 0.5 ~ "Significant & Large Effect",
      padj < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"))

  pl <- ggplot(statres_sig, aes(x = estimate, y = -log10(pval), color = sig_level)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_text_repel(data = statres_sig |> filter(padj < 0.05),
                    aes(label = gene_family), size = 3, max.overlaps = 20,
                    box.padding = 0.5, point.padding = 0.3, color = "black") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Significant & Large Effect" = "red",
                                   "Significant" = "orange",
                                   "Not Significant" = "gray70")) +
    theme_Publication() +
    labs(x = "Estimate (Male vs Female)",
         y = "-log10(p-value)",
         title = paste0("CAZymes: Sex diff in TDP43 - ", tp_label), color = "")

  volcano_plist[[as.character(tp)]] <- pl
}

# Save individual volcano plots
for (tp in names(volcano_plist)) {
  ggsave(file.path(res_dir, paste0("volcano_sexdiff_TDP43_", tp, "wk.pdf")),
         volcano_plist[[tp]], width = 7, height = 7)
}

# Combined volcano plot
n_volc <- length(volcano_plist)
n_cols_v <- min(3, n_volc)
n_rows_v <- ceiling(n_volc / n_cols_v)
(volcano_combined <- ggarrange(plotlist = volcano_plist, common.legend = TRUE, legend = "bottom",
                                labels = LETTERS[1:n_volc],
                                nrow = n_rows_v, ncol = n_cols_v))
ggsave(file.path(res_dir, "volcano_sexdiff_TDP43_all_timepoints.pdf"), volcano_combined,
       width = 7 * n_cols_v, height = 7 * n_rows_v)

# Save volcano stats
volcano_all <- bind_rows(volcano_results)
write.csv2(volcano_all, file.path(res_dir, "lm_cazyme_sexdiff_TDP43_per_timepoint.csv"), row.names = FALSE)

# ============================================================================
# 4a. BOXPLOTS FOR SIGNIFICANT FAMILIES FROM GENOTYPE VOLCANO (per timepoint)
# ============================================================================

sig_geno <- volcano_geno_all |> filter(padj < 0.05)
if (nrow(sig_geno) > 0) {
  sig_geno_split <- split(sig_geno, sig_geno$timepoint)
  for (tp_label in names(sig_geno_split)) {
    tp_int <- as.integer(str_extract(tp_label, "\\d+"))
    sig_fams <- sig_geno_split[[tp_label]]$gene_family
    df_tp <- df_tot |> filter(Age_ints == tp_int)

    bp_list <- list()
    for (gf in sig_fams) {
      df_tp$cazy_val <- clr_mat[df_tp$sampleID, gf]

      p_fem <- tryCatch(
        format.pval(wilcox.test(cazy_val ~ Genotype, data = df_tp |> filter(Sex == "Female"))$p.value, digits = 2),
        error = function(e) "NA")
      p_mal <- tryCatch(
        format.pval(wilcox.test(cazy_val ~ Genotype, data = df_tp |> filter(Sex == "Male"))$p.value, digits = 2),
        error = function(e) "NA")
      caption_text <- paste0("Genotype: p = ", p_fem, " (F); p = ", p_mal, " (M)")

      pl <- ggplot(df_tp, aes(x = Sex, y = cazy_val)) +
        geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.4) +
        geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
        labs(title = gf, caption = caption_text, y = "CLR(CPM)", x = "") +
        facet_wrap(~ Genotype) +
        theme_Publication() +
        theme(plot.caption = element_text(size = 8), legend.position = "none")

      bp_list[[length(bp_list) + 1]] <- pl
    }

    # Split into pages of max 9 plots each
    max_per_page <- 9
    n_pages <- ceiling(length(bp_list) / max_per_page)
    for (pg in seq_len(n_pages)) {
      idx <- ((pg - 1) * max_per_page + 1):min(pg * max_per_page, length(bp_list))
      bp_page <- bp_list[idx]
      n_bp <- length(bp_page)
      n_cols_bp <- min(3, n_bp)
      n_rows_bp <- ceiling(n_bp / n_cols_bp)
      suffix <- if (n_pages > 1) paste0("_p", pg) else ""
      (bp_combined <- ggarrange(plotlist = bp_page, labels = LETTERS[idx],
                                 nrow = n_rows_bp, ncol = n_cols_bp))
      ggsave(file.path(res_dir, paste0("boxplots_sig_genotype_", tp_int, "wk", suffix, ".pdf")),
             bp_combined, width = 4 * n_cols_bp, height = 4 * n_rows_bp)
    }
  }
}

# ============================================================================
# 4b. BOXPLOTS FOR SIGNIFICANT FAMILIES FROM SEX DIFF VOLCANO (per timepoint)
# ============================================================================

sig_sex <- volcano_all |> filter(padj < 0.05)
if (nrow(sig_sex) > 0) {
  sig_sex_split <- split(sig_sex, sig_sex$timepoint)
  for (tp_label in names(sig_sex_split)) {
    tp_int <- as.integer(str_extract(tp_label, "\\d+"))
    sig_fams <- sig_sex_split[[tp_label]]$gene_family
    df_tp <- df_tot |> filter(Age_ints == tp_int)

    bp_list <- list()
    for (gf in sig_fams) {
      df_tp$cazy_val <- clr_mat[df_tp$sampleID, gf]

      p_fem <- tryCatch(
        format.pval(wilcox.test(cazy_val ~ Genotype, data = df_tp |> filter(Sex == "Female"))$p.value, digits = 2),
        error = function(e) "NA")
      p_mal <- tryCatch(
        format.pval(wilcox.test(cazy_val ~ Genotype, data = df_tp |> filter(Sex == "Male"))$p.value, digits = 2),
        error = function(e) "NA")
      caption_text <- paste0("Genotype: p = ", p_fem, " (F); p = ", p_mal, " (M)")

      pl <- ggplot(df_tp, aes(x = Sex, y = cazy_val)) +
        geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.4) +
        geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
        labs(title = gf, caption = caption_text, y = "CLR(CPM)", x = "") +
        facet_wrap(~ Genotype) +
        theme_Publication() +
        theme(plot.caption = element_text(size = 8), legend.position = "none")

      bp_list[[length(bp_list) + 1]] <- pl
    }

    # Split into pages of max 9 plots each
    max_per_page <- 9
    n_pages <- ceiling(length(bp_list) / max_per_page)
    for (pg in seq_len(n_pages)) {
      idx <- ((pg - 1) * max_per_page + 1):min(pg * max_per_page, length(bp_list))
      bp_page <- bp_list[idx]
      n_bp <- length(bp_page)
      n_cols_bp <- min(3, n_bp)
      n_rows_bp <- ceiling(n_bp / n_cols_bp)
      suffix <- if (n_pages > 1) paste0("_p", pg) else ""
      (bp_combined <- ggarrange(plotlist = bp_page,
                                 nrow = n_rows_bp, ncol = n_cols_bp))
      ggsave(file.path(res_dir, paste0("boxplots_sig_sexdiff_TDP43_", tp_int, "wk", suffix, ".pdf")),
             bp_combined, width = 4 * n_cols_bp, height = 4 * n_rows_bp)
    }
  }
}

