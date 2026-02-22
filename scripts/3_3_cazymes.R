# CAZymes in ALS mouse experiments
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Load libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rio)
library(vegan)
library(ggrepel)

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

res_dir <- "results/microbiome/cazymes"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

## Load and prepare data
cayman <- rio::import("data/microbiome/families_cpm_table.tsv")
meta <- readRDS("data/microbiome/meta_microbiome.RDS")
stat <- rio::import("data/microbiome/sample_statistics.tsv")

head(cayman)[1:5,1:5]
rownames(cayman) <- cayman$family
cayman$family <- NULL
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
  filter(!is.na(Genotype)) |>
  filter(!sampleID %in% lowcountids) |>
  filter(Age_ints >= 6 & Age_ints <= 18)
message(nrow(df_tot), " samples after filtering")

# Check sample counts per group
print(table(interaction(df_tot$Sex, df_tot$Genotype), df_tot$Age_ints))

# Filter CAZyme families: present in >=30% of samples AND mean CPM > 20
prevalence <- colMeans(df_tot[, gene_families] > 0)
mean_abundance <- colMeans(df_tot[, gene_families])
gene_families <- gene_families[prevalence >= 0.30 & mean_abundance > 20]
message(length(gene_families), " CAZyme families retained after filtering")

# Timepoints for analyses (6, 8, 10, 12, 14, 16, 18 weeks)
timepoints <- c(6, 8, 10, 12, 14, 16, 18)
timepoints <- timepoints[timepoints %in% unique(df_tot$Age_ints)]
message("Timepoints: ", paste(timepoints, collapse = ", "))

# Determine which timepoints have enough samples for sex comparisons (min 3 per sex)
# within each genotype
min_per_group <- 3
tp_sufficient_tdp <- c()
tp_sufficient_wt <- c()
for (tp in timepoints) {
  counts_tdp <- df_tot |> filter(Genotype == "TDP43" & Age_ints == tp) |> count(Sex)
  if (all(counts_tdp$n >= min_per_group)) tp_sufficient_tdp <- c(tp_sufficient_tdp, tp)
  counts_wt <- df_tot |> filter(Genotype == "WT" & Age_ints == tp) |> count(Sex)
  if (all(counts_wt$n >= min_per_group)) tp_sufficient_wt <- c(tp_sufficient_wt, tp)
}
message("TDP43 timepoints with >= ", min_per_group, " per sex: ", paste(tp_sufficient_tdp, collapse = ", "))
message("WT timepoints with >= ", min_per_group, " per sex: ", paste(tp_sufficient_wt, collapse = ", "))

# ============================================================================
# 1a. BRAY-CURTIS PCoA PLOTS PER TIMEPOINT - Genotype (TDP43 vs WT)
# ============================================================================

## --- PCoA on all mice, faceted by timepoint ---
caz_mat_all <- log10(df_tot[, gene_families] + 1)
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
  mat <- log10(dfsel[, gene_families] + 1)
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
# 1b. BRAY-CURTIS PCoA PLOTS PER TIMEPOINT - Sex differences
# ============================================================================

## --- Bray-Curtis in TDP43 mice per timepoint ---
dfsel_tdp <- df_tot |> filter(Genotype == "TDP43" & Age_ints %in% tp_sufficient_tdp)
caz_mat_tdp <- dfsel_tdp[, gene_families]
caz_mat_tdp <- log10(caz_mat_tdp + 1)
rownames(caz_mat_tdp) <- dfsel_tdp$sampleID

bray_tdp <- vegan::vegdist(caz_mat_tdp, method = "bray")
pcoord_tdp <- ape::pcoa(bray_tdp, correction = "cailliez")
expl_variance_tdp <- pcoord_tdp$values$Rel_corr_eig * 100

dbray_tdp <- as.data.frame(pcoord_tdp$vectors[, c("Axis.1", "Axis.2")])
dbray_tdp$sampleID <- rownames(dbray_tdp)
dbray_tdp <- left_join(dbray_tdp, dfsel_tdp |> select(sampleID, Sex, Age_weeks, Age_ints), by = "sampleID")

# PERMANOVA per timepoint
bray_per_tp_tdp <- function(tp) {
  set.seed(14)
  dfsel <- dfsel_tdp |> filter(Age_ints == tp)
  mat <- log10(dfsel[, gene_families] + 1)
  rownames(mat) <- dfsel$sampleID
  bray <- vegan::vegdist(mat, method = "bray")
  res <- adonis2(bray ~ Sex, data = dfsel)
  return(res$`Pr(>F)`[1])
}

res_tdp <- data.frame(
  pvalue = sapply(tp_sufficient_tdp, bray_per_tp_tdp),
  Age_ints = tp_sufficient_tdp
) |> left_join(df_tot |> distinct(Age_ints, Age_weeks), by = "Age_ints")
res_tdp$Age_weeks <- factor(res_tdp$Age_weeks, levels = sort(unique(res_tdp$Age_weeks)))

dbray_tdp$Age_weeks <- factor(dbray_tdp$Age_weeks, levels = levels(res_tdp$Age_weeks))

(braycurt_tdp <- ggplot(dbray_tdp, aes(Axis.1, Axis.2)) +
    geom_point(aes(color = Sex), size = 2, alpha = 0.7) +
    xlab(paste0("PCo1 (", round(expl_variance_tdp[1], 1), "%)")) +
    ylab(paste0("PCo2 (", round(expl_variance_tdp[2], 1), "%)")) +
    scale_color_manual(values = pal_nejm()(2)) +
    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
    theme_Publication() +
    labs(color = "", fill = "", title = "CAZyme PCoA Bray-Curtis: TDP43") +
    stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
                 alpha = 0.1, linewidth = 1.0) +
    theme(legend.position = "top") +
    geom_text(data = res_tdp, aes(x = Inf, y = Inf, label = paste0("p = ", round(pvalue, 3))),
              hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
    facet_wrap(~ Age_weeks))
ggsave(file.path(res_dir, "PCoA_BrayCurtis_TDP43_sex.pdf"), braycurt_tdp, width = 14, height = 8)

## --- Bray-Curtis in WT mice per timepoint ---
dfsel_wt <- df_tot |> filter(Genotype == "WT" & Age_ints %in% tp_sufficient_wt)
caz_mat_wt <- dfsel_wt[, gene_families]
caz_mat_wt <- log10(caz_mat_wt + 1)
rownames(caz_mat_wt) <- dfsel_wt$sampleID

bray_wt <- vegan::vegdist(caz_mat_wt, method = "bray")
pcoord_wt <- ape::pcoa(bray_wt, correction = "cailliez")
expl_variance_wt <- pcoord_wt$values$Rel_corr_eig * 100

dbray_wt <- as.data.frame(pcoord_wt$vectors[, c("Axis.1", "Axis.2")])
dbray_wt$sampleID <- rownames(dbray_wt)
dbray_wt <- left_join(dbray_wt, dfsel_wt |> select(sampleID, Sex, Age_weeks, Age_ints), by = "sampleID")

bray_per_tp_wt <- function(tp) {
  set.seed(14)
  dfsel <- dfsel_wt |> filter(Age_ints == tp)
  mat <- log10(dfsel[, gene_families] + 1)
  rownames(mat) <- dfsel$sampleID
  bray <- vegan::vegdist(mat, method = "bray")
  res <- adonis2(bray ~ Sex, data = dfsel)
  return(res$`Pr(>F)`[1])
}

res_wt <- data.frame(
  pvalue = sapply(tp_sufficient_wt, bray_per_tp_wt),
  Age_ints = tp_sufficient_wt
) |> left_join(df_tot |> distinct(Age_ints, Age_weeks), by = "Age_ints")
res_wt$Age_weeks <- factor(res_wt$Age_weeks, levels = sort(unique(res_wt$Age_weeks)))

dbray_wt$Age_weeks <- factor(dbray_wt$Age_weeks, levels = levels(res_wt$Age_weeks))

(braycurt_wt <- ggplot(dbray_wt, aes(Axis.1, Axis.2)) +
    geom_point(aes(color = Sex), size = 2, alpha = 0.7) +
    xlab(paste0("PCo1 (", round(expl_variance_wt[1], 1), "%)")) +
    ylab(paste0("PCo2 (", round(expl_variance_wt[2], 1), "%)")) +
    scale_color_manual(values = pal_nejm()(2)) +
    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
    theme_Publication() +
    labs(color = "", fill = "", title = "CAZyme PCoA Bray-Curtis: WT") +
    stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
                 alpha = 0.1, linewidth = 1.0) +
    theme(legend.position = "top") +
    geom_text(data = res_wt, aes(x = Inf, y = Inf, label = paste0("p = ", round(pvalue, 3))),
              hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
    facet_wrap(~ Age_weeks))
ggsave(file.path(res_dir, "PCoA_BrayCurtis_WT_sex.pdf"), braycurt_wt, width = 14, height = 8)

# ============================================================================
# 2. GH78 AND GH106 BOXPLOTS AND LINE PLOTS OVER TIME
# ============================================================================

families_of_interest <- c("GH78", "GH106")
families_present <- families_of_interest[families_of_interest %in% names(df_tot)]
families_missing <- setdiff(families_of_interest, families_present)
if (length(families_missing) > 0) message("CAZy families not found: ", paste(families_missing, collapse = ", "))

## --- Boxplots per timepoint and genotype, sex-stratified ---
for (i in seq_along(families_present)) {
  gf <- families_present[i]
  df_tot$cazy_val <- log10(df_tot[[gf]] + 1)

  # Faceted plot: Age_weeks (columns) x Genotype (rows)
  pl <- ggplot(df_tot, aes(x = Sex, y = cazy_val)) +
    geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
    labs(title = gf, y = "log10(CPM + 1)", x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 2.5) +
    facet_grid(Genotype ~ Age_weeks) +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  ggsave(file.path(res_dir, paste0(gf, "_boxplots_per_timepoint.pdf")), pl,
         width = 2.5 * length(timepoints), height = 8)
}

## --- Line plots over time for GH78 and GH106 (like 2_7 TCAM lineplots) ---
line_plist <- list()
for (i in seq_along(families_present)) {
  gf <- families_present[i]
  df_tot$cazy_val <- log10(df_tot[[gf]] + 1)

  # Genotype line plot
  means_geno <- df_tot |>
    group_by(Genotype, Age_ints) |>
    summarise(mean = mean(cazy_val, na.rm = TRUE),
              se = sd(cazy_val, na.rm = TRUE) / sqrt(n()), .groups = "drop")

  pl_geno <- ggplot(means_geno, aes(x = Age_ints, y = mean, group = Genotype, color = Genotype)) +
    geom_line() +
    geom_point(size = 1) +
    geom_jitter(data = df_tot, aes(x = Age_ints, y = cazy_val, color = Genotype),
                width = 0.2, alpha = 0.5, inherit.aes = FALSE) +
    scale_color_manual(values = pal_nejm()(8)[c(3, 6)]) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3) +
    labs(title = paste0(gf, ": TDP43 vs WT"), x = "Age (weeks)", y = "log10(CPM + 1)") +
    scale_x_continuous(breaks = timepoints) +
    theme_Publication()

  line_plist[[length(line_plist) + 1]] <- pl_geno

  # Sex line plot within TDP43
  df_tdp <- df_tot |> filter(Genotype == "TDP43")
  means_sex <- df_tdp |>
    group_by(Sex, Age_ints) |>
    summarise(mean = mean(cazy_val, na.rm = TRUE),
              se = sd(cazy_val, na.rm = TRUE) / sqrt(n()), .groups = "drop")

  pl_sex <- ggplot(means_sex, aes(x = Age_ints, y = mean, group = Sex, color = Sex)) +
    geom_line() +
    geom_point(size = 1) +
    geom_jitter(data = df_tdp, aes(x = Age_ints, y = cazy_val, color = Sex),
                width = 0.2, alpha = 0.5, inherit.aes = FALSE) +
    scale_color_nejm() +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3) +
    labs(title = paste0(gf, ": sex diff in TDP43"), x = "Age (weeks)", y = "log10(CPM + 1)") +
    scale_x_continuous(breaks = timepoints) +
    theme_Publication()

  line_plist[[length(line_plist) + 1]] <- pl_sex
}

(line_plots <- ggarrange(plotlist = line_plist, ncol = 2, nrow = 2,
                          labels = LETTERS[1:length(line_plist)], common.legend = FALSE))
ggsave(file.path(res_dir, "GH78_GH106_lineplots.pdf"), line_plots, width = 12, height = 10)

## --- Interaction model: GH78 ~ Genotype * Sex ---
df_tot$GH78_log <- log10(df_tot$GH78 + 1)
model_gh78 <- lm(GH78_log ~ Genotype * Sex, data = df_tot)
summary(model_gh78)
confint(model_gh78)

## --- Interaction model: GH106 ~ Genotype * Sex ---
df_tot$GH106_log <- log10(df_tot$GH106 + 1)
model_gh106 <- lm(GH106_log ~ Genotype * Sex, data = df_tot)
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
    df_tp$mb <- log10(df_tp[[gf]] + 1)
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
    df_tp$mb <- log10(df_tp[[gf]] + 1)
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
      df_tp$cazy_val <- log10(df_tp[[gf]] + 1)

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
        labs(title = gf, caption = caption_text, y = "log10(CPM + 1)", x = "") +
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
      df_tp$cazy_val <- log10(df_tp[[gf]] + 1)

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
        labs(title = gf, caption = caption_text, y = "log10(CPM + 1)", x = "") +
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
      ggsave(file.path(res_dir, paste0("boxplots_sig_sexdiff_TDP43_", tp_int, "wk", suffix, ".pdf")),
             bp_combined, width = 4 * n_cols_bp, height = 4 * n_rows_bp)
    }
  }
}
