# Cayman analyses human cohort 2

# Library
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rio)
library(rstatix)
library(ggrepel)

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

stat <- rio::import("data/human_cohort2/cayman_sample_statistics.tsv") |> 
    mutate(Group = case_when(str_detect(sample, "AP") ~ "ALS", .default = "Control"))
cayman <- rio::import("data/human_cohort2/cayman_families_cpm_table.tsv")
head(cayman)[1:5,1:5]
rownames(cayman) <- cayman$family
cayman$family <- NULL
cayman <- as.data.frame(t(as.matrix(cayman)))
names(cayman)
dim(cayman)
head(cayman)[1:5,1:5]
cayman$sampleID <- rownames(cayman)
fam <- ncol(cayman)
clinical <- readRDS("data/human_cohort2/metadata.RDS")

# 1. Histogram of pct_cazy_reads
(gg_pct_cazy_hist <- gghistogram(stat, x = "pct_cazy_reads", bins = 10, fill = "#69b3a2", color = "black") +
  theme_Publication() +
  facet_wrap(~Group) +
  labs(title = "% CAZy Reads per Sample", x = "% CAZy Reads", y = "Count"))
#ggsave(file.path(qc_dir, "pct_cazy_reads_histogram.pdf"), gg_pct_cazy_hist, width = 6, height = 4)

# 2. Boxplot of pct_cazy_reads by timepoint (if timepoint column exists)
(gg_pct_cazy_box <- ggboxplot(stat, x = "Group", y = "pct_cazy_reads", fill = "Group", width = 0.4) +
    theme_Publication() +
    stat_compare_means() +
    scale_fill_jco() +
    labs(title = "% CAZy Reads by Group", x = "Group", y = "% CAZy Reads"))
#ggsave(file.path(qc_dir, "pct_cazy_reads_by_timepoint_boxplot.pdf"), gg_pct_cazy_box, width = 4, height = 6)

# 3. Scatterplot of pct_cazy_reads vs total_reads with correlation
  cor_pctcazy_total <- cor.test(stat$total_reads, stat$pct_cazy_reads, method = "spearman")
  subtitle_pctcazy_total <- paste0("Spearman's rho = ", round(cor_pctcazy_total$estimate, 3),
                                   ", p = ", formatC(cor_pctcazy_total$p.value, format = "e", digits = 2))
  (gg_pct_cazy_scatter <- ggplot(stat, aes(x = total_reads, y = pct_cazy_reads)) +
    geom_point(color = "#0072B2", alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    theme_Publication() +
    labs(title = "% CAZy Reads vs Total Reads", x = "Total Reads", y = "% CAZy Reads", subtitle = subtitle_pctcazy_total))
  #ggsave(file.path(qc_dir, "pct_cazy_reads_vs_total_reads.pdf"), gg_pct_cazy_scatter, width = 6, height = 4)

# 4. Histogram of richness
(gg_richness_hist <- gghistogram(stat, x = "richness", bins = 30, fill = "#E69F00", color = "black") +
  theme_Publication() +
  facet_wrap(~Group) +
  labs(title = "Sample Richness", x = "Richness", y = "Count"))
#ggsave(file.path(qc_dir, "richness_histogram.pdf"), gg_richness_hist, width = 6, height = 4)

# 5. Boxplot of richness by timepoint
(gg_richness_box <- ggboxplot(stat, x = "Group", y = "richness", fill = "Group", width = 0.4) +
  theme_Publication() +
  stat_compare_means() +
  scale_fill_jco() +
  labs(title = "Richness by Group", x = "Group", y = "Richness"))
# ggsave(file.path(qc_dir, "richness_by_timepoint_boxplot.pdf"), gg_richness_box, width = 4, height = 6)

(gg_richness_box <- ggboxplot(stat, x = "Group", y = "complexity", fill = "Group", width = 0.4) +
  theme_Publication() +
  stat_compare_means() +
  scale_fill_jco() +
  labs(title = "Complexity by Group", x = "Group", y = "Complexity"))
# ggsave(file.path(qc_dir, "richness_by_timepoint_boxplot.pdf"), gg_richness_box, width = 4, height = 6)

# 6. Scatterplot of richness vs complexity with correlation
cor_richness_complexity <- cor.test(stat$richness, stat$complexity, method = "spearman")
subtitle_richness_complexity <- paste0("Spearman's rho = ", round(cor_richness_complexity$estimate, 3),
                                        ", p = ", formatC(cor_richness_complexity$p.value, format = "e", digits = 2))
(gg_richness_complexity <- ggplot(stat, aes(x = richness, y = complexity)) +
  geom_point(color = "#D55E00", alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  theme_Publication() +
  labs(title = "Richness vs Complexity", x = "Richness", y = "Complexity", subtitle = subtitle_richness_complexity))
# ggsave(file.path(qc_dir, "richness_vs_complexity.pdf"), gg_richness_complexity, width = 6, height = 4)

# 8. Scatterplot of pct_cazy_reads vs richness with correlation
cor_pctcazy_richness <- cor.test(stat$richness, stat$pct_cazy_reads, method = "spearman")
subtitle_pctcazy_richness <- paste0("Spearman's rho = ", round(cor_pctcazy_richness$estimate, 3),
                                    ", p = ", formatC(cor_pctcazy_richness$p.value, format = "e", digits = 2))
(gg_pctcazy_richness <- ggplot(stat, aes(x = richness, y = pct_cazy_reads)) +
  geom_point(color = "#009E73", alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  theme_Publication() +
  labs(title = "% CAZy Reads vs Richness", x = "Richness", y = "% CAZy Reads", subtitle = subtitle_pctcazy_richness))
# ggsave(file.path(qc_dir, "pct_cazy_reads_vs_richness.pdf"), gg_pctcazy_richness, width = 6, height = 4)

# --- CAZy FAMILY ABUNDANCE DIFFERENCES: ALS vs Control ---
# Prepare: join cayman abundance with clinical metadata
gene_families <- setdiff(names(cayman), c("sampleID"))

# Filter CAZyme families: present in ≥10% of samples AND mean CPM > 1
prevalence <- colMeans(cayman[, gene_families] > 0)
mean_abundance <- colMeans(cayman[, gene_families])
gene_families <- gene_families[prevalence >= 0.30 & mean_abundance > 20]
message(length(gene_families), " CAZyme families retained after filtering (prevalence ≥10%, mean CPM >1)")

df_tot <- cayman |>
  left_join(clinical, by = c("sampleID" = "ID")) |>
  filter(!is.na(Group))

# Keep only families that survived the join as columns
gene_families <- gene_families[gene_families %in% names(df_tot)]

# Wilcoxon tests per CAZy family: ALS vs Control
wilcox_results <- purrr::map_dfr(gene_families, function(gf) {
  if (length(unique(df_tot[[gf]])) > 1 && length(unique(df_tot$Group)) > 1) {
    test <- wilcox.test(df_tot[[gf]] ~ df_tot$Group)
    tibble(gene_family = gf, pval_diag = test$p.value)
  } else {
    tibble(gene_family = gf, pval_diag = NA_real_)
  }
})

# Sex-stratified Wilcoxon tests
for (s in c("Female", "Male")) {
  df_sex <- df_tot |> filter(Sex == s)
  col_name <- paste0("pval_", tolower(s))
  sex_res <- purrr::map_dfr(gene_families, function(gf) {
    if (length(unique(df_sex[[gf]])) > 1 && length(unique(df_sex$Group)) > 1) {
      test <- wilcox.test(df_sex[[gf]] ~ df_sex$Group)
      tibble(gene_family = gf, pval = test$p.value)
    } else {
      tibble(gene_family = gf, pval = NA_real_)
    }
  })
  names(sex_res)[2] <- col_name
  wilcox_results <- wilcox_results |> left_join(sex_res, by = "gene_family")
}

# Sex differences within ALS and within Controls
for (g in c("ALS", "Control")) {
  df_grp <- df_tot |> filter(Group == g)
  col_name <- paste0("pval_sexdiff_", tolower(g))
  grp_res <- purrr::map_dfr(gene_families, function(gf) {
    if (length(unique(df_grp[[gf]])) > 1 && length(unique(df_grp$Sex)) > 1) {
      test <- wilcox.test(df_grp[[gf]] ~ df_grp$Sex)
      tibble(gene_family = gf, pval = test$p.value)
    } else {
      tibble(gene_family = gf, pval = NA_real_)
    }
  })
  names(grp_res)[2] <- col_name
  wilcox_results <- wilcox_results |> left_join(grp_res, by = "gene_family")
}

wilcox_results <- wilcox_results |>
  mutate(
    padj_diag = p.adjust(pval_diag, method = "fdr"),
    padj_female = p.adjust(pval_female, method = "fdr"),
    padj_male = p.adjust(pval_male, method = "fdr"),
    padj_sexdiff_als = p.adjust(pval_sexdiff_als, method = "fdr"),
    padj_sexdiff_control = p.adjust(pval_sexdiff_control, method = "fdr")
  ) |>
  arrange(pval_diag)

print(wilcox_results |> filter(padj_diag < 0.05))
print(wilcox_results |> filter(padj_sexdiff_control < 0.05))
write.csv2(wilcox_results, file.path(res_dir, "wilcoxons_cayman.csv"), row.names = FALSE)

# Linear model per CAZy family: Sex differences within ALS (log10 transformed)
df_als <- df_tot |> filter(Group == "ALS")
statres <- data.frame()
for (gf in gene_families) {
  df_als$mb <- log10(df_als[[gf]] + 1)
  tryCatch({
    model <- lm(mb ~ Sex, data = df_als)
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
  mutate(padj = p.adjust(pval, method = "fdr"))
write.csv2(statres, file.path(res_dir, "lm_cayman_sexdiff_als.csv"), row.names = FALSE)

# Volcano plot: sex differences in ALS
statres_sig <- statres |>
  mutate(sig_level = case_when(
    padj < 0.05 & abs(estimate) > 0.5 ~ "Significant & Large Effect",
    padj < 0.05 ~ "Significant",
    TRUE ~ "Not Significant"))
(p_volcano <- ggplot(statres_sig, aes(x = estimate, y = -log10(pval), color = sig_level)) +
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
       title = "CAZy Families: Sex Differences in ALS", color = ""))
ggsave(file.path(res_dir, "cayman_sexdiff_als_volcano.pdf"), p_volcano, width = 7, height = 7)

# --- Boxplots for specific CAZy families of interest ---
families_of_interest <- c("GH78", "GH106", "PL4", "PL11", "PL26", "GH43")

# Check which families are present in the data
families_present <- families_of_interest[families_of_interest %in% names(df_tot)]
families_missing <- setdiff(families_of_interest, families_present)
if (length(families_missing) > 0) message("CAZy families not found in data: ", paste(families_missing, collapse = ", "))

plist <- list()
for (i in seq_along(families_present)) {
  gf <- families_present[i]
  df_tot$cazy_val <- log10(df_tot[[gf]] + 1)

  # Sex-stratified Wilcoxon tests
  df_fem <- df_tot |> filter(Sex == "Female")
  df_mal <- df_tot |> filter(Sex == "Male")
  diagdiff_f <- wilcox.test(cazy_val ~ Group, data = df_fem)
  diagdiff_m <- wilcox.test(cazy_val ~ Group, data = df_mal)
  p_female <- format.pval(diagdiff_f$p.value, digits = 2)
  p_male <- format.pval(diagdiff_m$p.value, digits = 2)
  diag_diff_text <- paste0("ALS-Control p = ", p_female, " (women); p = ", p_male, " (men)")

  pl <- ggplot(df_tot, aes(x = Sex, y = cazy_val)) +
    geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
    labs(title = gf,
         caption = diag_diff_text,
         y = "log10(cpm + 1)",
         x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
          plot.caption = element_text(size = 10),
          legend.position = "none") +
    facet_wrap(~ Group)

  plist[[i]] <- pl
}

n_plots <- length(plist)
n_cols <- min(3, n_plots)
n_rows <- ceiling(n_plots / n_cols)

(plots_cazy <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
                         labels = LETTERS[1:n_plots],
                         nrow = n_rows, ncol = n_cols))
ggsave(file.path(res_dir, "cayman_families_of_interest_boxplots.pdf"), plots_cazy,
       width = 4 * n_cols, height = 5 * n_rows)

(plots_cazy_small <- ggarrange(plotlist = list(plist[[1]], plist[[2]]), common.legend = TRUE, legend = "bottom",
                         labels = LETTERS[1:2],
                         nrow = 1, ncol = 2))
ggsave(file.path(res_dir, "cayman_families_of_interest_boxplots.pdf"), plots_cazy_small,
       width = 8, height = 5)

# --- Linear model: GH78 ~ ALS * Sex + Age ---
df_tot$GH78_log <- log10(df_tot$GH78 + 1)
model_gh78 <- lm(GH78_log ~ Group * Sex, data = df_tot)
summary(model_gh78)
confint(model_gh78)

# --- Linear model: GH43 ~ ALS * Sex + Age ---
df_tot$GH43_log <- log10(df_tot$GH43 + 1)
model_gh43 <- lm(GH43_log ~ Group * Sex, data = df_tot)
summary(model_gh43)
confint(model_gh43)

# --- Linear model: GH43 ~ ALS * Sex + Age ---
df_tot$GH106_log <- log10(df_tot$GH106 + 1)
model_gh106 <- lm(GH106_log ~ Group * Sex, data = df_tot)
summary(model_gh106)
confint(model_gh106)

# --- Correlation plot: GH106 vs ALSFRS (ALS only) ---
df_als$GH106_log <- log10(df_als$GH106 + 1)
df_als$ALSFRS_num <- as.numeric(df_als$ALSFRS)
df_corr <- df_als |> filter(!is.na(ALSFRS_num))

cor_gh106 <- cor.test(df_corr$GH106_log, df_corr$ALSFRS_num, method = "spearman")
subtitle_gh106 <- paste0("Spearman's rho = ", round(cor_gh106$estimate, 3),
                         ", p = ", formatC(cor_gh106$p.value, format = "e", digits = 2))

(p_gh106_alsfrs <- ggplot(df_corr, aes(x = GH106_log, y = ALSFRS_num)) +
  geom_point(aes(color = Sex), alpha = 0.7, size = 2) +
    facet_wrap(~Sex) +
    stat_cor(method = "spearman") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(values = pal_nejm()(2)) +
  theme_Publication() +
  labs(x = "GH106 (log10 CPM + 1)", y = "ALSFRS",
       title = "GH106 vs ALSFRS in ALS", subtitle = subtitle_gh106))
ggsave(file.path(res_dir, "cayman_GH106_vs_ALSFRS.pdf"), p_gh106_alsfrs, width = 6, height = 5)

# --- Beta diversity: CAZyme families ---
# Prepare numeric matrix (all families, unfiltered)
caz_mat <- cayman[, setdiff(names(cayman), "sampleID")]
caz_mat <- caz_mat[rownames(caz_mat) %in% clinical$ID, ]
# caz_mat <- log10(caz_mat + 1)
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
        labs(color = "", fill = "", title = "Patients with ALS") +
        stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
            alpha = 0.1, linewidth = 1.0) +
        theme(legend.position = "top"))

(pl_caz1 <- pl_caz1 + annotate("text", x = Inf, y = Inf, label = paste0("Sex: p = ", round(res_als$`Pr(>F)`[1], 3)),
                    hjust = 1.1, vjust = 1.1, size = 3))

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
        labs(color = "", fill = "", title = "Controls") +
        stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
            alpha = 0.1, linewidth = 1.0) +
        theme(legend.position = "top"))

pl_caz2 <- pl_caz2 + annotate("text", x = Inf, y = Inf, label = paste0("Sex: p = ", round(res_ctrl$`Pr(>F)`[1], 3)),
                    hjust = 1.1, vjust = 1.1, size = 3)
pl_caz2

ggarrange(pl_caz1, pl_caz2, nrow = 1, labels = LETTERS[1:2], common.legend = TRUE, legend = "bottom")
ggsave(file.path(res_dir, "cayman_PCoA_BrayCurtis.pdf"), width = 12, height = 6)

# --- PCA on log10-transformed CAZyme family abundances ---
caz_log <- log10(cayman[, setdiff(names(cayman), "sampleID")] + 1)
caz_log <- caz_log[rownames(caz_log) %in% clinical$ID, ]

# PCA in ALS
caz_log_als <- caz_log[which(rownames(caz_log) %in% dfals$ID), ]
caz_log_als <- caz_log_als[, apply(caz_log_als, 2, var) > 0]
pca_als <- prcomp(caz_log_als, center = TRUE, scale. = TRUE)
pca_var_als <- summary(pca_als)$importance[2, ] * 100
pca_df_als <- as.data.frame(pca_als$x[, 1:2])
pca_df_als$ID <- rownames(pca_df_als)
pca_df_als <- left_join(pca_df_als, clinical, by = "ID")

set.seed(14)
euc_als <- vegan::vegdist(caz_log_als, method = "euclidean")
res_pca_als <- adonis2(euc_als ~ Sex, data = pca_df_als)
res_pca_als <- as.data.frame(res_pca_als)

(pl_pca1 <- ggplot(pca_df_als, aes(PC1, PC2)) +
    geom_point(aes(color = Sex), size = 3, alpha = 0.7) +
    xlab(paste0("PC1 (", round(pca_var_als[1], 1), "%)")) +
    ylab(paste0("PC2 (", round(pca_var_als[2], 1), "%)")) +
    scale_color_manual(values = pal_nejm()(2)) +
    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
    theme_Publication() +
    labs(color = "", fill = "", title = "Patients with ALS") +
    stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
                 alpha = 0.1, linewidth = 1.0) +
    theme(legend.position = "top") +
    annotate("text", x = Inf, y = Inf, label = paste0("Sex: p = ", round(res_pca_als$`Pr(>F)`[1], 3)),
             hjust = 1.1, vjust = 1.1, size = 3))

# PCA in Controls
caz_log_ctrl <- caz_log[which(rownames(caz_log) %in% dfctrl$ID), ]
caz_log_ctrl <- caz_log_ctrl[, apply(caz_log_ctrl, 2, var) > 0]
pca_ctrl <- prcomp(caz_log_ctrl, center = TRUE, scale. = TRUE)
pca_var_ctrl <- summary(pca_ctrl)$importance[2, ] * 100
pca_df_ctrl <- as.data.frame(pca_ctrl$x[, 1:2])
pca_df_ctrl$ID <- rownames(pca_df_ctrl)
pca_df_ctrl <- left_join(pca_df_ctrl, clinical, by = "ID")

set.seed(14)
euc_ctrl <- vegan::vegdist(caz_log_ctrl, method = "euclidean")
res_pca_ctrl <- adonis2(euc_ctrl ~ Sex, data = pca_df_ctrl)
res_pca_ctrl <- as.data.frame(res_pca_ctrl)

(pl_pca2 <- ggplot(pca_df_ctrl, aes(PC1, PC2)) +
    geom_point(aes(color = Sex), size = 3, alpha = 0.7) +
    xlab(paste0("PC1 (", round(pca_var_ctrl[1], 1), "%)")) +
    ylab(paste0("PC2 (", round(pca_var_ctrl[2], 1), "%)")) +
    scale_color_manual(values = pal_nejm()(2)) +
    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
    theme_Publication() +
    labs(color = "", fill = "", title = "Controls") +
    stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
                 alpha = 0.1, linewidth = 1.0) +
    theme(legend.position = "top") +
    annotate("text", x = Inf, y = Inf, label = paste0("Sex: p = ", round(res_pca_ctrl$`Pr(>F)`[1], 3)),
             hjust = 1.1, vjust = 1.1, size = 3))

ggarrange(pl_pca1, pl_pca2, nrow = 1, labels = LETTERS[1:2], common.legend = TRUE, legend = "bottom")
ggsave(file.path(res_dir, "cayman_PCA_log10.pdf"), width = 12, height = 6)
