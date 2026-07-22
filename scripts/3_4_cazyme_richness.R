# CAZyme richness and Shannon diversity over time
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Load libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rio)
library(vegan)     # For Shannon diversity
library(lme4)
library(lmerTest)  # For the linear mixed model (lmer) with p-values via Satterthwaite

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

# 1. Load data
# Richness = number of CAZyme families with CPM > 1 per sample (Ducarmon et al. definition,
#            applied to CPM here since gene-length-corrected RPKM isn't available at the
#            family level).
# Shannon  = Shannon diversity of the family-level CPM-normalized abundances per sample
#            (how evenly community abundance is spread across families, complementing
#            richness which only counts presence above the threshold).
# The sample_statistics.tsv "richness"/"complexity" columns are a different, unnormalized
# alignment-stage statistic (unique reference genes hit, pre-RPKM/CPM) and are not used here.
fam_cpm <- rio::import("data/microbiome/families_cpm_table.tsv")
rownames(fam_cpm) <- fam_cpm$family
fam_cpm$family <- NULL

rich    <- tibble(sample = colnames(fam_cpm), richness = colSums(fam_cpm > 1))
shannon_vec <- vegan::diversity(t(fam_cpm), index = "shannon")
shannon <- tibble(sample = names(shannon_vec), shannon = as.numeric(shannon_vec))

stat <- rio::import("data/microbiome/sample_statistics.tsv") |>
  select(sample, passed_reads)
meta <- readRDS("data/microbiome/meta_microbiome.RDS")

# 2. Join a metric with metadata
join_metadata <- function(metric_df) {
  metric_df |>
    left_join(stat, by = "sample") |>
    rename(sample_id = sample) |>
    mutate(sample_id = str_trim(as.character(sample_id))) |>
    left_join(
      meta |>
        mutate(ID = str_trim(as.character(ID))) |>
        select(ID, Sex, Genotype, MouseID, Age_weeks, Age_ints, GenotypePerSex, LowCount),
      by = c("sample_id" = "ID")
    ) |>
    filter(Age_ints %in% c(6, 8, 10, 12, 14, 16)) |>
    filter(is.na(LowCount) | LowCount == FALSE) |>
    mutate(Age_fac = factor(Age_ints))
}

rich_plot    <- join_metadata(rich)
shannon_plot <- join_metadata(shannon)

# 3. Summary statistics for plotting
summarise_metric <- function(data_plot, metric) {
  data_plot |>
    group_by(Genotype, Sex, Age_ints) |>
    summarise(
      mean = mean(.data[[metric]], na.rm = TRUE),
      se   = sd(.data[[metric]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

rich_means_ses    <- summarise_metric(rich_plot, "richness")
shannon_means_ses <- summarise_metric(shannon_plot, "shannon")

# 3b. Shared y-axis range across genotypes, so the two panels are directly comparable.
# Padding is a fraction of the observed range (not of the raw value), so it scales
# sensibly regardless of the metric's absolute magnitude.
# zero_anchor = TRUE starts the axis at 0 (appropriate for a count like richness);
# zero_anchor = FALSE zooms into the actual data range (appropriate for Shannon
# diversity, which only varies by a small amount around ~4.3).
shared_y_limits <- function(data_plot, means_ses, metric, zero_anchor = TRUE, pad_frac = 0.15) {
  y_ceiling <- max(means_ses$mean + means_ses$se, na.rm = TRUE)
  y_min_raw <- min(data_plot[[metric]], na.rm = TRUE)
  y_max_raw <- max(data_plot[[metric]], na.rm = TRUE)
  y_pad     <- (y_max_raw - y_min_raw) * pad_frac
  y_lower   <- if (zero_anchor) 0 else y_min_raw - y_pad
  list(
    y_ceiling = y_ceiling,
    y_pad     = y_pad,
    y_limits  = c(y_lower, max(y_max_raw, y_ceiling + y_pad) + y_pad)
  )
}

rich_range    <- shared_y_limits(rich_plot,    rich_means_ses,    "richness", zero_anchor = TRUE)
shannon_range <- shared_y_limits(shannon_plot, shannon_means_ses, "shannon",  zero_anchor = FALSE)

# 4. Plot per genotype with LMM baseline-referenced sex-difference contrasts
plot_metric_for_genotype <- function(genotype_name, data_plot, means_ses, metric, y_label,
                                      y_ceiling, y_pad, y_limits, title_prefix = "") {

  means_df <- means_ses |> filter(Genotype == genotype_name)
  raw_df   <- data_plot |> filter(Genotype == genotype_name)

  # Difference-in-differences: is the Female-Male difference at age X different from the
  # Female-Male difference at baseline (week 6)? Tested by refitting metric ~ Sex * Age_fac
  # on just that pair of timepoints and reading off the Sex:Age_fac interaction coefficient
  # directly - that coefficient *is* the DiD estimate. BH-adjusted across the non-baseline ages.
  # scale(passed_reads) adjusts for residual sequencing-depth differences between samples.
  baseline_age <- min(raw_df$Age_ints)
  other_ages   <- setdiff(sort(unique(raw_df$Age_ints)), baseline_age)

  contr_rows <- lapply(other_ages, function(age) {
    d2 <- raw_df |> filter(Age_ints %in% c(baseline_age, age))
    d2$Age_fac <- relevel(droplevels(factor(d2$Age_ints)), ref = as.character(baseline_age))
    f <- as.formula(paste(metric, "~ Sex * Age_fac + scale(passed_reads) + (1 | MouseID)"))
    m <- tryCatch(lmerTest::lmer(f, data = d2), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    coefs <- summary(m)$coefficients
    int_row <- grep("^Sex.*:Age_fac", rownames(coefs), value = TRUE)
    if (length(int_row) == 0) return(NULL)
    data.frame(Age_ints = age,
               estimate = coefs[int_row[1], "Estimate"],
               p.value  = coefs[int_row[1], "Pr(>|t|)"])
  })

  contr_df <- bind_rows(contr_rows)
  if (nrow(contr_df) > 0) {
    contr_df$q.value <- p.adjust(contr_df$p.value, method = "fdr")
    contr_df$sig <- case_when(
      contr_df$q.value < 0.001 ~ "***",
      contr_df$q.value < 0.01  ~ "**",
      contr_df$q.value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
    contr_df$y_pos <- y_ceiling + y_pad
  } else {
    contr_df <- NULL
  }

  p <- ggplot(means_df, aes(x = Age_ints, y = mean, group = Sex, color = Sex)) +
    geom_line() +
    geom_point(size = 0.5) +
    geom_jitter(data = raw_df, aes(x = Age_ints, y = .data[[metric]], color = Sex),
                width = 0.1, alpha = 0.4) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_color_nejm() +
    labs(title = paste0(title_prefix, genotype_name), x = "Age (weeks)", y = y_label) +
    scale_x_continuous(breaks = c(6, 8, 10, 12, 14, 16)) +
    coord_cartesian(ylim = y_limits) +
    theme_Publication()

  if (!is.null(contr_df)) {
    p <- p + geom_text(
      data = contr_df |> filter(sig != "ns"),
      aes(x = Age_ints, y = y_pos, label = sig),
      inherit.aes = FALSE, size = 5, vjust = 0, fontface = "bold", color = "black"
    )
  }

  list(plot = p, contrasts = contr_df)
}

# --- Richness ---
res_wt    <- plot_metric_for_genotype("WT", rich_plot, rich_means_ses, "richness",
                                       "Number of CAZyme families",
                                       rich_range$y_ceiling, rich_range$y_pad, rich_range$y_limits)
res_tdp43 <- plot_metric_for_genotype("TDP43", rich_plot, rich_means_ses, "richness",
                                       "Number of CAZyme families",
                                       rich_range$y_ceiling, rich_range$y_pad, rich_range$y_limits)

richness_combined <- ggarrange(
  res_tdp43$plot, res_wt$plot,
  ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "bottom"
)
richness_combined

ggsave(file.path(res_dir, "cazyme_richness.pdf"), richness_combined, width = 8, height = 5, units = "in")

# Save LMM contrasts: change in the sex difference vs baseline (week 6), per genotype
richness_contrasts_all <- bind_rows(
  if (!is.null(res_tdp43$contrasts)) data.frame(res_tdp43$contrasts, Genotype = "TDP43"),
  if (!is.null(res_wt$contrasts))    data.frame(res_wt$contrasts, Genotype = "WT")
)
write.csv2(richness_contrasts_all, file.path(res_dir, "cazyme_richness_sex_contrasts.csv"), row.names = FALSE)

# --- Shannon diversity ---
res_wt_shannon    <- plot_metric_for_genotype("WT", shannon_plot, shannon_means_ses, "shannon",
                                               "Shannon diversity",
                                               shannon_range$y_ceiling, shannon_range$y_pad, shannon_range$y_limits,
                                               title_prefix = "CAZyme diversity: ")
res_tdp43_shannon <- plot_metric_for_genotype("TDP43", shannon_plot, shannon_means_ses, "shannon",
                                               "Shannon diversity",
                                               shannon_range$y_ceiling, shannon_range$y_pad, shannon_range$y_limits,
                                               title_prefix = "CAZyme diversity: ")

shannon_combined <- ggarrange(
  res_tdp43_shannon$plot, res_wt_shannon$plot,
  ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "bottom"
)
shannon_combined

ggsave(file.path(res_dir, "cazyme_shannon.pdf"), shannon_combined, width = 8, height = 5, units = "in")

# Save LMM contrasts: change in the sex difference vs baseline (week 6), per genotype
shannon_contrasts_all <- bind_rows(
  if (!is.null(res_tdp43_shannon$contrasts)) data.frame(res_tdp43_shannon$contrasts, Genotype = "TDP43"),
  if (!is.null(res_wt_shannon$contrasts))    data.frame(res_wt_shannon$contrasts, Genotype = "WT")
)
write.csv2(shannon_contrasts_all, file.path(res_dir, "cazyme_shannon_sex_contrasts.csv"), row.names = FALSE)
