# CAZyme richness over time
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Load libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rio)
library(lme4)     # For Poisson GLMM (glmer)
library(emmeans)  # For FDR-adjusted longitudinal contrasts

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
rich <- rio::import("data/microbiome/sample_statistics.tsv")
meta <- readRDS("data/microbiome/meta_microbiome.RDS")

# 2. Join richness with metadata (raw unlogged richness)
rich_plot <- rich |>
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

# 3. Summary statistics for plotting
rich_means_ses <- rich_plot |>
  group_by(Genotype, Sex, Age_ints) |>
  summarise(
    mean = mean(richness, na.rm = TRUE),
    se   = sd(richness, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 3b. Shared y-axis range across genotypes, so the two panels are directly comparable
# Upper bound covers the raw (jittered) points, including outliers, not just mean +/- SE
y_ceiling_shared <- max(rich_means_ses$mean + rich_means_ses$se, na.rm = TRUE)
y_max_raw        <- max(rich_plot$richness, na.rm = TRUE)
y_pad            <- y_max_raw * 0.05
y_limits_shared  <- c(0, max(y_max_raw, y_ceiling_shared + y_pad) + y_pad)

# 4. Plot per genotype with Poisson GLMM sex contrasts
plot_richness_for_genotype <- function(genotype_name, y_ceiling, y_limits) {

  means_df <- rich_means_ses |> filter(Genotype == genotype_name)
  raw_df <- rich_plot |> filter(Genotype == genotype_name)

  # Richness ~ Sex * Age_fac with a random intercept for individual mice
  glmm_model <- tryCatch(
    lme4::glmer(richness ~ Sex * Age_fac + (1 | MouseID), data = raw_df, family = poisson(link = "log")),
    error = function(e) NULL
  )

  contr_df <- NULL
  if (!is.null(glmm_model)) {
    emm      <- emmeans::emmeans(glmm_model, ~ Sex | Age_fac, type = "response")
    contr_df <- as.data.frame(contrast(emm, method = "pairwise", adjust = "fdr"))

    contr_df$Age_ints <- as.integer(as.character(contr_df$Age_fac))
    contr_df$sig <- case_when(
      contr_df$p.value < 0.001 ~ "***",
      contr_df$p.value < 0.01  ~ "**",
      contr_df$p.value < 0.05  ~ "*",
      TRUE ~ "ns"
    )

    contr_df$y_pos <- y_ceiling + y_pad
  }

  p <- ggplot(means_df, aes(x = Age_ints, y = mean, group = Sex, color = Sex)) +
    geom_line() +
    geom_point(size = 0.5) +
    geom_jitter(data = raw_df, aes(x = Age_ints, y = richness, color = Sex),
                width = 0.1, alpha = 0.4) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_color_nejm() +
    labs(title = genotype_name, x = "Age (weeks)", y = "CAZyme richness") +
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

res_wt    <- plot_richness_for_genotype("WT", y_ceiling_shared, y_limits_shared)
res_tdp43 <- plot_richness_for_genotype("TDP43", y_ceiling_shared, y_limits_shared)

richness_combined <- ggarrange(
  res_tdp43$plot, res_wt$plot,
  ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "bottom"
)
richness_combined

ggsave(file.path(res_dir, "cazyme_richness.pdf"), richness_combined, width = 8, height = 5, units = "in")

# Save GLMM contrasts (sex differences per timepoint, per genotype)
contrasts_all <- bind_rows(
  if (!is.null(res_tdp43$contrasts)) data.frame(res_tdp43$contrasts, Genotype = "TDP43"),
  if (!is.null(res_wt$contrasts))    data.frame(res_wt$contrasts, Genotype = "WT")
)
write.csv2(contrasts_all, file.path(res_dir, "cazyme_richness_sex_contrasts.csv"), row.names = FALSE)
