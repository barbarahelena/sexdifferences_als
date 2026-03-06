# Alpha diversity plots ALS project - Line plots with LMM statistics
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ggsci)
library(lme4)
library(lmerTest)
library(emmeans)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(family = "Helvetica"),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.9)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.9)),
                axis.text.x = element_text(angle = 0),
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

res_dir <- "results/microbiome/alphadiversity"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

## Load data
mb <- readRDS("data/microbiome/microbiome.RDS")
df <- readRDS("data/microbiome/meta_microbiome.RDS")

## Shannon diversity
shannon <- vegan::diversity(mb, index = "shannon")
df_shan <- data.frame(ID = names(shannon), shannon = shannon)
df_shan <- left_join(df_shan, df, by = "ID")

# Filter: remove low-count samples, restrict to weeks 6-18
lowcountids <- df$ID[df$LowCount == TRUE]
message("Removing ", length(lowcountids), " low-count samples")

df_shan <- df_shan |>
  filter(!is.na(Genotype)) |>
  filter(!ID %in% lowcountids) |>
  filter(Age_ints >= 6 & Age_ints <= 18)
message(nrow(df_shan), " samples after filtering")

# Set reference levels and age factor
df_shan$Genotype <- relevel(factor(df_shan$Genotype), ref = "WT")
df_shan$Sex      <- relevel(factor(df_shan$Sex),      ref = "Female")
df_shan$Age_fac  <- factor(df_shan$Age_ints)

# Define timepoints
timepoints <- sort(unique(df_shan$Age_ints))

# Timepoints with sufficient TDP43 samples (>=3 per sex)
tp_sufficient_tdp <- timepoints[sapply(timepoints, function(tp) {
  d <- df_shan |> filter(Genotype == "TDP43" & Age_ints == tp)
  nrow(d) > 0 && all(table(d$Sex) >= 3)
})]

print(table(interaction(df_shan$Sex, df_shan$Genotype), df_shan$Age_ints))

# ============================================================================
# 1. LINE PLOT: Genotype (TDP43 vs WT) over time
# ============================================================================

means_geno <- df_shan |>
  group_by(Genotype, Age_ints) |>
  summarise(mean = mean(shannon, na.rm = TRUE),
            se = sd(shannon, na.rm = TRUE) / sqrt(n()), .groups = "drop")

lmm_geno <- tryCatch(
  lmerTest::lmer(shannon ~ Genotype * Age_fac + (1 | MouseID), data = df_shan, REML = FALSE),
  error = function(e) NULL
)

if (!is.null(lmm_geno)) {
  anova_geno <- anova(lmm_geno)
  emm_geno   <- emmeans(lmm_geno, ~ Genotype | Age_fac)
  contr_geno <- as.data.frame(contrast(emm_geno, method = "pairwise", adjust = "fdr"))
  contr_geno$Age_ints <- as.integer(as.character(contr_geno$Age_fac))
  contr_geno$sig <- case_when(
    contr_geno$p.value < 0.001 ~ "***",
    contr_geno$p.value < 0.01  ~ "**",
    contr_geno$p.value < 0.05  ~ "*",
    TRUE ~ "ns"
  )
  y_max_geno <- max(means_geno$mean + means_geno$se, na.rm = TRUE) + 0.1
  contr_geno$y_pos <- y_max_geno
} else {
  contr_geno <- NULL
}

pl_shan_geno <- ggplot(means_geno, aes(x = Age_ints, y = mean, group = Genotype, color = Genotype)) +
  geom_line() +
  geom_point(size = 1) +
  geom_jitter(data = df_shan, aes(x = Age_ints, y = shannon, color = Genotype),
              width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  scale_color_manual(values = pal_nejm()(8)[c(6, 3)]) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3) +
  labs(title = "Shannon index: TDP43 vs WT", x = "Age (weeks)", y = "Shannon index") +
  scale_x_continuous(breaks = timepoints) +
  theme_Publication()
if (!is.null(contr_geno)) {
  pl_shan_geno <- pl_shan_geno +
    geom_text(data = contr_geno |> filter(sig != "ns"),
              aes(x = Age_ints, y = y_pos, label = sig),
              inherit.aes = FALSE, size = 6, vjust = 0)
}
ggsave(file.path(res_dir, "shannon_lineplot_genotype.pdf"), pl_shan_geno, width = 8, height = 5)

# ============================================================================
# 2. LINE PLOT: Sex differences within TDP43
# ============================================================================

df_tdp <- df_shan |> filter(Genotype == "TDP43" & Age_ints %in% tp_sufficient_tdp)

means_sex_tdp <- df_tdp |>
  group_by(Sex, Age_ints) |>
  summarise(mean = mean(shannon, na.rm = TRUE),
            se = sd(shannon, na.rm = TRUE) / sqrt(n()), .groups = "drop")

lmm_sex_tdp <- tryCatch(
  lmerTest::lmer(shannon ~ Sex * Age_fac + (1 | MouseID), data = df_tdp, REML = FALSE),
  error = function(e) NULL
)

if (!is.null(lmm_sex_tdp)) {
  anova_sex_tdp <- anova(lmm_sex_tdp)
  emm_sex_tdp   <- emmeans(lmm_sex_tdp, ~ Sex | Age_fac)
  contr_sex_tdp <- as.data.frame(contrast(emm_sex_tdp, method = "pairwise", adjust = "fdr"))
  contr_sex_tdp$Age_ints <- as.integer(as.character(contr_sex_tdp$Age_fac))
  contr_sex_tdp$sig <- case_when(
    contr_sex_tdp$p.value < 0.001 ~ "***",
    contr_sex_tdp$p.value < 0.01  ~ "**",
    contr_sex_tdp$p.value < 0.05  ~ "*",
    TRUE ~ "ns"
  )
  y_max_sex_tdp <- max(means_sex_tdp$mean + means_sex_tdp$se, na.rm = TRUE) + 0.1
  contr_sex_tdp$y_pos <- y_max_sex_tdp
} else {
  contr_sex_tdp <- NULL
}

pl_shan_sex_tdp <- ggplot(means_sex_tdp, aes(x = Age_ints, y = mean, group = Sex, color = Sex)) +
  geom_line() +
  geom_point(size = 1) +
  geom_jitter(data = df_tdp, aes(x = Age_ints, y = shannon, color = Sex),
              width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  scale_color_nejm() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3) +
  labs(title = "Shannon index: TDP43 males vs females", x = "Age (weeks)", y = "Shannon index") +
  scale_x_continuous(breaks = tp_sufficient_tdp) +
  theme_Publication()
if (!is.null(contr_sex_tdp)) {
  pl_shan_sex_tdp <- pl_shan_sex_tdp +
    geom_text(data = contr_sex_tdp |> filter(sig != "ns"),
              aes(x = Age_ints, y = y_pos, label = sig),
              inherit.aes = FALSE, size = 6, vjust = 0)
}
ggsave(file.path(res_dir, "shannon_lineplot_sex_tdp43.pdf"), pl_shan_sex_tdp, width = 8, height = 5)

# ============================================================================
# 3. LINE PLOT: Sex differences within WT
# ============================================================================

df_wt <- df_shan |> filter(Genotype == "WT")

means_sex_wt <- df_wt |>
  group_by(Sex, Age_ints) |>
  summarise(mean = mean(shannon, na.rm = TRUE),
            se = sd(shannon, na.rm = TRUE) / sqrt(n()), .groups = "drop")

lmm_sex_wt <- tryCatch(
  lmerTest::lmer(shannon ~ Sex * Age_fac + (1 | MouseID), data = df_wt, REML = FALSE),
  error = function(e) NULL
)

if (!is.null(lmm_sex_wt)) {
  anova_sex_wt <- anova(lmm_sex_wt)
  emm_sex_wt   <- emmeans(lmm_sex_wt, ~ Sex | Age_fac)
  contr_sex_wt <- as.data.frame(contrast(emm_sex_wt, method = "pairwise", adjust = "fdr"))
  contr_sex_wt$Age_ints <- as.integer(as.character(contr_sex_wt$Age_fac))
  contr_sex_wt$sig <- case_when(
    contr_sex_wt$p.value < 0.001 ~ "***",
    contr_sex_wt$p.value < 0.01  ~ "**",
    contr_sex_wt$p.value < 0.05  ~ "*",
    TRUE ~ "ns"
  )
  y_max_sex_wt <- max(means_sex_wt$mean + means_sex_wt$se, na.rm = TRUE) + 0.1
  contr_sex_wt$y_pos <- y_max_sex_wt
} else {
  contr_sex_wt <- NULL
}

pl_shan_sex_wt <- ggplot(means_sex_wt, aes(x = Age_ints, y = mean, group = Sex, color = Sex)) +
  geom_line() +
  geom_point(size = 1) +
  geom_jitter(data = df_wt, aes(x = Age_ints, y = shannon, color = Sex),
              width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  scale_color_nejm() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3) +
  labs(title = "Shannon index: WT males vs females", x = "Age (weeks)", y = "Shannon index") +
  scale_x_continuous(breaks = timepoints) +
  theme_Publication()
if (!is.null(contr_sex_wt)) {
  pl_shan_sex_wt <- pl_shan_sex_wt +
    geom_text(data = contr_sex_wt |> filter(sig != "ns"),
              aes(x = Age_ints, y = y_pos, label = sig),
              inherit.aes = FALSE, size = 6, vjust = 0)
}
ggsave(file.path(res_dir, "shannon_lineplot_sex_wt.pdf"), pl_shan_sex_wt, width = 8, height = 5)

# ============================================================================
# 4. LINE PLOT + LMM: Genotype * Sex interaction
# ============================================================================

# Create a 4-group factor for colour mapping
df_shan$GenoSex <- interaction(df_shan$Genotype, df_shan$Sex, sep = " ")

means_int <- df_shan |>
  group_by(Genotype, Sex, GenoSex, Age_ints) |>
  summarise(mean = mean(shannon, na.rm = TRUE),
            se = sd(shannon, na.rm = TRUE) / sqrt(n()), .groups = "drop")

# LMM: Genotype * Sex + Age_fac + (1 | MouseID)
# Tests whether the Genotype effect differs between sexes (interaction term)
lmm_int <- tryCatch(
  lmerTest::lmer(shannon ~ Genotype * Sex + Age_fac + (1 | MouseID), data = df_shan, REML = FALSE),
  error = function(e) NULL
)

if (!is.null(lmm_int)) {
  anova_int <- anova(lmm_int)
  print(anova_int)
  # Pairwise contrasts for all 4 groups, FDR-corrected
  emm_int   <- emmeans(lmm_int, ~ Genotype * Sex)
  contr_int <- as.data.frame(contrast(emm_int, method = "pairwise", adjust = "fdr"))
  print(contr_int)
}

# Colours: WT Female, WT Male, TDP43 Female, TDP43 Male
int_cols <- pal_nejm()(8)[c(6, 6, 3, 3)]
int_ltys  <- c("solid", "dashed", "solid", "dashed")   # solid = Female, dashed = Male

pl_shan_int <- ggplot(means_int, aes(x = Age_ints, y = mean, group = GenoSex,
                                      color = GenoSex, linetype = Sex)) +
  geom_line() +
  geom_point(size = 1) +
  geom_jitter(data = df_shan, aes(x = Age_ints, y = shannon,
                                   color = GenoSex, shape = Sex),
              width = 0.2, alpha = 0.35, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3) +
  scale_color_manual(values = pal_nejm()(8)[c(6, 6, 3, 3)],
                     labels = levels(factor(means_int$GenoSex))) +
  scale_linetype_manual(values = c(Female = "solid", Male = "dashed")) +
  labs(title = "Shannon index: Genotype × Sex", x = "Age (weeks)", y = "Shannon index",
       color = "", linetype = "Sex") +
  scale_x_continuous(breaks = timepoints) +
  theme_Publication()
ggsave(file.path(res_dir, "shannon_lineplot_genotype_sex_interaction.pdf"), pl_shan_int, width = 8, height = 5)

# ============================================================================
# Save LMM results
# ============================================================================

lmm_anova_all <- bind_rows(
  if (!is.null(lmm_geno))    data.frame(anova_geno,    term = rownames(anova_geno),    model = "shannon_genotype"),
  if (!is.null(lmm_sex_tdp)) data.frame(anova_sex_tdp, term = rownames(anova_sex_tdp), model = "shannon_sex_TDP43"),
  if (!is.null(lmm_sex_wt))  data.frame(anova_sex_wt,  term = rownames(anova_sex_wt),  model = "shannon_sex_WT"),
  if (!is.null(lmm_int))     data.frame(anova_int,     term = rownames(anova_int),     model = "shannon_genotype_sex_interaction")
)
lmm_contr_all <- bind_rows(
  if (!is.null(contr_geno))    data.frame(contr_geno,    model = "shannon_genotype"),
  if (!is.null(contr_sex_tdp)) data.frame(contr_sex_tdp, model = "shannon_sex_TDP43"),
  if (!is.null(contr_sex_wt))  data.frame(contr_sex_wt,  model = "shannon_sex_WT"),
  if (!is.null(lmm_int))       data.frame(contr_int,     model = "shannon_genotype_sex_interaction")
)

write.csv2(lmm_anova_all, file.path(res_dir, "lmm_shannon_anova.csv"),     row.names = FALSE)
write.csv2(lmm_contr_all, file.path(res_dir, "lmm_shannon_contrasts.csv"), row.names = FALSE)
