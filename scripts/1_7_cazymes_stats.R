# Cayman analyses - exploratory

# Library
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rio)
library(rstatix)
library(ggrepel)
library(vegan)

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
cayman[is.na(cayman)] <- 0
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