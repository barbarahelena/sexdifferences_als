# Alpha diversity plots ALS project
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ggsci)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.5, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

## Load data
mb <- readRDS("data/microbiome/microbiome.RDS")
df <- readRDS("data/microbiome/meta_microbiome.RDS")
head(mb)[1:5,1:5]
rowSums(mb)

# Diversity metrics between 4 groups
## Shannon plot
shannon <- vegan::diversity(mb, index = 'shannon')
df_shan <- data.frame(ID = names(shannon), shannon = shannon)
df_shan <- left_join(df_shan, df, by = "ID")
(plshan <- ggplot(data = df_shan %>% filter(Age_ints < 20), aes(x = as.factor(Age_ints), 
                        y = shannon, fill = Genotype,
                        group = interaction(as.factor(Age_ints), Genotype))) +
    scale_fill_manual(values = pal_nejm()(6)[c(3,6)]) +
    stat_compare_means(aes(group = interaction(as.factor(Age_ints), Genotype)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_dodge(0.75)) +
    labs(title = "Shannon index", y = "Shannon index", x="Age (weeks)") +
    theme_Publication())
#ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.svg", width = 12, height = 5)
ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.pdf", width = 7, height = 5)

(plshan <- ggplot(data = df_shan %>% filter(Age_ints < 20), aes(x = as.factor(Age_ints), 
                        y = shannon, fill = Sex,
                        group = interaction(as.factor(Age_ints), Sex))) +
    scale_fill_manual(values = pal_nejm()(2)) +
    stat_compare_means(aes(group = interaction(as.factor(Age_ints), Sex)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_dodge(0.75)) +
    labs(title = "Shannon index", y = "Shannon index", x="Age (weeks)") +
    facet_wrap(~Genotype, nrow = 2) +
    theme_Publication())
#ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.svg", width = 12, height = 5)
ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon_sex.pdf", width = 7, height = 10)

## Simpsons
simpson <- vegan::diversity(mb, index = 'simpson')
df_simp <- data.frame(ID = names(simpson), simpson = simpson)
df_simp <- left_join(df_simp, df, by = "ID")
(plsimp <- ggplot(data = df_simp, aes(x = as.factor(Age_ints), 
                        y = simpson, fill = Genotype,
                        group = interaction(as.factor(Age_ints), Genotype))) +
    scale_fill_manual(values = pal_nejm()(4)[3:4]) +
    stat_compare_means(aes(group = interaction(as.factor(Age_ints), Genotype)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_dodge(0.75)) +
    labs(title = "Simpsons index", y = "Simpsons index", x="Age (weeks)") +
    facet_wrap(~ Sex) +
    theme_Publication())
#ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.svg", width = 12, height = 5)
ggsave(plsimp, filename = "results/microbiome/alphadiversity/simpson.pdf", width = 12, height = 5)

## Species richness
specrich <- specnumber(mb)
dfspec <- data.frame(ID = names(specrich), richness = specrich)
dfspec <- left_join(dfspec, df, by = "ID")
(plrich <- ggplot(data = dfspec, aes(x = as.factor(Age_ints), 
                        y = shannon, fill = Genotype,
                        group = interaction(as.factor(Age_ints), Genotype))) +
    scale_fill_manual(values = pal_nejm()(4)[3:4]) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(aes(group = interaction(as.factor(Age_ints), Genotype)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    geom_jitter(position = position_dodge(0.75)) +
    labs(title = "Species richness", y = "Number of species", x="Age (weeks)") +
    facet_wrap(~ Sex) +
    theme_Publication())
ggsave(plrich, filename = "results/microbiome/alphadiversity/richness.pdf", width = 4, height = 5)
#ggsave(plrich, filename = "results/microbiome/alphadiversity/richness.svg", width = 4, height = 5)

# Diversity metrics between genotype
## Shannon plot
(plshan <- ggplot(data = df_shan %>% filter(!Age_ints %in% c(20,26)), 
                        aes(x = as.factor(Age_ints), 
                        y = shannon, fill = Genotype,
                        group = interaction(as.factor(Age_ints), Genotype))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_dodge(0.75)) +
    scale_fill_manual(values = pal_nejm()(4)[3:4]) +
    stat_compare_means(aes(group = interaction(as.factor(Age_ints), Genotype)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    labs(title = "Shannon index", y = "Shannon index", x="") +
    theme_Publication())
#ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.svg", width = 12, height = 5)
ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.pdf", width = 6, height = 5)

## Simpsons
(plsimp <- ggplot(data = df_simp %>% filter(!Age_ints %in% c(20,26)), 
                        aes(x = as.factor(Age_ints), 
                        y = simpson, fill = Genotype,
                        group = interaction(as.factor(Age_ints), Genotype))) +
    geom_boxplot() +
    geom_jitter(position = position_dodge(0.75)) +
    stat_compare_means(aes(group = interaction(as.factor(Age_ints), Genotype)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    scale_fill_manual(values = pal_nejm()(4)[3:4]) +
    labs(title = "Simpsons index", y = "Simpsons index", x="") +
    theme_Publication())
#ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.svg", width = 12, height = 5)
ggsave(plsimp, filename = "results/microbiome/alphadiversity/simpson.pdf", width = 6, height = 5)

# Diversity metrics between sex - TDP43
## Shannon plot
(plshan <- ggplot(data = df_shan %>% filter(Genotype == "TDP43" &
                                                !Age_ints %in% c(20,26)), 
                        aes(x = as.factor(Age_ints), 
                        y = shannon, fill = Sex,
                        group = interaction(as.factor(Age_ints), Sex))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_dodge(0.75)) +
    scale_fill_manual(values = pal_nejm()(4)[3:4]) +
    stat_compare_means(aes(group = interaction(as.factor(Age_ints), Sex)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    labs(title = "Shannon index (TDP43)", y = "Shannon index", x="") +
    theme_Publication())
#ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.svg", width = 12, height = 5)
ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon_sex_tdp43.pdf", width = 6, height = 5)

## Simpsons
(plsimp <- ggplot(data = df_simp %>% filter(Genotype == "TDP43" &
                                                !Age_ints %in% c(20,26)), 
                        aes(x = as.factor(Age_ints), 
                        y = simpson, fill = Sex,
                        group = interaction(as.factor(Age_ints), Sex))) +
    geom_boxplot() +
    geom_jitter(position = position_dodge(0.75)) +
    stat_compare_means(aes(group = interaction(as.factor(Age_ints), Sex)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    scale_fill_manual(values = pal_nejm()(4)[3:4]) +
    labs(title = "Simpsons index (TDP43)", y = "Simpsons index", x="") +
    theme_Publication())
ggsave(plsimp, filename = "results/microbiome/alphadiversity/simpson_sex_tdp43.pdf", width = 6, height = 5)

# Diversity metrics between sex - ctrl
## Shannon plot
(plshan <- ggplot(data = df_shan %>% filter(Genotype == "WT" &
                                                !Age_ints %in% c(20,26)), 
                        aes(x = as.factor(Age_ints), 
                        y = shannon, fill = Sex,
                        group = interaction(as.factor(Age_ints), Sex))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_dodge(0.75)) +
    stat_compare_means(aes(group = interaction(as.factor(Age_ints), Sex)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    scale_fill_manual(values = pal_nejm()(2)) +
    labs(title = "Shannon index (WT)", y = "Shannon index", x="") +
    theme_Publication())
ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon_sex_wt.pdf", width = 6, height = 5)

## Simpsons
(plsimp <- ggplot(data = df_simp %>% filter(Genotype == "WT" &
                                                !Age_ints %in% c(20,26)), 
                        aes(x = as.factor(Age_ints), 
                        y = simpson, fill = Sex,
                        group = interaction(as.factor(Age_ints), Sex))) +
    geom_boxplot() +
    geom_jitter(position = position_dodge(0.75)) +
    stat_compare_means(aes(group = interaction(as.factor(Age_ints), Sex)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    scale_fill_manual(values = pal_nejm()(2)) +
    labs(title = "Simpsons index (WT)", y = "Simpsons index", x="") +
    theme_Publication())
#ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.svg", width = 12, height = 5)
ggsave(plsimp, filename = "results/microbiome/alphadiversity/simpson_sex_wt.pdf", width = 6, height = 5)
