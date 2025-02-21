# TCAM plots ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

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
tcam <- read.csv("results/microbiome/tcam/pythonoutput/df_plot.csv") %>% 
  select(1:17) %>% 
  mutate(Genotype = fct_relevel(Genotype, "WT", after = 0L),
          GenotypePerSex = fct_relevel(GenotypePerSex, "Female WT", after = 0L),
          GenotypePerSex = fct_relevel(GenotypePerSex, "Male WT", after = 1L))
head(tcam)
load <- read.csv("results/microbiome/tcam/pythonoutput/df_loadings.csv")
head(load)
mb <- readRDS("data/imputed_microbiome_data.RDS")
meta <- readRDS("data/meta_microbiome.RDS")

## TCAM plot
tcam <- tcam %>% filter(Age_ints == 6)
dim(tcam)
head(tcam)
names(tcam)
f1n <- names(tcam)[7]
f2n <- names(tcam)[8]
f3n <- names(tcam)[9]
f4n <- names(tcam)[10]

(pl <- ggplot(data = tcam, aes(x = .data[[f1n]], y = .data[[f2n]], 
                        color = Sex, fill = Sex, shape = Sex)) +
        stat_ellipse(geom = "polygon", alpha = 0.3) +
        geom_point(size = 3, alpha = 0.75) +
        ggtitle('TCAM: sex differences') +
        scale_color_manual(values = pal_nejm()(2)) +
        scale_fill_manual(values = pal_nejm()(2)) +
        theme_minimal() +
        labs(x=str_c(str_replace(f1n, "[.]", " "), "%"),
            y=str_c(str_replace(f2n, "[.]", " "), "%")) +
        theme_Publication() +
        facet_wrap(~Genotype) +
        theme(legend.title = element_blank()))
ggsave("results/microbiome/tcam/routput/f1f2_scatter.pdf", width = 6, height = 6)   

(pl <- ggplot(data = tcam, aes(x = .data[[f1n]], y = .data[[f2n]], 
                        color = Genotype, fill = Genotype, shape = Sex)) +
        #stat_ellipse(geom = "polygon", alpha = 0.3) +
        geom_point(size = 5) +
        ggtitle('TCAM') +
        scale_color_manual(values = pal_nejm()(6)[c(6,3)]) +
        scale_fill_manual(values = pal_nejm()(6)[c(6,3)]) +
        theme_minimal() +
        labs(x=str_c(str_replace(f1n, "[.]", " "), "%"),
            y=str_c(str_replace(f2n, "[.]", " "), "%")) +
        theme_Publication() +
        theme(legend.title = element_blank()))
ggsave("results/microbiome/tcam/routput/f1f2_scatter_genotype.pdf", width = 6, height = 6) 

## Loading of component plots
head(load)
names(load)[1:10]
last <- nrow(load)
min10 <- last - 9
f1 <- load %>% select(X, all_of(f1n)) %>% arrange(-.data[[f1n]]) %>% slice(1:10, min10:last)
f2 <- load %>% select(X, all_of(f2n)) %>% arrange(-.data[[f2n]]) %>% slice(1:10, min10:last)
f3 <- load %>% select(X, all_of(f3n)) %>% arrange(-.data[[f3n]]) %>% slice(1:10, min10:last)
f4 <- load %>% select(X, all_of(f4n)) %>% arrange(-.data[[f4n]]) %>% slice(1:10, min10:last)

plot_loadings <- function(data, title) {
  data$comp <- data[[2]]
  ggplot(data, aes(x = reorder(X, comp), y = comp)) +
    geom_bar(stat = "identity", 
        fill = c(rep("firebrick", 10), rep("royalblue",10))) +
    coord_flip() +
    labs(title = title, x = "Features", y = "Loadings") +
    theme_Publication()
}

plot_loadings(f1, "TCAM F1")
ggsave("results/microbiome/tcam/routput/loading_pc1.pdf", width = 8, height = 7)

plot_loadings(f2, "TCAM F2")
ggsave("results/microbiome/tcam/routput/loading_pc2.pdf", width = 8, height = 7)

plot_loadings(f3, "TCAM F3")
ggsave("results/microbiome/tcam/routput/loading_pc3.pdf", width = 8, height = 7)

plot_loadings(f4, "TCAM F4")
ggsave("results/microbiome/tcam/routput/loading_pc4.pdf", width = 8, height = 7)

# Lineplots
min5 <- last - 4
f1 <- load %>% select(X, all_of(f1n)) %>% arrange(-.data[[f1n]]) %>% slice(1:5, min5:last)

## Transform data
pseudocounts <- mb %>%  select(any_of(f1$X)) %>%
  summarise(across(everything(), ~ min(.x[.x > 0], na.rm = TRUE) / 2))
mbsel <- mb %>% ungroup(.) %>% select(any_of(f1$X), ID, Genotype, Age_ints, Sex, GenotypePerSex) %>%
  mutate( across(f1$X, ~log10(.x + pseudocounts[[cur_column()]]))
          ) %>%
  rename_at(c(f1$X), ~str_remove(.x, " \\(.*\\)$"))
head(mbsel)

## Calculate means and standard deviations per Genotype
means_ses <- mbsel %>%
  group_by(Genotype, Age_ints) %>%
  summarise(across(1:10, list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~(sd(.x, na.rm = TRUE) / sqrt(n())))))
means_ses

# Convert means_ses to long format for plotting
means_ses_long <- means_ses %>%
  pivot_longer(cols = c(-Genotype, -Age_ints), 
                        names_to = c("microbe", ".value"), names_sep = "_")
mbsel_long <- mbsel %>%
  pivot_longer(cols = c(-ID, -Genotype, -Age_ints, -GenotypePerSex, -Sex), 
                        names_to = "microbe", values_to = "value")

# Plot means with error bars representing standard deviations
ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                group = Genotype, color = Genotype)) +
  geom_line() +
  geom_point(size = 0.5) +
  geom_jitter(data = mbsel_long, aes(x = Age_ints, y = value, color = Genotype),
                width = 0.1, alpha = 0.7) +
  scale_color_manual(values = pal_nejm()(8)[c(3,6)]) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(title = "Abundance over time - TDP43 vs WT",
       x = "Age (weeks)",
       y = "log10(relative abundance)") +
  scale_x_continuous(breaks = c(6,8,10,12,14,16)) +
  facet_wrap(~ microbe, scales = "free", nrow = 2) +
  theme_Publication() +
  theme(strip.text = element_text(size = 8))
ggsave("results/microbiome/tcam/routput/lineplots_genotype_log10.pdf", width = 15, height = 8)

## Calculate means and standard deviations per Genotype
means_ses <- mbsel %>% filter(Genotype == "TDP43") %>%
  group_by(Sex, Age_ints) %>%
  summarise(across(1:10, list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~(sd(.x, na.rm = TRUE) / sqrt(n())))))
means_ses

# Convert means_ses to long format for plotting
means_ses_long <- means_ses %>%
  pivot_longer(cols = c(-Sex, -Age_ints), 
                        names_to = c("microbe", ".value"), names_sep = "_")
mbsel_long <- mbsel %>%
  pivot_longer(cols = c(-ID, -Genotype, -Age_ints, -GenotypePerSex, -Sex), 
                        names_to = "microbe", values_to = "value")

# Line plots TDP43 - sex differences
ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                group = Sex, color = Sex)) +
  geom_line() +
  geom_jitter(data = mbsel_long %>% filter(Genotype == "TDP43"), aes(x = Age_ints, y = value, 
                color = Sex),
                width = 0.1, alpha = 0.7) +
  scale_color_nejm() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(title = "Abundance over time for TDP43 - males and females",
       x = "Age (weeks)",
       y = "log10(relative abundance)") +
  scale_x_continuous(breaks = c(6,8,10,12,14,16)) +
  facet_wrap(~ microbe, scales = "free", nrow = 2) +
  theme_Publication() +
  theme(strip.text = element_text(size = 8))
ggsave("results/microbiome/tcam/routput/lineplots_tdp43_sex.pdf", width = 15, height = 8)

