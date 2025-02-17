# TCAM plots ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
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
tcam <- read.csv("results/pathways/tcam/df_plot_paths.csv") %>% 
  select(1:17) %>% 
  mutate(GenotypePerSex = fct_relevel(GenotypePerSex, "Female WT", after = 0L),
         GenotypePerSex = fct_relevel(GenotypePerSex, "Male WT", after = 1L),
         Genotype = fct_relevel(Genotype, "WT", after = 0L))
head(tcam)
load <- read.csv("results/pathways/tcam/df_loadings_paths.csv")
head(load)
paths <- readRDS("data/pathways.RDS")
meta <- readRDS("data/meta_microbiome.RDS")
keypath <- readRDS("data/pathwaykeys.RDS")

## TCAM plot
tcam <- tcam %>% filter(Age_ints == 6)
dim(tcam)
head(tcam)
names(tcam)

(pl <- ggplot(data = tcam, aes(x = `F1.13.4.`, y = `F2.6.13.`, 
                        color = GenotypePerSex, fill = GenotypePerSex)) +
        stat_ellipse(geom = "polygon", alpha = 0.3) +
        geom_point(size = 2) +
        ggtitle('TCAM PC1-PC2') +
        scale_color_manual(values = pal_nejm()(4)) +
        scale_fill_manual(values = pal_nejm()(4)) +
        theme_minimal() +
        labs(x=str_c('PC1 13.4%'), y=str_c('PC2 6.1%')) +
        theme_Publication() +
        facet_wrap(~ Genotype) +
        theme(legend.title = element_blank()))
ggsave("results/pathways/tcam/f1f2_scatter_mf.pdf", width = 9, height = 6)  

(pl <- ggplot(data = tcam, aes(x = `F3.4.38.`, y = `F4.3.74.`, 
                        color = GenotypePerSex, fill = GenotypePerSex)) +
        stat_ellipse(geom = "polygon", alpha = 0.3) +
        geom_point() +
        ggtitle('TCAM PC3-PC4') +
        scale_color_manual(values = pal_nejm()(4)) +
        scale_fill_manual(values = pal_nejm()(4)) +
        theme_minimal() +
        labs(x=str_c('PC3 4.4%'), y=str_c('PC4 3.7%')) +
        theme_Publication() +
        facet_wrap(~Genotype) +
        theme(legend.title = element_blank()))
ggsave("results/pathways/tcam/f3f4_scatter.pdf", width = 6, height = 6) 

## Loading of component plots
head(load)
names(load)
last <- nrow(load)
min10 <- last - 9
f1 <- load %>% select(X, F1.13.4.) %>% arrange(-F1.13.4.) %>% slice(1:10, min10:last)
f2 <- load %>% select(X, F2.6.13.) %>% arrange(-F2.6.13.) %>% slice(1:10, min10:last)
f3 <- load %>% select(X, F3.4.38.) %>% arrange(-F3.4.38.) %>% slice(1:10, min10:last)
f4 <- load %>% select(X, F4.3.74.) %>% arrange(-F4.3.74.) %>% slice(1:10, min10:last)

plot_loadings <- function(data, keypath, title) {
  df <- data %>% 
    left_join(keypath, by = c("X" = "keys")) %>%
    mutate(X = expl) %>%
    select(-expl)
  df$comp <- df[[2]]
  ggplot(df, aes(x = reorder(X, comp), y = comp)) +
    geom_bar(stat = "identity", 
        fill = c(rep("firebrick", 10), rep("royalblue",10))) +
    coord_flip() +
    labs(title = title, x = "Features", y = "Loadings") +
    theme_Publication()
}

plot_loadings(f1, keypath, "TCAM PC1")
ggsave("results/pathways/tcam/loading_pc1.pdf", width = 10, height = 7)

plot_loadings(f2, keypath, "TCAM PC2")
ggsave("results/pathways/tcam/loading_pc2.pdf", width = 10, height = 7)

plot_loadings(f3, keypath, "TCAM PC3")
ggsave("results/pathways/tcam/loading_pc3.pdf", width = 10, height = 7)

plot_loadings(f4, keypath, "TCAM PC4")
ggsave("results/pathways/tcam/loading_pc4.pdf", width = 10, height = 7)

## Lineplots
pathsel <- paths %>% select(any_of(f1$X))
head(pathsel)
pathsel$ID <- rownames(pathsel)
pathmerge <- left_join(pathsel, meta, by = "ID")

## Calculate means and standard deviations per GenotypePerSex
means_ses <- pathmerge %>%
  group_by(GenotypePerSex, Age_ints) %>%
  summarise(across(1:20, list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~sd(.x, na.rm = TRUE) / sqrt(n()))))
means_ses

# Convert means_ses to long format for plotting
means_ses_long <- means_ses %>%
  pivot_longer(cols = c(-GenotypePerSex, -Age_ints), 
                        names_to = c("path", ".value"), names_sep = "_") %>%
  filter(Age_ints < 20)

# Full names for the pathways instead of abbrev
means_ses_long <- means_ses_long %>%
  left_join(keypath, by = c("path" = "keys")) %>%
  mutate(path = expl) %>%
  select(-expl)

# Plot means with error bars representing standard deviations
ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                group = GenotypePerSex, color = GenotypePerSex)) +
  geom_line() +
  geom_point() +
  scale_color_nejm() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(title = "Means and Standard Deviations per GroupPerSex",
       x = "Age (weeks)",
       y = "Relative abundance") +
  facet_wrap(~ path, scales = "free") +
  scale_x_continuous(breaks = c(6,8,10,12,14,16,18))+
  theme_Publication() +
  theme(strip.text = element_text(size = 8))
ggsave("results/pathways/tcam/lineplots.pdf", width = 20, height = 15)
