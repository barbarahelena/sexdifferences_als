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
tcam <- read.csv("results/microbiome/tcam/df_plot.csv") %>% 
  select(1:17) %>% 
  mutate(Genotype = fct_relevel(Genotype, "WT", after = 0L),
          GenotypePerSex = fct_relevel(GenotypePerSex, "Female WT", after = 0L),
          GenotypePerSex = fct_relevel(GenotypePerSex, "Male WT", after = 1L))
head(tcam)
load <- read.csv("results/microbiome/tcam/df_loadings.csv")
head(load)
mb <- readRDS("data/microbiome.RDS")
meta <- readRDS("data/meta_microbiome.RDS")

## TCAM plot
tcam <- tcam %>% filter(Age_ints == 6)
dim(tcam)
head(tcam)
names(tcam)

(pl <- ggplot(data = tcam, aes(x = `F1.19.38.`, y = `F2.8.21.`, 
                        color = GenotypePerSex, fill = GenotypePerSex)) +
        stat_ellipse(geom = "polygon", alpha = 0.3) +
        geom_point() +
        ggtitle('TCAM') +
        scale_color_manual(values = pal_nejm()(4)) +
        scale_fill_manual(values = pal_nejm()(4)) +
        theme_minimal() +
        labs(x=str_c('F1 19.4%'), y=str_c('F2 8.2%')) +
        theme_Publication() +
        facet_wrap(~Genotype) +
        theme(legend.title = element_blank()))
ggsave("results/microbiome/tcam/f1f2_scatter.pdf", width = 6, height = 6)      

(pl <- ggplot(data = tcam, aes(x = `F3.6.22.`, y = `F4.5.53.`, 
                        color = GenotypePerSex, fill = GenotypePerSex)) +
        stat_ellipse(geom = "polygon", alpha = 0.3) +
        geom_point() +
        ggtitle('TCAM') +
        scale_color_manual(values = pal_nejm()(4)) +
        scale_fill_manual(values = pal_nejm()(4)) +
        theme_minimal() +
        labs(x=str_c('F3 6.2%'), y=str_c('F4 5.5%')) +
        theme_Publication() +
        facet_wrap(~Genotype) +
        theme(legend.title = element_blank()))
ggsave("results/microbiome/tcam/f3f4_scatter.pdf", width = 6, height = 6) 


## Loading of component plots
head(load)
names(load)[1:10]
last <- nrow(load)
min10 <- last - 9
f1 <- load %>% select(X, F1.19.38.) %>% arrange(-F1.19.38.) %>% slice(1:10, min10:last)
f2 <- load %>% select(X, F2.8.21.) %>% arrange(-F2.8.21.) %>% slice(1:10, min10:last)
f3 <- load %>% select(X, F3.6.22.) %>% arrange(-F3.6.22.) %>% slice(1:10, min10:last)
f4 <- load %>% select(X, F4.5.53.) %>% arrange(-F4.5.53.) %>% slice(1:10, min10:last)

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
ggsave("results/microbiome/tcam/loading_pc1.pdf", width = 6, height = 7)

plot_loadings(f2, "TCAM F2")
ggsave("results/microbiome/tcam/loading_pc2.pdf", width = 6, height = 7)

plot_loadings(f3, "TCAM F3")
ggsave("results/microbiome/tcam/loading_pc3.pdf", width = 6, height = 7)

plot_loadings(f4, "TCAM F4")
ggsave("results/microbiome/tcam/loading_pc4.pdf", width = 6, height = 7)

## Lineplots
mbsel <- mb %>% select(any_of(f1$X))
head(mbsel)
mbsel$ID <- rownames(mbsel)
mbmerge <- left_join(mbsel, meta, by = "ID")

## Calculate means and standard deviations per GenotypePerSex
means_ses <- mbmerge %>%
  group_by(GenotypePerSex, Age_ints) %>%
  summarise(across(1:20, list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~sd(.x, na.rm = TRUE) / sqrt(n()))))
means_ses

# Convert means_ses to long format for plotting
means_ses_long <- means_ses %>%
  pivot_longer(cols = c(-GenotypePerSex, -Age_ints), 
                        names_to = c("microbe", ".value"), names_sep = "_") %>%
  filter(Age_ints < 20)

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
  facet_wrap(~ microbe, scales = "free") +
  theme_Publication() +
  theme(strip.text = element_text(size = 8))
ggsave("results/microbiome/tcam/lineplots.pdf", width = 20, height = 15)

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
  facet_wrap(~ microbe, scales = "free") +
  theme_Publication() +
  theme(strip.text = element_text(size = 8))
ggsave("results/microbiome/tcam/lineplots_genotype.pdf", width = 20, height = 15)

# Create a summary table counting the number of samples per group per timepoint
sample_counts <- meta %>%
  group_by(GenotypePerSex, Age_ints) %>%
  summarise(count = n()) %>%
  ungroup()

# Plot the bar plot
ggplot(sample_counts, aes(x = Age_ints, y = count, fill = GenotypePerSex)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of samples",
       x = "Age (weeks)",
       y = "Number of samples") +
  scale_fill_nejm() +
  scale_x_continuous(n.breaks = 16) +
  theme_Publication() +
  theme(legend.title = element_blank())
ggsave("results/microbiome/sample_counts_barplot.pdf", width = 10, height = 6)
