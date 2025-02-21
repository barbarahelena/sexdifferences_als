# Exploratory plots pathways ALS
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)
library(rio)
library(forcats)

theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(family = 'Helvetica'),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
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

## Data
dim(meta)
mb <- readRDS("data/imputed_microbiome_data.RDS") %>%
    mutate(across(where(is.numeric), log10))
df <- readRDS("data/pathways.RDS") %>% filter(rownames(.) %in% mb$ID)
head(df)
keypath <- readRDS("data/pathwaykeys.RDS")
meta <- readRDS("data/meta_microbiome_run1.RDS")
load <- read.csv("results/microbiome/tcam/pythonoutput/df_loadings.csv")
min5 <- ncol(load) - 5
f1 <- load %>% select(X, F1.19.99.) %>% arrange(-F1.19.99.) %>% slice(c(1:5,130:134))
f1
mbsel <- mb %>% select(all_of(f1$X), ID, Age_weeks, Genotype, Sex)

bugs <- f1$X
pathways <- colnames(df)[1:ncol(df)-1]

# Merge the mb dataframe with the pathways dataframe
df$ID <- rownames(df)
merged_df <- mbsel %>%
  select(ID, all_of(f1$X), Age_weeks) %>%
  right_join(df, by = "ID") %>%
  filter(Age_weeks == "14 weeks") %>%
    select(-ID, -Age_weeks)

# Function to calculate correlations and identify significant ones
correlate_bugs_pathways <- function(merged_df, bugs, pathways) {
  results <- data.frame()
  for (bug in bugs) {
    for (pathway in pathways) {
      cor_test <- cor.test(merged_df[[bug]], merged_df[[pathway]], method = "spearman")
      results <- rbind(results, data.frame(Bug = bug, Pathway = pathway, 
                  Correlation = cor_test$estimate, P.value = cor_test$p.value))
    }
  }
  return(results)
}

# Calculate correlations and identify significant ones
correlations <- correlate_bugs_pathways(merged_df, bugs, pathways) %>%
    mutate(q.value = p.adjust(P.value, method = "fdr")) %>%
    left_join(keypath, by = c("Pathway" = "keys")) %>%
    filter(q.value < 0.05)

corr_long <- correlations %>% pivot_longer(cols = Correlation, names_to = "Metric", values_to = "value")
ggplot(corr_long, aes(x = Bug, y = expl, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "royalblue", high = "firebrick", mid = "white", midpoint = 0, 
        limit = c(-1, 1), space = "Lab", name = "Correlation") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1) 
          ) +
    labs(title = "Correlations between microbes and pathways", x = "Microbe", y = "Pathway")
ggsave("results/pathways/heatmap_microbes_pathways.pdf", width = 10, height = 12)

# Akkermansia
akk <- corr_long %>% filter(str_detect(Bug, "Akkermansia"))
akk
head(merged_df)
akk_df <- mbsel %>%
  select(ID, Genotype, Sex, all_of(f1$X), Age_weeks) %>%
  right_join(df, by = "ID") %>%
  filter(Age_weeks == "14 weeks") %>%
  select(ID, Genotype, Sex, `Akkermansia muciniphila`, all_of(akk$Pathway))
pathway_expl <- setNames(keypath$expl, keypath$keys)

# Function for scatterplots
plot_correlation <- function(data, microbe, pathway) {
  microbe_expl <- strtrim(microbe, 23)
  pathway_expl_name <- pathway_expl[[pathway]]
  ggplot(data, aes(x = .data[[microbe]], y = .data[[pathway]],)) +
    geom_point(size = 3, alpha = 0.8, aes(color = Genotype, shape = Sex)) +
    scale_color_manual(values = pal_nejm()(6)[c(3,6)]) +
    geom_smooth(method = "lm", se = FALSE, color = "grey5") +
    stat_cor(method = "spearman", label.y = max(data[[pathway]])*1.1) +
    labs(x = str_c("log10(rel abundance) ", microbe_expl), 
          y = pathway_expl_name) +
    theme_Publication()
}

plots <- list()
for (pathway in akk$Pathway) {
  p <- plot_correlation(akk_df, "Akkermansia muciniphila", paste0(pathway))
  plots[[pathway]] <- p
}
ggarrange(plotlist = plots, ncol = 4, nrow = 5)
ggsave("results/pathways/scatterplots_akkermansia.pdf", width = 20, height = 26)

# Porphyromonadaceae
por <- corr_long %>% filter(str_detect(Bug, "Porphyromonadaceae"))
por
por_df <- mbsel %>%
  select(ID, Genotype, Sex, all_of(f1$X), Age_weeks) %>%
  right_join(df, by = "ID") %>%
  filter(Age_weeks == "14 weeks") %>%
  select(ID, Genotype, Sex, `Porphyromonadaceae bacterium UBA7141`, all_of(por$Pathway))
pathway_expl <- setNames(keypath$expl, keypath$keys)

plots <- list()
for (pathway in por$Pathway) {
  p <- plot_correlation(por_df, "Porphyromonadaceae bacterium UBA7141", paste0(pathway))
  plots[[pathway]] <- p
}
ggarrange(plotlist = plots, ncol = 4, nrow = 2)
ggsave("results/pathways/scatterplots_porphyromonadaceae.pdf", width = 20, height = 12)

# Group differences
res_box <- list()
for(a in akk$Pathway){
    pathway_expl_name <- pathway_expl[[a]]
    print(pathway_expl_name)
    print(a)
    dfpath <- akk_df %>% select(Sex, Genotype, all_of(a))
    dfpath$path_y <- dfpath[[3]]
    comp <- list(c("Female", "Male"))
    (pl <- ggplot(data = dfpath, aes(x = Sex, y = path_y)) + 
            ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comp,
                                       hide.ns = TRUE, size = 5, #step.increase = 0.1,
                                       var.equal = FALSE) +
            geom_boxplot(aes(fill = Sex), outlier.shape = NA, 
                         width = 0.5, alpha = 0.9) +
            geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
            scale_fill_manual(guide = "none", values = pal_nejm()(2)) +
            labs(y=paste0(pathway_expl_name), x = "") +
            facet_wrap(~Genotype) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    ggsave(str_c("results/pathways/boxplots/", a, ".pdf"), width = 4, height = 5.5, device = "pdf")
    ggsave(str_c("results/pathways/boxplots/", a, ".svg"), width = 4, height = 5.5, device = "svg")
    res_box[[a]] <- pl
}

res_box <- list()
for(a in por$Pathway){
    pathway_expl_name <- pathway_expl[[a]]
    print(pathway_expl_name)
    print(a)
    dfpath <- por_df %>% select(Sex, Genotype, all_of(a))
    dfpath$path_y <- dfpath[[3]]
    comp <- list(c("Female", "Male"))
    (pl <- ggplot(data = dfpath, aes(x = Sex, y = path_y)) + 
            ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comp,
                                       hide.ns = TRUE, size = 5, #step.increase = 0.1,
                                       var.equal = FALSE) +
            geom_boxplot(aes(fill = Sex), outlier.shape = NA, 
                         width = 0.5, alpha = 0.9) +
            geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
            scale_fill_manual(guide = "none", values = pal_nejm()(2)) +
            labs(y=paste0(pathway_expl_name), x = "") +
            facet_wrap(~Genotype) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    ggsave(str_c("results/pathways/boxplots/", a, ".pdf"), width = 4, height = 5.5, device = "pdf")
    ggsave(str_c("results/pathways/boxplots/", a, ".svg"), width = 4, height = 5.5, device = "svg")
    res_box[[a]] <- pl
}
