# Exploratory plots pathways ALS
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rio)
library(forcats)
library(ComplexHeatmap)

# Plot theme
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

## Data
mb <- readRDS("data/microbiome/imputed_microbiome_data.RDS") %>%
    mutate(across(where(is.numeric), log10))
df <- readRDS("data/microbiome/pathways.RDS") %>% filter(rownames(.) %in% mb$ID) # only samples that are in the microbiome data
head(df)
keypath <- readRDS("data/microbiome/pathwaykeys.RDS")
meta <- readRDS("data/microbiome/meta_microbiome_run1.RDS")
load <- read.csv("results/microbiome/tcam/pythonoutput/df_loadings.csv")
min5 <- ncol(load) - 5
f1 <- load %>% select(X, F1.19.99.) %>% arrange(-F1.19.99.) 
f1 <- f1[c(1:10,126:134),]
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

### Correlations and heatmap ###
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

# Calculate correlations and identify significant ones - KEEP THIS FILTER
correlations <- correlate_bugs_pathways(merged_df, bugs, pathways) %>%
    mutate(q.value = p.adjust(P.value, method = "fdr")) %>%
    left_join(keypath, by = c("Pathway" = "keys"))
sig_pathways <- correlations %>% filter(q.value < 0.05) %>% pull(Pathway) %>% unique()
correlations_filtered <- correlations %>% filter(Pathway %in% sig_pathways)

# Transform correlation data to matrix format - TRANSPOSED
corr_matrix <- correlations_filtered %>%
  select(Bug, Pathway, Correlation) %>%
  pivot_wider(names_from = Bug, values_from = Correlation) %>%
  column_to_rownames("Pathway")

# Store q-values in a matrix format for cell annotation - TRANSPOSED
qval_matrix <- correlations_filtered %>%
  select(Bug, Pathway, q.value) %>%
  pivot_wider(names_from = Bug, values_from = q.value) %>%
  column_to_rownames("Pathway")

# Get pathway explanations for better labels
pathway_labels <- correlations_filtered %>%
  select(Pathway, expl) %>%
  mutate(expl = str_remove(expl, "\\s*\\([^)]*\\)"), 
        expl = as.factor(expl)) %>% 
  mutate(expl = gsub("&beta;", "β", expl, fixed=TRUE)) %>%
  distinct() %>%
  arrange(Pathway)

# Define color scheme
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c(pal_nejm()(6)[6], "white", pal_nejm()(3)[3]))

# Create the ComplexHeatmap with transposed matrix
heatmap_bugs_pathways <- Heatmap(
  as.matrix(corr_matrix),
  name = "Correlation",
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 2),
  na_col = "grey95",  # Color for NA values
  cell_fun = function(j, i, x, y, width, height, fill) {   # Add significance markers to cells
    # Add asterisks based on q-values
    if (!is.na(corr_matrix[i, j])) {
      if (qval_matrix[i, j] < 0.001) {
        grid.text("***", x, y, gp = gpar(fontsize = 11), vjust = 0.75)
      } else if (qval_matrix[i, j] < 0.01) {
        grid.text("**", x, y, gp = gpar(fontsize = 11), vjust = 0.75)
      } else if (qval_matrix[i, j] < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 11), vjust = 0.75)
      }
    }
  },
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_dend_side = "right",
  row_labels = pathway_labels$expl[match(rownames(corr_matrix), pathway_labels$Pathway)],
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),
  row_names_rot = 0,
  column_names_rot = 45,
  # width = unit(6, "cm"), 
  # height = unit(0.25*31, "cm"),
  heatmap_legend_param = list(
    title = "Spearman\nCorrelation",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1")
  )
)
heatmap_bugs_pathways

lgd_sig = Legend(pch = c("*","**","***","****"), type = "points", 
                    labels = c("0.05", "0.01", "0.001", "<0.001"),
                    legend_gp = gpar(fontsize = 8))
heatmap_anno <- draw(heatmap_bugs_pathways, annotation_legend_list = list(lgd_sig))

# Save
# Replace your PDF save code with this:
cairo_pdf("results/pathways/heatmap_microbes_pathways_complex.pdf", width = 12, height = 10,
          family = "Helvetica")
draw(heatmap_anno)
dev.off()


### Scatterplots for Akkermansia and Porphyromonadaceae ###
merged_df <- mbsel %>%
  select(ID, Sex, Genotype, all_of(f1$X), Age_weeks) %>%
  right_join(df, by = "ID") %>%
  filter(Age_weeks == "14 weeks")
pathway_expl <- setNames(str_remove(keypath$expl, "\\s*\\([^)]*\\)"), keypath$keys)

# Akkermansia
# akk <- corr_long %>% filter(str_detect(Bug, "Akkermansia"))
# akk_df <- merged_df %>% select(ID, Genotype, Sex, `Akkermansia muciniphila`, all_of(akk$Pathway))

# plots <- list()
# for (pathway in akk$Pathway) {
#   p <- plot_correlation(akk_df, "Akkermansia muciniphila", paste0(pathway))
#   plots[[pathway]] <- p
# }
# ggarrange(plotlist = plots, ncol = 4, nrow = 5)
# ggsave("results/pathways/scatterplots_akkermansia.pdf", width = 20, height = 26)

# # Porphyromonadaceae
# por <- corr_long %>% filter(str_detect(Bug, "Porphyromonadaceae"))
# por_df <- merged_df %>% select(ID, Genotype, Sex, `Porphyromonadaceae bacterium UBA7141`, all_of(por$Pathway))

# plots <- list()
# for (pathway in por$Pathway) {
#   p <- plot_correlation(por_df, "Porphyromonadaceae bacterium UBA7141", paste0(pathway))
#   plots[[pathway]] <- p
# }
# ggarrange(plotlist = plots, ncol = 4, nrow = 2)
# ggsave("results/pathways/scatterplots_porphyromonadaceae.pdf", width = 20, height = 12)

### Sex differences in significant pathways ###
pathway_expl_wrapped <- setNames(str_wrap(pathway_expl, width = 35), names(pathway_expl))
res_box <- list()
for(a in sig_pathways){
    pathway_expl_name <- pathway_expl_wrapped[[a]]
    print(pathway_expl_name)
    print(a)
    dfpath <- merged_df %>% select(Sex, Genotype, all_of(a))
    dftdp <- dfpath %>% filter(Genotype == "TDP43")
    dfpath$path_y <- dfpath[[3]]
    comp <- list(c("Female", "Male"))
    (pl <- ggplot(data = dfpath, aes(x = Sex, y = path_y)) + 
            ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comp,
                                       hide.ns = TRUE, size = 5, #step.increase = 0.1,
                                       var.equal = FALSE) +
            geom_boxplot(aes(fill = Sex), outlier.shape = NA, 
                         width = 0.5, alpha = 0.9) +
            geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
            scale_fill_manual(guide = "none", values = ggsci::pal_nejm()(2)) +
            labs(y="log10(cpm)", x = "", title = paste0(pathway_expl_name)) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
            facet_wrap(~Genotype) +
            theme(plot.title = element_text(size = 14)) +  # Set fontsize here
            theme_Publication())
    test <- wilcox.test(dftdp[[3]] ~ dftdp$Sex)
    if(test$p.value < 0.05){
      res_box[[a]] <- pl
    } else {
      next
    }
}
length(res_box)
(paths_sexdiff <- ggarrange(plotlist = res_box, ncol = 3, nrow = 3, labels = LETTERS[2:9]))
ggsave("results/pathways/boxplots/all_pathways_sexdifferences.pdf", width = 16, height = 12)
# ggsave("results/pathways/boxplots/all_pathways_sexdifferences.svg", width = 12, height = 12)

### Assembled plot ###
# Convert the heatmap to a grob object with minimal margins
heatmap_grob <- grid.grabExpr(draw(
      heatmap_bugs_pathways, 
      annotation_legend_list = list(lgd_sig),
      padding = unit(c(15, 65, 10, 2), "mm"), # bottom, left, top, right padding
  )
)
ggsave("results/pathways/heatmap_microbes_pathways_complex.pdf", plot = as_ggplot(heatmap_grob),
       width = 12, height = 10)
ggarrange(as_ggplot(heatmap_grob), paths_sexdiff, ncol = 1, nrow = 2, heights = c(0.3, 0.35), labels = c("A", ""))
ggsave("results/pathways/assembled_plot.pdf", width = 12, height = 22)
ggsave("results/pathways/assembled_plot.pdf", width = 12, height = 22,
       device = cairo_pdf, family = "Helvetica")
# ggsave("results/pathways/assembled_plot.svg", width = 20, height = 12)
