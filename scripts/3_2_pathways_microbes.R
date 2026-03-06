# Functional pathway analyses ALS
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

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
        grid.text("***", x, y, gp = gpar(fontsize = 16), vjust = 0.75)
      } else if (qval_matrix[i, j] < 0.01) {
        grid.text("**", x, y, gp = gpar(fontsize = 16), vjust = 0.75)
      } else if (qval_matrix[i, j] < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 16), vjust = 0.75)
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

lgd_sig = Legend(title = "q-value", pch = c("*","**","***","****"), type = "points", 
                    labels = c("<0.05", "<0.01", "<0.001", "<0.001"),
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

### Sex differences in significant pathways ###
pathway_expl_wrapped <- setNames(str_wrap(pathway_expl, width = 35), names(pathway_expl))

# Pass 1: collect Sex x Genotype interaction p-values across all 43 pathways
sexdiff_pvals <- sapply(sig_pathways, function(a) {
    dfpath <- merged_df %>% select(Sex, Genotype, all_of(a))
    dfpath$path_y <- dfpath[[3]]
    fit <- lm(path_y ~ Sex * Genotype, data = dfpath)
    anova(fit)["Sex:Genotype", "Pr(>F)"]
})
sexdiff_qvals <- p.adjust(sexdiff_pvals, method = "fdr")

# Pass 2: build plots only for pathways with BH q < 0.05
res_box <- list()
for(a in sig_pathways){
    if(sexdiff_qvals[a] >= 0.05) next

    pathway_expl_name <- pathway_expl_wrapped[[a]]
    print(pathway_expl_name)
    print(a)
    dfpath <- merged_df %>% select(Sex, Genotype, all_of(a))
    dfpath$path_y <- dfpath[[3]]

    fit <- lm(path_y ~ Sex * Genotype, data = dfpath)
    p_interaction <- anova(fit)["Sex:Genotype", "Pr(>F)"]
    q_interaction <- sexdiff_qvals[a]
    sex_diff_text <- paste0("Sex \u00d7 Genotype: p = ", format.pval(p_interaction, digits = 2),
                            ", q = ", format.pval(q_interaction, digits = 2), " (BH)")

    (pl <- ggplot(data = dfpath, aes(x = Sex, y = path_y)) +
            geom_boxplot(aes(fill = Sex), outlier.shape = NA,
                         width = 0.5, alpha = 0.9) +
            geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
            scale_fill_manual(guide = "none", values = ggsci::pal_nejm()(2)) +
            labs(y="log10(cpm)", x = "", title = paste0(pathway_expl_name), caption = sex_diff_text) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
            facet_wrap(~Genotype) +
            theme(plot.title = element_text(size = 14), plot.caption = element_text(size = 12)) +
            theme_Publication())
    res_box[[a]] <- pl
}
length(res_box)
res_box[[1]]
(paths_sexdiff <- ggarrange(plotlist = res_box, ncol = 3, nrow = 3, labels = LETTERS[2:10]))
ggsave("results/pathways/sexdiff_qval_interactions.pdf", width = 5, height = 4)
# ggsave("results/pathways/boxplots/all_pathways_sexdifferences.svg", width = 12, height = 12)

# Supplement: boxplots for heatmap pathways with nominal Sex x Genotype interaction (p < 0.05)
res_box_all <- list()
for(a in sig_pathways){
    pathway_expl_name <- pathway_expl_wrapped[[a]]
    dfpath <- merged_df %>% select(Sex, Genotype, all_of(a))
    dfpath$path_y <- dfpath[[3]]

    fit <- lm(path_y ~ Sex * Genotype, data = dfpath)
    p_interaction <- anova(fit)["Sex:Genotype", "Pr(>F)"]
    if(p_interaction >= 0.05) next
    q_interaction <- sexdiff_qvals[a]
    sex_diff_text <- paste0("Sex \u00d7 Genotype: p = ", format.pval(p_interaction, digits = 2),
                            ", q = ", format.pval(q_interaction, digits = 2), " (BH)")

    pl <- ggplot(data = dfpath, aes(x = Sex, y = path_y)) +
            geom_boxplot(aes(fill = Sex), outlier.shape = NA,
                         width = 0.5, alpha = 0.9) +
            geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
            scale_fill_manual(guide = "none", values = ggsci::pal_nejm()(2)) +
            labs(y="log10(cpm)", x = "", title = pathway_expl_name, caption = sex_diff_text) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
            facet_wrap(~Genotype) +
            theme(plot.title = element_text(size = 14), plot.caption = element_text(size = 12)) +
            theme_Publication()
    res_box_all[[a]] <- pl
}
length(res_box_all)
ggexport(plotlist = res_box_all,
         filename = "results/pathways/supplement_all_heatmap_pathways_boxplots.pdf",
         nrow = 4, ncol = 3, width = 10, height = 14)

### Assembled plot ###
# Convert the heatmap to a grob object with minimal margins
heatmap_grob_pathway <- grid.grabExpr(draw(
      heatmap_bugs_pathways, 
      annotation_legend_list = list(lgd_sig),
      padding = unit(c(15, 65, 10, 2), "mm"), # bottom, left, top, right padding
  )
)
ggsave("results/pathways/heatmap_microbes_pathways_complex.pdf", plot = as_ggplot(heatmap_grob_pathway),
       width = 12, height = 10)

# ggsave("results/pathways/assembled_plot.svg", width = 20, height = 12)

### Top 15 strongest correlations heatmap ###
n_top <- 15
top_pathways <- correlations_filtered %>%
  group_by(Pathway) %>%
  summarise(max_abs_cor = max(abs(Correlation), na.rm = TRUE)) %>%
  arrange(desc(max_abs_cor)) %>%
  slice(1:n_top) %>%
  pull(Pathway)

corr_top <- correlations_filtered %>% filter(Pathway %in% top_pathways)

corr_matrix_top <- corr_top %>%
  select(Bug, Pathway, Correlation) %>%
  pivot_wider(names_from = Bug, values_from = Correlation) %>%
  column_to_rownames("Pathway")

qval_matrix_top <- corr_top %>%
  select(Bug, Pathway, q.value) %>%
  pivot_wider(names_from = Bug, values_from = q.value) %>%
  column_to_rownames("Pathway")

pathway_labels_top <- corr_top %>%
  select(Pathway, expl) %>%
  mutate(expl = str_remove(expl, "\\s*\\([^)]*\\)"),
         expl = as.factor(expl)) %>%
  mutate(expl = gsub("&beta;", "\u03b2", expl, fixed = TRUE)) %>%
  distinct() %>%
  arrange(Pathway)

wrapped_labels_top <- pathway_labels_top$expl[match(rownames(corr_matrix_top), pathway_labels_top$Pathway)]

heatmap_top <- Heatmap(
  as.matrix(corr_matrix_top),
  name = "Correlation",
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 2),
  na_col = "grey95",
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(corr_matrix_top[i, j])) {
      if (qval_matrix_top[i, j] < 0.001) {
        grid.text("***", x, y, gp = gpar(fontsize = 16), vjust = 0.75)
      } else if (qval_matrix_top[i, j] < 0.01) {
        grid.text("**", x, y, gp = gpar(fontsize = 16), vjust = 0.75)
      } else if (qval_matrix_top[i, j] < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 16), vjust = 0.75)
      }
    }
  },
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_dend_side = "right",
  row_labels = wrapped_labels_top,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),
  row_names_rot = 0,
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = "Spearman\nCorrelation",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1")
  )
)

heatmap_top_grob <- grid.grabExpr(draw(
  heatmap_top,
  annotation_legend_list = list(lgd_sig),
  padding = unit(c(15, 65, 10, 2), "mm")
))
ggsave("results/pathways/heatmap_microbes_pathways_top15.pdf", plot = as_ggplot(heatmap_top_grob),
       width = 12, height = 9)

