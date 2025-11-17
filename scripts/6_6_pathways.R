## Pathways - microbes in human hosts
## Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Libraries
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

wrap_label <- function(label, width = 20) {
  # Insert \n every `width` characters at spaces if possible
  sapply(label, function(x) {
    if (is.na(x)) return(NA)
    paste(strwrap(x, width = width), collapse = "\n")
  })
}


## Data
mb <- readRDS("data/human_cohort/microbiome.RDS") %>% mutate(across(where(is.numeric), log10))
microb <- read.csv("results/humancohort/microbes_validate.csv")
mb_cols <- colnames(mb)
mb_subset <- mb[, microb$x, drop = FALSE]
bugs <- names(mb_subset)
mb_subset$ID <- rownames(mb_subset)

df <- readRDS("data/human_cohort/human_als_pathways.RDS") %>% filter(rownames(.) %in% rownames(mb_subset)) # only samples that are in the microbiome data
head(df)
dim(df)
pathways <- colnames(df)[1:ncol(df)-1]
pathways2 <- c("DTDPRHAMSYN-PWY", "FUC-RHAMCAT-PWY", "RHAMCAT-PWY")
keypath <- readRDS("data/human_cohort/pathwaykeys.RDS")
meta <- readRDS("data/human_cohort/metadata.RDS")

# Merge the mb dataframe with the pathways dataframe
df$ID <- rownames(df)
merged_df <- mb_subset %>%
  right_join(df, by = "ID") %>%
  select(-ID)

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

# Calculate correlations and identify significant ones
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
wrapped_labels <- wrap_label(pathway_labels$expl[match(rownames(corr_matrix), pathway_labels$Pathway)], width = 100)

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
  row_labels = wrapped_labels,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),
  row_names_rot = 0,
  column_names_rot = 45,
  width = unit(100, "mm"),
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

heatmap_grob <- grid.grabExpr(draw(
      heatmap_bugs_pathways, 
      annotation_legend_list = list(lgd_sig),
      padding = unit(c(15, 65, 10, 2), "mm"), # bottom, left, top, right padding
  )
)
ggsave("results/humancohort/heatmap_microbes_pathways_complex.pdf", plot = as_ggplot(heatmap_grob),
       width = 12, height = 10)

# Save
# Replace your PDF save code with this:
cairo_pdf("results/humancohort/heatmap_microbes_pathways.pdf", width = 15, height = 45,
          family = "Helvetica")
draw(heatmap_anno, padding = unit(c(5,5,5,20), "mm"))
dev.off()

### Rhamnose pathways
# Calculate correlations and identify significant ones
correlations <- correlate_bugs_pathways(merged_df, bugs, pathways2) %>%
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
wrapped_labels <- wrap_label(pathway_labels$expl[match(rownames(corr_matrix), pathway_labels$Pathway)], width = 35)

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
  row_labels = unname(wrapped_labels),
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

heatmap_bugs_pathways

lgd_sig = Legend(pch = c("*","**","***","****"), type = "points", 
                    labels = c("0.05", "0.01", "0.001", "<0.001"),
                    legend_gp = gpar(fontsize = 8))
heatmap_anno <- draw(heatmap_bugs_pathways, annotation_legend_list = list(lgd_sig))

# Save
cairo_pdf("results/humancohort/heatmap_microbes_pathways_rhamnose.pdf", width = 12, height = 4,
          family = "Helvetica")
draw(heatmap_anno)
dev.off()

svg("results/humancohort/heatmap_microbes_pathways_rhamnose.svg", width = 12, height = 4,
          family = "Helvetica")
draw(heatmap_anno)
dev.off()

png("results/humancohort/heatmap_microbes_pathways_rhamnose.png", width = "1200", height = "400",
          family = "Helvetica")
draw(heatmap_anno)
dev.off()

### Sex differences in significant pathways 2: facet per sex ###
merged_df <- df %>% right_join(meta, by = "ID")
pathway_expl_wrapped <- setNames(str_wrap(pathway_labels$expl, width = 35), names(pathway_labels$expl))
res_box <- list()
for(a in 1:length(sig_pathways)){
    print(a)
    pw <- sig_pathways[[a]]
    pathway_expl_name <- pathway_expl_wrapped[[a]]
    print(pathway_expl_name)
    dfpath <- merged_df %>% select(Sex, Group, all_of(pw))
    dfals <- dfpath %>% filter(Group == "ALS")
    dfpath$path_y <- dfpath[[3]]
  
    (sexdiff1 <- wilcox.test(path_y ~ Group, data = dfpath %>% filter(Sex == "Female")))
    (sexdiff2 <- wilcox.test(path_y ~ Group, data = dfpath %>% filter(Sex == "Male")))
    p_female <- format.pval(sexdiff1$p.value, digits = 2)
    p_male <- format.pval(sexdiff2$p.value, digits = 2)
    sex_diff_text <- paste0("ALS-Healthy p = ", p_female, " (women); p = ", p_male, " (men)")
  
    (pl <- ggplot(data = dfpath, aes(x = Sex, y = path_y)) + 
            ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
            geom_boxplot(aes(fill = Sex), outlier.shape = NA, 
                         width = 0.5, alpha = 0.9) +
            geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
            scale_fill_manual(guide = "none", values = ggsci::pal_nejm()(2)) +
            labs(y="log10(cpm)", x = "", title = pathway_expl_name, caption = sex_diff_text) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
            facet_wrap(~Group) +
            theme(plot.title = element_text(size = 14), plot.caption = element_text(size = 12)) +
            theme_Publication())    
    res_box[[a]] <- pl
}
length(res_box)
(paths_sexdiff <- ggarrange(plotlist = res_box, ncol = 3, nrow = 1, labels = LETTERS[2:4]))
ggsave("results/humancohort/rhamnose_pathways_sexdifferences_pval.pdf", width = 12, height = 4)

cairo_pdf("results/humancohort/rhamnose_pathways_sexdifferences_anno.pdf", width = 12, height = 4)
(paths_sexdiff <- ggarrange(plotlist = res_box, ncol = 3, nrow = 1))
dev.off()
svg(filename = "results/humancohort/rhamnose_pathways_sexdifferences_wilcox_anno.svg", 
        width = 12, height = 4)
paths_sexdiff
dev.off()
ggsave(paths_sexdiff, filename = "results/pathways/boxplots/all_pathways_sexdifferences_pval.png", 
        width = 10, height = 4)

### Sex differences in significant pathways ###
merged_df <- df %>% right_join(meta, by = "ID")
pathway_expl_wrapped <- setNames(str_wrap(pathway_labels$expl, width = 35), names(pathway_labels$expl))
res_box <- list()
for(a in 1:length(sig_pathways)){
    print(a)
    pw <- sig_pathways[[a]]
    pathway_expl_name <- pathway_expl_wrapped[[a]]
    dfpath <- merged_df %>% select(Sex, Group, all_of(pw))
    dftdp <- dfpath %>% filter(Group == "ALS")
    dfpath$path_y <- dfpath[[3]]
    comp <- list(c("Female", "Male"))
    (pl <- ggplot(data = dfpath, aes(x = Sex, y = path_y)) + 
            ggpubr::stat_compare_means(method = "t.test", label = "p.format", comparisons = comp,
                                       hide.ns = TRUE, size = 5, #step.increase = 0.1,
                                       var.equal = FALSE) +
            geom_boxplot(aes(fill = Sex), outlier.shape = NA, 
                         width = 0.5, alpha = 0.9) +
            geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
            scale_fill_manual(guide = "none", values = ggsci::pal_nejm()(2)) +
            labs(y="log10(cpm)", x = "", title = paste0(pathway_expl_name)) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
            facet_wrap(~Group) +
            theme(plot.title = element_text(size = 14)) +  # Set fontsize here
            theme_Publication())
    print(pl)
    test <- wilcox.test(dftdp[[3]] ~ dftdp$Sex)
    print(test)
  res_box[[a]] <- pl
}
length(res_box)
(paths_sexdiff <- ggarrange(plotlist = res_box, ncol = 3, nrow = 1, labels = LETTERS[2:4]))
ggsave("results/humancohort/rhamnose_pathways_sexdifferences_pval.pdf", width = 10, height = 4)
(paths_sexdiff <- ggarrange(plotlist = res_box, ncol = 3, nrow = 1))
svg(filename = "results/humancohort/rhamnose_pathways_sexdifferences_ttest.svg", 
        width = 12, height = 4)
paths_sexdiff
dev.off()
ggsave(paths_sexdiff, filename = "results/pathways/boxplots/all_pathways_sexdifferences_pval.png", 
        width = 10, height = 4)

### Assembled plot ###
# Convert the heatmap to a grob object with minimal margins
heatmap_grob <- grid.grabExpr(draw(
      heatmap_bugs_pathways, 
      annotation_legend_list = list(lgd_sig),
      padding = unit(c(15, 10, 10, 2), "mm"), # bottom, left, top, right padding
  )
)
ggarrange(as_ggplot(heatmap_grob), paths_sexdiff, ncol = 1, nrow = 2, heights = c(0.2, 0.2), labels = c("A", ""))
ggsave("results/humancohort/assembled_plot.pdf", width = 13, height = 12,
       device = cairo_pdf, family = "Helvetica")
