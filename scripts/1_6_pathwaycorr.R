## Pathways - microbes in Calgary cohort
## Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Libraries
## Load libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rio)
library(forcats)
library(ComplexHeatmap)
library(Cairo)

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

wrap_label <- function(label, width = 20) {
  # Insert \n every `width` characters at spaces if possible
  sapply(label, function(x) {
    if (is.na(x)) return(NA)
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

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


## Data
mb <- readRDS("data/human_cohort2/microbiome_pruned.RDS") |> 
  mutate(across(where(is.numeric), log10))
dim(mb)
microb <- read.csv("results/humancohort2/wilcoxons_microbes.csv") |> 
  filter(pval_als < 0.05) |> 
  arrange(pval_als) |> 
  slice(1:15)
microb
mb_cols <- colnames(mb)
mb_subset <- mb[, microb$mbname, drop = FALSE]
bugs <- names(mb_subset)
mb_subset$ID <- rownames(mb_subset)

df <- readRDS("data/human_cohort2/human_als_pathways.RDS") %>%
  filter(rownames(.) %in% rownames(mb_subset)) # only samples that are in the microbiome data
head(df)
dim(df)
pathways <- colnames(df)[1:ncol(df)-1]
pathways2 <- c("DTDPRHAMSYN-PWY", "RHAMCAT-PWY")
keypath <- readRDS("data/human_cohort2/pathwaykeys.RDS")
meta <- readRDS("data/human_cohort2/metadata.RDS")

# Merge the mb dataframe with the pathways dataframe
df$ID <- rownames(df)
merged_df <- mb_subset %>%
  right_join(df, by = "ID") %>%
  select(-ID)

### Correlations and heatmap ###
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
ggsave("results/humancohort2/heatmap_microbes_pathways_complex.pdf", plot = as_ggplot(heatmap_grob),
       width = 12, height = 30)

# Save with less margin
cairo_pdf("results/humancohort/heatmap_microbes_pathways.pdf", width = 15, height = 45,
          family = "Helvetica")
draw(heatmap_anno, padding = unit(c(5,5,5,20), "mm"))
dev.off()

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

wrapped_labels_top <- wrap_label(
  pathway_labels_top$expl[match(rownames(corr_matrix_top), pathway_labels_top$Pathway)],
  width = 100
)

heatmap_top <- Heatmap(
  as.matrix(corr_matrix_top),
  name = "Correlation",
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 2),
  na_col = "grey95",
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(corr_matrix_top[i, j])) {
      if (qval_matrix_top[i, j] < 0.001) {
        grid.text("***", x, y, gp = gpar(fontsize = 11), vjust = 0.75)
      } else if (qval_matrix_top[i, j] < 0.01) {
        grid.text("**", x, y, gp = gpar(fontsize = 11), vjust = 0.75)
      } else if (qval_matrix_top[i, j] < 0.05) {
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
  row_labels = wrapped_labels_top,
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

heatmap_top_grob <- grid.grabExpr(draw(
  heatmap_top,
  annotation_legend_list = list(lgd_sig),
  padding = unit(c(15, 65, 10, 2), "mm")
))
ggsave("results/humancohort2/heatmap_microbes_pathways_top15.pdf", plot = as_ggplot(heatmap_top_grob),
       width = 12, height = 9)

### Sex differences in significant pathways: Wilcoxon test loop ###
merged_df <- df %>% right_join(meta, by = "ID")

# Wilcoxon test loop for sex differences in ALS and controls
sex_diff_results <- data.frame()
for(a in 1:length(sig_pathways)){
    pw <- sig_pathways[[a]]
    dfpath <- merged_df %>% select(Sex, Group, all_of(pw))
    dfpath$path_y <- dfpath[[3]]

    # Test sex difference in ALS
    test_als <- wilcox.test(path_y ~ Sex, data = dfpath %>% filter(Group == "ALS"))

    # Test sex difference in controls
    test_control <- wilcox.test(path_y ~ Sex, data = dfpath %>% filter(Group == "Control"))

    # Store results
    sex_diff_results <- rbind(sex_diff_results, data.frame(
        Pathway = pw,
        p_als = test_als$p.value,
        p_control = test_control$p.value
    ))
}

# Apply FDR correction
sex_diff_results <- sex_diff_results %>%
    mutate(q_als = p.adjust(p_als, method = "fdr"),
           q_control = p.adjust(p_control, method = "fdr"))

# View results
print(sex_diff_results)

# Filter pathways significant in ALS after FDR correction
sig_als_sexdiff <- sex_diff_results %>%
    filter(q_als < 0.05) %>%
    pull(Pathway)

print(paste0("Number of pathways with significant sex differences in ALS (FDR < 0.05): ", length(sig_als_sexdiff)))

# Save results
write.csv(sex_diff_results, "results/humancohort2/pathways_sexdiff_wilcoxon.csv", row.names = FALSE)

### Create plots only for pathways significant in ALS after FDR correction ###
pathway_expl_wrapped <- setNames(str_wrap(pathway_labels$expl, width = 35), pathway_labels$Pathway)
res_box <- list()
plot_counter <- 1

for(pw in sig_als_sexdiff){
    print(paste0("Plotting pathway: ", pw))
    pathway_expl_name <- pathway_expl_wrapped[[pw]]

    dfpath <- merged_df %>% select(Sex, Group, all_of(pw))
    dfpath$path_y <- dfpath[[3]]

    # Get p-values from results
    pval_als <- sex_diff_results %>% filter(Pathway == pw) %>% pull(p_als)
    pval_control <- sex_diff_results %>% filter(Pathway == pw) %>% pull(p_control)
    qval_als <- sex_diff_results %>% filter(Pathway == pw) %>% pull(q_als)

    caption_text <- paste0("ALS sex diff: p = ", format.pval(pval_als, digits = 2),
                          ", q = ", format.pval(qval_als, digits = 2),
                          "\nControl sex diff: p = ", format.pval(pval_control, digits = 2))

    (pl <- ggplot(data = dfpath, aes(x = Sex, y = path_y)) +
            ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
            geom_boxplot(aes(fill = Sex), outlier.shape = NA,
                         width = 0.5, alpha = 0.9) +
            geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
            scale_fill_manual(guide = "none", values = ggsci::pal_nejm()(2)) +
            labs(y="log10(cpm)", x = "", title = pathway_expl_name, caption = caption_text) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
            facet_wrap(~Group) +
            theme(plot.title = element_text(size = 14), plot.caption = element_text(size = 10)) +
            theme_Publication())
    print(pl)
    res_box[[plot_counter]] <- pl
    plot_counter <- plot_counter + 1
}

# Save plots if there are any
if(length(res_box) > 0){
    # Determine layout based on number of plots
    ncol_plot <- min(3, length(res_box))
    nrow_plot <- ceiling(length(res_box) / ncol_plot)

    (paths_sexdiff <- ggarrange(plotlist = res_box, ncol = ncol_plot, nrow = nrow_plot,
                                labels = LETTERS[1:length(res_box)]))

    ggsave("results/humancohort2/rhamnose_pathways_sexdifferences_fdr_sig.pdf",
           width = 4*ncol_plot, height = 4*nrow_plot)

    cairo_pdf("results/humancohort2/rhamnose_pathways_sexdifferences_fdr_sig_cairo.pdf",
              width = 4*ncol_plot, height = 4*nrow_plot, family = "Helvetica")
    print(paths_sexdiff)
    dev.off()

    svg(filename = "results/humancohort2/rhamnose_pathways_sexdifferences_fdr_sig.svg",
        width = 4*ncol_plot, height = 4*nrow_plot)
    print(paths_sexdiff)
    dev.off()
} else {
    print("No pathways significant after FDR correction in ALS")
}



### Rhamnose pathways
# Recreate merged_df with bugs + pathways for correlations
merged_df_full <- mb_subset %>%
  right_join(df, by = "ID") %>%
  select(-ID)

# Calculate correlations and identify significant ones
correlations <- correlate_bugs_pathways(merged_df_full, bugs, pathways2) %>%
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
cairo_pdf("results/humancohort2/heatmap_microbes_pathways_rhamnose.pdf", width = 12, height = 4,
          family = "Helvetica")
draw(heatmap_anno)
dev.off()

# svg("results/humancohort/heatmap_microbes_pathways_rhamnose.svg", width = 12, height = 4,
#           family = "Helvetica")
# draw(heatmap_anno)
# dev.off()

# png("results/humancohort/heatmap_microbes_pathways_rhamnose.png", width = "1200", height = "400",
#           family = "Helvetica")
# draw(heatmap_anno)
# dev.off()

### Sex differences in rhamnose pathways ###
merged_df_rham <- df %>% right_join(meta, by = "ID")

sex_diff_rhamnose <- data.frame()
for (pw in pathways2) {
  dfpath <- merged_df_rham %>% select(Sex, Group, all_of(pw))
  dfpath$path_y <- dfpath[[3]]

  test_als <- wilcox.test(path_y ~ Sex, data = dfpath %>% filter(Group == "ALS"))
  test_control <- wilcox.test(path_y ~ Sex, data = dfpath %>% filter(Group == "Control"))

  sex_diff_rhamnose <- rbind(sex_diff_rhamnose, data.frame(
    Pathway = pw,
    p_als = test_als$p.value,
    p_control = test_control$p.value
  ))
}

sex_diff_rhamnose <- sex_diff_rhamnose %>%
  left_join(keypath, by = c("Pathway" = "keys"))
print(sex_diff_rhamnose)
write.csv(sex_diff_rhamnose, "results/humancohort2/rhamnose_pathways_sexdiff.csv", row.names = FALSE)

# Boxplots for rhamnose pathways
rham_plots <- list()
for (i in seq_along(pathways2)) {
  pw <- pathways2[i]
  dfpath <- merged_df_rham %>% select(Sex, Group, all_of(pw))
  dfpath$path_y <- dfpath[[3]]

  pval_als <- sex_diff_rhamnose %>% filter(Pathway == pw) %>% pull(p_als)
  pval_control <- sex_diff_rhamnose %>% filter(Pathway == pw) %>% pull(p_control)
  pw_label <- sex_diff_rhamnose %>% filter(Pathway == pw) %>% pull(expl)
  if (length(pw_label) == 0 || is.na(pw_label)) pw_label <- pw

  caption_text <- paste0("ALS sex diff: p = ", format.pval(pval_als, digits = 2),
                         "\nControl sex diff: p = ", format.pval(pval_control, digits = 2))

  rham_plots[[i]] <- ggplot(data = dfpath, aes(x = Sex, y = path_y)) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
    geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.5, alpha = 0.9) +
    geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
    scale_fill_manual(guide = "none", values = ggsci::pal_nejm()(2)) +
    labs(y = "log10(cpm)", x = "", title = str_wrap(pw_label, width = 35), caption = caption_text) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    facet_wrap(~Group) +
    theme_Publication() +
    theme(plot.caption = element_text(size = 10))
}

(rham_sexdiff_plot <- ggarrange(plotlist = rham_plots, ncol = 2, nrow = 1, labels = LETTERS[1:length(rham_plots)]))
ggsave("results/humancohort2/rhamnose_pathways_sexdifferences.pdf", rham_sexdiff_plot,
       width = 8, height = 5)