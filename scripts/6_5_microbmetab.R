## Welch t tests for metabolomics sex differences
## Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Libraries
library(tidyverse)
library(ggpubr)
library(ggsci)
library(Hmisc)
library(ComplexHeatmap)
library(circlize)

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

# Function to calculate correlation matrices
calculate_corr_matrices <- function(met_data, mb_data, met_names, mb_names) {
  corr_mat <- matrix(NA, nrow = length(met_names), ncol = length(mb_names))
  pval_mat <- matrix(NA, nrow = length(met_names), ncol = length(mb_names))
  
  for (i in 1:length(met_names)) {
    for (j in 1:length(mb_names)) {
      test <- cor.test(met_data[, met_names[i]], mb_data[, mb_names[j]], 
                      method = "spearman")
      corr_mat[i, j] <- test$estimate
      pval_mat[i, j] <- test$p.value
    }
  }
  
  rownames(corr_mat) <- met_names
  colnames(corr_mat) <- mb_names
  rownames(pval_mat) <- met_names
  colnames(pval_mat) <- mb_names
  
  pval_adj_mat <- matrix(p.adjust(pval_mat, method = "BH"), 
                        nrow = nrow(pval_mat), ncol = ncol(pval_mat))
  rownames(pval_adj_mat) <- rownames(pval_mat)
  colnames(pval_adj_mat) <- colnames(pval_mat)
  
  return(list(corr = corr_mat, pval = pval_mat, pval_adj = pval_adj_mat))
}

# Update the group-specific heatmap function to use clustering
create_group_heatmap <- function(corr_matrix, pval_adj_matrix, group_name) {
  Heatmap(
    corr_matrix,
    name = "Correlation",
    col = col_fun,  # Using the same NEJM color palette
    rect_gp = gpar(col = "white", lwd = 1),
    cell_fun = function(j, i, x, y, width, height, fill) {
      corr_val <- corr_matrix[i, j]
      p_val <- pval_adj_matrix[i, j]
      
      grid.text(sprintf("%.2f", corr_val), x, y, gp = gpar(fontsize = 10))
      
      if (p_val < 0.05) {
        stars <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", "*"))
        grid.text(stars, x, y + unit(0.3, "cm"), gp = gpar(fontsize = 12, fontface = "bold"))
      }
    },
    # Use clustering instead of alphabetical ordering
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    row_names_gp = gpar(fontsize = 10),
    show_column_dend = TRUE,
    show_row_dend = TRUE,
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title = paste(group_name, "Microbe-Metabolite Correlations"),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_title = "Metabolites",
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Correlation",
      at = c(-1, -0.5, 0, 0.5, 1),
      labels = c("-1", "-0.5", "0", "0.5", "1")
    )
  )
}

# Create a function to generate scatter plots for metabolite-microbe pairs
create_correlation_plot <- function(metabolite_idx, microbe_idx, 
                                   met_data, mb_data, corr_val, p_val,
                                   meta_data) {
  metabolite <- rownames(corr_matrix)[metabolite_idx]
  microbe <- colnames(corr_matrix)[microbe_idx]
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    Metabolite = met_data[, metabolite],
    Microbe = mb_data[, microbe],
    Group = meta_data$Group
  )
  
  # Format the correlation and p-value for the plot
  corr_text <- sprintf("Spearman rho = %.2f", corr_val)
  p_text <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
  
  # Create the scatter plot with regression line
  p <- ggplot(plot_data, aes(x = log10(Metabolite), y = log10(Microbe+0.01))) +
    geom_point(aes(color = Group), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    stat_cor(method = "spearman",
         label.x.npc = "left", 
         label.y.npc = "top",
         p.accuracy = 0.001,
         r.accuracy = 0.01,
         alternative = "two.sided") +
    scale_color_manual(values = pal_nejm()(6)[c(3,6)]) +
    labs(
      title = paste(microbe, "-", metabolite),
      x = paste0("log10(",metabolite,")"),
      y = paste0("log10(",microbe,")")
    ) +
    theme_Publication() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 10, face = "italic")
    )
  
  return(p)
}

## Data
meta <- readRDS("data/human_cohort/metadata.RDS")
met <- readRDS("data/human_cohort/metabolomics_human_cohort.RDS")
mb <- readRDS("data/human_cohort/microbiome.RDS")

selmet <- read.csv("results/humancohort/metabolites_validation.csv")
microb <- read.csv("results/humancohort/microbes_validate.csv")
selected_metabolites <- selmet$metabolite
selected_microbes <- microb$x
met_cols <- colnames(met)
mb_cols <- colnames(mb)
met_subset <- met[, selected_metabolites, drop = FALSE]
mb_subset <- mb[, selected_microbes, drop = FALSE]

# Find common samples between metabolome and microbiome data
common_samples <- intersect(rownames(met_subset), rownames(mb_subset))
met_subset <- met_subset[common_samples, ]
mb_subset <- mb_subset[common_samples, ]

# Add metadata for these samples
meta_subset <- meta %>% filter(ID %in% common_samples)

# Calculate Spearman correlations
corr_matrix <- matrix(NA, nrow = length(selected_metabolites), 
                     ncol = length(selected_microbes))
pval_matrix <- matrix(NA, nrow = length(selected_metabolites), 
                     ncol = length(selected_microbes))

# Calculate correlations and p-values
for (i in 1:length(selected_metabolites)) {
  for (j in 1:length(selected_microbes)) {
    test <- cor.test(met_subset[, selected_metabolites[i]], 
                    mb_subset[, selected_microbes[j]], 
                    method = "spearman")
    corr_matrix[i, j] <- test$estimate
    pval_matrix[i, j] <- test$p.value
  }
}

# Set row and column names
rownames(corr_matrix) <- selected_metabolites
colnames(corr_matrix) <- selected_microbes
rownames(pval_matrix) <- selected_metabolites
colnames(pval_matrix) <- selected_microbes

# Adjust p-values for multiple testing
pval_adj_matrix <- matrix(p.adjust(pval_matrix, method = "BH"), 
                         nrow = nrow(pval_matrix),
                         ncol = ncol(pval_matrix))
rownames(pval_adj_matrix) <- rownames(pval_matrix)
colnames(pval_adj_matrix) <- colnames(pval_matrix)

# Use NEJM color palette instead of blue-white-red
# Extract colors from NEJM palette for negative and positive correlations
nejm_colors <- pal_nejm()(8)
neg_color <- nejm_colors[6]  # Dark blue 
pos_color <- nejm_colors[3]  # Orange/red
col_fun <- colorRamp2(c(-1, 0, 1), c(neg_color, "white", pos_color))

# Create the updated heatmap with clustering
heatmap <- Heatmap(
  corr_matrix,
  name = "Correlation",
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    corr_val <- corr_matrix[i, j]
    p_val <- pval_adj_matrix[i, j]
    
    # Add correlation value
    grid.text(sprintf("%.2f", corr_val), x, y, gp = gpar(fontsize = 10))
    
    # Add asterisks for significance
    if (p_val < 0.05) {
      stars <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", "*"))
      grid.text(stars, x, y + unit(0.3, "cm"), gp = gpar(fontsize = 12, fontface = "bold"))
    }
  },
  # Use clustering instead of alphabetical ordering
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  row_names_gp = gpar(fontsize = 10),
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45,
  column_title = "Microbes",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title = "Metabolites",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Correlation",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1")
  )
)

## Create a legend for significance
significance_legend <- Legend(pch = c("*", "**", "***"), type = "points", 
                             labels = c("p<0.05", "p<0.01", "p<0.001"),
                             legend_gp = gpar(fontsize = 8), 
                             title = "Significance")

# Save the heatmap with the significance legend
pdf("results/humancohort/metab_microb_correlations.pdf", width = 8, height = 12)
draw(heatmap, annotation_legend_list = list(significance_legend))
dev.off()

# Also create separate correlation heatmaps for ALS and Healthy groups
# For ALS patients
als_samples <- meta_subset$ID[meta_subset$Group == "ALS"]
als_met <- met_subset[als_samples, ]
als_mb <- mb_subset[als_samples, ]

# For Healthy controls
healthy_samples <- meta_subset$ID[meta_subset$Group == "Healthy"]
healthy_met <- met_subset[healthy_samples, ]
healthy_mb <- mb_subset[healthy_samples, ]

# Calculate group-specific correlations
als_results <- calculate_corr_matrices(als_met, als_mb, selected_metabolites, selected_microbes)
healthy_results <- calculate_corr_matrices(healthy_met, healthy_mb, selected_metabolites, selected_microbes)

# Save group-specific heatmaps with significance legend
pdf("results/humancohort/als_metab_microb_correlations.pdf", width = 10, height = 8)
draw(create_group_heatmap(als_results$corr, als_results$pval_adj, "ALS"), 
     annotation_legend_list = list(significance_legend))
dev.off()

pdf("results/humancohort/healthy_metab_microb_correlations.pdf", width = 10, height = 8)
draw(create_group_heatmap(healthy_results$corr, healthy_results$pval_adj, "Healthy"), 
     annotation_legend_list = list(significance_legend))
dev.off()

# Export correlation results as tables
write.csv(corr_matrix, "results/humancohort/metab_microb_correlation_values.csv")
write.csv(pval_adj_matrix, "results/humancohort/metab_microb_correlation_pvalues.csv")


# Create scatter plots for significant correlations from the overall correlation analysis
dir.create("results/humancohort/correlation_plots", showWarnings = FALSE, recursive = TRUE)
sig_correlations <- which(pval_matrix < 0.05, arr.ind = TRUE)

# Generate and save individual plots
for (i in 1:nrow(sig_correlations)) {
  met_idx <- sig_correlations[i, 1]
  mb_idx <- sig_correlations[i, 2]
  
  corr_val <- corr_matrix[met_idx, mb_idx]
  p_val <- pval_matrix[met_idx, mb_idx]
  
  p <- create_correlation_plot(met_idx, mb_idx, met_subset, mb_subset, 
                              corr_val, p_val, meta_subset)
  
  # Save individual plot
  metabolite <- rownames(corr_matrix)[met_idx]
  microbe <- colnames(corr_matrix)[mb_idx]
  filename <- paste0("results/humancohort/correlations_metabmicrob/", gsub(" ", "_", metabolite), "_vs_", 
                    gsub(" ", "_", microbe), ".pdf")
  ggsave(filename, p, width = 8, height = 6)
}

# Create a combined figure with multiple plots
# Determine grid dimensions based on number of significant correlations
n_plots <- nrow(sig_correlations)
n_cols <- min(3, n_plots)
n_rows <- ceiling(n_plots / n_cols)

# Create a list of plots
plot_list <- list()
for (i in 1:nrow(sig_correlations)) {
  met_idx <- sig_correlations[i, 1]
  mb_idx <- sig_correlations[i, 2]
  
  corr_val <- corr_matrix[met_idx, mb_idx]
  p_val <- pval_matrix[met_idx, mb_idx]
  
  p <- create_correlation_plot(met_idx, mb_idx, met_subset, mb_subset, 
                              corr_val, p_val, meta_subset)
  plot_list[[i]] <- p
}

# Combine plots
(combined_plot <- ggarrange(plotlist = plot_list, 
                          ncol = n_cols, nrow = n_rows,
                          common.legend = TRUE, legend = "bottom"))
ggsave("results/humancohort/significant_metab_microb_correlations.pdf", 
       width = n_cols * 5, height = n_rows * 4)

# Arrange
heatmap_grob2 <- grid.grabExpr(draw(
      heatmap, 
      annotation_legend_list = list(significance_legend),
      annotation_legend_side = "right",
      padding = unit(c(2, 20, 10, 2), "mm"), # top, right, bottom, left padding
  )
)

heatmap_grob <- grid.grabExpr(draw(
      heatmap, 
      annotation_legend_list = list(lgd_sig),
      annotation_legend_side = "right",
      padding = unit(c(2, 30, 10, 2), "mm"), # top, right, bottom, left padding
  )
)
hm <- as_ggplot(heatmap_grob)
ggsave(hm, filename = "results/humancohort/metab_microb_correlations.pdf", 
        width = 8, height = 12)

(corrs <- ggarrange(plot_list[[3]], plot_list[[4]], plot_list[[5]], nrow = 3, labels = c(LETTERS[12:14])))

block1 <- ggarrange(boxplots, 
                    ggarrange(heatmap_grob, NULL, heights = c(1, 0.6), nrow = 2), 
                    ggarrange(corr_pls, NULL, nrow = 2, heights = c(1, 0.3)),
                    ncol = 3, widths = c(1.8, 1, 1), labels = c("", "H", ""))
(block3 <- ggarrange(as_ggplot(heatmap_grob2), 
                      corrs, 
                      nrow = 1, widths = c(1, 0.60),
                      labels = c(LETTERS[11], "", "")
                      )
                )
ggarrange(boxplots, block2, block3, nrow = 3, heights = c(1.4, 1.0, 1.9))
ggsave("results/humancohort/assembled_microbes_metabolites.pdf", width = 17, height = 28)
