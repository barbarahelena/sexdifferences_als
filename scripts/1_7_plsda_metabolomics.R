## Volcano plot metabolomics sex differences
## Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

# Libraries
library(mixOmics)
library(ggplot2)
library(tidyverse)
library(ggsci)

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
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
}

# Prepare data
met <- readRDS("data/metabolomics.RDS")
meta <- readRDS("data/metadata.RDS")
groups <- meta$GroupPerSex  # group labels

# Perform PLS-DA
set.seed(123)
plsda_model <- plsda(met, groups, ncomp = 2, scale = TRUE)

perf_model <- perf(plsda_model, 
                   validation = "Mfold",  # Cross-validation method
                   folds = 5,             # Number of folds
                   nrepeat = 100,          # Number of repeats (adjust this)
                   progressBar = TRUE)

# Model performance
perf_model$error.rate  # Classification error rates per component
perf_model$Q2.total    # Q² values (predictive ability)
perf_model$R2          # R² values (variance explained)

# Extract PLS-DA scores
scores_df <- data.frame(plsda_model$variates$X, Group = as.factor(groups))
pc_ev1 <- round(plsda_model$explained_variance$X[[1]] * 100, 2)
pc_ev2 <- round(plsda_model$explained_variance$X[[2]] * 100, 2)

ggplot(scores_df, aes(x = comp.1, y = comp.2, color = Group)) +
    stat_ellipse(geom = "polygon", aes(fill = Group, color = Group),
                        alpha = 0.3) +
    scale_color_nejm() +
    scale_fill_nejm() +
    geom_point(size = 3, alpha = 0.8) +  # Customize point size & transparency
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    labs(title = "PLS-DA", x = str_c('Component 1 (', pc_ev1, '%)'), 
        y = str_c('Component 2 (', pc_ev2, '%)'),
        color = "", fill = "") +
    theme_Publication()

(pca <- ggplot(df, aes(x=PC1, y=PC2, color=Group)) +
        stat_ellipse(geom = "polygon", aes(fill = Group, color = Group),
                     alpha = 0.3) +
        geom_point() +
        ggtitle('PCA plasma metabolites') +
        scale_color_bmj() +
        scale_fill_bmj() +
        theme_minimal() +
        labs(x=str_c('PC1 ', pc1_ev, '%'), y=str_c('PC2 ', pc2_ev, '%')) +
        theme_Publication() +
        theme(legend.title = element_blank()))

# Save the plot
ggsave("results/plsda_plot.pdf", width = 8, height = 6, device = "pdf")

# Optional: Evaluate the model
perf_plsda <- perf(plsda_result, validation = "Mfold", folds = 5, progressBar = TRUE, nrepeat = 10)
print(perf_plsda)

# Save the performance results
saveRDS(perf_plsda, file = "results/plsda_performance.rds")