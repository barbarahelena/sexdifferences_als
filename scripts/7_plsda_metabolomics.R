## Volcano plot metabolomics sex differences
## Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

# Libraries
library(mixOmics)
library(ggplot2)

# Prepare data
X <- df_sum[, -which(names(df_sum) == "group")]  # Features
Y <- df_sum$group  # Class labels

# Perform PLS-DA
plsda_result <- plsda(X, Y, ncomp = 2)

# Plot PLS-DA results
plotIndiv(plsda_result, comp = 1:2, group = Y, 
          legend = TRUE, title = "PLS-DA on df_sum")

# Save the plot
ggsave("r_results/plsda_plot.pdf", width = 8, height = 6, device = "pdf")

# Optional: Evaluate the model
perf_plsda <- perf(plsda_result, validation = "Mfold", folds = 5, progressBar = TRUE, nrepeat = 10)
print(perf_plsda)

# Save the performance results
saveRDS(perf_plsda, file = "r_results/plsda_performance.rds")