## Pathways - microbes in human hosts
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

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


## Data
mb <- readRDS("data/human_cohort2/microbiome.RDS") %>% mutate(across(where(is.numeric), log10))
microb <- read.csv("results/humancohort2/boxplots_validation/microbes_validate.csv")
mb_cols <- colnames(mb)
mb_subset <- mb[, microb$x, drop = FALSE]
bugs <- names(mb_subset)
mb_subset$ID <- rownames(mb_subset)

df <- readRDS("data/human_cohort2/human_als_pathways.RDS") %>% filter(rownames(.) %in% rownames(mb_subset)) # only samples that are in the microbiome data
head(df)
dim(df)
pathways <- colnames(df)[1:ncol(df)-1]
length(pathways)
pathways[which(str_detect(pathways, "RHAM"))]
pathways2 <- c("DTDPRHAMSYN-PWY", "RHAMCAT-PWY")
keypath <- readRDS("data/human_cohort2/pathwaykeys.RDS")
meta <- readRDS("data/human_cohort2/metadata.RDS")

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
    print(bug)
    for (pathway in pathways) {
      print(pathway)
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
ggsave("results/humancohort2/heatmap_microbes_pathways_complex.pdf", plot = as_ggplot(heatmap_grob),
       width = 12, height = 10)

# Save
# Replace your PDF save code with this:
# cairo_pdf("results/humancohort2/heatmap_microbes_pathways.pdf", width = 15, height = 45,
#           family = "Helvetica")
# draw(heatmap_anno, padding = unit(c(5,5,5,20), "mm"))
# dev.off()

### Rhamnose pathways
# Calculate correlations and identify significant ones
correlations <- correlate_bugs_pathways(merged_df, bugs, pathways2) %>% #pathwys2 contains rhamnose pathways
    mutate(q.value = p.adjust(P.value, method = "fdr")) %>%
    left_join(keypath, by = c("Pathway" = "keys"))
sig_pathways <- correlations %>% # filter(P.value < 0.05) %>% 
    pull(Pathway) %>% unique()
correlations_filtered <- correlations %>% filter(Pathway %in% sig_pathways)

# Transform correlation data to matrix format - TRANSPOSED
corr_matrix <- correlations %>%
  select(Bug, Pathway, Correlation) %>%
  pivot_wider(names_from = Bug, values_from = Correlation) %>%
  column_to_rownames("Pathway")

# Store q-values in a matrix format for cell annotation - TRANSPOSED
qval_matrix <- correlations %>%
  select(Bug, Pathway, q.value) %>%
  pivot_wider(names_from = Bug, values_from = q.value) %>%
  column_to_rownames("Pathway")

# Get pathway explanations for better labels
pathway_labels <- correlations %>%
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
# heatmap_anno <- draw(heatmap_bugs_pathways, annotation_legend_list = list(lgd_sig)) # no significant pathways
heatmap_anno <- draw(heatmap_bugs_pathways)

# Save
cairo_pdf("results/humancohort2/heatmap_microbes_pathways_rhamnose.pdf", width = 8, height = 4,
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
(paths_sexdiff <- ggarrange(plotlist = res_box, ncol = 2, nrow = 1))
ggsave("results/humancohort2/rhamnose_pathways_sexdifferences_pval.pdf", width = 9, height = 4)

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

## Broader approach: testing differences ALL pathways
df$ID <- rownames(df)
df_tot <- left_join(meta, df, by = c("ID"))
names(df_tot)

statres <- c()
for(i in pathways) {
    df_tot$pathway_val <- df_tot[,i]
    pwname <- i
    test <- wilcox.test(pathway_val ~ Group, data = df_tot)
    pval_diag <- as.numeric(test$p.value)
    pval_diag <- as.numeric(format(round(pval_diag, 3), nsmall = 3))
    sig_diag <- case_when(
        pval_diag < 0.0001 ~ paste0("****"),
        pval_diag < 0.001 ~paste0("***"),
        pval_diag < 0.01 ~paste0("**"),
        pval_diag <= 0.05 ~paste0("*"),
        pval_diag > 0.05 ~paste0("")
    )
    statres_line <- cbind(pwname, pval_diag, sig_diag)
    statres <- rbind(statres, statres_line)
}
statres <- as.data.frame(statres)
statres <- statres %>% arrange(pval_diag) %>%
    mutate(qval_diag = p.adjust(pval_diag, method = "fdr")) |>
    mutate(sigq_diag = case_when(
        qval_diag < 0.0001 ~ paste0("****"),
        qval_diag < 0.001 ~paste0("***"),
        qval_diag < 0.01 ~paste0("**"),
        qval_diag <= 0.05 ~paste0("*"),
        qval_diag > 0.05 ~paste0("")
    ))

statresfem <- c()
for(i in pathways) {
    df_tot2 <- df_tot |> filter(Sex == "Female")
    df_tot2$pathway_val <- df_tot2[,i]
    pwname <- i
    test <- wilcox.test(pathway_val ~ Group, data = df_tot2)
    pval_fem <- as.numeric(test$p.value)
    pval_fem <- as.numeric(format(round(pval_fem, 3), nsmall = 3))
    sig_fem <- case_when(
        pval_fem < 0.0001 ~ paste0("****"),
        pval_fem < 0.001 ~paste0("***"),
        pval_fem < 0.01 ~paste0("**"),
        pval_fem <= 0.05 ~paste0("*"),
        pval_fem > 0.05 ~paste0("")
    )
    statres_line <- cbind(pwname, pval_fem, sig_fem)
    statresfem <- rbind(statresfem, statres_line)
}
statresfem <- as.data.frame(statresfem)
statresfem <- statresfem %>% arrange(pval_fem) %>%
    mutate(qval_fem = p.adjust(pval_fem, method = "fdr")) |>
    mutate(sigq_fem = case_when(
        qval_fem < 0.0001 ~ paste0("****"),
        qval_fem < 0.001 ~paste0("***"),
        qval_fem < 0.01 ~paste0("**"),
        qval_fem <= 0.05 ~paste0("*"),
        qval_fem > 0.05 ~paste0("")
    ))

statresmale <- c()
for(i in pathways) {
    df_tot2 <- df_tot |> filter(Sex == "Male")
    df_tot2$pathway_val <- df_tot2[,i]
    pwname <- i
    test <- wilcox.test(pathway_val ~ Group, data = df_tot2)
    pval_male <- as.numeric(test$p.value)
    pval_male <- as.numeric(format(round(pval_male, 3), nsmall = 3))
    sig_male <- case_when(
        pval_male < 0.0001 ~ paste0("****"),
        pval_male < 0.001 ~paste0("***"),
        pval_male < 0.01 ~paste0("**"),
        pval_male <= 0.05 ~paste0("*"),
        pval_male > 0.05 ~paste0("")
    )
    statres_line <- cbind(pwname, pval_male, sig_male)
    statresmale <- rbind(statresmale, statres_line)
}
statresmale <- as.data.frame(statresmale)
statresmale <- statresmale %>% arrange(pval_male) %>%
    mutate(qval_male = p.adjust(pval_male, method = "fdr")) |>
    mutate(sigq_male = case_when(
        qval_male < 0.0001 ~ paste0("****"),
        qval_male < 0.001 ~paste0("***"),
        qval_male < 0.01 ~paste0("**"),
        qval_male <= 0.05 ~paste0("*"),
        qval_male > 0.05 ~paste0("")
    ))

statresals <- c()
for(i in pathways) {
    df_tot3 <- df_tot |> filter(Group == "ALS")
    df_tot3$pathway_val <- df_tot3[,i]
    pwname <- i
    test <- wilcox.test(pathway_val ~ Sex, data = df_tot3)
    pval_als <- as.numeric(test$p.value)
    pval_als <- as.numeric(format(round(pval_als, 3), nsmall = 3))
    sig_als <- case_when(
        pval_als < 0.0001 ~ paste0("****"),
        pval_als < 0.001 ~paste0("***"),
        pval_als < 0.01 ~paste0("**"),
        pval_als <= 0.05 ~paste0("*"),
        pval_als > 0.05 ~paste0("")
    )
    statres_line <- cbind(pwname, pval_als, sig_als)
    statresals <- rbind(statresals, statres_line)
}
statresals <- as.data.frame(statresals)
statresals <- statresals %>% arrange(pval_als) %>%
    mutate(qval_als = p.adjust(pval_als, method = "fdr")) |>
    mutate(sigq_als = case_when(
        qval_als < 0.0001 ~ paste0("****"),
        qval_als < 0.001 ~paste0("***"),
        qval_als < 0.01 ~paste0("**"),
        qval_als <= 0.05 ~paste0("*"),
        qval_als > 0.05 ~paste0("")
    ))

statresctrl <- c()
for(i in pathways) {
    df_tot3 <- df_tot |> filter(Group == "Control")
    df_tot3$pathway_val <- df_tot3[,i]
    pwname <- i
    test <- wilcox.test(pathway_val ~ Sex, data = df_tot3)
    pval_ctrl <- as.numeric(test$p.value)
    pval_ctrl <- as.numeric(format(round(pval_ctrl, 3), nsmall = 3))
    sig_ctrl <- case_when(
        pval_ctrl < 0.0001 ~ paste0("****"),
        pval_ctrl < 0.001 ~paste0("***"),
        pval_ctrl < 0.01 ~paste0("**"),
        pval_ctrl <= 0.05 ~paste0("*"),
        pval_ctrl > 0.05 ~paste0("")
    )
    statres_line <- cbind(pwname, pval_ctrl, sig_ctrl)
    statresctrl <- rbind(statresctrl, statres_line)
}
statresctrl <- as.data.frame(statresctrl)
statresctrl <- statresctrl %>% arrange(pval_ctrl) %>%
    mutate(qval_ctrl = p.adjust(pval_ctrl, method = "fdr")) |>
    mutate(sigq_ctrl = case_when(
        qval_ctrl < 0.0001 ~ paste0("****"),
        qval_ctrl < 0.001 ~paste0("***"),
        qval_ctrl < 0.01 ~paste0("**"),
        qval_ctrl <= 0.05 ~paste0("*"),
        qval_ctrl > 0.05 ~paste0("")
    ))

restot <- right_join(statres, statresfem) |> right_join(statresmale) |> right_join(statresals) |> right_join(statresctrl)
head(restot)
restot |> arrange(qval_diag)
restot |> arrange(qval_fem)
restot |> arrange(qval_male)
restot |> arrange(qval_als)
restot |> arrange(qval_ctrl)

res_sigdiag <- restot |> filter(pval_als < 0.05)
dim(res_sigdiag)

plist <- list()
for(i in 1:nrow(res_sigdiag)){
    nm <- res_sigdiag$pwname[i]
    pw_label <- keypath$expl[keypath$keys == nm]
    if(length(pw_label) == 0) pw_label <- nm
    df_tot$pathway_val <- df_tot[[nm]]
    df_fem <- df_tot |> filter(Sex == "Female")
    df_fem$pathway_val <- df_fem[[nm]]
    df_mal <- df_tot |> filter(Sex == "Male")
    df_mal$pathway_val <- df_mal[[nm]]

    diagdiff1 <- wilcox.test(pathway_val ~ Group, data = df_fem)
    diagdiff2 <- wilcox.test(pathway_val ~ Group, data = df_mal)

    # Format p-values
    p_female <- format.pval(diagdiff1$p.value, digits = 2)
    p_male <- format.pval(diagdiff2$p.value, digits = 2)

    # Create subtitle with sex difference p-values
    diag_diff_text <- paste0("ALS-Control p = ", p_female, " (women); p = ", p_male, " (men)")

    pl1 <- ggplot(df_tot, aes(x = Sex, y = pathway_val)) +
        geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.4) +
        geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        labs(title = str_wrap(pw_label, width = 40),
             caption = diag_diff_text,
             y = "log10(cpm)",
             x = "") +
        stat_compare_means(method = "wilcox.test", label = "p.format") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
              plot.caption = element_text(size = 12),
              legend.position = "none") +
        facet_wrap(~ Group)

    plist[[i]] <- pl1
}

n_plots_diag <- nrow(res_sigdiag)
n_cols_diag <- min(4, n_plots_diag)
n_rows_diag <- ceiling(n_plots_diag / n_cols_diag)

(plots <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:n_plots_diag],
          nrow = n_rows_diag, ncol = n_cols_diag))
ggsave(plots, filename = "results/humancohort2/diffpathways_boxplots.pdf",
       width = 4 * n_cols_diag, height = 5 * n_rows_diag)
write.csv2(restot, "results/humancohort2/wilcoxons_pathways.csv")

# Boxplots for pathways with pval_als < 0.05 (sex differences in ALS)
res_sigals <- statresals |> filter(pval_als < 0.05)
dim(res_sigals)

plist_als <- list()
for(i in 1:nrow(res_sigals)){
    nm <- res_sigals$pwname[i]
    pw_label <- keypath$expl[keypath$keys == nm]
    if(length(pw_label) == 0) pw_label <- nm
    df_tot$pathway_val <- df_tot[[nm]]
    df_fem <- df_tot |> filter(Sex == "Female")
    df_fem$pathway_val <- df_fem[[nm]]
    df_mal <- df_tot |> filter(Sex == "Male")
    df_mal$pathway_val <- df_mal[[nm]]

    # Test diagnosis difference in females and males
    diagdiff1 <- wilcox.test(pathway_val ~ Group, data = df_fem)
    diagdiff2 <- wilcox.test(pathway_val ~ Group, data = df_mal)

    # Format p-values
    p_female <- format.pval(diagdiff1$p.value, digits = 2)
    p_male <- format.pval(diagdiff2$p.value, digits = 2)

    # Create subtitle with diagnosis difference p-values
    diag_diff_text <- paste0("ALS-Control p = ", p_female, " (women); p = ", p_male, " (men)")

    pl_als <- ggplot(df_tot, aes(x = Sex, y = pathway_val)) +
        geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.4) +
        geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        labs(title = str_wrap(pw_label, width = 40),
             caption = diag_diff_text,
             y = "log10(cpm)",
             x = "") +
        stat_compare_means(method = "wilcox.test", label = "p.format") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
              plot.caption = element_text(size = 12),
              legend.position = "none") +
        facet_wrap(~ Group)

    plist_als[[i]] <- pl_als
}

n_plots <- length(plist_als)
n_cols <- ceiling(sqrt(n_plots))
n_rows <- ceiling(n_plots / n_cols)

(plots_als <- ggarrange(plotlist = plist_als, common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:n_plots],
          nrow = n_rows, ncol = n_cols))
ggsave(plots_als, filename = "results/humancohort2/sexdiff_als_pathways_boxplots.pdf",
       width = 4 * n_cols, height = 5 * n_rows)
