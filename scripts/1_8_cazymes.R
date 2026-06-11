# Cayman analyses human cohort Calgary

# Library
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rio)
library(rstatix)
library(ggrepel)
library(vegan)

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

res_dir <- "results/humancohort2"

# Data
cayman <- readRDS("data/human_cohort2/cayman_families_merged.RDS")
head(cayman)[1:5,1:5]
rownames(cayman) <- cayman$family
cayman$family <- NULL
cayman <- as.data.frame(t(as.matrix(cayman)))
cayman[is.na(cayman)] <- 0
names(cayman)
dim(cayman)
head(cayman)[1:5,1:5]
cayman$sampleID <- rownames(cayman)
fam <- ncol(cayman)
clinical <- readRDS("data/human_cohort2/metadata.RDS")

# --- Beta diversity: CAZyme families ---
# Prepare numeric matrix (all families, unfiltered)
caz_mat <- cayman[, setdiff(names(cayman), "sampleID")]
caz_mat <- caz_mat[rownames(caz_mat) %in% clinical$ID, ]
caz_mat <- log10(caz_mat + 1)

# Beta diversity in ALS
dfals <- clinical |> filter(Group == "ALS")
caz_als <- caz_mat[which(rownames(caz_mat) %in% dfals$ID), ]
bray_als <- vegan::vegdist(caz_als, method = 'bray')
set.seed(14)
pcoord_als <- ape::pcoa(bray_als, correction = "cailliez")
str(pcoord_als$values)
expl_variance_als <- pcoord_als$values$Relative_eig * 100
head(expl_variance_als)
dbray_als <- pcoord_als$vectors[, c('Axis.1', 'Axis.2')]
dbray_als <- as.data.frame(dbray_als)
dbray_als$ID <- rownames(dbray_als)
dbray_als <- left_join(dbray_als, clinical, by = 'ID')
set.seed(14)
res_als <- adonis2(bray_als ~ Sex, data = dbray_als)
res_als <- as.data.frame(res_als)

(pl_caz1 <- ggplot(data = dbray_als, aes(Axis.1, Axis.2)) +
        geom_point(aes(color = Sex), size = 3, alpha = 0.7) +
        xlab(paste0('PCo1 (', round(expl_variance_als[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_als[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_nejm()(2)) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        theme_Publication() +
        labs(color = "", fill = "", title = "CAZymes: Patients with ALS") +
        stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
            alpha = 0.1, linewidth = 1.0) +
        theme(legend.position = "top"))

(pl_caz1 <- pl_caz1 + annotate("text", x = Inf, y = Inf,
                    label = paste0("Sex: p = ", round(res_als$`Pr(>F)`[1], 3), ", R\u00b2 = ", round(res_als$R2[1], 3)),
                    hjust = 1.1, vjust = 1.1, size = 5))

# Beta diversity in controls
dfctrl <- clinical |> filter(Group == "Control")
caz_ctrl <- caz_mat[which(rownames(caz_mat) %in% dfctrl$ID), ]
bray_ctrl <- vegan::vegdist(caz_ctrl, method = 'bray')
pcoord_ctrl <- ape::pcoa(bray_ctrl, correction = "cailliez")
str(pcoord_ctrl$values)
expl_variance_ctrl <- pcoord_ctrl$values$Relative_eig * 100
head(expl_variance_ctrl)
dbray_ctrl <- pcoord_ctrl$vectors[, c('Axis.1', 'Axis.2')]
dbray_ctrl <- as.data.frame(dbray_ctrl)
dbray_ctrl$ID <- rownames(dbray_ctrl)
dbray_ctrl <- left_join(dbray_ctrl, clinical, by = 'ID')
res_ctrl <- adonis2(bray_ctrl ~ Sex, data = dbray_ctrl)
res_ctrl <- as.data.frame(res_ctrl)

(pl_caz2 <- ggplot(data = dbray_ctrl, aes(Axis.1, Axis.2)) +
        geom_point(aes(color = Sex), size = 3, alpha = 0.7) +
        xlab(paste0('PCo1 (', round(expl_variance_ctrl[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_ctrl[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_nejm()(2)) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        theme_Publication() +
        labs(color = "", fill = "", title = "CAZymes: Controls") +
        stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
            alpha = 0.1, linewidth = 1.0) +
        theme(legend.position = "top"))

pl_caz2 <- pl_caz2 + annotate("text", x = Inf, y = Inf,
                    label = paste0("Sex: p = ", round(res_ctrl$`Pr(>F)`[1], 3), ", R\u00b2 = ", round(res_ctrl$R2[1], 3)),
                    hjust = 1.1, vjust = 1.1, size = 5)
pl_caz2
pca_combined <- ggarrange(pl_caz1, pl_caz2, nrow = 1, labels = LETTERS[1:2], common.legend = TRUE, legend = "bottom")
ggsave(file.path(res_dir, "cayman_PCoA_BrayCurtis.pdf"), plot = pca_combined, width = 12, height = 6)

# --- CAZy FAMILY ABUNDANCE DIFFERENCES: ALS vs Control ---
# Prepare: join cayman abundance with clinical metadata
gene_families <- setdiff(names(cayman), c("sampleID"))

# Filter CAZyme families: present in ≥10% of samples AND mean CPM > 1
prevalence <- colMeans(cayman[, gene_families] > 0)
mean_abundance <- colMeans(cayman[, gene_families])
gene_families <- gene_families[prevalence >= 0.30 & mean_abundance > 20]
message(length(gene_families), " CAZyme families retained after filtering (prevalence ≥10%, mean CPM >1)")

df_tot <- cayman |>
  left_join(clinical, by = c("sampleID" = "ID")) |>
  filter(!is.na(Group))

# Keep only families that survived the join as columns
gene_families <- gene_families[gene_families %in% names(df_tot)]

# --- Boxplots for specific CAZy families of interest ---
families_of_interest <- c("GH78", "GH106")
families_present <- families_of_interest[families_of_interest %in% names(df_tot)]

plist <- list()
for (i in seq_along(families_present)) {
  gf <- families_present[i]
  df_tot$cazy_val <- log10(df_tot[[gf]] + 1)

  # LM interaction test: Sex * Group
  model_int <- lm(cazy_val ~ Sex * Group, data = df_tot)
  int_row <- grep(":", rownames(summary(model_int)$coefficients), value = TRUE)
  p_int <- if (length(int_row) > 0) summary(model_int)$coefficients[int_row[1], 4] else NA
  int_text <- paste0("Sex \u00d7 Group interaction: p = ", format.pval(p_int, digits = 2))

  pl <- ggplot(df_tot, aes(x = Sex, y = cazy_val)) +
    geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
    labs(title = gf,
         caption = int_text,
         y = "log10(cpm + 1)",
         x = "") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
          plot.caption = element_text(size = 14),
          legend.position = "none") +
    facet_wrap(~ Group)

  plist[[i]] <- pl
}

(plots_cazy <- ggarrange(plotlist = list(plist[[1]], plist[[2]]), common.legend = TRUE, legend = "bottom",
                         labels = LETTERS[1:2],
                         nrow = 1, ncol = 2))
ggsave(file.path(res_dir, "cayman_families_of_interest_boxplots.pdf"), plots_cazy,
       width = 8, height = 5)
