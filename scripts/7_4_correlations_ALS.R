# Correlations with ALS scores: corrplots and heatmap
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Library
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)

# Functions
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

calculate_correlations <- function(data, microb) {
  corr_results <- numeric(length(microb))
  pval_results <- numeric(length(microb))
  for (i in seq_along(microb)) {
    bact <- microb[i]    
    min_value <- min(tot[[bact]][which(tot[[bact]] > 0)])
    x <- log10(data[[bact]] + min_value)
    y <- as.numeric(data$ALS)

    cor_test <- cor.test(x, y, method = "spearman")
    corr_results[i] <- cor_test$estimate
    pval_results[i] <- cor_test$p.value
  }
  return(list(correlation = corr_results, pvalue = pval_results))
}


# Data
mb <- readRDS("data/human_cohort2/microbiome.RDS")
df <- readRDS("data/human_cohort2/metadata.RDS")
mb <- mb[rownames(mb) %in% df$ID, ]
dim(mb)

# Select bacteria
names(mb)
names_of_interest <- c("Akkermansia muciniphila", "bacterium 1xD42-87", "Phocaeicola sartorii", 
                        "Parabacteroides distasonis", "Porphyromonadaceae bacterium UBA7141", 
                        "Lachnospiraceae bacterium UBA7143", "Erysipelotrichaceae bacterium",
                        "Bacteroidales bacterium M5", "Parabacteroides goldsteinii", "uncultured Prevotella sp.",
                        "Clostridiales bacterum VE202-03", "Bacteroides acidifaciens", "Flavonifractor plautii",
                        "Bifidobacterium pseudolongum", "Lactobacillus intestinalis")

mb <- mb[, names(mb) %in% names_of_interest]
dim(mb)
rownames(mb)
mb$ID <- rownames(mb)
tot <- left_join(mb, df, by = "ID")
head(tot)
microb <- names(tot[1:14,])[which(colSums(tot[1:7]) / nrow(tot) > 0.01)]
tot <- tot %>% select(all_of(microb), ID:ALSFRS)
head(tot)
write.csv(microb, "results/humancohort2/microbes_validate.csv", row.names = FALSE)

# Correlation plots with ALS score
## 1: Flavonifractor plautii
ggplot(tot %>% filter(Group == "ALS"), 
    aes(x = log10(`Flavonifractor plautii` + min(tot$`Flavonifractor plautii`[which(tot$`Flavonifractor plautii` > 0)])), 
        y = as.numeric(ALSFRS), color = Sex)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_cor(method = "spearman") +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Flavonifractor plautii", x = "log10(Flavonifractor plautii)",
    y = "ALS Score", color = "Sex") +
  theme_Publication() +
  scale_color_nejm()
ggsave("results/humancohort2/als_flavonifractor.pdf", width = 7, height = 7)

## 4: Parabacteroides distasonis
ggplot(tot %>% filter(Group == "ALS"), 
    aes(x = log10(`Parabacteroides distasonis` + min(tot$`Parabacteroides distasonis`[which(tot$`Parabacteroides distasonis` > 0)])), 
      y = as.numeric(ALSFRS), color = Sex)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_cor(method = "spearman") +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Parabacteroides distasonis", x = "log10(Parabacteroides distasonis)",
    y = "ALS-FRS Score", color = "Sex") +
  theme_Publication() +
  scale_color_nejm()
ggsave("results/humancohort2/als_parabacteroides.pdf", width = 7, height = 7)

## 5: Parabacteroides goldsteinii
ggplot(tot %>% filter(Group == "ALS"), 
    aes(x = log10(`Parabacteroides goldsteinii` + min(tot$`Parabacteroides goldsteinii`[which(tot$`Parabacteroides goldsteinii` > 0)])), 
      y = as.numeric(ALSFRS), color = Sex)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_cor(method = "spearman") +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Parabacteroides goldsteinii", x = "log10(Parabacteroides goldsteinii)",
    y = "ALS-FRS Score", color = "Sex") +
  theme_Publication() +
  scale_color_nejm()
ggsave("results/humancohort2/als_parabacteroidesgold.pdf", width = 7, height = 7)

## 7: Akkermansia muciniphila
(pl9 <- ggplot(tot %>% filter(Group == "ALS"), 
    aes(x = log10(`Akkermansia muciniphila` + min(tot$`Akkermansia muciniphila`[which(tot$`Akkermansia muciniphila` > 0)])), 
          y = as.numeric(ALSFRS), color = Sex)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_cor(method = "spearman") +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Akkermansia muciniphila", x = "log10(Akkermansia muciniphila)",
    y = "ALS-FRS Score", color = "Sex") +
  theme_Publication() +
  scale_color_nejm())
ggsave("results/humancohort2/als_akkermansia.pdf", width = 7, height = 7)

# Define als_data for heatmap correlation analysis
als_data <- tot %>% filter(Group == "ALS")

# Heatmap with correlations with ALS score
## Correlation matrix
male_data <- als_data %>% filter(Sex == "Male")
male_corr <- calculate_correlations(male_data, microb)
female_data <- als_data %>% filter(Sex == "Female")
female_corr <- calculate_correlations(female_data, microb)
corr_matrix <- cbind(Male = male_corr$correlation, Female = female_corr$correlation)
rownames(corr_matrix) <- microb

## Create p-value matrix
pval_matrix <- cbind(Male = male_corr$pvalue, Female = female_corr$pvalue)
rownames(pval_matrix) <- microb

col_fun <- colorRamp2(c(-1, 0, 1), c(pal_nejm()(3)[3], "white", pal_nejm()(6)[6]))

## Create the heatmap
heatmap <- Heatmap(
  t(corr_matrix),  # Transpose the correlation matrix
  name = "Correlation",
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    # After transposition, indices are flipped
    corr_val <- t(corr_matrix)[i, j]
    p_val <- t(pval_matrix)[i, j]
    
    # Add correlation value
    grid.text(sprintf("%.2f", corr_val), x, y, gp = gpar(fontsize = 10))
    
    # Add asterisks for significance
    if (p_val < 0.05) {
      stars <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", "*"))
      grid.text(stars, x, y + unit(0.3, "cm"), gp = gpar(fontsize = 12, fontface = "bold"))
    }
  },
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title = "Sex",  
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "Correlation of Microbes with ALS-FRS",  
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Correlation",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1")
  )
)

## Create a legend for significance
significance_legend <- Legend(pch = c("*","**"), type = "points", 
                    labels = c("p<0.05", "p<0.01"),
                    legend_gp = gpar(fontsize = 8), title = "Significance",)

lgd_sig <- significance_legend

## Save the heatmap
# pdf("results/humancohort/als_microbe_correlation_heatmap.pdf", width = 6, height = 6)
# draw(
#   heatmap,
#   annotation_legend_list = list(significance_legend),
#   annotation_legend_side = "right"
# )
# dev.off()

heatmap_grob <- grid.grabExpr(draw(
      heatmap, 
      annotation_legend_list = list(lgd_sig),
      annotation_legend_side = "right",
      padding = unit(c(2, 20, 10, 2), "mm"), # top, right, bottom, left padding
  )
)
hm <- as_ggplot(heatmap_grob)
ggsave(hm, filename = "results/humancohort2/als_microbe_correlation_heatmap.pdf", 
        width = 6, height = 3)

corr_pls <- ggarrange(pl9, nrow = 1, labels = c(LETTERS[9]))

block2 <- ggarrange(ggarrange(as_ggplot(heatmap_grob), NULL, nrow =2, heights = c(1, 0.2)),
          ggarrange(pl9), widths = c(2, 1.2), ncol = 2,
          labels = c(LETTERS[8], LETTERS[9]))
ggsave("results/humancohort2/als_microbes_correlations.pdf", width = 9, height = 3.5)
