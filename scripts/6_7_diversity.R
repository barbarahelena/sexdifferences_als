# Diversity plots Human cohort ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ggsci)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
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
                plot.caption = element_text(size = rel(0.7), face = "italic")
        ))
    
} 

## Load data
mb <- readRDS("data/human_cohort/microbiome.RDS")
df <- readRDS("data/human_cohort/metadata.RDS")
head(mb)[1:5,1:5]
rowSums(mb)
mb <- mb[rownames(mb) %in% df$ID, ]

# Alpha diversity
## Shannon plot
shannon <- vegan::diversity(mb, index = 'shannon')
df_shan <- data.frame(ID = names(shannon), shannon = shannon)
df_shan <- left_join(df_shan, df, by = "ID")
p_shan_female <- round(wilcox.test(shannon ~ Group, data = df_shan |> filter(Sex == "Female"))$p.value, 3)
p_shan_male <- round(wilcox.test(shannon ~ Group, data = df_shan |> filter(Sex == "Male"))$p.value, 3)
shan_caption <- paste0("ALS-Control p = ", p_shan_female, " (women); p = ", p_shan_male, " (men)")
(plshan <- ggplot(data = df_shan, aes(x = Sex, y = shannon, fill = Sex)) +
    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
    stat_compare_means(method = "wilcox.test") +
    facet_wrap(~Group) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    labs(title = "Shannon index", y = "Shannon index", x="Sex", caption = shan_caption) +
    theme_Publication())
# ggsave(plshan, filename = "results/humancohort2/shannon.pdf", width = 4, height = 5)

## Species richness
specrich <- specnumber(mb)
dfspec <- data.frame(ID = names(specrich), richness = specrich)
dfspec <- left_join(dfspec, df, by = "ID")
p_rich_female <- round(wilcox.test(richness ~ Group, data = dfspec |> filter(Sex == "Female"))$p.value, 3)
p_rich_male <- round(wilcox.test(richness ~ Group, data = dfspec |> filter(Sex == "Male"))$p.value, 3)
rich_caption <- paste0("ALS-Control p = ", p_rich_female, " (women); p = ", p_rich_male, " (men)")
(plrich <- ggplot(data = dfspec, aes(x = Sex, y = richness, fill = Sex)) +
    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
    stat_compare_means(method = "wilcox.test") +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_wrap(~Group) +
    labs(title = "Richness", y = "Number of species", x="Sex", caption = rich_caption) +
    theme_Publication())
# ggsave(plrich, filename = "results/humancohort2/richness.pdf", width = 4, height = 5)
ggarrange(plshan, plrich, nrow = 1)
ggsave("results/humancohort/alphadiversity.pdf", width = 10, height = 6)

# Beta diversity in ALS
dfals <- df |> filter(Group == "ALS") |> filter(!is.na(Sex))
mb2 <- mb[which(rownames(mb) %in% dfals$ID),]
bray <- vegan::vegdist(mb2, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
head(expl_variance_bray)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$ID <- rownames(dbray)
dbray <- left_join(dbray, df, by = 'ID') # add metadata / covariates
res <- adonis2(bray ~ Sex, data = dbray, by = "terms")
res <- as.data.frame(res)

## Plots
(pl1 <- ggplot(data = dbray, aes(Axis.1, Axis.2)) +
        geom_point(aes(color = Sex), size = 3, alpha = 0.7) +
        xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_nejm()(2)) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        theme_Publication() +
        labs(color = "", fill = "", title = "Patients with ALS") +
        stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
            alpha = 0.1, linewidth = 1.0) +
        theme(legend.position = "top"))


pl1 <- pl1 + annotate("text", x = Inf, y = Inf, label = paste0("Sex: p = ", round(res$`Pr(>F)`[1], 3)),
                    hjust = 1.1, vjust = 1.1, size = 3)
pl1

# Beta diversity in controls
dfctrl <- df |> filter(Group == "Healthy")
mb2 <- mb[which(rownames(mb) %in% dfctrl$ID),]
bray <- vegan::vegdist(mb2, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Relative_eig * 100
head(expl_variance_bray)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$ID <- rownames(dbray)
dbray <- left_join(dbray, df, by = 'ID') # add metadata / covariates
res <- adonis2(bray ~ Sex, data = dbray)
res <- as.data.frame(res)

## Plots
(pl2 <- ggplot(data = dbray, aes(Axis.1, Axis.2)) +
        geom_point(aes(color = Sex), size = 3, alpha = 0.7) +
        xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_nejm()(2)) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        theme_Publication() +
        labs(color = "", fill = "", title = "Controls") +
        stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
            alpha = 0.1, linewidth = 1.0) +
        theme(legend.position = "top"))


pl2 <- pl2 + annotate("text", x = Inf, y = Inf, label = paste0("Sex: p = ", round(res$`Pr(>F)`[1], 3)),
                    hjust = 1.1, vjust = 1.1, size = 3)
pl2

ggarrange(pl1, pl2, nrow = 1, labels = LETTERS[1:2], common.legend = TRUE, legend = "bottom")
ggsave("results/humancohort/PCoA_BrayCurtis.pdf", width = 12, height = 6)
