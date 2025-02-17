# Beta diversity plots ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(aplot)
library(doParallel)
registerDoParallel(8)

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
        legend.spacing = unit(0, "cm"),
        # legend.title = element_text(face="italic"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"),
        plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
}

## Load data
mb <- readRDS("data/microbiome.RDS")
df <- readRDS("data/meta_microbiome.RDS")
head(mb)[1:5,1:5]
rowSums(mb)

#### Bray-Curtis distance ####
bray <- vegan::vegdist(mb, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
head(expl_variance_bray)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$ID <- rownames(dbray)
dbray <- left_join(dbray, df, by = 'ID') # add metadata / covariates

## PERMANOVAs
dfanova <- df %>% slice(match(rownames(mb), ID)) # distance matrix and metadata must have the same sample order
all(dfanova$ID == rownames(mb)) # TRUE
dim(dfanova)
set.seed(14)
res1 <- adonis2(bray ~ Age_weeks, data = dfanova) # PERMANOVA
print(res1)
set.seed(14)
res2 <- adonis2(bray ~ Genotype, data = dfanova) # PERMANOVA
print(res2)
set.seed(14)
res3 <- adonis2(bray ~ Sex, data = dfanova) # PERMANOVA
print(res3)

#### Bray-Curtis per time point ####
braypertimepoint_mice <- function(timepoint, df = dfanova, tab = mb) {
    set.seed(14)
    dfsel <- df %>% filter(Age_weeks == timepoint)
    bray <- vegan::vegdist(tab[rownames(tab) %in% dfsel$ID,], method = 'bray')
    print(adonis2(bray ~ Genotype, data = dfsel))
}

braypertimepoint_mice(timepoint = "8 weeks")
braypertimepoint_mice(timepoint = "10 weeks")
braypertimepoint_mice(timepoint = "12 weeks")
braypertimepoint_mice(timepoint = "14 weeks")
braypertimepoint_mice(timepoint = "16 weeks")
braypertimepoint_mice(timepoint = "18 weeks")
braypertimepoint_mice(timepoint = "20 weeks")
braypertimepoint_mice(timepoint = "26 weeks")

braypertimepoint_sex <- function(timepoint, mouse, df = dfanova, tab = mb) {
    set.seed(14)
    dfsel <- df %>% filter(Age_weeks == timepoint & Genotype == mouse)
    bray <- vegan::vegdist(tab[rownames(tab) %in% dfsel$ID,], method = 'bray')
    print(adonis2(bray ~ Sex, data = dfsel))
}

braypertimepoint_sex(timepoint = "8 weeks", mouse = "TDP43")
braypertimepoint_sex(timepoint = "10 weeks", mouse = "TDP43")
braypertimepoint_sex(timepoint = "12 weeks", mouse = "TDP43")
braypertimepoint_sex(timepoint = "14 weeks", mouse = "TDP43")
braypertimepoint_sex(timepoint = "16 weeks", mouse = "TDP43")
braypertimepoint_sex(timepoint = "18 weeks", mouse = "TDP43")

braypertimepoint_sex(timepoint = "8 weeks", mouse = "WT")
braypertimepoint_sex(timepoint = "10 weeks", mouse = "WT")
braypertimepoint_sex(timepoint = "12 weeks", mouse = "WT")
braypertimepoint_sex(timepoint = "14 weeks", mouse = "WT")
braypertimepoint_sex(timepoint = "16 weeks", mouse = "WT")
braypertimepoint_sex(timepoint = "18 weeks", mouse = "WT")

## Plots
(braycurt <- dbray %>%
                ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(aes(color = Genotype, shape = Sex), size = 2, alpha = 0.7) +
                    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
                    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
                    scale_color_manual(values = pal_nejm()(4)[c(3:4)]) +
                    scale_fill_manual(values = pal_nejm()(4)[c(3:4)], guide = "none") +
                    theme_Publication() +
                    labs(color = "", fill = "", title = "PCoA Bray-Curtis distance") +
                    stat_ellipse(geom = "polygon", aes(color = Genotype, fill = Genotype), type = "norm",
                        alpha = 0.1, linewidth = 1.0) +
                    theme(legend.position = "top") +
                    facet_wrap(~ Age_weeks))
ggsave(braycurt, filename = "results/microbiome/betadiversity/PCoA_BrayCurtis_genotype.pdf", width = 10, height = 10)


# Sex differences in TDP43
dfsel <- df %>% filter(Genotype == "TDP43" & Age_ints %in% c(6,10,12,14,16))
bray <- vegan::vegdist(mb[rownames(mb) %in% dfsel$ID,], method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
head(expl_variance_bray)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$ID <- rownames(dbray)
dbray <- left_join(dbray, dfsel, by = 'ID') # add metadata / covariates

(braycurt <- dbray %>%
                ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
                    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
                    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
                    scale_color_manual(values = pal_nejm()(4)[3:4]) +
                    scale_fill_manual(values = pal_nejm()(4)[3:4], guide = "none") +
                    theme_Publication() +
                    labs(color = "", fill = "", title = "PCoA Bray-Curtis distance TDP43") +
                    stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
                    alpha = 0.1, linewidth = 1.0) +
                    theme(legend.position = "top") +
                    facet_wrap(~ Age_weeks, nrow = 1))
ggsave(braycurt, filename = "results/microbiome/betadiversity/PCoA_BrayCurtis_tdp43_sex.pdf", 
    width = 14, height = 5)

# Sex differences in WT
dfsel <- df %>% filter(Genotype == "WT" & Age_ints %in% c(6,10,12,14,16))
bray <- vegan::vegdist(mb[rownames(mb) %in% dfsel$ID,], method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
head(expl_variance_bray)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$ID <- rownames(dbray)
dbray <- left_join(dbray, dfsel, by = 'ID') # add metadata / covariates
(braycurt <- dbray %>%
                ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
                    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
                    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
                    scale_color_manual(values = pal_nejm()(2)) +
                    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
                    theme_Publication() +
                    labs(color = "", fill = "", title = "PCoA Bray-Curtis distance WT") +
                    stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
                    alpha = 0.1, linewidth = 1.0) +
                    theme(legend.position = "top") +
                    facet_wrap(~ Age_weeks, nrow = 1))
ggsave(braycurt, filename = "results/microbiome/betadiversity/PCoA_BrayCurtis_ctrl_sex.pdf", 
    width = 14, height = 5)
