# Beta diversity plots ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
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

#### Bray-Curtis per time point ####
braypertimepoint_mice <- function(timepoint, df = dbray, tab = mb) {
    set.seed(14)
    dfsel <- df %>% filter(Age_weeks == timepoint)
    bray <- vegan::vegdist(tab[rownames(tab) %in% dfsel$ID,], method = 'bray')
    return(adonis2(bray ~ Genotype, data = dfsel))
}

wk6 <- c(braypertimepoint_mice(timepoint = "6 weeks")['Pr(>F)'][[1]][1], "6 weeks")
wk8 <- c(braypertimepoint_mice(timepoint = "8 weeks")['Pr(>F)'][[1]][1], "8 weeks")
wk10 <- c(braypertimepoint_mice(timepoint = "10 weeks")['Pr(>F)'][[1]][1], "10 weeks")
wk12 <- c(braypertimepoint_mice(timepoint = "12 weeks")['Pr(>F)'][[1]][1], "12 weeks")
wk14 <- c(braypertimepoint_mice(timepoint = "14 weeks")['Pr(>F)'][[1]][1], "14 weeks")
wk16 <- c(braypertimepoint_mice(timepoint = "16 weeks")['Pr(>F)'][[1]][1], "16 weeks")
wk18 <- c(braypertimepoint_mice(timepoint = "18 weeks")['Pr(>F)'][[1]][1], "18 weeks")
res <- rbind(wk6, wk8, wk10, wk12, wk14, wk16, wk18)
colnames(res) <- c("pvalue", "Age_weeks")
res <- as.data.frame(res)
res$Age_weeks <- factor(res$Age_weeks, levels = c("6 weeks", "8 weeks", "10 weeks", "12 weeks", "14 weeks", "16 weeks", "18 weeks"))
res$pvalue <- as.numeric(res$pvalue)

## Plots
(braycurt <- dbray %>% filter(!Age_ints %in% c(20,22,26)) %>% 
                ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(aes(color = Genotype, shape = Sex), size = 3, alpha = 0.7) +
                    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
                    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
                    scale_color_manual(values = pal_nejm()(6)[c(3,6)]) +
                    scale_fill_manual(values = pal_nejm()(6)[c(3,6)], guide = "none") +
                    theme_Publication() +
                    labs(color = "", fill = "", title = "PCoA Bray-Curtis distance") +
                    stat_ellipse(geom = "polygon", aes(color = Genotype, fill = Genotype), type = "norm",
                        alpha = 0.1, linewidth = 1.0) +
                    theme(legend.position = "top") +
                    geom_text(data = res, aes(x = Inf, y = Inf, label = paste0("p = ", round(pvalue, 3))),
                              hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
                    facet_wrap(~ Age_weeks))
ggsave(braycurt, filename = "results/microbiome/betadiversity/PCoA_BrayCurtis_genotype.pdf", width = 8, height = 8)


### Bray-Curtis sex differences ###
braypertimepoint_sex <- function(timepoint, mouse, df = dfanova, tab = mb) {
    set.seed(14)
    dfsel <- df %>% filter(Age_weeks == timepoint & Genotype == mouse)
    bray <- vegan::vegdist(tab[rownames(tab) %in% dfsel$ID,], method = 'bray')
    print(adonis2(bray ~ Sex, data = dfsel))
}

braypertimepoint_sex(timepoint = "8 weeks", mouse = "WT")
braypertimepoint_sex(timepoint = "10 weeks", mouse = "WT")
braypertimepoint_sex(timepoint = "12 weeks", mouse = "WT")
braypertimepoint_sex(timepoint = "14 weeks", mouse = "WT")
braypertimepoint_sex(timepoint = "16 weeks", mouse = "WT")
braypertimepoint_sex(timepoint = "18 weeks", mouse = "WT")

# Sex differences in TDP43
dfsel <- df %>% filter(Genotype == "TDP43" & Age_ints %in% c(6,10,12,14,16,18))
bray <- vegan::vegdist(mb[rownames(mb) %in% dfsel$ID,], method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
head(expl_variance_bray)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$ID <- rownames(dbray)
dbray <- left_join(dbray, dfsel, by = 'ID') # add metadata / covariates

braypertimepoint_sex <- function(timepoint, mouse, df = dfanova, tab = mb) {
    set.seed(14)
    dfsel <- df %>% filter(Age_weeks == timepoint & Genotype == mouse)
    bray <- vegan::vegdist(tab[rownames(tab) %in% dfsel$ID,], method = 'bray')
    return(adonis2(bray ~ Sex, data = dfsel))
}

tdp43_wk6 <- c(braypertimepoint_sex(timepoint = "6 weeks", mouse = "TDP43")['Pr(>F)'][[1]][1], "6 weeks")
tdp43_wk8 <- c(braypertimepoint_sex(timepoint = "8 weeks", mouse = "TDP43")['Pr(>F)'][[1]][1], "8 weeks")
tdp43_wk10 <- c(braypertimepoint_sex(timepoint = "10 weeks", mouse = "TDP43")['Pr(>F)'][[1]][1], "10 weeks")
tdp43_wk12 <- c(braypertimepoint_sex(timepoint = "12 weeks", mouse = "TDP43")['Pr(>F)'][[1]][1], "12 weeks")
tdp43_wk14 <- c(braypertimepoint_sex(timepoint = "14 weeks", mouse = "TDP43")['Pr(>F)'][[1]][1], "14 weeks")
tdp43_wk16 <- c(braypertimepoint_sex(timepoint = "16 weeks", mouse = "TDP43")['Pr(>F)'][[1]][1], "16 weeks")
tdp43_wk18 <- c(braypertimepoint_sex(timepoint = "18 weeks", mouse = "TDP43")['Pr(>F)'][[1]][1], "18 weeks")
res_tdp43 <- rbind(tdp43_wk6, tdp43_wk8, tdp43_wk10, tdp43_wk12, tdp43_wk14, tdp43_wk16, tdp43_wk18)
colnames(res_tdp43) <- c("pvalue", "Age_weeks")
res_tdp43 <- as.data.frame(res_tdp43)
res_tdp43$Age_weeks <- factor(res_tdp43$Age_weeks, levels = c("6 weeks", "8 weeks", "10 weeks", "12 weeks", "14 weeks", "16 weeks", "18 weeks"))
res_tdp43$pvalue <- as.numeric(res_tdp43$pvalue)

(braycurt <- dbray %>%
                ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
                    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
                    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
                    scale_color_manual(values = pal_nejm()(2)) +
                    scale_fill_manual(values = pal_nejm()(2), guide = "none") +
                    theme_Publication() +
                    labs(color = "", fill = "", title = "PCoA Bray-Curtis distance TDP43") +
                    stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), type = "norm",
                    alpha = 0.1, linewidth = 1.0) +
                    theme(legend.position = "top") +
                    geom_text(data = res_tdp43, aes(x = Inf, y = Inf, label = paste0("p = ", round(pvalue, 3))),
                              hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
                    facet_wrap(~ Age_weeks))
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
