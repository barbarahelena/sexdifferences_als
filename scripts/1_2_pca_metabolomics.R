# Descriptives metabolomics
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

#### Libraries ####
library(tidyverse)
library(ggsci)
library(ggpubr)

##### Functions ####
theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                #family = 'Helvetica'
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line.x = element_line(colour="black"),
                axis.line.y = element_line(colour="black"),
                axis.ticks.x = element_line(),
                axis.ticks.y = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(5,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(face = "italic", size=rel(0.6))
        ))
} 

#### Opening data files ####
als <- readRDS("data/metadata.RDS") %>% mutate(Sex = fct_rev(Sex))
met <- readRDS("data/metabolome/metabolomics.RDS")

#### PCA plot metabolites ####
head(met)[1:5,1:5]
dim(met)
met <- as.matrix(met)
set.seed(112)
pca <- prcomp(met, center = TRUE, scale = TRUE)
df <- as.data.frame(pca$x[, 1:5])
summary(pca)
pc1_ev <- round(summary(pca)$importance[2,1] * 100, 2)
pc2_ev <- round(summary(pca)$importance[2,2] * 100, 2)
head(df)
df$ID <- rownames(df)
df$Group <- als$Intervention[match(df$ID, als$ID)]
df$Sex <- als$Sex[match(df$ID, als$ID)]
df$GroupPerSex <- als$GroupPerSex[match(df$ID, als$ID)]
df$ID <-als$ID[match(df$ID, als$ID)]
df <- df %>% arrange(Group)

(pca1 <- ggplot(df, aes(x=PC1, y=PC2, color=Group)) +
        stat_ellipse(geom = "polygon", aes(fill = Group, color = Group),
                     alpha = 0.3) +
        geom_point(size = 3, aes(shape = Sex)) +
        ggtitle('PCA plasma metabolites') +
        scale_color_manual(values = pal_nejm()(6)[c(6,3)]) +
        scale_fill_manual(values = pal_nejm()(6)[c(6,3)]) +
        theme_minimal() +
        labs(x=str_c('PC1 ', pc1_ev, '%'), y=str_c('PC2 ', pc2_ev, '%')) +
        theme_Publication() +
        theme(legend.title = element_blank()))
ggsave("results/metabolomics/pca/pca_groups.pdf", width = 5, height = 4.5)

(pca2 <- ggplot(df, aes(x=PC1, y=PC2, color=Sex)) +
        stat_ellipse(geom = "polygon", aes(fill = Sex, color = Sex),
                     alpha = 0.3) +
        geom_point(size = 3, aes(shape = Sex)) +
        ggtitle('PCA: sex differences') +
        scale_color_manual(values = pal_nejm()(2)) +
        scale_fill_manual(values = pal_nejm()(2)) +
        theme_minimal() +
        labs(x=str_c('PC1 ', pc1_ev, '%'), y=str_c('PC2 ', pc2_ev, '%')) +
        facet_wrap(~Group, nrow = 1) +
        theme_Publication() +
        theme(legend.title = element_blank()))
ggsave("results/metabolomics/pca/pca_sex.pdf", width = 8, height = 4.5)
