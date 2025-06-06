# Boxplots microbes in human cohort
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

# Library
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(car)

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

# Data
mb <- readRDS("data/human_cohort/microbiome.RDS")
df <- readRDS("data/human_cohort/metadata.RDS")  %>%
    mutate(ALS = as.numeric(ALS),
           Group = as.factor(Group),
           ALScat = cut(ALS, breaks = c(0, 25, 35, 45, Inf), labels = c("0-25", "25-35", "35-45", ">45"), right = FALSE),
           ALScat = as.factor(ALScat),
           Age = as.numeric(Age)
    )
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
microb <- names(tot[1:14,])[which(colSums(tot[1:14]) / nrow(tot) > 0.01)]
tot <- tot %>% select(all_of(microb), ID:ALScat)
head(tot)

# Boxplots for each bacterium
## 1: Flavonifractor plautii
(sexdiff1 <- wilcox.test(`Flavonifractor plautii` ~ Sex, data = tot %>% filter(Group == "ALS")))
(sexdiff2 <- wilcox.test(`Flavonifractor plautii` ~ Sex, data = tot %>% filter(Group == "Healthy")))

# Format p-values
p_als <- format.pval(sexdiff1$p.value, digits = 2)
p_healthy <- format.pval(sexdiff2$p.value, digits = 2)

# Create subtitle with sex difference p-values
sex_diff_text <- paste0("Sex diff ALS p = ", p_als, "; Sex diff Healthy p = ", p_healthy)

(pl1 <- ggplot(tot, aes(x = Group, y = log10(`Flavonifractor plautii` + 0.5))) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
    scale_fill_manual(values = pal_nejm()(8)[c(3,6)], guide = "none") +
    labs(title = "Flavonifractor plautii", 
         caption = sex_diff_text,
         y = "Flavonifractor plautii", 
         x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
            plot.caption = element_text(size = 12),
            legend.position = "none") + 
    facet_wrap(~ Sex))
ggsave("results/humancohort/flavonifractor.pdf", width = 5, height = 5.5)

## 2: Bacteroides acidifaciens
(sexdiff1 <- wilcox.test(`Bacteroides acidifaciens` ~ Sex, data = tot %>% filter(Group == "ALS")))
(sexdiff2 <- wilcox.test(`Bacteroides acidifaciens` ~ Sex, data = tot %>% filter(Group == "Healthy")))

# Format p-values
p_als <- format.pval(sexdiff1$p.value, digits = 2)
p_healthy <- format.pval(sexdiff2$p.value, digits = 2)

# Create caption with sex difference p-values
sex_diff_text <- paste0("Sex diff ALS p = ", p_als, "; Sex diff Healthy p = ", p_healthy)

(pl2 <- ggplot(tot, aes(x = Group, y = log10(`Bacteroides acidifaciens` + 0.5))) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
    scale_fill_manual(values = pal_nejm()(8)[c(3,6)], guide = "none") +
    labs(title = "Bacteroides acidifaciens", 
         caption = sex_diff_text,
         y = "Bacteroides acidifaciens", 
         x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
            plot.caption = element_text(size = 12),
            legend.position = "none") + 
    facet_wrap(~ Sex))

ggsave("results/humancohort/bacteroidesacidi.pdf", width = 5, height = 5.5)

## 3: Phocaeicola sartorii
(sexdiff1 <- wilcox.test(`Phocaeicola sartorii` ~ Sex, data = tot %>% filter(Group == "ALS")))
(sexdiff2 <- wilcox.test(`Phocaeicola sartorii` ~ Sex, data = tot %>% filter(Group == "Healthy")))

# Format p-values
p_als <- format.pval(sexdiff1$p.value, digits = 2)
p_healthy <- format.pval(sexdiff2$p.value, digits = 2)

# Create caption with sex difference p-values
sex_diff_text <- paste0("Sex diff ALS p = ", p_als, "; Sex diff Healthy p = ", p_healthy)

(pl3 <- ggplot(tot, aes(x = Group, y = log10(`Phocaeicola sartorii` + 0.2))) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
    scale_fill_manual(values = pal_nejm()(8)[c(3,6)], guide = "none") +
    labs(title = "Phocaeicola sartorii", 
         caption = sex_diff_text,
         y = "Phocaeicola sartorii", 
         x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
            plot.caption = element_text(size = 12),
            legend.position = "none") + 
    facet_wrap(~ Sex))
ggsave("results/humancohort/phocaeicola.pdf", width = 5, height = 5.5)

## 4: Parabacteroides distasonis
(sexdiff1 <- wilcox.test(`Parabacteroides distasonis` ~ Sex, data = tot %>% filter(Group == "ALS")))
(sexdiff2 <- wilcox.test(`Parabacteroides distasonis` ~ Sex, data = tot %>% filter(Group == "Healthy")))

# Format p-values
p_als <- format.pval(sexdiff1$p.value, digits = 2)
p_healthy <- format.pval(sexdiff2$p.value, digits = 2)

# Create caption with sex difference p-values
sex_diff_text <- paste0("Sex diff ALS p = ", p_als, "; Sex diff Healthy p = ", p_healthy)

(pl4 <- ggplot(tot, aes(x = Group, y = log10(`Parabacteroides distasonis` + 0.1))) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
    scale_fill_manual(values = pal_nejm()(8)[c(3,6)], guide = "none") +
    labs(title = "Parabacteroides distasonis", 
         caption = sex_diff_text,
         y = "Parabacteroides distasonis", 
         x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
            plot.caption = element_text(size = 12),
            legend.position = "none") + 
    facet_wrap(~ Sex))
ggsave("results/humancohort/parabacteroides.pdf", width = 5, height = 5.5)

## 5: Parabacteroides goldsteinii
(sexdiff1 <- wilcox.test(`Parabacteroides goldsteinii` ~ Sex, data = tot %>% filter(Group == "ALS")))
(sexdiff2 <- wilcox.test(`Parabacteroides goldsteinii` ~ Sex, data = tot %>% filter(Group == "Healthy")))

# Format p-values
p_als <- format.pval(sexdiff1$p.value, digits = 2)
p_healthy <- format.pval(sexdiff2$p.value, digits = 2)

# Create caption with sex difference p-values
sex_diff_text <- paste0("Sex diff ALS p = ", p_als, "; Sex diff Healthy p = ", p_healthy)

(pl5 <- ggplot(tot, aes(x = Group, y = log10(`Parabacteroides goldsteinii` + 0.1))) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
    scale_fill_manual(values = pal_nejm()(8)[c(3,6)], guide = "none") +
    labs(title = "Parabacteroides goldsteinii", 
         caption = sex_diff_text,
         y = "Parabacteroides goldsteinii", 
         x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
            plot.caption = element_text(size = 12),
            legend.position = "none") + 
    facet_wrap(~ Sex))
ggsave("results/humancohort/parabacteroidesgold.pdf", width = 5, height = 5.5)

## 6: uncultured Prevotella sp.
(sexdiff1 <- wilcox.test(`uncultured Prevotella sp.` ~ Sex, data = tot %>% filter(Group == "ALS")))
(sexdiff2 <- wilcox.test(`uncultured Prevotella sp.` ~ Sex, data = tot %>% filter(Group == "Healthy")))

# Format p-values
p_als <- format.pval(sexdiff1$p.value, digits = 2)
p_healthy <- format.pval(sexdiff2$p.value, digits = 2)

# Create caption with sex difference p-values
sex_diff_text <- paste0("Sex diff ALS p = ", p_als, "; Sex diff Healthy p = ", p_healthy)

(pl6 <- ggplot(tot, aes(x = Group, y = log10(`uncultured Prevotella sp.` + 0.1))) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
    scale_fill_manual(values = pal_nejm()(8)[c(3,6)], guide = "none") +
    labs(title = "uncultured Prevotella sp.", 
         caption = sex_diff_text,
         y = "uncultured Prevotella sp.", 
         x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
            plot.caption = element_text(size = 12),
            legend.position = "none") + 
    facet_wrap(~ Sex))
ggsave("results/humancohort/prevotella.pdf", width = 5, height = 5.5)

## 7: Akkermansia muciniphila
(sexdiff1 <- wilcox.test(`Akkermansia muciniphila` ~ Sex, data = tot %>% filter(Group == "ALS")))
(sexdiff2 <- wilcox.test(`Akkermansia muciniphila` ~ Sex, data = tot %>% filter(Group == "Healthy")))

# Format p-values
p_als <- format.pval(sexdiff1$p.value, digits = 2)
p_healthy <- format.pval(sexdiff2$p.value, digits = 2)

# Create caption with sex difference p-values
sex_diff_text <- paste0("Sex diff ALS p = ", p_als, "; Sex diff Healthy p = ", p_healthy)

(pl7 <- ggplot(tot, aes(x = Group, y = log10(`Akkermansia muciniphila` + 0.5))) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
    scale_fill_manual(values = pal_nejm()(8)[c(3,6)], guide = "none") +
    labs(title = "Akkermansia muciniphila", 
         caption = sex_diff_text,
         y = "Akkermansia muciniphila", 
         x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
            plot.caption = element_text(size = 12),
            legend.position = "none") + 
    facet_wrap(~ Sex))
ggsave("results/humancohort/akkermansia.pdf", width = 5, height = 5.5)

# (boxplots <- ggarrange(pl7, pl2, pl1, pl4, pl5, pl6, pl3, labels = LETTERS[1:7],
#           ncol = 2, nrow = 4,
#           common.legend = TRUE, legend = "bottom"))
# ggsave("results/humancohort/boxplots_microbes.pdf", width = 10, height = 18)

(boxplots <- ggarrange(pl7, pl2, pl1, pl4, pl5, pl6, pl3, labels = LETTERS[1:7],
          ncol = 4, nrow = 2,
          common.legend = TRUE, legend = "bottom"))
ggsave("results/humancohort/boxplots_microbes.pdf", width = 18, height = 10)
