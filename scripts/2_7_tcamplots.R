# TCAM plots ALS project
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

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
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
} 


plot_loadings <- function(data, title) {
  data$comp <- data[[2]]
  ggplot(data, aes(x = reorder(X, comp), y = comp)) +
    geom_bar(stat = "identity", 
        fill = c(rep("firebrick", 10), rep("royalblue",10))) +
    coord_flip() +
    labs(title = title, x = "", y = "Loading") +
    theme_Publication()
}

## Load data
tcam <- read.csv("results/microbiome/tcam/pythonoutput/df_plot.csv") %>% 
  select(1:17) %>% 
  mutate(Genotype = fct_relevel(Genotype, "WT", after = 0L),
          GenotypePerSex = fct_relevel(GenotypePerSex, "Female WT", after = 0L),
          GenotypePerSex = fct_relevel(GenotypePerSex, "Male WT", after = 1L))
head(tcam)
load <- read.csv("results/microbiome/tcam/pythonoutput/df_loadings.csv")
head(load)
mb <- readRDS("data/microbiome/imputed_microbiome_data.RDS")
meta <- readRDS("data/microbiome/meta_microbiome.RDS")

## PERMANOVAs ##
tcam <- tcam %>% filter(Age_ints == 6) 
f1n <- names(tcam)[7]
f2n <- names(tcam)[8]

set.seed(1234)
# 1. PERMANOVA for overall genotype effect
tcam_dist <- tcam %>% select(all_of(c(f1n, f2n))) %>% dist(method = "euclidean")
(permanova_genotype <- adonis2(tcam_dist ~ Genotype, data = tcam, permutations = 999, method = "euclidean", by = "term"))

# 2. PERMANOVA for genotype in females
tcam_fem <- tcam %>% filter(Sex == "Female")
tcam_fem_dist <- tcam_fem %>% select(all_of(c(f1n, f2n))) %>% dist(method = "euclidean")
(permanova_fem <- adonis2(tcam_fem_dist ~ Genotype, data = tcam_fem, permutations = 999, method = "euclidean", by = "term"))

# 3. PERMANOVA for genotype in males
tcam_male <- tcam %>% filter(Sex == "Male")
tcam_male_dist <- tcam_male %>% select(all_of(c(f1n, f2n))) %>% dist(method = "euclidean")
(permanova_male <- adonis2(tcam_male_dist ~ Genotype, data = tcam_male, permutations = 999, method = "euclidean", by = "term"))

# 3. PERMANOVA for genotype in males
tcam_dist <- tcam %>% select(all_of(c(f1n, f2n))) %>% dist(method = "euclidean")
(permanova_male <- adonis2(tcam_dist ~ Genotype * Sex, data = tcam, permutations = 999, method = "euclidean", by = "term"))

annotation_text <- str_c("all: p = ", permanova_genotype$`Pr(>F)`[1], "\n",
                    "females: p = ", permanova_fem$`Pr(>F)`[1], "\n",
                    "males: p = ", permanova_male$`Pr(>F)`[1])

### TCAM plot ###
(tcampl <- ggplot(data = tcam, aes(x = .data[[f1n]], y = .data[[f2n]], 
                        color = Genotype, fill = Genotype, shape = Sex)) +
        #stat_ellipse(geom = "polygon", alpha = 0.3) +
        geom_point(size = 5) +
        ggtitle('TCAM') +
        scale_color_manual(values = pal_nejm()(6)[c(6,3)]) +
        scale_fill_manual(values = pal_nejm()(6)[c(6,3)]) +
        theme_minimal() +
        labs(x=str_c(str_replace(f1n, "[.]", " "), "%"),
            y=str_c(str_replace(f2n, "[.]", " "), "%")) +
        annotate("text", x = Inf, y = Inf, 
                label = annotation_text,
                hjust = 1, vjust = 1, 
                size = 4) +
        theme_Publication() +
        theme(legend.title = element_blank()))
ggsave("results/microbiome/tcam/routput/f1f2_scatter.pdf", width = 6, height = 6) 


### PERMANOVA ###
# 1. PERMANOVA for Sex effect within WT
tcam_wt <- tcam %>% filter(Genotype == "WT")
tcam_wt_dist <- tcam_wt %>% select(all_of(c(f1n, f2n))) %>%dist(method = "euclidean")
(permanova_sex_wt <- adonis2(tcam_wt_dist ~ Sex, data = tcam_wt, permutations = 999, method = "euclidean"))

# 2. PERMANOVA for Sex effect within TDP-43
tcam_tdp <- tcam %>% filter(Genotype == "TDP43") 
tcam_tdp_dist <- tcam_tdp %>% select(all_of(c(f1n, f2n))) %>% dist(method = "euclidean")
(permanova_sex_tdp <- adonis2(tcam_tdp_dist ~ Sex, data = tcam_tdp, permutations = 999, method = "euclidean"))

pvalues <- c(permanova_sex_wt$`Pr(>F)`[1], permanova_sex_tdp$`Pr(>F)`[1])
genotype <- c("WT", "TDP43")
res <- data.frame(Genotype = genotype, pvalue = pvalues)

## TCAM plot ##
(tcamgenotype <- ggplot(data = tcam, aes(x = .data[[f1n]], y = .data[[f2n]], 
                        color = Sex, fill = Sex, shape = Sex)) +
        stat_ellipse(geom = "polygon", alpha = 0.3) +
        geom_point(size = 3, alpha = 0.75) +
        ggtitle('TCAM: sex differences') +
        scale_color_manual(values = pal_nejm()(2)) +
        scale_fill_manual(values = pal_nejm()(2)) +
        theme_minimal() +
        labs(x=str_c(str_replace(f1n, "[.]", " "), "%"),
            y=str_c(str_replace(f2n, "[.]", " "), "%")) +
        geom_text(data = res, aes(x = Inf, y = Inf, label = paste0("p = ", round(pvalue, 3))),
                              hjust = 1.1, vjust = 1.1, size = 4, inherit.aes = FALSE) +
        theme_Publication() +
        facet_wrap(~Genotype) +
        theme(legend.title = element_blank()))
ggsave(tcamgenotype, filename = "results/microbiome/tcam/routput/f1f2_scatter_pergenotype.pdf", width = 6, height = 6)   

## Loading of component plots ##
head(load)
names(load)[1:10]
last <- nrow(load)
min10 <- last - 9
f1 <- load %>% select(X, all_of(f1n)) %>% arrange(-.data[[f1n]])
f1 <- f1[c(1:10, min10:last),]
f2 <- load %>% select(X, all_of(f2n)) %>% arrange(-.data[[f2n]])
f2 <- f2[c(1:10, min10:last),]

(f1plot <- plot_loadings(f1, "TCAM F1"))
ggsave("results/microbiome/tcam/routput/loading_pc1.pdf", width = 8, height = 7)

plot_loadings(f2, "TCAM F2")
ggsave("results/microbiome/tcam/routput/loading_pc2.pdf", width = 8, height = 7)

## Lineplots ##
f1 <- load %>% select(X, all_of(f1n)) %>% arrange(-.data[[f1n]])
f1 <- f1[c(1:3, (last-1):last),]
pseudocounts <- mb %>%  select(any_of(f1$X)) %>%
  summarise(across(everything(), ~ min(.x[.x > 0], na.rm = TRUE) / 2))
mbsel <- mb %>% ungroup(.) %>% select(any_of(f1$X), ID, Genotype, Age_ints, Sex, GenotypePerSex) %>%
  mutate( across(f1$X, ~log10(.x + pseudocounts[[cur_column()]]))
          ) %>%
  rename_at(c(f1$X), ~str_remove(.x, " \\(.*\\)$"))

means_ses <- mbsel %>% # Calculate means and standard deviations per Genotype
  group_by(Genotype, Age_ints) %>%
  summarise(across(1:5, list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~(sd(.x, na.rm = TRUE) / sqrt(n())))))
means_ses

means_ses_long <- means_ses %>% # to long format for plotting
  pivot_longer(cols = c(-Genotype, -Age_ints), 
                        names_to = c("microbe", ".value"), names_sep = "_")
mbsel_long <- mbsel %>%
  pivot_longer(cols = c(-ID, -Genotype, -Age_ints, -GenotypePerSex, -Sex), 
                        names_to = "microbe", values_to = "value")

genotype_lineplots <- ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                                    group = Genotype, color = Genotype)) +
                      geom_line() +
                      geom_point(size = 0.5) +
                      geom_jitter(data = mbsel_long, aes(x = Age_ints, y = value, color = Genotype),
                                    width = 0.1, alpha = 0.7) +
                      scale_color_manual(values = pal_nejm()(8)[c(3,6)]) +
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                      labs(title = "Abundance over time - TDP43 vs WT",
                          x = "Age (weeks)",
                          y = "log10(relative abundance)") +
                      scale_x_continuous(breaks = c(6,8,10,12,14,16)) +
                      facet_wrap(~ microbe, scales = "free", nrow = 1) +
                      theme_Publication() +
                      theme(strip.text = element_text(size = 8))
ggsave("results/microbiome/tcam/routput/lineplots_genotype_log10.pdf", width = 15, height = 5)

# The same but then for sex differences within TDP43
means_ses <- mbsel %>% filter(Genotype == "TDP43") %>% ## Calculate means and standard deviations
  group_by(Sex, Age_ints) %>%
  summarise(across(1:5, list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~(sd(.x, na.rm = TRUE) / sqrt(n())))))
means_ses

means_ses_long <- means_ses %>% # to long
  pivot_longer(cols = c(-Sex, -Age_ints), 
                        names_to = c("microbe", ".value"), names_sep = "_")
mbsel_long <- mbsel %>%
  pivot_longer(cols = c(-ID, -Genotype, -Age_ints, -GenotypePerSex, -Sex), 
                        names_to = "microbe", values_to = "value")

(sexdiff_lineplots <- ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                                    group = Sex, color = Sex)) +
                      geom_line() +
                      geom_jitter(data = mbsel_long %>% filter(Genotype == "TDP43"), aes(x = Age_ints, y = value, 
                                    color = Sex),
                                    width = 0.1, alpha = 0.7) +
                      scale_color_nejm() +
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                      labs(title = "Abundance over time for TDP43 - males vs females",
                          x = "Age (weeks)",
                          y = "log10(relative abundance)") +
                      scale_x_continuous(breaks = c(6,8,10,12,14,16)) +
                      facet_wrap(~ microbe, scales = "free", nrow = 1) +
                      theme_Publication() +
                      theme(strip.text = element_text(size = 8)))
ggsave("results/microbiome/tcam/routput/lineplots_tdp43_sex.pdf", width = 15, height = 5)

# The same but then for sex differences within TDP43
means_ses <- mbsel %>% filter(Genotype == "WT") %>% ## Calculate means and standard deviations
  group_by(Sex, Age_ints) %>%
  summarise(across(1:5, list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~(sd(.x, na.rm = TRUE) / sqrt(n())))))
means_ses

means_ses_long <- means_ses %>% # to long
  pivot_longer(cols = c(-Sex, -Age_ints), 
                        names_to = c("microbe", ".value"), names_sep = "_")
mbsel_long <- mbsel %>%
  pivot_longer(cols = c(-ID, -Genotype, -Age_ints, -GenotypePerSex, -Sex), 
                        names_to = "microbe", values_to = "value")

(sexdiff_lineplots <- ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                                    group = Sex, color = Sex)) +
                      geom_line() +
                      geom_jitter(data = mbsel_long %>% filter(Genotype == "WT"), aes(x = Age_ints, y = value, 
                                    color = Sex),
                                    width = 0.1, alpha = 0.7) +
                      scale_color_nejm() +
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                      labs(title = "Abundance over time for WT - males vs females",
                          x = "Age (weeks)",
                          y = "log10(relative abundance)") +
                      scale_x_continuous(breaks = c(6,8,10,12,14,16)) +
                      facet_wrap(~ microbe, scales = "free", nrow = 1) +
                      theme_Publication() +
                      theme(strip.text = element_text(size = 8)))
ggsave("results/microbiome/tcam/routput/lineplots_wt_sex.pdf", width = 15, height = 5)

### Ggarrange of microbiome plots ###
(pl <- ggarrange(ggarrange(comp_species, braycurt, nrow = 1, labels = LETTERS[1:2]),
                  ggarrange(tcampl, tcamgenotype, f1plot, nrow = 1, widths = c(1, 1, 1.6), labels = LETTERS[3:5]),
                  genotype_lineplots,
                  sexdiff_lineplots,
                  ncol = 1, heights = c(1.3, 0.8, 0.8, 0.8), labels = c("", "", LETTERS[6:7])))
ggsave("results/microbiome/assembled_microbiome.pdf", width = 16, height = 20)

# 20 microbes with highest loadings
## Lineplots ##
f1 <- load %>% select(X, all_of(f1n)) %>% arrange(-.data[[f1n]])
f1 <- f1[c(1:10, (last-9):last),]
pseudocounts <- mb %>%  select(any_of(f1$X)) %>%
  summarise(across(everything(), ~ min(.x[.x > 0], na.rm = TRUE) / 2))
mbsel <- mb %>% ungroup(.) %>% select(any_of(f1$X), ID, Genotype, Age_ints, Sex, GenotypePerSex) %>%
  mutate( across(f1$X, ~log10(.x + pseudocounts[[cur_column()]]))
          ) %>%
  rename_at(c(f1$X), ~str_remove(.x, " \\(.*\\)$"))

means_ses <- mbsel %>% # Calculate means and standard deviations per Genotype
  group_by(Genotype, Age_ints) %>%
  summarise(across(1:20, list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~(sd(.x, na.rm = TRUE) / sqrt(n())))))
means_ses

means_ses_long <- means_ses %>% # to long format for plotting
  pivot_longer(cols = c(-Genotype, -Age_ints), 
                        names_to = c("microbe", ".value"), names_sep = "_")
mbsel_long <- mbsel %>%
  pivot_longer(cols = c(-ID, -Genotype, -Age_ints, -GenotypePerSex, -Sex), 
                        names_to = "microbe", values_to = "value")

genotype_lineplots20 <- ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                                    group = Genotype, color = Genotype)) +
                      geom_line() +
                      geom_point(size = 0.5) +
                      geom_jitter(data = mbsel_long, aes(x = Age_ints, y = value, color = Genotype),
                                    width = 0.1, alpha = 0.7) +
                      scale_color_manual(values = pal_nejm()(8)[c(3,6)]) +
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                      labs(title = "Abundance over time - TDP43 vs WT",
                          x = "Age (weeks)",
                          y = "log10(relative abundance)") +
                      scale_x_continuous(breaks = c(6,8,10,12,14,16)) +
                      facet_wrap(~ microbe, scales = "free", nrow = 4) +
                      theme_Publication() +
                      theme(strip.text = element_text(size = 8))
ggsave("results/microbiome/tcam/routput/20lineplots_genotype_log10.pdf", width = 18, height = 15)

# The same but then for sex differences within TDP43
means_ses <- mbsel %>% filter(Genotype == "TDP43") %>% ## Calculate means and standard deviations
  group_by(Sex, Age_ints) %>%
  summarise(across(1:20, list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~(sd(.x, na.rm = TRUE) / sqrt(n())))))
means_ses

means_ses_long <- means_ses %>% # to long
  pivot_longer(cols = c(-Sex, -Age_ints), 
                        names_to = c("microbe", ".value"), names_sep = "_")
mbsel_long <- mbsel %>%
  pivot_longer(cols = c(-ID, -Genotype, -Age_ints, -GenotypePerSex, -Sex), 
                        names_to = "microbe", values_to = "value")

(sexdiff_lineplots <- ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                                    group = Sex, color = Sex)) +
                      geom_line() +
                      geom_jitter(data = mbsel_long %>% filter(Genotype == "TDP43"), aes(x = Age_ints, y = value, 
                                    color = Sex),
                                    width = 0.1, alpha = 0.7) +
                      scale_color_nejm() +
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                      labs(title = "Abundance over time for TDP43 - males vs females",
                          x = "Age (weeks)",
                          y = "log10(relative abundance)") +
                      scale_x_continuous(breaks = c(6,8,10,12,14,16)) +
                      facet_wrap(~ microbe, scales = "free", nrow = 4) +
                      theme_Publication() +
                      theme(strip.text = element_text(size = 8)))
ggsave("results/microbiome/tcam/routput/20lineplots_tdp43_sex.pdf", width = 18, height = 15)

# The same but then for sex differences within WT
means_ses <- mbsel %>% filter(Genotype == "WT") %>% ## Calculate means and standard deviations
  group_by(Sex, Age_ints) %>%
  summarise(across(1:20, list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~(sd(.x, na.rm = TRUE) / sqrt(n())))))
means_ses

means_ses_long <- means_ses %>% # to long
  pivot_longer(cols = c(-Sex, -Age_ints), 
                        names_to = c("microbe", ".value"), names_sep = "_")
mbsel_long <- mbsel %>%
  pivot_longer(cols = c(-ID, -Genotype, -Age_ints, -GenotypePerSex, -Sex), 
                        names_to = "microbe", values_to = "value")

(sexdiff_lineplots <- ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                                    group = Sex, color = Sex)) +
                      geom_line() +
                      geom_jitter(data = mbsel_long %>% filter(Genotype == "WT"), aes(x = Age_ints, y = value, 
                                    color = Sex),
                                    width = 0.1, alpha = 0.7) +
                      scale_color_nejm() +
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                      labs(title = "Abundance over time for WT - males vs females",
                          x = "Age (weeks)",
                          y = "log10(relative abundance)") +
                      scale_x_continuous(breaks = c(6,8,10,12,14,16)) +
                      facet_wrap(~ microbe, scales = "free", nrow = 4) +
                      theme_Publication() +
                      theme(strip.text = element_text(size = 8)))
ggsave("results/microbiome/tcam/routput/20lineplots_wt_sex.pdf", width = 18, height = 15)

## Boxplots top 10 F1 microbes ##
f1 <- load %>% select(X, all_of(f1n)) %>% arrange(-.data[[f1n]])
f1_top10 <- f1[c(1:10, (last-9):last),]
pseudocounts_box <- mb %>% select(any_of(f1_top10$X)) %>%
  summarise(across(everything(), ~ min(.x[.x > 0], na.rm = TRUE) / 2))
mbsel_box <- mb %>% ungroup(.) %>% select(any_of(f1_top10$X), ID, Genotype, Age_ints, Sex, GenotypePerSex) %>%
  filter(Age_ints == 12) |> 
  mutate(across(f1_top10$X, ~log10(.x + pseudocounts_box[[cur_column()]]))) %>%
  rename_at(c(f1_top10$X), ~str_remove(.x, " \\(.*\\)$"))

microbe_names <- str_remove(f1_top10$X, " \\(.*\\)$")

res_box <- list()
for(i in seq_along(microbe_names)) {
  mic <- microbe_names[i]
  print(mic)
  dfmic <- mbsel_box %>% select(Sex, Genotype, all_of(mic)) %>%
    rename(value = all_of(mic))

  # Genotype difference per sex
  sexdiff_f <- wilcox.test(value ~ Genotype, data = dfmic %>% filter(Sex == "Female"), var.equal = FALSE)
  sexdiff_m <- wilcox.test(value ~ Genotype, data = dfmic %>% filter(Sex == "Male"), var.equal = FALSE)
  p_female <- format.pval(sexdiff_f$p.value, digits = 2)
  p_male <- format.pval(sexdiff_m$p.value, digits = 2)
  geno_diff_text <- paste0("Genotype: p = ", p_female, " (F); p = ", p_male, " (M)")

  labely <- max(dfmic$value, na.rm = TRUE) * 1.05

  pl <- ggplot(data = dfmic, aes(x = Sex, y = value)) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 2.8, label.x = 1,
                       label.y = labely) +
    geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.5, alpha = 0.9) +
    geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
    scale_fill_manual(guide = "none", values = pal_nejm()(2)) +
    labs(title = str_wrap(mic, width = 25), y = "log10(relative abundance)", x = "",
         caption = geno_diff_text) +
    facet_wrap(~Genotype) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.caption = element_text(size = 7))

  res_box[[i]] <- pl
}

pdf("results/microbiome/tcam/routput/boxplots_top10_f1.pdf", width = 20, height = 25)
gridExtra::grid.arrange(grobs = res_box, ncol = 5)
dev.off()

## Lactobacillus rhamnosus
df <- readRDS("data/microbiome/microbiome_run1.RDS")
lac <- df %>% select(contains("Lactobacillus"))
lac <- lac[,apply(lac, 2, mean) > 0.001]
dim(lac)
names(lac)
lac$ID <- rownames(lac)
lac2 <- left_join(lac, meta, by = "ID")
head(lac2)

lac2 <- lac2 %>% filter(Genotype == "TDP43") %>% filter(Age_ints <= 18)
plist <- list()
for(i in 1:7){
  lac2$var <- log10(lac2[[i]] + 0.001)
  print(names(lac2)[i])
  means <- lac2 %>% group_by(Sex, Age_ints) %>% summarise(mean = mean(var), sd = sd(var), se = sd / sqrt(nrow(.)))
  pl <- ggplot(means, aes(x = Age_ints, y = mean, group = Sex, color = Sex)) +
                      geom_line() +
                      geom_jitter(data = lac2, aes(x = Age_ints, y = var, 
                                    color = Sex),
                                    width = 0.1, alpha = 0.7) +
                      scale_color_nejm() +
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                      labs(title = names(lac2)[i], x = "Age (weeks)",
                          y = "log10(relative abundance)") +
                      scale_x_continuous(breaks = c(6,8,10,12,14,16,18)) +
                      theme_Publication() +
                      theme(strip.text = element_text(size = 8))
  plist[[i]] <- pl
}

ggarrange(plotlist = plist, ncol = 4, nrow = 2, labels = LETTERS[1:length(plist)], common.legend = TRUE)
ggsave("results/microbiome/lactobacilli.pdf", width = 15, height = 15)
