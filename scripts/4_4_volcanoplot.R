## Volcano plot metabolomics sex differences
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)

## Functions
calculate_log2_fold_change <- function(df) {
  # Ensure the data frame is grouped by metabolite
  df <- df %>% group_by(metabolite) %>%
        mutate(log2fold = mean[Sex == "Female"] - mean[Sex == "Male"])
  return(df)
}

calculate_log2_fold_change2 <- function(df) {
  # Ensure the data frame is grouped by metabolite
  df <- df %>% group_by(metabolite) %>%
        mutate(log2fold = mean[Intervention == "TDP43"] - mean[Intervention == "WT"])
  return(df)
}

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
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
} 

## Female-male ctrl
res1 <- rio::import("results/metabolomics/ttests/metabolites_welcht_ctrl_sex_diff.csv")
res1 <- calculate_log2_fold_change(res1)
res1
df_sum <- res1 %>% filter(Sex == "Female") # otherwise everything double
head(df_sum)

## Volcano plot
df_sum <- df_sum %>% 
    mutate(group = case_when(
    log2fold > 0 ~ paste0("higher in female"),
    log2fold < 0 ~ paste0("lower in female")
), 
sigdir = case_when(
    sig != "" & group == "higher in female" ~ paste0("up"),
    sig != "" & group == "lower in female" ~ paste0("down"),
    sig == "" ~ paste0("no")
), 
    sigdir = as.factor(sigdir),
    delabel = case_when(
       log2fold > 2 | log2fold < 2 | sig != "" ~ paste0(metabolite),
       sig == "" ~ paste0("")
    )
    ) %>% 
    ungroup(.) 
df_sum

set.seed(1234)
(plwt <- ggplot(df_sum, aes(x = log2fold, y = -log10(p.value), color = group, label = delabel)) +
    theme_Publication() +
    theme(axis.title = element_text(size = rel(0.8))) +
    geom_hline(aes(yintercept = -log10(0.05)), color = "darkgrey", linetype = "dashed") +
    geom_vline(aes(xintercept = -2), color = "darkgrey", linetype = "dashed") +
    geom_vline(aes(xintercept = 2), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 0.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey70",
                    force = 1.5, max.overlaps = 10) +
    scale_color_manual(values = rev(c(ggsci::pal_lancet()(2))), guide = "none") +
    #scale_x_continuous(limits = c(-0.75, 0.75), n.breaks = 6) +
    labs(x = "Log2 fold change (Female - Male)",
         y = "-log10(p-value)",
         title = "WT mice: sex difference"))

ggsave("results/metabolomics/volcanoplots/volcanoplot_ctrl.pdf", width = 5, height = 5, device = "pdf")
ggsave("results/metabolomics/volcanoplots/volcanoplot_ctrl.svg", width = 5, height = 5, device = "svg")

## Female-male TDP43
res2 <- rio::import("results/metabolomics/ttests/metabolites_welcht_tdp_sex_diff.csv")
res2 <- calculate_log2_fold_change(res2)
res2
df_sum <- res2 %>% filter(Sex == "Female") # to deduplicate
head(df_sum)

## Volcano plot
df_sum <- df_sum %>% 
    mutate(group = case_when(
    log2fold > 0 ~ paste0("higher in female"),
    log2fold < 0 ~ paste0("lower in female")
), 
sigdir = case_when(
    sig != "" & group == "higher in female" ~ paste0("up"),
    sig != "" & group == "lower in female" ~ paste0("down"),
    sig == "" ~ paste0("no")
), 
    sigdir = as.factor(sigdir),
    delabel = case_when(
       sig != "" ~ paste0(metabolite),
       sig == "" ~ paste0("")
    )
    ) %>% 
    ungroup(.) 
df_sum

set.seed(1234)
(pltdp <- ggplot(df_sum, aes(x = log2fold, y = -log10(p.value), color = group, label = delabel)) +
    theme_Publication() +
    theme(axis.title = element_text(size = rel(0.8))) +
    geom_hline(aes(yintercept = -log10(0.05)), color = "darkgrey", linetype = "dashed") +
    # geom_vline(aes(xintercept = -2), color = "darkgrey", linetype = "dashed") +
    # geom_vline(aes(xintercept = 2), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 0.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey70",
                    force = 1.5, max.overlaps = 10) +
    scale_color_manual(values = rev(c(ggsci::pal_lancet()(2))), guide = "none") +
    #scale_x_continuous(limits = c(-1, 1), breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0)) +
    labs(x = "Log2 fold change (Female - Male)",
         y = "-log10(p-value)",
         title = "TDP43 mice: sex difference"))

ggsave("results/metabolomics/volcanoplots/volcanoplot_tdp43.pdf", width = 5, height = 5, device = "pdf")
ggsave("results/metabolomics/volcanoplots/volcanoplot_tdp43.svg", width = 5, height = 5, device = "svg")

## TDP43 - Control (male and female together)
res3 <- rio::import("results/metabolomics/ttests/metabolites_welcht_genotype_diff.csv")
res3 <- calculate_log2_fold_change2(res3)
res3
df_sum <- res3 %>% filter(Intervention == "TDP43") # to deduplicate
head(df_sum)

## Volcano plot
df_sum <- df_sum %>% 
    mutate(group = case_when(
    log2fold > 0 ~ paste0("higher in TDP43"),
    log2fold < 0 ~ paste0("lower in TDP43")
), 
sigdir = case_when(
    sig != "" & group == "higher in TDP43" ~ paste0("up"),
    sig != "" & group == "lower in TDP43" ~ paste0("down"),
    sig == "" ~ paste0("no")
), 
    sigdir = as.factor(sigdir),
    delabel = case_when(
       sig != "" ~ paste0(metabolite),
       sig == "" ~ paste0("")
    )
    ) %>% 
    ungroup(.) 
df_sum

set.seed(1234)
(plgeno <- ggplot(df_sum, aes(x = log2fold, y = -log10(p.value), color = group, label = delabel)) +
    theme_Publication() +
    theme(axis.title = element_text(size = rel(0.8))) +
    geom_hline(aes(yintercept = -log10(0.05)), color = "darkgrey", linetype = "dashed") +
    # geom_vline(aes(xintercept = -2), color = "darkgrey", linetype = "dashed") +
    # geom_vline(aes(xintercept = 2), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 0.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey70",
                    force = 1.5, max.overlaps = 10) +
    scale_color_manual(values = rev(c(ggsci::pal_lancet()(2))), guide = "none") +
    #scale_x_continuous(limits = c(-1, 1), breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0)) +
    labs(x = "Log2 fold change (TDP43 - Control)",
         y = "-log10(p-value)",
         title = "All mice: genotype"))

ggsave("results/metabolomics/volcanoplots/volcanoplot_all.pdf", width = 5, height = 5, device = "pdf")
ggsave("results/metabolomics/volcanoplots/volcanoplot_all.svg", width = 5, height = 5, device = "svg")

## TDP43 - Control (female)
res4 <- rio::import("results/metabolomics/ttests/metabolites_welcht_fem_diff.csv")
res4 <- calculate_log2_fold_change2(res4)
res4
df_sum <- res4 %>% filter(Intervention == "TDP43")
head(df_sum)

## Volcano plot
df_sum <- df_sum %>% 
    mutate(group = case_when(
    log2fold > 0 ~ paste0("higher in TDP43"),
    log2fold < 0 ~ paste0("lower in TDP43")
), 
sigdir = case_when(
    sig != "" & group == "higher in TDP43" ~ paste0("up"),
    sig != "" & group == "lower in TDP43" ~ paste0("down"),
    sig == "" ~ paste0("no")
), 
    sigdir = as.factor(sigdir),
    delabel = case_when(
       sig != "" ~ paste0(metabolite),
       sig == "" ~ paste0("")
    )
    ) %>% 
    ungroup(.) 
df_sum

set.seed(1234)
ggplot(df_sum, aes(x = log2fold, y = -log10(p.value), color = group, label = delabel)) +
    theme_Publication() +
    theme(axis.title = element_text(size = rel(0.8))) +
    geom_hline(aes(yintercept = -log10(0.05)), color = "darkgrey", linetype = "dashed") +
    # geom_vline(aes(xintercept = -2), color = "darkgrey", linetype = "dashed") +
    # geom_vline(aes(xintercept = 2), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 0.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey70",
                    force = 1.5, max.overlaps = 10) +
    scale_color_manual(values = c(ggsci::pal_lancet()(2)), guide = "none") +
    #scale_x_continuous(limits = c(-1, 1), breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0)) +
    labs(x = "Log2 fold change (TDP43 - Control)",
         y = "-log10(p-value)")

ggsave("results/metabolomics/volcanoplots/volcanoplot_tdp43_fem.pdf", width = 5, height = 5, device = "pdf")
ggsave("results/metabolomics/volcanoplots/volcanoplot_tdp43_fem.svg", width = 5, height = 5, device = "svg")

## TDP43 - Control (male)
res5 <- rio::import("results/metabolomics/ttests/metabolites_welcht_male_diff.csv")
res5 <- calculate_log2_fold_change2(res5)
res5
df_sum <- res5 %>% filter(Intervention == "TDP43")
head(df_sum)

## Volcano plot
df_sum <- df_sum %>% 
    mutate(group = case_when(
    log2fold > 0 ~ paste0("higher in TDP43"),
    log2fold < 0 ~ paste0("lower in TDP43")
), 
sigdir = case_when(
    sig != "" & group == "higher in TDP43" ~ paste0("up"),
    sig != "" & group == "lower in TDP43" ~ paste0("down"),
    sig == "" ~ paste0("no")
), 
    sigdir = as.factor(sigdir),
    delabel = case_when(
       sig != "" ~ paste0(metabolite),
       sig == "" ~ paste0("")
    )
    ) %>% 
    ungroup(.) 
df_sum

set.seed(1234)
ggplot(df_sum, aes(x = log2fold, y = -log10(p.value), color = group, label = delabel)) +
    theme_Publication() +
    theme(axis.title = element_text(size = rel(0.8))) +
    geom_hline(aes(yintercept = -log10(0.05)), color = "darkgrey", linetype = "dashed") +
    # geom_vline(aes(xintercept = -2), color = "darkgrey", linetype = "dashed") +
    # geom_vline(aes(xintercept = 2), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 0.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey70",
                    force = 1.5, max.overlaps = 10) +
    scale_color_manual(values = c(ggsci::pal_lancet()(2)), guide = "none") +
    #scale_x_continuous(limits = c(-1, 1), breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0)) +
    labs(x = "Log2 fold change (TDP43 - Control)",
         y = "-log10(p-value)")

ggsave("results/metabolomics/volcanoplots/volcanoplot_tdp43_male.pdf", width = 5, height = 5, device = "pdf")
ggsave("results/metabolomics/volcanoplots/volcanoplot_tdp43_male.svg", width = 5, height = 5, device = "svg")
