# Compositional plots ALS cohort 2
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Load libraries
library(tidyverse)
library(ggpubr)
library(circlize)
library(ggsci)

## Functions
cols <- colorRampPalette(c(ggsci::pal_nejm()(8)))

theme_composition <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                #axis.text.x =  element_text(angle = 45, hjust = 1),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                legend.key.size= unit(0.4, "cm"),
                legend.spacing  = unit(0, "cm"),
                legend.text = element_text(size = rel(0.7)),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
}

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
} 

# Data
mb <- readRDS("data/human_cohort2/microbiome.RDS")
head(mb)[1:5,1:5]
mean(rowSums(mb)); sd(rowSums(mb)) # adds up to 1, with 0 var
df <- readRDS("data/human_cohort2/metadata.RDS")
head(df)

# Species level summarised per group
mb2 <- mb %>% # convert to long format per species per sample
    rownames_to_column(var = 'Sample') %>% 
    pivot_longer(-Sample, names_to = 'Tax', values_to = 'Abundance') %>% 
    mutate(sampleID = Sample)

top_taxa <- mb2 %>% # summarise data and select top 20 taxa
    group_by(Tax, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>% 
    group_by(Tax) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    arrange(-Abundance) %>% 
    dplyr::select(Tax) %>% 
    head(20)
top_taxa

mb3 <- mb2 %>% 
    mutate(
        Tax2 = case_when(
            Tax %in% top_taxa$Tax ~ paste(Tax),
            .default = "Other species"
        ),
        Tax2 = as.factor(Tax2)
    ) %>% 
    group_by(Tax2, Sample) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    mutate(Group = df$Group[match(Sample, df$ID)],
           Sex = df$Sex[match(Sample, df$ID)],
           Age = df$Age[match(Sample, df$ID)] ) %>% 
    group_by(Tax2, Group, Sex) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    ungroup(.) %>% 
    mutate(
           Tax2 = fct_reorder(Tax2, Abundance),
           Tax2 = fct_relevel(Tax2, "Other species", after = 0L)
    )
lev <- levels(mb3$Tax2)

set.seed(1234)
(comp_species <- mb3 %>% 
    ggplot(aes(x = Group, y = Abundance, fill = Tax2)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = rev(c(sample(cols(20)), 
        "grey90")), labels = lev) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(y="Composition (relative abundances)", x = "", 
        title = "Microbiota composition", fill = "") +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~Sex) +
    theme_composition())
ggsave("results/humancohort2/compositionalplots/speciescomp.pdf", width = 6, height = 7)
