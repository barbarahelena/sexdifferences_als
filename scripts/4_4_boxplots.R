## Boxplots metabolites - Sex x Genotype interaction
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(ggsci)

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
                legend.key.size= unit(0.4, "cm"),
                legend.spacing  = unit(0, "cm"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))

}

### Data ###
met <- readRDS("data/metabolome/metabolomics.RDS")
meta <- readRDS("data/metadata.RDS")

# Load interaction results from heatmap script
interaction_results <- rio::import("results/metabolomics/ttests/interaction_lm_results.csv")

# Metabolites with significant Sex x Genotype interaction (p < 0.05)
sig_metabolites <- interaction_results %>%
  filter(p_interaction < 0.05) %>%
  pull(metabolite)

# Merge metabolomics with metadata
df <- met %>%
  rownames_to_column("ID") %>%
  left_join(meta %>% select(ID, Sex, Intervention), by = "ID") %>%
  mutate(Sex = fct_rev(Sex),
         Intervention = factor(Intervention, levels = c("WT", "TDP43")))

### Boxplots for interaction metabolites ###
res_box <- list()
for (a in seq_along(sig_metabolites)) {
  metname <- sig_metabolites[a]
  print(metname)

  p_int <- interaction_results$p_interaction[interaction_results$metabolite == metname]
  p_int_label <- paste0("p(interaction) = ", format.pval(p_int, digits = 2))

  dfmet <- df %>%
    select(ID, Sex, Intervention, value = all_of(metname)) %>%
    mutate(value_z = as.numeric(scale(value)))

  labely <- max(dfmet$value_z, na.rm = TRUE) * 1.1

  pl <- ggplot(dfmet, aes(x = Intervention, y = value_z, fill = Sex)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.9,
                 position = position_dodge(0.6)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.6),
                alpha = 0.75, size = 1.5) +
    scale_fill_manual(values = pal_nejm()(2)) +
    scale_color_manual(values = pal_nejm()(2), guide = "none") +
    labs(title = str_wrap(metname, width = 22),
         y = "Metabolite (z-score)", x = "",
         caption = p_int_label) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_Publication() +
    theme(plot.caption = element_text(size = 8),
          legend.title = element_text(face = "bold"))

  ggsave(str_c("results/metabolomics/boxplots/", metname, "_interaction.pdf"),
         plot = pl, width = 3.5, height = 5, device = "pdf")

  res_box[[a]] <- pl
}

(boxplots <- ggarrange(plotlist = res_box, nrow = 1, ncol = 5,
                       common.legend = TRUE, legend = "bottom"))
ggsave("results/metabolomics/boxplots/boxplots_interaction_met.pdf", width = 10, height = 4)
