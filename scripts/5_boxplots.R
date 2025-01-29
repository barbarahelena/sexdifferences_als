## Boxplots metabolites
## Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Libraries
library(tidyverse)
library(ggpubr)
library(gridExtra)

theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(family = 'Helvetica'),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
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
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
} 

### Data ###
df3_scale <- readRDS("data/proteins_expression_scaled_data.RDS") # all protein expression data (scaled to mean 0 and sd 1)
df5 <- readRDS("data/proteins_ttest_diff.RDS")
qval_sig <- "Gcg"
pval_sig <- unique(df5_sel$protein)

### Boxplot differences ###
res_box <- list()
for(a in 1:length(pval_sig)){
    protein <- pval_sig[a]
    proteinname <- paste0(pval_sig[a])
    print(proteinname)
    df_protein <- df3_scale %>% select(Treatment, all_of(protein))
    df_protein$protein_y <- df_protein[,2]
    comp <- list(c("Control + CI", "Pneumonia + CI"))
    (pl <- ggplot(data = df_protein, aes(x = Treatment, y = protein_y)) + 
            ggpubr::stat_compare_means(method = "t.test", label = "p.signif", comparisons = comp,
                                       hide.ns = TRUE, bracket.size = 0.5, size = 5) +
            geom_boxplot(aes(color = Treatment), outlier.shape = NA, 
                         width = 0.5, alpha = 0.9) +
            geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
            scale_color_manual(guide = "none", values = c("black","red")) +
            ylim(NA, max(df_protein$protein_y)*1.3) +
            labs(title=proteinname, y="Protein expression (z-score)") +
            theme_Publication())
    ggsave(str_c("results/pdf/proteins/", proteinname, ".pdf"), width = 3, height = 4, device = "pdf")
    ggsave(str_c("results/svg/proteins/", proteinname, ".svg"), width = 3, height = 4, device = "svg")
    res_box[[a]] <- pl
}

pdf("results/pdf/boxplots_proteins.pdf", width = 12, height = 10)
gridExtra::grid.arrange(grobs=res_box, ncol=4)
dev.off()

svg("results/svg/boxplots_proteins.svg", width = 12, height = 10)
gridExtra::grid.arrange(grobs=res_box, ncol=5)
dev.off()

