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
df <- readRDS("data/metabolomics.RDS") 
allmets <- names(df)
df$ID <- rownames(df)
meta <- readRDS("data/metadata.RDS")
df <- right_join(df, meta, by = "ID") %>% mutate(Sex = fct_rev(Sex))
tests <- rio::import("results/metabolomics/ttests/metabolites_welcht_diff.csv")
tests2 <- rio::import("results/metabolomics/ttests/metabolites_welcht_ctrl_diff.csv")
qval_sig <- unique(c(tests$metabolite[which(tests$q.value < 0.05)], tests2$metabolite[which(tests2$q.value < 0.05)]))
pval_sig <- unique(c(tests$metabolite[which(tests$p.value < 0.05)], tests2$metabolite[which(tests2$p.value < 0.05)]))

### Boxplot differences ###
res_box <- list()
for(a in 1:length(pval_sig)){
    metab <- pval_sig[a]
    metname <- paste0(pval_sig[a])
    print(metname)
    dfmet <- df %>% select(Sex, Intervention, all_of(metname))
    dfmet$met_y <- scale(dfmet[,3])
    # comp <- list(c("Female Control", "Male Control"), c("Female TDP43", "Male TDP43"))
    comp <- list(c("Female", "Male"))
    (pl <- ggplot(data = dfmet, aes(x = Sex, y = met_y)) + 
            ggpubr::stat_compare_means(method = "t.test", label = "p.signif", comparisons = comp,
                                       hide.ns = TRUE, size = 5, #step.increase = 0.1,
                                       var.equal = FALSE) +
            geom_boxplot(aes(fill = Sex), outlier.shape = NA, 
                         width = 0.5, alpha = 0.9) +
            geom_jitter(color = "grey5", height = 0, width = 0.1, alpha = 0.75) +
            scale_fill_manual(guide = "none", values = pal_nejm()(2)) +
            # ylim(NA, max(dfmet$met_y)*2) +
            labs(title=metname, y="Metabolite (z-score)", x = "") +
            facet_wrap(~Intervention) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    ggsave(str_c("results/metabolomics/boxplots/", metname, ".pdf"), width = 4, height = 5.5, device = "pdf")
    ggsave(str_c("results/metabolomics/boxplots/", metname, ".svg"), width = 4, height = 5.5, device = "svg")
    res_box[[a]] <- pl
}

pdf("results/metabolomics/boxplots/boxplots_met.pdf", width = 15, height = 22)
gridExtra::grid.arrange(grobs=res_box, ncol=5)
dev.off()

svg("results/metabolomics/boxplots/boxplots_met.svg", width = 15, height = 22)
gridExtra::grid.arrange(grobs=res_box, ncol=5)
dev.off()
