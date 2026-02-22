# all microbes
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
                                          size = rel(0.8), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
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
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
} 

# Data
mb <- readRDS("data/human_cohort/microbiome_pruned.RDS")
df <- readRDS("data/human_cohort/metadata.RDS")


dim(mb)
mb <- apply(mb, 2, function(x) log10(x + 0.01))
mean(mb[,3])
mb <- as.data.frame(mb)
mb$ID <- rownames(mb)

# Metadata
df_tot <- left_join(mb, df, by = c("ID"))

statres <- c()
for(i in c(1:(ncol(mb)-1))) {
    df_tot$microbe <- df_tot[,i]
    mbname <- colnames(df_tot)[i]
    test <- wilcox.test(microbe ~ Group, data = df_tot)
    pval_diag <- as.numeric(test$p.value)
    pval_diag <- as.numeric(format(round(pval_diag, 3), nsmall = 3))
    sig_diag <- case_when(
        pval_diag < 0.0001 ~ paste0("****"),
        pval_diag < 0.001 ~paste0("***"),
        pval_diag < 0.01 ~paste0("**"),
        pval_diag <= 0.05 ~paste0("*"),
        pval_diag > 0.05 ~paste0("")
    )
    statres_line <- cbind(mbname, pval_diag, sig_diag)
    statres <- rbind(statres, statres_line)
}
statres <- as.data.frame(statres)
statres <- statres %>% arrange(pval_diag) %>% 
    mutate(qval_diag = p.adjust(pval_diag, method = "fdr")) |> 
    mutate(sigq_diag = case_when(
        qval_diag < 0.0001 ~ paste0("****"),
        qval_diag < 0.001 ~paste0("***"),
        qval_diag < 0.01 ~paste0("**"),
        qval_diag <= 0.05 ~paste0("*"),
        qval_diag > 0.05 ~paste0("")
    ))
names(statres)

statresfem <- c()
for(i in c(1:(ncol(mb)-1))) {
    df_tot2 <- df_tot |> filter(Sex == "Female")
    df_tot2$microbe <- df_tot2[,i]
    mbname <- colnames(df_tot2)[i]
    test <- wilcox.test(microbe ~ Group, data = df_tot2)
    pval_fem <- as.numeric(test$p.value)
    pval_fem <- as.numeric(format(round(pval_fem, 3), nsmall = 3))
    sig_fem <- case_when(
        pval_fem < 0.0001 ~ paste0("****"),
        pval_fem < 0.001 ~paste0("***"),
        pval_fem < 0.01 ~paste0("**"),
        pval_fem <= 0.05 ~paste0("*"),
        pval_fem > 0.05 ~paste0("")
    )
    statres_line <- cbind(mbname, pval_fem, sig_fem)
    statresfem <- rbind(statresfem, statres_line)
}
statresfem <- as.data.frame(statresfem)
statresfem <- statresfem %>% arrange(pval_fem) %>% 
    mutate(qval_fem = p.adjust(pval_fem, method = "fdr")) |> 
    mutate(sigq_fem = case_when(
        qval_fem < 0.0001 ~ paste0("****"),
        qval_fem < 0.001 ~paste0("***"),
        qval_fem < 0.01 ~paste0("**"),
        qval_fem <= 0.05 ~paste0("*"),
        qval_fem > 0.05 ~paste0("")
    ))

statresmale <- c()
for(i in c(1:(ncol(mb)-1))) {
    df_tot2 <- df_tot |> filter(Sex == "Male")
    df_tot2$microbe <- df_tot2[,i]
    mbname <- colnames(df_tot2)[i]
    test <- wilcox.test(microbe ~ Group, data = df_tot2)
    pval_male <- as.numeric(test$p.value)
    pval_male <- as.numeric(format(round(pval_male, 3), nsmall = 3))
    sig_male <- case_when(
        pval_male < 0.0001 ~ paste0("****"),
        pval_male < 0.001 ~paste0("***"),
        pval_male < 0.01 ~paste0("**"),
        pval_male <= 0.05 ~paste0("*"),
        pval_male > 0.05 ~paste0("")
    )
    statres_line <- cbind(mbname, pval_male, sig_male)
    statresmale <- rbind(statresmale, statres_line)
}
statresmale <- as.data.frame(statresmale)
statresmale <- statresmale %>% arrange(pval_male) %>% 
    mutate(qval_male = p.adjust(pval_male, method = "fdr")) |> 
    mutate(sigq_fem = case_when(
        qval_male < 0.0001 ~ paste0("****"),
        qval_male < 0.001 ~paste0("***"),
        qval_male < 0.01 ~paste0("**"),
        qval_male <= 0.05 ~paste0("*"),
        qval_male > 0.05 ~paste0("")
    ))

statresals <- c()
for(i in c(1:(ncol(mb)-1))) {
    df_tot3 <- df_tot |> filter(Group == "ALS")
    df_tot3$microbe <- df_tot3[,i]
    mbname <- colnames(df_tot3)[i]
    test <- wilcox.test(microbe ~ Sex, data = df_tot3)
    pval_als <- as.numeric(test$p.value)
    pval_als <- as.numeric(format(round(pval_als, 3), nsmall = 3))
    sig_als <- case_when(
        pval_als < 0.0001 ~ paste0("****"),
        pval_als < 0.001 ~paste0("***"),
        pval_als < 0.01 ~paste0("**"),
        pval_als <= 0.05 ~paste0("*"),
        pval_als > 0.05 ~paste0("")
    )
    statres_line <- cbind(mbname, pval_als, sig_als)
    statresals <- rbind(statresals, statres_line)
}
statresals <- as.data.frame(statresals)
statresals <- statresals %>% arrange(pval_als) %>%
    mutate(qval_als = p.adjust(pval_als, method = "fdr")) |>
    mutate(sigq_als = case_when(
        qval_als < 0.0001 ~ paste0("****"),
        qval_als < 0.001 ~paste0("***"),
        qval_als < 0.01 ~paste0("**"),
        qval_als <= 0.05 ~paste0("*"),
        qval_als > 0.05 ~paste0("")
    ))

statresctrl <- c()
for(i in c(1:(ncol(mb)-1))) {
    df_tot3 <- df_tot |> filter(Group == "Healthy")
    df_tot3$microbe <- df_tot3[,i]
    mbname <- colnames(df_tot3)[i]
    test <- wilcox.test(microbe ~ Sex, data = df_tot3)
    pval_ctrl <- as.numeric(test$p.value)
    pval_ctrl <- as.numeric(format(round(pval_ctrl, 3), nsmall = 3))
    sig_ctrl <- case_when(
        pval_ctrl < 0.0001 ~ paste0("****"),
        pval_ctrl < 0.001 ~paste0("***"),
        pval_ctrl < 0.01 ~paste0("**"),
        pval_ctrl <= 0.05 ~paste0("*"),
        pval_ctrl > 0.05 ~paste0("")
    )
    statres_line <- cbind(mbname, pval_ctrl, sig_ctrl)
    statresctrl <- rbind(statresctrl, statres_line)
}
statresctrl <- as.data.frame(statresctrl)
statresctrl <- statresctrl %>% arrange(pval_ctrl) %>%
    mutate(qval_ctrl = p.adjust(pval_ctrl, method = "fdr")) |>
    mutate(sigq_ctrl = case_when(
        qval_ctrl < 0.0001 ~ paste0("****"),
        qval_ctrl < 0.001 ~paste0("***"),
        qval_ctrl < 0.01 ~paste0("**"),
        qval_ctrl <= 0.05 ~paste0("*"),
        qval_ctrl > 0.05 ~paste0("")
    ))

restot <- right_join(statres, statresfem) |> right_join(statresmale) |> right_join(statresals) |> right_join(statresctrl)
head(restot)
restot |> arrange(qval_diag)
restot |> arrange(qval_fem)
restot |> arrange(qval_male)
restot |> arrange(qval_als)
restot |> arrange(qval_ctrl)
write.csv(restot, "results/humancohort/wilcoxons_microbes.csv")

# Boxplots for microbes with pval_als < 0.05 (sex differences in ALS)
res_sigals <- statresals |> filter(pval_als < 0.05)
dim(res_sigals)

plist_als <- list()
for(i in 1:nrow(res_sigals)){
    nm <- res_sigals$mbname[i]
    df_tot$mb <- df_tot[, nm]

    # Test diagnosis difference in females and males
    diagdiff1 <- wilcox.test(mb ~ Group, data = df_tot %>% filter(Sex == "Female"))
    diagdiff2 <- wilcox.test(mb ~ Group, data = df_tot %>% filter(Sex == "Male"))

    # Format p-values
    p_female <- format.pval(diagdiff1$p.value, digits = 2)
    p_male <- format.pval(diagdiff2$p.value, digits = 2)

    # Create subtitle with diagnosis difference p-values
    diag_diff_text <- paste0("ALS-Control p = ", p_female, " (women); p = ", p_male, " (men)")

    pl_als <- ggplot(df_tot, aes(x = Sex, y = mb)) +
        geom_boxplot(aes(fill = Sex), outlier.shape = NA, width = 0.4) +
        geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        labs(title = nm,
             caption = diag_diff_text,
             y = "log10(abundance + 0.01)",
             x = "") +
        stat_compare_means(method = "wilcox.test", label = "p.format") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
              plot.caption = element_text(size = 12),
              legend.position = "none") +
        facet_wrap(~ Group)

    plist_als[[i]] <- pl_als
}

n_plots <- length(plist_als)
n_cols <- ceiling(sqrt(n_plots))
n_rows <- ceiling(n_plots / n_cols)

(plots_als <- ggarrange(plotlist = plist_als, common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:n_plots],
          nrow = n_rows, ncol = n_cols))
ggsave(plots_als, filename = "results/humancohort/sexdiff_als_boxplots.pdf",
       width = 4 * n_cols, height = 5 * n_rows)
