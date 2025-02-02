## Welch t tests for metabolomics sex differences
## Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Libraries
library(tidyverse)

## Data
als <- readRDS("data/metadata.RDS")
met <- readRDS("data/metabolomics.RDS")
met$ID <- rownames(met)
met <- met %>% # needs to be long for join with metadata
    pivot_longer(., cols = 1:ncol(.)-1, names_to = "metabolite", values_to = "value")
dftot <- full_join(als, met, by = "ID") # join
dftdp <- dftot %>% # filter and pivot back so that metabolites are vars
    filter(Intervention == "TDP43") %>%
    select(ID, Sex, metabolite, value) %>%
    pivot_wider(names_from = "metabolite", values_from = "value")
head(dftdp)

## Welch t tests between male and female TDP43 mice
res <- c()
for(a in 3:ncol(dftdp)){
    dftdp$met <- dftdp[[a]]
    welch.res <- t.test(dftdp$met ~ dftdp$Sex, var.equal = FALSE) # var unequal = Welch
    res.row <- c(colnames(dftdp)[a], welch.res$p.value)
    names(res.row) <- c("metabolite", "p.value")
    res <- rbind(res, res.row)
    dftdp$met <- NULL
}

res <- as.data.frame(res)
rownames(res) <- res$metabolite
res <- res %>% mutate(metabolite = as.factor(metabolite), p.value = as.numeric(p.value))
res$q.value <- p.adjust(res$p.value, method = "fdr")
nrow(res %>% filter(p.value < 0.05))
res %>% filter(p.value < 0.05)
nrow(res %>% filter(q.value < 0.05))
res %>% filter(q.value < 0.05)

res <- res %>% mutate(sig = case_when(
    p.value >= 0.05 ~ paste0(""),
    p.value < 0.0001 ~ paste0("****"),
    p.value < 0.001 ~ paste0("***"),
    p.value < 0.01 ~ paste0("**"),
    p.value < 0.05 ~ paste0("*")
), sigq = case_when(
    q.value >= 0.05 ~ paste0(""),
    q.value < 0.0001 ~ paste0("****"),
    q.value < 0.001 ~ paste0("***"),
    q.value < 0.01 ~ paste0("**"),
    q.value < 0.05 ~ paste0("*")
))

df.sum <- dftdp %>% group_by(Sex) %>% summarise(across(2:(ncol(.)-1), 
                                                    list(mean = ~mean(.x, na.rm = TRUE),
                                                         sd = ~sd(.x, na.rm = TRUE),
                                                         n = ~length(.x)), 
                                                    .names = "{.col}__{.fn}"))

df.sumlong <- df.sum %>% pivot_longer(2:ncol(.), names_to = c("metabolite", "stat"),
                                      names_sep = "__",
                                      values_to = "number") %>%
    pivot_wider(id_cols = 1:2, values_from = "number", names_from = "stat")

restot <- left_join(res, df.sumlong, by = "metabolite") 
ressel <- restot %>% filter(p.value < 0.05) %>% arrange(p.value)
pval_sig <- unique(ressel$metabolite)
pval_sig

# Save tables
write_csv(restot, "r_results/metabolites_welcht_diff.csv")
write_csv(ressel, "r_results/metabolites_welcht_sig.csv")

## Welch t tests between male and female control mice
dfctrl <- dftot %>% # filter and pivot back so that metabolites are vars
    filter(Intervention == "Control") %>%
    select(ID, Sex, metabolite, value) %>%
    pivot_wider(names_from = "metabolite", values_from = "value")
head(dfctrl)

res <- c()
for(a in 3:ncol(dfctrl)){
    dfctrl$met <- dfctrl[[a]]
    welch.res <- t.test(dfctrl$met ~ dfctrl$Sex, var.equal = FALSE)
    res.row <- c(colnames(dfctrl)[a], welch.res$p.value)
    names(res.row) <- c("metabolite", "p.value")
    res <- rbind(res, res.row)
    dfctrl$met <- NULL
}

res <- as.data.frame(res)
rownames(res) <- res$metabolite
res <- res %>% mutate(metabolite = as.factor(metabolite), p.value = as.numeric(p.value))
res$q.value <- p.adjust(res$p.value, method = "fdr")
nrow(res %>% filter(p.value < 0.05))
res %>% filter(p.value < 0.05)
nrow(res %>% filter(q.value < 0.05))
res %>% filter(q.value < 0.05)

res <- res %>% mutate(sig = case_when(
    p.value >= 0.05 ~ paste0(""),
    p.value < 0.0001 ~ paste0("****"),
    p.value < 0.001 ~ paste0("***"),
    p.value < 0.01 ~ paste0("**"),
    p.value < 0.05 ~ paste0("*")
), sigq = case_when(
    q.value >= 0.05 ~ paste0(""),
    q.value < 0.0001 ~ paste0("****"),
    q.value < 0.001 ~ paste0("***"),
    q.value < 0.01 ~ paste0("**"),
    q.value < 0.05 ~ paste0("*")
))

df.sum <- dfctrl %>% group_by(Sex) %>% summarise(across(2:(ncol(.)-1), 
                                                    list(mean = ~mean(.x, na.rm = TRUE),
                                                         sd = ~sd(.x, na.rm = TRUE),
                                                         n = ~length(.x)), 
                                                    .names = "{.col}__{.fn}"))

df.sumlong <- df.sum %>% pivot_longer(2:ncol(.), names_to = c("metabolite", "stat"),
                                      names_sep = "__",
                                      values_to = "number") %>%
    pivot_wider(id_cols = 1:2, values_from = "number", names_from = "stat")

restot <- left_join(res, df.sumlong, by = "metabolite") 
ressel <- restot %>% filter(p.value < 0.05) %>% arrange(p.value)
pval_sig <- unique(ressel$metabolite)
pval_sig

# Save tables
write_csv(restot, "r_results/metabolites_welcht_ctrl_diff.csv")
write_csv(ressel, "r_results/metabolites_welcht_ctrl_sig.csv")

## Welch t tests between control and TDP43 mice
dftot2 <- dftot %>% # filter and pivot back so that metabolites are vars
    select(ID, Intervention, metabolite, value) %>%
    pivot_wider(names_from = "metabolite", values_from = "value")
head(dftot2)

res <- c()
for(a in 3:ncol(dftot2)){
    dftot2$met <- dftot2[[a]]
    welch.res <- t.test(dftot2$met ~ dftot2$Intervention, var.equal = FALSE)
    res.row <- c(colnames(dftot2)[a], welch.res$p.value)
    names(res.row) <- c("metabolite", "p.value")
    res <- rbind(res, res.row)
    dftot2$met <- NULL
}

res <- as.data.frame(res)
rownames(res) <- res$metabolite
res <- res %>% mutate(metabolite = as.factor(metabolite), p.value = as.numeric(p.value))
res$q.value <- p.adjust(res$p.value, method = "fdr")
nrow(res %>% filter(p.value < 0.05))
res %>% filter(p.value < 0.05)
nrow(res %>% filter(q.value < 0.05))
res %>% filter(q.value < 0.05)

res <- res %>% mutate(sig = case_when(
    p.value >= 0.05 ~ paste0(""),
    p.value < 0.0001 ~ paste0("****"),
    p.value < 0.001 ~ paste0("***"),
    p.value < 0.01 ~ paste0("**"),
    p.value < 0.05 ~ paste0("*")
), sigq = case_when(
    q.value >= 0.05 ~ paste0(""),
    q.value < 0.0001 ~ paste0("****"),
    q.value < 0.001 ~ paste0("***"),
    q.value < 0.01 ~ paste0("**"),
    q.value < 0.05 ~ paste0("*")
))

df.sum <- dftot2 %>% group_by(Intervention) %>% summarise(across(2:(ncol(.)-1), 
                                                    list(mean = ~mean(.x, na.rm = TRUE),
                                                         sd = ~sd(.x, na.rm = TRUE),
                                                         n = ~length(.x)), 
                                                    .names = "{.col}__{.fn}"))

df.sumlong <- df.sum %>% pivot_longer(2:ncol(.), names_to = c("metabolite", "stat"),
                                      names_sep = "__",
                                      values_to = "number") %>%
    pivot_wider(id_cols = 1:2, values_from = "number", names_from = "stat")

restot <- left_join(res, df.sumlong, by = "metabolite") 
ressel <- restot %>% filter(p.value < 0.05) %>% arrange(p.value)
pval_sig <- unique(ressel$metabolite)
pval_sig

# Save tables
write_csv(restot, "r_results/metabolites_welcht_mice_diff.csv")
write_csv(ressel, "r_results/metabolites_welcht_mice_sig.csv")

## Welch t tests between control and TDP43 mice - females
dftot4 <- dftot %>% # filter and pivot back so that metabolites are vars
    filter(Sex == "Female") %>%
    select(ID, Intervention, metabolite, value) %>%
    pivot_wider(names_from = "metabolite", values_from = "value")
head(dftot4)

res <- c()
for(a in 3:ncol(dftot4)){
    dftot4$met <- dftot4[[a]]
    welch.res <- t.test(dftot4$met ~ dftot4$Intervention, var.equal = FALSE)
    res.row <- c(colnames(dftot4)[a], welch.res$p.value)
    names(res.row) <- c("metabolite", "p.value")
    res <- rbind(res, res.row)
    dftot4$met <- NULL
}

res <- as.data.frame(res)
rownames(res) <- res$metabolite
res <- res %>% mutate(metabolite = as.factor(metabolite), p.value = as.numeric(p.value))
res$q.value <- p.adjust(res$p.value, method = "fdr")
nrow(res %>% filter(p.value < 0.05))
res %>% filter(p.value < 0.05)
nrow(res %>% filter(q.value < 0.05))
res %>% filter(q.value < 0.05)

res <- res %>% mutate(sig = case_when(
    p.value >= 0.05 ~ paste0(""),
    p.value < 0.0001 ~ paste0("****"),
    p.value < 0.001 ~ paste0("***"),
    p.value < 0.01 ~ paste0("**"),
    p.value < 0.05 ~ paste0("*")
), sigq = case_when(
    q.value >= 0.05 ~ paste0(""),
    q.value < 0.0001 ~ paste0("****"),
    q.value < 0.001 ~ paste0("***"),
    q.value < 0.01 ~ paste0("**"),
    q.value < 0.05 ~ paste0("*")
))

df.sum <- dftot2 %>% group_by(Intervention) %>% summarise(across(2:(ncol(.)-1), 
                                                    list(mean = ~mean(.x, na.rm = TRUE),
                                                         sd = ~sd(.x, na.rm = TRUE),
                                                         n = ~length(.x)), 
                                                    .names = "{.col}__{.fn}"))

df.sumlong <- df.sum %>% pivot_longer(2:ncol(.), names_to = c("metabolite", "stat"),
                                      names_sep = "__",
                                      values_to = "number") %>%
    pivot_wider(id_cols = 1:2, values_from = "number", names_from = "stat")

restot <- left_join(res, df.sumlong, by = "metabolite") 
ressel <- restot %>% filter(p.value < 0.05) %>% arrange(p.value)
pval_sig <- unique(ressel$metabolite)
pval_sig

# Save tables
write_csv(restot, "r_results/metabolites_welcht_mice_fem_diff.csv")
write_csv(ressel, "r_results/metabolites_welcht_mice_fem_sig.csv")

## Welch t tests between control and TDP43 mice - males
dftot3 <- dftot %>% # filter and pivot back so that metabolites are vars
    filter(Sex == "Male") %>%
    select(ID, Intervention, metabolite, value) %>%
    pivot_wider(names_from = "metabolite", values_from = "value")
head(dftot3)

res <- c()
for(a in 3:ncol(dftot3)){
    dftot3$met <- dftot3[[a]]
    welch.res <- t.test(dftot3$met ~ dftot3$Intervention, var.equal = FALSE)
    res.row <- c(colnames(dftot3)[a], welch.res$p.value)
    names(res.row) <- c("metabolite", "p.value")
    res <- rbind(res, res.row)
    dftot3$met <- NULL
}

res <- as.data.frame(res)
rownames(res) <- res$metabolite
res <- res %>% mutate(metabolite = as.factor(metabolite), p.value = as.numeric(p.value))
res$q.value <- p.adjust(res$p.value, method = "fdr")
nrow(res %>% filter(p.value < 0.05))
res %>% filter(p.value < 0.05)
nrow(res %>% filter(q.value < 0.05))
res %>% filter(q.value < 0.05)

res <- res %>% mutate(sig = case_when(
    p.value >= 0.05 ~ paste0(""),
    p.value < 0.0001 ~ paste0("****"),
    p.value < 0.001 ~ paste0("***"),
    p.value < 0.01 ~ paste0("**"),
    p.value < 0.05 ~ paste0("*")
), sigq = case_when(
    q.value >= 0.05 ~ paste0(""),
    q.value < 0.0001 ~ paste0("****"),
    q.value < 0.001 ~ paste0("***"),
    q.value < 0.01 ~ paste0("**"),
    q.value < 0.05 ~ paste0("*")
))

df.sum <- dftot2 %>% group_by(Intervention) %>% summarise(across(2:(ncol(.)-1), 
                                                    list(mean = ~mean(.x, na.rm = TRUE),
                                                         sd = ~sd(.x, na.rm = TRUE),
                                                         n = ~length(.x)), 
                                                    .names = "{.col}__{.fn}"))

df.sumlong <- df.sum %>% pivot_longer(2:ncol(.), names_to = c("metabolite", "stat"),
                                      names_sep = "__",
                                      values_to = "number") %>%
    pivot_wider(id_cols = 1:2, values_from = "number", names_from = "stat")

restot <- left_join(res, df.sumlong, by = "metabolite") 
ressel <- restot %>% filter(p.value < 0.05) %>% arrange(p.value)
pval_sig <- unique(ressel$metabolite)
pval_sig

# Save tables
write_csv(restot, "r_results/metabolites_welcht_mice_male_diff.csv")
write_csv(ressel, "r_results/metabolites_welcht_mice_male_sig.csv")
