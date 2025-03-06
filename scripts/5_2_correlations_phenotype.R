# Weights
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

# Libraries
library(readxl)
library(tidyverse)
library(lubridate)

# Data
weight <- readRDS("data/weights_metadata.RDS")
onset <- readRDS("data/onset_metadata.RDS")
microbiome <- readRDS("data/microbiome_run1.RDS")

akk <- microbiome %>% select(`Akkermansia muciniphila`)
head(akk)
akk$ID <- rownames(akk)
head(weight)
tot <- left_join(weight, akk) %>% left_join(., onset)
head(tot)

akk6 <- tot %>% filter(Age_ints == "10", ) %>% select(MouseID, Akk_baseline = `Akkermansia muciniphila`)
tot2 <- tot %>% left_join(., akk6, relationship = "many-to-one") %>% filter(Age_ints == 12)

ggplot(data = tot2, aes(x = log10(Akk_baseline + 1), y = weight, color = Sex)) + 
  theme_minimal() + 
  ggsci::scale_color_nejm() +
  geom_jitter() +
  geom_smooth(method = "lm") + 
  ggpubr::stat_cor(method = "spearman")

ggplot(data = tot2, aes(x = log10(Akk_baseline + 1), y = Age_onset_weeks, color = Sex)) + 
  theme_minimal() + 
  ggsci::scale_color_nejm() +
  geom_jitter() +
  geom_smooth(method = "lm") + 
  ggpubr::stat_cor(method = "spearman")

por <- microbiome %>% select(`Porphyromonadaceae bacterium UBA7141`)
head(por)
por$ID <- rownames(por)
head(weight)
tot <- left_join(weight, por) %>% filter(!is.na(`Porphyromonadaceae bacterium UBA7141`))
head(tot)
ggplot(data = tot %>% filter(Age_ints == "6"), 
          aes(x = log10(`Porphyromonadaceae bacterium UBA7141` + 1), y = weight, color = Sex)) + 
  theme_minimal() + 
  ggsci::scale_color_nejm() +
  geom_jitter() +
  geom_smooth(method = "lm") + 
  ggpubr::stat_cor(method = "spearman")
