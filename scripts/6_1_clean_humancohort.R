# Cleaning human cohort data
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

# Libraries
library(tidyverse)
library(ggpubr)

# Data
meta <- rio::import("data/human_cohort/human_cohort_metadata.xlsx") %>%
    mutate(Sex = case_when(Sex == "F" ~ "Female", Sex == "M" ~ "Male"),
           across(c("Height", "Weight", "Age", "ALS"), as.numeric),
           BMI = Weight / (Height/100)^2, # BMI in kg/m^2
           )
saveRDS(meta, "data/human_cohort/metadata.RDS")

# Metabolomics
met <- rio::import("data/human_cohort/cleaned_metabolomics.xlsx")
names(met)
dim(met)
metamet <- rio::import("data/human_cohort/metadata_metabolites.xlsx", header = FALSE)
metamet <- t(as.matrix(metamet))
colnames(metamet) <- metamet[1,]
metamet <- metamet[-1,]
head(metamet)
metamet <- as.data.frame(metamet)
metamet$GROUP <- case_when(
    metamet$GROUP == "Control" ~ "C",
    metamet$GROUP == "AL" ~ "ALS"
)
metamet$ID <- str_c(metamet$GROUP, "_", metamet$SUBJECT_ID) # same format as metadata
head(metamet); tail(metamet)
head(met)
rownames(met) <- met$BIOCHEMICAL
met$BIOCHEMICAL <- NULL
head(met)
colnames(met) <- metamet$ID[which(metamet$SAMPLE_ID %in% colnames(met))]
head(met)
met <- as.data.frame(t(as.matrix(met)))
head(met)
filtr <- sapply(met, function(x) sd(x, na.rm = TRUE) == 0) 
met <- met[, !filtr]
dim(met)
rownames(met)
saveRDS(met, "data/human_cohort/metabolomics_human_cohort.RDS")

# Microbiome
tab <- rio::import("data/human_cohort/combined_brackenoutput_species.txt")
dim(tab) 
names(tab) 
tab <- tab %>% dplyr::select(name, taxonomy_lvl, contains("_frac")) %>% 
                    dplyr::select(-taxonomy_lvl) %>% 
                    mutate(across(contains("_frac"), ~ .x * 100)) %>% # make %
                    rename_with(~ str_remove(., ".tsv_frac"), 2:ncol(.)) %>% 
                    rename_with(~ str_remove(., "_species"), 2:ncol(.)) 
dim(tab) # 67 samples

readcounts <- rio::import("data/human_cohort/readcounts_table.csv")
gghistogram(log10(readcounts$ReadCount)) 
min(readcounts$ReadCount) # 8.6 million reads - none below 7.5 million
rownames(tab) <- tab$name
tab$name <- NULL
df <- as.data.frame(t(as.matrix(tab)))
head(df)[1:5,1:5]
summary(meta$ID %in% rownames(df))
meta$ID[which(!meta$ID %in% rownames(df))] # 9 samples missing
rownames(df)[which(!rownames(df) %in% meta$ID)] # no samples missing
dim(df)
summary(as.factor(meta$Sex)) # missing sex for 5 samples
table(meta$Sex, meta$Group) # in ALS group 9 vs 30
df <- df[rownames(df) %in% meta$ID[which(!is.na(meta$Sex))],]
dim(df) # 61 samples

# Pruning of species
threshold <- 0.05 # 0.05% threshold
min_samples <- 0.1 * nrow(df) # in 10% of samples
species_filter <- apply(df, 2, function(x) sum(x >= threshold) >= min_samples)
df_filtered <- df[, species_filter]
dim(df_filtered) # 340 species and 57 samples

saveRDS(df_filtered, "data/human_cohort/microbiome_pruned.RDS")
saveRDS(df, "data/human_cohort/microbiome.RDS")
