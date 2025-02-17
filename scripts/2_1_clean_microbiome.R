# Clean microbiome data ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)
library(rio)
library(forcats)

# Metadata table
meta <- rio::import("data/sample_info_valeriia.csv") %>% 
                    filter(str_detect(ID, "L")) %>% 
                    mutate(GenotypePerSex = str_c(Sex, " ", Genotype),
                           GenotypePerSex = fct_relevel(GenotypePerSex, "Female WT", after = 0L),
                           GenotypePerSex = fct_relevel(GenotypePerSex, "Male WT", after = 1L),
                            Age_weeks = fct_relevel(Age_weeks, "8 weeks", after = 0L),
                            Age_weeks = fct_relevel(Age_weeks, "6 weeks", after = 0L))
saveRDS(meta, "data/meta_microbiome.RDS")
write.csv(meta, "data/meta_microbiome.csv")

# Tidy bracken table
tab <- rio::import("data/combined_brackenoutput.txt")
dim(tab)
tab <- tab %>% dplyr::select(name, taxonomy_lvl, contains("_frac")) %>% 
                    filter(taxonomy_lvl == "S") %>% 
                    dplyr::select(-taxonomy_lvl) %>% 
                    rename_with(~ str_remove(., ".tsv_frac"), 2:ncol(.))
names(tab)
rownames(tab) <- tab$name
tab$name <- NULL
df <- as.data.frame(t(as.matrix(tab)))
head(df)[1:5,1:5]
df <- df %>% filter(rownames(.) %in% meta$ID)
dim(df)

threshold <- 0.0005
min_samples <- 0.1 * nrow(df)
species_filter <- apply(df, 2, function(x) sum(x >= threshold) >= min_samples)
df_filtered <- df[, species_filter]
dim(df_filtered)

saveRDS(df, "data/microbiome.RDS")
write.csv(df, "data/microbiome.csv")
saveRDS(df_filtered, "data/microbiome_filtered.RDS")
write.csv(df_filtered, "data/microbiome_filtered.csv")
