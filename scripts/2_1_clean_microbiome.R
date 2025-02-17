# Clean microbiome data ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)
library(rio)
library(forcats)

# Metadata table
meta <- rio::import("data/sample_info_valeriia.csv") %>% 
                    filter(str_detect(ID, "L") | str_detect(ID, "S")) %>% 
                    mutate(GenotypePerSex = str_c(Sex, " ", Genotype),
                           GenotypePerSex = fct_relevel(GenotypePerSex, "Female WT", after = 0L),
                           GenotypePerSex = fct_relevel(GenotypePerSex, "Male WT", after = 1L),
                            Age_weeks = fct_relevel(Age_weeks, "8 weeks", after = 0L),
                            Age_weeks = fct_relevel(Age_weeks, "6 weeks", after = 0L),
                            LowCount = case_when(counts < 1000000 ~ TRUE,
                                                 .default = FALSE))
lowcountids <- meta$ID[which(meta$LowCount == TRUE)]

saveRDS(meta, "data/meta_microbiome.RDS")
write.csv(meta, "data/meta_microbiome.csv")
meta2 <- meta %>% filter(str_detect(ID, "L")) # Run 1
saveRDS(meta, "data/meta_microbiome_run1.RDS")
write.csv(meta, "data/meta_microbiome_run1.csv")
meta3 <- meta %>% filter(str_detect(ID, "S")) # Run 2
saveRDS(meta, "data/meta_microbiome_run2.RDS")
write.csv(meta, "data/meta_microbiome_run2.csv")

# Tidy bracken table
tab <- rio::import("data/combined_brackenoutput.txt")
dim(tab) 
names(tab) 
tab <- tab %>% dplyr::select(name, taxonomy_lvl, contains("_frac")) %>% 
                    filter(taxonomy_lvl == "S") %>% 
                    dplyr::select(-taxonomy_lvl) %>% 
                    rename_with(~ str_remove(., ".tsv_frac"), 2:ncol(.)) 
dim(tab) # 323 samples
tab <- tab %>% dplyr::select(-all_of(lowcountids))
names(tab)
dim(tab) # 318 samples
rownames(tab) <- tab$name
tab$name <- NULL
df <- as.data.frame(t(as.matrix(tab)))
head(df)[1:5,1:5]
df2 <- df %>% filter(rownames(.) %in% meta2$ID) # Run 1
dim(df2) # 137 samples
df3 <- df %>% filter(rownames(.) %in% meta3$ID) # Run 2
dim(df3) # 144 samples

saveRDS(df2, "data/microbiome_run1.RDS")
write.csv(df2, "data/microbiome_run1.csv")
saveRDS(df3, "data/microbiome_run2.RDS")
write.csv(df3, "data/microbiome_run2.csv")

# Pruning of species #
## Run 1
threshold <- 0.0005 # 0.05% threshold
min_samples <- 0.1 * nrow(df2) # in 10% of samples
species_filter <- apply(df2, 2, function(x) sum(x >= threshold) >= min_samples)
df_filtered2 <- df2[, species_filter]
dim(df_filtered2) # 134 species

## Run 2
threshold <- 0.0005 # 0.05% threshold
min_samples <- 0.1 * nrow(df3) # in 10% of samples
species_filter <- apply(df3, 2, function(x) sum(x >= threshold) >= min_samples)
df_filtered3 <- df3[, species_filter]
dim(df_filtered3) # 169 species

saveRDS(df_filtered2, "data/microbiome_filtered_run1.RDS")
write.csv(df_filtered2, "data/microbiome_filtered_run1.csv")
saveRDS(df_filtered3, "data/microbiome_filtered_run2.RDS")
write.csv(df_filtered3, "data/microbiome_filtered_run2.csv")
