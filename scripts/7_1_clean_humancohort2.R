# Cleaning human cohort data - Calgary cohort
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggpubr)

# Data
meta <- rio::import("data/human_cohort2/metadata.xlsx") |> 
          mutate(
                  #ID = str_c("X", fecal_ID),
                  Group = case_when(
                      diagnosis == "control" ~ "Control",
                      diagnosis == "als" ~ "ALS"
                  ),
                  Group = as.factor(Group),
                  sex = fct_recode(sex, "Male" = "male", "Female" = "female")
                ) |> 
          select(ID = fecal_ID, Household = household_ID, Group, Age = age_yrs,
                  Sex = sex, Duration = als_symptom_duration_months, 
                  Familial = sporadic_or_genetic, Mutation = genetic_mutation,
                  ALSFRS = alsfrs)
saveRDS(meta, "data/human_cohort2/metadata.RDS")

# Microbiome
tab <- rio::import("data/human_cohort2/combined_brackenoutput.txt")
dim(tab) 
names(tab) 
tab <- tab %>% dplyr::select(name, taxonomy_lvl, contains("_frac")) %>% 
                    dplyr::select(-taxonomy_lvl) %>% 
                    mutate(across(contains("_frac"), ~ .x * 100)) %>% # make %
                    rename_with(~ str_remove(., ".tsv_frac"), 2:ncol(.)) %>% 
                    rename_with(~ str_remove(., "_species"), 2:ncol(.)) 
dim(tab) # 39 samples

readcounts <- rio::import("data/human_cohort2/readcounts_table.csv")
gghistogram(readcounts$ReadCount)
min(readcounts$ReadCount) # 7.7 million reads
rownames(tab) <- tab$name
tab$name <- NULL
df <- as.data.frame(t(as.matrix(tab)))
head(df)[1:5,1:5]
summary(meta$ID %in% rownames(df))
meta$ID[which(!meta$ID %in% rownames(df))] # 2 samples missing
rownames(df)[which(!rownames(df) %in% meta$ID)] # no samples missing
dim(df)
summary(as.factor(meta$Sex))
table(meta$Sex, meta$Group) # in ALS group 9 vs 30
df <- df[rownames(df) %in% meta$ID[which(!is.na(meta$Sex))],]
dim(df) # 38 samples

# Pruning of species
threshold <- 0.05 # 0.05% threshold
min_samples <- 0.1 * nrow(df) # in 10% of samples
species_filter <- apply(df, 2, function(x) sum(x >= threshold) >= min_samples)
df_filtered <- df[, species_filter]
dim(df_filtered) # 210 species and 38 samples

saveRDS(df_filtered, "data/human_cohort2/microbiome_pruned.RDS")
saveRDS(df, "data/human_cohort2/microbiome.RDS")

# Pathways
pathw <- rio::import("data/human_cohort2/pathway_abundance_cpm_unstratified.txt") |> 
                select(Pathways = `# Pathway HUMAnN v4.0.0.alpha.1`, everything()) |> 
                filter(Pathways != "UNMAPPED" & Pathways != "UNINTEGRATED") %>%
                rename_with(~ str_remove(., "_Abundance"), 1:ncol(.))
dim(pathw)
head(pathw)

paths <- pathw$Pathways
keys <- str_extract(paths, "[A-Z0-9-]+(?=:)")
expl <- trimws(str_extract(paths, "(?<=:).*"))
lib <- data.frame(paths, keys, expl)
pathw <- pathw |> mutate(
    Pathways = str_extract(Pathways, "[A-Z0-9-]+(?=:)"))
rownames(pathw) <- pathw$Pathways
pathw$Pathways <- NULL
pathdf <- as.data.frame(t(as.matrix(pathw)))

# Filter
pathway_filter <- apply(pathdf, 2, function(x) sum(x > 5) >= (0.5 * nrow(pathdf)))
pw_filtered <- pathdf[, pathway_filter]

# Log transform filtered pathways
pw <- log10(pw_filtered + 1)

# Print dimensions
dim(pw)
head(pw)

saveRDS(pw, "data/human_cohort2/human_als_pathways.RDS")
write.csv(pw, "data/human_cohort2/human_als_pathways.csv")
saveRDS(lib, "data/human_cohort2/pathwaykeys.RDS")
write.csv(lib, "data/human_cohort2/pathwaykeys.csv", row.names = FALSE)

lib %>% filter(str_detect(expl, "rhamnose"))
#                                                              paths            keys                                            expl
# 1             DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis DTDPRHAMSYN-PWY             dTDP-&beta;-L-rhamnose biosynthesis
# 2 FUC-RHAMCAT-PWY: superpathway of fucose and rhamnose degradation FUC-RHAMCAT-PWY superpathway of fucose and rhamnose degradation
# 3                            RHAMCAT-PWY: L-rhamnose degradation I     RHAMCAT-PWY                        L-rhamnose degradation I