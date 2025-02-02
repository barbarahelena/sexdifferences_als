# Clean microbiome data ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)
library(rio)
library(forcats)

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
saveRDS(df, "data/microbiome.RDS")

# Tidy pathway table
pathw <- rio::import("data/pathway_abundance_cpm_unstratified.txt") %>% 
                filter(`# Pathway` != "UNMAPPED" & `# Pathway` != "UNINTEGRATED")
head(pathw)
names(pathw)
paths <- pathw$`# Pathway`
keys <- str_extract(paths, "[A-Z0-9-]+(?=:)")
expl <- trimws(str_extract(paths, "(?<=:).*"))
lib <- data.frame(paths, keys, expl)
pathw <- pathw %>% rename_with(~ str_remove(., "_Abundance"), 2:ncol(.)) %>% 
                    mutate(`# Pathway` = str_extract(`# Pathway`, "[A-Z0-9-]+(?=:)"))
rownames(pathw) <- pathw$`# Pathway`
pathw$`# Pathway` <- NULL
pathdf <- as.data.frame(t(as.matrix(pathw)))
head(pathdf)[1:5,1:5]
dim(pathdf)
variances <- sapply(pathdf, sd)
summary(variances < 5)
means <- sapply(pathdf, mean)
summary(means < 40)
medians <- sapply(pathdf, median)
summary(medians < 40)
saveRDS(pathdf, "data/pathways.RDS")