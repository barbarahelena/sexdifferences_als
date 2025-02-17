# Clean microbiome data ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)
library(rio)
library(forcats)

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

# Explore
head(pathdf)[1:5,1:5]
dim(pathdf) # 234 pathways
variances <- sapply(pathdf, sd)
summary(variances < 5)
pw <- pathdf[,variances > 5] # 21 pathways removed
means <- sapply(pathdf, mean)
summary(means < 0.1)
medians <- sapply(pathdf, median)
summary(medians < 40)
zero_counts <- apply(pw, 2, function(x) sum(x == 0))
summary(zero_counts > (0.8*140))
pw <- pw[,zero_counts<(0.8*140)] # 59 pathways removed

# Log-transform the pathways
pw <- log10(pw + 1)
dim(pw)
head(pw)

saveRDS(pw, "data/pathways.RDS")
write.csv(pw, "data/pathways.csv")
saveRDS(lib, "data/pathwaykeys.RDS")
write.csv(lib, "data/pathwaykeys.csv", row.names = FALSE)
