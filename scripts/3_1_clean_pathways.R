# Clean microbiome data ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)
library(rio)
library(forcats)

# Tidy pathway table
pathw <- rio::import("data/pathway_abundance_cpm_unstratified.txt") %>% 
                filter(`# Pathway` != "UNMAPPED" & `# Pathway` != "UNINTEGRATED") %>%
                select(`# Pathway`, contains("L")) %>%
                rename_with(~ str_remove(., "_Abundance"), 1:ncol(.))

head(pathw)
names(pathw)
paths <- pathw$`# Pathway`
keys <- str_extract(paths, "[A-Z0-9-]+(?=:)")
expl <- trimws(str_extract(paths, "(?<=:).*"))
lib <- data.frame(paths, keys, expl)
pathw <- pathw %>% mutate(`# Pathway` = str_extract(`# Pathway`, "[A-Z0-9-]+(?=:)"))
rownames(pathw) <- pathw$`# Pathway`
pathw$`# Pathway` <- NULL
pathdf <- as.data.frame(t(as.matrix(pathw)))

# Filter
pathway_filter <- apply(pathdf, 2, function(x) sum(x > 5) >= (0.5 * nrow(pathdf)))
pw_filtered <- pathdf[, pathway_filter]

# Log transform filtered pathways
pw <- log10(pw_filtered + 1)

# Print dimensions
dim(pw)
head(pw)

saveRDS(pw, "data/pathways.RDS")
write.csv(pw, "data/pathways.csv")
saveRDS(lib, "data/pathwaykeys.RDS")
write.csv(lib, "data/pathwaykeys.csv", row.names = FALSE)
