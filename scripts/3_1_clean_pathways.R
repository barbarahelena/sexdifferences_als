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

# Tidy stratified pathway table
# pathw <- rio::import("data/pathway_abundance_cpm_stratified.txt") %>%
#                 filter(!(str_detect(`# Pathway`, "UNMAPPED")) & !(str_detect(`# Pathway`, "UNINTEGRATED"))) %>%
#                 select(contains("L"), `# Pathway`) %>%
#                 rename(pathway = `# Pathway`) %>%
#                 rename_all(~ str_remove(., "_Abundance"))
# head(pathw)
# names(pathw)
# load <- read.csv("results/microbiome/tcam/pythonoutput/df_loadings.csv")
# f1 <- load %>% select(X, F1.19.99.) %>% arrange(-F1.19.99.) %>% slice(1:10, min10:last)
# f1$X <- str_replace_all(f1$X, " ", "_")
# res <- c()
# for (i in f1$X){
#     pathfiltered <- pathw %>% filter(str_detect(pathway, i))
#     res <- rbind(res, pathfiltered)
# }
# res <- as.data.frame(res)
# head(res)
# dim(res) # Ok so this only works for Akkermansia.. skip

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
summary(zero_counts > (0.5*140))
pw <- pw[,zero_counts<(0.5*140)] # 105 pathways removed

# Log-transform the pathways
pw <- log10(pw + 1)
dim(pw) # 108 pathways
head(pw)

saveRDS(pw, "data/pathways.RDS")
write.csv(pw, "data/pathways.csv")
saveRDS(lib, "data/pathwaykeys.RDS")
write.csv(lib, "data/pathwaykeys.csv", row.names = FALSE)
