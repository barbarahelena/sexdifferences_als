# Clean metabolomics data ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)
library(rio)
library(forcats)

## Open metadata
df <- import("data/metadata.xlsx") %>% select(ID = Sample, Group) %>%
    mutate(Intervention = case_when(str_detect(Group, "Disease") ~ "TDP43",
                           str_detect(Group, "Control") ~ "Control",
                           .default = NA),
            Sex = case_when(str_detect(Group, "Female") ~ "Female",
                            str_detect(Group, "Male") ~ "Male",
                            .default = NA),
            GroupPerSex = str_c(Sex, " ", Intervention),
            across(c(Sex, Intervention, GroupPerSex), as.factor),
            Sex = fct_rev(Sex))
summary(df$Sex); summary(df$Intervention)
saveRDS(df, "data/metadata.RDS")

# Open metabolomics data
met <- import("data/metabolomics.xlsx") %>% 
    select(!contains("Std") & !contains("Blank")) %>%
    mutate(compoundId = str_to_lower(compoundId))
IDs <- colnames(met)[2:ncol(met)]
rownames(met) <- met$compoundId
met$compoundId <- NULL
met <- apply(t(met), 2, function(x) log2(x + 0.01))
rownames(met) <- IDs
met <- as.data.frame(met)
head(met)
dim(met)
any(sapply(met, function(x) sum(is.na(x))))
sapply(met, function(x) sd(x, na.rm = TRUE))
saveRDS(met, "data/metabolomics.RDS")
