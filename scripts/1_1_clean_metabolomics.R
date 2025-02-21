# Clean metabolomics data ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)
library(rio)
library(forcats)

## Open metadata
df <- import("data/metadata.xlsx") %>% select(ID = Sample, Group) %>%
    mutate(Intervention = case_when(str_detect(Group, "Disease") ~ "TDP43",
                           str_detect(Group, "Control") ~ "WT",
                           .default = NA),
            Intervention = fct_relevel(Intervention, "WT", after = 0L),
            Sex = case_when(str_detect(Group, "Female") ~ "Female",
                            str_detect(Group, "Male") ~ "Male",
                            .default = NA),
            GroupPerSex = str_c(Sex, " ", Intervention),
            GroupPerSex = fct_relevel(GroupPerSex, "Female WT", after = 0L),
            GroupPerSex = fct_relevel(GroupPerSex, "Male WT", after = 1L),
            GroupPerSex = fct_relevel(GroupPerSex, "Female TDP43", after = 2L),
            across(c(Sex, Intervention, GroupPerSex), as.factor),
            Sex = fct_rev(Sex))
summary(df$Sex); summary(df$Intervention); summary(df$GroupPerSex)
saveRDS(df, "data/metadata.RDS")

# Open metabolomics data
met <- import("data/metabolomics.xlsx") %>% 
    select(!contains("Std") & !contains("Blank")) %>%
    mutate(compoundId = str_to_lower(compoundId))
hmdb <- met[,1:2]
hmdb <- hmdb %>% rename(metabolite = compoundId)
saveRDS(hmdb, "data/hmdb_ids.RDS")

met <- met[,c(1,3:ncol(met))]
IDs <- colnames(met)[2:ncol(met)]
rownames(met) <- met$compoundId
met$compoundId <- NULL

met2 <- t(met)
rownames(met2) <- IDs
met2 <- as.data.frame(met2)
colnames(met2) <- hmdb$HMDB[which(hmdb$metabolite == colnames(met2))]
met2 <- met2[,1:ncol(met2)-1]
met2$group <- meta$Intervention[which(meta$ID == rownames(met2))]
met2 <- met2[,c("group", colnames(met2)[1:ncol(met2)-1])]
write.csv(met2, "data/metabolomics_hmdb.csv")

met3 <- met2 %>% filter(group == "TDP43")
met3$group <- meta$Sex[which(meta$ID == rownames(met3))]
write.csv(met3, "data/metabolomics_hmdb_tdp43.csv")

met <- apply(t(met), 2, function(x) log10(x + 0.01))
rownames(met) <- IDs
met <- as.data.frame(met)
head(met)
dim(met)
any(sapply(met, function(x) sum(is.na(x))))
sapply(met, function(x) sd(x))
saveRDS(met, "data/metabolomics.RDS")
