# Table 1
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

# Libraries
library(tidyverse)
library(tableone)

# Data
meta <- readRDS("data/human_cohort/metadata.RDS")

# Table 1
table_one <- CreateTableOne(
  vars = c("Sex", "Age", "ALS", "BMI"),
  strata = "Group",
  data = meta,
  factorVars = "Sex",
  test = FALSE
)
table_one_csv <- print(table_one, nonnormal = c("ALS"), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table_one_csv, "results/humancohort/table_one_group.csv", row.names = TRUE)

table_one <- CreateTableOne(
  vars = c("Age", "ALS", "BMI"),
  strata = "Sex",
  data = meta %>% filter(Group == "ALS"),
  factorVars = "Sex",
  test = FALSE
)
table_one_csv <- print(table_one, nonnormal = c("ALS"), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table_one_csv, "results/humancohort/table_one_onlyals.csv", row.names = TRUE)
