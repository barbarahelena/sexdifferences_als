# Table 1
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(tableone)

# Data
meta <- readRDS("data/human_cohort2/metadata.RDS")
head(meta)
dim(meta)

# Table 1
table_one <- CreateTableOne(
  vars = c("Sex", "Age", "ALSFRS", "Mutation"),
  strata = "Group",
  data = meta,
  factorVars = c("Sex", "Mutation"),
  test = FALSE
)
table_one_csv <- print(table_one, nonnormal = c("ALSFRS"), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table_one_csv, "results/humancohort2/table_one_group.csv", row.names = TRUE)