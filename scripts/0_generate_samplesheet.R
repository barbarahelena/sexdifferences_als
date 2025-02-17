# Generate sample sheet
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

# Library
library(tidyverse)

# Data
df <- rio::import("data/sample_info_valeriia.csv")
path <- "/omics/groups/OE0554/internal_temp/barbara/projects/als/rawdata/"
ext <- ".fastq.gz"

# Sample sheet
head(df)
sample <- df$ID
fastq_1 <- str_c(path, df$ID, ext)
sample; fastq_1
samplesheet <- data.frame(sample, fastq_1)
head(samplesheet)
write.csv(samplesheet, "data/samplesheet2.csv", row.names = FALSE)
