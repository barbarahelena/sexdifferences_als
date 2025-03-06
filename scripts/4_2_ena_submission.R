# ENA submission
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Library
library(tidyverse)

# Data
meta <- readRDS("data/meta_microbiome_run1.RDS")

templ <- meta %>% mutate(tax_id = "10090",
                            scientific_name = "Mus musculus", 
                            sample_alias = ID,
                            sample_title = ID,
                            sample_description = Genotype,
                            `collection date` = "2020",
                            sex = Sex,
                            `geographic location (country and/or sea)` = "Canada") %>%
                    select(tax_id, scientific_name, sample_alias, sample_title, sample_description,
                            `collection date`, sex, `geographic location (country and/or sea)`)
head(templ)
write.table(templ, "results/ENA/samplesheet.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

fastqs <- meta %>% mutate(
        sample = ID,
        study = "PRJEB86491",
        instrument_model = "Illumina NovaSeq 6000",
        library_name = "",
        library_source = "METAGENOMIC",
        library_selection = "PCR",
        library_strategy = "WGS",
        library_layout = "SINGLE",
        file_name = str_c(ID, '.fastq.gz'),
        file_md5 = sapply(ID, function(x) {
            md5_file <- paste0("rawdata/ena/", x, ".fastq.gz.md5")
            md5_content <- readLines(md5_file)
            strsplit(md5_content, " ")[[1]][1]
        })
        ) %>%
        select(sample, study, instrument_model, library_name,
                library_source, library_selection, library_strategy,
                library_layout, file_name, 
                file_md5
                )
write.table(fastqs, "results/ENA/fastqsheet.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
