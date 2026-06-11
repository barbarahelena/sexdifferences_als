# Cleaning human cohort data - Calgary cohort
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggpubr)

# Data
meta <- rio::import("data/human_cohort2/metadata.xlsx") |> 
          mutate(
                  #ID = str_c("X", fecal_ID),
                  Group = case_when(
                      diagnosis == "control" ~ "Control",
                      diagnosis == "als" ~ "ALS"
                  ),
                  Group = as.factor(Group),
                  sex = fct_recode(sex, "Male" = "male", "Female" = "female")
                ) |> 
          dplyr::select(ID = fecal_ID, Household = household_ID, Group, Age = age_yrs,
                  Sex = sex, Duration = als_symptom_duration_months, 
                  Familial = sporadic_or_genetic, Mutation = genetic_mutation,
                  ALSFRS = alsfrs, Onset = onset_site) |> 
            mutate(Mutation = fct_relevel(Mutation, "none", after = 2L))
saveRDS(meta, "data/human_cohort2/metadata.RDS")
dim(meta)

# Microbiome - merge existing and new bracken profiles
tab_old <- rio::import("data/human_cohort2/combined_brackenoutput.txt")
tab_new <- rio::import("data/human_cohort2/combined_brackenoutput_newfiles.txt")

old_samples <- names(tab_old) |> str_subset("_frac$") |> str_remove("\\.tsv_frac")
new_samples <- names(tab_new) |> str_subset("_frac$") |> str_remove("\\.tsv_frac")
overlap_samp <- intersect(old_samples, new_samples)
if (length(overlap_samp) > 0) stop("Sample overlap in bracken files: ", paste(overlap_samp, collapse = ", "))
message("Bracken merge: ", length(old_samples), " existing + ", length(new_samples),
        " new = ", length(old_samples) + length(new_samples), " total samples")

tab <- full_join(tab_old, tab_new, by = c("name", "taxonomy_id", "taxonomy_lvl")) |>
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

stopifnot(ncol(tab |> dplyr::select(contains("_frac"))) == length(old_samples) + length(new_samples))
message("Bracken row count after merge: ", nrow(tab))

tab <- tab %>% dplyr::select(name, taxonomy_lvl, contains("_frac")) %>%
                    filter(name != "Homo sapiens") |>
                    dplyr::select(-taxonomy_lvl) %>%
                    mutate(across(contains("_frac"), ~ .x * 100)) %>% # make %
                    rename_with(~ str_remove(., ".tsv_frac"), 2:ncol(.)) %>%
                    rename_with(~ str_remove(., "_species"), 2:ncol(.))
dim(tab) # 40 samples
head(tab)[1:5,1:5]

readcounts <- rio::import("data/human_cohort2/readcounts_table.csv")
gghistogram(readcounts$ReadCount)
min(readcounts$ReadCount) # 7.7 million reads
rownames(tab) <- tab$name
tab$name <- NULL
df <- as.data.frame(t(as.matrix(tab)))
head(df)[1:5,1:5]
meta_missing <- meta$ID[which(!meta$ID %in% rownames(df))]
extra_in_df  <- rownames(df)[which(!rownames(df) %in% meta$ID)]
message("Samples in metadata but missing from bracken: ",
        if (length(meta_missing) == 0) "none" else paste(meta_missing, collapse = ", "))
message("Samples in bracken but missing from metadata: ",
        if (length(extra_in_df) == 0) "none" else paste(extra_in_df, collapse = ", "))
stopifnot(length(extra_in_df) == 0)
dim(df)
summary(as.factor(meta$Sex))
table(meta$Sex, meta$Group) # in ALS group 9 vs 30
df <- df[rownames(df) %in% meta$ID[which(!is.na(meta$Sex))],]
dim(df) # 40 samples

# Pruning of species - stratified by patient group - done in Linda so not needed here
# threshold <- 0.05 # 0.05% threshold
# min_proportion <- 0.3 # in 20% of samples within group

# tk <- apply(df, 2, function(x) sum(x >= threshold) > (min_proportion * nrow(df)))
# summary(tk)

saveRDS(df, "data/human_cohort2/microbiome.RDS")

# Pathways - merge existing and new HUMAnN tables
pathw_old <- rio::import("data/human_cohort2/pathway_abundance_cpm_unstratified.txt") |>
    dplyr::select(Pathways = `# Pathway HUMAnN v4.0.0.alpha.1`, everything())
pathw_new <- rio::import("data/human_cohort2/pathway_abundance_cpm_unstratified_newfiles.txt") |>
    dplyr::select(Pathways = `# Pathway HUMAnN v4.0.0.alpha.1`, everything())

old_pw_samps <- names(pathw_old)[-1] |> str_remove("_Abundance")
new_pw_samps <- names(pathw_new)[-1] |> str_remove("_Abundance")
overlap_pw <- intersect(old_pw_samps, new_pw_samps)
if (length(overlap_pw) > 0) stop("Sample overlap in pathway files: ", paste(overlap_pw, collapse = ", "))
message("Pathway merge: ", length(old_pw_samps), " existing + ", length(new_pw_samps),
        " new = ", length(old_pw_samps) + length(new_pw_samps), " total samples")

pathw <- full_join(pathw_old, pathw_new, by = "Pathways") |>
    mutate(across(where(is.numeric), ~ replace_na(.x, 0))) |>
    filter(Pathways != "UNMAPPED" & Pathways != "UNINTEGRATED") |>
    rename_with(~ str_remove(., "_Abundance"), contains("_Abundance"))

stopifnot((ncol(pathw) - 1) == length(old_pw_samps) + length(new_pw_samps))
message("Pathway rows after merge: ", nrow(pathw))
dim(pathw)
head(pathw)

paths <- pathw$Pathways
keys <- str_extract(paths, "[A-Z0-9-]+(?=:)")
expl <- trimws(str_extract(paths, "(?<=:).*"))
lib <- data.frame(paths, keys, expl)
pathw <- pathw |> mutate(
    Pathways = str_extract(Pathways, "[A-Z0-9-]+(?=:)"))
rownames(pathw) <- pathw$Pathways
pathw$Pathways <- NULL
pathdf <- as.data.frame(t(as.matrix(pathw)))

# Filter
pathway_filter <- apply(pathdf, 2, function(x) sum(x > 5) >= (0.5 * nrow(pathdf)))
pw_filtered <- pathdf[, pathway_filter]

# Log transform filtered pathways
pw <- log10(pw_filtered + 1)

# Print dimensions
dim(pw)
head(pw)

saveRDS(pw, "data/human_cohort2/human_als_pathways.RDS")
write.csv(pw, "data/human_cohort2/human_als_pathways.csv")
saveRDS(lib, "data/human_cohort2/pathwaykeys.RDS")
write.csv(lib, "data/human_cohort2/pathwaykeys.csv", row.names = FALSE)

lib %>% filter(str_detect(expl, "rhamnose"))
#                                                              paths            keys                                            expl
# 1             DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis DTDPRHAMSYN-PWY             dTDP-&beta;-L-rhamnose biosynthesis
# 2                            RHAMCAT-PWY: L-rhamnose degradation I     RHAMCAT-PWY                        L-rhamnose degradation I

# CAZyme sample statistics - merge existing and new
cay_stat_old <- rio::import("data/human_cohort2/cayman_sample_statistics.tsv")
cay_stat_new <- rio::import("data/human_cohort2/sample_statistics_newfiles.tsv")
overlap_stat <- intersect(cay_stat_old$sample, cay_stat_new$sample)
if (length(overlap_stat) > 0) stop("Sample overlap in cayman stats: ", paste(overlap_stat, collapse = ", "))
message("CAZyme stats merge: ", nrow(cay_stat_old), " existing + ", nrow(cay_stat_new),
        " new = ", nrow(cay_stat_old) + nrow(cay_stat_new), " total samples")
cayman_stat <- bind_rows(cay_stat_old, cay_stat_new)
stopifnot(nrow(cayman_stat) == nrow(cay_stat_old) + nrow(cay_stat_new))
stopifnot(all(names(cay_stat_old) == names(cay_stat_new)))
saveRDS(cayman_stat, "data/human_cohort2/cayman_sample_statistics.RDS")

# CAZyme families - merge existing and new
cay_old <- rio::import("data/human_cohort2/cayman_families_cpm_table.tsv")
cay_new <- rio::import("data/human_cohort2/families_cpm_table_newfiles.tsv")
old_cay_samps <- names(cay_old)[-1]
new_cay_samps <- names(cay_new)[-1]
overlap_cay <- intersect(old_cay_samps, new_cay_samps)
if (length(overlap_cay) > 0) stop("Sample overlap in cayman families: ", paste(overlap_cay, collapse = ", "))
message("CAZyme families merge: ", length(old_cay_samps), " existing + ", length(new_cay_samps),
        " new = ", length(old_cay_samps) + length(new_cay_samps), " total samples")
cayman_fam <- full_join(cay_old, cay_new, by = "family") |>
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
stopifnot((ncol(cayman_fam) - 1) == length(old_cay_samps) + length(new_cay_samps))
message("CAZyme families rows after merge: ", nrow(cayman_fam))
saveRDS(cayman_fam, "data/human_cohort2/cayman_families_merged.RDS")