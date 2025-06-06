# Clean microbiome data ALS project
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)
library(rio)
library(forcats)

# Function for plots
theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.5, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

# Function to add GenotypePerSex, Sex, and Genotype to the missing samples
add_metadata <- function(mbdata, metadata, metadata_cols) {
  # Filter the original data to get the metadata columns
  metadata <- meta %>%
    filter(Age_weeks == "6 weeks") %>%
    select(MouseID, all_of(metadata_cols)) %>%
    filter(!duplicated(MouseID))
  # Merge the metadata with the microbiome data, redefine metadata columns
  mbdata_with_metadat <- mbdata %>% select(-Sex, -Genotype, -GenotypePerSex) %>%
    left_join(., metadata, by = "MouseID")
  return(mbdata_with_metadat)
}

# Metadata table
meta <- rio::import("data/sample_info_valeriia.csv") %>% 
                    filter(Run == "Run1" | Run == "Run2") %>% 
                    mutate(GenotypePerSex = str_c(Sex, " ", Genotype),
                           GenotypePerSex = fct_relevel(GenotypePerSex, "Female WT", after = 0L),
                           GenotypePerSex = fct_relevel(GenotypePerSex, "Male WT", after = 1L),
                            Age_weeks = fct_relevel(Age_weeks, "8 weeks", after = 0L),
                            Age_weeks = fct_relevel(Age_weeks, "6 weeks", after = 0L),
                            LowCount = case_when(counts < 1000000 ~ TRUE,
                                                 .default = FALSE))
lowcountids <- meta$ID[which(meta$LowCount == TRUE)]

saveRDS(meta, "data/meta_microbiome.RDS")
write.csv(meta, "data/meta_microbiome.csv")
meta2 <- meta %>% filter(str_detect(ID, "L")) %>% droplevels(.) # Run 1
saveRDS(meta2, "data/meta_microbiome_run1.RDS")
write.csv(meta2, "data/meta_microbiome_run1.csv")
meta3 <- meta %>% filter(str_detect(ID, "S")) %>% droplevels(.) # Run 2
saveRDS(meta3, "data/meta_microbiome_run2.RDS")
write.csv(meta3, "data/meta_microbiome_run2.csv")

# Tidy bracken table
tab <- rio::import("data/combined_brackenoutput.txt")
dim(tab) 
names(tab) 
tab <- tab %>% dplyr::select(name, taxonomy_lvl, contains("_frac")) %>% 
                    filter(taxonomy_lvl == "S") %>% 
                    dplyr::select(-taxonomy_lvl) %>% 
                    mutate(across(contains("_frac"), ~ .x * 100)) %>% # make %
                    rename_with(~ str_remove(., ".tsv_frac"), 2:ncol(.))
dim(tab) # 323 samples
tab <- tab %>% dplyr::select(-all_of(lowcountids))
names(tab)
dim(tab) # 318 samples
rownames(tab) <- tab$name
tab$name <- NULL
df <- as.data.frame(t(as.matrix(tab)))
head(df)[1:5,1:5]
df2 <- df %>% filter(rownames(.) %in% meta2$ID) # Run 1
dim(df2) # 137 samples
df3 <- df %>% filter(rownames(.) %in% meta3$ID) # Run 2
dim(df3) # 144 samples

saveRDS(df2, "data/microbiome_run1.RDS")
write.csv(df2, "data/microbiome_run1.csv")
saveRDS(df3, "data/microbiome_run2.RDS")
write.csv(df3, "data/microbiome_run2.csv")

# Pruning of species
## Run 1
threshold <- 0.05 # 0.05% threshold
min_samples <- 0.1 * nrow(df2) # in 10% of samples
species_filter <- apply(df2, 2, function(x) sum(x >= threshold) >= min_samples)
df_filtered2 <- df2[, species_filter]
dim(df_filtered2) # 134 species and 137 samples

## Run 2
threshold <- 0.05 # 0.05% threshold
min_samples <- 0.1 * nrow(df3) # in 10% of samples
species_filter <- apply(df3, 2, function(x) sum(x >= threshold) >= min_samples)
df_filtered3 <- df3[, species_filter]
dim(df_filtered3) # 169 species and 144 samples

saveRDS(df_filtered2, "data/microbiome_filtered_run1.RDS")
write.csv(df_filtered2, "data/microbiome_filtered_run1.csv")
saveRDS(df_filtered3, "data/microbiome_filtered_run2.RDS")
write.csv(df_filtered3, "data/microbiome_filtered_run2.csv")

# Plot a bar plot of available samples
## In total per timepoint per group
sample_counts <- meta3 %>%
  group_by(GenotypePerSex, Age_ints) %>%
  summarise(count = n()) %>%
  ungroup()
ggplot(sample_counts, aes(x = Age_ints, y = count, fill = GenotypePerSex)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of samples",
       x = "Age (weeks)",
       y = "Number of samples") +
  ggsci::scale_fill_nejm() +
  scale_x_continuous(n.breaks = 8) +
  theme_Publication() +
  theme(legend.title = element_blank())
ggsave("results/microbiome/sample_counts_barplot_run2.pdf", width = 10, height = 6)
## Per mouse
sample_counts <- meta3 %>%
  group_by(MouseID, Age_ints) %>%
  summarise(count = n(), GenotypePerSex = GenotypePerSex) %>%
  ungroup()
ggplot(sample_counts, aes(x = Age_ints, y = count, fill = GenotypePerSex)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of samples",
       x = "Age (weeks)",
       y = "Number of samples") +
  ggsci::scale_fill_nejm() +
  # scale_x_continuous(breaks = c(6,8,10,12,14,16,18,20,22,24)) +
  facet_wrap(~MouseID) +
  theme_Publication() +
  theme(legend.title = element_blank())
ggsave("results/microbiome/sample_counts_barplot_id_run2.pdf", width = 15, height = 15)

# Filter mice that are missing more than 2 timepoints
## First make a complete dataset: microbiome and metadata
df_filtered2$ID <- rownames(df_filtered2)
meta2 <- meta2 %>% filter(Age_ints < 20) %>% droplevels(.) # from week 20, a lot of samples missing
all_levels <- levels(meta2$Age_weeks)
table(meta2$Age_weeks)
dftot <- left_join(meta2, df_filtered2, by = "ID") # merge metadata with abundance table

## Add missing timepoints to table and fill in metadata for these samples
metadata_cols <- c("GenotypePerSex", "Sex", "Genotype") # to add to missing samples below
dftot_compl <- dftot %>% 
  group_by(MouseID) %>% 
  complete(Age_weeks = all_levels) %>% # complete levels to make missing timepoints explicit
  add_metadata(dftot_compl2, metadata_cols) %>% # see function at the top
  mutate(Age_ints = as.integer(str_remove(Age_weeks, " weeks"))) %>% # redefine for missing samples
  arrange(Age_ints) %>%
  mutate(Age_weeks = fct_inorder(as.factor(Age_weeks)), # so that it sorts nicely later
         MouseID = fct_inorder(as.factor(MouseID)))

dim(dftot_compl) # 140 samples (3 samples added, by filling in the missing timepoints)
table(dftot_compl$Age_weeks) # now each timepoint has 20 samples
table(dftot_compl$GenotypePerSex, dftot_compl$Age_weeks) 

## Filter out mice with more than 2 missing timepoints
missingdf <- dftot_compl %>% filter(!Age_ints %in% c(18)) # at week 18, many samples are missing
missingids <- missingdf %>% mutate(missing = sum(is.na(ID))) %>% filter(missing > 1) # get mice with >1 missing
unique(missingids$MouseID) # IDs 92 and 94 have 2 or more missing timepoints out of 5 total
all_levels <- levels(as.factor(dftot_compl$Age_weeks))
all_levels
dftot_compl2 <- dftot_compl %>%
       filter(!MouseID %in% missingids$MouseID) # based on missingdf, filter mice 92 and 94 out

dim(dftot_compl2) # 126 samples left (14 samples filtered out)
table(dftot_compl2$Age_weeks) # all timepoints have 18 samples
table(dftot_compl2$GenotypePerSex, dftot_compl2$Age_weeks) # the sample no is constant across timepoints

## Save the filtered data (no mice with > 2 missing timepoints, ignoring wk 18)
saveRDS(dftot_compl2, "data/filtered_microbiome_data.RDS")
write.csv(dftot_compl2, "data/filtered_microbiome_data.csv", row.names = FALSE)
