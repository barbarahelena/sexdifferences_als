# Impute microbiome data
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)

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

# Function to impute missing values using the average of the week before and after
impute_missing <- function(x) {
  n <- length(x)
  if (n < 3) return(x)  # If there are less than 3 elements, return as is
  
  for (i in 2:(n - 1)) {
    if (is.na(x[i]) & !is.na(x[i - 1]) & !is.na(x[i + 1])) {
      x[i] <- mean(c(x[i - 1], x[i + 1]), na.rm = TRUE)
    } else if (is.na(x[i]) & !is.na(x[i - 1])){
      x[i] <- x[i - 1] + (x[i - 1] - x[i - 2])
      if(x[i] < 0) x[i] <- x[i - 1]
    }
  }
  return(x)
}

# Open data
df <- readRDS("data/filtered_microbiome_data.RDS")
names(df)
microbiome_cols <- names(df)[c(15:ncol(df)-3)] # last three cols are metadata
head(microbiome_cols); tail(microbiome_cols)

# Pivot to long format for imputation
dftot_long <- df %>%
  pivot_longer(., cols = all_of(microbiome_cols), names_to = "microbe", values_to = "abundance") %>%
  arrange(MouseID, Age_ints, microbe) %>%
  mutate(Age_weeks = fct_inorder(as.factor(Age_weeks)), 
         MouseID = fct_inorder(as.factor(MouseID)),MouseID) 

# Apply the imputation function to each microbiome variable
dftot_long_imputed <- dftot_long %>%
  group_by(MouseID, microbe) %>%
  arrange(Age_weeks) %>%
  mutate(abundance = impute_missing(abundance)) %>%
  ungroup()

# And pivot back to wide format
df_imputed_wide <- dftot_long_imputed %>%
  pivot_wider(names_from = microbe, values_from = abundance) %>%
  select(ID, MouseID, Age_weeks, Age_ints, Genotype, Sex, GenotypePerSex, all_of(microbiome_cols)) %>%
  mutate(ID = case_when(is.na(ID) ~ str_c("I", row_number()), 
            .default = ID))

# Filter out timepoint 18 weeks - too many missings (see script 2_1)
df_imputed_wide <- df_imputed_wide %>% filter(!Age_ints %in% c(18))
## This filter is applied last so that the timepoints could still be used to impute neigbouring timepoints ##
any(is.na(df_imputed_wide)) # FALSE - There's no missing values left

# Print the resulting data frame and check for 2 mice
df %>% # BEFORE IMPUTATION, including week 8 and 18
       filter(MouseID == "80") %>% 
       select(microbiome_cols[15], Age_weeks) %>% head(n = 30)
df_imputed_wide %>% # AFTER IMPUTATION
       filter(MouseID == "80") %>% 
       select(microbiome_cols[15], Age_weeks, GenotypePerSex) %>% 
       head(n = 8)

df %>% # BEFORE IMPUTATION, including week 8 and 18
       filter(MouseID == "79") %>% 
       select(microbiome_cols[15], Age_weeks) %>% head(n = 30)
df_imputed_wide %>% # AFTER IMPUTATION
       filter(MouseID == "79") %>% 
       select(microbiome_cols[15], Age_weeks, GenotypePerSex) %>% 
       head(n = 8)

# Save the imputed df
write.csv(df_imputed_wide, "data/imputed_microbiome_data.csv", row.names = FALSE)
saveRDS(df_imputed_wide, "data/imputed_microbiome_data.RDS")

# Transform imputed df
log_and_normalize <- function(x, age_ints) {
  week6_value <- x[age_ints == 6]
  log2((x + 0.001) / (week6_value + 0.001))
}
rowSums(df_imputed_wide[, 8:141])

# Apply the log2 transformation and normalization
data_normalized <- df_imputed_wide %>%
  arrange(MouseID, Age_ints) %>%
  group_by(MouseID) %>%
  mutate(across(c("Muribaculaceae bacterium Isolate-037 (Harlan)":"Synergistes sp. Zagget9"), 
          ~ log_and_normalize(.x, Age_ints))) %>%
  ungroup()
head(data_normalized)[1:5, c(2,3,9)]
head(df_imputed_wide %>% arrange(MouseID))[1:5, c(2,3,9)]

# Save the normalized-to-6wk df
write.csv(data_normalized, "data/imputed_norm_microbiome_data.csv", row.names = FALSE)
saveRDS(data_normalized, "data/imputed_norm_microbiome_data.RDS")

# Lineplot before and after imp as example
microbiome_var <- microbiome_cols[1]  # Select the first microbiome variable for plotting
data_before_imputation <- df %>%
  arrange(Age_ints) %>%
  mutate(Age_weeks = fct_inorder(as.factor(Age_weeks))) %>%
  filter(MouseID == "79") %>%
  select(MouseID, Age_weeks, all_of(microbiome_var), Age_ints) %>%
  filter(!Age_ints %in% c(18))
  
data_after_imputation <- df_imputed_wide %>% filter(MouseID == "79") %>% 
       select(MouseID, Age_weeks, any_of(microbiome_var), Age_ints)

(pl1 <- ggplot(data = data_before_imputation, 
       aes(x = Age_ints, y = .data[[microbiome_var]], na.rm = FALSE)) +
  geom_line(aes(group = MouseID)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  labs(title = "Sample 79 before imputation",
       x = "Age (weeks)",
       y = str_c("Relative abundance ", microbiome_var)) +
  theme_Publication())

(pl2 <- ggplot(data = data_after_imputation, 
        aes(x = Age_ints, y = .data[[microbiome_var]], na.rm = FALSE)) +
  geom_line(aes(group = MouseID)) +
  geom_point() +
  labs(title = "Sample 79 after imputation",
       x = "Age (weeks)",
       y = str_c("Relative abundance ", microbiome_var)) +
  theme_Publication())

ggpubr::ggarrange(pl1, pl2, labels = c("A", "B"), ncol = 2, nrow = 1)

# Save the plot
ggsave("results/microbiome/sample_line_plot.pdf", width = 10, height = 6)

