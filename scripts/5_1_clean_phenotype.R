# Weights
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

# Libraries
library(readxl)
library(tidyverse)
library(lubridate)

# Data
meta <- readRDS("data/meta_microbiome_run1.RDS") %>% mutate(MouseID = as.numeric(MouseID))
df <- read_xlsx("data/weight.xlsx", skip = 1)
names(df)
colnames(df) <- str_remove(colnames(df), "Weight_")

df <- df %>% 
    filter(ID %in% meta$MouseID) %>%
    pivot_longer(., 4:ncol(.), names_to = "date", values_to = "weight")
head(df)
df2 <- df %>% 
  mutate(
    # Clean up date string
    date = str_replace(date, "([A-Za-z]+)(\\d+)_(\\d+)", "\\1 \\2 \\3"),
    # Parse as date using mdy
    date = mdy(date),
    # Convert birthdate to date
    birthdate = ymd(Birthdate),
    # Calculate age in weeks
    age_weeks = as.numeric(difftime(date, birthdate, units = "weeks"))
  )
head(df2)

ggplot(data = df2, aes(x = age_weeks, y = weight)) + 
    geom_jitter() + 
    theme_minimal()

df3 <- df2 %>%
  mutate(age_weeks_rounded = round(age_weeks)) %>%
  group_by(ID, age_weeks_rounded) %>%
  summarise(
    weight = mean(weight, na.rm = TRUE),
    n_measurements = n()
  ) %>%
  ungroup() %>%
  mutate(MouseID = ID,) %>%
  rename(Age_ints = age_weeks_rounded) %>%
  select(-ID)
saveRDS(df3, "data/weights.RDS")

# Merge weights with metadata
head(df3)
tot <- left_join(df3, meta, by = c("MouseID", "Age_ints"))
head(tot)
saveRDS(tot, "data/weights_metadata.RDS")

# Open onset data
df_onset <- read_xlsx("data/data_onset_death.xlsx")
head(df_onset)
tot <- left_join(meta, df_onset, by = "MouseID", relationship = "many-to-one")
head(tot)
saveRDS(tot, "data/onset_metadata.RDS")
