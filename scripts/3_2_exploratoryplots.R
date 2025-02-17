# Exploratory plots pathways ALS
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

## Load libraries
library(tidyverse)
library(rio)
library(forcats)

theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(family = 'Helvetica'),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.9)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text.x = element_text(size = rel(0.9)),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
} 

## Data
df <- readRDS("data/pathways.RDS")
head(df)
keypath <- readRDS("data/pathwaykeys.RDS")
meta <- readRDS("data/metadata.RDS")

## Bar plots pathways
df_mean <- df %>% summarise_all(mean, na.rm = TRUE)
df_sd <- df %>% summarise_all(sd, na.rm = TRUE)

df_summary <- bind_rows(
  df_mean %>% mutate(stat = "mean"),
  df_sd %>% mutate(stat = "sd")
)
dim(df_summary) # 235 pathways

# Reshape df_summary to long format for plotting
df_long <- df_summary %>%
  pivot_longer(-stat, names_to = "pathway", values_to = "value")

# Merge with keypath to replace pathway names with explanations
df_long <- df_long %>%
  left_join(keypath, by = c("pathway" = "keys")) %>%
  mutate(pathway = expl) %>%
  select(-expl)

# Filter top 20 pathways by mean value
top_20_pathways <- df_long %>%
  filter(stat == "mean") %>%
  arrange(desc(value)) %>%
  slice(1:20) %>%
  pull(pathway)

df_long2 <- df_long %>% 
  filter(stat == "mean" & pathway %in% top_20_pathways) %>% 
  mutate(pathway = fct_reorder(pathway, value))

# Create bar plots arranged by mean pathways
ggplot(df_long2, aes(x = pathway, y = value)) +
  geom_bar(stat = "identity", fill = "royalblue") +
  coord_flip() +
  labs(title = "Top 20 Pathways", x = "Pathway", y = "log10(cpm)") +
  theme_Publication()
ggsave("results/pathways/top20.pdf", width = 7, height = 6)
