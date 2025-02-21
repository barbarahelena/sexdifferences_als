## Enrichment analyses metabolomics
## Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

# Libraries
library(ggplot2)
library(tidyverse)
library(ggsci)

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

# Load metabolomics data
met <- readRDS("data/metabolomics.RDS")
meta <- readRDS("data/metadata.RDS")
hmdb <- readRDS("data/hmdb_ids.RDS") %>% mutate(metabolite = str_to_lower(metabolite))   # hmdb IDs per metabolite
colnames(met) <- hmdb$HMDB[which(hmdb$metabolite == colnames(met))]

# # Test RaMP #
# libary(RaMP)
# # Load differential metabolite results (assuming it has metabolite IDs and p-values)
# diff_metabolites_genotype <- rio::import("results/metabolomics/ttests/metabolites_welcht_mice_diff.csv")
# diff_metabolites_wt_sex <- rio::import("results/metabolomics/ttests/metabolites_welcht_diff.csv")
# diff_metabolites_tdp_sex <- rio::import("results/metabolomics/ttests/metabolites_welcht_ctrl_diff.csv")

# hmdb <- hmdb %>% 
#     filter(metabolite %in% c(diff_metabolites_genotype$metabolite, 
#                                 diff_metabolites_tdp_sex$metabolite, 
#                                 diff_metabolites_wt_sex$metabolite))

# # Connect to RaMP
# db_path <- "RaMP.sqlite"  # Adjust if needed
# con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

# # Perform pathway enrichment analysis
# enrichment_results <- RaMP::getPathwayFromAnalyte(analytes = str_c("hmdb:",hmdb$HMDB), NameOrIds = "ids")

# fisher.results <- runCombinedFisherTest(analytes = str_c("hmdb:",hmdb$HMDB))
# filtered.fisher.results <- FilterFishersResults(fisher.results, pval_type = 'holm', pval_cutoff=0.05)
# clusters <- RaMP::findCluster(fishers_df = filtered.fisher.results,
#   perc_analyte_overlap = 0.2,
#   min_pathway_tocluster = 2, perc_pathway_overlap = 0.2
# )
# data.table(clusters$fishresults %>% mutate_if(is.numeric, ~ round(., 8)),
#   rownames = FALSE
# )
# pathwayResultsPlot(pathwaysSig = filtered.fisher.results, text_size = 8, perc_analyte_overlap = 0.2, 
# 	min_pathway_tocluster = 2, perc_pathway_overlap = 0.2, interactive = FALSE)

# # Disconnect database
# DBI::dbDisconnect(con)

### Try MetaboAnalyst ### cannot install...
library(MetaboAnalystR)
# mSet <- InitDataObjects("conc", "pathora", FALSE)

# Read in data table
mSet<-Read.TextData(mset, filePath = "data/metabolomics_hmdb.csv", "rowu", "disc");

# Perform cross-referencing of compound names
mSet <- CrossReferencing(mSet, "hmdb")
mSet <- SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);

# Perform no normalization
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, rowNorm="NULL", transNorm="NULL", scaleNorm="NULL", ref="WT")
# mSet <- PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
# mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
# Set the metabolome filter
mSet<-SetMetabolomeFilter(mSet, F);
# Set the metabolite set library to pathway
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway")
# Calculate the global test score
mSet<-CalculateGlobalTestScore(mSet)
# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_genotype_", "bar", "png", 72, width=NA)

# Read in data table
mSet<-Read.TextData(mset, filePath = "data/metabolomics_hmdb_tdp43.csv", "rowu", "disc");

# Perform cross-referencing of compound names
mSet <- CrossReferencing(mSet, "hmdb")
mSet <- SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);

# Perform no normalization
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, rowNorm="NULL", transNorm="NULL", scaleNorm="NULL", ref="Female")
# mSet <- PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
# mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
# Set the metabolome filter
mSet<-SetMetabolomeFilter(mSet, F);
# Set the metabolite set library to pathway
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway")
# Calculate the global test score
mSet<-CalculateGlobalTestScore(mSet)
# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_tdp43_", "bar", "png", 72, width=NA)