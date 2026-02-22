# Gut Microbiota and Metabolome Analysis in TDP43 Mice

## Introduction
In this project, we investigated the role of the gut microbiota and metabolome in driving sexual dimorphism in ALS-linked TDP43 mice. Using shotgun metagenomics and semi-targeted metabolomics, we analyzed the microbiome and metabolome profiles of WT and TDP43 mice over time, focusing on genotype and sex differences. The repository also includes scripts for the validation human ALS cohort (Calgary cohort). Scripts that are no longer part of the main analysis are stored in `scripts/archive/`.

## Repository Structure

#### Chapter 1: Human Cohort 2 (Calgary ALS Cohort)
1. **`1_1_clean_humancohort2.R`**
   - Cleans and processes metadata for the Calgary human ALS cohort.
   - Outputs: `data/human_cohort2/metadata.RDS`.
2. **`1_2_tableone.R`**
   - Generates Table 1 with demographic and clinical characteristics.
   - Outputs: Summary table.
3. **`1_3_composition.R`**
   - Analyzes microbiome composition and generates compositional plots for cohort 2.
   - Outputs: Bar plots of microbiome composition.
4. **`1_4_diversity.R`**
   - Calculates alpha and beta diversity metrics and tests group differences.
   - Outputs: Diversity plots.
5. **`1_4_pathways.R`**
   - Analyzes microbial pathways in human hosts.
   - Outputs: Pathway plots.
6. **`1_5_microbes.R`**
   - Runs linear mixed models (LMMs) for microbial species differences.
   - Outputs: LMM results and plots.
7. **`1_6_pathwaycorr.R`**
   - Correlates microbial pathways with clinical variables in the Calgary cohort.
   - Outputs: Pathway correlation plots.
8. **`1_7_cayman.R`**
   - Analyzes Cayman metabolomics data for human cohort 2.
   - Outputs: Metabolomics plots and statistics.

#### Chapter 2: Mouse Microbiome Analysis
1. **`2_1_clean_microbiome.R`**
   - Cleans and processes microbiome data from TDP43 mouse experiments.
   - Outputs: `data/microbiome_filtered_run1.RDS`, `data/microbiome_filtered_run2.RDS`.
2. **`2_2_composition.R`**
   - Analyzes microbiome composition and generates compositional plots.
   - Outputs: Bar plots of microbiome composition.
3. **`2_3_alphadiversity.R`**
   - Calculates alpha diversity metrics (Shannon, Simpson) and tests group differences.
   - Outputs: Alpha diversity plots.
4. **`2_4_betadiversity.R`**
   - Calculates Bray-Curtis distances and performs PCoA.
   - Outputs: Beta diversity plots.
5. **`2_5_tcam_prep.R`**
   - Imputes and prepares microbiome data for TCAM analysis.
   - Outputs: Imputed microbiome data matrices.
6. **`2_6_tcam.py`**
   - Performs TCAM (Tensor Component Analysis of Microbiome) analysis.
   - Outputs: TCAM scatterplots, loadings, and PERMANOVA results.
7. **`2_7_tcamplots.R`**
   - Visualizes TCAM results and performs PERMANOVA for genotype and sex differences.
   - Outputs: TCAM plots, PERMANOVA results.

#### Chapter 3: Mouse Pathway and CAZyme Analysis
1. **`3_1_clean_pathways.R`**
   - Cleans and processes HUMAnN pathway abundance data.
   - Outputs: `data/pathways.RDS`, `data/pathwaykeys.RDS`.
2. **`3_2_pathways_microbes.R`**
   - Correlates microbial species with pathways and visualizes results.
   - Outputs: Correlation heatmaps, pathway boxplots.
3. **`3_3_cazymes.R`**
   - Analyzes carbohydrate-active enzymes (CAZymes) in TDP43 mouse experiments.
   - Outputs: CAZyme abundance plots and statistics.

#### Chapter 4: Mouse Metabolomics Analysis
1. **`4_1_clean_metabolomics.R`**
   - Cleans and processes metabolomics data.
   - Outputs: `data/metabolomics.RDS`, `data/metabolomics_hmdb.csv`.
2. **`4_2_pca_metabolomics.R`**
   - Performs PCA on metabolomics data.
   - Outputs: PCA plots of metabolite profiles.
3. **`4_3_welch_metabolomics.R`**
   - Performs Welch's t-tests for metabolite differences by genotype and sex.
   - Outputs: Statistical results tables.
4. **`4_4_volcanoplot.R`**
   - Generates volcano plots for metabolomics comparisons.
   - Outputs: Volcano plots for genotype and sex differences.
5. **`4_5_boxplots.R`**
   - Creates boxplots for metabolites showing significant differences.
   - Outputs: Boxplots for selected metabolites.
6. **`4_6_heatmap.R`**
   - Generates heatmaps for metabolomics data.
   - Outputs: Heatmaps of metabolite concentrations.

#### Chapter 5: ENA Submission
1. **`5_1_ena_checksum.sh`**
   - Computes MD5 checksums for raw FASTQ files prior to ENA upload.
2. **`5_2_ena_submission.R`**
   - Prepares submission metadata for upload to the European Nucleotide Archive (ENA).

#### Pipeline Scripts
- **`run_pipeline.sh`** — Runs the metagenomics processing pipeline (Nextflow) for TDP43 mouse samples.
- **`run_pipeline_human.sh`** — Runs the metagenomics processing pipeline for the human ALS cohort.
- **`check_metagenomic_ids.sh`** — Checks and lists unique sample IDs in the metagenomic data subset.

### Outputs
- **Human Cohort 2**: Table 1, compositional plots, diversity plots, pathway plots, LMM results, metabolomics plots.
- **Microbiome Analysis**: Compositional plots, alpha/beta diversity plots, TCAM scatterplots and loadings.
- **Pathway Analysis**: Correlation heatmaps, pathway boxplots, CAZyme plots.
- **Metabolomics Analysis**: PCA plots, volcano plots, heatmaps, boxplots.
- **Statistical Results**: PERMANOVA results, Welch's t-test results, LMM results.

## Requirements
### R Packages
- `tidyverse`
- `ggplot2`
- `ggpubr`
- `ComplexHeatmap`
- `circlize`
- `vegan`
- `ggsci`
- `rio`
- `gridExtra`
- `ape`
- `lme4`
- `afex`
- `tableone`
- `ggrepel`

### Python Packages
- `numpy`
- `pandas`
- `seaborn`
- `matplotlib`
- `scipy`
- `mprod`

## Code Availability
All scripts and workflows are available in this repository.
The metagenomics processing workflow is available at: [https://github.com/barbarahelena/metagenomicspipeline](https://github.com/barbarahelena/metagenomicspipeline).

## Data Availability
The raw metagenomics sequencing data is deposited in the European Nucleotide Archive (ENA) under accession number [PRJEB86491](https://www.ebi.ac.uk/ena/browser/view/PRJEB86491).
