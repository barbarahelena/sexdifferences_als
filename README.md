# Gut Microbiota and Metabolome Analysis in TDP43 Mice

## Introduction
In this project, we investigated the role of the gut microbiota and metabolome in driving sexual dimorphism in ALS-linked TDP43 mice. Using shotgun metagenomics and semi-targeted metabolomics, we analyzed the microbiome and metabolome profiles of WT and TDP43 mice over time, focusing on genotype and sex differences. The repository includes scripts for microbiome composition analysis, pathway analysis and metabolomics analysis.

## Repository Structure
#### Chapter 1: Metabolomics Analysis
1. **`1_1_clean_metabolomics.R`**  
   - Cleans and processes metabolomics data.  
   - Outputs: `data/metabolomics.RDS`, `data/metabolomics_hmdb.csv`.
2. **`1_2_pca_metabolomics.R`**  
   - Performs PCA on metabolomics data.  
   - Outputs: PCA plots of metabolite profiles.
3. **`1_3_welch_metabolomics.R`**  
   - Performs Welch's t-tests for metabolite differences by genotype and sex.  
   - Outputs: Volcano plots, significant metabolite tables.
4. **`1_4_volcanoplot.R`**  
   - Generates volcano plots for metabolomics data.  
   - Outputs: Volcano plots for genotype and sex differences.
5. **`1_5_boxplots.R`**  
   - Creates boxplots for metabolites showing significant differences.  
   - Outputs: Boxplots for selected metabolites.
6. **`1_6_heatmap.R`**  
   - Generates heatmaps for metabolomics data.  
   - Outputs: Heatmaps of metabolite concentrations.

#### Chapter 2: Microbiome Analysis
1. **`2_1_clean_microbiome.R`**  
   - Cleans and processes microbiome data.  
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
   - Prepares microbiome data for TCAM analysis.  
   - Outputs: Imputed microbiome data.
6. **`2_6_tcam.py`**  
   - Performs TCAM analysis and generates scatterplots.  
   - Outputs: TCAM scatterplots, loadings, and PERMANOVA results.
7. **`2_7_tcamplots.R`**  
   - Visualizes TCAM results and performs PERMANOVA for genotype and sex differences.  
   - Outputs: TCAM plots, PERMANOVA results.

#### Chapter 3: Pathway Analysis
1. **`3_1_clean_pathways.R`**  
   - Cleans and processes pathway data.  
   - Outputs: `data/pathways.RDS`, `data/pathwaykeys.RDS`.
2. **`3_2_pathways_microbes.R`**  
   - Correlates microbial species with pathways and visualizes results.  
   - Outputs: Correlation heatmaps, pathway boxplots.

### Outputs
- **Microbiome Analysis**:  
  - Compositional plots, alpha diversity plots, beta diversity plots, TCAM scatterplots, and loadings.
- **Pathway Analysis**:  
  - Correlation heatmaps, pathway boxplots.
- **Metabolomics Analysis**:  
  - PCA plots, volcano plots, heatmaps, and boxplots.
- **Statistical Results**:  
  - PERMANOVA results, Welch's t-test results, significant metabolite/pathway tables.

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