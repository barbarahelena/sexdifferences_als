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
   - Calculates beta diversity (Bray-Curtis PCoA) and tests group differences by sex in ALS and control groups.
   - Outputs: PCoA plots (`pl1`, `pl2`).
5. **`1_5_diffabundance.R`**
   - Runs LinDA (Linear models for Differential Abundance) to identify microbial species differing by sex in ALS and control groups.
   - Outputs: Volcano plots (`p_li_als_sex`, `p_li_ctrl_sex`).
6. **`1_6_pathwaycorr.R`**
   - Correlates microbial pathways with clinical variables in the Calgary cohort and generates a heatmap.
   - Outputs: `heatmap_top`, `lgd_sig_path`.
7. **`1_7_cazymes_stats.R`**
   - Exploratory statistical analysis of carbohydrate-active enzymes (CAZymes) in the Calgary cohort.
   - Outputs: CAZyme statistics and exploratory plots.
8. **`1_8_cazymes.R`**
   - Analyzes CAZyme profiles in the Calgary cohort; generates PCoA plots and boxplots for key CAZyme families (GH78, GH106).
   - Outputs: `pl_caz1`, `pl_caz2`, `plist`, `families_present`.
9. **`1_9_assembleplot.R`**
   - Assembles all human cohort 2 analysis outputs into a single multi-panel figure.
   - Outputs: `results/humancohort2/assembled_figure.pdf`.

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
   - Outputs: PCA plots of metabolite profiles by sex and genotype.
3. **`4_3_heatmap.R`**
   - Generates a heatmap of metabolites with a significant Sex × Genotype interaction (two-way ANOVA).
   - Outputs: Metabolomics heatmaps.
4. **`4_4_boxplots.R`**
   - Creates boxplots for metabolites showing significant Sex × Genotype interactions.
   - Outputs: Boxplots for selected metabolites.
5. **`4_5_assembleplot.R`**
   - Assembles a cross-dataset figure combining mouse pathway/CAZyme and human metabolomics results.
   - Outputs: Combined multi-panel figure.

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
- **Human Cohort 2**: Table 1, compositional plots, beta diversity PCoA plots, LinDA volcano plots, pathway correlation heatmap, CAZyme plots, assembled figure.
- **Microbiome Analysis**: Compositional plots, alpha/beta diversity plots, TCAM scatterplots and loadings.
- **Pathway Analysis**: Correlation heatmaps, pathway boxplots, CAZyme plots.
- **Metabolomics Analysis**: PCA plots, Sex × Genotype interaction heatmaps and boxplots.
- **Statistical Results**: PERMANOVA results, LinDA results, LMM results.

## Requirements
### Environment
Dependencies are managed with [pixi](https://pixi.sh). To set up the environment:

```bash
bash setup_pixi.sh
```

This installs all R and Python dependencies defined in `pixi.toml`, plus additional Bioconductor packages (TreeSummarizedExperiment, mia, phyloseq, microbiome).

### R Packages (key)
- `tidyverse`, `ggplot2`, `ggpubr`, `ggrepel`, `ggthemes`, `ggsci`
- `ComplexHeatmap`, `circlize`
- `vegan`, `ape`, `MicrobiomeStat`
- `rio`, `tableone`, `gt`
- `Cairo`
- `rstatix`, `gridExtra`

### Python Packages (key)
- `numpy`, `pandas`, `scipy`
- `matplotlib`, `seaborn`
- `mprod` (for TCAM analysis)

## Code Availability
All scripts and workflows are available in this repository.
The metagenomics processing workflow is available at: [https://github.com/barbarahelena/metagenomicspipeline](https://github.com/barbarahelena/metagenomicspipeline).

## Data Availability
The raw metagenomics sequencing data is deposited in the European Nucleotide Archive (ENA) under accession number [PRJEB86491](https://www.ebi.ac.uk/ena/browser/view/PRJEB86491).
