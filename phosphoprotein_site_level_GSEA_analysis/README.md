# Phosphoprotein Site-Level GSEA Analysis

This folder contains the pipeline for phosphosite-level Gene Set Enrichment Analysis (GSEA) using PTMsigDB across CPTAC proteomics datasets.

## Overview

The analysis examines associations between COP9 Signalosome (CSN) subunit expression and phosphosite-level gene sets, stratified by TP53 mutation status.

## Files

### Main Script
- **`phosphosite_GSEA_main.R`** - Main analysis pipeline that orchestrates all steps

### Utility Functions (R/)
| File | Description |
|------|-------------|
| `utils_helpers.R` | General helper functions (z-score, normalization, logging) |
| `utils_geneset_loading.R` | PTMsigDB gene set loading functions |
| `utils_data_loading.R` | Data loading and phosphosite matrix functions |
| `utils_csn_score.R` | CSN composite score calculation |
| `utils_gsea.R` | fgsea wrappers and result processing |
| `utils_limma_analysis.R` | limma alignment and analysis functions |
| `utils_covariates.R` | Batch and covariate handling functions |
| `utils_tp53.R` | TP53 mutation status functions |
| `utils_meta_analysis.R` | Stouffer meta-analysis across datasets |

## Data Download

Before running, download the following CPTAC datasets from [cBioPortal](https://www.cbioportal.org/) and place them in the CSN_CPTAC folder:

- `brca_cptac_2020`
- `luad_cptac_2020`
- `lusc_cptac_2021`
- `ucec_cptac_2020`
- `coad_cptac_2019`
- `gbm_cptac_2021`
- `paad_cptac_2021`

Additionally, download the PTMsigDB database file (`data_PTMsigDB_all_sites_v2.0.0.xlsx`) from https://proteogenomics.shinyapps.io/ptmsigdb/ptmsigdb.Rmd and place it in the CSN_CPTAC root folder.

## Usage

```r
# Set working directory to CSN_CPTAC root
setwd("/path/to/CSN_CPTAC")

# Run the main analysis
source("phosphoprotein_site_level_GSEA_analysis/phosphosite_GSEA_main.R")
```

## Output Structure

```
phosphoproteomic_site_level_GSEA/
├── PTMsigDB/BatchAdj/{dataset}/{stratum}/{subunit}/
│   ├── GSEA_limma_t_cont.csv
│   └── GSEA_limma_interaction.csv
└── meta_fdr/
    ├── {stratum}/BatchAdj/{subunit}/PTMsigDB/
    │   └── GSEA_limma_t_cont_meta_fdr.csv
    └── summary/{stratum}/BatchAdj/PTMsigDB/
        └── Summary_PTMsigDB_GSEA_limma_t_cont_meta_fdr_*.csv
```

## Requirements

- R >= 4.0
- Required packages: tidyverse, data.table, limma, fgsea, imputeLCMD, openxlsx, yaml

## Datasets

Analyzes the following CPTAC datasets:
- BRCA (2020), LUAD (2020), UCEC (2020), COAD (2019), GBM (2021), PAAD (2021)

> **Note:** LUSC (2021) is excluded from this analysis because its phosphoproteome data format is incompatible with the other datasets (different column structure in the phosphoprotein quantification file).
