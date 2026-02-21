# mRNA GSEA Analysis Pipeline

CSN Subunits → mRNA-based GSEA Analysis for CPTAC Datasets

## Overview

This pipeline performs Gene Set Enrichment Analysis (GSEA) using protein-level CSN subunit expression as predictors for mRNA-level pathway enrichment across 7 CPTAC cancer datasets.

## Data Download

Before running, download the following CPTAC datasets from [cBioPortal](https://www.cbioportal.org/) and place them in the CSN_CPTAC folder:

- `brca_cptac_2020`
- `luad_cptac_2020`
- `lusc_cptac_2021`
- `ucec_cptac_2020`
- `coad_cptac_2019`
- `gbm_cptac_2021`
- `paad_cptac_2021`

## Project Structure

```
mRNA_GSEA_analysis/
├── main_mRNA_GSEA_pipeline.R    # Main orchestrating script
├── README.md                     # This file
└── R/                            # Utility modules
    ├── 01_config.R               # Configuration and parameters
    ├── 02_utils_logging.R        # Logging functions
    ├── 03_utils_io.R             # File I/O utilities
    ├── 04_utils_geneset.R        # MSigDB gene set preparation
    ├── 05_utils_csn_score.R      # CSN score calculation 
    ├── 06_utils_gsea.R           # GSEA execution functions
    ├── 07_utils_limma.R          # Limma analysis utilities
    ├── 08_utils_tp53.R           # TP53 stratification
    ├── 09_utils_batch.R          # Batch effect handling
    ├── 10_utils_meta_analysis.R  # Meta-analysis (Stouffer's z)
    ├── 11_utils_covariates.R     # Covariate auditing
    ├── 12_utils_correlation.R    # CSN subunit correlations
    ├── 13_utils_tpm.R            # TPM calculations
    ├── 14_utils_covariates_ext.R # Covariate extraction
    ├── 15_core_analysis_functions.R # Core GSEA analysis
    └── 16_utils_meta_fdr_summary.R  # Meta-FDR summary
```


## Analysis Strata

- **ALL**: All samples
- **TP53_mutant**: Samples with TP53 mutations
- **TP53_wild_type**: Samples without TP53 mutations

## Output Structure

```
mRNA_GSEA/
├── <GENESET_GROUP>/CovarAdj/<DATASET>/<STRATUM>/<PREDICTOR>/
│   ├── GSEA_limma_t_cont.csv
│   └── GSEA_limma_interaction.csv
├── csn_gsea_pan_summary_TP53/meta_fdr/
│   ├── <STRATUM>/<RAW|CovarAdj>/<SUBUNIT>/<GROUP>/
│   └── summary/
├── RNA_CSN_subunits_correlation_coefficient/
└── TPM/
    ├── CSN_subunits_TPM_median_TPM__RAW.csv
    └── CSN_subunits_TPM_median_TPM__AdjCovars.csv
```

## Usage

```r
# Set working directory to project root
setwd()

# Run the complete pipeline
source("mRNA_GSEA_analysis/main_mRNA_GSEA_pipeline.R")
```

## Dependencies

- tidyverse, data.table, janitor, glue
- SummarizedExperiment, MultiAssayExperiment, S4Vectors
- limma, imputeLCMD, msigdbr, fgsea
- openxlsx, ComplexHeatmap, cowplot, matrixStats
- yaml, readr
