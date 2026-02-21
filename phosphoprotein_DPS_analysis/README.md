# Phosphoprotein DPS Analysis

## Directory Structure

```
phosphoprotein_DPS_analysis/
├── DPS_main.R              # Main workflow script
├── README.md               # This file
└── R/                      # Utility functions
    ├── 00_packages.R       # Package loading
    ├── 01_utils_general.R  # General utilities
    ├── 02_data_loading.R   # Data loading functions
    ├── 03_batch_processing.R # Batch effect handling
    ├── 04_covariates.R     # Covariate processing
    ├── 05_tp53_status.R    # TP53 mutation status
    ├── 06_csn_score.R      # CSN complex score calculation
    ├── 07_dps_analysis.R   # Core DPS analysis functions
    ├── 08_meta_analysis.R  # Stouffer meta-analysis
    └── 09_summary_output.R # Summary table generation
```

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

>  **Note:** LUSC (2021) is excluded from this analysis because its phosphoproteome data format is incompatible with the other datasets (different column structure in the phosphoprotein quantification file).

## Usage

1. **Set working directory** to the CSN_CPTAC project root
2. **Run the main script**:
   ```r
   setwd("/path/to/CSN_CPTAC")
   source("phosphoprotein_DPS_analysis/DPS_main.R")
   ```

## Analysis Overview

### Predictors Analyzed
- **CSN Subunits**: GPS1, COPS2-9, COPS7A, COPS7B
- **CSN_SCORE**: PC1 of CSN subunit expression (complex-level)
- **RESIDUAL_\<SU\>**: Subunit effect adjusted for CSN_SCORE 

### Stratification
- **ALL**: All samples
- **TP53_mutant**: TP53-mutated samples only
- **TP53_wild_type**: TP53-wild-type samples only

### Statistical Methods
- **Limma**: Linear modeling with empirical Bayes moderation
- **Covariates**: Batch, sex, age 
- **Meta-analysis**: Stouffer's method across datasets
- **Multiple testing**: Benjamini-Hochberg FDR correction

## Output Structure

```
phosphoproteomic_DPS/
├── <dataset>/
│   ├── ALL/DPS/BatchAdj/
│   │   ├── SU__GPS1/DPS_results.csv
│   │   ├── SU__CSN_SCORE/DPS_results.csv
│   │   ├── SU__RESIDUAL_GPS1/DPS_results.csv
│   │   └── DPS_summary_wide.csv
│   ├── TP53_mutant/DPS/BatchAdj/...
│   └── TP53_wild_type/DPS/BatchAdj/...
└── phospho_DPS_meta/
    ├── ALL/BatchAdj/<subunit>/DPS_meta_stouffer.csv
    └── summary/<stratum>/BatchAdj/DPS/...
```

## Key Output Files
| File | Description |
|------|-------------|
| `DPS_results.csv` | Per-predictor limma results (logFC, p-value, adj.P.Val) |
| `DPS_summary_wide.csv` | Wide-format summary of all predictors |
| `DPS_meta_stouffer.csv` | Meta-analysis results across datasets |

## Dependencies
- tidyverse, data.table, limma, yaml, vroom, readr
- imputeLCMD, openxlsx, ComplexHeatmap, cowplot, matrixStats

