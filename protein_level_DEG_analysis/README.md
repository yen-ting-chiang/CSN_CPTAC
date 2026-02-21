# Proteomic DEG Analysis

## Data Download

Before running, download the following CPTAC datasets from [cBioPortal](https://www.cbioportal.org/) and place them in the CSN_CPTAC folder:

- `brca_cptac_2020`
- `luad_cptac_2020`
- `lusc_cptac_2021`
- `ucec_cptac_2020`
- `coad_cptac_2019`
- `gbm_cptac_2021`
- `paad_cptac_2021`


###  1: Data Processing Modules
- **config.R** : Global parameters and settings
- **utils_data.R** : Data loading and filtering
- **utils_tp53.R** : TP53 mutation status determination

###  2: Advanced Processing Modules
- **utils_batch.R** : Batch effect detection and handling
- **utils_csn.R** : CSN complex score calculation 
- **utils_stats.R** : Statistical utilities and covariate processing

###  3: Core Analysis Modules
- **analysis_deg.R** : DEG analysis core functions
- **analysis_stratum.R** : Stratum-level workflow 

###  4: Integration & Documentation
- **main_analysis.R**: Complete pipeline orchestration
- **README.md**: This file 


##  Quick Start

### Running the Complete Analysis

```r
# 1. Set working directory to CSN_CPTAC project root
setwd("/path/to/CSN_CPTAC")

# 2. Run main analysis pipeline
source("protein_level_DEG_analysis/main_analysis.R")
```

The script will verify you are in the correct directory and automatically handle all analysis steps.


##  Output Structure

```
{OUTPUT_PREFIX}/
├── {dataset}/
│   ├── ALL/DEG/BatchAdj/
│   │   ├── GPS1/DEG_limma_cont.csv
│   │   ├── CSN_SCORE/DEG_limma_cont.csv
│   │   ├── GPS1_adj_CSN/DEG_limma_cont.csv        # Interaction model
│   │   └── Summary_DEG_wide.csv                   # Cross-predictor summary
│   ├── TP53_mutant/DEG/BatchAdj/...
│   └── TP53_wild_type/DEG/BatchAdj/...
└── ANALYSIS_SUMMARY.md                            # Execution summary
```

## Dependencies

Required R packages: tidyverse, data.table, limma, msigdbr, fgsea, imputeLCMD, openxlsx, ComplexHeatmap

