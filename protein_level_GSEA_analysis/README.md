# CSN CPTAC Proteomic GSEA Analysis

## Data Download

Before running, download the following CPTAC datasets from [cBioPortal](https://www.cbioportal.org/) and place them in the CSN_CPTAC folder:

- `brca_cptac_2020`
- `luad_cptac_2020`
- `lusc_cptac_2021`
- `ucec_cptac_2020`
- `coad_cptac_2019`
- `gbm_cptac_2021`
- `paad_cptac_2021`

## Directory Structure

```
protein_level_GSEA_analysis/
├── R/                              # Modular function libraries
│   ├── 01_utility_functions.R      # Logging, naming, helper utilities
│   ├── 02_data_loading.R           # Matrix loading, clinical data extraction
│   ├── 03_csn_score.R              # CSN complex score calculation
│   ├── 04_batch_detection.R        # Batch effect identification
│   ├── 05_covariate_selection.R    # Covariate selection
│   ├── 06_purity_extraction.R      # Tumor purity data extraction
│   ├── 07_msigdb_utils.R           # MSigDB gene set utilities
│   ├── 08_meta_analysis.R          # Pan-cancer meta-analysis
│   ├── 09_meta_analysis_summary.R  # Meta-analysis summary across subunits
│   └── 10_pairwise_correlation.R   # CSN subunit correlation analysis
├── config/
│   └── analysis_parameters.yaml    # Analysis configuration
├── main_analysis.R                 # Main execution script
└── README.md                       # This file
```

## Quick Start

```r
# Set working directory to CSN_CPTAC root
setwd("/path/to/CSN_CPTAC")

# Run the main analysis
source("protein_level_GSEA_analysis/main_analysis.R")
```

## Configuration

Edit `config/analysis_parameters.yaml` to customize analysis parameters.

## Output

Results are saved in:
- `proteomic_GSEA/` - GSEA results by collection
- `csn_gsea_pan_summary_TP53/meta_fdr/` - Meta-analysis results
- `csn_gsea_pan_summary_TP53/meta_fdr/summary/` - Meta-analysis summary
- `CSN_subunits_correlation_coefficient/` - Pairwise correlations
- `run_info/` - Audit and metadata

## Requirements

R packages: tidyverse, data.table, limma, fgsea, msigdbr, imputeLCMD, openxlsx, yaml, janitor, glue, cowplot
