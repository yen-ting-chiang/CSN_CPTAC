# Phosphoprotein Gene-Level GSEA Analysis

## Project Structure

```
phosphoprotein_gene_level_GSEA_analysis/
├── phosphoprotein_gene_level_main.R    # Main analysis script
├── R/                                   # Utility modules
│   ├── 01_config.R                     # Configuration and parameters
│   ├── 02_utils_logging.R              # Logging utilities
│   ├── 03_utils_io.R                   # Data I/O and imputation
│   ├── 04_utils_geneset.R              # Gene set building (MSigDB)
│   ├── 05_utils_csn_score.R            # CSN score calculation
│   ├── 06_utils_gsea.R                 # GSEA analysis (fgsea)
│   ├── 07_utils_limma.R                # limma regression analysis
│   ├── 08_utils_tp53.R                 # Tp53 stratification
│   ├── 09_utils_batch.R                # Batch correction handling
│   ├── 10_utils_meta_analysis.R        # Cross-dataset meta-analysis
│   ├── 11_utils_covariates.R           # Covariate handling (purity, sex, age)
│   └── 12_utils_audit.R                # Dataset auditing utilities
└── README.md                           # This file
```

## Data Download

Before running, download the following CPTAC datasets from [cBioPortal](https://www.cbioportal.org/) and place them in the `CSN_CPTAC` folder:

- `brca_cptac_2020`
- `luad_cptac_2020`
- `lusc_cptac_2021`
- `ucec_cptac_2020`
- `coad_cptac_2019`
- `gbm_cptac_2021`
- `paad_cptac_2021`

Each dataset folder should contain the `data_phosphoprotein_quantification.txt` file required for this analysis.

## Usage

```r
# Set working directory to the CSN_CPTAC project root
setwd("/path/to/CSN_CPTAC")

# Run the complete analysis
source("phosphoprotein_gene_level_GSEA_analysis/phosphoprotein_gene_level_main.R")
```

## Output Structure

```
phosphoproteomic_gene_level_GSEA/
├── {GROUP}/                            # Gene set group (H, C2, etc.)
│   ├── RAW/{dataset}/{stratum}/{subunit}/
│   │   ├── GSEA_limma_t_cont.csv
│   │   └── GSEA_limma_interaction.csv
│   └── BatchAdj/{dataset}/{stratum}/{subunit}/
│       └── ...
└── phospho_csn_gsea_pan_summary_TP53/
    └── meta_fdr/
        ├── {stratum}/{version}/{subunit}/{group}/
        │   └── GSEA_limma_t_cont_meta_fdr.csv
        └── summary/
            └── {stratum}/{version}/{group}/{stat_tag}/
                ├── Summary_*_ALL.csv
                ├── Summary_*_padjLT0.05.csv
                └── Summary_*.xlsx
```


## Dependencies

Required R packages:
- `tidyverse`, `data.table`, `readr`, `dplyr`
- `limma` (Bioconductor)
- `fgsea` (Bioconductor)
- `msigdbr` (MSigDB gene sets)
- `imputeLCMD` (imputation)
- `openxlsx` (Excel output)

## Analysis Pipeline

1. **Load Data**: Phosphoprotein quantification matrix (protein-adjusted)
2. **Imputation**: MinProb imputation for missing values
3. **Covariate Extraction**: Purity, sex, age from clinical data
4. **CSN Score**: PCA-based complex activity score
5. **limma Analysis**: Continuous and interaction models
6. **GSEA**: fgseaMultilevel on limma t-statistics
7. **Per-Dataset Summary**: Merge results across subunits
8. **Meta-Analysis**: Stouffer's Z across datasets with BH FDR
