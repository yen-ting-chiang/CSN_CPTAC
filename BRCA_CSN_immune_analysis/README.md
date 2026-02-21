# BRCA CSN Immune Analysis

CSN subunits vs. immune/stromal score correlation analysis for BRCA CPTAC dataset.

## Data Download

Before running, download the following CPTAC dataset from [cBioPortal](https://www.cbioportal.org/) and place it in the CSN_CPTAC folder:

- `brca_cptac_2020`

## Quick Start

```r
# run:
source("BRCA_CSN_immune_analysis/BRCA_immune_main.R")
```

## Project Structure

```
BRCA_CSN_immune_analysis/
├── BRCA_immune_main.R          # Main workflow script
├── README.md                   # This file
└── R/
    ├── 01_config.R             # Configuration and parameters
    ├── 02_utils_logging.R      # Logging utilities
    ├── 03_utils_io.R           # File I/O utilities
    ├── 04_utils_csn_score.R    # CSN score calculation
    ├── 05_utils_batch.R        # Batch effect handling
    ├── 06_utils_tp53.R         # TP53 mutation status
    ├── 07_utils_covariates.R   # Purity, sex, age covariates
    ├── 08_utils_correlation.R  # Correlation helpers
    ├── 09_utils_plotting.R     # Plotting functions
    └── 10_utils_audit.R        # Audit utilities
```

## Analysis Overview

**Predictors:**
- Individual CSN subunits (GPS1, COPS2-9)
- CSN_SCORE (PC1 of subunits)
- RESIDUAL_\<subunit\> (regressed on CSN_SCORE)

**Outcomes (from data_clinical_patient.txt):**
- ESTIMATE_IMMUNE_SCORE
- ESTIMATE_STROMAL_SCORE
- XCELL_IMMUNE_SCORE
- XCELL_STROMAL_SCORE
- CIBERSORT_ABSOLUTE_SCORE


**Covariates:**
- Sex, Age, Tumor purity
- Batch (BatchAdj version only)

**Stratification:**
- ALL, TP53_mutant, TP53_wild_type

## Output

```
BRCA_immune_analysis/
├── ALL/
│   ├── RAW/
│   │   ├── immune_vs_predictor_correlations.csv
│   │   ├── immune_vs_predictor_interaction.csv
│   │   ├── immune_vs_predictor_heatmap_PEARSON.tiff
│   │   └── plots/
│   └── BatchAdj/
│       └── ...
├── TP53_mutant/
│   └── ...
└── TP53_wild_type/
    └── ...
```

## Configuration

Edit `BRCA_immune_main.R` to customize:

### Scatter Plot Colors
```r
BRCA_PLOT_COLORSET <- "CELL_BLUE"  # Options:
# CELL_BLUE, CELL_TEAL, CELL_ORANGE, CELL_PURPLE
# CELL_BRIGHT_BLUE, CELL_AQUA, CELL_MINT
# CELL_MAGENTA, CELL_SCARLET, CELL_GOLD
```

### Heatmap Colors
```r
BRCA_HM_COLORSET <- "BLUE_RED_CRIMSON"  # Options:
# GSEA_DEFAULT, BLUE_RED, BLUE_ORANGE, GREEN_MAGENTA
# BLUE_RED_DEEP, BLUE_RED_BRIGHT, BLUE_RED_LIGHT
# BLUE_RED_TWILIGHT, BLUE_RED_CRIMSON
```

### Plot Filters
```r
BRCA_PLOT_FILTERS <- list(
    stratum   = c("ALL"),
    version   = c("RAW", "BatchAdj"),
    score     = c("CIBERSORT_ABSOLUTE_SCORE"),
    predictor = c("CSN_SCORE", "COPS7A", "COPS7B"),
    method    = c("pearson")
)
```

## Requirements

R packages:
- tidyverse, data.table, janitor, glue
- SummarizedExperiment, MultiAssayExperiment, S4Vectors
- limma, imputeLCMD, msigdbr
- openxlsx, ComplexHeatmap, cowplot
- matrixStats, yaml, ragg, scales

