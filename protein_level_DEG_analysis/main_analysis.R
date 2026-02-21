## ============================================================================
## Main Analysis Pipeline - Proteomic DEG Analysis
## ============================================================================

# --- Working Directory Setup ---
# NOTE: This script should be run from the CSN_CPTAC project root directory.
# Please refer to README.md for data download instructions and setup guide.
#
# Example:
#   setwd("/path/to/CSN_CPTAC")
#   source("protein_level_DEG_analysis/main_analysis.R")

message("====================================================================")
message("CSN CPTAC Proteomic DEG Analysis")
message("====================================================================\n")

# Verify we are in the correct directory
if (!file.exists("protein_level_DEG_analysis/main_analysis.R")) {
    stop(
        "ERROR: Script must be run from CSN_CPTAC project root directory.\n",
        "Current directory: ", getwd(), "\n",
        "Please setwd() to the CSN_CPTAC root directory."
    )
}

message("Working directory: ", getwd(), "\n")

# Load All Modules --------------------------------------------------------
cat("[Main] Loading configuration and utility modules...\n")

source("protein_level_DEG_analysis/R/config.R")
source("protein_level_DEG_analysis/R/utils_data.R")
source("protein_level_DEG_analysis/R/utils_tp53.R")
source("protein_level_DEG_analysis/R/utils_batch.R")
source("protein_level_DEG_analysis/R/utils_csn.R")
source("protein_level_DEG_analysis/R/utils_stats.R")
source("protein_level_DEG_analysis/R/analysis_deg.R")
source("protein_level_DEG_analysis/R/analysis_stratum.R")
source("protein_level_DEG_analysis/R/analysis_meta.R") # Meta-analysis functions
source("protein_level_DEG_analysis/R/analysis_summary.R") # Summary aggregation

log_msg("========================================")
log_msg("Proteomic DEG Analysis Pipeline")
log_msg("CSN Subunits -> Differential Expression")
log_msg("========================================")

# Gene Set Preparation ----------------------------------------------------
log_msg("\n[1/5] Preparing gene sets")
log_msg("Gene set collections to use: %s", paste(GENESET_GROUPS_TO_RUN, collapse = ", "))

# Helper functions for gene set loading
.fetch_msig <- function(cat) {
    try(msigdbr::msigdbr(species = "Homo sapiens", category = toupper(cat)), silent = TRUE)
}

.parse_group_token <- function(tok) {
    t <- toupper(trimws(tok))
    sp <- unlist(strsplit(t, ":", fixed = TRUE))
    list(cat = sp[1], sub = if (length(sp) >= 2) paste(sp[-1], collapse = ":") else NULL, raw = t)
}

.get_subcat_col <- function(df) {
    cand <- c("gs_subcollection", "gs_subcat", "gs_subcategory", "subcategory", "sub_category")
    hit <- cand[cand %in% names(df)]
    if (length(hit)) hit[1] else NULL
}

.make_group_label <- function(cat, df, want_sub) {
    subcol <- .get_subcat_col(df)
    if (!is.null(want_sub)) {
        if (!is.null(subcol)) {
            uniq <- unique(df[[subcol]])
            hit <- uniq[grepl(want_sub, uniq, ignore.case = TRUE)]
            if (length(hit) >= 1) {
                return(sprintf("%s__%s", cat, hit[1]))
            }
        }
        return(sprintf("%s__%s", cat, safe_fs_name(want_sub)))
    }
    return(cat)
}

# Build gene sets
genesets_by_group <- list()

for (tok in GENESET_GROUPS_TO_RUN) {
    pp <- .parse_group_token(tok)
    cat <- pp$cat
    df <- suppressMessages(.fetch_msig(cat))

    if (inherits(df, "try-error") || is.null(df) || !nrow(df)) {
        log_msg("  Note: msigdbr cannot provide category %s; skipping %s", cat, pp$raw)
        next
    }

    if (!is.null(pp$sub)) {
        subcol <- .get_subcat_col(df)
        subtxt <- if (!is.null(subcol)) as.character(df[[subcol]]) else ""
        esc <- function(s) gsub("([\\^\\$\\.|\\(\\)\\[\\]\\{\\}\\*\\+\\?\\\\])", "\\\\\\1", s)
        pattern <- paste0("^", esc(pp$sub), "($|:)")
        keep <- grepl(pattern, subtxt, ignore.case = TRUE)
        df <- df[keep, , drop = FALSE]

        if (!nrow(df)) {
            log_msg("  Note: %s cannot find any subset of %s; skipping", cat, pp$sub)
            next
        }
        grp <- .make_group_label(cat, df, pp$sub)
    } else {
        grp <- cat
    }

    genesets_by_group[[grp]] <- lapply(split(df$gene_symbol, df$gs_name), unique)
    log_msg("  Got %s: %d sets", grp, length(genesets_by_group[[grp]]))
}

if (!length(genesets_by_group)) {
    stop("No available gene-set. Please check GENESET_GROUPS_TO_RUN setting.")
}

log_msg("Final gene-set groups: %s", paste(names(genesets_by_group), collapse = ", "))

# Dataset Preparation -----------------------------------------------------
log_msg("\n[2/5] Preparing datasets")

dataset_dirs <- get_dataset_dirs(dataset_ids, getwd())
log_msg("Datasets to analyze: %d", length(dataset_dirs))
for (ds_name in names(dataset_dirs)) {
    log_msg("  - %s", ds_name)
}

# Optional: Run TP53 audit for all datasets
log_msg("\nRunning TP53 status audit...")
summarize_tp53_counts_all_datasets(dataset_dirs)

# Main Analysis Loop ------------------------------------------------------
log_msg("\n[3/5] Running main analysis")
log_msg("========================================")

for (ds_id in names(dataset_dirs)) {
    ds_dir <- dataset_dirs[ds_id]

    log_msg("\n========== Dataset: %s ==========", ds_id)

    # Load full matrix
    mat0_full <- load_matrix_from_dataset_dir(ds_dir)
    log_msg("Loaded matrix: %d genes × %d samples", nrow(mat0_full), ncol(mat0_full))

    # Get TP53 status for stratification
    tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full))

    # Prepare strata
    strata_info <- list(
        ALL = list(
            name = "ALL",
            samples = colnames(mat0_full),
            description = "All samples"
        ),
        TP53_mutant = list(
            name = "TP53_mutant",
            samples = names(tp53_status)[tp53_status == "TP53_mutant"],
            description = "TP53 mutant samples"
        ),
        TP53_wild_type = list(
            name = "TP53_wild_type",
            samples = names(tp53_status)[tp53_status == "TP53_wild_type"],
            description = "TP53 wild-type samples"
        )
    )

    # Loop through strata
    for (stratum_id in names(strata_info)) {
        stratum <- strata_info[[stratum_id]]

        log_msg("\n--- Stratum: %s ---", stratum$name)
        log_msg("Description: %s", stratum$description)
        log_msg("Sample count: %d", length(stratum$samples))

        # Skip if too few samples
        if (length(stratum$samples) < 4) {
            log_msg("  [Skip] Insufficient samples (< 4)")
            next
        }

        # Define output root for this stratum (MATCH original script structure)
        out_root <- file.path(OUTPUT_PREFIX, ds_id, "csn_gsea_results_TP53", stratum$name)

        # Run complete stratum analysis
        run_one_stratum(
            ds_id = ds_id,
            ds_dir = ds_dir,
            mat0_full = mat0_full,
            sample_keep = stratum$samples,
            out_root = out_root,
            genesets_by_group = genesets_by_group
        )

        # Summarize DEG across predictors for this stratum
        log_msg("Creating cross-predictor summary for %s...", stratum$name)
        summarize_deg_across_predictors(out_root)
    }

    log_msg("\n========== Completed: %s ==========", ds_id)
}

# Meta-Analysis Across Datasets ------------------------------------------
log_msg("\n[4/5] Running meta-analysis across datasets")
log_msg("========================================")

# Run meta-FDR using Stouffer method (matches original script structure)
if (length(dataset_dirs) >= 2) {
    meta_fdr_stouffer_DEG(
        dataset_dirs = dataset_dirs,
        strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
        versions = c("BatchAdj"),
        output_prefix = OUTPUT_PREFIX
    )
    log_msg("\n[Meta-DEG] Meta-analysis complete")
    log_msg("Results saved to: %s/csn_deg_pan_summary_TP53/meta_fdr/", OUTPUT_PREFIX)

    # Summarize meta-FDR results across subunits
    log_msg("\n[Meta-Summary] Aggregating results across subunits...")
    summarize_meta_fdr_across_subunits_DEG(
        meta_root = file.path(OUTPUT_PREFIX, "csn_deg_pan_summary_TP53/meta_fdr"),
        strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
        versions = c("BatchAdj")
    )
    log_msg("[Meta-Summary] Summary tables created")
} else {
    log_msg("Skipping meta-analysis (requires >= 2 datasets)")
}

# Post-Analysis Summary ---------------------------------------------------
log_msg("\n[5/5] Creating final summaries")

log_msg("\nAll datasets processed!")
log_msg("\nOutput structure:")
log_msg("  %s/", OUTPUT_PREFIX)
log_msg("    ├── {dataset}/")
log_msg("    │   ├── ALL/DEG/BatchAdj/")
log_msg("    │   │   ├── CSN_SCORE/DEG_limma_cont.csv")
log_msg("    │   │   ├── GPS1/DEG_limma_cont.csv")
log_msg("    │   │   ├── GPS1_adj_CSN/DEG_limma_cont.csv")
log_msg("    │   │   └── Summary_DEG_wide.csv")
log_msg("    │   ├── TP53_mutant/DEG/BatchAdj/...")
log_msg("    │   └── TP53_wild_type/DEG/BatchAdj/...")
log_msg("    └── ...")

log_msg("\n========================================")
log_msg("Analysis Pipeline Complete!")
log_msg("========================================")

# Generate summary report
log_msg("\nGenerating analysis summary report...")

summary_lines <- c(
    "# Proteomic DEG Analysis - Execution Summary",
    "",
    sprintf("**Analysis Date**: %s", Sys.Date()),
    sprintf("**Output Directory**: %s", OUTPUT_PREFIX),
    "",
    "## Datasets Analyzed",
    ""
)

for (ds_id in names(dataset_dirs)) {
    summary_lines <- c(
        summary_lines,
        sprintf("### %s", ds_id)
    )

    for (stratum_name in c("ALL", "TP53_mutant", "TP53_wild_type")) {
        out_root <- file.path(OUTPUT_PREFIX, ds_id, stratum_name)
        deg_dir <- file.path(out_root, "DEG", "BatchAdj")

        if (dir.exists(deg_dir)) {
            predictors <- list.dirs(deg_dir, full.names = FALSE, recursive = FALSE)
            n_predictors <- length(predictors)

            summary_lines <- c(
                summary_lines,
                sprintf("- **%s**: %d predictors analyzed", stratum_name, n_predictors)
            )

            # Check for summary file
            summary_file <- file.path(deg_dir, "Summary_DEG_wide.csv")
            if (file.exists(summary_file)) {
                summary_df <- data.table::fread(summary_file, nrows = 1)
                summary_lines <- c(
                    summary_lines,
                    sprintf(
                        "  - Summary file: %d genes",
                        nrow(data.table::fread(summary_file))
                    )
                )
            }
        }
    }

    summary_lines <- c(summary_lines, "")
}

summary_lines <- c(
    summary_lines,
    "## Analysis Parameters",
    "",
    sprintf("- Minimum fraction complete: %.2f", min_frac_complete),
    sprintf("- Minimum samples per group: %d", min_per_group),
    sprintf("- CSN subunits: %s", paste(csn_subunits, collapse = ", ")),
    sprintf("- Gene set collections: %s", paste(names(genesets_by_group), collapse = ", ")),
    "",
    "## Output Files",
    "",
    "For each dataset and stratum:",
    "- `DEG/BatchAdj/{predictor}/DEG_limma_cont.csv`: Individual DEG results",
    "- `DEG/BatchAdj/Summary_DEG_wide.csv`: Cross-predictor summary",
    ""
)

# Write summary report
summary_file <- file.path(OUTPUT_PREFIX, "ANALYSIS_SUMMARY.md")
writeLines(summary_lines, summary_file)
log_msg("Summary report written to: %s", summary_file)

log_msg("\n All analyses complete!")
log_msg("Check %s for results", OUTPUT_PREFIX)
