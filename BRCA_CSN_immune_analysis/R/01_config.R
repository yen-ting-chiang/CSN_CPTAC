# =============================================================================
# 01_config.R
# Configuration and Global Parameters for BRCA Immune Analysis
# =============================================================================

# -----------------------------------------------------------------------------
# Random Seed
# -----------------------------------------------------------------------------
set.seed(1234)

# -----------------------------------------------------------------------------
# CSN Complex Subunits
# -----------------------------------------------------------------------------
csn_subunits <- c(
    "GPS1", "COPS2", "COPS3", "COPS4", "COPS5",
    "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9"
)

# -----------------------------------------------------------------------------
# Batch Cleaning Settings
# -----------------------------------------------------------------------------
# BATCH_PIPE_POLICY: "NA" or "b_small"
#   - "NA": set values containing '|' to NA
#   - "b_small": merge values containing '|' to b_small
BATCH_PIPE_POLICY <- "NA"

# BATCH_MIN_PER_LEVEL: Levels with fewer samples will be merged to b_small
# (to avoid non-estimable coefficients in linear models)
BATCH_MIN_PER_LEVEL <- 2

# -----------------------------------------------------------------------------
# TP53 Variant Classification Constants
# -----------------------------------------------------------------------------
# Only count protein-altering variants as TP53-mutant
# Others (Silent, UTR, Intron, IGR, RNA, lincRNA, Flank...) are treated as WT
TP53_KEEP_CLASSES <- c(
    "MISSENSE_MUTATION", "NONSENSE_MUTATION",
    "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
    "IN_FRAME_DEL", "IN_FRAME_INS",
    "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
)

# -----------------------------------------------------------------------------
# Age Missingness Handling
# -----------------------------------------------------------------------------
# If TRUE, use missing-indicator approach for age
# If FALSE (default), keep NA only
USE_AGE_MISSING_INDICATOR <- FALSE

# -----------------------------------------------------------------------------
# Global Helper Functions
# -----------------------------------------------------------------------------

#' Get option value with default fallback
#' @param nm Option name
#' @param default Default value if option not found
#' @return Option value or default
opt <- function(nm, default) {
    if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
}

#' Null coalescing operator
#' @param a Primary value
#' @param b Fallback value if a is NULL
#' @return a if not NULL, otherwise b
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Create safe filesystem name
#' @param s Input string
#' @return Sanitized string safe for filesystem use
safe_fs_name <- function(s) {
    s <- gsub('[<>:"/\\\\|?*]', "_", s)
    s <- gsub("\\s+", "_", s)
    s
}
