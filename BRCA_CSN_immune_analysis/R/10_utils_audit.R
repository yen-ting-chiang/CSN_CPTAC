# =============================================================================
# 10_utils_audit.R
# Audit and Summary Utilities
# =============================================================================

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Format batch level sample sizes into string
#'
#' @param fac Factor with batch levels
#' @return Formatted string like "batch1=10; batch2=8"
.format_batch_sizes <- function(fac) {
    if (is.null(fac)) {
        return(NA_character_)
    }
    tb <- sort(table(fac), decreasing = TRUE)
    paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
}

#' Convert table to string representation
#'
#' @param tb Table object
#' @return Formatted string
.tbl_str <- function(tb) {
    paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
}
