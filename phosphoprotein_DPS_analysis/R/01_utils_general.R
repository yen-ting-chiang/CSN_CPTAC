## =========================================================
## 01_utils_general.R
## General Utility Functions for DPS Analysis
## =========================================================

## ===== Global helper: opt (retrieve options with default) =====
opt <- function(nm, default) {
    if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
}

## ===== Null coalescing operator =====
`%||%` <- function(a, b) if (!is.null(a)) a else b

## ===== Timestamped logging function =====
log_msg <- function(text, ..., .envir = parent.frame()) {
    ts <- format(Sys.time(), "%H:%M:%S")
    msg <- tryCatch(
        {
            if (grepl("\\{[^}]+\\}", text)) {
                # Has { } -> use glue
                glue::glue(text, ..., .envir = .envir)
            } else if (grepl("%", text)) {
                # Has % -> use sprintf
                do.call(sprintf, c(list(fmt = text), list(...)))
            } else {
                # Plain text; if there are extra arguments, treat as sprintf format
                if (nargs() > 1L) do.call(sprintf, c(list(fmt = text), list(...))) else text
            }
        },
        error = function(e) text
    )
    cat(sprintf("[%s] %s\n", ts, msg))
}

## ===== Filesystem-safe name conversion =====
safe_fs_name <- function(s) {
    s <- gsub('[<>:"/\\\\|?*]', "_", s)
    s <- gsub("\\s+", "_", s)
    s
}

## ===== Helper: coerce_covariates_safely =====
coerce_covariates_safely <- function(df) {
    df <- as.data.frame(df, check.names = FALSE)
    keep <- rep(TRUE, ncol(df))
    names(keep) <- colnames(df)

    for (cn in colnames(df)) {
        v <- df[[cn]]
        if (is.factor(v) || is.character(v) || is.logical(v)) {
            v <- factor(v)
            lv <- levels(droplevels(v[!is.na(v)]))
            if (length(lv) <= 1) { # Single-level factor -> drop to avoid contrasts error
                keep[cn] <- FALSE
                if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
            } else {
                df[[cn]] <- v
            }
        } else {
            df[[cn]] <- suppressWarnings(as.numeric(v)) # Keep numeric as numeric
            
        }
    }
    df <- df[, keep, drop = FALSE]
    df
}

## ===== Force single-thread execution =====
.force_serial_execution <- function() {
    # base R / data.table
    options(mc.cores = 1L)
    if ("package:data.table" %in% search() &&
        exists("setDTthreads", asNamespace("data.table"))) {
        data.table::setDTthreads(1L)
    }

    # BiocParallel
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
        BiocParallel::register(BiocParallel::SerialParam())
    }

    # foreach / doParallel
    if (requireNamespace("foreach", quietly = TRUE)) {
        foreach::registerDoSEQ()
    }

    # future
    if (requireNamespace("future", quietly = TRUE)) {
        future::plan(future::sequential)
    }
    Sys.setenv("R_FUTURE_FORK_ENABLE" = "FALSE")
}

## ===== Clean matrix for PCA: remove Inf->NA, impute with row median =====
.clean_for_pca <- function(X, min_samples = 10L, min_genes = 5L) {
    X <- as.matrix(X)
    X[!is.finite(X)] <- NA
    keep_rows <- rowSums(is.finite(X)) >= min_samples
    if (!any(keep_rows)) {
        return(NULL)
    }
    X <- X[keep_rows, , drop = FALSE]
    if (nrow(X) < min_genes) {
        return(NULL)
    }
    if (anyNA(X)) {
        med <- apply(X, 1, function(r) median(r[is.finite(r)], na.rm = TRUE))
        for (i in seq_len(nrow(X))) {
            xi <- X[i, ]
            xi[!is.finite(xi)] <- med[i]
            X[i, ] <- xi
        }
    }
    X
}

## ===== Impute and filter matrix =====
impute_and_filter <- function(mat, min_frac = 0.75) {
    keep <- rowMeans(!is.na(mat)) >= min_frac
    m <- mat[keep, , drop = FALSE]
    if (any(is.na(m))) {
        set.seed(1234)
        m <- imputeLCMD::impute.MinProb(m, q = 0.01)
    }
    m
}

## ===== Orthogonalize matrix to nuisance variables =====
orthogonalize_to <- function(mat, nuisance) {
    if (is.null(mat) || is.null(nuisance)) {
        return(mat)
    }
    Y <- as.matrix(mat)
    X <- as.matrix(nuisance)

    # 1) Add intercept
    X <- cbind("(Intercept)" = rep(1, nrow(X)), X)

    # 2) Handle non-finite values + remove zero-variance columns
    for (j in seq_len(ncol(X))) {
        v <- X[, j]
        if (any(!is.finite(v))) {
            mu <- mean(v[is.finite(v)], na.rm = TRUE)
            if (!is.finite(mu)) mu <- 0
            v[!is.finite(v)] <- mu
            X[, j] <- v
        }
    }
    var_ok <- apply(X, 2, function(z) {
        z <- as.numeric(z)
        is.finite(var(z)) && var(z) > 0
    })
    if (!all(var_ok)) X <- X[, var_ok, drop = FALSE]

    beta <- tryCatch(qr.coef(qr(X), Y), error = function(e) NULL)
    if (is.null(beta)) {
        return(as.matrix(mat))
    }

    resid <- Y - X %*% beta
    rownames(resid) <- rownames(Y)
    colnames(resid) <- colnames(Y)
    return(resid)
}

## ===== Align limma inputs (matrix and design) =====
.align_limma_inputs <- function(M, design_df, rhs_terms, label = "limma") {
    stopifnot(!is.null(M), is.matrix(M), ncol(M) > 0, !is.null(colnames(M)))
    # First convert design table to data.frame, ensure rownames are sample IDs
    design_df <- as.data.frame(design_df, check.names = FALSE)
    stopifnot(!is.null(rownames(design_df)))
    samp_M <- as.character(colnames(M))
    samp_des <- as.character(rownames(design_df))

    # First remove NA samples from design table
    na_mask <- rep(FALSE, nrow(design_df))
    rn <- rownames(design_df)

    # predictor must exist
    if (!("predictor" %in% colnames(design_df))) {
        stop(sprintf("[%s] design_df missing predictor column", label))
    }

    # predictor not finite
    v <- design_df$predictor
    na_mask <- na_mask | is.na(v) | !is.finite(v)

    # Other numeric columns
    num_cols <- setdiff(colnames(design_df), c("predictor"))
    for (cn in num_cols) {
        v <- design_df[[cn]]
        if (is.numeric(v)) {
            na_mask <- na_mask | is.na(v) | !is.finite(v)
        } else {
            # Factor/string/logical: remove if NA
            na_mask <- na_mask | is.na(v)
        }
    }
    design_df <- design_df[!na_mask, , drop = FALSE]

    # Get intersection and sort
    common <- intersect(as.character(colnames(M)), as.character(rownames(design_df)))
    if (!length(common)) {
        return(NULL)
    }
    common <- sort(common)

    M2 <- M[, common, drop = FALSE]
    design_df2 <- design_df[common, , drop = FALSE]

    # Build design matrix
    rhs <- unique(rhs_terms)
    des2 <- stats::model.matrix(stats::reformulate(rhs), data = design_df2)

    # Final check again
    if (nrow(des2) != ncol(M2)) {
        msg <- sprintf(
            "[%s] nrow(des)=%d != ncol(M)=%d (still inconsistent after alignment)",
            label, nrow(des2), ncol(M2)
        )
        stop(msg)
    }
    list(M = M2, des = des2, design_df = design_df2)
}

## ===== Audit design alignment =====
audit_design_alignment <- function(tag, samples, mod_interest, mod_nuisance, out_dir = NULL) {
    safe_rn <- function(X) if (is.null(X)) rep(NA_character_, length(samples)) else rownames(as.matrix(X))
    df <- data.frame(
        pos = seq_along(samples),
        M_sample = samples,
        interest_rn = safe_rn(mod_interest),
        nuisance_rn = safe_rn(mod_nuisance),
        interest_ok = if (!is.null(mod_interest)) safe_rn(mod_interest) == samples else NA,
        nuisance_ok = if (!is.null(mod_nuisance)) safe_rn(mod_nuisance) == samples else NA,
        stringsAsFactors = FALSE
    )
    n_i <- sum(df$interest_ok, na.rm = TRUE)
    n_n <- sum(df$nuisance_ok, na.rm = TRUE)
    log_msg(
        "[align:%s] interest match = %d/%d, nuisance match = %d/%d, both = %d",
        tag, n_i, length(samples), n_n, length(samples), sum(df$interest_ok & df$nuisance_ok, na.rm = TRUE)
    )
    if (!is.null(out_dir)) {
        fn <- file.path(out_dir, sprintf("ALIGN_%s.csv", gsub("[^A-Za-z0-9]+", "_", tag)))
        utils::write.csv(df, fn, row.names = FALSE)
        log_msg("[align:%s] wrote %s", tag, fn)
    }
    # Additional: list first 6 rows for quick log viewing
    head_show <- utils::head(df[, c("pos", "M_sample", "interest_rn", "nuisance_rn", "interest_ok", "nuisance_ok")], 6)
    log_msg("[align:%s] head:\n%s", tag, utils::capture.output(print(head_show)) |> paste(collapse = "\n"))
    invisible(df)
}

## ===== Ensure matrix or NULL helper =====
.ensure_mat_or_null <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    m <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (is.null(m) || length(m) == 0) {
        return(NULL)
    }
    if (is.null(dim(m))) m <- matrix(m, ncol = 1)
    if (nrow(m) == 0 || ncol(m) == 0) {
        return(NULL)
    }
    m[!is.finite(m)] <- NA_real_
    m
}

## ===== Helper functions for alignment =====
.take_first <- function(x) sub("\\|.*$", "", as.character(x))

.sanitize_levels <- function(f) {
    if (is.null(f)) {
        return(NULL)
    }
    f <- droplevels(factor(f))
    lv <- levels(f)
    lv2 <- make.names(lv)
    lv2 <- paste0("b_", lv2) # Avoid starting with number
    levels(f) <- lv2
    f
}

.align_by_colnames <- function(vec_or_df, target_names) {
    if (is.null(vec_or_df)) {
        return(NULL)
    }
    if (is.vector(vec_or_df) || is.factor(vec_or_df)) {
        nm <- names(vec_or_df)
        if (!is.null(nm) && length(nm) == length(vec_or_df)) {
            return(vec_or_df[match(target_names, nm)])
        } else {
            return(vec_or_df)
        }
    } else {
        rn <- rownames(vec_or_df)
        if (!is.null(rn)) {
            return(vec_or_df[match(target_names, rn), , drop = FALSE])
        }
        return(vec_or_df)
    }
}

.mk_batch_factor <- function(batch, sample_order, min_count_collapse = 1L) {
    if (is.null(batch)) {
        return(NULL)
    }
    b <- .align_by_colnames(batch, sample_order)
    b <- .take_first(b)
    b[!nzchar(b) | is.na(b)] <- "unknown"
    b <- .sanitize_levels(b)
    # Optional: merge very small batches
    if (min_count_collapse > 0L) {
        tb <- table(b)
        small <- names(tb)[tb <= min_count_collapse]
        if (length(small)) {
            b[b %in% small] <- factor("b_small")
            b <- droplevels(factor(b))
        }
    }
    b
}
