setwd("C:/Users/danny/Documents/R_project/CSN_CPTAC") ## YTC laptop
## setwd("C:/Users/cmuh/Documents/YenTing_document/CSN_CPTAC") ## lab computer
getwd()


suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(janitor)
  library(glue)
  library(SummarizedExperiment)
  library(MultiAssayExperiment)
  library(S4Vectors)
  library(limma)
  library(imputeLCMD)
  library(msigdbr)
  library(openxlsx)
  library(ComplexHeatmap)
  library(cowplot)
  library(matrixStats)
  library(yaml)
})

## ===== global helper: opt (must exist before any writer uses it) =====
if (!exists("opt", mode = "function")) {
  opt <- function(nm, default) {
    if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
  }
}


## ===== helper: coerce_covariates_safely (must be global, defined before use) =====
if (!exists("coerce_covariates_safely", mode = "function")) {
  coerce_covariates_safely <- function(df) {
    df <- as.data.frame(df, check.names = FALSE)
    keep <- rep(TRUE, ncol(df))
    names(keep) <- colnames(df)

    for (cn in colnames(df)) {
      v <- df[[cn]]
      if (is.factor(v) || is.character(v) || is.logical(v)) {
        v <- factor(v)
        lv <- levels(droplevels(v[!is.na(v)]))
        if (length(lv) <= 1) { # single-level factor -> drop to avoid contrasts error
          keep[cn] <- FALSE
          if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
        } else {
          df[[cn]] <- v
        }
      } else {
        df[[cn]] <- suppressWarnings(as.numeric(v)) # keep numeric as numeric
      }
    }
    df <- df[, keep, drop = FALSE]
    df
  }
}


## ==== Force single-threaded, close all parallel backends (across common frameworks) ====
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
  if (requireNamespace("doParallel", quietly = TRUE)) {}

  # future
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
  Sys.setenv("R_FUTURE_FORK_ENABLE" = "FALSE")
}
.force_serial_execution()

## Create run_info directory (before writing any run_info/* files)
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)


## ===== Parameters =====
set.seed(1234)

csn_subunits <- c("GPS1", "COPS2", "COPS3", "COPS4", "COPS5", "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9")


log_msg <- function(text, ..., .envir = parent.frame()) {
  ts <- format(Sys.time(), "%H:%M:%S")
  msg <- tryCatch(
    {
      if (grepl("\\{[^}]+\\}", text)) {
        # has { } -> use glue
        glue::glue(text, ..., .envir = .envir)
      } else if (grepl("%", text)) {
        # has % -> use sprintf
        do.call(sprintf, c(list(fmt = text), list(...)))
      } else {
        # plain text; if extra args, treat as sprintf format
        if (nargs() > 1L) do.call(sprintf, c(list(fmt = text), list(...))) else text
      }
    },
    error = function(e) text
  )
  cat(sprintf("[%s] %s\n", ts, msg))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_fs_name <- function(s) {
  s <- gsub('[<>:"/\\\\|?*]', "_", s)
  s <- gsub("\\s+", "_", s)
  s
}



## ===== File reading and common utilities =====
read_case_list <- function(path_file) {
  if (!file.exists(path_file)) {
    return(character(0))
  }
  x <- readLines(path_file, warn = FALSE, encoding = "UTF-8")
  line <- x[grepl("^case_list_ids:", x)]
  if (!length(line)) {
    return(character(0))
  }
  ids <- sub("^case_list_ids:\\s*", "", line[1])
  ids <- unlist(strsplit(ids, "[,\\s]+"))
  unique(ids[nchar(ids) > 0])
}

load_matrix_from_dataset_dir <- function(dir) {
  fp <- file.path(dir, "data_protein_quantification.txt")
  if (!file.exists(fp)) stop(glue::glue("File not found: {fp}"))
  log_msg("Reading protein matrix: {basename(fp)}")
  dat <- suppressMessages(readr::read_tsv(fp, guess_max = 200000, show_col_types = FALSE))
  gene_cols <- c("Hugo_Symbol", "hugo_symbol", "Gene", "Gene_Symbol", "HugoSymbol", "GENE_SYMBOL", "gene", "gene_symbol")
  gcol <- intersect(gene_cols, names(dat))
  if (!length(gcol)) gcol <- names(dat)[1]
  dat <- dplyr::rename(dat, Gene = !!gcol[1])
  dat$Gene <- sub("\\|.*$", "", dat$Gene)
  not_sample <- c("Gene", "Entrez_Gene_Id", "Entrez_Gene_Id.", "ENTREZ_GENE_ID", "Description", "Gene_Name", "GeneName", "Gene_Symbol")
  sample_cols_all <- setdiff(names(dat), not_sample)
  # If folder has case_list, try to match
  case_file <- file.path(dir, "case_lists", "cases_protein_quantification.txt")
  keep_ids <- read_case_list(case_file)
  if (length(keep_ids)) {
    inter <- intersect(sample_cols_all, keep_ids)
    sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
    if (length(inter) < 10) log_msg("Note: case_list intersection too small ({length(inter)}), using all sample columns")
  } else {
    sample_cols <- sample_cols_all
  }
  m <- dat %>%
    dplyr::select(Gene, dplyr::all_of(sample_cols)) %>%
    janitor::remove_empty("cols")
  rn <- m$Gene
  m <- as.matrix(m[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- rn
  if (ncol(m) == 0) stop("Read 0 sample columns")
  if (anyDuplicated(rownames(m))) {
    log_msg("Detected duplicate genes, averaging duplicate rows")
    m <- rowsum(m, group = rownames(m), reorder = FALSE) / as.vector(table(rownames(m)))
  }
  log_msg("Matrix dimensions: {nrow(m)} genes x {ncol(m)} samples")
  m
}


# For PCA pre-cleaning: remove Inf->NA, keep rows with at least min_samples finite values; row-wise median imputation for NA
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

# Safe version of CSN SCORE (try original method first; fallback if failed)
build_csn_score_safe <- function(mat0, subunits, combine_7AB = TRUE,
                                 min_members = 5L, pca_min_samples = 10L) {
  sub <- intersect(subunits, rownames(mat0))
  out_na <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
  if (length(sub) < min_members) {
    return(out_na)
  }

  cs_try <- try(
    {
      build_csn_score(mat0, subunits = sub, combine_7AB = combine_7AB, min_members = min_members)
    },
    silent = TRUE
  )

  if (!inherits(cs_try, "try-error") && sum(is.finite(cs_try)) >= pca_min_samples) {
    return(cs_try)
  }

  X <- .clean_for_pca(mat0[sub, , drop = FALSE], min_samples = pca_min_samples, min_genes = min_members)
  if (is.null(X)) {
    if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] Insufficient subunits or samples -> all NA")
    return(out_na)
  }
  pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) {
    if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] prcomp failed -> all NA")
    return(out_na)
  }
  sc <- pc$x[, 1]
  names(sc) <- rownames(pc$x)
  ref <- colMeans(X, na.rm = TRUE)
  rr <- suppressWarnings(stats::cor(sc, ref, use = "pairwise.complete.obs"))
  if (is.finite(rr) && rr < 0) sc <- -sc
  out <- out_na
  out[names(sc)] <- as.numeric(sc)
  if (exists("log_msg", mode = "function")) {
    varpc1 <- if (!is.null(pc$sdev)) (pc$sdev[1]^2) / sum(pc$sdev^2) else NA_real_
    log_msg(
      "[CSN_SCORE-safe] fallback: genes=%d; PC1%%=%.1f; nonNA=%d/%d",
      nrow(X), 100 * varpc1, sum(is.finite(out)), length(out)
    )
  }
  out
}


## ==== PATCH: helpers ====

## === NEW: CSN complex score (PC1) & residualization helpers ==================

# Use z-scores of subunits across samples for PCA, take PC1 as CSN score; correct direction to match mean z sign
build_csn_score <- function(mat0,
                            subunits = csn_subunits,
                            combine_7AB = TRUE,
                            min_members = 5L) {
  present <- intersect(subunits, rownames(mat0))
  # Pre-build return skeleton (must have names)
  s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
  if (!length(present)) {
    return(s)
  }

  # z-score (preserve sample names)
  get_z <- function(v) {
    nm <- names(v) # save sample names first
    v <- as.numeric(v)
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    out <- (v - mu) / sdv
    names(out) <- nm # restore sample names
    out
  }

  X <- do.call(rbind, lapply(present, function(g) get_z(mat0[g, ])))
  rownames(X) <- present
  colnames(X) <- colnames(mat0) # This line is important: ensure sample names exist

  # Optional: merge COPS7A/7B
  if (combine_7AB && all(c("COPS7A", "COPS7B") %in% rownames(X))) {
    Z7 <- colMeans(X[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
    X <- rbind(X[setdiff(rownames(X), c("COPS7A", "COPS7B")), , drop = FALSE],
      "COPS7*" = Z7
    )
  }

  # Use non-NA count from original mat0 to determine coverage
  enough <- colSums(is.finite(mat0[present, , drop = FALSE])) >= min_members
  keep_sam <- names(s)[enough]

  if (length(keep_sam) >= 10) {
    pc <- stats::prcomp(t(X[, keep_sam, drop = FALSE]), center = TRUE, scale. = FALSE)
    sc <- pc$x[, 1]
    # Direction correction: same sign as subunit mean z
    mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)
    if (suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs")) < 0) sc <- -sc
    s[keep_sam] <- sc # Sample name alignment must succeed here
  }

  s
}

# Regress vector y (a subunit) on CSN score + batch + other covariates, then take residuals
# Replaces original residualize_vector()
residualize_vector <- function(y, csn_score, batch = NULL, covars = NULL, min_n = 8L) {
  # 1) Align sample names
  if (is.null(names(y))) stop("[residualize_vector] y must be a named numeric vector")
  sam <- names(y)
  common <- sam
  if (!is.null(csn_score)) common <- intersect(common, names(csn_score))
  if (!is.null(batch)) common <- intersect(common, names(batch))
  if (!is.null(covars)) {
    rn <- rownames(as.data.frame(covars, check.names = FALSE))
    if (is.null(rn)) stop("[residualize_vector] covars must have rownames=samples")
    common <- intersect(common, rn)
  }
  if (length(common) < min_n) {
    return(setNames(rep(NA_real_, length(y)), names(y)))
  }

  # 2) Build data frame -> single model.matrix expansion
  DF <- data.frame(
    csn = as.numeric(csn_score[common]),
    row.names = common, check.names = FALSE
  )
  if (!is.null(batch)) {
    bf <- droplevels(batch[common])
    DF[["batch"]] <- bf
  }
  if (!is.null(covars)) {
    C <- as.data.frame(covars[common, , drop = FALSE], check.names = FALSE)
    C <- coerce_covariates_safely(C) # <- added
    for (cn in colnames(C)) DF[[cn]] <- C[[cn]] # <- put directly back to DF
  }
  DF_y <- suppressWarnings(as.numeric(y[common]))

  ## === NEW: First remove columns that would make complete.cases all FALSE ===
  # 2a) Drop all-NA columns
  all_na_col <- vapply(DF, function(v) all(is.na(v)), logical(1))
  if (any(all_na_col)) {
    if (exists("log_msg", mode = "function")) {
      try(log_msg(
        "  [resid] drop all-NA columns: %s",
        paste(names(all_na_col)[all_na_col], collapse = ",")
      ), silent = TRUE)
    }
    DF <- DF[, !all_na_col, drop = FALSE]
  }

  # 2b) Drop single-level factors or constant numeric columns (to avoid singular design matrix)
  is_const <- vapply(DF, function(v) {
    vv <- v[!is.na(v)]
    if (!length(vv)) {
      TRUE
    } else {
      if (is.factor(v)) nlevels(droplevels(vv)) <= 1 else stats::var(as.numeric(vv)) == 0
    }
  }, logical(1))
  if (any(is_const)) {
    if (exists("log_msg", mode = "function")) {
      try(log_msg(
        "  [resid] drop constant/single-level columns: %s",
        paste(names(is_const)[is_const], collapse = ",")
      ), silent = TRUE)
    }
    DF <- DF[, !is_const, drop = FALSE]
  }
  ## === /NEW ===

  ok <- is.finite(DF_y) & stats::complete.cases(DF)
  if (sum(ok) < min_n) {
    out <- setNames(rep(NA_real_, length(y)), names(y))
    return(out)
  }
  DF <- DF[ok, , drop = FALSE]
  y_ok <- DF_y[ok]

  # 3) Expand design matrix
  # - batch if exists and is factor, auto dummy encode
  des <- stats::model.matrix(~ 1 + ., data = DF, na.action = stats::na.fail)
  stopifnot(nrow(des) == length(y_ok))

  # 4) Linear regression y ~ csn + batch + covars
  fit <- lm.fit(x = des, y = y_ok)
  res <- rep(NA_real_, length(DF_y))
  res[ok] <- fit$residuals

  # 5) Map back to full samples and (optionally) z-score
  out <- setNames(rep(NA_real_, length(y)), names(y))
  out[common] <- res
  fin <- is.finite(out)
  if (sum(fin) >= 3L) {
    mu <- mean(out[fin])
    sdv <- stats::sd(out[fin])
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    out[fin] <- (out[fin] - mu) / sdv
  }
  out
}


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


## ===== Batch cleaning settings (you can adjust as needed) =====
BATCH_PIPE_POLICY <- "NA" # "NA" or "b_small": set values containing '|' to NA or merge to b_small
BATCH_MIN_PER_LEVEL <- 2 # Levels below this threshold will be merged to b_small (to avoid non-estimable coefficients)

## ---- Batch value cleaning: handle '|', sanitize names, merge sparse levels ----
sanitize_batch_levels <- function(x,
                                  pipe_policy = BATCH_PIPE_POLICY,
                                  min_per_level = BATCH_MIN_PER_LEVEL) {
  x0 <- as.character(x)
  has_pipe <- grepl("\\|", x0 %||% "")
  if (any(has_pipe)) {
    log_msg("  [batch] Detected {sum(has_pipe)} values containing '|' -> processing by policy {pipe_policy}")
    x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
  }
  # Sanitize names (to avoid eBayes/design matrix column name issues)
  fac <- factor(make.names(x0))
  fac <- droplevels(fac)

  # Merge sparse levels (e.g. those with only 1 sample)
  if (!is.null(min_per_level) && min_per_level > 1) {
    tab <- table(fac, useNA = "no")
    small <- names(tab)[tab < min_per_level]
    if (length(small)) {
      log_msg(
        "  [batch] Merging sparse levels to 'b_small': %s",
        paste(sprintf("%s(n=%d)", small, as.integer(tab[small])), collapse = ", ")
      )
      fac_chr <- as.character(fac)
      fac_chr[fac_chr %in% small] <- "b_small"
      fac <- droplevels(factor(fac_chr))
    }
  }
  fac
}

## ---- Auto-detect batch column (with cleaning pipeline) ----
detect_batch_column <- function(meta,
                                pipe_policy = BATCH_PIPE_POLICY,
                                min_per_level = BATCH_MIN_PER_LEVEL) {
  cand <- c("TMT_PLEX", "EXPERIMENT")
  hit <- intersect(cand, colnames(meta))
  if (!length(hit)) {
    return(NULL)
  }

  for (cn in hit) {
    fac <- sanitize_batch_levels(meta[[cn]],
      pipe_policy   = pipe_policy,
      min_per_level = min_per_level
    )
    # Must have at least 2 valid levels and valid sample count >= 3
    if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
      return(list(name = cn, fac = fac))
    }
  }
  NULL
}

## ---- Get batch factor by sample order (already cleaned) ----
get_batch_factor <- function(ds_dir, sample_ids,
                             pipe_policy = BATCH_PIPE_POLICY,
                             min_per_level = BATCH_MIN_PER_LEVEL) {
  meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  if (!file.exists(meta_fp)) {
    return(NULL)
  }

  meta <- suppressMessages(
    readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")
  ) |> as.data.frame()

  id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"), names(meta))
  if (!length(id_cols)) {
    return(NULL)
  }

  meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
  meta$SAMPLE_ID <- as.character(meta$SAMPLE_ID)

  # First align by sample_ids, ensure detection and return vector length consistent
  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
  rownames(meta) <- sample_ids

  det <- detect_batch_column(meta,
    pipe_policy   = pipe_policy,
    min_per_level = min_per_level
  )
  if (is.null(det)) {
    ## === [NEW] Fallback: derive TMT-plex as batch from TMT_protein.csv ===
    tmt_fp <- file.path(ds_dir, "TMT_protein.csv")
    if (file.exists(tmt_fp)) {
      tmt <- suppressMessages(readr::read_csv(tmt_fp, show_col_types = FALSE)) |> as.data.frame()

      ## Normalize column names (remove leading/trailing spaces; handle " Run Metadata ID")
      cn <- names(tmt)
      cn_trim <- trimws(cn)
      names(tmt) <- cn_trim
      run_hits <- grep("^Run\\s*Metadata\\s*ID$", cn_trim, ignore.case = TRUE, value = TRUE)
      tmt_cols <- grep("^tmt_", cn_trim, ignore.case = TRUE, value = TRUE)

      if (length(run_hits) >= 1 && length(tmt_cols) >= 1) {
        run_col <- run_hits[1]

        ## Build sample_id -> plex mapping (one plex per row; all samples in tmt_* columns belong to this plex)
        plex_by_sample <- list()
        nR <- nrow(tmt)
        for (i in seq_len(nR)) {
          run_id <- as.character(tmt[[run_col]][i])
          if (!nzchar(run_id) || is.na(run_id)) next
          plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id) # take first two digits
          if (!nzchar(plex2) || is.na(plex2)) next

          for (tc in tmt_cols) {
            cell <- tmt[[tc]][i]
            if (is.na(cell)) next
            cell <- as.character(cell)
            if (!nzchar(cell)) next
            sid <- sub("\\r?\\n.*$", "", cell) # take content before newline
            sid <- trimws(sid)
            if (!nzchar(sid)) next
            if (is.null(plex_by_sample[[sid]])) plex_by_sample[[sid]] <- plex2
          }
        }

        if (length(plex_by_sample)) {
          v <- rep(NA_character_, length(sample_ids))
          names(v) <- sample_ids
          mm <- intersect(names(plex_by_sample), sample_ids)
          if (length(mm)) v[mm] <- unlist(plex_by_sample[mm], use.names = FALSE)

          fac2 <- sanitize_batch_levels(v,
            pipe_policy   = pipe_policy,
            min_per_level = min_per_level
          )
          ## Success condition: at least 2 valid levels and valid sample count >= 3
          if (nlevels(fac2) >= 2 && sum(!is.na(fac2)) >= 3) {
            names(fac2) <- sample_ids
            return(list(name = "TMT_protein.csv:RunMetadataID", fac = fac2))
          }
        }
      }
    }
    ## fallback also failed -> return NULL (maintain original semantics)
    return(NULL)
  }

  fac <- det$fac # already aligned with sample_ids
  names(fac) <- sample_ids
  list(name = det$name, fac = fac)
}


datasets_root <- getwd()

dataset_ids <- c(
  "brca_cptac_2020"
)


dataset_dirs <- setNames(file.path(datasets_root, dataset_ids), dataset_ids)
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")
message("datasets_root = ", datasets_root)

## Old version was stopifnot(all(dir.exists(dataset_dirs))), would abort entire batch when some directories are missing.
## Changed to: list missing ones, filter to "actually runnable" set using dataset_dirs_run.
missing_dirs <- names(dataset_dirs)[!dir.exists(dataset_dirs)]
if (length(missing_dirs)) {
  log_msg("Detected %d missing directories, will skip: %s", length(missing_dirs), paste(missing_dirs, collapse = ", "))
}


## List of datasets actually runnable this round (directory exists and contains data_protein_quantification.txt)
dataset_dirs_run <- dataset_dirs[
  dir.exists(dataset_dirs) &
    file.exists(file.path(dataset_dirs, "data_protein_quantification.txt"))
]

if (!length(dataset_dirs_run)) stop("dataset_dirs_run is empty, please check if directories and data_protein_quantification.txt exist")

log_msg("Available datasets this round (%d): %s", length(dataset_dirs_run), paste(names(dataset_dirs_run), collapse = ", "))


## ===== Get TP53 mutation status =====
## ===== TP53 status (protein-altering baseline) =====
## Only count protein-altering variants as TP53-mutant; others (Silent, UTR, Intron, IGR, RNA, lincRNA, Flank...) are treated as wild type

TP53_KEEP_CLASSES <- c(
  "MISSENSE_MUTATION", "NONSENSE_MUTATION",
  "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
  "IN_FRAME_DEL", "IN_FRAME_INS",
  "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
)

normalize_vc <- function(x) {
  # Standardize Variant_Classification: uppercase, unify various delimiters to underscore
  x <- toupper(trimws(as.character(x)))
  gsub("[^A-Z0-9]+", "_", x)
}

get_tp53_status <- function(ds_dir, sample_ids) {
  # Default all to wild-type
  status <- setNames(rep("TP53_wild_type", length(sample_ids)), sample_ids)

  mut_fp <- file.path(ds_dir, "data_mutations.txt")
  if (!file.exists(mut_fp)) {
    log_msg("  [TP53] data_mutations.txt not found, treating all as wild-type/ALL available")
    return(status)
  }

  mutation_df <- tryCatch(
    readr::read_tsv(mut_fp, comment = "#", show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(mutation_df)) {
    return(status)
  }

  req <- c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode")
  miss <- setdiff(req, colnames(mutation_df))
  if (length(miss)) {
    log_msg("  [TP53] data_mutations.txt missing columns: %s -> treating as wild-type", paste(miss, collapse = ", "))
    return(status)
  }

  tp53_df <- subset(mutation_df, Hugo_Symbol == "TP53",
    select = c("Variant_Classification", "Tumor_Sample_Barcode")
  )
  if (!nrow(tp53_df)) {
    return(status)
  }

  vc_norm <- normalize_vc(tp53_df$Variant_Classification)

  # Strictly use protein-altering classes; also tolerate FRAME_SHIFT_* / IN_FRAME_* prefixes
  keep <- vc_norm %in% TP53_KEEP_CLASSES |
    grepl("^FRAME_SHIFT_", vc_norm) |
    grepl("^IN_FRAME_", vc_norm)

  if (!any(keep)) {
    return(status)
  }

  tp53_samples <- toupper(unique(tp53_df$Tumor_Sample_Barcode[keep]))
  sid_up <- toupper(sample_ids)
  status[sid_up %in% tp53_samples] <- "TP53_mutant"
  status
}

## ===== TP53 status auditing: per-dataset sample counts =====
## Output:
##   - run_info/tp53_status/<dataset>_tp53_class_sample_counts.csv
##      (original Variant_Classification -> "sample count"; also includes Any_TP53_mutation and protein_altering)
##   - run_info/tp53_status/tp53_class_sample_counts_long.csv (aggregated long format)
##   - run_info/tp53_status/tp53_binary_counts_by_dataset.csv (wild type / mutant summary table)

dir.create(file.path("run_info", "tp53_status"), recursive = TRUE, showWarnings = FALSE)

# Safety: if not defined earlier, add here
if (!exists("normalize_vc", mode = "function")) {
  normalize_vc <- function(x) {
    x <- toupper(trimws(as.character(x)))
    gsub("[^A-Z0-9]+", "_", x)
  }
}
if (!exists("TP53_KEEP_CLASSES")) {
  TP53_KEEP_CLASSES <- c(
    "MISSENSE_MUTATION", "NONSENSE_MUTATION",
    "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
    "IN_FRAME_DEL", "IN_FRAME_INS",
    "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
  )
}


## ===== Missingness sensitivity（AGE）=====
USE_AGE_MISSING_INDICATOR <- FALSE # Main analysis default: no missing-indicator, keep NA only

## Safe z-score: estimate mean/sd using only finite values, keep NA, no mean imputation
.z_no_impute <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  fin <- is.finite(v)
  mu <- mean(v[fin], na.rm = TRUE)
  sdv <- stats::sd(v[fin], na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) sdv <- 1
  out <- (v - mu) / sdv
  out[!fin] <- NA_real_
  out
}

## Utilities: column name normalization, z-score, 0~1
.norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))
.zscore <- function(v) {
  v <- as.numeric(v)
  mu <- mean(v[is.finite(v)], na.rm = TRUE)
  sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) sdv <- 1
  v[!is.finite(v)] <- mu
  (v - mu) / sdv
}
.to01 <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  if (sum(is.finite(v) & v > 1, na.rm = TRUE) > sum(is.finite(v) & v <= 1, na.rm = TRUE)) v <- v / 100
  pmin(pmax(v, 0), 1)
}

## For PAAD: take median of semicolon-separated values
.median_from_semicolon <- function(x_chr) {
  vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
  vv <- vv[is.finite(vv)]
  if (!length(vv)) {
    return(NA_real_)
  }
  stats::median(vv)
}


## 3) Get purity by dataset (named numeric, names=sample_ids, 0~1)
get_purity_covariate <- function(ds_id, ds_dir, sample_ids) {
  samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  samp <- if (file.exists(samp_fp)) suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) else NULL
  if (!is.null(samp)) {
    samp <- as.data.frame(samp)
    names(samp) <- .norm_names(names(samp))
  }
  purity <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)

  if (ds_id == "brca_cptac_2020") {
    pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
    if (file.exists(pat_fp) && !is.null(samp)) {
      pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
      names(pat) <- .norm_names(names(pat))
      sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
      pid_in_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
      pid_in_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
      if (!is.na(sid) && !is.na(pid_in_samp) && !is.na(pid_in_pat) && "ESTIMATE_TUMORPURITY" %in% names(pat)) {
        map_pt <- setNames(as.character(samp[[pid_in_samp]]), samp[[sid]])
        pt_pur <- setNames(.to01(pat[["ESTIMATE_TUMORPURITY"]]), as.character(pat[[pid_in_pat]]))
        purity[] <- unname(pt_pur[map_pt[sample_ids]])
      }
    }
  } else if (ds_id == "luad_cptac_2020") {
    if (!is.null(samp)) {
      sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
      if (!is.na(sid) && "TUMOR_PURITY_BYESTIMATE_RNASEQ" %in% names(samp)) {
        purity[] <- .to01(samp[["TUMOR_PURITY_BYESTIMATE_RNASEQ"]][match(sample_ids, samp[[sid]])])
      }
    }
  } else if (ds_id == "lusc_cptac_2021") {
    if (!is.null(samp)) {
      sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
      if (!is.na(sid) && "ESTIMATE_TUMORPURITY" %in% names(samp)) {
        purity[] <- .to01(samp[["ESTIMATE_TUMORPURITY"]][match(sample_ids, samp[[sid]])])
      }
    }
  } else if (ds_id == "paad_cptac_2021") {
    if (!is.null(samp)) {
      sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
      if (!is.na(sid) && "NEOPLASTIC_CELLULARITY" %in% names(samp)) {
        medv <- vapply(as.character(samp[["NEOPLASTIC_CELLULARITY"]]), .median_from_semicolon, numeric(1))
        purity[] <- .to01(medv[match(sample_ids, samp[[sid]])])
      }
    }
  } else if (ds_id == "ucec_cptac_2020") {
    if (!is.null(samp)) {
      sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
      if (!is.na(sid)) {
        pc <- suppressWarnings(as.numeric(samp[["PURITY_CANCER"]]))
        pi <- suppressWarnings(as.numeric(samp[["PURITY_IMMUNE"]]))
        ps <- suppressWarnings(as.numeric(samp[["PURITY_STROMA"]]))
        idx <- match(sample_ids, samp[[sid]])
        purity_calc <- ifelse(is.finite(pc[idx]), pc[idx], 1 - (pi[idx] %||% 0) - (ps[idx] %||% 0))
        purity[] <- pmin(pmax(purity_calc, 0), 1)
      }
    }
  }
  purity
}


## ========================
## Get sex/age (aligned to sample; sex: 0/1; age: z-score)
## ========================

if (!exists(".norm_names", inherits = FALSE)) {
  .norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))
}
if (!exists(".zscore", inherits = FALSE)) {
  .zscore <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    (v - mu) / sdv
  }
}

## Main process: strictly use patient file's SEX and AGE to generate covariates (sex: 0/1; age: z-score)
get_sex_age_covariates <- function(ds_dir, sample_ids) {
  samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")

  if (!file.exists(samp_fp) || !file.exists(pat_fp)) {
    log_msg("[covars] Missing sample/patient file, using NA for sex/age")
    out <- cbind(
      sex = rep(NA_real_, length(sample_ids)),
      age = rep(NA_real_, length(sample_ids))
    )
    rownames(out) <- sample_ids
    return(out)
  }

  .NN <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))
  .z_no_impute <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    fin <- is.finite(v)
    mu <- mean(v[fin], na.rm = TRUE)
    sdv <- stats::sd(v[fin], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    out <- (v - mu) / sdv
    out[!fin] <- NA_real_
    out
  }
  # If .zscore (mean-imputed then z) not defined externally, provide backup version to maintain consistent logic
  if (!exists(".zscore", mode = "function")) {
    .zscore <- function(v) {
      v <- suppressWarnings(as.numeric(v))
      mu <- mean(v[is.finite(v)], na.rm = TRUE)
      sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      v[!is.finite(v)] <- mu
      (v - mu) / sdv
    }
  }

  samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
  pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
  names(samp) <- .NN(names(samp))
  names(pat) <- .NN(names(pat))

  sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
  pid_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
  pid_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
  if (is.na(sid) || is.na(pid_samp) || is.na(pid_pat)) {
    stop("[get_sex_age_covariates] Cannot build sample<->patient mapping (missing SAMPLE_ID/PATIENT_ID)")
  }
  map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])

  ## SEX: Male=1, Female=0
  sex_col <- intersect(c("SEX", "GENDER"), names(pat))[1]
  sex <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  if (!is.na(sex_col)) {
    raw <- toupper(as.character(pat[[sex_col]]))
    val <- ifelse(grepl("^M", raw), 1, ifelse(grepl("^F", raw), 0, NA_real_))
    names(val) <- as.character(pat[[pid_pat]])
    sex[] <- unname(val[map_pt[sample_ids]])
  }

  ## AGE: main analysis no imputation; sensitivity can optionally use missing-indicator + mean-imputed z
  age_col <- intersect(c("AGE", "AGE_AT_DIAGNOSIS", "AGE_AT_INDEX", "AGE_YEARS"), names(pat))[1]
  age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  age_missing <- age_z_imputed <- NULL
  if (!is.na(age_col)) {
    v <- suppressWarnings(as.numeric(pat[[age_col]]))
    names(v) <- as.character(pat[[pid_pat]])
    age_raw <- unname(v[map_pt[sample_ids]])
    age[] <- .z_no_impute(age_raw) # NA preserved (main analysis)
    if (!exists("USE_AGE_MISSING_INDICATOR", inherits = FALSE)) USE_AGE_MISSING_INDICATOR <- FALSE
    if (isTRUE(USE_AGE_MISSING_INDICATOR)) {
      age_missing <- as.numeric(is.na(age_raw))
      age_z_imputed <- .zscore(age_raw) # mean-imputed z (sensitivity)
    }
  }

  cov_sex <- mean(is.finite(sex)) * 100
  cov_age <- mean(is.finite(age)) * 100
  if (isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))) {
    cov_age_imp <- mean(is.finite(age_z_imputed)) * 100
    cov_age_mis <- mean(is.finite(age_missing)) * 100
    log_msg(
      "    covariates coverage：sex %.1f%%, age(NA-as-NA) %.1f%%, age_z_imputed %.1f%%, age_missing %.1f%%",
      cov_sex, cov_age, cov_age_imp, cov_age_mis
    )
  } else {
    log_msg("    covariates coverage：sex %.1f%%, age(NA-as-NA) %.1f%%", cov_sex, cov_age)
  }

  out <- data.frame(
    sex = as.numeric(sex), age = as.numeric(age),
    row.names = sample_ids, check.names = FALSE
  )
  out$age_z <- out$age
  if (isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))) {
    out$age_missing <- as.numeric(age_missing)
    out$age_z_imputed <- as.numeric(age_z_imputed)
  }
  out
}


## ==============================
## Mini audit: sex/age + purity + batch overall check
## ==============================

# Safety: utility functions (if not yet defined)
if (!exists(".norm_names", inherits = FALSE)) {
  .norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))
}
if (!exists("%||%", inherits = FALSE)) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
if (!exists(".to01", inherits = FALSE)) {
  .to01 <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    # If points >1 outnumber points <=1, treat as 0~100 percentage, convert to 0~1
    if (sum(is.finite(v) & v > 1, na.rm = TRUE) > sum(is.finite(v) & v <= 1, na.rm = TRUE)) v <- v / 100
    pmin(pmax(v, 0), 1)
  }
}
.median_from_semicolon <- function(x_chr) {
  vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
  vv <- vv[is.finite(vv)]
  if (!length(vv)) {
    return(NA_real_)
  }
  stats::median(vv)
}

# Utility: format batch level sample sizes into string
.format_batch_sizes <- function(fac) {
  if (is.null(fac)) {
    return(NA_character_)
  }
  tb <- sort(table(fac), decreasing = TRUE)
  paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
}


## =========================================================
## [NEW] BRCA_immune_analysis: immune/stromal score x predictors correlation for brca_cptac_2020 only
## - Source: data_clinical_patient.txt (5 columns)
## - Predictors: CSN subunits protein abundance + CSN_SCORE + RESIDUAL_<subunit>
## - Methods: Pearson/Spearman; RAW and BatchAdj covariate strategy same as protein GSEA version
## - Stratification: ALL / TP53_mutant / TP53_wild_type; also do TP53 interaction for ALL
## - Output: <ds_dir>/BRCA_immune_analysis/<stratum>/<version>/*.csv
## =========================================================

.run_brca_immune_cor <- function(dataset_dirs_map = NULL, make_plots = FALSE, min_pairs = 8L) {
  ds_id <- "brca_cptac_2020"
  # Parse ds_dir
  ds_dir <- NULL
  if (!is.null(dataset_dirs_map) && ds_id %in% names(dataset_dirs_map)) ds_dir <- dataset_dirs_map[[ds_id]]
  if (is.null(ds_dir) && exists("dataset_dirs_run")) {
    dsm <- get("dataset_dirs_run", inherits = TRUE)
    if (ds_id %in% names(dsm)) ds_dir <- dsm[[ds_id]]
  }
  if (is.null(ds_dir) && exists("dataset_dirs")) {
    dsm <- get("dataset_dirs", inherits = TRUE)
    if (ds_id %in% names(dsm)) ds_dir <- dsm[[ds_id]]
  }
  if (is.null(ds_dir) || !dir.exists(ds_dir)) {
    message("[BRCA-immune] Cannot find brca_cptac_2020 directory, skipping")
    return(invisible(NULL))
  }

  out_root <- "BRCA_immune_analysis"
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  # [NEW] Plot filters and output formats (can be overridden by global variables or options)
  plot_filters <- get0("BRCA_PLOT_FILTERS", ifnotfound = list(
    stratum   = NULL, # e.g. c("ALL","TP53_mutant","TP53_wild_type")
    version   = NULL, # e.g. c("RAW","BatchAdj")
    score     = NULL, # e.g. c("ESTIMATE_IMMUNE_SCORE","XCELL_IMMUNE_SCORE", ...)
    predictor = NULL, # e.g. c("CSN_SCORE","GPS1","RESIDUAL_COPS5", ...)
    method    = NULL # e.g. c("pearson","spearman")
  ))
  # Default output only .tiff; can change to c("TIFF","PNG","JPG","PDF")
  plot_formats <- toupper(get0("BRCA_PLOT_FORMATS", ifnotfound = c("TIFF")))
  plot_dpi <- get0("BRCA_PLOT_DPI", ifnotfound = 600L) # commonly used 600 dpi for submission
  plot_w <- get0("BRCA_PLOT_WIDTH", ifnotfound = 4.0) # inch
  plot_h <- get0("BRCA_PLOT_HEIGHT", ifnotfound = 3.0) # inch

  # ===== [NEW] Cell-style theme and color palette =====
  CELL_PALETTES <- list(
    # Original four sets
    CELL_BLUE        = list(pt = "#2D2D2D", ln = "#1F77B4"), # dark gray point / blue line
    CELL_TEAL        = list(pt = "#2D2D2D", ln = "#1B9E77"), # dark gray point / teal line
    CELL_ORANGE      = list(pt = "#2D2D2D", ln = "#D55E00"), # dark gray point / orange line
    CELL_PURPLE      = list(pt = "#2D2D2D", ln = "#6A3D9A"), # dark gray point / purple line
    # New: brighter series (for regression line color ln; point default still dark gray pt, can override separately)
    CELL_BRIGHT_BLUE = list(pt = "#2D2D2D", ln = "#0077FF"),
    CELL_AQUA        = list(pt = "#2D2D2D", ln = "#00D5C6"),
    CELL_MINT        = list(pt = "#2D2D2D", ln = "#00C853"),
    CELL_MAGENTA     = list(pt = "#2D2D2D", ln = "#E91E63"),
    CELL_SCARLET     = list(pt = "#2D2D2D", ln = "#FF3B30"),
    CELL_GOLD        = list(pt = "#2D2D2D", ln = "#F6B400")
  )

  # Can use BRCA_PLOT_COLORSET global variable to override (e.g. "CELL_TEAL")
  plot_colorset <- toupper(get0("BRCA_PLOT_COLORSET", ifnotfound = "CELL_BLUE"))
  if (!plot_colorset %in% names(CELL_PALETTES)) plot_colorset <- "CELL_BLUE"

  # Lighter CI fill color
  .ci_fill <- function(col) tryCatch(grDevices::adjustcolor(col, alpha.f = 0.25), error = function(e) col)

  # Cell-like theme: no outer frame, thin black axis lines, small text; all title text smaller and positioned outside chart
  .cell_theme <- function(base_size = 10) {
    # Can override font size via global variables; if not set, use default ratios
    title_sz <- get0("BRCA_TITLE_SIZE", ifnotfound = base_size * 1.0)
    subtitle_sz <- get0("BRCA_SUBTITLE_SIZE", ifnotfound = base_size * 0.9)
    axis_title_sz <- get0("BRCA_AXIS_TITLE_SIZE", ifnotfound = base_size * 0.9)
    axis_text_sz <- get0("BRCA_AXIS_TEXT_SIZE", ifnotfound = base_size * 0.85)

    ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(), # no outer frame
        axis.line = ggplot2::element_line(linewidth = 0.6, colour = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.6, colour = "black"),
        plot.title = ggplot2::element_text(
          face = "bold", size = title_sz,
          margin = ggplot2::margin(b = 2)
        ),
        plot.subtitle = ggplot2::element_text(
          size = subtitle_sz,
          margin = ggplot2::margin(b = 4)
        ),
        axis.title = ggplot2::element_text(size = axis_title_sz),
        axis.text = ggplot2::element_text(size = axis_text_sz)
      )
  }

  # Read protein matrix & TP53
  mat0_full <- load_matrix_from_dataset_dir(ds_dir)
  tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full))

  # Only keep samples where CSN subunits exist, and prepare score data
  present_sub <- intersect(csn_subunits, rownames(mat0_full))
  if (!length(present_sub)) {
    message("[BRCA-immune] This dataset lacks CSN subunits, skipping")
    return(invisible(NULL))
  }
  scores_all <- .read_brca_scores_from_patient(ds_dir, colnames(mat0_full))

  strata_list <- list(
    ALL = colnames(mat0_full),
    TP53_mutant = names(tp53_status)[tp53_status == "TP53_mutant"],
    TP53_wild_type = names(tp53_status)[tp53_status == "TP53_wild_type"]
  )
  versions <- c("RAW", "BatchAdj")

  # Calculate by stratum x version
  for (st_name in names(strata_list)) {
    keep <- intersect(colnames(mat0_full), strata_list[[st_name]])
    if (length(keep) < 4L) {
      message(sprintf("[BRCA-immune][%s] Too few samples (%d), skipping", st_name, length(keep)))
      next
    }
    mat0 <- mat0_full[, keep, drop = FALSE]
    # CSN_SCORE (calculated using stratified samples)
    csn <- build_csn_score_safe(mat0, present_sub, combine_7AB = TRUE)
    csn <- setNames(as.numeric(csn[colnames(mat0)]), colnames(mat0))

    # outcome (five scores)
    Y <- scores_all[keep, , drop = FALSE]

    # Covariate preparation (aligned to samples)
    purity <- get_purity_covariate(ds_id, ds_dir, keep)
    sa <- get_sex_age_covariates(ds_dir, keep) # sex, age
    cov0 <- cbind(sex = sa[, "sex"], age = sa[, "age"], purity = purity)
    rownames(cov0) <- keep
    bi <- get_batch_factor(ds_dir, keep)
    batch0 <- if (!is.null(bi)) droplevels(bi$fac[keep]) else NULL
    if (!is.null(batch0)) names(batch0) <- keep

    for (ver in versions) {
      if (ver == "RAW") {
        batch_use <- NULL
        cov_use <- cov0
      } else { # BatchAdj
        batch_use <- batch0
        cov_use <- cov0
      }

      # Prepare predictors: each subunit, CSN_SCORE, RESIDUAL_<subunit>
      predictors <- list()
      for (su in present_sub) {
        v <- as.numeric(mat0[su, ])
        names(v) <- keep
        predictors[[su]] <- v
      }
      predictors[["CSN_SCORE"]] <- csn
      for (su in present_sub) {
        vv <- as.numeric(mat0[su, ])
        names(vv) <- keep
        predictors[[paste0("RESIDUAL_", su)]] <- residualize_vector(
          y = vv, csn_score = csn, batch = batch_use, covars = cov_use
        )
      }

      rows <- list()
      # (default no plotting) Prepare plot directory
      plot_dir <- file.path(out_root, st_name, ver, "plots")
      if (isTRUE(make_plots)) dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
      # [NEW] Collect plot jobs; wait until BH adjustment is done before actually plotting (to correctly label padj)
      plot_jobs <- list()

      for (sc in colnames(Y)) {
        y0 <- suppressWarnings(as.numeric(Y[[sc]]))
        names(y0) <- keep
        for (pn in names(predictors)) {
          x0 <- predictors[[pn]]
          ok <- is.finite(x0) & is.finite(y0)
          if (sum(ok) < min_pairs) next

          
          # Prepare data for LM
          # X: x0[ok], Y: y0[ok]
          # Covariates: cov_use[ok, ], Batch: batch_use[ok]
          df_model <- data.frame(
            y = y0[ok],
            x = x0[ok],
            check.names = FALSE
          )
          # Add covariates
          if (!is.null(cov_use)) {
            cov_sub <- as.data.frame(cov_use[ok, , drop = FALSE])
            df_model <- cbind(df_model, cov_sub)
          }
          # Add batch
          if (!is.null(batch_use)) {
            df_model$batch <- droplevels(batch_use[ok])
          }

          # Helper to extract partial correlation r and p from lm fit
          .extract_partial_stats <- function(fit, term = "x") {
            # If fit failed or singular
            if (inherits(fit, "try-error") || is.na(coef(fit)[term])) {
              return(list(r = NA_real_, p = NA_real_))
            }
            sm <- summary(fit)
            # t-statistic and degrees of freedom
            t_val <- sm$coefficients[term, "t value"]
            df_resid <- fit$df.residual
            p_val <- sm$coefficients[term, "Pr(>|t|)"]

            # Partial correlation from t-statistic: r = t / sqrt(t^2 + df)
            # Sign of r is same as sign of t
            r_part <- t_val / sqrt(t_val^2 + df_resid)

            list(r = r_part, p = p_val)
          }

          # 1) PEARSON (Linear Model on raw values)
          # Model: y ~ x + cov1 + cov2 + ...
          form_p <- as.formula(paste("y ~ .")) # y against all other columns in df_model
          fit_p <- try(stats::lm(form_p, data = df_model), silent = TRUE)
          stats_p <- .extract_partial_stats(fit_p, "x")

          # For equation string (optional, just for display), we can still show the simple relationship
          # from the 'residualized' space (Partial Regression Plot line), or the main effect.
          # Here we keep the original visible equation style (y_resid ~ x_resid),
          # because the plot will show residuals.
          # So we still need .residualize_pair for the PLOT and EQUATION display.
          xy <- .residualize_pair(x0[ok], y0[ok], batch = batch_use[ok], covars = as.data.frame(cov_use[ok, , drop = FALSE]))

          # Refit simple lm on residuals JUST for the plot line equation / R2 display
          lm_resid_p <- try(stats::lm(xy$y ~ xy$x), silent = TRUE)
          intercept_p <- if (!inherits(lm_resid_p, "try-error")) unname(coef(lm_resid_p)[1]) else NA_real_
          slope_p <- if (!inherits(lm_resid_p, "try-error")) unname(coef(lm_resid_p)[2]) else NA_real_
          # Note: The R2 of this residual plot is the squared partial correlation
          r2_p <- if (is.finite(stats_p$r)) stats_p$r^2 else NA_real_
          eq_p <- sprintf("y = %.6f + %.6f * x", intercept_p, slope_p)


          # 2) SPEARMAN (Linear Model on Ranks)
          # We rank x and y, then fit the SAME covariate model
          df_rank <- df_model
          df_rank$y <- rank(df_model$y)
          df_rank$x <- rank(df_model$x)
          

          fit_s <- try(stats::lm(form_p, data = df_rank), silent = TRUE)
          stats_s <- .extract_partial_stats(fit_s, "x")

          
          xy_rank <- .residualize_pair(rank(x0[ok]), rank(y0[ok]), batch = batch_use[ok], covars = as.data.frame(cov_use[ok, , drop = FALSE]))

          lm_resid_s <- try(stats::lm(xy_rank$y ~ xy_rank$x), silent = TRUE)
          intercept_s <- if (!inherits(lm_resid_s, "try-error")) unname(coef(lm_resid_s)[1]) else NA_real_
          slope_s <- if (!inherits(lm_resid_s, "try-error")) unname(coef(lm_resid_s)[2]) else NA_real_
          r2_s <- if (is.finite(stats_s$r)) stats_s$r^2 else NA_real_
          eq_s <- sprintf("rank(y) = %.6f + %.6f * rank(x)", intercept_s, slope_s)

          rows[[length(rows) + 1L]] <- data.frame(
            dataset = ds_id, stratum = st_name, version = ver,
            score = sc, predictor = pn,
            method = "pearson", n = sum(ok),
            r = as.numeric(stats_p$r), p = as.numeric(stats_p$p),
            padj = NA_real_, R2 = r2_p, equation = eq_p,
            stringsAsFactors = FALSE, check.names = FALSE
          )
          rows[[length(rows) + 1L]] <- data.frame(
            dataset = ds_id, stratum = st_name, version = ver,
            score = sc, predictor = pn,
            method = "spearman", n = sum(ok),
            r = as.numeric(stats_s$r), p = as.numeric(stats_s$p),
            padj = NA_real_, R2 = r2_s, equation = eq_s,
            stringsAsFactors = FALSE, check.names = FALSE
          )

          # (Optional) Scatter plot: now collect plot jobs (since padj is not yet available)
          if (isTRUE(make_plots)) {
            # Decide whether to collect based on user-specified filter conditions
            want_stratum <- is.null(plot_filters$stratum) || st_name %in% plot_filters$stratum
            want_version <- is.null(plot_filters$version) || ver %in% plot_filters$version
            want_score <- is.null(plot_filters$score) || sc %in% plot_filters$score
            want_predictor <- is.null(plot_filters$predictor) || pn %in% plot_filters$predictor
            if (want_stratum && want_version && want_score && want_predictor) {
              # Visual: always use residuals (Partial Regression plot style)
              dfp <- data.frame(x = xy$x, y = xy$y)

              methods_want <- c("pearson", "spearman")
              if (!is.null(plot_filters$method)) {
                methods_want <- intersect(methods_want, tolower(plot_filters$method))
              }
              if (length(methods_want)) {
                for (mm in methods_want) {
                  
                  rr <- if (mm == "pearson") as.numeric(stats_p$r) else as.numeric(stats_s$r)
                  pp <- if (mm == "pearson") as.numeric(stats_p$p) else as.numeric(stats_s$p)

                  
                  r2_val <- if (mm == "pearson") r2_p else r2_s
                  eq_val <- if (mm == "pearson") eq_p else eq_s

                  plot_jobs[[length(plot_jobs) + 1L]] <- list(
                    dataset = ds_id, stratum = st_name, version = ver,
                    score = sc, predictor = pn, method = mm,
                    r = rr, p = pp,
                    R2 = r2_val,
                    equation = eq_val,
                    df = dfp, n = sum(ok)
                  )
                }
              }
            }
          }
        }
      } # end scores

      out_dir <- file.path(out_root, st_name, ver)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      if (length(rows)) {
        out_tbl <- do.call(rbind, rows)
        # Group by (dataset x stratum x version x method x score) for BH adjustment
        out_tbl <- split(
          out_tbl,
          list(out_tbl$dataset, out_tbl$stratum, out_tbl$version, out_tbl$method, out_tbl$score),
          drop = TRUE
        )
        out_tbl <- lapply(out_tbl, function(d) {
          d$padj <- stats::p.adjust(d$p, method = "BH")
          d
        })
        # [NEW] Plot based on BH-adjusted results from out_tbl (Cell style; with CI; subtitle shows r and padj; no overlap with chart)
        if (isTRUE(make_plots) && length(plot_jobs)) {
          ot <- do.call(rbind, out_tbl)
          pal <- CELL_PALETTES[[plot_colorset]]
          ln_col <- pal$ln
          pt_col <- get0("BRCA_POINT_COLOR", ifnotfound = pal$pt)
          pt_size <- as.numeric(get0("BRCA_POINT_SIZE", ifnotfound = 1.6))
          ci_col <- .ci_fill(ln_col)

          for (job in plot_jobs) {
            # Re-check method filter (in case external changed)
            if (!is.null(plot_filters$method) &&
              !(tolower(job$method) %in% tolower(plot_filters$method))) {
              next
            }

            # Get corresponding padj
            idx <- which(ot$dataset == job$dataset & ot$stratum == job$stratum &
              ot$version == job$version & ot$method == job$method &
              ot$score == job$score & ot$predictor == job$predictor)
            padj_txt <- "NA"
            if (length(idx)) {
              padj_val <- ot$padj[idx[1]]
              if (is.finite(padj_val)) padj_txt <- format(padj_val, digits = 3, scientific = TRUE)
            }

            title_txt <- sprintf("%s vs %s", job$score, job$predictor)
            sub_txt <- sprintf(
              "%s | %s / %s / %s | r=%.3f, padj=%s, n=%d",
              toupper(job$method), job$dataset, job$stratum, job$version,
              job$r, padj_txt, job$n
            )

            p1 <- ggplot2::ggplot(job$df, ggplot2::aes(x = x, y = y)) +
              ggplot2::geom_point(size = pt_size, alpha = 0.85, colour = pt_col) +
              ggplot2::geom_smooth(
                method = "lm", se = TRUE, linewidth = 0.8,
                colour = ln_col, fill = ci_col
              ) +
              ggplot2::labs(
                title = title_txt,
                subtitle = sub_txt, # r / padj in subtitle -> no overlap with chart
                x = sprintf("%s (residualized)", job$predictor),
                y = sprintf("%s (residualized)", job$score)
              ) +
              .cell_theme(base_size = 10)

            # Output multiple formats (default only .tiff; can override with BRCA_PLOT_FORMATS)
            for (fmt in plot_formats) {
              ext <- switch(toupper(fmt),
                "TIFF" = "tiff",
                "TIF" = "tiff",
                "PNG" = "png",
                "JPG" = "jpg",
                "JPEG" = "jpg",
                "PDF" = "pdf",
                tolower(fmt)
              )
              fn <- sprintf(
                "%s_vs_%s_%s_%s_%s.%s",
                job$score, job$predictor, job$stratum, job$version, job$method, ext
              )
              ggplot2::ggsave(
                filename = file.path(plot_dir, fn),
                plot = p1,
                width = plot_w, height = plot_h, dpi = plot_dpi,
                bg = "white", device = ext, limitsize = FALSE
              )
            }
          }
        }
        out_tbl <- do.call(rbind, out_tbl)
        data.table::fwrite(out_tbl, file.path(out_dir, "immune_vs_predictor_correlations.csv"))
        message(sprintf(
          "[BRCA-immune][%s/%s] Output %d results -> %s",
          st_name, ver, nrow(out_tbl),
          file.path(out_dir, "immune_vs_predictor_correlations.csv")
        ))
      } else {
        # Empty file also output for downstream review
        header <- data.frame(
          dataset = character(0), stratum = character(0), version = character(0),
          score = character(0), predictor = character(0), method = character(0),
          correlation = numeric(0), p = numeric(0), padj = numeric(0),
          R2 = numeric(0), equation = character(0),
          stringsAsFactors = FALSE, check.names = FALSE
        )
        data.table::fwrite(header, file.path(out_dir, "immune_vs_predictor_correlations.csv"))
        message(sprintf(
          "[BRCA-immune][%s/%s] No qualified pairs (rows=0); empty file output.",
          st_name, ver
        ))
      }
    }

    ## ---- TP53 interaction (ALL only; follows GSEA workflow spirit) ----
    if (identical(st_name, "ALL")) {
      tp53_num <- as.numeric(tp53_status[keep] == "TP53_mutant")
      names(tp53_num) <- keep
      for (ver in versions) {
        if (ver == "RAW") {
          batch_use <- NULL
          cov_use <- cov0
        } else {
          batch_use <- batch0
          cov_use <- cov0
        }
        rows_i <- list()
        for (sc in colnames(Y)) {
          y0 <- suppressWarnings(as.numeric(Y[[sc]]))
          names(y0) <- keep
          for (pn in names(predictors)) {
            x0 <- predictors[[pn]]
            names(x0) <- keep
            ok <- is.finite(x0) & is.finite(y0) & is.finite(tp53_num)
            if (sum(ok) < min_pairs) next
            dfm <- data.frame(
              y = y0[ok], x = x0[ok], tp53 = tp53_num[ok],
              cov_use[ok, , drop = FALSE]
            )
            if (!is.null(batch_use)) dfm$batch <- droplevels(batch_use[ok])
            # y ~ x * tp53 + covariates (+ batch)
            form <- as.formula(paste0(
              "y ~ x * tp53 + ",
              paste(colnames(cov_use), collapse = " + "),
              if (!is.null(batch_use)) " + batch" else ""
            ))
            fit <- try(stats::lm(form, data = dfm), silent = TRUE)
            if (inherits(fit, "try-error")) next
            sm <- summary(fit)
            co <- coef(sm)
            beta_int <- if ("x:tp53" %in% rownames(co)) co["x:tp53", "Estimate"] else NA_real_
            p_int <- if ("x:tp53" %in% rownames(co)) co["x:tp53", "Pr(>|t|)"] else NA_real_
            r2m <- as.numeric(sm$r.squared)
            eq_line <- paste0(
              "y = ", sprintf("%.6f", unname(coef(fit)[1])),
              " + ", sprintf("%.6f", unname(coef(fit)["x"])), " * x",
              " + ", sprintf("%.6f", unname(coef(fit)["tp53"])), " * tp53",
              " + ", sprintf("%.6f", beta_int), " * x:tp53 + ..."
            )
            rows_i[[length(rows_i) + 1L]] <- data.frame(
              dataset = ds_id, stratum = "ALL", version = ver,
              score = sc, predictor = pn,
              method = "lm_interaction",
              beta_interaction = beta_int, p = p_int, padj = NA_real_,
              R2 = r2m, equation = eq_line,
              stringsAsFactors = FALSE, check.names = FALSE
            )
          }
        }
        out_dir <- file.path(out_root, "ALL", ver)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        if (length(rows_i)) {
          out_i <- do.call(rbind, rows_i)
          out_i$padj <- ave(out_i$p, out_i$dataset, out_i$stratum, out_i$version, out_i$score,
            FUN = function(v) stats::p.adjust(v, method = "BH")
          )
          data.table::fwrite(out_i, file.path(out_dir, "immune_vs_predictor_interaction.csv"))
          message(sprintf(
            "[BRCA-immune][ALL/%s] Interaction: output %d rows -> %s",
            ver, nrow(out_i),
            file.path(out_dir, "immune_vs_predictor_interaction.csv")
          ))
        } else {
          header_i <- data.frame(
            dataset = character(0), stratum = character(0), version = character(0),
            score = character(0), predictor = character(0), method = character(0),
            beta_interaction = numeric(0), p = numeric(0), padj = numeric(0),
            R2 = numeric(0), equation = character(0),
            stringsAsFactors = FALSE, check.names = FALSE
          )
          data.table::fwrite(header_i, file.path(out_dir, "immune_vs_predictor_interaction.csv"))
          message(sprintf("[BRCA-immune][ALL/%s] Interaction: no qualified results, empty file output.", ver))
        }
      }
    } # end interaction
  } # end strata
  message("[BRCA-immune] Complete: output -> ", out_root)
}

.read_brca_scores_from_patient <- function(ds_dir, sample_ids) {
  samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
  if (!file.exists(pat_fp)) {
    message("[BRCA-immune] Cannot find data_clinical_patient.txt: ", ds_dir)
    out <- as.data.frame(matrix(NA_real_, nrow = length(sample_ids), ncol = 0),
      stringsAsFactors = FALSE
    )
    rownames(out) <- sample_ids
    return(out)
  }
  .NN <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))
  samp <- if (file.exists(samp_fp)) suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame() else NULL
  pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
  if (!is.null(samp)) names(samp) <- .NN(names(samp))
  names(pat) <- .NN(names(pat))

  sid <- if (!is.null(samp)) intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1] else NA
  pid_samp <- if (!is.null(samp)) intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1] else NA
  pid_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
  map_pt <- NULL
  if (!is.na(sid) && !is.na(pid_samp)) map_pt <- setNames(as.character(samp[[pid_samp]]), as.character(samp[[sid]]))

  keep_cols <- c(
    "ESTIMATE_IMMUNE_SCORE", "ESTIMATE_STROMAL_SCORE",
    "XCELL_IMMUNE_SCORE", "XCELL_STROMAL_SCORE", "CIBERSORT_ABSOLUTE_SCORE"
  )
  out <- matrix(NA_real_,
    nrow = length(sample_ids), ncol = length(keep_cols),
    dimnames = list(sample_ids, keep_cols)
  )
  if (!is.null(pat) && !is.na(pid_pat)) {
    for (nm in intersect(keep_cols, names(pat))) {
      vec_pt <- suppressWarnings(as.numeric(pat[[nm]]))
      names(vec_pt) <- as.character(pat[[pid_pat]])
      if (!is.null(map_pt)) {
        out[, nm] <- vec_pt[map_pt[sample_ids]]
      } else {
        # Fallback: align using patient ID with suffix removed
        pids <- sub("([A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+).*", "\\1", toupper(sample_ids))
        out[, nm] <- vec_pt[pids]
      }
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

.residualize_pair <- function(x, y, batch = NULL, covars = NULL, min_n = 8L) {
  # If residualize_to_covars is defined globally, use it; otherwise use orthogonalize_to as fallback
  if (exists("residualize_to_covars", mode = "function")) {
    xr <- residualize_to_covars(x, batch = batch, covars = covars)
    yr <- residualize_to_covars(y, batch = batch, covars = covars)
  } else {
    X <- NULL
    if (!is.null(batch)) X <- cbind(X, stats::model.matrix(~ 0 + batch))
    if (!is.null(covars) && ncol(covars) > 0) X <- cbind(X, as.matrix(covars))
    xr <- orthogonalize_to(x, X)
    yr <- orthogonalize_to(y, X)
  }
  list(x = xr, y = yr)
}


BRCA_TITLE_SIZE <- 8 # main title
BRCA_SUBTITLE_SIZE <- 5 # subtitle (for r and padj)
BRCA_AXIS_TITLE_SIZE <- 8 # x / y axis title
BRCA_AXIS_TEXT_SIZE <- 8 # tick numbers


BRCA_PLOT_FILTERS <- list(
  stratum   = c("ALL"), # choose one or more
  version   = c("RAW", "BatchAdj"),
  score     = c("CIBERSORT_ABSOLUTE_SCORE"),
  predictor = c("CSN_SCORE", "COPS7A", "COPS7B"), # example
  method    = c("pearson")
)

# score     = c("ESTIMATE_IMMUNE_SCORE","ESTIMATE_STROMAL_SCORE","XCELL_IMMUNE_SCORE",
#              "XCELL_STROMAL_SCORE","CIBERSORT_ABSOLUTE_SCORE")

# BRCA_PLOT_COLORSET <- "CELL_BLUE"
# BRCA_PLOT_COLORSET <- "CELL_TEAL"
# BRCA_PLOT_COLORSET <- "CELL_ORANGE"
# BRCA_PLOT_COLORSET <- "CELL_PURPLE"

# Choose color scheme (line color): brighter
# BRCA_PLOT_COLORSET <- "CELL_BRIGHT_BLUE"   # or CELL_BRIGHT_BLUE / CELL_AQUA / CELL_SCARLET / CELL_GOLD / CELL_MINT
# BRCA_PLOT_COLORSET <- "CELL_AQUA"
# BRCA_PLOT_COLORSET <- "CELL_SCARLET"
# BRCA_PLOT_COLORSET <- "CELL_MAGENTA"
# BRCA_PLOT_COLORSET <- "CELL_GOLD"
BRCA_PLOT_COLORSET <- "CELL_MINT"

"#0077FF" # CELL_BRIGHT_BLUE
"#00D5C6" # CELL_AQUA
"#00C853" # CELL_MINT
"#E91E63" # CELL_MAGENTA
"#FF3B30" # CELL_SCARLET
"#F6B400" # CELL_GOLD


# Adjust scatter point size and color (optional)
BRCA_POINT_SIZE <- 0.9
BRCA_POINT_COLOR <- "#111111" # If not specified, defaults to palette pt (dark gray)


BRCA_PLOT_FORMATS <- c("TIFF") # default is tiff can be omitted; for multiple formats: c("TIFF","PNG","JPG","PDF")
# BRCA_PLOT_FORMATS <- c("TIFF","PNG","JPG","PDF")  # for all formats like this
BRCA_PLOT_DPI <- 600

tryCatch(
  {
    .run_brca_immune_cor(
      dataset_dirs_map = get0("dataset_dirs_run", ifnotfound = NULL),
      make_plots = TRUE
    )
  },
  error = function(e) {
    message("[BRCA-immune] Execution failed: ", conditionMessage(e))
    # Simple call stack (avoid R's traceback not working well in tryCatch)
    calls <- sys.calls()
    message(
      "[BRCA-immune] Call stack: ",
      paste(utils::capture.output(print(calls)), collapse = " | ")
    )
    stop(e) # Explicitly throw to let outer workflow see the error
  }
)


## ===== [NEW] Immune heatmap color options (follows GSEA default, can be switched) =====
BRCA_HM_COLORSET <- get0("BRCA_HM_COLORSET", ifnotfound = "GSEA_DEFAULT")
.BRCA_HM_PALETTES <- list(
  GSEA_DEFAULT      = c(neg = "#053061", mid = "#FFFFFF", pos = "#67001F"), # blue-white-burgundy (original default)
  BLUE_RED          = c(neg = "#2166AC", mid = "#F7F7F7", pos = "#B2182B"), # similar to RdBu
  BLUE_ORANGE       = c(neg = "#2B8CBE", mid = "#F7F7F7", pos = "#E34A33"),
  GREEN_MAGENTA     = c(neg = "#1B7837", mid = "#F7F7F7", pos = "#762A83"),

  # --- Added red-blue series (more vivid/different shades) ---
  BLUE_RED_DEEP     = c(neg = "#08306B", mid = "#F7F7F7", pos = "#7F0000"), # deep ocean blue <-> deep burgundy
  BLUE_RED_BRIGHT   = c(neg = "#1F78B4", mid = "#FFFFFF", pos = "#E31A1C"), # bright blue <-> bright red
  BLUE_RED_LIGHT    = c(neg = "#6BAED6", mid = "#F7F7F7", pos = "#FB6A4A"), # light blue <-> light red
  BLUE_RED_TWILIGHT = c(neg = "#2C7FB8", mid = "#EEEEEE", pos = "#D7301F"), # soft blue <-> warm red
  BLUE_RED_CRIMSON  = c(neg = "#0B3C5D", mid = "#FAFAFA", pos = "#B80C09") # navy blue <-> crimson
)


# Then call as in your original workflow (if plotting, set make_plots=TRUE and specify desired conditions in BRCA_PLOT_FILTERS)


## ===== [NEW] Immediately execute BRCA_immune_analysis (only brca_cptac_2020; default no plotting) =====

## ===== [NEW] Immune score x predictors correlation heatmap (follows GSEA single dataset style) =====
# Description:
# - Read BRCA_immune_analysis/<dataset>/<stratum>/<version>/immune_vs_predictor_correlations.csv
# - x-axis: predictor ordering/spacing same as  GSEA single dataset heatmap (uses .make_heatmap_plot_with_yorder_colored)
#   approach: temporarily add "NES_" prefix to predictor names to apply existing .pred_order_all and .x_build_positions rules
# - y-axis (top to bottom): CIBERSORT_ABSOLUTE_SCORE, XCELL_IMMUNE_SCORE, ESTIMATE_IMMUNE_SCORE,
#                            XCELL_STROMAL_SCORE, ESTIMATE_STROMAL_SCORE
# - Color: by r (correlation coefficient); padj < 0.05 marked with dot
# - Options: stratum / version / method / colorset / output format (default only .tiff)

.plot_brca_immune_heatmap <- function(
  dataset_id = "brca_cptac_2020",
  dataset_dirs_map = get0("dataset_dirs_run", ifnotfound = NULL),
  stratum = c("ALL", "TP53_mutant", "TP53_wild_type"),
  version = c("RAW", "BatchAdj"),
  method = c("pearson", "spearman"),
  hm_colorset = get0("BRCA_HM_COLORSET", ifnotfound = "GSEA_DEFAULT"),
  out_formats = get0("BRCA_PLOT_FORMATS", ifnotfound = c("TIFF")),
  width = 7.0, height = 3.2, dpi = get0("BRCA_PLOT_DPI", ifnotfound = 600)
) {
  stratum <- match.arg(stratum)
  version <- match.arg(version)
  method <- match.arg(method)

  # 1) Find immune_vs_predictor_correlations.csv
  ds_root <- if (!is.null(dataset_dirs_map) && dataset_id %in% names(dataset_dirs_map)) {
    dataset_dirs_map[[dataset_id]]
  } else {
    NULL
  }
  if (is.null(ds_root) || !dir.exists(ds_root)) {
    stop("Cannot find dataset directory: ", dataset_id, " (dataset_dirs_map not set or path does not exist)")
  }

  # Try several candidate paths based on actual output location
  candidates <- c(
    file.path("BRCA_immune_analysis", stratum, version, "immune_vs_predictor_correlations.csv"),
    file.path(ds_root, "BRCA_immune_analysis", stratum, version, "immune_vs_predictor_correlations.csv")
  )
  csv_file <- NA_character_
  for (p in candidates) {
    if (file.exists(p)) {
      csv_file <- p
      break
    }
  }
  if (!is.character(csv_file) || is.na(csv_file) || !nzchar(csv_file)) {
    stop("Cannot find immune correlation CSV: ", paste(candidates, collapse = " | "))
  }

  dt <- data.table::fread(csv_file, na.strings = c("NA", "NaN", ""))
  if (!nrow(dt)) stop("File is empty: ", csv_file)

  # 2) Filter by method (PEARSON / SPEARMAN)
  keep <- grepl(method, dt$method, ignore.case = TRUE)
  dt <- dt[keep, ]
  if (!nrow(dt)) stop("No rows matching method: ", method, " (check CSV's method column)")

  # 3) padj column selection: prefer method-specific padj, otherwise fall back to generic padj
  padj_col <- if (tolower(method) == "pearson" && "padj_pearson" %in% names(dt)) {
    "padj_pearson"
  } else if (tolower(method) == "spearman" && "padj_spearman" %in% names(dt)) {
    "padj_spearman"
  } else if ("padj" %in% names(dt)) "padj" else stop("Cannot find padj column in CSV")

  # 4) Prepare long format for .make_heatmap_plot_with_yorder_colored
  #    To use GSEA heatmap's same x-axis ordering/spacing, temporarily add 'NES_' prefix to predictor
  df_long <- data.frame(
    pathway = dt$score,
    predictor = paste0("NES_", dt$predictor),
    NES = dt$r,
    padj = dt[[padj_col]],
    stringsAsFactors = FALSE
  )

  # 5) y-axis order (top -> bottom)
  y_order <- c(
    "CIBERSORT_ABSOLUTE_SCORE",
    "XCELL_IMMUNE_SCORE",
    "ESTIMATE_IMMUNE_SCORE",
    "XCELL_STROMAL_SCORE",
    "ESTIMATE_STROMAL_SCORE"
  )

  # 6) Get palette
  pal <- .BRCA_HM_PALETTES[[hm_colorset]]
  if (is.null(pal)) pal <- .BRCA_HM_PALETTES[["GSEA_DEFAULT"]]
  # If read an unnamed three-color vector, add neg/mid/pos names to match internal function requirements
  if (is.character(pal) && is.null(names(pal)) && length(pal) >= 3) {
    pal <- c(neg = pal[1], mid = pal[2], pos = pal[3])
  }

  # 7) 繪圖（完全沿用你 GSEA 單一 dataset 熱圖的內部函式）
  title0 <- paste0("Immune-score vs Predictors (", toupper(method), ")")
  subt0 <- paste(dataset_id, stratum, version, sep = " / ")
  # 7) Plot (use existing helper: the function takes csv_file, not df_long and extra parameters)
  title0 <- paste0("Immune-score vs Predictors (", toupper(method), ")")
  subt0 <- paste(dataset_id, stratum, version, sep = " / ")

  # Write temporary long table CSV for .make_heatmap_plot_with_yorder_colored to read
  tmp_long_csv <- tempfile(pattern = "immune_long_", fileext = ".csv")
  data.table::fwrite(df_long, tmp_long_csv)

  # Get heatmap: helper may return ggplot, or list(plot=..., width=..., height=...)
  g_res <- .make_heatmap_plot_with_yorder_colored(
    csv_file = tmp_long_csv,
    y_order  = y_order,
    palette  = pal
  )

  # Decompose return object
  g <- if (is.list(g_res) && !is.null(g_res$plot)) g_res$plot else g_res
  if (is.list(g_res) && !is.null(g_res$width) && !is.null(g_res$height)) {
    width <- g_res$width
    height <- g_res$height
  }

  # Add title/subtitle
  g <- g + ggplot2::labs(title = title0, subtitle = subt0)

  # ---- sanity check：g 必須是 ggplot 物件，否則視為生成失敗 ----
  if (is.null(g) || !inherits(g, "ggplot")) {
    stop("[BRCA-immune][heatmap] Failed to generate heatmap object: g is NULL or not ggplot.")
  }


  # 8) Output (default only .tiff; other formats can be specified via out_formats)
  out_dir <- dirname(csv_file)
  base_bn <- paste0("immune_vs_predictor_heatmap_", toupper(method))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fmts <- toupper(out_formats)
  fmts <- if (length(fmts)) fmts else c("TIFF")
  any_ok <- FALSE
  for (f in fmts) {
    target <- file.path(out_dir, paste0(base_bn, ".", tolower(f)))
    if (f == "TIFF") {
      tryCatch(
        {
          ragg::agg_tiff(
            filename = target, width = width, height = height,
            units = "in", res = dpi, compression = "lzw"
          )
          # Don't show "NULL" in console
          invisible(print(g))
          grDevices::dev.off()
          message(sprintf("[BRCA-immune][heatmap] Output: %s", target))
        },
        error = function(e) {
          message(sprintf(
            "[BRCA-immune][heatmap] Failed to save TIFF: %s -> %s",
            conditionMessage(e), target
          ))
          stop(e)
        }
      )
      # Check file after tryCatch to avoid scoping issues with <<-
      if (file.exists(target) && file.size(target) > 0) any_ok <- TRUE
    } else if (f %in% c("PNG", "JPG", "JPEG", "PDF")) {
      tryCatch(
        {
          ggplot2::ggsave(
            filename = target, plot = g,
            width = width, height = height, units = "in", dpi = dpi,
            bg = "white", limitsize = FALSE
          )
          message(sprintf("[BRCA-immune][heatmap] Output: %s", target))
        },
        error = function(e) {
          message(sprintf(
            "[BRCA-immune][heatmap] Failed to save %s: %s -> %s",
            f, conditionMessage(e), target
          ))
          stop(e)
        }
      )
      # Check file after tryCatch to avoid scoping issues with <<-
      if (file.exists(target) && file.size(target) > 0) any_ok <- TRUE
    } else {
      stop("[BRCA-immune][heatmap] Unsupported format: ", f)
    }
  }
  if (!any_ok) {
    stop("[BRCA-immune][heatmap] No output generated, please check if out_dir is writable: ", out_dir)
  }

  invisible(list(plot = g, file_base = file.path(out_dir, base_bn)))
}

.make_heatmap_plot_with_yorder_colored <- function(csv_file, y_order, palette = .get_interaction_palette()) {
  df_raw <- suppressMessages(.read_csv_safe(csv_file))
  path_candidates <- c("pathway", "Pathway", "term", "Term", "gs_name", "NAME", "set", "Set")
  path_col <- intersect(path_candidates, names(df_raw))[1]
  if (is.na(path_col)) stop("Cannot find pathway column: ", paste(path_candidates, collapse = ", "))

  has_long <- all(c("predictor", "NES", "padj") %in% names(df_raw))
  if (has_long) {
    df_long <- df_raw %>%
      dplyr::rename(pathway = dplyr::all_of(path_col)) %>%
      dplyr::mutate(
        predictor = as.character(.data$predictor),
        NES = as.numeric(.data$NES),
        padj = as.numeric(.data$padj)
      )
  } else {
    nes_cols <- grep("^NES_", names(df_raw), value = TRUE)
    padj_cols <- grep("^padj_", names(df_raw), value = TRUE)
    if (!length(nes_cols)) stop("Wide table format: cannot find NES_* columns.")
    df_wide <- df_raw %>% dplyr::rename(pathway = dplyr::all_of(path_col))
    nes_map <- tibble::tibble(nes_col = nes_cols, key = sub("^NES_", "", nes_cols))
    padj_map <- tibble::tibble(padj_col = padj_cols, key = sub("^padj_", "", padj_cols))
    pair_map <- dplyr::left_join(nes_map, padj_map, by = "key")
    df_list <- lapply(seq_len(nrow(pair_map)), function(i) {
      nes_c <- pair_map$nes_col[i]
      padj_c <- pair_map$padj_col[i]
      tibble::tibble(
        pathway   = df_wide$pathway,
        predictor = paste0("NES_", pair_map$key[i]),
        NES       = suppressWarnings(as.numeric(df_wide[[nes_c]])),
        padj      = if (!is.na(padj_c)) suppressWarnings(as.numeric(df_wide[[padj_c]])) else NA_real_
      )
    })
    df_long <- dplyr::bind_rows(df_list)
  }

  present_preds <- intersect(.pred_order_all, unique(df_long$predictor))
  df_long <- df_long %>% dplyr::filter(.data$predictor %in% present_preds)

  # y-axis uses ALL's order; only draw ALL's genesets and fully follow its order
  df_long <- df_long %>% dplyr::filter(.data$pathway %in% y_order)
  y_levels <- y_order[y_order %in% df_long$pathway]
  if (!length(y_levels)) stop("Ordered version: no intersection with ALL's pathway order.")
  df_long <- df_long %>% dplyr::mutate(pathway = factor(.data$pathway, levels = rev(y_levels)))


  # gap rules and positions
  gap <- 0.4
  needs_gap1 <- all(c("NES_CSN_SCORE", "NES_GPS1") %in% present_preds)
  needs_gap2 <- "NES_RESIDUAL_GPS1" %in% present_preds
  pos_map <- list()
  pos <- 0
  for (p in present_preds) {
    if (p == "NES_GPS1" && needs_gap1) pos <- pos + gap
    if (p == "NES_RESIDUAL_GPS1" && needs_gap2) pos <- pos + gap
    pos <- pos + 1
    pos_map[[p]] <- pos
  }
  pos_map <- unlist(pos_map)
  df_long <- df_long %>% dplyr::mutate(xpos = unname(pos_map[.data$predictor]))

  # Symmetric color scale range
  L <- max(abs(df_long$NES), na.rm = TRUE)
  if (!is.finite(L) || L == 0) L <- 1

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = pathway, fill = NES)) +
    ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
    ggplot2::geom_point(
      data = df_long %>% dplyr::filter(is.finite(.data$padj), .data$padj < 0.05),
      ggplot2::aes(x = xpos, y = pathway),
      shape = 16, size = 1.6, color = "black", inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = palette[["neg"]], mid = palette[["mid"]], high = palette[["pos"]],
      limits = c(-L, L), midpoint = 0, oob = scales::squish, name = "NES"
    ) +
    ggplot2::scale_x_continuous(
      breaks = unname(pos_map[present_preds]),
      labels = present_preds,
      expand = ggplot2::expansion(mult = c(0.01, 0.01)),
      position = "top"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0, size = 9, margin = ggplot2::margin(b = 8)),
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = "right",
      plot.margin = ggplot2::margin(6, 12, 6, 6),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  n_path <- nlevels(df_long$pathway)
  W <- max(8, length(present_preds) * 0.45)
  H <- max(6, n_path * 0.22)
  list(plot = p, width = W, height = H)
}

.read_csv_safe <- function(path) {
  p <- path

  # Windows: if original string doesn't exist, try replacing / with \
  if (!file.exists(p) && .Platform$OS.type == "windows") {
    p2 <- gsub("/", "\\\\", p, fixed = TRUE)
    if (file.exists(p2)) p <- p2
  }
  if (!file.exists(p)) stop("File does not exist: ", p)

  fi <- suppressWarnings(file.info(p))
  if (isTRUE(fi$isdir)) stop("Target is not a file: ", p)
  if (!is.na(fi$size) && fi$size == 0) stop("Empty file (size = 0): ", p)

  # 1) readr default
  out <- try(suppressMessages(readr::read_csv(p, show_col_types = FALSE, progress = FALSE)), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(out)
  }

  # 2) readr + UTF-8-BOM
  out <- try(suppressMessages(readr::read_csv(
    p,
    locale = readr::locale(encoding = "UTF-8-BOM"),
    show_col_types = FALSE, progress = FALSE
  )), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(out)
  }

  # 3) Don't manually open connection: read as full text string then parse
  txt <- try(readr::read_file(p), silent = TRUE)
  if (!inherits(txt, "try-error")) {
    out <- try(suppressMessages(readr::read_csv(
      readr::I(txt),
      show_col_types = FALSE, progress = FALSE
    )), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
  }

  # 4) Last resort: base R
  out <- try(utils::read.csv(p, check.names = FALSE), silent = TRUE)
  if (inherits(out, "try-error")) stop(out)
  tibble::as_tibble(out)
}

.pred_order_all <- c(
  "NES_CSN_SCORE", "NES_GPS1",
  "NES_COPS2", "NES_COPS3", "NES_COPS4", "NES_COPS5", "NES_COPS6",
  "NES_COPS7A", "NES_COPS7B", "NES_COPS8", "NES_COPS9",
  "NES_RESIDUAL_GPS1",
  "NES_RESIDUAL_COPS2", "NES_RESIDUAL_COPS3", "NES_RESIDUAL_COPS4", "NES_RESIDUAL_COPS5", "NES_RESIDUAL_COPS6",
  "NES_RESIDUAL_COPS7A", "NES_RESIDUAL_COPS7B", "NES_RESIDUAL_COPS8", "NES_RESIDUAL_COPS9"
)


.plot_brca_immune_heatmap(
  stratum = "ALL", version = "BatchAdj", method = "pearson",
  hm_colorset = "BLUE_RED_CRIMSON",
  out_formats = c("TIFF") # for multiple formats: c("TIFF","PNG","PDF","JPG")
)
