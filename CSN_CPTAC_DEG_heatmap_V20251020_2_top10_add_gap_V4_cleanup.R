## =========================================================
## CSN subunits -> proteomic GSEA (CPTAC/TCGA, add TP53 stratification)
##  - By strata = c("ALL","TP53_mutant","TP53_wild_type") output to
##    <dataset>/csn_gsea_results_TP53/<STRATUM>/...
##  - Also do pan-cancer integration and two summary tables to csn_gsea_pan_summary_TP53/
## =========================================================

setwd("C:/Users/danny/Documents/R_project/CSN_CPTAC") ## YTC laptop
## setwd("C:/Users/cmuh/Documents/YenTing_document/CSN_CPTAC") ## lab computer
getwd()


## ===== Add yaml to "Packages" block (1) =====
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
  library(fgsea)
  library(openxlsx)
  library(ComplexHeatmap)
  library(cowplot)
  library(matrixStats)
  library(yaml) # <— New: Used to write run_manifest.yml
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
        if (length(lv) <= 1) { # Single-level factor -> drop, avoid contrasts error
          keep[cn] <- FALSE
          if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
        } else {
          df[[cn]] <- v
        }
      } else {
        df[[cn]] <- suppressWarnings(as.numeric(v)) # Values remain numeric
        # If need to filter out "completely constant" numeric columns together, uncomment below:
        # vv <- df[[cn]]
        # if (sum(!is.na(vv)) >= 2 && isTRUE(all(vv[!is.na(vv)] == vv[which(!is.na(vv))[1]]))) {
        #   keep[cn] <- FALSE
        #   if (exists("logf")) try(logf("  [covars] drop constant numeric: %s", cn), silent = TRUE)
        # }
      }
    }
    df <- df[, keep, drop = FALSE]
    df
  }
}


## ==== Force single thread, close all parallel backends (across common frameworks) ====
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
  if (requireNamespace("doParallel", quietly = TRUE)) {
    # Try to close residual cluster (if you stored cluster in .GlobalEnv, stopCluster yourself)
    # Here we don't hard stop unknown objects, just use DoSEQ.
  }

  # future
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
  Sys.setenv("R_FUTURE_FORK_ENABLE" = "FALSE")

  # fgsea: Some versions support nproc, here only set an option for wrapper to read
  options(.fgsea_nproc = 1L)
}
.force_serial_execution()

## Create run_info directory (before writing any run_info/* files)
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)


## Attach multiple hypothesis level explanation (Update: add meta-FDR)
cat(
  paste(
    "Multiple-testing policy:",
    " - Within each gene-set group and per statistic, we control FDR using fgsea padj (per-dataset/stratum).",
    " - Pan-cancer summaries aggregate directions and counts across datasets/strata (descriptive).",
    " - Additionally, we provide a simple meta-analysis across datasets using Stouffer's z to obtain meta p-values,",
    "   followed by Benjamini-Hochberg to report meta-level FDR (see csn_gsea_pan_summary_TP53/meta_fdr).",
    sep = "\n"
  ),
  file = file.path("run_info", "analysis_notes.txt")
)



## ===== Parameters =====
set.seed(1234)

csn_subunits <- c("GPS1", "COPS2", "COPS3", "COPS4", "COPS5", "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9")

min_frac_complete <- 0.75 # At least 75% non-NA per gene
hi_lo_quantile <- 0.25 # High/Low each 25%


min_per_group <- 8 # High/Low at least how many samples each (skip limma if insufficient)
MAKE_PLOTS <- FALSE

## ---- Global switch: Whether to run High/Low sensitivity (Default FALSE; main analysis uses continuous predictor) ----
RUN_HILO_SENSITIVITY <- FALSE




## ---- DEG Analysis Mode ----
## This script performs DEG (Differential Expression Gene) analysis using limma continuous model
.RUN_DEG <- TRUE


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
        # Plain text; if extra arguments, treat as sprintf format
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

## ===== Gene set selection (Supports all MSigDB collections/subcollections) =====
## Optional tokens (case insensitive):
##   "H",
##   "C1",
##   "C2", "C2:REACTOME", "C2:BIOCARTA", "C2:CGP", "C2:KEGG"  # Note: Most new versions have no KEGG
##   "C3", "C3:TFT", "C3:MIR",
##   "C4", "C4:CGN", "C4:CM",
##   "C5", "C5:BP", "C5:CC", "C5:MF",
##   "C6", "C7", "C8"


# GENESET_GROUPS_TO_RUN <- c("H")  # ← Change as needed, e.g., c("H","C5:BP","C2:REACTOME","C7")
# GENESET_GROUPS_TO_RUN <- c("H", "C5:BP", "C2:REACTOME", "C3:TFT", "C7")
# GENESET_GROUPS_TO_RUN <- c("H","C1","C2","C3","C4","C5","C6","C7","C8")
# GENESET_GROUPS_TO_RUN <- c("C5:BP","C5:CC","C5:MF")
# GENESET_GROUPS_TO_RUN <- c("H","C3:MIR","C2:CP:REACTOME","C5:GO:BP")

GENESET_GROUPS_TO_RUN <- c("H")



## ===== Gene sets (Build only selected groups) =====
log_msg(
  "Preparing MSigDB gene sets; will build according to GENESET_GROUPS_TO_RUN: %s",
  paste(GENESET_GROUPS_TO_RUN, collapse = ", ")
)
# (Optional) View and output currently available collections / subcollections in msigdbr (including geneset counts)
# Placement: After log_msg(...) above, before genesets_by_group <- list()
df0 <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens"), error = function(e) NULL)

if (!is.null(df0) && NROW(df0) > 0) {
  df0 <- as.data.frame(df0) # 避免 tibble 子集嚴格行為

  ## 1) Detect actual column names (2025.1 is gs_collection / gs_subcollection)
  cat_candidates <- c("gs_collection", "gs_cat", "gs_category", "category", "collection")
  sub_candidates <- c("gs_subcollection", "gs_subcat", "gs_subcategory", "subcategory", "sub_category")
  name_candidates <- c("gs_name", "geneset_name", "set_name")

  cat_hits <- intersect(cat_candidates, names(df0))
  sub_hits <- intersect(sub_candidates, names(df0))
  name_hits <- intersect(name_candidates, names(df0))

  if (length(cat_hits) == 0 || length(name_hits) == 0) {
    log_msg("(Skipped) Necessary columns (collection or gs_name) not found; collections list not output.")
  } else {
    catcol <- cat_hits[1]
    subcol <- if (length(sub_hits) >= 1) sub_hits[1] else NULL
    namecol <- name_hits[1]

    ## 2) Assemble data frame and count "unique geneset (gs_name) count"
    tmp <- data.frame(
      gs_cat = as.character(df0[[catcol]]),
      gs_subcat = if (!is.null(subcol)) as.character(df0[[subcol]]) else NA_character_,
      gs_name = as.character(df0[[namecol]]),
      stringsAsFactors = FALSE
    )

    counts <- stats::aggregate(gs_name ~ gs_cat + gs_subcat,
      data = tmp,
      FUN = function(x) length(unique(x))
    )
    names(counts)[3] <- "n_genesets"

    ## 3) Sort (Treat NA subclass as empty string to avoid order encountering NA)
    ord <- order(counts$gs_cat, ifelse(is.na(counts$gs_subcat), "", counts$gs_subcat))
    counts <- counts[ord, , drop = FALSE]

    ## 4) Output list to run_info
    dir.create("run_info", showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(counts, file.path("run_info", "msigdb_available_collections.csv"))
    log_msg("Output msigdb available collections (including geneset counts) to run_info/msigdb_available_collections.csv")
  }
} else {
  log_msg("(Skipped) msigdbr() returned empty or call failed; collections list not output.")
}

genesets_by_group <- list()

.get_subcat_col <- function(df) {
  cand <- c("gs_subcollection", "gs_subcat", "gs_subcategory", "subcategory", "sub_category")
  hit <- cand[cand %in% names(df)]
  if (length(hit)) hit[1] else NULL
}


.fetch_msig <- function(cat) {
  try(msigdbr(species = "Homo sapiens", category = toupper(cat)), silent = TRUE)
}

.parse_group_token <- function(tok) {
  t <- toupper(trimws(tok))
  sp <- unlist(strsplit(t, ":", fixed = TRUE))
  list(cat = sp[1], sub = if (length(sp) >= 2) paste(sp[-1], collapse = ":") else NULL, raw = t)
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
    # Match only at the "beginning" of the subclass string (e.g., MIR, CP:REACTOME, GO:BP)
    # And allow following ":" or end
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


## ===== Read file and shared tools =====
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
    if (length(inter) < 10) log_msg("Hint: case_list intersection too small ({length(inter)}), using all sample columns instead")
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


## ===== [REMOVED] write_geneset_manifest - only needed for GSEA =====


## ===== [REMOVED] GSEA functions - only for GSEA pathway enrichment =====
## Removed: GSEA ranking and execution functions


# Clean for PCA: Remove Inf->NA, keep rows with at least min_samples finite values; row-wise median imputation for NA
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

# Safe version CSN SCORE (Try original method first; fallback if failed)
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
    if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] Available subunits or samples insufficient -> All NA")
    return(out_na)
  }
  pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) {
    if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] prcomp failed -> All NA")
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

# Safe version audit: Try original audit first; if failed, use clean matrix to run prcomp for record only, do not stop process
audit_csn_score_feasibility_safe <- function(ds_id, stratum, mat0, present_sub,
                                             min_members = 5L, pca_min_samples = 10L,
                                             min_per_group = NULL,
                                             out_dir = file.path("run_info", "csn_score_audit")) {
  ## [FIX] resolve min_per_group safely to avoid recursive default evaluation
  if (is.null(min_per_group)) {
    if (exists("min_per_group", envir = .GlobalEnv, inherits = FALSE)) {
      min_per_group <- get("min_per_group", envir = .GlobalEnv, inherits = FALSE)
    } else {
      min_per_group <- 8L
    }
  }
  ok <- try(
    {
      audit_csn_score_feasibility(
        ds_id = ds_id,
        stratum = stratum,
        mat0 = mat0,
        present_sub = present_sub,
        min_members = min_members,
        pca_min_samples = pca_min_samples,
        min_per_group = min_per_group,
        out_dir = out_dir
      )
      TRUE
    },
    silent = TRUE
  )
  if (!inherits(ok, "try-error")) {
    return(invisible(TRUE))
  }

  X <- .clean_for_pca(mat0[intersect(present_sub, rownames(mat0)), , drop = FALSE],
    min_samples = pca_min_samples, min_genes = min_members
  )
  if (is.null(X)) {
    if (exists("log_msg", mode = "function")) {
      log_msg("[CSN-audit-safe] %s | %s: Available subunits or samples insufficient, skip audit (do not stop)", ds_id, stratum)
    }
    return(invisible(FALSE))
  }
  pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) {
    if (exists("log_msg", mode = "function")) {
      log_msg("[CSN-audit-safe] %s | %s: fallback prcomp still failed, skip (do not stop)", ds_id, stratum)
    }
    return(invisible(FALSE))
  }
  varpc1 <- (pc$sdev[1]^2) / sum(pc$sdev^2)
  if (exists("log_msg", mode = "function")) {
    log_msg(
      "[CSN-audit-safe] %s | %s: fallback OK; genes=%d; PC1%%=%.1f",
      ds_id, stratum, nrow(X), 100 * varpc1
    )
  }
  invisible(TRUE)
}

# --- helper: Ensure limma t-stat has correct gene names ---
._ensure_stats_names <- function(stats, gene_names, label = NULL) {
  # stats: numeric vector (limma t etc.)
  # gene_names: rownames corresponding to matrix (Gene ID)
  if (is.null(stats)) {
    return(NULL)
  }
  v <- suppressWarnings(as.numeric(stats))
  if (length(v) != length(gene_names)) {
    stop(sprintf(
      "[ensure-names%s] stats(%d) and gene_names(%d) length mismatch",
      if (!is.null(label)) paste0("-", label) else "",
      length(v), length(gene_names)
    ))
  }
  names(v) <- as.character(gene_names)
  # Do not filter yet, leave to next step cleaning (will record log)
  v
}



## ===== [REMOVED] GSEA helper functions =====
## Removed: .gsea_from_ranks, ._finite_rank_stats, ._intersect_and_filter_pathways


.align_limma_inputs <- function(M, design_df, rhs_terms, label = "limma") {
  stopifnot(!is.null(M), is.matrix(M), ncol(M) > 0, !is.null(colnames(M)))
  # Convert design table to data.frame first, ensure rownames are sample IDs
  design_df <- as.data.frame(design_df, check.names = FALSE)
  stopifnot(!is.null(rownames(design_df)))
  samp_M <- as.character(colnames(M))
  samp_des <- as.character(rownames(design_df))

  # Kick out NA samples in design table (predictor / numeric covariates / factor NA are all considered unavailable)
  na_mask <- rep(FALSE, nrow(design_df))
  rn <- rownames(design_df)

  # predictor must exist
  if (!("predictor" %in% colnames(design_df))) {
    stop(sprintf("[%s] design_df missing predictor column", label))
  }

  # predictor non-finite
  v <- design_df$predictor
  na_mask <- na_mask | is.na(v) | !is.finite(v)

  # Other numeric columns
  num_cols <- setdiff(colnames(design_df), c("predictor"))
  for (cn in num_cols) {
    v <- design_df[[cn]]
    if (is.numeric(v)) {
      na_mask <- na_mask | is.na(v) | !is.finite(v)
    } else {
      # Factor/String/Logical: Remove if NA
      na_mask <- na_mask | is.na(v)
    }
  }
  design_df <- design_df[!na_mask, , drop = FALSE]

  # Intersect and same order
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

  # Check once more
  if (nrow(des2) != ncol(M2)) {
    msg <- sprintf(
      "[%s] nrow(des)=%d != ncol(M)=%d (Still inconsistent after alignment)",
      label, nrow(des2), ncol(M2)
    )
    stop(msg)
  }
  list(M = M2, des = des2, design_df = design_df2)
}

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
  # Extra: List top 6 comparisons, convenient for quick look in log
  head_show <- utils::head(df[, c("pos", "M_sample", "interest_rn", "nuisance_rn", "interest_ok", "nuisance_ok")], 6)
  log_msg("[align:%s] head:\n%s", tag, utils::capture.output(print(head_show)) |> paste(collapse = "\n"))
  invisible(df)
}


## ===== impute_and_filter: Fix seed before imputation (5) =====
impute_and_filter <- function(mat, min_frac = 0.75) {
  keep <- rowMeans(!is.na(mat)) >= min_frac
  m <- mat[keep, , drop = FALSE]
  if (any(is.na(m))) {
    set.seed(1234)
    m <- imputeLCMD::impute.MinProb(m, q = 0.01)
  }
  m
}




## ===== Coverage and covariates audit: audit_covars_coverage() =====
## Suggested placement before run_one_stratum()
## Will append summary to run_info/covars_audit/audit_rows.csv, and print a log line
if (!exists(".ensure_mat_or_null", mode = "function")) {
  .ensure_mat_or_null <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    m <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (is.null(m)) {
      return(NULL)
    }
    if (is.null(dim(m))) { # Vector -> Single column matrix
      m <- matrix(m, ncol = 1)
    }
    if (ncol(m) == 0) {
      return(NULL)
    }
    m
  }
}




## ===== [REMOVED] GSEA ranking and execution functions =====
## Removed: GSEA ranking and execution functions


## ===== [REMOVED] GSEA summarization functions =====
## Removed: read_gsea_table, merge_subunit_tables, add_sig_counts, write_summary_outputs, summarize_all_groups

## ===== [KEPT] DEG summarization (for cross-predictor summary) =====
summarize_deg_across_predictors <- function(out_root) {
  log_msg("[Summary-DEG] Root directory: %s", out_root)
  sfn <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x)) # Safe column name

  ver <- "BatchAdj"
  base_dir <- file.path(out_root, "DEG", ver)
  if (!dir.exists(base_dir)) {
    log_msg("  [DEG-summary] %s | Directory not found (%s), skip", ver, base_dir)
    return(invisible(NULL))
  }

  preds <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
  if (!length(preds)) {
    log_msg("  [DEG-summary] %s | No predictor subdirectories, skip", ver)
    return(invisible(NULL))
  }

  tabs <- list()
  for (p in preds) {
    fp <- file.path(base_dir, p, "DEG_limma_cont.csv")
    if (!file.exists(fp)) next
    tb <- tryCatch(data.table::fread(fp, showProgress = FALSE), error = function(e) NULL)
    if (is.null(tb) || !nrow(tb)) next

    need <- c("gene", "logFC", "adj.P.Val")
    if (!all(need %in% names(tb))) next

    tb_slim <- tb[, ..need]
    setnames(tb_slim,
      old = c("gene", "logFC", "adj.P.Val"),
      new = c(
        "gene",
        paste0("logFC_", sfn(p)),
        paste0("padj_", sfn(p))
      )
    )
    tabs[[p]] <- tb_slim
  }

  if (!length(tabs)) {
    log_msg("  [DEG-summary] %s | Cannot find any available DEG_limma_cont.csv, skip", ver)
    return(invisible(NULL))
  }

  df <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), tabs)
  padj_cols <- grep("^padj_", names(df), value = TRUE)
  if (length(padj_cols)) {
    suppressWarnings({
      df$min_padj <- apply(as.data.frame(df[, ..padj_cols]), 1, function(z) {
        z <- suppressWarnings(as.numeric(z))
        if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE)
      })
    })
    data.table::setcolorder(df, c(
      "gene", "min_padj",
      setdiff(names(df), c("gene", "min_padj"))
    ))
    df <- df[order(df$min_padj), ]
  }

  out_dir <- file.path(out_root, "summary", "DEG", ver)
  out_dir <- file.path(out_root, "DEG", ver)
  base <- file.path(out_dir, "Summary_DEG_wide")

  data.table::fwrite(df, paste0(base, ".csv"))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "DEG_wide")
  openxlsx::writeData(wb, "DEG_wide", df)
  openxlsx::saveWorkbook(wb, paste0(base, ".xlsx"), overwrite = TRUE)

  log_msg(
    "[DEG-summary] %s | predictors=%d | genes=%d: Wrote %s",
    ver, length(tabs), nrow(df), paste0(base, ".csv")
  )
}


## ===== [REMOVED] Meta-analysis functions for GSEA =====
## Removed: safe_read_gsea, .z_from_p_signed, meta_fdr_stouffer


## ==== PATCH: helpers ====

## === NEW: CSN complex score (PC1) & residualization helpers ==================

# Perform PCA on "cross-sample z-scores" of subunits, take PC1 as CSN score; correct direction to same sign as mean z
build_csn_score <- function(mat0,
                            subunits = csn_subunits,
                            combine_7AB = TRUE,
                            min_members = 5L) {
  present <- intersect(subunits, rownames(mat0))
  # Pre-build return skeleton (Must have names)
  s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
  if (!length(present)) {
    return(s)
  }

  # z-score (Keep sample names)
  get_z <- function(v) {
    nm <- names(v) # Save sample names first
    v <- as.numeric(v)
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    out <- (v - mu) / sdv
    names(out) <- nm # Put back sample names
    out
  }

  X <- do.call(rbind, lapply(present, function(g) get_z(mat0[g, ])))
  rownames(X) <- present
  colnames(X) <- colnames(mat0) # This line is important: Ensure sample names exist

  # Optional: Combine COPS7A/7B
  if (combine_7AB && all(c("COPS7A", "COPS7B") %in% rownames(X))) {
    Z7 <- colMeans(X[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
    X <- rbind(X[setdiff(rownames(X), c("COPS7A", "COPS7B")), , drop = FALSE],
      "COPS7*" = Z7
    )
  }

  # Use non-NA count of "original mat0" to judge coverage
  enough <- colSums(is.finite(mat0[present, , drop = FALSE])) >= min_members
  keep_sam <- names(s)[enough]

  if (length(keep_sam) >= 10) {
    pc <- stats::prcomp(t(X[, keep_sam, drop = FALSE]), center = TRUE, scale. = FALSE)
    sc <- pc$x[, 1]

    # Direction correction: Ensure PC1 has same sign as subunit mean z-score
    mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)

    # Critical fix: Check if correlation is valid before using it
    rho <- suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs"))

    if (is.finite(rho) && rho < 0) {
      sc <- -sc
    } else if (!is.finite(rho)) {
      # If correlation cannot be computed (e.g., all subunits NA for some samples),
      # keep original PC1 direction but log a warning
      if (exists("log_msg", mode = "function")) {
        try(log_msg(
          "[CSN_SCORE] Warning: Cannot compute cor(PC1, mean_z) for direction correction (rho=%s). Using uncorrected PC1.",
          as.character(rho)
        ), silent = TRUE)
      }
    }

    s[keep_sam] <- sc # Here alignment by sample name will definitely succeed
  }

  s
}

# Regress vector y (a subunit) on CSN score + batch + other covariates, then take residuals
# Replace original residualize_vector()
residualize_vector <- function(y, csn_score, batch = NULL, covars = NULL, min_n = 8L) {
  # 1) Align sample names
  if (is.null(names(y))) stop("[residualize_vector] y must be a numeric vector with sample names")
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

  # 2) Build data frame -> Single model.matrix expansion
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
    C <- coerce_covariates_safely(C) # <- New
    for (cn in colnames(C)) DF[[cn]] <- C[[cn]] # <- Put back to DF directly
  }
  DF_y <- suppressWarnings(as.numeric(y[common]))

  ## === NEW: Clear columns that would wipe out complete.cases first ===
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

  # 2b) Drop single-level factor or constant numeric columns (Avoid design matrix singularity)
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
  # - If batch exists and is factor, auto dummy encode
  des <- stats::model.matrix(~ 1 + ., data = DF, na.action = stats::na.fail)
  stopifnot(nrow(des) == length(y_ok))

  # 4) Linear regression y ~ csn + batch + covars
  fit <- lm.fit(x = des, y = y_ok)

  # Critical fix: Ensure res has sample names for safe alignment
  # Create named vector aligned with 'common' samples
  res <- setNames(rep(NA_real_, length(DF_y)), common)

  # Extract sample names from DF (which were filtered by ok)
  kept_samples <- rownames(DF)
  res[kept_samples] <- fit$residuals

  # 5) Map back to full samples and (optional) z-score
  out <- setNames(rep(NA_real_, length(y)), names(y))
  out[common] <- res # Now uses name-based matching, not positional
  fin <- is.finite(out)
  if (sum(fin) >= 3L) {
    mu <- mean(out[fin])
    sdv <- stats::sd(out[fin])
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    out[fin] <- (out[fin] - mu) / sdv
  }
  out
}

audit_covars_coverage <- function(tag, ds_id, stratum, su,
                                  sample_ids,
                                  batch = NULL,
                                  covars = NULL,
                                  sv = NULL,
                                  tech = NULL) {
  dir.create(file.path("run_info", "covars_audit"), recursive = TRUE, showWarnings = FALSE)

  cov_df <- if (is.null(covars)) NULL else as.data.frame(covars, check.names = FALSE)
  cov_nms <- if (!is.null(cov_df)) colnames(cov_df) else character(0)

  get_cov <- function(k) {
    if (is.null(cov_df) || !(k %in% names(cov_df))) {
      return(NA_real_)
    }
    v <- cov_df[[k]]
    if (is.numeric(v)) {
      mean(is.finite(v)) * 100
    } else {
      mean(!is.na(v)) * 100 # Factor/String -> Count by non-NA ratio
    }
  }
  cov_purity <- get_cov("purity")
  cov_sex <- get_cov("sex")
  cov_age <- get_cov("age")

  sv_nms <- if (!is.null(sv)) colnames(.ensure_mat_or_null(sv)) else NULL
  tech_nms <- if (!is.null(tech)) colnames(.ensure_mat_or_null(tech)) else NULL

  line <- sprintf(
    "  [audit:%s] covars={%s} | coverage: purity=%.1f%%, sex=%.1f%%, age=%.1f%% | batch_levels=%s | SV={%s} | tech={%s}",
    tag,
    if (length(cov_nms)) paste(cov_nms, collapse = ",") else "NULL",
    cov_purity %||% NaN, cov_sex %||% NaN, cov_age %||% NaN,
    if (!is.null(batch)) nlevels(batch) else 0L,
    if (length(sv_nms)) paste(sv_nms, collapse = ",") else "none",
    if (length(tech_nms)) paste(tech_nms, collapse = ",") else "none"
  )
  log_msg(line)

  # Append a line to summary table
  df <- data.frame(
    dataset = ds_id, stratum = stratum, subunit = su, tag = tag,
    batch_levels = if (!is.null(batch)) nlevels(batch) else 0L,
    batch_sizes = if (!is.null(batch)) {
      paste(sprintf(
        "%s=%d", names(sort(table(batch), decreasing = TRUE)),
        as.integer(sort(table(batch), decreasing = TRUE))
      ), collapse = "; ")
    } else {
      NA_character_
    },
    covars_cols = if (length(cov_nms)) paste(cov_nms, collapse = ";") else "NULL",
    sv_cols = if (length(sv_nms)) paste(sv_nms, collapse = ";") else "NULL",
    tech_cols = if (length(tech_nms)) paste(tech_nms, collapse = ";") else "NULL",
    purity_cov = cov_purity, sex_cov = cov_sex, age_cov = cov_age,
    stringsAsFactors = FALSE
  )
  fp <- file.path("run_info", "covars_audit", "audit_rows.csv")
  data.table::fwrite(df, fp, append = file.exists(fp))
}

orthogonalize_to <- function(mat, nuisance) {
  if (is.null(mat) || is.null(nuisance)) {
    return(mat)
  }
  Y <- as.matrix(mat)
  X <- as.matrix(nuisance)

  # 1) Add intercept
  X <- cbind("(Intercept)" = rep(1, nrow(X)), X)

  # 2) Non-finite value handling + Remove zero variance columns
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


## ---- Auto detect batch column & read (Expand candidate columns) ----
## Small tool (Only used in this file)
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
  # If named vector/factor/data.frame is provided, align by colnames(mat); otherwise treat as aligned
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
  # Optional: Collapse extremely small batches (Default collapse singletons, set 0 to not collapse)
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

## ===== Batch cleaning settings (You can adjust yourself) =====
BATCH_PIPE_POLICY <- "NA" # "NA" or "b_small": Should values containing '|' be set to NA or merged to b_small
BATCH_MIN_PER_LEVEL <- 2 # Levels smaller than this threshold will be merged to b_small (Avoid inestimable coefficients)

## ---- Batch value cleaning: Handle '|', legalize names, merge sparse levels ----
sanitize_batch_levels <- function(x,
                                  pipe_policy = BATCH_PIPE_POLICY,
                                  min_per_level = BATCH_MIN_PER_LEVEL) {
  x0 <- as.character(x)
  has_pipe <- grepl("\\|", x0 %||% "")
  if (any(has_pipe)) {
    log_msg("  [batch] Detected {sum(has_pipe)} values containing '|' -> Handle by policy {pipe_policy}")
    x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
  }
  # Legalize names (Avoid eBayes/design matrix column name issues)
  fac <- factor(make.names(x0))
  fac <- droplevels(fac)

  # Merge sparse levels (e.g., only 1 sample)
  if (!is.null(min_per_level) && min_per_level > 1) {
    tab <- table(fac, useNA = "no")
    small <- names(tab)[tab < min_per_level]
    if (length(small)) {
      log_msg(
        "  [batch] Merge sparse level to 'b_small': %s",
        paste(sprintf("%s(n=%d)", small, as.integer(tab[small])), collapse = ", ")
      )
      fac_chr <- as.character(fac)
      fac_chr[fac_chr %in% small] <- "b_small"
      fac <- droplevels(factor(fac_chr))
    }
  }
  fac
}

## ---- Auto detect batch column (Add cleaning pipeline) ----
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

## ---- Get batch factor by sample order (Cleaned) ----
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

  # Align by sample_ids first, ensure detected and returned vector lengths are consistent
  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
  rownames(meta) <- sample_ids


  det <- detect_batch_column(meta,
    pipe_policy   = pipe_policy,
    min_per_level = min_per_level
  )
  if (is.null(det)) {
    ## === [NEW] Fallback: Infer TMT-plex from TMT_protein.csv as batch ===
    tmt_fp <- file.path(ds_dir, "TMT_protein.csv")
    if (file.exists(tmt_fp)) {
      tmt <- suppressMessages(readr::read_csv(tmt_fp, show_col_types = FALSE)) |> as.data.frame()

      ## Normalize column names (Trim whitespace; handle " Run Metadata ID")
      cn <- names(tmt)
      cn_trim <- trimws(cn)
      names(tmt) <- cn_trim
      run_hits <- grep("^Run\\s*Metadata\\s*ID$", cn_trim, ignore.case = TRUE, value = TRUE)
      tmt_cols <- grep("^tmt_", cn_trim, ignore.case = TRUE, value = TRUE)

      if (length(run_hits) >= 1 && length(tmt_cols) >= 1) {
        run_col <- run_hits[1]

        ## Build sample_id -> plex mapping (One plex per row; all samples in tmt_* columns of that row belong to this plex)
        plex_by_sample <- list()
        nR <- nrow(tmt)
        for (i in seq_len(nR)) {
          run_id <- as.character(tmt[[run_col]][i])
          if (!nzchar(run_id) || is.na(run_id)) next
          plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id) # Take first two digits
          if (!nzchar(plex2) || is.na(plex2)) next

          for (tc in tmt_cols) {
            cell <- tmt[[tc]][i]
            if (is.na(cell)) next
            cell <- as.character(cell)
            if (!nzchar(cell)) next
            sid <- sub("\\r?\\n.*$", "", cell) # Take before newline
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
          ## Success condition: At least two valid levels and valid sample count >= 3
          if (nlevels(fac2) >= 2 && sum(!is.na(fac2)) >= 3) {
            names(fac2) <- sample_ids
            return(list(name = "TMT_protein.csv:RunMetadataID", fac = fac2))
          }
        }
      }
    }
    ## Fallback also failed -> Return NULL (Maintain original semantics)
    return(NULL)
  }

  fac <- det$fac # Aligned with sample_ids
  names(fac) <- sample_ids
  list(name = det$name, fac = fac)
}

## ---- Batch requirement screening (Follow your logic) ----
screen_batch_need <- function(ds_dir, min_frac_complete = 0.75) {
  log_msg("== Batch Check: %s ==", basename(ds_dir))
  mat0 <- load_matrix_from_dataset_dir(ds_dir)
  mx <- suppressWarnings(max(mat0, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)
  m <- impute_and_filter(mat0, min_frac = min_frac_complete)

  bi <- get_batch_factor(ds_dir, colnames(m))
  if (is.null(bi)) {
    log_msg("  [batch] No clear batch column found -> Do not adjust for now")
    return(invisible(list(
      dataset = basename(ds_dir),
      batch_col = NA, n_levels = 0, sizes = NA,
      pc_R2 = NA, pc_p = NA, frac_FDR05 = NA, frac_FDR25 = NA,
      recommend = FALSE
    )))
  }

  batch <- droplevels(bi$fac[colnames(m)])
  tab <- sort(table(batch), decreasing = TRUE)
  log_msg(
    "  [batch] Column: %s; levels=%d; n per level: %s",
    bi$name, nlevels(batch),
    paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")
  )

  X <- t(scale(t(m), center = TRUE, scale = TRUE))
  pc <- prcomp(t(X), scale. = FALSE)
  K <- min(5, ncol(pc$x))
  r2 <- vapply(seq_len(K), function(i) summary(lm(pc$x[, i] ~ batch))$r.squared, numeric(1))
  pv <- vapply(seq_len(K), function(i) {
    a <- anova(lm(pc$x[, i] ~ batch))
    as.numeric(a$`Pr(>F)`[1])
  }, numeric(1))
  log_msg(
    "  PCA (by batch) R2: %s ; p: %s", paste(round(r2, 3), collapse = ", "),
    paste(signif(pv, 3), collapse = ", ")
  )

  design <- model.matrix(~batch)
  fit <- limma::lmFit(m, design)
  fit <- limma::eBayes(fit)
  padj <- p.adjust(fit$F.p.value, "BH")
  prop05 <- mean(padj < 0.05, na.rm = TRUE)
  prop25 <- mean(padj < 0.25, na.rm = TRUE)
  log_msg(
    "  Gene-level F-test: FDR<0.05 prop = %.1f%%; FDR<0.25 = %.1f%%",
    100 * prop05, 100 * prop25
  )

  recommend <- (any(r2 >= 0.10 & pv[seq_along(r2)] < 0.01) || prop05 >= 0.05)
  if (recommend) {
    log_msg("  **Suggest Correction**: R2 or gene proportion met threshold (>=10%% R2 and p<0.01, or FDR<0.05 genes >=5%%)")
  } else {
    log_msg("  **Can Skip Correction**: No obvious batch effect (Record kept)")
  }

  invisible(list(
    dataset = basename(ds_dir),
    batch_col = bi$name, n_levels = nlevels(batch),
    sizes = tab, pc_R2 = r2, pc_p = pv,
    frac_FDR05 = prop05, frac_FDR25 = prop25,
    recommend = recommend
  ))
}

datasets_root <- getwd()

# Dataset IDs for proteomic DEG analysis (7 CPTAC datasets)
dataset_ids <- c(
  "brca_cptac_2020", "luad_cptac_2020", "lusc_cptac_2021", "ucec_cptac_2020",
  "coad_cptac_2019", "gbm_cptac_2021", "paad_cptac_2021"
)

# Output prefix for all results
OUTPUT_PREFIX <- "proteomic_DEG"

dataset_dirs <- setNames(file.path(datasets_root, dataset_ids), dataset_ids)
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")
message("datasets_root = ", datasets_root)

## Old version was stopifnot(all(dir.exists(dataset_dirs))), which would stop batch if some directories missing.
## Changed to: List missing ones, filter "actually runnable this round" set with dataset_dirs_run.
missing_dirs <- names(dataset_dirs)[!dir.exists(dataset_dirs)]
if (length(missing_dirs)) {
  log_msg("Detected %d missing folders, skipping: %s", length(missing_dirs), paste(missing_dirs, collapse = ", "))
}


## List of datasets actually to run and available this round (Folder exists and contains data_protein_quantification.txt)
dataset_dirs_run <- dataset_dirs[
  dir.exists(dataset_dirs) &
    file.exists(file.path(dataset_dirs, "data_protein_quantification.txt"))
]

if (!length(dataset_dirs_run)) stop("dataset_dirs_run is empty, please check if folders and data_protein_quantification.txt exist")

log_msg("Available datasets this round (%d): %s", length(dataset_dirs_run), paste(names(dataset_dirs_run), collapse = ", "))


## ---- Helper: Align any vector/data frame to specified sample order ----
.align_to_samples <- function(x, sam, what = "covariate") {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.vector(x) || is.factor(x)) {
    # Allow length exactly equal to sam and no names; otherwise prioritize alignment by names
    if (is.null(names(x))) {
      if (length(x) != length(sam)) stop(sprintf("[%s] Length mismatch: %d vs %d", what, length(x), length(sam)))
      names(x) <- sam
    }
    if (!all(sam %in% names(x))) stop(sprintf("[%s] Missing samples: %s", what, paste(setdiff(sam, names(x)), collapse = ", ")))
    out <- x[sam]
    return(out)
  } else {
    x <- as.data.frame(x)
    # Expect rownames(x) to be samples
    if (is.null(rownames(x))) {
      if (nrow(x) != length(sam)) stop(sprintf("[%s] Row count mismatch: %d vs %d and no rownames for alignment", what, nrow(x), length(sam)))
      rownames(x) <- sam
    }
    if (!all(sam %in% rownames(x))) stop(sprintf("[%s] Missing samples: %s", what, paste(setdiff(sam, rownames(x)), collapse = ", ")))
    out <- x[sam, , drop = FALSE]
    return(out)
  }
}

# Any object -> Matrix; Force set rownames; Return NULL if 0 columns or row count mismatch specified rows
.mk_mat_or_null <- function(x, rows) {
  if (is.null(x)) {
    return(NULL)
  }
  m <- as.matrix(x)
  if (is.null(rownames(m))) rownames(m) <- rows
  if (nrow(m) != length(rows) || ncol(m) == 0) {
    return(NULL)
  }
  m
}
# More generic safe matrix conversion (Allow vector, 1 column, or zero length; no specified rows)
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


## Safely convert vector to factor, and explicitly convert NA to "NA" (then add prefix)
.factorize_with_explicit_NA <- function(x) {
  x_chr <- as.character(x)
  x_chr[!nzchar(x_chr) | is.na(x_chr)] <- "NA" # Empty string/NA treated as "NA"
  factor(x_chr)
}

## ===== Get TP53 mutation status =====
## ===== TP53 status (protein-altering baseline) =====
## Only count mutations that alter protein as TP53-mutant; others (Silent, UTR, Intron, IGR, RNA, lincRNA, Flank...) treated as wild type

TP53_KEEP_CLASSES <- c(
  "MISSENSE_MUTATION", "NONSENSE_MUTATION",
  "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
  "IN_FRAME_DEL", "IN_FRAME_INS",
  "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
)

normalize_vc <- function(x) {
  # Normalize Variant_Classification: To upper case, unify various separators to underscore
  x <- toupper(trimws(as.character(x)))
  gsub("[^A-Z0-9]+", "_", x)
}

get_tp53_status <- function(ds_dir, sample_ids) {
  # Default all wild-type
  status <- setNames(rep("TP53_wild_type", length(sample_ids)), sample_ids)

  mut_fp <- file.path(ds_dir, "data_mutations.txt")
  if (!file.exists(mut_fp)) {
    log_msg("  [TP53] data_mutations.txt not found, treat all as wild-type/ALL available")
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
    log_msg("  [TP53] data_mutations.txt missing columns: %s -> Treat as wild-type", paste(miss, collapse = ", "))
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
##      (Original Variant_Classification -> "Sample count"; also includes Any_TP53_mutation and protein_altering)
##   - run_info/tp53_status/tp53_class_sample_counts_long.csv (Aggregated long format)
##   - run_info/tp53_status/tp53_binary_counts_by_dataset.csv (wild type / mutant summary table)

dir.create(file.path("run_info", "tp53_status"), recursive = TRUE, showWarnings = FALSE)

# Insurance: If not defined earlier, add a copy here
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

# TP53 sample count summary for a single dataset
summarize_tp53_counts_for_dataset <- function(ds_dir) {
  ds_id <- basename(ds_dir)

  # Use protein matrix samples as "population"
  M <- try(load_matrix_from_dataset_dir(ds_dir), silent = TRUE)
  if (inherits(M, "try-error")) {
    log_msg("[TP53-audit] %s: Cannot read matrix, skip", ds_id)
    return(NULL)
  }
  sample_ids <- colnames(M)
  sid_up <- toupper(sample_ids)
  n_all <- length(sample_ids)

  # Get mutant/wild type (protein-altering definition)
  status <- get_tp53_status(ds_dir, sample_ids)
  tb_bin <- table(factor(status, levels = c("TP53_wild_type", "TP53_mutant")))
  bin_row <- data.frame(
    dataset = ds_id,
    in_matrix_n = n_all,
    WT_n = as.integer(tb_bin["TP53_wild_type"]),
    MUT_n = as.integer(tb_bin["TP53_mutant"]),
    stringsAsFactors = FALSE
  )

  # Read MAF (Original Variant_Classification distribution -> "Sample count")
  mut_fp <- file.path(ds_dir, "data_mutations.txt")
  class_df <- NULL
  if (file.exists(mut_fp)) {
    mutation_df <- try(mutation_df <- readr::read_tsv(mut_fp, comment = "#", show_col_types = FALSE), silent = TRUE)
    if (!inherits(mutation_df, "try-error") &&
      all(c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode") %in% colnames(mutation_df))) {
      tp53 <- subset(mutation_df, Hugo_Symbol == "TP53",
        select = c("Variant_Classification", "Tumor_Sample_Barcode")
      )

      if (nrow(tp53) > 0) {
        tp53$Variant_Classification <- normalize_vc(tp53$Variant_Classification)
        tp53$Tumor_Sample_Barcode <- toupper(tp53$Tumor_Sample_Barcode)

        # Only count samples "in protein matrix"
        tp53 <- tp53[tp53$Tumor_Sample_Barcode %in% sid_up, , drop = FALSE]

        # Each classification -> How many "unique samples" hit
        class_counts <- tp53 |>
          dplyr::group_by(Variant_Classification) |>
          dplyr::summarise(sample_n = dplyr::n_distinct(Tumor_Sample_Barcode), .groups = "drop") |>
          dplyr::arrange(dplyr::desc(sample_n), Variant_Classification)

        # Add "Any_TP53_mutation" and "protein_altering" summary rows (sample level)
        any_samples <- unique(tp53$Tumor_Sample_Barcode)

        keep_flag <- tp53$Variant_Classification %in% TP53_KEEP_CLASSES |
          grepl("^FRAME_SHIFT_", tp53$Variant_Classification) |
          grepl("^IN_FRAME_", tp53$Variant_Classification)
        prot_alt_samples <- unique(tp53$Tumor_Sample_Barcode[keep_flag])

        add_rows <- data.frame(
          Variant_Classification = c("Any_TP53_mutation", "protein_altering"),
          sample_n = c(length(any_samples), length(prot_alt_samples)),
          stringsAsFactors = FALSE
        )

        class_df <- dplyr::bind_rows(class_counts, add_rows)
        class_df$dataset <- ds_id
        class_df <- class_df[, c("dataset", "Variant_Classification", "sample_n")]
      }
    }
  }

  # Write files individually
  if (!is.null(class_df)) {
    data.table::fwrite(
      class_df,
      file.path("run_info", "tp53_status", paste0(ds_id, "_tp53_class_sample_counts.csv"))
    )
  } else {
    # If no MAF or no TP53 record, output empty shell
    data.table::fwrite(
      data.frame(dataset = ds_id, Variant_Classification = NA_character_, sample_n = NA_integer_),
      file.path("run_info", "tp53_status", paste0(ds_id, "_tp53_class_sample_counts.csv"))
    )
  }

  list(binary = bin_row, class_long = class_df)
}

# Run all datasets, aggregate summary table
summarize_tp53_counts_all_datasets <- function(dataset_dirs) {
  all_bin <- list()
  all_class <- list()
  k <- 1L
  j <- 1L
  for (ds in names(dataset_dirs)) {
    ds_dir <- dataset_dirs[[ds]]
    if (!dir.exists(ds_dir)) next
    log_msg("[TP53-audit] Start: %s", ds)
    res <- summarize_tp53_counts_for_dataset(ds_dir)
    if (is.null(res)) next
    all_bin[[k]] <- res$binary
    k <- k + 1L
    if (!is.null(res$class_long)) {
      all_class[[j]] <- res$class_long
      j <- j + 1L
    }
  }

  if (length(all_bin)) {
    bin_df <- dplyr::bind_rows(all_bin) |>
      dplyr::mutate(Any_TP53_mutation_n = in_matrix_n - WT_n)
    data.table::fwrite(bin_df, file.path("run_info", "tp53_status", "tp53_binary_counts_by_dataset.csv"))
    log_msg("[TP53-audit] Wrote: tp53_binary_counts_by_dataset.csv")
  }
  if (length(all_class)) {
    class_df <- dplyr::bind_rows(all_class)
    data.table::fwrite(class_df, file.path("run_info", "tp53_status", "tp53_class_sample_counts_long.csv"))
    log_msg("[TP53-audit] Wrote: tp53_class_sample_counts_long.csv")
  }
}

## === Call (Can be placed before or after main flow; does not affect downstream analysis) ===
## summarize_tp53_counts_all_datasets(dataset_dirs)


## === run_deg_for_predictor (limma continuous; no imputation; batch/covariate consistent with GSEA) ===
run_deg_for_predictor <- function(
  predictor_name,
  predictor_vec, # named numeric by sample (can contain NA)
  ds_id, ds_dir,
  mat0, # RAW genes x samples; no imputation
  out_root,
  batch_all = NULL,
  purity_all = NULL,
  sa_all = NULL, # may contain age_missing, age_z_imputed
  min_per_group = NULL
) {
  log_msg <- get0("log_msg", ifnotfound = function(.) cat(sprintf("[LOG] %s\n", sprintf(.))))
  # [FIX] Backfill with global settings to avoid hitting local NULL (inherits=TRUE will hit parameter layer first)
  if (is.null(min_per_group) || length(min_per_group) == 0L || !is.finite(min_per_group)) {
    if (exists("min_per_group", envir = .GlobalEnv, inherits = FALSE)) {
      min_per_group <- get("min_per_group", envir = .GlobalEnv, inherits = FALSE)
    } else {
      min_per_group <- 8L
    }
  }
  min_per_group <- as.integer(min_per_group)[1L]
  if (!is.finite(min_per_group) || min_per_group < 1L) min_per_group <- 8L

  ## 1) Align samples (Based on predictor
  if (is.null(names(predictor_vec))) names(predictor_vec) <- colnames(mat0)
  keep <- intersect(colnames(mat0), names(predictor_vec))
  pred_all <- suppressWarnings(as.numeric(predictor_vec[keep]))
  names(pred_all) <- keep
  fin <- is.finite(pred_all)
  if (sum(fin) < (2L * min_per_group)) {
    log_msg("  [DEG-%s] predictor non-NA samples insufficient (%d < %d), skip", predictor_name, sum(fin), 2L * min_per_group)
    return(invisible(NULL))
  }
  sample_order <- keep[fin]
  pred <- pred_all[fin]
  M0 <- as.matrix(mat0[, sample_order, drop = FALSE])

  ## 2) Covariates (Consistent with GSEA; no forced imputation)
  USE_AGE_MISSING_INDICATOR <- isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))
  build_covars_df <- function(so) {
    df <- data.frame(row.names = so, check.names = FALSE)
    if (!is.null(purity_all)) df$purity <- suppressWarnings(as.numeric(purity_all[so]))
    if (!is.null(sa_all)) {
      if ("sex" %in% colnames(sa_all)) df$sex <- factor(sa_all[so, "sex"])
      if ("age" %in% colnames(sa_all)) df$age <- suppressWarnings(as.numeric(sa_all[so, "age"]))
      if (USE_AGE_MISSING_INDICATOR) {
        if ("age_missing" %in% colnames(sa_all)) df$age_missing <- suppressWarnings(as.numeric(sa_all[so, "age_missing"]))
        if ("age_z_imputed" %in% colnames(sa_all)) df$age_z_imputed <- suppressWarnings(as.numeric(sa_all[so, "age_z_imputed"]))
      }
    }
    df
  }
  df_covars0 <- build_covars_df(sample_order)

  ## 3) BatchAdj version: Add batch/covars to design formula
  DF_ba <- df_covars0
  if (!is.null(batch_all)) DF_ba$batch <- factor(batch_all[sample_order])

  ok <- intersect(colnames(M0), rownames(DF_ba))
  DF_ba <- DF_ba[ok, , drop = FALSE]
  ## Clean covariates: Normalize types + Remove single-level factors (Avoid contrasts error)
  DF_ba <- coerce_covariates_safely(DF_ba)
  pred_ba <- pred[ok]
  names(pred_ba) <- ok
  M <- M0[, ok, drop = FALSE]

  DF_ba$predictor <- pred_ba - mean(pred_ba, na.rm = TRUE)
  ## Remove all-NA covariate columns, and filter out single-level factors (Avoid model.matrix dropping all samples)
  if (ncol(DF_ba)) {
    all_na <- vapply(DF_ba, function(z) all(is.na(z)), logical(1))
    if (any(all_na)) DF_ba <- DF_ba[, !all_na, drop = FALSE]
  }
  if (ncol(DF_ba)) {
    DF_ba <- coerce_covariates_safely(DF_ba)
  }
  # Ensure predictor is still there (coerce won't touch numeric columns, but just in case)
  if (!("predictor" %in% colnames(DF_ba))) {
    DF_ba$predictor <- pred_ba - mean(pred_ba, na.rm = TRUE)
  }
  des <- stats::model.matrix(~ 1 + ., data = DF_ba) # Auto expand batch/covars
  ## Align expression matrix to samples actually kept in design matrix (model.matrix might drop rows due to NA)
  use <- rownames(des)
  M <- M[, use, drop = FALSE]
  pass_label <- "BatchAdj"

  ## DEG safety: Filter out genes with insufficient degrees of freedom based on current design matrix (Avoid NA in eBayes weights)
  rnk <- qr(des)$rank
  need <- rnk + 1L # At least 1 residual degree of freedom
  nobs <- rowSums(is.finite(M))
  keep_rows <- nobs >= need

  if (!all(keep_rows)) {
    log_msg(
      "  [DEG|%s] Remove %d/%d genes (Valid samples < rank(des)+1 = %d; avoid NA df/weights)",
      pass_label, sum(!keep_rows), nrow(M), need
    )
    M <- M[keep_rows, , drop = FALSE]
  }
  if (nrow(M) == 0L) {
    log_msg("  [DEG|%s] No available genes after filtering, skip", pass_label)
    return(invisible(NULL))
  }

  fit <- limma::lmFit(M, design = des) # Allow NA
  eb <- limma::eBayes(fit, trend = TRUE)
  coef_name <- if ("predictor" %in% colnames(des)) "predictor" else tail(colnames(des), 1)
  tbl <- limma::topTable(eb, coef = coef_name, number = Inf, sort.by = "P")
  tbl$gene <- rownames(tbl)
  tbl <- tbl[, c("gene", "logFC", "t", "P.Value", "adj.P.Val", "B")]

  out_dir <- file.path(out_root, "DEG", pass_label, predictor_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(tbl, file.path(out_dir, "DEG_limma_cont.csv"))
  log_msg(
    "[DEG] %s | %s | %s: Wrote %s",
    ds_id, basename(out_root), predictor_name,
    file.path(out_dir, "DEG_limma_cont.csv")
  )
  invisible(NULL)
}


## === run_deg_interaction_model: One-step model for subunit controlling for CSN_SCORE ===
## This function implements the statistically correct one-step approach instead of
## the problematic two-step residualization. It includes both subunit and CSN_SCORE
## in the same limma model, allowing us to test the subunit's partial effect while
## controlling for complex activity.
##
## Biological interpretation:
##   The coefficient for 'subunit' represents the effect of individual subunit
##   variation that is independent of CSN complex activity (same as residualization
##   but statistically correct).
##
run_deg_interaction_model <- function(
  subunit_name,
  subunit_vec, # Subunit expression (named numeric by sample)
  csn_score_vec, # CSN complex score (named numeric by sample)
  ds_id, ds_dir,
  mat0, # RAW genes x samples; no imputation
  out_root,
  batch_all = NULL,
  purity_all = NULL,
  sa_all = NULL,
  min_per_group = NULL
) {
  log_msg <- get0("log_msg", ifnotfound = function(.) cat(sprintf("[LOG] %s\n", sprintf(.))))

  # Backfill min_per_group
  if (is.null(min_per_group) || length(min_per_group) == 0L || !is.finite(min_per_group)) {
    if (exists("min_per_group", envir = .GlobalEnv, inherits = FALSE)) {
      min_per_group <- get("min_per_group", envir = .GlobalEnv, inherits = FALSE)
    } else {
      min_per_group <- 8L
    }
  }
  min_per_group <- as.integer(min_per_group)[1L]
  if (!is.finite(min_per_group) || min_per_group < 1L) min_per_group <- 8L

  ## 1) Align samples (Based on both subunit and CSN_SCORE)
  if (is.null(names(subunit_vec))) names(subunit_vec) <- colnames(mat0)
  if (is.null(names(csn_score_vec))) names(csn_score_vec) <- colnames(mat0)

  # Samples must have both subunit and CSN_SCORE
  common <- intersect(colnames(mat0), names(subunit_vec))
  common <- intersect(common, names(csn_score_vec))

  subunit_all <- suppressWarnings(as.numeric(subunit_vec[common]))
  csn_all <- suppressWarnings(as.numeric(csn_score_vec[common]))
  names(subunit_all) <- common
  names(csn_all) <- common

  # Both must be finite
  fin <- is.finite(subunit_all) & is.finite(csn_all)
  if (sum(fin) < (2L * min_per_group)) {
    log_msg(
      "  [DEG-INTERACTION-%s] Insufficient samples with both subunit and CSN_SCORE (%d < %d), skip",
      subunit_name, sum(fin), 2L * min_per_group
    )
    return(invisible(NULL))
  }

  sample_order <- common[fin]
  subunit <- subunit_all[fin]
  csn_score <- csn_all[fin]
  M0 <- as.matrix(mat0[, sample_order, drop = FALSE])

  ## 2) Build covariates
  USE_AGE_MISSING_INDICATOR <- isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))
  build_covars_df <- function(so) {
    df <- data.frame(row.names = so, check.names = FALSE)
    if (!is.null(purity_all)) df$purity <- suppressWarnings(as.numeric(purity_all[so]))
    if (!is.null(sa_all)) {
      if ("sex" %in% colnames(sa_all)) df$sex <- factor(sa_all[so, "sex"])
      if ("age" %in% colnames(sa_all)) df$age <- suppressWarnings(as.numeric(sa_all[so, "age"]))
      if (USE_AGE_MISSING_INDICATOR) {
        if ("age_missing" %in% colnames(sa_all)) df$age_missing <- suppressWarnings(as.numeric(sa_all[so, "age_missing"]))
        if ("age_z_imputed" %in% colnames(sa_all)) df$age_z_imputed <- suppressWarnings(as.numeric(sa_all[so, "age_z_imputed"]))
      }
    }
    df
  }
  df_covars0 <- build_covars_df(sample_order)

  ## 3) ONE-STEP MODEL: Include both subunit and CSN_SCORE
  DF_ba <- df_covars0
  if (!is.null(batch_all)) DF_ba$batch <- factor(batch_all[sample_order])

  ok <- intersect(colnames(M0), rownames(DF_ba))
  DF_ba <- DF_ba[ok, , drop = FALSE]
  DF_ba <- coerce_covariates_safely(DF_ba)

  subunit_use <- subunit[ok]
  csn_use <- csn_score[ok]
  M <- M0[, ok, drop = FALSE]

  # Add both predictors (centered)
  DF_ba$subunit <- subunit_use - mean(subunit_use, na.rm = TRUE)
  DF_ba$CSN_SCORE <- csn_use - mean(csn_use, na.rm = TRUE)

  # Remove all-NA and single-level columns
  if (ncol(DF_ba)) {
    all_na <- vapply(DF_ba, function(z) all(is.na(z)), logical(1))
    if (any(all_na)) DF_ba <- DF_ba[, !all_na, drop = FALSE]
  }
  if (ncol(DF_ba)) {
    DF_ba <- coerce_covariates_safely(DF_ba)
  }

  # Ensure both predictors are still there
  if (!("subunit" %in% colnames(DF_ba))) {
    DF_ba$subunit <- subunit_use - mean(subunit_use, na.rm = TRUE)
  }
  if (!("CSN_SCORE" %in% colnames(DF_ba))) {
    DF_ba$CSN_SCORE <- csn_use - mean(csn_use, na.rm = TRUE)
  }

  # Build design matrix with BOTH subunit and CSN_SCORE
  des <- stats::model.matrix(~ 1 + ., data = DF_ba)
  use <- rownames(des)
  M <- M[, use, drop = FALSE]
  pass_label <- "BatchAdj"

  ## DEG safety: Filter genes with insufficient df
  rnk <- qr(des)$rank
  need <- rnk + 1L
  nobs <- rowSums(is.finite(M))
  keep_rows <- nobs >= need

  if (!all(keep_rows)) {
    log_msg(
      "  [DEG-INTERACTION|%s] Remove %d/%d genes (df insufficient)",
      pass_label, sum(!keep_rows), nrow(M)
    )
    M <- M[keep_rows, , drop = FALSE]
  }
  if (nrow(M) == 0L) {
    log_msg("  [DEG-INTERACTION|%s] No available genes, skip", pass_label)
    return(invisible(NULL))
  }

  ## Run limma
  fit <- limma::lmFit(M, design = des)
  eb <- limma::eBayes(fit, trend = TRUE)

  # Test the subunit coefficient (controlling for CSN_SCORE)
  coef_name <- if ("subunit" %in% colnames(des)) "subunit" else "subunit"
  tbl <- limma::topTable(eb, coef = coef_name, number = Inf, sort.by = "P")
  tbl$gene <- rownames(tbl)
  tbl <- tbl[, c("gene", "logFC", "t", "P.Value", "adj.P.Val", "B")]

  # Save with suffix indicating CSN_SCORE adjusted
  out_dir <- file.path(out_root, "DEG", pass_label, paste0(subunit_name, "_adj_CSN"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(tbl, file.path(out_dir, "DEG_limma_cont.csv"))

  log_msg(
    "[DEG-INTERACTION] %s | %s | %s (controlling for CSN_SCORE): Wrote %s",
    ds_id, basename(out_root), subunit_name,
    file.path(out_dir, "DEG_limma_cont.csv")
  )
  invisible(NULL)
}


## === robust run_predictor_analyses (center predictors for limma DEG models) ===
run_predictor_analyses <- function(
  predictor_name,
  predictor_vec, # named numeric by sample (can contain NA)
  exclude_genes = NULL, # Genes to exclude before ranking (e.g., self/whole complex)
  ds_id, ds_dir,
  mat0, # Genes x samples; no imputation (This function handles limma continuous DEG)
  mat, # After impute+filter (For limma)
  out_root,
  genesets_by_group, # list(group -> pathways)
  batch_all = NULL, # factor by sample or NULL
  purity_all = NULL, # named numeric by sample or NULL
  sa_all = NULL, # data.frame(sex, age[, age_missing, age_z_imputed]); rownames=sample
  tp53_num_all = NULL, # named numeric (0/1) or NULL
  is_ALL = FALSE # Only possible to add TP53_mutant when TRUE
) {
  ## -------- Safety parameters --------
  opt <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
  min_per_group <- opt("min_per_group", 8L)
  minSize <- opt("minSize", 15L)
  maxSize <- opt("maxSize", 500L)
  fgsea_eps <- opt("fgsea_eps", 0)
  USE_AGE_MISSING_INDICATOR <- isTRUE(opt("USE_AGE_MISSING_INDICATOR", FALSE))

  logf <- function(...) if (exists("log_msg", mode = "function")) try(log_msg(...), silent = TRUE)
  sfn <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

  stopifnot(!is.null(genesets_by_group), length(genesets_by_group) > 0)

  ## -------- Align samples (Based on predictor) --------
  if (is.null(names(predictor_vec))) names(predictor_vec) <- colnames(mat)
  keep <- intersect(colnames(mat), names(predictor_vec))
  pred_all <- suppressWarnings(as.numeric(predictor_vec[keep]))
  names(pred_all) <- keep
  fin <- is.finite(pred_all)
  if (sum(fin) < (2L * min_per_group)) {
    logf("  [%s] predictor non-NA samples insufficient (%d < %d) -> skip", predictor_name, sum(fin), 2L * min_per_group)
    return(invisible(NULL))
  }
  sample_order <- keep[fin]
  pred <- pred_all[fin] # Original (uncentered) predictor
  stopifnot(identical(names(pred), sample_order))

  ## -------- Build covariates data frame (No forced NA imputation; align sample order only) --------
  build_covars_df <- function(so) {
    df <- data.frame(row.names = so, check.names = FALSE)
    if (!is.null(purity_all)) df$purity <- suppressWarnings(as.numeric(purity_all[so]))
    if (!is.null(sa_all)) {
      if ("sex" %in% colnames(sa_all)) df$sex <- factor(sa_all[so, "sex"])
      if ("age" %in% colnames(sa_all)) df$age <- suppressWarnings(as.numeric(sa_all[so, "age"]))
      if (USE_AGE_MISSING_INDICATOR) {
        if ("age_missing" %in% colnames(sa_all)) df$age_missing <- suppressWarnings(as.numeric(sa_all[so, "age_missing"]))
        if ("age_z_imputed" %in% colnames(sa_all)) df$age_z_imputed <- suppressWarnings(as.numeric(sa_all[so, "age_z_imputed"]))
      }
    }
    # [POLICY] TP53 covariate disabled in ALL strata
    df
  }
  df_covars0 <- build_covars_df(sample_order)

  coerce_covariates_safely <- function(df) {
    df <- as.data.frame(df, check.names = FALSE)
    keep <- rep(TRUE, ncol(df))
    names(keep) <- colnames(df)

    for (cn in colnames(df)) {
      v <- df[[cn]]

      if (is.factor(v) || is.character(v) || is.logical(v)) {
        v <- factor(v) # Keep as factor
        # Check actual levels using only non-NA samples
        lv <- levels(droplevels(v[!is.na(v)]))
        if (length(lv) <= 1) {
          # Single-level factor -> Drop, avoid contrasts error
          keep[cn] <- FALSE
          if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
        } else {
          df[[cn]] <- v
        }
      } else {
        # Keep numeric as numeric; do not force convert to factor
        df[[cn]] <- suppressWarnings(as.numeric(v))
        # (Optional) If you want to remove "completely constant" numeric columns together, uncomment below:
        # vv <- df[[cn]]; if (sum(!is.na(vv)) >= 2 && isTRUE(all(vv[!is.na(vv)] == vv[which(!is.na(vv))[1]]))) {
        #   keep[cn] <- FALSE
        #   if (exists("logf")) try(logf("  [covars] drop constant numeric: %s", cn), silent = TRUE)
        # }
      }
    }

    df <- df[, keep, drop = FALSE]
    df
  }


  ## -------- coverage threshold --------
  base_thr <- c(purity = 0.60, sex = 0.80, age = 0.80)
  present <- colnames(df_covars0)
  extra <- setdiff(present, names(base_thr))
  if (length(extra)) base_thr <- c(base_thr, stats::setNames(rep(min(base_thr), length(extra)), extra))

  ## -------- Select covariates (Return numeric data.frame aligned with sample_order) --------
  pick_covars_df <- function(label) {
    ## NEW: Determine covariate source (fallback to df_covars0 if covars_all not in scope)
    cov_src <- tryCatch(get("covars_all", inherits = TRUE), error = function(e) NULL)
    if (is.null(cov_src)) {
      cov_src <- tryCatch(get("df_covars0", inherits = TRUE), error = function(e) NULL)
    }
    if (is.null(cov_src)) {
      logf("  [covars-%s] NO covariate table in scope -> skip", label)
      return(NULL)
    }
    sel <- tryCatch(
      select_covars_safely(
        df           = cov_src,
        sample_order = sample_order,
        label        = label,
        y            = pred
      ),
      error = function(e) {
        logf("  [covars-%s] SELECT-FAIL: %s", label, conditionMessage(e))
        ## Diagnostic info: Use cov_src again, avoid error due to non-existent object
        rns <- tryCatch(rownames(cov_src), error = function(e) NULL)
        nrowC <- tryCatch(NROW(cov_src), error = function(e) -1L)
        ncolC <- tryCatch(NCOL(cov_src), error = function(e) -1L)
        hasrn <- is.character(rns) || is.factor(rns)
        logf("  [covars-%s] diag: C_dim = %d x %d; has_rownames=%s", label, nrowC, ncolC, hasrn)
        miss <- tryCatch(setdiff(sample_order, rns), error = function(e) character(0))
        logf(
          "  [covars-%s] samples_not_in_covars = %d / %d (eg: %s)",
          label, length(miss), length(sample_order), paste(head(miss, 5), collapse = ",")
        )
        NULL
      }
    )
    if (is.null(sel)) {
      return(NULL)
    }
    # === NEW: Extract drop reasons returned by select_covars_safely() (Save variable first) ===
    dl_sel <- tryCatch(attr(sel, "drop_log"), error = function(e) NULL)
    ## === /NEW ===
    X <- NULL
    if (is.character(sel) && is.null(dim(sel))) {
      nm <- intersect(sel, colnames(cov_src))
      if (length(nm)) X <- cov_src[, nm, drop = FALSE]
    } else {
      X <- as.data.frame(sel, stringsAsFactors = FALSE)
      if (is.null(rownames(X))) {
        if (nrow(X) == length(sample_order)) {
          rownames(X) <- sample_order
        } else if (!is.null(colnames(X)) && ncol(X) == length(sample_order) && all(colnames(X) == sample_order)) {
          X <- as.data.frame(t(as.matrix(X)))
        } else {
          logf("  [covars-%s] Return shape/alignment unknown, discard", label)
          return(NULL)
        }
      }
      X <- X[sample_order, , drop = FALSE]
    }
    if (is.null(X) || !ncol(X)) {
      return(NULL)
    }
    X <- coerce_covariates_safely(X) # Let character/logical -> factor; numeric remains numeric
    # Factor-friendly minimum quality threshold: At least 3 non-NA
    good <- vapply(X, function(v) sum(!is.na(v)) >= 3, logical(1))
    if (any(!good)) {
      logf("  [covars-%s] drop low-coverage columns: %s", label, paste(colnames(X)[!good], collapse = ","))
      X <- X[, good, drop = FALSE]
    }
    if (!ncol(X)) {
      return(NULL)
    }
    logf("  [covars-picked:%s] kept = {%s}", label, if (!is.null(X)) paste(colnames(X), collapse = ",") else "NULL")
    ## === NEW: Save dropped covariates and reasons to CSV (Append mode) ===
    if (!is.null(dl_sel) && nrow(dl_sel)) {
      dl_sel$dataset <- ds_id
      dl_sel$stratum <- basename(out_root)
      dl_sel$subunit <- predictor_name
      dl_sel$pass <- label # "limma-cont:RAW" or "limma-cont:base"
      dl_sel$samples <- length(sample_order)

      out_dir <- file.path("run_info", "covars_audit")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      fp <- file.path(out_dir, "covariate_drop_reasons.csv")

      data.table::fwrite(dl_sel, fp, append = file.exists(fp))
    }
    ## === /NEW ===
    X
  }
  X_raw_cov <- pick_covars_df("limma-cont:RAW")
  if (!is.null(X_raw_cov)) X_raw_cov <- coerce_covariates_safely(X_raw_cov)
  X_ba_cov <- pick_covars_df("limma-cont:base")
  if (!is.null(X_ba_cov)) X_ba_cov <- coerce_covariates_safely(X_ba_cov)
  ## -------- RAW: Build DF -> Drop NA -> **Center predictor** -> model.matrix --------
  DF_raw <- data.frame(predictor = as.numeric(pred), row.names = sample_order, check.names = FALSE)
  if (!is.null(X_raw_cov)) for (nm in colnames(X_raw_cov)) DF_raw[[nm]] <- X_raw_cov[[nm]]
  ok_raw <- stats::complete.cases(DF_raw)
  if (sum(ok_raw) < (2L * min_per_group)) {
    logf("  [%s|RAW] Insufficient available samples (%d < %d) -> skip", predictor_name, sum(ok_raw), 2L * min_per_group)
    return(invisible(NULL))
  }
  DF_raw <- DF_raw[ok_raw, , drop = FALSE]
  ## -- Centering (Shift mean only, no scaling): Does not affect coefficients/tests, only intercept interpretation differs
  DF_raw$predictor <- DF_raw$predictor - mean(DF_raw$predictor, na.rm = TRUE)
  M_raw <- as.matrix(mat[, rownames(DF_raw), drop = FALSE])
  des_raw <- stats::model.matrix(~ 1 + ., data = DF_raw, na.action = stats::na.fail)
  stopifnot(nrow(des_raw) == ncol(M_raw))
  # Audit: Coverage/Covariates summary + Alignment
  audit_covars_coverage(
    tag = sprintf("%s_%s_%s|RAW", ds_id, basename(out_root), predictor_name),
    ds_id = ds_id,
    stratum = basename(out_root),
    su = predictor_name,
    sample_ids = rownames(DF_raw),
    batch = if (!is.null(batch_all)) droplevels(as.factor(batch_all[rownames(DF_raw)])) else NULL,
    covars = DF_raw[, setdiff(colnames(DF_raw), "predictor"), drop = FALSE],
    sv = NULL,
    tech = NULL
  ) # -> Will append to run_info/covars_audit/audit_rows.csv  :contentReference[oaicite:17]{index=17}

  audit_design_alignment(
    tag = sprintf("%s_%s_%s|RAW", ds_id, basename(out_root), predictor_name),
    samples = colnames(M_raw),
    mod_interest = des_raw,
    mod_nuisance = NULL,
    out_dir = file.path("run_info", "covars_audit")
  ) # -> Will write ALIGN_RAW.csv type filename  :contentReference[oaicite:18]{index=18}

  ## -------- BatchAdj: Same as above; batch into design table -> **Center predictor** --------
  DF_ba <- data.frame(predictor = as.numeric(pred), row.names = sample_order, check.names = FALSE)
  if (!is.null(X_ba_cov)) for (nm in colnames(X_ba_cov)) DF_ba[[nm]] <- X_ba_cov[[nm]]
  if (!is.null(batch_all)) DF_ba$batch <- droplevels(as.factor(batch_all[sample_order]))

  ## NEW: Remove all-NA columns (e.g., purity version dropped by coverage but remaining)
  all_na_col <- vapply(DF_ba, function(v) all(is.na(v)), logical(1))
  if (any(all_na_col)) {
    logf(
      "[BatchAdj] drop all-NA columns: %s",
      paste(names(all_na_col)[all_na_col], collapse = ",")
    )
    DF_ba <- DF_ba[, !all_na_col, drop = FALSE]
  }

  ## Select columns actually used in this round (predictor + selected covars)
  use_cols <- unique(c(
    "predictor",
    intersect(c("purity", "sex", "age", "batch", "TP53_mutant"), colnames(DF_ba))
  ))
  use_cols <- use_cols[use_cols %in% colnames(DF_ba)] # Guard

  ## Check complete cases using "columns to be used"
  ok_rows <- complete.cases(DF_ba[, use_cols, drop = FALSE])
  n_ok <- sum(ok_rows)

  if (n_ok < 16L) {
    logf("  [%s|BatchAdj] Insufficient available samples (%d < 16) -> skip", predictor_name, n_ok)
    return(invisible(NULL)) # <- Consistent with your original control flow
  }

  ## Align to samples with complete cases (This step is critical)
  DF_ba <- DF_ba[ok_rows, , drop = FALSE]

  ## -- Centering (Keep your original design)
  DF_ba$predictor <- DF_ba$predictor - mean(DF_ba$predictor, na.rm = TRUE)

  ## Generate M_ba using final sample subset (Keep your original writing)
  M_ba <- as.matrix(mat[, rownames(DF_ba), drop = FALSE])

  form_ba <- if (ncol(DF_ba) > 0) {
    as.formula(paste("~ 1 +", paste(colnames(DF_ba), collapse = " + ")))
  } else {
    ~1
  }
  des_ba <- stats::model.matrix(form_ba, data = DF_ba, na.action = stats::na.fail)
  stopifnot(nrow(des_ba) == ncol(M_ba))
  audit_covars_coverage(
    tag        = sprintf("%s|BatchAdj", predictor_name),
    ds_id      = ds_id,
    stratum    = basename(out_root),
    su         = predictor_name,
    sample_ids = rownames(DF_ba),
    batch      = if (!is.null(batch_all)) droplevels(as.factor(batch_all[rownames(DF_ba)])) else NULL,
    covars     = DF_ba[, intersect(c("purity", "sex", "age", "batch", "TP53_mutant"), colnames(DF_ba)), drop = FALSE]
  )

  ## --- NEW: Add dataset/stratum prefix only to RESIDUAL_* (BatchAdj); keep others as is ---
  tag_ba_align <- if (grepl("^RESIDUAL_", predictor_name)) {
    sprintf("%s_%s_%s|BatchAdj", ds_id, basename(out_root), predictor_name)
  } else {
    sprintf("%s|BatchAdj", predictor_name)
  }
  audit_design_alignment(
    tag = tag_ba_align,
    samples = colnames(M_ba),
    mod_interest = des_ba,
    mod_nuisance = NULL,
    out_dir = file.path("run_info", "covars_audit")
  )


  ## -------- Helper tools --------
  coef_name <- function(des) if ("predictor" %in% colnames(des)) "predictor" else colnames(des)[ncol(des)]
  rank_finite <- function(v, exclude = NULL, tag = NULL) {
    keep <- is.finite(v)
    if (any(!keep)) logf("[gsea-%s] drop non-finite: %d", ifelse(is.null(tag), "NA", tag), sum(!keep))
    v <- v[keep]
    if (!is.null(exclude) && length(exclude)) v <- v[setdiff(names(v), exclude)]
    if (!length(v)) {
      return(NULL)
    }
    v
  }
  gsea_from <- function(pw, stats) {
    if (is.null(pw) || !length(pw) || is.null(stats) || !length(stats)) {
      return(NULL)
    }
    # Read global thresholds
    minSize <- opt("minSize", 15L)
    maxSize <- opt("maxSize", 500L)
    # Explicit intersection + size filter
    orig_n <- length(pw)
    pw_use <- ._intersect_and_filter_pathways(pw, names(stats), minSize = minSize, maxSize = maxSize)
    message(sprintf("[gsea] (gsea_from) |H| Orig=%d -> Used=%d", orig_n, length(pw_use)))
    if (!length(pw_use)) {
      return(NULL)
    }

    set.seed(1L)
    res <- tryCatch(
      {
        suppressWarnings(fgsea::fgseaMultilevel(
          pathways = pw_use, stats = stats,
          minSize = minSize, maxSize = maxSize,
          eps = 1e-10
        ))
      },
      error = function(e) {
        message(sprintf("[gsea] Multilevel failed: %s -> Switch to fgseaSimple", conditionMessage(e)))
        suppressWarnings(fgsea::fgseaSimple(
          pathways = pw_use, stats = stats,
          nperm = 10000, minSize = minSize, maxSize = maxSize
        ))
      }
    )
    res
  }
}


## ===== CSN subunits coverage & CSN_SCORE (PC1) feasibility audit (Robust version) =====
audit_csn_score_feasibility <- function(ds_id, stratum, mat0, present_sub,
                                        min_members = 5L,
                                        pca_min_samples = 10L, # Consistent with threshold in build_csn_score
                                        min_per_group = 8L, # Your limma threshold
                                        out_dir = file.path("run_info", "csn_score_audit")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  n_samples <- ncol(mat0)
  if (!is.matrix(mat0) || is.null(rownames(mat0)) || is.null(colnames(mat0))) {
    log_msg("[CSN-audit] %s | %s: mat0 structure incomplete, skip", ds_id, stratum)
    return(invisible(NULL))
  }
  present_sub <- intersect(present_sub, rownames(mat0))
  if (!length(present_sub) || n_samples == 0) {
    log_msg("[CSN-audit] %s | %s: No available CSN subunit or samples is 0, skip", ds_id, stratum)
    return(invisible(NULL))
  }

  ## 1) Coverage of each subunit
  sub_cov <- vapply(present_sub, function(g) mean(is.finite(mat0[g, ])) * 100, numeric(1))
  sub_tbl <- data.frame(
    dataset = ds_id,
    stratum = stratum,
    subunit = present_sub,
    nonNA_pct = round(sub_cov, 1),
    nonNA_n = vapply(present_sub, function(g) sum(is.finite(mat0[g, ])), integer(1)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  cov_min <- if (length(sub_cov)) round(min(sub_cov, na.rm = TRUE), 1) else NA_real_
  cov_med <- if (length(sub_cov)) round(stats::median(sub_cov, na.rm = TRUE), 1) else NA_real_
  cov_max <- if (length(sub_cov)) round(max(sub_cov, na.rm = TRUE), 1) else NA_real_

  ## 2) How many subunits non-NA per sample; sample count satisfying >= min_members
  sample_counts <- colSums(is.finite(mat0[present_sub, , drop = FALSE]))
  enough <- sample_counts >= min_members
  n_enough <- sum(enough)

  ## 3) Try calculating CSN_SCORE (PC1), and record feasibility
  csn_score <- build_csn_score(mat0,
    subunits = present_sub,
    combine_7AB = TRUE, min_members = min_members
  )
  csn_nonNA <- sum(is.finite(csn_score))
  csn_can_pca <- (n_enough >= pca_min_samples) && (csn_nonNA >= pca_min_samples)

  ## Extra: PC1 explained variance (when feasible)
  pc1_var_pct <- NA_real_
  if (isTRUE(csn_can_pca)) {
    ok_sam <- names(sample_counts)[sample_counts >= min_members]
    get_z <- function(v) {
      v <- as.numeric(v)
      mu <- mean(v[is.finite(v)], na.rm = TRUE)
      sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      v[!is.finite(v)] <- mu
      (v - mu) / sdv
    }
    Z <- do.call(rbind, lapply(present_sub, function(g) get_z(mat0[g, ok_sam, drop = FALSE])))
    rownames(Z) <- present_sub
    if (all(c("COPS7A", "COPS7B") %in% rownames(Z))) {
      Z7 <- colMeans(Z[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
      Z <- rbind(Z[setdiff(rownames(Z), c("COPS7A", "COPS7B")), , drop = FALSE],
        "COPS7*" = Z7
      )
    }
    pc <- tryCatch(stats::prcomp(t(Z), center = TRUE, scale. = FALSE), error = function(e) NULL)
    if (!is.null(pc)) {
      ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
      pc1_var_pct <- round(ve[1], 1)
    }
  }

  ## 4) Write files
  # 4a) Summary table (append)
  sum_row <- data.frame(
    dataset = ds_id, stratum = stratum,
    n_samples = n_samples,
    n_subunits_present = length(present_sub),
    subunit_cov_min = cov_min, subunit_cov_median = cov_med, subunit_cov_max = cov_max,
    min_members = min_members,
    samples_with_ge_min_members = n_enough,
    pca_min_samples = pca_min_samples,
    csn_score_nonNA = csn_nonNA,
    csn_pc1_feasible = csn_can_pca,
    csn_pc1_var_pct = pc1_var_pct,
    can_form_high_low_groups = (csn_nonNA >= 2 * min_per_group),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  fp_sum <- file.path(out_dir, "csn_score_feasibility_summary.csv")
  data.table::fwrite(sum_row, fp_sum, append = file.exists(fp_sum))

  # 4b) subunit coverage (one file per stratum)
  tag <- paste(ds_id, stratum, sep = "_")
  fp_sub <- file.path(out_dir, sprintf("%s_subunit_coverage.csv", tag))
  data.table::fwrite(sub_tbl, fp_sub)

  # 4c) non-NA subunit count per sample (one file per stratum)
  sample_tbl <- data.frame(
    dataset = ds_id, stratum = stratum,
    sample_id = colnames(mat0),
    nonNA_subunits_n = as.integer(sample_counts),
    ge_min_members = as.logical(enough),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  fp_sam <- file.path(out_dir, sprintf("%s_sample_subunit_counts.csv", tag))
  data.table::fwrite(sample_tbl, fp_sam)

  log_msg(
    "[CSN-audit] %s | %s: subunits=%d; min/median/max coverage=%.1f/%.1f/%.1f%%; eligible_samples=%d; CSN_SCORE nonNA=%d; PC1 feasible=%s (PC1%%=%.1f)",
    ds_id, stratum, length(present_sub), cov_min %||% NaN, cov_med %||% NaN, cov_max %||% NaN,
    n_enough, csn_nonNA, as.character(csn_can_pca), pc1_var_pct %||% NaN
  )

  invisible(list(summary = sum_row, per_subunit = sub_tbl, per_sample = sample_tbl))
}


## ===== Small function to run complete GSEA for "specified sample set" =====
## =========================================================
## Each stratum: Produce RAW and BatchAdj simultaneously (New version)
## =========================================================
## ===== run_one_stratum: Unify Spearman filename, write sample_counts, enhance batch notes (2,3,7) =====
## ===== Small function to run complete GSEA for "specified sample set" (Fixed version) =====
## ===== Small function to run complete GSEA for "specified sample set" (Fixed version; BatchAdj can also take age_missing/age_z_imputed) =====
run_one_stratum <- function(ds_id, ds_dir, mat0_full, sample_keep, out_root, genesets_by_group) {
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  log_msg("  -- stratum：%s | N(sample_keep)=%d", basename(out_root), length(sample_keep))
  ## [DEG mode] Set flags for DEG analysis only
  .RUN_DEG <<- TRUE


  log_msg("[FLAGS] .RUN_DEG=%s", as.character(.RUN_DEG))


  ## 1) Subset samples + Size check + Log scale
  keep <- intersect(colnames(mat0_full), sample_keep)
  if (length(keep) < 4) {
    log_msg("  [Skip] Samples too few: %d", length(keep))
    return(invisible(NULL))
  }
  mat0 <- mat0_full[, keep, drop = FALSE]
  mx <- suppressWarnings(max(mat0, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)

  ## 2) Imputation+Filter matrix for limma (Spearman still uses raw mat0) + CSN_SCORE feasibility audit
  mat <- impute_and_filter(mat0, min_frac = min_frac_complete)
  present_sub <- intersect(csn_subunits, rownames(mat0))
  if (!length(present_sub)) {
    log_msg("  [Skip] No CSN subunit in this stratum")
    return(invisible(NULL))
  }
  audit_csn_score_feasibility_safe(
    ds_id = ds_id,
    stratum = basename(out_root),
    mat0 = mat0,
    present_sub = present_sub,
    min_members = 5L,
    pca_min_samples = 10L,
    min_per_group = min_per_group,
    out_dir = file.path("run_info", "csn_score_audit")
  )

  ## === In ALL analysis: Prepare TP53 status (0/1; mutant=1) ===
  is_ALL <- identical(basename(out_root), "ALL")
  tp53_num_all <- NULL
  if (is_ALL) {
    tp53_status_all <- get_tp53_status(ds_dir, colnames(mat0))
    tp53_num_all <- as.numeric(tp53_status_all == "TP53_mutant")
    names(tp53_num_all) <- colnames(mat0)
  }

  ## 3) Covariates / Batch (Align samples) -- Build covariates for limma (Optionally take missing-indicator)
  bi_all <- get_batch_factor(ds_dir, colnames(mat0))
  batch_all <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL
  purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
  sa_all <- get_sex_age_covariates(ds_dir, colnames(mat0)) # data.frame(sex, age [, age_missing, age_z_imputed])
  sa_all_limma <- sa_all
  ## ---- [NEW2] DEG_EAGER: Run early (Avoid DEG being skipped by any early-return after step 4) ----
  if (isTRUE(get0(".RUN_DEG", ifnotfound = FALSE)) ||
    isTRUE(get0("RUN_DEG", ifnotfound = FALSE))) {
    ## 4a-1) Each subunit (predictor = protein ratio of that subunit)
    for (su in present_sub) {
      run_deg_for_predictor(
        predictor_name = su,
        predictor_vec = mat0[su, ],
        ds_id = ds_id, ds_dir = ds_dir,
        mat0 = mat0,
        out_root = out_root,
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma
      )
    }

    ## 4b-1) CSN complex score (PC1; direction corrected); exclude CSN members when ranking
    csn_score <- build_csn_score_safe(
      mat0,
      subunits = present_sub, combine_7AB = TRUE,
      min_members = 5L, pca_min_samples = 10L
    )
    if (sum(is.finite(csn_score)) >= (2 * min_per_group)) {
      run_deg_for_predictor(
        predictor_name = "CSN_SCORE",
        predictor_vec = csn_score,
        ds_id = ds_id, ds_dir = ds_dir,
        mat0 = mat0,
        out_root = out_root,
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma
      )
    } else {
      log_msg("  [DEG-CSN_SCORE] Insufficient non-NA samples, skip")
    }

    ## 4c-1) residual subunits (Require base_covars_all)
    base_covars_all <- data.frame(
      purity = as.numeric(purity_all[colnames(mat0)]),
      sex = suppressWarnings(as.numeric(sa_all[colnames(mat0), "sex"])),
      age = suppressWarnings(as.numeric(sa_all[colnames(mat0), "age"])),
      row.names = colnames(mat0), check.names = FALSE
    )
    if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing", "age_z_imputed") %in% colnames(sa_all))) {
      base_covars_all$age_missing <- as.numeric(sa_all[colnames(mat0), "age_missing"])
      base_covars_all$age_z_imputed <- as.numeric(sa_all[colnames(mat0), "age_z_imputed"])
    }
    is_ALL <- identical(basename(out_root), "ALL")
    tp53_num_all <- if (is_ALL) {
      status_all <- get_tp53_status(ds_dir, colnames(mat0))
      setNames(as.numeric(status_all == "TP53_mutant"), colnames(mat0))
    } else {
      NULL
    }
    if (is_ALL && !is.null(tp53_num_all)) {
      base_covars_all$TP53_mutant <- as.numeric(tp53_num_all[colnames(mat0)])
    }

    min_n_resid <- opt("min_n_resid", 10L)
    csn_nonNA <- sum(is.finite(csn_score))
    if (csn_nonNA < min_n_resid) {
      log_msg("  [DEG-INTERACTION] CSN score non-NA samples insufficient (%d < %d), skip all subunit interaction models", csn_nonNA, min_n_resid)
    } else {
      ## NEW APPROACH: One-step interaction model
      ## Instead of two-step residualization, use single limma model with both
      ## subunit and CSN_SCORE. This is statistically correct and tests the same
      ## biological hypothesis: "subunit effect independent of CSN complex activity"

      log_msg("  [DEG-INTERACTION] Running one-step models (subunit + CSN_SCORE + covariates)")

      for (su in present_sub) {
        # Check if subunit has sufficient non-NA samples
        su_nonNA <- sum(is.finite(mat0[su, ]))
        if (su_nonNA < (2 * min_per_group)) {
          log_msg("  [%s_adj_CSN] Insufficient non-NA samples, skip", su)
          next
        }

        # Run ONE-STEP interaction model
        run_deg_interaction_model(
          subunit_name = su,
          subunit_vec = mat0[su, ], # Subunit expression
          csn_score_vec = csn_score, # CSN complex score
          ds_id = ds_id, ds_dir = ds_dir,
          mat0 = mat0,
          out_root = out_root,
          batch_all = batch_all,
          purity_all = purity_all,
          sa_all = sa_all_limma,
          min_per_group = min_per_group
        )
      }
    }
  }
  if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing", "age_z_imputed") %in% colnames(sa_all))) {
    keep_cols <- intersect(c("sex", "age", "age_missing", "age_z_imputed"), colnames(sa_all))
    sa_all_limma <- sa_all[, keep_cols, drop = FALSE]
  }


  ## 5) Aggregate each version independently
  present_sub <- intersect(csn_subunits, rownames(mat0))
  # Updated: Use new interaction model naming (*_adj_CSN instead of RESIDUAL_*)
  sum_units <- c(present_sub, "CSN_SCORE", paste0(present_sub, "_adj_CSN"))

  # RAW: limma_t_cont + interaction (Aggregate existing CSVs only, do not rerun)
  passes <- getOption("csn.run_passes", default = c("BatchAdj"))
  if (("RAW" %in% passes) && (isTRUE(get0(".RUN_LIMMA", ifnotfound = FALSE)) || isTRUE(get0(".RUN_SPEARMAN", ifnotfound = FALSE)))) {
    for (grp_name in names(genesets_by_group)) {
      out_root_coll <- file.path(OUTPUT_PREFIX, "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), ds_id, basename(out_root))
      ver_root <- file.path(OUTPUT_PREFIX, "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), "RAW", ds_id, basename(out_root))
      summarize_all_groups(
        out_root = ver_root,
        csn_subunits = sum_units,
        genesets_by_group = setNames(list(genesets_by_group[[grp_name]]), grp_name),
        stat_tags = c("GSEA_limma_t_cont", "GSEA_limma_interaction")
      )
    }
  }
  # BatchAdj: limma_t_cont + spearman + interaction (Aggregate existing CSVs only, do not rerun)
  if (("BatchAdj" %in% passes) && (isTRUE(get0(".RUN_LIMMA", ifnotfound = FALSE)) || isTRUE(get0(".RUN_SPEARMAN", ifnotfound = FALSE)))) {
    for (grp_name in names(genesets_by_group)) {
      out_root_coll <- file.path(OUTPUT_PREFIX, "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), ds_id, basename(out_root))
      ver_root <- file.path(OUTPUT_PREFIX, "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root))
      summarize_all_groups(
        out_root = ver_root,
        csn_subunits = sum_units,
        genesets_by_group = setNames(list(genesets_by_group[[grp_name]]), grp_name),
        stat_tags = {
          st <- c("GSEA_limma_t_cont", "GSEA_spearman", "GSEA_limma_interaction")
          if (isTRUE(RUN_HILO_SENSITIVITY)) st <- c(st, "GSEA_limma_t_hilo")
          st
        }
      )
    }
  }
  invisible(NULL)
}


## Inspect which PLEX are currently fetched (Original vs Cleaned)
inspect_plex <- function(ds_dir, col = "TMT_PLEX",
                         pipe_policy = BATCH_PIPE_POLICY,
                         min_per_level = BATCH_MIN_PER_LEVEL) {
  # Read samples actually used (Based on protein matrix)
  mat0 <- load_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(mat0)

  # Read clinical table and align samples
  meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  stopifnot(file.exists(meta_fp))
  meta <- suppressMessages(readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
  id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"), names(meta))
  stopifnot(length(id_cols) > 0)
  meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
  # Get original and cleaned PLEX
  raw <- as.character(meta[[col]])
  clean <- sanitize_batch_levels(raw, pipe_policy = pipe_policy, min_per_level = min_per_level)

  # Summary
  cat("\n==== ", basename(ds_dir), " | Column: ", col, " ====\n", sep = "")
  cat("[Original PLEX levels]:\n")
  print(sort(table(raw), decreasing = TRUE))
  cat("\n[Original values containing '|' (first few)]:\n")
  print(utils::head(unique(raw[grepl("\\|", raw %||% "")]), 10))

  cat("\n[Cleaned PLEX levels] (pipe_policy = ", pipe_policy,
    ", min_per_level = ", min_per_level, "):\n",
    sep = ""
  )
  print(sort(table(clean, useNA = "ifany"), decreasing = TRUE))

  # Output mapping table (First 10 rows demo)
  df_map <- data.frame(SAMPLE_ID = sample_ids, raw_plex = raw, clean_plex = as.character(clean))
  cat("\n[Sample mapping table (First 10 rows)]\n")
  print(utils::head(df_map, 10))

  invisible(list(
    raw_counts = sort(table(raw), decreasing = TRUE),
    clean_counts = sort(table(clean, useNA = "ifany"), decreasing = TRUE),
    map = df_map
  ))
}


## ===== Missingness sensitivity（AGE）=====
USE_AGE_MISSING_INDICATOR <- FALSE # Main analysis default: Do not use missing-indicator, keep NA only

## Safe z-score: Estimate mean/sd using finite values only, keep NA, no mean imputation
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

## Tools: Column name normalization, z-score, 0~1
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

## For PAAD: Take median of values split by semicolon
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
## Get sex/age (Sample aligned; sex: 0/1; age: z-score)
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

## For main pipeline: Strictly generate covariates using SEX and AGE from patient file (sex: 0/1; age: z-score)
get_sex_age_covariates <- function(ds_dir, sample_ids) {
  samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")

  if (!file.exists(samp_fp) || !file.exists(pat_fp)) {
    log_msg("[covars] Missing sample/patient file, replace sex/age with NA")
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
  # If .zscore (z after mean imputation) is not defined externally, provide fallback to maintain consistent logic
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
    stop("[get_sex_age_covariates] Cannot establish sample<->patient mapping (Missing SAMPLE_ID/PATIENT_ID)")
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

  ## AGE: No imputation for main analysis; Sensitivity can choose missing-indicator + z after mean imputation
  age_col <- intersect(c("AGE", "AGE_AT_DIAGNOSIS", "AGE_AT_INDEX", "AGE_YEARS"), names(pat))[1]
  age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  age_missing <- age_z_imputed <- NULL
  if (!is.na(age_col)) {
    v <- suppressWarnings(as.numeric(pat[[age_col]]))
    names(v) <- as.character(pat[[pid_pat]])
    age_raw <- unname(v[map_pt[sample_ids]])
    age[] <- .z_no_impute(age_raw) # Keep NA (Main analysis)
    if (!exists("USE_AGE_MISSING_INDICATOR", inherits = FALSE)) USE_AGE_MISSING_INDICATOR <- FALSE
    if (isTRUE(USE_AGE_MISSING_INDICATOR)) {
      age_missing <- as.numeric(is.na(age_raw))
      age_z_imputed <- .zscore(age_raw) # z after mean imputation (Sensitivity)
    }
  }

  cov_sex <- mean(is.finite(sex)) * 100
  cov_age <- mean(is.finite(age)) * 100
  if (isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))) {
    cov_age_imp <- mean(is.finite(age_z_imputed)) * 100
    cov_age_mis <- mean(is.finite(age_missing)) * 100
    log_msg(
      "    covariates coverage: sex %.1f%%, age(NA-as-NA) %.1f%%, age_z_imputed %.1f%%, age_missing %.1f%%",
      cov_sex, cov_age, cov_age_imp, cov_age_mis
    )
  } else {
    log_msg("    covariates coverage: sex %.1f%%, age(NA-as-NA) %.1f%%", cov_sex, cov_age)
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


select_covars_safely <- function(
  df, # Candidate covariates; rownames=sample_order
  sample_order,
  label = "covars",
  y = NULL, # Continuous predictor; can be NULL
  min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
  max_abs_cor = 0.30,
  min_pairs = 20L
) {
  logf <- function(...) if (exists("log_msg", mode = "function")) try(log_msg(...), silent = TRUE)
  ## NEW: Collector for drop reasons per column
  drop_log <- list()
  .append_drop <- function(step, col, reason, metric = NA_real_, threshold = NA_real_, extra = NA_character_) {
    drop_log[[length(drop_log) + 1L]] <<- data.frame(
      step = step,
      covariate = col,
      reason = reason,
      metric = metric,
      threshold = threshold,
      extra = extra,
      stringsAsFactors = FALSE
    )
  }
  if (is.null(df) || !nrow(df)) {
    return(NULL)
  }


  # 1) Align sample order
  if (is.null(rownames(df))) {
    logf("  [covars-%s] df has no rownames -> skip", label)
    return(NULL)
  }
  so <- as.character(sample_order)
  logf("  [covars-%s] before align: C_dim=%d x %d; has_rownames=%s", label, NROW(df), NCOL(df), !is.null(rownames(df)))
  df <- df[so, , drop = FALSE]

  ## [POLICY] Never include TP53 as a covariate (ALL / MT / WT strata)
  tp_cols_idx <- grep("^TP53($|_|)", colnames(df), ignore.case = FALSE)
  if (length(tp_cols_idx) > 0L) {
    tp_cols <- colnames(df)[tp_cols_idx]
    for (cc in tp_cols) {
      .append_drop(
        step = "policy", col = cc, reason = "exclude_TP53_as_covariate",
        metric = NA_real_, threshold = NA_real_, extra = "global policy"
      )
    }
    df <- df[, setdiff(colnames(df), tp_cols), drop = FALSE]
    logf(sprintf("  [covars-%s] policy: drop TP53 columns from covariates -> %s", label, paste(tp_cols, collapse = ",")))
  }

  logf("  [covars-%s] after  align: C_dim=%d x %d", label, NROW(df), NCOL(df))


  # 2) Always use data.frame; keep factors
  df <- as.data.frame(df, stringsAsFactors = TRUE, check.names = FALSE)

  # 3) coverage: Use global minimum threshold for unnamed columns (Safe indexing, do not use [[cn]])
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  global_min <- if (length(min_cov_named)) min(min_cov_named, na.rm = TRUE) else 0
  keep_cov <- rep(TRUE, ncol(df))
  names(keep_cov) <- colnames(df)

  for (cn in colnames(df)) {
    v <- df[[cn]]
    # Unified coverage definition: is.finite for numeric; !is.na for non-numeric
    cover <- if (is.numeric(v)) mean(is.finite(v)) else mean(!is.na(v))

    thr <- if (!is.null(names(min_cov_named)) && (cn %in% names(min_cov_named))) {
      min_cov_named[[cn]] # * Safe: Only [[cn]] when name exists
    } else {
      global_min
    }

    if (is.na(cover) || cover < thr) {
      logf(
        "  [covars-%s] drop %s (coverage=%s < %s)", label, cn,
        sprintf("%.0f%%", 100 * ifelse(is.na(cover), 0, cover)),
        sprintf("%.0f%%", 100 * thr)
      )
      keep_cov[cn] <- FALSE
      if (is.na(cover) || cover < thr) {
        logf(
          "  [covars-%s] drop %s (coverage=%s < %s)", label, cn,
          sprintf("%.0f%%", 100 * ifelse(is.na(cover), 0, cover)),
          sprintf("%.0f%%", 100 * thr)
        )
        keep_cov[cn] <- FALSE
        ## NEW:
        .append_drop("coverage", cn, "low_coverage",
          metric = ifelse(is.na(cover), 0, cover), threshold = thr
        )
      }
    }
  }
  df <- df[, keep_cov, drop = FALSE]
  if (!ncol(df)) {
    return(NULL)
  }

  ## NEW: Re-initialize keep_cov, make it same length and column names as "shrunk df"
  keep_cov <- rep(TRUE, ncol(df))
  names(keep_cov) <- colnames(df)
  ## /NEW

  # 4) Filter by correlation with biology (Only for numeric columns)
  if (!is.null(y)) {
    y <- suppressWarnings(as.numeric(y))
    ## Do not do rho gate for these covariates
    skip_rho_gate <- c("purity", "age", "sex")

    for (cn in colnames(df)) {
      v <- df[[cn]]
      ## Only perform |rho| filter when "numeric and not in skip list"
      if (is.numeric(v) && !(cn %in% skip_rho_gate)) {
        fin <- is.finite(v) & is.finite(y)
        if (sum(fin) >= min_pairs) {
          r <- suppressWarnings(stats::cor(v[fin], y[fin], method = "spearman"))
          if (is.finite(r) && abs(r) >= max_abs_cor) {
            logf("  [covars-%s] drop %s (|rho|=%.3f ≥ %.2f vs biology)", label, cn, abs(r), max_abs_cor)
            keep_cov[cn] <- FALSE
          }
        }
      }
    }
    df <- df[, keep_cov, drop = FALSE]
    if (!ncol(df)) {
      return(NULL)
    }
  }


  # 5) Return: Row order = sample_order; for model.matrix expansion
  stopifnot(nrow(df) == length(so), identical(rownames(df), so))
  ## NEW: Attach drop reasons per column to attributes (Caller decides whether to write CSV)
  attr(df, "drop_log") <- if (length(drop_log)) do.call(rbind, drop_log) else NULL
  df
}


## ==============================
## Small Audit: sex/age + purity + batch overall check
## ==============================

# Safety: Helper functions (if not defined)
if (!exists(".norm_names", inherits = FALSE)) {
  .norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))
}
if (!exists("%||%", inherits = FALSE)) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
if (!exists(".to01", inherits = FALSE)) {
  .to01 <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    # If points > 1 are more than <= 1, treat as 0~100 percentage, convert to 0~1
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

# Helper: Format sample counts per batch level into string
.format_batch_sizes <- function(fac) {
  if (is.null(fac)) {
    return(NA_character_)
  }
  tb <- sort(table(fac), decreasing = TRUE)
  paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
}

# --- New: Audit purity source and coverage by dataset rules ---
.audit_purity_for_dataset <- function(ds_id, ds_dir, sample_ids) {
  samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
  has_samp <- file.exists(samp_fp)
  has_pat <- file.exists(pat_fp)

  purity_col <- "NONE"
  purity_vec <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)

  samp <- pat <- NULL
  sid <- pid_samp <- pid_pat <- NA

  if (has_samp) {
    samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    names(samp) <- .norm_names(names(samp))
    sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
    pid_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
  }
  if (has_pat) {
    pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    names(pat) <- .norm_names(names(pat))
    pid_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
  }

  # Build sample->patient mapping (when needed)
  map_pt <- NULL
  if (!is.na(sid) && !is.na(pid_samp)) map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])

  ds <- ds_id

  if (ds == "brca_cptac_2020" && has_pat && !is.na(pid_pat) && "ESTIMATE_TUMORPURITY" %in% names(pat) && !is.null(map_pt)) {
    purity_col <- "PATIENT:ESTIMATE_TUMORPURITY"
    pt_pur <- setNames(.to01(pat[["ESTIMATE_TUMORPURITY"]]), as.character(pat[[pid_pat]]))
    purity_vec[] <- unname(pt_pur[map_pt[sample_ids]])
  } else if (ds == "luad_cptac_2020" && has_samp && !is.na(sid) && "TUMOR_PURITY_BYESTIMATE_RNASEQ" %in% names(samp)) {
    purity_col <- "SAMPLE:TUMOR_PURITY_BYESTIMATE_RNASEQ"
    purity_vec[] <- .to01(samp[["TUMOR_PURITY_BYESTIMATE_RNASEQ"]][match(sample_ids, samp[[sid]])])
  } else if (ds == "lusc_cptac_2021" && has_samp && !is.na(sid) && "ESTIMATE_TUMORPURITY" %in% names(samp)) {
    purity_col <- "SAMPLE:ESTIMATE_TUMORPURITY"
    purity_vec[] <- .to01(samp[["ESTIMATE_TUMORPURITY"]][match(sample_ids, samp[[sid]])])
  } else if (ds == "paad_cptac_2021" && has_samp && !is.na(sid) && "NEOPLASTIC_CELLULARITY" %in% names(samp)) {
    purity_col <- "SAMPLE:NEOPLASTIC_CELLULARITY(median;)"
    medv <- vapply(as.character(samp[["NEOPLASTIC_CELLULARITY"]]), .median_from_semicolon, numeric(1))
    purity_vec[] <- .to01(medv[match(sample_ids, samp[[sid]])])
  } else if (ds == "ucec_cptac_2020" && has_samp && !is.na(sid)) {
    if ("PURITY_CANCER" %in% names(samp)) {
      purity_col <- "SAMPLE:PURITY_CANCER"
      purity_vec[] <- .to01(samp[["PURITY_CANCER"]][match(sample_ids, samp[[sid]])])
    } else {
      # Try inferring from immune/stroma (if columns exist)
      pi <- if ("PURITY_IMMUNE" %in% names(samp)) suppressWarnings(as.numeric(samp[["PURITY_IMMUNE"]])) else NULL
      ps <- if ("PURITY_STROMA" %in% names(samp)) suppressWarnings(as.numeric(samp[["PURITY_STROMA"]])) else NULL
      if (!is.null(pi) || !is.null(ps)) {
        idx <- match(sample_ids, samp[[sid]])
        purity_est <- 1 - (pi[idx] %||% 0) - (ps[idx] %||% 0)
        purity_vec[] <- .to01(purity_est)
        purity_col <- "SAMPLE:1-(PURITY_IMMUNE+PURITY_STROMA)"
      }
    }
  } else {
    purity_col <- "NONE"
  }

  nonNA <- mean(is.finite(purity_vec)) * 100
  pur_min <- if (any(is.finite(purity_vec))) min(purity_vec, na.rm = TRUE) else NA_real_
  pur_med <- if (any(is.finite(purity_vec))) stats::median(purity_vec, na.rm = TRUE) else NA_real_
  pur_max <- if (any(is.finite(purity_vec))) max(purity_vec, na.rm = TRUE) else NA_real_

  list(
    column = purity_col,
    nonNA_pct = nonNA,
    min = pur_min, median = pur_med, max = pur_max,
    vec = purity_vec
  )
}

## ==============================
## Small Audit: sex/age + purity + batch overall check (Improved version)
## -- Added: policy/batch_min/limma_min in log and output columns
## ==============================
audit_one_dataset_sa_batch <- function(ds_dir, pipe_policy = "NA", min_per_level = 2) {
  ds_id <- basename(ds_dir)
  # 1) Get protein matrix samples
  m <- load_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(m)

  # 2) Fetch SEX/AGE from patient file only
  samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
  has_samp <- file.exists(samp_fp)
  has_pat <- file.exists(pat_fp)

  sex_nonNA <- age_nonNA <- NA_real_
  sex_M <- sex_F <- NA_integer_
  age_min <- age_med <- age_max <- NA_real_
  sex_col <- age_col <- "MISSING"
  sex_src <- age_src <- if (has_pat) "patient" else "NONE"

  if (has_samp && has_pat) {
    samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    names(samp) <- .norm_names(names(samp))
    names(pat) <- .norm_names(names(pat))

    sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
    pid_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
    pid_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]

    if (!is.na(sid) && !is.na(pid_samp) && !is.na(pid_pat)) {
      map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])

      # SEX
      if ("SEX" %in% names(pat)) {
        sex_col <- "SEX"
        raw <- toupper(as.character(pat$SEX))
        val <- ifelse(grepl("^M", raw), 1, ifelse(grepl("^F", raw), 0, NA_real_))
        names(val) <- as.character(pat[[pid_pat]])
        sex_vec <- unname(val[map_pt[sample_ids]])
        sex_nonNA <- mean(is.finite(sex_vec)) * 100
        sex_M <- sum(sex_vec == 1, na.rm = TRUE)
        sex_F <- sum(sex_vec == 0, na.rm = TRUE)
      }

      # AGE (Original value stats; z-score done in main pipeline)
      if ("AGE" %in% names(pat)) {
        age_col <- "AGE"
        v <- suppressWarnings(as.numeric(pat$AGE))
        names(v) <- as.character(pat[[pid_pat]])
        age_vec <- unname(v[map_pt[sample_ids]])
        age_nonNA <- mean(is.finite(age_vec)) * 100
        if (any(is.finite(age_vec))) {
          age_min <- min(age_vec, na.rm = TRUE)
          age_med <- stats::median(age_vec, na.rm = TRUE)
          age_max <- max(age_vec, na.rm = TRUE)
        }
      }
    }
  }

  # 3) purity audit
  pur <- .audit_purity_for_dataset(ds_id, ds_dir, sample_ids)

  # 4) batch check (Use your existing get_batch_factor)
  bi <- get_batch_factor(ds_dir, sample_ids,
    pipe_policy   = pipe_policy,
    min_per_level = min_per_level
  )
  batch_col <- if (!is.null(bi)) bi$name else "NONE"
  batch_levels <- if (!is.null(bi)) nlevels(bi$fac) else 0L
  batch_nonNA <- if (!is.null(bi)) mean(!is.na(bi$fac)) * 100 else NA_real_
  batch_sizes <- if (!is.null(bi)) .format_batch_sizes(bi$fac) else NA_character_

  # Print summary line (Append policy/batch_min/limma_min)
  line <- sprintf(
    "[Audit] %-24s | SEX:%-7s nonNA=%5.1f%% (M=%d,F=%d) | AGE:%-7s nonNA=%5.1f%% (min/med/max=%s/%s/%s) | PURITY:%-35s nonNA=%5.1f%% (min/med/max=%s/%s/%s) | BATCH:%-18s levels=%-2s nonNA=%5.1f%% [%s] | policy=%s | batch_min=%d | limma_min=%d",
    ds_id,
    sex_col, sex_nonNA %||% NaN, sex_M %||% NA_integer_, sex_F %||% NA_integer_,
    age_col, age_nonNA %||% NaN,
    ifelse(is.finite(age_min), round(age_min, 1), "NA"),
    ifelse(is.finite(age_med), round(age_med, 1), "NA"),
    ifelse(is.finite(age_max), round(age_max, 1), "NA"),
    pur$column, pur$nonNA_pct %||% NaN,
    ifelse(is.finite(pur$min), round(pur$min, 3), "NA"),
    ifelse(is.finite(pur$median), round(pur$median, 3), "NA"),
    ifelse(is.finite(pur$max), round(pur$max, 3), "NA"),
    batch_col, as.integer(batch_levels), batch_nonNA %||% NaN, batch_sizes %||% "NA",
    pipe_policy, min_per_level, min_per_group
  )
  if (exists("log_msg", mode = "function")) log_msg(line) else cat(line, "\n")

  data.frame(
    dataset = ds_id,
    sex_column = sex_col, sex_nonNA_pct = round(sex_nonNA, 1), sex_M = sex_M, sex_F = sex_F,
    age_column = age_col, age_nonNA_pct = round(age_nonNA, 1),
    age_min = ifelse(is.finite(age_min), round(age_min, 1), NA),
    age_median = ifelse(is.finite(age_med), round(age_med, 1), NA),
    age_max = ifelse(is.finite(age_max), round(age_max, 1), NA),
    purity_column = pur$column, purity_nonNA_pct = round(pur$nonNA_pct, 1),
    purity_min = ifelse(is.finite(pur$min), round(pur$min, 3), NA),
    purity_median = ifelse(is.finite(pur$median), round(pur$median, 3), NA),
    purity_max = ifelse(is.finite(pur$max), round(pur$max, 3), NA),
    batch_column = batch_col, batch_levels = as.integer(batch_levels),
    batch_nonNA_pct = round(batch_nonNA, 1), batch_sizes = batch_sizes,
    batch_policy = pipe_policy,
    batch_min_per_level = as.integer(min_per_level),
    limma_min_per_group = as.integer(min_per_group),
    stringsAsFactors = FALSE
  )
}


audit_all_datasets_sa_batch <- function(dataset_dirs, pipe_policy = "NA", min_per_level = 2) {
  out <- lapply(dataset_dirs, function(p) {
    tryCatch(audit_one_dataset_sa_batch(p, pipe_policy, min_per_level),
      error = function(e) {
        if (exists("log_msg", mode = "function")) log_msg("[Audit] %s ERROR: %s", basename(p), e$message)
        data.frame(
          dataset = basename(p),
          sex_column = "ERROR", sex_nonNA_pct = NA, sex_M = NA, sex_F = NA,
          age_column = "ERROR", age_nonNA_pct = NA, age_min = NA, age_median = NA, age_max = NA,
          purity_column = "ERROR", purity_nonNA_pct = NA, purity_min = NA, purity_median = NA, purity_max = NA,
          batch_column = "ERROR", batch_levels = NA, batch_nonNA_pct = NA, batch_sizes = NA,
          stringsAsFactors = FALSE
        )
      }
    )
  })
  dplyr::bind_rows(out)
}




## ===== DO_AUDIT switch for QC/Exploration (4) =====

DO_AUDIT <- FALSE # <— Default OFF; change to TRUE when needed

if (DO_AUDIT) {
  results <- list()
  for (ds in names(dataset_dirs)) {
    ds_dir <- dataset_dirs[[ds]]
    if (!dir.exists(ds_dir)) {
      log_msg("Skip: Folder %s not found", ds_dir)
      next
    }
    results[[ds]] <- screen_batch_need(ds_dir)
  }

  ds_dir <- file.path(getwd(), "brca_cptac_2020")
  if (dir.exists(ds_dir)) try(inspect_plex(ds_dir), silent = TRUE)

  sex_age_batch_audit <- audit_all_datasets_sa_batch(dataset_dirs)
  print(sex_age_batch_audit)
  readr::write_csv(sex_age_batch_audit, file.path(getwd(), "sex_age_batch_audit_summary.csv"))
}


yaml::write_yaml(list(
  seed = 1234,
  min_frac_complete = min_frac_complete,
  hi_lo_quantile = hi_lo_quantile,
  min_per_group = min_per_group,
  MAKE_PLOTS = MAKE_PLOTS,
  datasets_root = datasets_root,
  dataset_ids = dataset_ids,
  strata = strata,
  geneset_groups_selected = GENESET_GROUPS_TO_RUN
), file.path("run_info", "run_manifest.yml"))


## Start rerun from which dataset (Changeable)
start_from <- "brca_cptac_2020"
ord <- dataset_ids
ix <- match(start_from, ord)
if (is.na(ix)) stop(sprintf("'%s' not in dataset_ids", start_from))
dataset_dirs_run <- dataset_dirs[ord[ix:length(ord)]]

exists("coerce_covariates_safely")
getAnywhere("coerce_covariates_safely")
exists("opt") # Should return TRUE
getAnywhere("opt") # Should show in .GlobalEnv


## ---- Optional: For testing, run specific strata only ----
# only_strata <- c("TP53_mutant")   # To run multiple, use c("ALL","TP53_mutant") etc.

## Helper: If only_strata not set, default run all
.should_run <- function(tag) {
  if (!exists("only_strata") || is.null(only_strata) || !length(only_strata)) {
    return(TRUE)
  }
  tag %in% only_strata
}


## =========================================================
## Process datasets to run in this round sequentially (Complete main loop + Define sample sets for each stratum)
## =========================================================
for (ds in names(dataset_dirs_run)) {
  ds_dir <- dataset_dirs_run[[ds]]
  log_msg("== Start Dataset: %s ==", ds)


  ## Protein Matrix & TP53 Status
  mat0_full <- load_matrix_from_dataset_dir(ds_dir)
  tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full))

  ## Sample sets for three strata
  samples_ALL <- colnames(mat0_full)
  samples_MUT <- names(tp53_status)[tp53_status == "TP53_mutant"]
  samples_WT <- names(tp53_status)[tp53_status == "TP53_wild_type"]

  ## Output root directory for each stratum
  base_tp53_root <- file.path(OUTPUT_PREFIX, ds, "csn_gsea_results_TP53")


  ## Run three strata sequentially (RAW & BatchAdj will be written in run_one_stratum)
  if (.should_run("ALL")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full,
      sample_keep = samples_ALL,
      out_root = file.path(base_tp53_root, "ALL"),
      genesets_by_group = genesets_by_group
    )
  }
  ## [NEW] ALL: Aggregate DEG across predictors
  try(summarize_deg_across_predictors(file.path(base_tp53_root, "ALL")), silent = TRUE)

  if (.should_run("TP53_mutant")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full,
      sample_keep = samples_MUT,
      out_root = file.path(base_tp53_root, "TP53_mutant"),
      genesets_by_group = genesets_by_group
    )
  }
  ## [NEW] TP53_mutant: Aggregate DEG across predictors
  try(summarize_deg_across_predictors(file.path(base_tp53_root, "TP53_mutant")), silent = TRUE)

  if (.should_run("TP53_wild_type")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full,
      sample_keep = samples_WT,
      out_root = file.path(base_tp53_root, "TP53_wild_type"),
      genesets_by_group = genesets_by_group
    )
  }
  ## [NEW] TP53_wild_type: Aggregate DEG across predictors
  try(summarize_deg_across_predictors(file.path(base_tp53_root, "TP53_wild_type")), silent = TRUE)

  log_msg("== Completed Dataset: %s (TP53 Stratum Output -> %s) ==", ds, base_tp53_root)
}


## ---------- [NEW] WT vs MT Delta NES Aggregation (Aggregate from existing outputs) ----------
tp53_delta_nes_aggregate <- function(datasets_root,
                                     dataset_ids = NULL,
                                     versions = c("BatchAdj"),
                                     groups = names(genesets_by_group),
                                     stat_tags = c("GSEA_limma_t_cont", "GSEA_spearman")) {
  if (is.null(dataset_ids)) {
    dataset_ids <- list.dirs(datasets_root, full.names = FALSE, recursive = FALSE)
  }
  sfn <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  for (ds in dataset_ids) {
    base_dir <- file.path(datasets_root, ds, "csn_gsea_results_TP53")
    for (ver in versions) {
      # Auto-detect subunits: Read subdirectory names under WT/MT folders
      subunits <- unique(unlist(lapply(groups, function(g) {
        b1 <- file.path(OUTPUT_PREFIX, "csn_gsea_results_TP53_by_collection", safe_fs_name(g), ver, ds, "TP53_wild_type")
        b2 <- file.path(OUTPUT_PREFIX, "csn_gsea_results_TP53_by_collection", safe_fs_name(g), ver, ds, "TP53_mutant")
        c(
          list.dirs(b1, full.names = FALSE, recursive = FALSE),
          list.dirs(b2, full.names = FALSE, recursive = FALSE)
        )
      })))
      subunits <- subunits[nzchar(subunits)]
      subunits <- subunits[nzchar(subunits)]
      if (!length(subunits)) next
      for (su in subunits) {
        for (grp in groups) {
          grp_safe <- safe_fs_name(grp)
          for (st in stat_tags) {
            fp_wt <- file.path(OUTPUT_PREFIX, "csn_gsea_results_TP53_by_collection", grp_safe, ver, ds, "TP53_wild_type", su, paste0(st, ".csv"))
            fp_mt <- file.path(OUTPUT_PREFIX, "csn_gsea_results_TP53_by_collection", grp_safe, ver, ds, "TP53_mutant", su, paste0(st, ".csv"))
            if (!file.exists(fp_wt) || !file.exists(fp_mt)) next
            wt <- tryCatch(data.table::fread(fp_wt, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
            mt <- tryCatch(data.table::fread(fp_mt, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
            if (is.null(wt) || is.null(mt) || !nrow(wt) || !nrow(mt)) next
            wt <- setNames(as.data.frame(wt), tolower(names(wt)))
            mt <- setNames(as.data.frame(mt), tolower(names(mt)))
            # Normalize column names
            if (!("pathway" %in% names(wt)) && ("term" %in% names(wt))) names(wt)[names(wt) == "term"] <- "pathway"
            if (!("pathway" %in% names(mt)) && ("term" %in% names(mt))) names(mt)[names(mt) == "term"] <- "pathway"
            if (!("nes" %in% names(wt)) && ("enrichment_score" %in% names(wt))) names(wt)[names(wt) == "enrichment_score"] <- "nes"
            if (!("nes" %in% names(mt)) && ("enrichment_score" %in% names(mt))) names(mt)[names(mt) == "enrichment_score"] <- "nes"
            for (nm in c("nes", "padj", "pval")) if (!nm %in% names(wt)) wt[[nm]] <- NA_real_
            for (nm in c("nes", "padj", "pval")) if (!nm %in% names(mt)) mt[[nm]] <- NA_real_
            df <- dplyr::full_join(wt[, c("pathway", "nes", "padj", "pval")],
              mt[, c("pathway", "nes", "padj", "pval")],
              by = "pathway", suffix = c("_WT", "_MT")
            )
            if (!nrow(df)) next
            df <- df |>
              dplyr::mutate(
                NES_WT = as.numeric(nes_WT),
                NES_MT = as.numeric(nes_MT),
                padj_WT = as.numeric(padj_WT),
                padj_MT = as.numeric(padj_MT),
                pval_WT = as.numeric(pval_WT),
                pval_MT = as.numeric(pval_MT),
                delta_NES = NES_MT - NES_WT,
                direction_consistency = dplyr::case_when(
                  is.finite(NES_WT) & is.finite(NES_MT) & NES_MT > 0 & NES_WT > 0 ~ "BOTH_POS",
                  is.finite(NES_WT) & is.finite(NES_MT) & NES_MT < 0 & NES_WT < 0 ~ "BOTH_NEG",
                  is.finite(NES_WT) & is.finite(NES_MT) ~ "MIXED",
                  TRUE ~ "NA"
                )
              ) |>
              dplyr::select(pathway, NES_WT, NES_MT, delta_NES, padj_WT, padj_MT, pval_WT, pval_MT, direction_consistency) |>
              dplyr::arrange(dplyr::desc(abs(delta_NES)), pathway)
            out_dir <- file.path(OUTPUT_PREFIX, "csn_gsea_results_TP53_by_collection", grp_safe, ver, ds, "DeltaWT_MT", su)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            out_csv <- file.path(out_dir, paste0(st, "_deltaWTMT.csv"))
            data.table::fwrite(df, out_csv)
          }
        }
      }
      message(sprintf("[DeltaNES] %s/%s Done", ds, ver))
    }
  }
  invisible(TRUE)
}



setwd("C:/Users/danny/Documents/R_project/CSN_CPTAC") ## YTC laptop
## setwd("C:/Users/cmuh/Documents/YenTing_document/CSN_CPTAC") ## lab computer
getwd()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
})


.read_csv_safe <- function(path) {
  p <- path

  # Windows: If original string not found, try replacing / with \
  if (!file.exists(p) && .Platform$OS.type == "windows") {
    p2 <- gsub("/", "\\\\", p, fixed = TRUE)
    if (file.exists(p2)) p <- p2
  }
  if (!file.exists(p)) stop("File not found: ", p)

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

  # 3) Do not open connection manually: Read as full text string then parse
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




.safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

# Target X-axis order (Inherit your original script)
.pred_order_all <- c(
  "NES_CSN_SCORE", "NES_GPS1",
  "NES_COPS2", "NES_COPS3", "NES_COPS4", "NES_COPS5", "NES_COPS6",
  "NES_COPS7A", "NES_COPS7B", "NES_COPS8", "NES_COPS9",
  "NES_RESIDUAL_GPS1",
  "NES_RESIDUAL_COPS2", "NES_RESIDUAL_COPS3", "NES_RESIDUAL_COPS4", "NES_RESIDUAL_COPS5", "NES_RESIDUAL_COPS6",
  "NES_RESIDUAL_COPS7A", "NES_RESIDUAL_COPS7B", "NES_RESIDUAL_COPS8", "NES_RESIDUAL_COPS9"
)

# Single CSV -> Long table -> ggplot object (No file output)
.make_heatmap_plot <- function(csv_file) {
  df_raw <- suppressMessages(.read_csv_safe(csv_file))

  # Auto-detect pathway column (Inherit your original script)
  path_candidates <- c("pathway", "Pathway", "term", "Term", "gs_name", "NAME", "set", "Set")
  path_col <- intersect(path_candidates, names(df_raw))[1]
  if (is.na(path_col)) stop("Pathway column not found: ", paste(path_candidates, collapse = ", "))

  # Long/Wide table determination + Convert to long table (Inherit)
  has_long <- all(c("predictor", "NES", "padj") %in% names(df_raw))
  if (has_long) {
    df_long <- df_raw %>%
      rename(pathway = all_of(path_col)) %>%
      mutate(
        predictor = as.character(predictor),
        NES = as.numeric(NES),
        padj = as.numeric(padj)
      )
  } else {
    nes_cols <- grep("^NES_", names(df_raw), value = TRUE)
    padj_cols <- grep("^padj_", names(df_raw), value = TRUE)
    if (!length(nes_cols)) stop("Wide format: NES_* columns not found.")

    df_wide <- df_raw %>% rename(pathway = all_of(path_col))
    nes_map <- tibble(nes_col = nes_cols, key = sub("^NES_", "", nes_cols))
    padj_map <- tibble(padj_col = padj_cols, key = sub("^padj_", "", padj_cols))
    pair_map <- nes_map %>% left_join(padj_map, by = "key")

    df_list <- lapply(seq_len(nrow(pair_map)), function(i) {
      nes_c <- pair_map$nes_col[i]
      padj_c <- pair_map$padj_col[i] # May be NA
      tibble(
        pathway   = df_wide$pathway,
        predictor = paste0("NES_", pair_map$key[i]),
        NES       = suppressWarnings(as.numeric(df_wide[[nes_c]])),
        padj      = if (!is.na(padj_c)) suppressWarnings(as.numeric(df_wide[[padj_c]])) else NA_real_
      )
    })
    df_long <- bind_rows(df_list)
  }

  # Keep only predictors in specified order, and sort by order
  present_preds <- intersect(.pred_order_all, unique(df_long$predictor))
  if (!length(present_preds)) stop("No specified predictors found in data.")
  df_long <- df_long %>% filter(predictor %in% present_preds)
  ## [NEW] Non-H collection: Filter top/bottom N by CSN_SCORE NES
  .coll_tok <- sub("^Summary_(.+?)_GSEA.*$", "\\1", basename(csv_file))
  if (!identical(toupper(.coll_tok), "H") && "NES_CSN_SCORE" %in% present_preds) {
    .top_n <- get0("DATASET_HEATMAP_TOP_N", ifnotfound = 25)
    .bot_n <- get0("DATASET_HEATMAP_BOTTOM_N", ifnotfound = 25)
    csn_tbl <- df_long %>%
      dplyr::filter(.data$predictor == "NES_CSN_SCORE")
    keep_up <- csn_tbl %>%
      dplyr::arrange(dplyr::desc(.data$NES)) %>%
      dplyr::slice_head(n = .top_n) %>%
      dplyr::pull(.data$pathway)
    keep_dn <- csn_tbl %>%
      dplyr::arrange(.data$NES) %>%
      dplyr::slice_head(n = .bot_n) %>%
      dplyr::pull(.data$pathway)
    keep <- unique(c(keep_up, keep_dn))
    if (length(keep)) {
      df_long <- df_long %>% dplyr::filter(.data$pathway %in% keep)
    }
  }

  # Y-axis sort: Prioritize NES_CSN_SCORE descending; otherwise by row mean
  if ("NES_CSN_SCORE" %in% present_preds) {
    nes_csn <- df_long %>%
      filter(predictor == "NES_CSN_SCORE") %>%
      select(pathway, NES) %>%
      distinct()
    path_order <- nes_csn %>%
      arrange(desc(NES)) %>%
      pull(pathway)
  } else {
    path_order <- df_long %>%
      group_by(pathway) %>%
      summarise(m = mean(NES, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(m)) %>%
      pull(pathway)
  }
  df_long <- df_long %>% mutate(pathway = factor(pathway, levels = rev(path_order)))


  # X-axis unequal spacing and two gaps
  # Rule adjustment: Even if "No NES_COPS9", keep second gap as long as NES_RESIDUAL_GPS1 exists
  gap <- 0.4
  needs_gap1 <- all(c("NES_CSN_SCORE", "NES_GPS1") %in% present_preds)
  needs_gap2 <- "NES_RESIDUAL_GPS1" %in% present_preds # Modified: No longer require COPS9 simultaneously

  pos_map <- list()
  pos <- 0
  for (p in present_preds) {
    if (p == "NES_GPS1" && needs_gap1) pos <- pos + gap
    if (p == "NES_RESIDUAL_GPS1" && needs_gap2) pos <- pos + gap
    pos <- pos + 1
    pos_map[[p]] <- pos
  }
  pos_map <- unlist(pos_map)
  df_long <- df_long %>% mutate(xpos = unname(pos_map[predictor]))

  # Color range (Symmetric 0) and cell style color
  L <- max(abs(df_long$NES), na.rm = TRUE)
  if (!is.finite(L) || L == 0) L <- 1
  col_low <- "#053061"
  col_mid <- "#FFFFFF"
  col_high <- "#67001F"

  # Plotting (Black dot: padj < 0.05)
  p <- ggplot(df_long, aes(x = xpos, y = pathway, fill = NES)) +
    geom_tile(width = 1, height = 0.9, color = NA) +
    geom_point(
      data = df_long %>% filter(is.finite(padj), padj < 0.05),
      aes(x = xpos, y = pathway),
      shape = 16, size = 1.6, color = "black", inherit.aes = FALSE
    ) +
    scale_fill_gradient2(
      low = col_low, mid = col_mid, high = col_high,
      limits = c(-L, L), oob = scales::squish, name = "NES"
    ) +
    scale_x_continuous(
      breaks = unname(pos_map[present_preds]),
      labels = present_preds,
      expand = expansion(mult = c(0.01, 0.01)),
      position = "top"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 9, margin = margin(b = 8)),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.margin = margin(6, 12, 6, 6),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    ) +
    coord_cartesian(clip = "off")

  # Size: Adjust with pathway count
  n_path <- nlevels(df_long$pathway)
  W <- max(8, length(present_preds) * 0.45) # inches
  H <- max(6, n_path * 0.22) # inches
  list(plot = p, width = W, height = H)
}


# Infer dataset / variant / stratum from path
.parse_meta_from_path <- function(csv_file) {
  parts <- strsplit(normalizePath(csv_file, winslash = "/"), "/")[[1]]
  parts_lc <- tolower(parts)
  bn <- basename(csv_file)

  # Variant (RAW / BatchAdj): Grab any layer
  vidx <- which(parts_lc %in% c("raw", "batchadj"))
  variant <- if (length(vidx)) parts[vidx[1]] else NA_character_

  # Stratum real location: "One layer above" variant
  stratum <- NA_character_
  if (length(vidx) && vidx[1] - 1 >= 1) {
    cand <- parts[vidx[1] - 1]
    if (tolower(cand) %in% c("all", "tp53_mutant", "tp53_wild_type")) {
      stratum <- cand
    }
  }
  # Fallback: If not found in previous step, search once in any layer
  if (is.na(stratum) || !nzchar(stratum)) {
    hit <- which(parts_lc %in% c("all", "tp53_mutant", "tp53_wild_type"))
    if (length(hit)) stratum <- parts[hit[1]]
  }

  # dataset: If path contains csn_gsea_results_TP53, take its parent; otherwise fallback to caller
  dataset <- NA_character_
  hit2 <- which(parts_lc == "csn_gsea_results_tp53")
  if (length(hit2) && hit2[1] > 1) dataset <- parts[hit2[1] - 1]

  list(dataset = dataset, variant = variant, stratum = stratum)
}


# Actual output: Four formats (In-place & Aggregated)
.save_both_places <- function(p, width, height, csv_file, meta,
                              collect_root = "single_dataset_GSEA_heatmap", dpi = 600) {
  bn <- basename(csv_file) # e.g. Summary_H_GSEA_limma_t_cont_ALL.csv
  dir_near <- dirname(csv_file)

  # 1) In-place filename
  near_prefix <- paste0(
    "heatmap_", .safe_fs(meta$dataset), "_", .safe_fs(meta$stratum),
    "_", .safe_fs(meta$variant), "_"
  )
  near_base <- file.path(dir_near, paste0(near_prefix, bn))

  # 2) Aggregated path
  out_dir2 <- file.path(collect_root, .safe_fs(meta$dataset))
  dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)
  collect_base <- file.path(out_dir2, paste0("heatmap_", bn))

  # Write four formats
  ggsave(paste0(near_base, ".tiff"), p,
    width = width, height = height, units = "in",
    dpi = dpi, bg = "white", compression = "lzw"
  )
  # ggsave(paste0(near_base, ".png"),  p, width = width, height = height, units = "in",
  #       dpi = dpi, bg = "white")
  # ggsave(paste0(near_base, ".jpg"),  p, width = width, height = height, units = "in",
  #       dpi = dpi, bg = "white", quality = 100)
  # ggsave(paste0(near_base, ".pdf"),  p, width = width, height = height, units = "in",
  #       bg = "white")

  ggsave(paste0(collect_base, ".tiff"), p,
    width = width, height = height, units = "in",
    dpi = dpi, bg = "white", compression = "lzw"
  )
  # ggsave(paste0(collect_base, ".png"),  p, width = width, height = height, units = "in",
  #       dpi = dpi, bg = "white")
  # ggsave(paste0(collect_base, ".jpg"),  p, width = width, height = height, units = "in",
  #       dpi = dpi, bg = "white", quality = 100)
  # ggsave(paste0(collect_base, ".pdf"),  p, width = width, height = height, units = "in",
  #       bg = "white")
}




if (!exists(".read_csv_safe", mode = "function")) {
  .read_csv_safe <- function(path) {
    p <- path
    # Windows: If original string not found, try / -> \
    if (!file.exists(p) && .Platform$OS.type == "windows") {
      p2 <- gsub("/", "\\\\", p, fixed = TRUE)
      if (file.exists(p2)) p <- p2
    }
    if (!file.exists(p)) stop("File not found: ", p)

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

    # 3) Read as full text string then parse
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
}

if (!exists(".safe_fs", mode = "function")) {
  .safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

# Target X-axis order (Updated: Use *_adj_CSN instead of RESIDUAL_*)
# IMPORTANT: Always redefine to ensure latest naming is used
.pred_order_all <- c(
  "NES_CSN_SCORE", "NES_GPS1",
  "NES_COPS2", "NES_COPS3", "NES_COPS4", "NES_COPS5", "NES_COPS6",
  "NES_COPS7A", "NES_COPS7B", "NES_COPS8", "NES_COPS9",
  "NES_GPS1_adj_CSN",
  "NES_COPS2_adj_CSN", "NES_COPS3_adj_CSN", "NES_COPS4_adj_CSN", "NES_COPS5_adj_CSN", "NES_COPS6_adj_CSN",
  "NES_COPS7A_adj_CSN", "NES_COPS7B_adj_CSN", "NES_COPS8_adj_CSN", "NES_COPS9_adj_CSN"
)


# Z / padj_meta predictor order (Inherit your original NES_* order)
.pred_order_meta <- gsub("^NES_", "Z_", .pred_order_all)

# Read Summary_*_meta_fdr_ALL.csv as "Long Table" (pathway, predictor, Z, padj_meta)
.meta_read_summary_long <- function(csv_file) {
  df <- suppressMessages(.read_csv_safe(csv_file))
  # Auto find pathway column
  path_candidates <- c("pathway", "Pathway", "term", "Term", "gs_name", "NAME", "set", "Set")
  path_col <- intersect(path_candidates, names(df))[1]
  if (is.na(path_col)) stop("Pathway column not found in: ", csv_file)

  z_cols <- grep("^Z_", names(df), value = TRUE)
  padj_cols <- grep("^padj_meta_", names(df), value = TRUE)
  if (!length(z_cols)) stop("Z_* columns not found in: ", csv_file)

  df <- dplyr::rename(df, pathway = dplyr::all_of(path_col))
  z_map <- tibble::tibble(z_col = z_cols, key = sub("^Z_", "", z_cols))
  padj_map <- tibble::tibble(p_col = padj_cols, key = sub("^padj_meta_", "", padj_cols))
  pair_map <- dplyr::left_join(z_map, padj_map, by = "key")

  lst <- lapply(seq_len(nrow(pair_map)), function(i) {
    zc <- pair_map$z_col[i]
    pc <- pair_map$p_col[i]
    tibble::tibble(
      pathway   = df$pathway,
      predictor = paste0("Z_", pair_map$key[i]),
      Z         = suppressWarnings(as.numeric(df[[zc]])),
      padj_meta = if (!is.na(pc)) suppressWarnings(as.numeric(df[[pc]])) else NA_real_
    )
  })
  dplyr::bind_rows(lst)
}

# Calculate y-axis order using meta file (Z): Prioritize Z_CSN_SCORE descending, otherwise column mean Z
.meta_compute_y_order <- function(csv_file) {
  df_long <- .meta_read_summary_long(csv_file)
  present <- intersect(.pred_order_meta, unique(df_long$predictor))
  if (!length(present)) stop("No available predictors in: ", csv_file)
  df_long <- dplyr::filter(df_long, .data$predictor %in% present)

  ## -- Key: Use "Top 25 + Bottom 25 (by Z_CSN_SCORE)" of ALL layer as unique geneset list -- ##
  .coll_tok <- sub("^Summary_(.+?)_GSEA.*$", "\\1", basename(csv_file))
  if ("Z_CSN_SCORE" %in% present) {
    z_tbl <- df_long %>%
      dplyr::filter(.data$predictor == "Z_CSN_SCORE") %>%
      dplyr::select(.data$pathway, .data$Z) %>%
      dplyr::distinct()
    .top_n <- get0("PAN_HEATMAP_TOP_N", ifnotfound = 25L)
    .bot_n <- get0("PAN_HEATMAP_BOTTOM_N", ifnotfound = 25L)
    if (!identical(toupper(.coll_tok), "H")) {
      keep_up <- z_tbl %>%
        dplyr::arrange(dplyr::desc(.data$Z)) %>%
        dplyr::slice_head(n = .top_n) %>%
        dplyr::pull(.data$pathway)
      keep_dn <- z_tbl %>%
        dplyr::arrange(.data$Z) %>%
        dplyr::slice_head(n = .bot_n) %>%
        dplyr::pull(.data$pathway)
      keep <- unique(c(keep_up, keep_dn))
      if (length(keep)) {
        df_long <- df_long %>% dplyr::filter(.data$pathway %in% keep)
        z_tbl <- z_tbl %>% dplyr::filter(.data$pathway %in% keep)
      }
    } else {
      keep <- z_tbl$pathway
      df_long <- df_long %>% dplyr::filter(.data$pathway %in% keep)
      z_tbl <- z_tbl %>% dplyr::filter(.data$pathway %in% keep)
    }
  }

  # Give "Order containing only keep" by Z_CSN_SCORE descending (or row mean Z)
  if ("Z_CSN_SCORE" %in% present) {
    ord <- df_long %>%
      dplyr::filter(.data$predictor == "Z_CSN_SCORE") %>%
      dplyr::select(.data$pathway, .data$Z) %>%
      dplyr::distinct() %>%
      dplyr::arrange(dplyr::desc(.data$Z)) %>%
      dplyr::pull(.data$pathway)
  } else {
    ord <- df_long %>%
      dplyr::group_by(.data$pathway) %>%
      dplyr::summarise(m = mean(.data$Z, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(.data$m)) %>%
      dplyr::pull(.data$pathway)
  }
  ord
}


# Build meta-FDR heatmap (Z as color scale; padj_meta<0.05 black dot)
.meta_make_heatmap_plot <- function(csv_file, y_order = NULL,
                                    palette = c(low = "#053061", mid = "#FFFFFF", high = "#67001F")) {
  df_long <- .meta_read_summary_long(csv_file)

  present <- intersect(.pred_order_meta, unique(df_long$predictor))
  if (!length(present)) stop("No specified predictors found in data: ", csv_file)
  df_long <- dplyr::filter(df_long, .data$predictor %in% present)
  ## [NEW] Non-H collection: Filter top/bottom N by Z of CSN_SCORE
  ##      **Important**: If y_order is not NULL (i.e. ordered version), do not perform any top/bottom filtering here,
  ##      only use y_order to precisely specify 50 genesets later.
  .coll_tok <- sub("^Summary_(.+?)_GSEA.*$", "\\1", basename(csv_file))
  if (is.null(y_order)) {
    if (!identical(toupper(.coll_tok), "H") && "Z_CSN_SCORE" %in% present) {
      .top_n <- get0("PAN_HEATMAP_TOP_N", ifnotfound = 25)
      .bot_n <- get0("PAN_HEATMAP_BOTTOM_N", ifnotfound = 25)
      csn_tbl <- df_long %>%
        dplyr::filter(.data$predictor == "Z_CSN_SCORE")
      keep_up <- csn_tbl %>%
        dplyr::arrange(dplyr::desc(.data$Z)) %>%
        dplyr::slice_head(n = .top_n) %>%
        dplyr::pull(.data$pathway)
      keep_dn <- csn_tbl %>%
        dplyr::arrange(.data$Z) %>%
        dplyr::slice_head(n = .bot_n) %>%
        dplyr::pull(.data$pathway)
      keep <- unique(c(keep_up, keep_dn))
      if (length(keep)) {
        df_long <- dplyr::filter(df_long, .data$pathway %in% keep)
      }
    }
  }

  # y-order: Use if provided externally, otherwise by Z_CSN_SCORE or mean Z
  if (is.null(y_order)) {
    y_order <- .meta_compute_y_order(csv_file)
  }
  # Use ALL order; plot only ALL genesets and fully inherit its order (Same as single dataset ordered)
  df_long <- df_long %>% dplyr::filter(.data$pathway %in% y_order)
  y_levels <- y_order[y_order %in% df_long$pathway]
  if (!length(y_levels)) stop("y-order has no intersection with data: ", csv_file)
  df_long <- dplyr::mutate(df_long, pathway = factor(.data$pathway, levels = rev(y_levels)))


  # X-axis position and gap rules (Consistent with your original rules; change to Z_* names)
  gap <- 0.4
  needs_gap1 <- all(c("Z_CSN_SCORE", "Z_GPS1") %in% present)
  needs_gap2 <- "Z_RESIDUAL_GPS1" %in% present

  pos_map <- list()
  pos <- 0
  for (p in present) {
    if (p == "Z_GPS1" && needs_gap1) pos <- pos + gap
    if (p == "Z_RESIDUAL_GPS1" && needs_gap2) pos <- pos + gap
    pos <- pos + 1
    pos_map[[p]] <- pos
  }
  pos_map <- unlist(pos_map)
  df_long <- dplyr::mutate(df_long, xpos = unname(pos_map[.data$predictor]))

  L <- max(abs(df_long$Z), na.rm = TRUE)
  if (!is.finite(L) || L == 0) L <- 1

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = pathway, fill = Z)) +
    ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
    ggplot2::geom_point(
      data = dplyr::filter(df_long, is.finite(.data$padj_meta), .data$padj_meta < 0.05),
      ggplot2::aes(x = xpos, y = pathway),
      shape = 16, size = 1.6, color = "black", inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = palette[["low"]], mid = palette[["mid"]], high = palette[["high"]],
      limits = c(-L, L), midpoint = 0, oob = scales::squish, name = "Z"
    ) +
    ggplot2::scale_x_continuous(
      breaks = unname(pos_map[present]),
      labels = present,
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
  W <- max(8, length(present) * 0.45)
  H <- max(6, n_path * 0.22)
  list(plot = p, width = W, height = H)
}



## =========================================================
## [NEW-DEG] WT vs MT Delta logFC Aggregation (Aggregate from existing DEG outputs)
##   - DEG Input Location (Per dataset / Per version / Per predictor):
##     proteomic_DEG/<dataset>/csn_gsea_results_TP53/<STRATUM>/DEG/<RAW|BatchAdj>/<PREDICTOR>/DEG_limma_cont.csv
##   - Output Location (Same hierarchy logic as GSEA Delta NES aggregation; but root changed to csn_deg_results_TP53):
##     proteomic_DEG/csn_deg_results_TP53/<RAW|BatchAdj>/<dataset>/DeltaWT_MT/<PREDICTOR>/DEG_limma_cont_deltaWTMT.csv
##   - Mapping: gene <=> pathway, logFC <=> NES
## =========================================================
`%||%` <- function(a, b) if (is.null(a)) b else a

deg_root_path <- function(ds, ver, stratum) {
  file.path(
    OUTPUT_PREFIX,
    ds, "csn_gsea_results_TP53", stratum, "DEG", ver
  )
}

tp53_delta_logfc_aggregate_DEG <- function(datasets_root = NULL,
                                           dataset_ids = NULL,
                                           versions = c("BatchAdj")) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
  if (is.null(dataset_ids)) {
    ## Prioritize using dataset_dirs_run already used in your main pipeline
    if (exists("dataset_dirs_run")) {
      dataset_ids <- names(dataset_dirs_run)
    } else if (!is.null(datasets_root)) {
      dataset_ids <- list.dirs(datasets_root, FALSE, FALSE)
    } else {
      stop("Please provide dataset_ids or datasets_root or existing dataset_dirs_run")
    }
  }

  for (ds in dataset_ids) {
    for (ver in versions) {
      ## Auto-detect predictor (subunit): Check subdirectories under WT folder
      wt_dir <- deg_root_path(ds, ver, "TP53_wild_type")
      mt_dir <- deg_root_path(ds, ver, "TP53_mutant")
      preds <- unique(c(list.dirs(wt_dir, FALSE, FALSE), list.dirs(mt_dir, FALSE, FALSE)))
      preds <- preds[nzchar(preds)]
      if (!length(preds)) next

      for (su in preds) {
        fp_wt <- file.path(wt_dir, su, "DEG_limma_cont.csv")
        fp_mt <- file.path(mt_dir, su, "DEG_limma_cont.csv")
        if (!file.exists(fp_wt) || !file.exists(fp_mt)) next

        wt <- tryCatch(data.table::fread(fp_wt, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
        mt <- tryCatch(data.table::fread(fp_mt, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
        if (is.null(wt) || is.null(mt) || !nrow(wt) || !nrow(mt)) next

        ## Normalize column names and map to GSEA conventions: gene->pathway; logFC->nes; Keep P.Value/adj.P.Val
        wt <- setNames(as.data.frame(wt), tolower(names(wt)))
        mt <- setNames(as.data.frame(mt), tolower(names(mt)))
        if (!("gene" %in% names(wt)) && "pathway" %in% names(wt)) names(wt)[names(wt) == "pathway"] <- "gene"
        if (!("gene" %in% names(mt)) && "pathway" %in% names(mt)) names(mt)[names(mt) == "pathway"] <- "gene"
        for (nm in c("logfc", "p.value", "adj.p.val")) if (!nm %in% names(wt)) wt[[nm]] <- NA_real_
        for (nm in c("logfc", "p.value", "adj.p.val")) if (!nm %in% names(mt)) mt[[nm]] <- NA_real_

        df <- dplyr::full_join(
          wt[, c("gene", "logfc", "adj.p.val", "p.value")],
          mt[, c("gene", "logfc", "adj.p.val", "p.value")],
          by = "gene", suffix = c("_WT", "_MT")
        )
        if (!nrow(df)) next

        df <- df |>
          dplyr::mutate(
            pathway = gene, # Align with GSEA aggregation interface
            NES_WT = as.numeric(logfc_WT), # logFC <=> NES
            NES_MT = as.numeric(logfc_MT),
            padj_WT = as.numeric(adj.p.val_WT),
            padj_MT = as.numeric(adj.p.val_MT),
            pval_WT = as.numeric(p.value_WT),
            pval_MT = as.numeric(p.value_MT),
            delta_NES = NES_MT - NES_WT,
            direction_consistency = dplyr::case_when(
              is.finite(NES_WT) & is.finite(NES_MT) & NES_MT > 0 & NES_WT > 0 ~ "BOTH_POS",
              is.finite(NES_WT) & is.finite(NES_MT) & NES_MT < 0 & NES_WT < 0 ~ "BOTH_NEG",
              is.finite(NES_WT) & is.finite(NES_MT) ~ "MIXED",
              TRUE ~ "NA"
            )
          ) |>
          dplyr::select(pathway, NES_WT, NES_MT, delta_NES, padj_WT, padj_MT, pval_WT, pval_MT, direction_consistency) |>
          dplyr::arrange(dplyr::desc(abs(delta_NES)), pathway)

        out_dir <- file.path(
          file.path(OUTPUT_PREFIX, "csn_deg_results_TP53"),
          ver, ds, "DeltaWT_MT", su
        )
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        data.table::fwrite(df, file.path(out_dir, "DEG_limma_cont_deltaWTMT.csv"))
        if (exists("log_msg", mode = "function")) {
          try(log_msg(
            "[Delta logFC-DEG] %s/%s/%s Output %s", ds, ver, su,
            file.path(out_dir, "DEG_limma_cont_deltaWTMT.csv")
          ), silent = TRUE)
        }
      }
      if (exists("log_msg", mode = "function")) {
        try(log_msg("[Delta logFC-DEG] %s/%s Done", ds, ver), silent = TRUE)
      }
    }
  }
  invisible(TRUE)
}


## =========================================================
## [NEW-DEG] meta-FDR (Stouffer -> BH) for DEG across datasets
##   - Read DEG_limma_cont.csv for each dataset, use P.Value as p, logFC sign determines Z direction
##   - Write to:
##     proteomic_DEG/csn_deg_pan_summary_TP53/meta_fdr/<STRATUM>/<RAW|BatchAdj>/<PREDICTOR>/DEG_meta_fdr.csv
##   - Columns follow GSEA aggregation interface: pathway(=gene), Z, padj_meta
## =========================================================
.meta_from_p_and_sign <- function(p, signv) {
  # Convert p-values to Z-scores with directional signs
  # Important: p must be valid (0 < p < 1), signv determines direction (+/-)
  p <- pmin(pmax(as.numeric(p), .Machine$double.eps), 1 - .Machine$double.eps)
  s <- sign(as.numeric(signv))

  # Critical fix: When sign is NA/NaN (unknown direction),
  # return NA instead of forcing positive to avoid systematic bias
  # Old behavior: s[!is.finite(s)] <- 1  # ← WRONG: forces NA → positive
  # New behavior: Keep NA as NA to exclude uncertain directions from meta-analysis

  z <- s * stats::qnorm(1 - p / 2) # two-sided p -> |Z|, then multiply by direction

  # Explicitly set non-finite signs to NA (defensive programming)
  z[!is.finite(s)] <- NA_real_

  z
}

meta_fdr_stouffer_DEG <- function(dataset_dirs,
                                  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
                                  versions = c("BatchAdj")) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
  if (is.null(dataset_dirs) || !length(dataset_dirs)) stop("dataset_dirs cannot be empty")

  ds_ids <- names(dataset_dirs)
  for (st in strata) {
    for (ver in versions) {
      ## Aggregate available predictor list from all datasets (Union) - Take subdirectory basename, remove empty strings and '.'
      preds_all <- unique(unlist(lapply(ds_ids, function(ds) {
        base <- deg_root_path(ds, ver, st)
        if (!dir.exists(base)) {
          return(character(0))
        }
        subdirs <- list.dirs(base, full.names = TRUE, recursive = FALSE)
        basename(subdirs)
      })))
      preds_all <- preds_all[nzchar(preds_all) & preds_all != "."]
      if (!length(preds_all)) next


      for (su in preds_all) {
        ## Read DEG table for this predictor per dataset, extract gene, logFC, P.Value
        lst <- lapply(ds_ids, function(ds) {
          fp <- file.path(deg_root_path(ds, ver, st), su, "DEG_limma_cont.csv")
          if (!file.exists(fp)) {
            return(NULL)
          }
          dt <- tryCatch(data.table::fread(fp, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
          if (is.null(dt) || !nrow(dt)) {
            return(NULL)
          }
          nm <- tolower(names(dt))
          names(dt) <- nm
          if (!("gene" %in% nm)) {
            if ("pathway" %in% nm) names(dt)[match("pathway", nm)] <- "gene" else return(NULL)
          }
          need <- c("gene", "logfc", "p.value")
          miss <- setdiff(need, names(dt))
          if (length(miss)) {
            return(NULL)
          }
          data.frame(
            gene = as.character(dt$gene),
            logfc = suppressWarnings(as.numeric(dt$logfc)),
            pval = suppressWarnings(as.numeric(dt$`p.value`)),
            stringsAsFactors = FALSE
          )
        })
        keep <- lst[!vapply(lst, is.null, logical(1))]
        if (!length(keep)) next

        ## Merge all datasets, calculate Stouffer Z and meta-BH for each gene
        ds_ids_keep <- ds_ids[!vapply(lst, is.null, logical(1))]
        all <- do.call(rbind, Map(function(df, id) {
          df$dataset <- id
          as.data.frame(df, stringsAsFactors = FALSE)
        }, keep, ds_ids_keep))
        if (is.null(all) || !ncol(all)) next
        if (!all(c("pval", "logfc") %in% names(all))) next
        all <- all[is.finite(all$pval) & is.finite(all$logfc), , drop = FALSE]
        if (!nrow(all)) next


        all$Z <- .meta_from_p_and_sign(all$pval, all$logfc)

        z_tab <- split(all, all$gene)
        z_sum <- lapply(names(z_tab), function(g) {
          z <- z_tab[[g]]$Z
          k <- sum(is.finite(z))
          Zmeta <- if (k) sum(z, na.rm = TRUE) / sqrt(k) else NA_real_
          pmeta <- if (is.finite(Zmeta)) 2 * pnorm(-abs(Zmeta)) else NA_real_
          c(gene = g, Z = Zmeta, p_meta = pmeta, n_ds = k)
        })
        Zdf <- data.frame(do.call(rbind, z_sum), stringsAsFactors = FALSE)
        Zdf$Z <- suppressWarnings(as.numeric(Zdf$Z))
        Zdf$p_meta <- suppressWarnings(as.numeric(Zdf$p_meta))
        Zdf$padj_meta <- p.adjust(Zdf$p_meta, method = "BH")
        Zdf$pathway <- Zdf$gene # Consistent with GSEA summary interface (pathway column exists)
        Zdf <- Zdf[, c("pathway", "Z", "p_meta", "padj_meta", "n_ds")]

        out_dir <- file.path(
          file.path(OUTPUT_PREFIX, "csn_deg_pan_summary_TP53/meta_fdr"),
          st, ver, su
        )
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        data.table::fwrite(Zdf, file.path(out_dir, "DEG_meta_fdr.csv"))
        if (exists("log_msg", mode = "function")) {
          try(log_msg(
            "[meta-DEG] %s | %s | %s: Write %s",
            st, ver, su, file.path(out_dir, "DEG_meta_fdr.csv")
          ), silent = TRUE)
        }
      }
    }
  }
  invisible(TRUE)
}

## ---- Summary Tool (DEG Version): Reuse your GSEA wide table + count output logic, but change stat_tag and filename to DEG ----
.read_meta_fdr_table_DEG <- function(meta_root, stratum, version, subunit,
                                     group_name = "GENE", stat_tag = "DEG_meta_fdr") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
  grp <- safe_fs_name(group_name)
  fp <- file.path(meta_root, stratum, version, subunit, grp, paste0(stat_tag, ".csv"))
  ## If no group subdirectory, also tolerate: .../<subunit>/DEG_meta_fdr.csv
  if (!file.exists(fp)) fp <- file.path(meta_root, stratum, version, subunit, paste0(stat_tag, ".csv"))
  if (!file.exists(fp)) {
    return(NULL)
  }
  dt <- tryCatch(data.table::fread(fp, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
  if (is.null(dt)) {
    return(NULL)
  }
  nms <- tolower(names(dt))
  names(dt) <- nms # Lowercase first, then rename to 'pathway'/'Z' as needed
  if (!("pathway" %in% tolower(names(dt)))) {
    return(NULL)
  }
  if (!("z" %in% tolower(names(dt)))) {
    return(NULL)
  }
  if (!("padj_meta" %in% tolower(names(dt)))) {
    return(NULL)
  }
  # Unify column names (Keep GSEA interface)
  if (!"pathway" %in% names(dt) && "pathway" %in% nms) names(dt)[match("pathway", nms)] <- "pathway"
  if (!"Z" %in% names(dt) && "z" %in% nms) names(dt)[match("z", nms)] <- "Z"
  if (!"padj_meta" %in% names(dt) && "padj_meta" %in% nms) names(dt)[match("padj_meta", nms)] <- "padj_meta"
  out <- as.data.frame(dt[, c("pathway", "Z", "padj_meta")])
  names(out)[names(out) == "Z"] <- paste0("Z_", subunit)
  names(out)[names(out) == "padj_meta"] <- paste0("padj_meta_", subunit)
  out
}

.merge_subunit_tables_meta_DEG <- function(tbl_list) {
  keep <- tbl_list[!vapply(tbl_list, is.null, logical(1))]
  if (!length(keep)) {
    return(NULL)
  }
  out <- Reduce(function(x, y) dplyr::full_join(x, y, by = "pathway"), keep)
  num_cols <- setdiff(names(out), "pathway")
  out[num_cols] <- lapply(out[num_cols], function(z) suppressWarnings(as.numeric(z)))
  out
}

.add_sig_counts_meta_DEG <- function(df, alphas = c(0.05, 0.25)) {
  if (is.null(df) || !nrow(df)) {
    return(df)
  }
  padj_cols <- grep("^padj_meta_", names(df), value = TRUE)
  subs <- sub("^padj_meta_", "", padj_cols)
  z_cols <- paste0("Z_", subs)
  keep_idx <- z_cols %in% names(df)
  padj_cols <- padj_cols[keep_idx]
  z_cols <- z_cols[keep_idx]
  if (!length(padj_cols)) {
    return(df)
  }
  for (a in alphas) {
    a_tag <- gsub("\\.", "_", sprintf("%.2f", a))
    sig_mat <- sapply(padj_cols, function(p) {
      pv <- df[[p]]
      as.integer(is.finite(pv) & pv < a)
    })
    if (is.null(dim(sig_mat))) sig_mat <- matrix(sig_mat, ncol = 1)
    df[[sprintf("sig_n_padj_meta_%s", a_tag)]] <- rowSums(sig_mat, na.rm = TRUE)

    pos_mat <- mapply(function(p, zc) {
      pv <- df[[p]]
      zv <- df[[zc]]
      as.integer(is.finite(pv) & pv < a & is.finite(zv) & zv > 0)
    }, padj_cols, z_cols, SIMPLIFY = TRUE)
    if (is.null(dim(pos_mat))) pos_mat <- matrix(pos_mat, ncol = 1)
    df[[sprintf("pos_n_padj_meta_%s", a_tag)]] <- rowSums(pos_mat, na.rm = TRUE)

    neg_mat <- mapply(function(p, zc) {
      pv <- df[[p]]
      zv <- df[[zc]]
      as.integer(is.finite(pv) & pv < a & is.finite(zv) & zv < 0)
    }, padj_cols, z_cols, SIMPLIFY = TRUE)
    if (is.null(dim(neg_mat))) neg_mat <- matrix(neg_mat, ncol = 1)
    df[[sprintf("neg_n_padj_meta_%s", a_tag)]] <- rowSums(neg_mat, na.rm = TRUE)
  }
  df
}

.write_summary_outputs_meta_csv_DEG <- function(df, out_dir, tag = "DEG_meta_fdr") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install openxlsx first")
  if (is.null(df) || !nrow(df)) {
    return(invisible(NULL))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(out_dir, paste0("Summary_GENE_", tag))

  ## Sort: Significant count (0.05) -> Positive/Negative count -> pathway
  ord_keys <- c("sig_n_padj_meta_0_05", "pos_n_padj_meta_0_05", "neg_n_padj_meta_0_05")
  ord_keys <- intersect(ord_keys, names(df))
  if (length(ord_keys)) {
    df <- df |>
      dplyr::arrange(
        dplyr::desc(.data[[ord_keys[1]]]),
        dplyr::desc(ifelse(length(ord_keys) > 1, .data[[ord_keys[2]]], 0)),
        dplyr::desc(ifelse(length(ord_keys) > 2, .data[[ord_keys[3]]], 0)),
        .data[["pathway"]]
      )
  } else {
    df <- df |> dplyr::arrange(.data[["pathway"]])
  }

  data.table::fwrite(df, paste0(base, "_ALL.csv"))
  if ("sig_n_padj_meta_0_05" %in% names(df)) {
    data.table::fwrite(dplyr::filter(df, .data[["sig_n_padj_meta_0_05"]] > 0), paste0(base, "_padjLT0.05.csv"))
  }
  if ("sig_n_padj_meta_0_25" %in% names(df)) {
    data.table::fwrite(dplyr::filter(df, .data[["sig_n_padj_meta_0_25"]] > 0), paste0(base, "_padjLT0.25.csv"))
  }

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "ALL")
  openxlsx::writeData(wb, "ALL", df)
  openxlsx::freezePane(wb, "ALL", firstRow = TRUE)
  openxlsx::setColWidths(wb, "ALL", cols = 1:ncol(df), widths = "auto")

  if ("sig_n_padj_meta_0_05" %in% names(df)) {
    openxlsx::addWorksheet(wb, "padjLT0.05")
    openxlsx::writeData(wb, "padjLT0.05", dplyr::filter(df, .data[["sig_n_padj_meta_0_05"]] > 0))
    openxlsx::freezePane(wb, "padjLT0.05", firstRow = TRUE)
    openxlsx::setColWidths(wb, "padjLT0.05", cols = 1:ncol(df), widths = "auto")
  }
  if ("sig_n_padj_meta_0_25" %in% names(df)) {
    openxlsx::addWorksheet(wb, "padjLT0.25")
    openxlsx::writeData(wb, "padjLT0.25", dplyr::filter(df, .data[["sig_n_padj_meta_0_25"]] > 0))
    openxlsx::freezePane(wb, "padjLT0.25", firstRow = TRUE)
    openxlsx::setColWidths(wb, "padjLT0.25", cols = 1:ncol(df), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, paste0(base, ".xlsx"), overwrite = TRUE)
  invisible(NULL)
}

summarize_meta_fdr_across_subunits_DEG <- function(
  meta_root = file.path(OUTPUT_PREFIX, "csn_deg_pan_summary_TP53/meta_fdr"),
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  versions = c("RAW", "BatchAdj")
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
  for (st in strata) {
    for (ver in versions) {
      base_dir <- file.path(meta_root, st, ver)
      if (!dir.exists(base_dir)) next
      subs <- basename(list.dirs(base_dir, full.names = TRUE, recursive = FALSE))
      subs <- subs[nzchar(subs) & subs != "."]
      if (!length(subs)) next

      ## Here group is fixed to "GENE" (DEG has no geneset grouping)
      for (grp in "GENE") {
        if (exists("log_msg", mode = "function")) try(log_msg("== [meta-summary-DEG] stratum=%s | version=%s ==", st, ver), silent = TRUE)
        lst <- setNames(vector("list", length(subs)), subs)
        for (su in subs) lst[[su]] <- .read_meta_fdr_table_DEG(meta_root, st, ver, su, grp, "DEG_meta_fdr")
        wide <- .merge_subunit_tables_meta_DEG(lst)
        wide <- .add_sig_counts_meta_DEG(wide, alphas = c(0.05, 0.25))
        out_dir <- file.path(meta_root, "summary", st, ver, "GENE", "DEG_meta_fdr")
        .write_summary_outputs_meta_csv_DEG(wide, out_dir, "DEG_meta_fdr")
      }
    }
  }
  invisible(TRUE)
}

## ---- One-click DEG Post-processing (Read existing output; do not rerun any models) ----
posthoc_summary_meta_fdr_DEG <- function() {
  ## 1) First do Delta logFC (WT vs MT)
  try(tp53_delta_logfc_aggregate_DEG(
    dataset_ids = names(dataset_dirs_run),
    versions = c("RAW", "BatchAdj")
  ), silent = TRUE)
  ## 2) Do DEG version meta-FDR & Summary
  meta_fdr_stouffer_DEG(
    dataset_dirs = dataset_dirs_run,
    strata       = c("ALL", "TP53_mutant", "TP53_wild_type"),
    versions     = c("BatchAdj")
  )
  summarize_meta_fdr_across_subunits_DEG(
    meta_root = file.path(OUTPUT_PREFIX, "csn_deg_pan_summary_TP53/meta_fdr"),
    strata    = c("ALL", "TP53_mutant", "TP53_wild_type"),
    versions  = c("BatchAdj")
  )
  invisible(TRUE)
}

## ---- Run directly once (Can be re-executed; overwrite output only) ----
posthoc_summary_meta_fdr_DEG()



## =========================================================
## [NEW-DEG-HEATMAP-GENESET] "Gene x Predictors" Pan-Heatmap based on DEG meta-FDR
##  - Source CSV: .../GENE/DEG_meta_fdr/Summary_GENE_DEG_meta_fdr_ALL.csv
##  - y-axis: Gene (= pathway column of CSV; actually gene)
##  - x-axis: Predictors (Z_* columns; order follows .pred_order_meta)
##  - ALL: Sort by Z_CSN_SCORE descending, and can take Top/Bottom (Default 25/25)
##  - MUT / WT: Each produces two versions (Self-sorted; Follow ALL list + order)
##  - Color: Customizable up/down; default consistent with GSEA Blue-White-Red
##  - Output: Default .tiff (600 dpi LZW); can add png/jpg/pdf
## =========================================================

# -- Get gene list corresponding to geneset names specified by user (Collect all groups from genesets_by_group)
.get_genes_from_gs_names <- function(gs_names) {
  stopifnot(exists("genesets_by_group", inherits = TRUE))
  gs_names <- unique(as.character(gs_names))
  out <- list()
  for (grp in names(genesets_by_group)) {
    lst <- genesets_by_group[[grp]]
    hits <- intersect(names(lst), gs_names)
    if (length(hits)) {
      for (nm in hits) out[[nm]] <- sort(unique(lst[[nm]]))
    }
  }
  if (!length(out)) {
    warning(
      "[deg-heatmap] Specified gs_name not found in genesets_by_group: ",
      paste(gs_names, collapse = ", ")
    )
  }
  out
}

# -- Find three layers of Summary_GENE_DEG_meta_fdr_ALL.csv (Restrict RAW / BatchAdj)
.deg_meta_find_csvs <- function(root, versions = c("BatchAdj")) {
  files <- list.files(root,
    pattern = "^Summary_GENE_DEG_meta_fdr_ALL\\.csv$",
    recursive = TRUE, full.names = TRUE
  )
  if (!length(files)) {
    return(character(0))
  }
  vv <- tolower(versions)
  keep_ver <- Reduce(`|`, lapply(vv, function(v) {
    grepl(paste0("(/|\\\\)", v, "(/|\\\\)"), tolower(files))
  }))
  keep_str <- grepl("(/|\\\\)(all|tp53_mutant|tp53_wild_type)(/|\\\\)", tolower(files))
  files <- files[keep_ver & keep_str]
  # Exclude size=0
  if (length(files)) {
    finfo <- suppressWarnings(file.info(files))
    files <- files[is.finite(finfo$size) & finfo$size > 0 & !finfo$isdir]
  }
  unique(files)
}

# -- Determine y-axis "Selection List + Order" based on ALL layer (Sort by Z_CSN_SCORE descending; then cut top/bottom)
.deg_meta_y_order_from_ALL <- function(csv_all, gene_set, top_n = 25L, bottom_n = 25L) {
  df_long <- .meta_read_summary_long(csv_all) # Reuse existing long table tool
  present <- intersect(.pred_order_meta, unique(df_long$predictor))
  if (!length(present)) stop("[deg-heatmap] File missing expected predictors: ", csv_all)
  df_long <- dplyr::filter(
    df_long, .data$predictor %in% present,
    .data$pathway %in% gene_set
  )
  # Sort by Z_CSN_SCORE; if not present, use column mean Z
  if ("Z_CSN_SCORE" %in% present) {
    z_tbl <- df_long %>%
      dplyr::filter(.data$predictor == "Z_CSN_SCORE") %>%
      dplyr::select(.data$pathway, .data$Z) %>%
      dplyr::distinct()
    keep_up <- z_tbl %>%
      dplyr::arrange(dplyr::desc(.data$Z)) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull(.data$pathway)
    keep_dn <- z_tbl %>%
      dplyr::arrange(.data$Z) %>%
      dplyr::slice_head(n = bottom_n) %>%
      dplyr::pull(.data$pathway)
    keep <- unique(c(keep_up, keep_dn))
    ord <- z_tbl %>%
      dplyr::arrange(dplyr::desc(.data$Z)) %>%
      dplyr::pull(.data$pathway)
    ord <- ord[ord %in% keep]
  } else {
    avg <- df_long %>%
      dplyr::group_by(.data$pathway) %>%
      dplyr::summarise(m = mean(.data$Z, na.rm = TRUE), .groups = "drop")
    keep_up <- avg %>%
      dplyr::arrange(dplyr::desc(.data$m)) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull(.data$pathway)
    keep_dn <- avg %>%
      dplyr::arrange(.data$m) %>%
      dplyr::slice_head(n = bottom_n) %>%
      dplyr::pull(.data$pathway)
    keep <- unique(c(keep_up, keep_dn))
    ord <- avg %>%
      dplyr::arrange(dplyr::desc(.data$m)) %>%
      dplyr::pull(.data$pathway)
    ord <- ord[ord %in% keep]
  }
  list(keep = keep, order = ord)
}

# -- Build DEG heatmap with "Gene List (Can include top/bottom filter)" (Style follows .meta_make_heatmap_plot)
.deg_meta_make_heatmap_plot_genes <- function(csv_file,
                                              gene_keep = NULL, y_order = NULL,
                                              palette = c(low = "#053061", mid = "#FFFFFF", high = "#67001F")) {
  df_long <- .meta_read_summary_long(csv_file)
  present <- intersect(.pred_order_meta, unique(df_long$predictor))
  if (!length(present)) stop("[deg-heatmap] No available predictors: ", csv_file)
  df_long <- dplyr::filter(df_long, .data$predictor %in% present)
  if (!is.null(gene_keep) && length(gene_keep)) {
    df_long <- dplyr::filter(df_long, .data$pathway %in% unique(gene_keep))
  }
  # y order: If given, follow strictly; otherwise sort by Z_CSN_SCORE (or mean Z)
  if (!is.null(y_order) && length(y_order)) {
    df_long$pathway <- factor(df_long$pathway, levels = rev(y_order))
  } else {
    if ("Z_CSN_SCORE" %in% present) {
      ord <- df_long %>%
        dplyr::filter(.data$predictor == "Z_CSN_SCORE") %>%
        dplyr::select(.data$pathway, .data$Z) %>%
        dplyr::distinct() %>%
        dplyr::arrange(dplyr::desc(.data$Z)) %>%
        dplyr::pull(.data$pathway)
    } else {
      ord <- df_long %>%
        dplyr::group_by(.data$pathway) %>%
        dplyr::summarise(m = mean(.data$Z, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(.data$m)) %>%
        dplyr::pull(.data$pathway)
    }
    df_long$pathway <- factor(df_long$pathway, levels = rev(ord))
  }
  # x-axis position and gap rules (Leave white space between column groups, consistent with proteomic GSEA meta heatmap)
  # Rules:
  #   1. Leave a small gap between Z_CSN_SCORE and Z_GPS1
  #   2. Leave another small gap before Z_GPS1_adj_CSN (interaction models)
  # These two gaps actually push the x-coordinate to the right, creating visual vertical white space.

  gap <- 0.4
  needs_gap1 <- all(c("Z_CSN_SCORE", "Z_GPS1") %in% present)
  needs_gap2 <- "Z_GPS1_adj_CSN" %in% present # Updated: use new naming

  pos_map <- list()
  pos <- 0
  for (pname in present) {
    if (pname == "Z_GPS1" && needs_gap1) {
      pos <- pos + gap
    }
    if (pname == "Z_GPS1_adj_CSN" && needs_gap2) { # Updated: use new naming
      pos <- pos + gap
    }
    pos <- pos + 1
    pos_map[[pname]] <- pos
  }
  pos_map <- unlist(pos_map)

  df_long$xpos <- unname(pos_map[df_long$predictor])

  L <- ceiling(max(abs(df_long$Z[is.finite(df_long$Z)]), na.rm = TRUE))
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = pathway, fill = Z)) +
    ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
    ggplot2::geom_point(
      data = dplyr::filter(df_long, is.finite(.data$padj_meta), .data$padj_meta < 0.05),
      ggplot2::aes(x = xpos, y = pathway),
      shape = 16, size = 1.6, color = "black", inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = palette[["low"]], mid = palette[["mid"]], high = palette[["high"]],
      limits = c(-L, L), midpoint = 0, oob = scales::squish, name = "Z"
    ) +
    ggplot2::scale_x_continuous(
      breaks = unname(pos_map[present]),
      labels = present, expand = ggplot2::expansion(mult = c(0.01, 0.01)), position = "top"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(
        angle = 45, hjust = 0, vjust = 0, size = 9,
        margin = ggplot2::margin(b = 8)
      ),
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = "right",
      plot.margin = ggplot2::margin(6, 12, 6, 6),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::coord_cartesian(clip = "off")
  n_path <- nlevels(df_long$pathway)
  W <- max(8, length(present) * 0.45)
  H <- max(6, n_path * 0.22)
  list(plot = p, width = W, height = H)
}

# -- Multi-format output (Default tiff; can add png/jpg/pdf)
.save_multi_formats <- function(p, base, width, height, formats = c("tiff"), dpi = 600) {
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  fmts <- unique(tolower(formats))
  for (f in fmts) {
    target <- paste0(base, ".", f)
    if (f %in% c("tiff", "tif")) {
      ggplot2::ggsave(target, p,
        width = width, height = height, units = "in",
        dpi = dpi, bg = "white", compression = "lzw"
      )
    } else {
      ggplot2::ggsave(target, p,
        width = width, height = height, units = "in",
        dpi = dpi, bg = "white"
      )
    }
  }
}

# -- Create output filename (Save nearby to csv folder)
.deg_meta_out_base <- function(csv_file, gs_name, stratum, version, suffix = "") {
  safe <- function(s) gsub('[<>:"/\\\\|?*\\s]+', "_", s)
  paste0(
    file.path(
      dirname(csv_file),
      paste0("heatmap_DEG_", safe(stratum), "_", safe(version), "_", safe(gs_name))
    ),
    suffix
  )
}


# [NEW] Export complete data for heatmap (Wide table: Row=gene; Col=Z_* and padj_meta_* for each predictor)
.deg_meta_write_data_csv <- function(csv_file, gene_set, out_csv, predictor_order = .pred_order_meta) {
  df <- suppressMessages(.read_csv_safe(csv_file))

  # Auto-detect pathway column and unify as pathway
  path_candidates <- c("pathway", "Pathway", "term", "Term", "gs_name", "NAME", "set", "Set")
  path_col <- intersect(path_candidates, names(df))[1]
  if (is.na(path_col)) stop("[deg-heatmap] Data export failed: Pathway column not found in: ", csv_file)
  df <- dplyr::rename(df, pathway = dplyr::all_of(path_col))

  # Keep only genes intersecting with geneset ("Intersection of all geneset genes", not just heatmap y-axis)
  df <- dplyr::filter(df, .data$pathway %in% unique(gene_set))
  if (!nrow(df)) stop("[deg-heatmap] Data export failed: Geneset has no intersection with file: ", csv_file)

  # Determine order based on actual Z_* columns; fill missing padj_meta_* columns with NA
  z_present <- intersect(predictor_order, grep("^Z_", names(df), value = TRUE))
  keys <- sub("^Z_", "", z_present)
  padj_need <- paste0("padj_meta_", keys)
  missing_p <- setdiff(padj_need, names(df))
  for (nm in missing_p) df[[nm]] <- NA_real_

  # Export columns: pathway + all Z_* (in order) + corresponding padj_meta_*
  out <- dplyr::select(df, dplyr::all_of(c("pathway", z_present, padj_need)))
  readr::write_csv(out, out_csv)
  message(sprintf("[deg-heatmap] Data exported: %s (rows=%d, cols=%d)", out_csv, nrow(out), ncol(out)))
}


# -- Main Program: Plot ALL/MUT/WT (including ALL-ordered) based on specified geneset name list and pass version

# -- Auto-detect DEG pan-summary root directory
.deg_meta_autoroot <- function(root_hint) {
  has_target <- function(r) {
    dir.exists(r) && length(list.files(
      r,
      pattern = "^Summary_GENE_DEG_meta_fdr_ALL\\.csv$", recursive = TRUE
    )) > 0
  }
  if (has_target(root_hint)) {
    return(normalizePath(root_hint, winslash = "/", mustWork = FALSE))
  }

  cand <- c(file.path(OUTPUT_PREFIX, "csn_deg_pan_summary_TP53/meta_fdr/summary"), "csn_deg_pan_summary_TP53/meta_fdr/summary")

  dirs <- try(list.dirs(".", full.names = TRUE, recursive = TRUE), silent = TRUE)
  if (!inherits(dirs, "try-error") && length(dirs)) {
    hit <- dirs[grepl("(/|\\\\)csn_deg_pan_summary_TP53(/|\\\\)meta_fdr(/|\\\\)summary$",
      dirs,
      ignore.case = TRUE
    )]
    cand <- c(cand, hit)
  }
  cand <- unique(cand)
  cand <- cand[dir.exists(cand)]

  if (length(cand)) {
    scores <- vapply(
      cand, function(r) {
        length(list.files(r, pattern = "^Summary_GENE_DEG_meta_fdr_ALL\\.csv$", recursive = TRUE))
      },
      integer(1)
    )
    if (any(scores > 0L)) {
      chosen <- cand[which.max(scores)]
      message(sprintf(
        "[deg-heatmap] Auto-selected root = %s (Detected %d target files)",
        normalizePath(chosen, winslash = "/", mustWork = FALSE),
        max(scores)
      ))
      return(normalizePath(chosen, winslash = "/", mustWork = FALSE))
    }
  }
  stop(sprintf(
    "[deg-heatmap] DEG pan-summary root directory not found. Tried root='%s' and candidates: %s",
    root_hint, paste(cand, collapse = "; ")
  ))
}


run_deg_geneset_heatmaps <- function(
  gs_names,
  root = file.path(OUTPUT_PREFIX, "csn_deg_pan_summary_TP53/meta_fdr/summary"),
  versions = getOption("csn.run_passes", c("BatchAdj")),
  top_n = get0("PAN_HEATMAP_TOP_N", ifnotfound = 25L),
  bottom_n = get0("PAN_HEATMAP_BOTTOM_N", ifnotfound = 25L),
  up_color = "#67001F", mid_color = "#FFFFFF", down_color = "#053061",
  formats = c("tiff") # <- Default output tiff only
) {
  gs_map <- .get_genes_from_gs_names(gs_names)
  if (!length(gs_map)) {
    stop("[deg-heatmap] No specified geneset found: ", paste(gs_names, collapse = ", "))
  }
  root <- .deg_meta_autoroot(root)
  csvs <- .deg_meta_find_csvs(root, versions = versions)
  if (!length(csvs)) {
    stop(sprintf(
      "[deg-heatmap] No Summary_GENE_DEG_meta_fdr_ALL.csv found in root='%s' (including auto-detect); please check: %s/<STRATUM>/%s/GENE/DEG_meta_fdr/",
      root, "csn_deg_pan_summary_TP53/meta_fdr/summary", paste(versions, collapse = "|")
    ))
  }
  tolc <- function(x) tolower(normalizePath(x, winslash = "/", mustWork = FALSE))
  for (gs in names(gs_map)) {
    genes <- unique(gs_map[[gs]])
    message(sprintf("[deg-heatmap] geneset = %s (n=%d)", gs, length(genes)))
    for (ver in versions) {
      # Grab CSV for each of the three layers
      pick <- function(st) {
        # 1) Combine standard paths directly (Most common; best performance)
        cand1 <- file.path(root, st, ver, "GENE", "DEG_meta_fdr", "Summary_GENE_DEG_meta_fdr_ALL.csv")
        if (file.exists(cand1)) {
          return(normalizePath(cand1, winslash = "/", mustWork = FALSE))
        }
        # 2) Tolerate configuration without GENE subdirectory
        cand2 <- file.path(root, st, ver, "DEG_meta_fdr", "Summary_GENE_DEG_meta_fdr_ALL.csv")
        if (file.exists(cand2)) {
          return(normalizePath(cand2, winslash = "/", mustWork = FALSE))
        }
        # 3) Finally search from scanned csvs with looser matching conditions
        patt <- paste0(
          "(/|\\\\)", tolower(st), "(/|\\\\).*?(/|\\\\)", tolower(ver),
          "(/|\\\\).*?(/|\\\\)gene(/|\\\\)deg_meta_fdr(/|\\\\)summary_gene_deg_meta_fdr_all\\.csv$"
        )
        hits <- csvs[grepl(patt, tolower(normalizePath(csvs, winslash = "/", mustWork = FALSE)))]
        if (length(hits) >= 1) hits[1] else NA_character_
      }
      csv_all <- pick("ALL")
      csv_mut <- pick("TP53_mutant")
      csv_wt <- pick("TP53_wild_type")
      if (is.na(csv_all)) {
        stop(sprintf(
          "[deg-heatmap] ALL layer not found (%s, %s). root='%s'; please check if Summary_GENE_DEG_meta_fdr_ALL.csv exists and is not empty in that layer.",
          gs, ver, normalizePath(root, winslash = "/", mustWork = FALSE)
        ))
      }

      pal <- c(low = down_color, mid = mid_color, high = up_color)


      ## --- A) ALL: Select and sort top/bottom by ALL's own Z_CSN_SCORE ---
      yinfo <- .deg_meta_y_order_from_ALL(csv_all,
        gene_set = genes,
        top_n = top_n, bottom_n = bottom_n
      )
      h_all <- .deg_meta_make_heatmap_plot_genes(csv_all,
        gene_keep = yinfo$keep, y_order = yinfo$order, palette = pal
      )
      base_all <- .deg_meta_out_base(csv_all, gs, "ALL", ver, suffix = "")
      .save_multi_formats(h_all$plot, base_all, h_all$width, h_all$height, formats = formats)
      message(sprintf("[deg-heatmap] Exported: %s.{%s}", base_all, paste(unique(tolower(formats)), collapse = ",")))
      # [NEW] Export complete data corresponding to ALL heatmap (All genes x All predictors)
      .deg_meta_write_data_csv(
        csv_file = csv_all, gene_set = genes,
        out_csv = paste0(base_all, "_data.csv"),
        predictor_order = .pred_order_meta
      )

      ## --- B) TP53_mutant: Self-sorted version + ALL-ordered version ---
      if (!is.na(csv_mut)) {
        # Self-sorted (Intersection of genes only)
        h_mut_self <- .deg_meta_make_heatmap_plot_genes(csv_mut,
          gene_keep = yinfo$keep, y_order = NULL, palette = pal
        )
        base_mut_self <- .deg_meta_out_base(csv_mut, gs, "TP53_mutant", ver, suffix = "")
        .save_multi_formats(h_mut_self$plot, base_mut_self, h_mut_self$width, h_mut_self$height,
          formats = formats
        )
        message(sprintf("[deg-heatmap] Exported: %s.{%s}", base_mut_self, paste(unique(tolower(formats)), collapse = ",")))
        # [NEW] Export complete data corresponding to TP53_mutant (Self-sorted) heatmap
        .deg_meta_write_data_csv(
          csv_file = csv_mut, gene_set = genes,
          out_csv = paste0(base_mut_self, "_data.csv"),
          predictor_order = .pred_order_meta
        )

        # ALL-ordered (Follow ALL selection list + order)
        h_mut_allord <- .deg_meta_make_heatmap_plot_genes(csv_mut,
          gene_keep = yinfo$keep, y_order = yinfo$order, palette = pal
        )
        base_mut_allord <- .deg_meta_out_base(csv_mut, gs, "TP53_mutant", ver, suffix = "_ALLordered")
        .save_multi_formats(h_mut_allord$plot, base_mut_allord, h_mut_allord$width, h_mut_allord$height,
          formats = formats
        )
        message(sprintf("[deg-heatmap] Exported: %s.{%s}", base_mut_allord, paste(unique(tolower(formats)), collapse = ",")))
        # [NEW] Export complete data corresponding to TP53_mutant (ALL-ordered) heatmap
        .deg_meta_write_data_csv(
          csv_file = csv_mut, gene_set = genes,
          out_csv = paste0(base_mut_allord, "_data.csv"),
          predictor_order = .pred_order_meta
        )
      }

      ## --- C) TP53_wild_type: Self-sorted version + ALL-ordered version ---
      if (!is.na(csv_wt)) {
        h_wt_self <- .deg_meta_make_heatmap_plot_genes(csv_wt,
          gene_keep = yinfo$keep, y_order = NULL, palette = pal
        )
        base_wt_self <- .deg_meta_out_base(csv_wt, gs, "TP53_wild_type", ver, suffix = "")
        .save_multi_formats(h_wt_self$plot, base_wt_self, h_wt_self$width, h_wt_self$height,
          formats = formats
        )
        message(sprintf("[deg-heatmap] Exported: %s.{%s}", base_wt_self, paste(unique(tolower(formats)), collapse = ",")))
        # [NEW] Export complete data corresponding to TP53_wild_type (Self-sorted) heatmap
        .deg_meta_write_data_csv(
          csv_file = csv_wt, gene_set = genes,
          out_csv = paste0(base_wt_self, "_data.csv"),
          predictor_order = .pred_order_meta
        )
        h_wt_allord <- .deg_meta_make_heatmap_plot_genes(csv_wt,
          gene_keep = yinfo$keep, y_order = yinfo$order, palette = pal
        )
        base_wt_allord <- .deg_meta_out_base(csv_wt, gs, "TP53_wild_type", ver, suffix = "_ALLordered")
        .save_multi_formats(h_wt_allord$plot, base_wt_allord, h_wt_allord$width, h_wt_allord$height,
          formats = formats
        )
        message(sprintf("[deg-heatmap] Exported: %s.{%s}", base_wt_allord, paste(unique(tolower(formats)), collapse = ",")))
        # [NEW] Export complete data corresponding to TP53_wild_type (ALL-ordered) heatmap
        .deg_meta_write_data_csv(
          csv_file = csv_wt, gene_set = genes,
          out_csv = paste0(base_wt_allord, "_data.csv"),
          predictor_order = .pred_order_meta
        )
      }
    }
  }
  invisible(TRUE)
}

## ---- [USER CONFIG | heatmap top/bottom limits] ----
## Non-H collection heatmap filtering threshold (Single dataset uses NES, PAN uses Z; both use CSN_SCORE as filtering basis)
## - Single dataset: Keep only top N and bottom M NES with padj<0.05 (by CSN_SCORE)
## - PAN (meta-FDR): Keep only top N and bottom M Z with padj_meta<0.05 (by CSN_SCORE)
DATASET_HEATMAP_TOP_N <- 25 # Default top 25 (NES)
DATASET_HEATMAP_BOTTOM_N <- 25 # Default bottom 25 (NES)
PAN_HEATMAP_TOP_N <- 10 # Default top 25 (Z)
PAN_HEATMAP_BOTTOM_N <- 5 # Default bottom 25 (Z)
## ---- [END USER CONFIG | heatmap top/bottom limits] ----


## need: .meta_read_summary_long, .read_csv_safe, .pred_order_all, .pred_order_all
## --------- [DEMO] Call Example (Uncomment to execute directly) ---------
run_deg_geneset_heatmaps(
  gs_names = c(  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
     "HALLMARK_HYPOXIA",
     "HALLMARK_GLYCOLYSIS", "HALLMARK_FATTY_ACID_METABOLISM"
     #"HALLMARK_INTERFERON_GAMMA_RESPONSE",
     #"HALLMARK_INTERFERON_ALPHA_RESPONSE"
  ),
  root = file.path(OUTPUT_PREFIX, "csn_deg_pan_summary_TP53/meta_fdr/summary"),
  versions = getOption("csn.run_passes", c("BatchAdj")), # or c("RAW","BatchAdj")
  top_n = get0("PAN_HEATMAP_TOP_N", ifnotfound = 25L),
  bottom_n = get0("PAN_HEATMAP_BOTTOM_N", ifnotfound = 25L),
  up_color = "#67001F", mid_color = "#FFFFFF", down_color = "#053061",
  formats = c("tiff") # Can be changed to c("tiff","png","jpg","pdf")
)
