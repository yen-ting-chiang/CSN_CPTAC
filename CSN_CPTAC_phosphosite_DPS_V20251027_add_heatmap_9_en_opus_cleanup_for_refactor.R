## =========================================================
## CSN subunits → Differential Phosphosite (DPS) Analysis
## (CPTAC, with TP53 stratification)
## Output to strata = c("ALL","TP53_mutant","TP53_wild_type")
## =========================================================

setwd("C:/Users/danny/Documents/R_project/CSN_CPTAC") ## YTC laptop
## setwd("C:/Users/cmuh/Documents/YenTing_document/CSN_CPTAC") ## lab computer
getwd()


## ===== Add yaml package in the package block (1) =====
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

  library(openxlsx)
  library(ComplexHeatmap)
  library(cowplot)
  library(matrixStats)
  library(yaml) # <— Added: for writing run_manifest.yml
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
        if (length(lv) <= 1) { # Single-level factor -> drop to avoid contrasts error
          keep[cn] <- FALSE
          if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
        } else {
          df[[cn]] <- v
        }
      } else {
        df[[cn]] <- suppressWarnings(as.numeric(v)) # Keep numeric as numeric
        # To also filter out constant numeric columns, uncomment the following:
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


## ==== Force single-thread execution, disable all parallel backends (across common frameworks) ====
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
    # Try to close residual cluster (if you stored cluster in .GlobalEnv, call stopCluster yourself)
    # Here we don't forcibly stop unknown objects, just use DoSEQ.
  }

  # future
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
  Sys.setenv("R_FUTURE_FORK_ENABLE" = "FALSE")
}
.force_serial_execution()

## Create run_info directory (before writing any run_info/* files)
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)


## Attach DPS analysis description
cat(
  paste(
    "DPS Analysis (Differential Phosphosite):",
    " - Uses limma for differential phosphosite analysis.",
    " - Meta-analysis across datasets using Stouffer's z to obtain meta p-values,",
    "   followed by Benjamini-Hochberg to report meta-level FDR.",
    sep = "\n"
  ),
  file = file.path("run_info", "analysis_notes.txt")
)


## ===== Parameters =====
set.seed(1234)

csn_subunits <- c("GPS1", "COPS2", "COPS3", "COPS4", "COPS5", "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9")

min_frac_complete <- 0.75 # At least 75% non-NA per gene
hi_lo_quantile <- 0.25 # Top and bottom 25% each


min_per_group <- 8 # Minimum samples for High/Low each (skip limma if insufficient)
MAKE_PLOTS <- FALSE

## ---- Global switch: whether to run High/Low sensitivity (default FALSE; main analysis uses continuous predictor) ----
RUN_HILO_SENSITIVITY <- FALSE

## ---- [USER CONFIG | Inline settings] ----
## Analysis version: BatchAdj (batch-adjusted)
RUN_PASSES <- c("BatchAdj")

## Write inline settings to options
options(csn.run_passes = RUN_PASSES)
## ---- [END USER CONFIG] ----

## ---- [DPS config] ----
## Set TRUE to run DPS only (skip PTM-SEA)
RUN_DPS_ONLY <- TRUE
RUN_DPS_ONLY <- get0("RUN_DPS_ONLY", ifnotfound = FALSE)

## DPS output prefix (hardcoded for consistent output paths)
COMBO_PREFIX_DPS <- "phosphoproteomic_DPS"
## ---- [END DPS config] ----


## ---- [USER CONFIG | heatmap top/bottom limits] ----
## Non-H collection heatmap filtering threshold (single dataset uses NES, PAN uses Z; both filter by CSN_SCORE)
## - Single dataset: keep only top N and bottom M NES with padj<0.05 (by CSN_SCORE)
## - PAN (meta-FDR): keep only top N and bottom M Z with padj_meta<0.05 (by CSN_SCORE)
DATASET_HEATMAP_TOP_N <- 25 # Default top 25 (NES)
DATASET_HEATMAP_BOTTOM_N <- 25 # Default bottom 25 (NES)
PAN_HEATMAP_TOP_N <- 25 # Default top 25 (Z)
PAN_HEATMAP_BOTTOM_N <- 25 # Default bottom 25 (Z)
## ---- [END USER CONFIG | heatmap top/bottom limits] ----


## ---- [USER CONFIG | heatmap collections] ----
## Which collections to plot heatmaps for (character vector; e.g., c("H", "C6", "C2:CP:BIOCARTA"))
## - Single dataset heatmap: leave empty or NULL for "plot all available collections"
## - PAN heatmap: leave empty or NULL for "plot all available collections"
PLOT_DATASET_COLLECTIONS <- NULL
PLOT_PAN_COLLECTIONS <- NULL
## ---- [END USER CONFIG | heatmap collections] ----


## ---- Pipeline: limma t-statistic only ----
PIPELINES_TO_RUN <- c("limma_t")
.RUN_LIMMA <- TRUE

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

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_fs_name <- function(s) {
  s <- gsub('[<>:"/\\\\|?*]', "_", s)
  s <- gsub("\\s+", "_", s)
  s
}




## ---- [NEW | PTMsigDB helper for site-level PTM-SEA] ----
.std_site_id <- function(x) {
  # Standardize to GENE_S123 / GENE_T123 / GENE_Y123
  x <- toupper(as.character(x))
  x <- gsub("[:\\s]", "_", x) # Convert to underscore
  x <- gsub("_([STY])(\\d+)[A-Z]*$", "_\\1\\2", x) # Remove trailing lowercase or other suffixes
  x
}

# Generate site ID from CPTAC phospho annotation; prioritize NAME, otherwise use GENE + PHOSPHOSITES
.make_site_id <- function(gene, name, phosphosites) {
  gene <- as.character(gene)
  name <- as.character(name)
  phosphosites <- as.character(phosphosites)

  row_fun <- function(g, n, p) {
    cand <- NA_character_
    if (!is.na(n) && nzchar(n) && grepl("_[STY]\\d+", n, ignore.case = TRUE)) cand <- n
    if (is.na(cand) || !nzchar(cand)) {
      ps <- NA_character_
      if (!is.na(p) && nzchar(p)) {
        m <- regmatches(p, regexpr("([STY])(\\d+)", p, ignore.case = TRUE))
        if (length(m) == 1 && nzchar(m)) ps <- m
      }
      if (!is.na(g) && nzchar(g) && !is.na(ps) && nzchar(ps)) cand <- paste0(g, "_", ps)
    }
    .std_site_id(cand)
  }

  vapply(seq_along(gene), function(i) row_fun(gene[i], name[i], phosphosites[i]), character(1))
}

# Automatically infer GENE_SITE using available columns:
# 1) If NAME contains GENE_S/T/Y### (allowing :isoform suffix) -> use NAME (remove after colon)
# 2) Otherwise if GENE_SYMBOL + PHOSPHOSITES/PHOSPHOSITE exist -> combine into GENE_SITE
# 3) Otherwise parse from ENTITY_STABLE_ID (supports AAAS_pS495, AAAS_S495:NP_... etc.)
.infer_site_id <- function(df) {
  g <- if ("GENE_SYMBOL" %in% names(df)) as.character(df$GENE_SYMBOL) else rep(NA_character_, nrow(df))
  nm <- if ("NAME" %in% names(df)) as.character(df$NAME) else rep(NA_character_, nrow(df))
  p1 <- if ("PHOSPHOSITES" %in% names(df)) as.character(df$PHOSPHOSITES) else rep(NA_character_, nrow(df))
  p2 <- if ("PHOSPHOSITE" %in% names(df)) as.character(df$PHOSPHOSITE) else rep(NA_character_, nrow(df))
  ent <- if ("ENTITY_STABLE_ID" %in% names(df)) as.character(df$ENTITY_STABLE_ID) else rep(NA_character_, nrow(df))

  # 1) from NAME
  id_from_name <- ifelse(!is.na(nm) & grepl("_[STY]\\d+", nm, ignore.case = TRUE),
    sub(":.*$", "", nm), NA_character_
  )

  # helper: extract first site token (S/T/Y + number)
  site_token <- function(x) {
    y <- NA_character_
    if (!is.na(x) && nzchar(x)) {
      m <- regmatches(x, regexpr("([STY])(\\d+)", x, ignore.case = TRUE))
      if (length(m) == 1 && nzchar(m)) y <- toupper(m)
    }
    y
  }

  # 2) from GENE + (PHOSPHOSITES/PHOSPHOSITE)
  ps <- vapply(seq_along(p1), function(i) {
    s <- if (!is.na(p1[i]) && nzchar(p1[i])) p1[i] else p2[i]
    site_token(s)
  }, character(1))
  id_from_gp <- ifelse(!is.na(g) & nzchar(g) & !is.na(ps) & nzchar(ps),
    paste0(g, "_", ps), NA_character_
  )

  # 3) from ENTITY_STABLE_ID
  ent_base <- toupper(sub(":.*$", "", ent))
  ent_gene <- ifelse(!is.na(g) & nzchar(g), toupper(g), sub("_.*$", "", ent_base))
  ent_site <- vapply(ent_base, function(x) {
    m <- regmatches(x, regexpr("([STY])(\\d+)", x, ignore.case = TRUE))
    if (length(m) == 1 && nzchar(m)) toupper(m) else NA_character_
  }, character(1))
  id_from_ent <- ifelse(!is.na(ent_gene) & nzchar(ent_gene) & !is.na(ent_site) & nzchar(ent_site),
    paste0(ent_gene, "_", ent_site), NA_character_
  )

  # Backfill in order
  cand <- id_from_name
  need <- is.na(cand) | !nzchar(cand)
  cand[need] <- id_from_gp[need]
  need <- is.na(cand) | !nzchar(cand)
  cand[need] <- id_from_ent[need]

  .std_site_id(cand)
}

# Read PTMsigDB GMT (each line: <set>\t<desc>\t<site1>\t<site2>...)
read_ptmsigdb_gmt <- function(fp) {
  if (missing(fp) || !nzchar(fp) || !file.exists(fp)) {
    stop("[PTMsigDB] File does not exist: ", fp %||% "<empty>")
  }

  std_id <- function(x) {
    # Standardization consistent with .std_site_id: keep GENE_SITE format
    x <- toupper(trimws(x))
    x <- sub(":.*$", "", x) # Remove PMID etc. annotations (after :)
    x <- gsub("[^A-Z0-9_]", "", x) # Conservative cleanup
    x
  }

  # ---- A) If .xlsx (with site.annotation column) -> use gene-symbol_site directly ----
  if (grepl("\\.xlsx$", fp, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("[PTMsigDB] Package readxl required to read .xlsx")
    }
    df <- readxl::read_xlsx(fp)
    req <- c("signature", "site.annotation")
    if (!all(req %in% names(df))) {
      stop("[PTMsigDB] .xlsx missing required columns: ", paste(setdiff(req, names(df)), collapse = ", "))
    }
    # site.annotation format: PPP1R12A_T696:15226371;20801872
    df$gene_site <- std_id(df$site.annotation)
    by_sig <- split(df$gene_site, df$signature)
    out <- lapply(by_sig, function(v) unique(v[nzchar(v)]))
    out[lengths(out) == 0] <- NULL
    return(out)
  }

  # ---- B) If .gmt -> compatible with two common formats ----
  lines <- readr::read_lines(fp)
  res <- vector("list", length(lines))
  nm <- character(length(lines))

  for (i in seq_along(lines)) {
    fields <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    if (length(fields) < 2) next
    sig <- fields[1]
    nm[i] <- sig

    # Normal GMT: sites start from column 3; some variants list gene-symbol in column 2, separated by '|'
    toks <- if (length(fields) >= 3) fields[-(1:2)] else character(0)
    if (length(toks) == 0 && grepl("\\|", fields[2], fixed = TRUE)) {
      t2 <- strsplit(fields[2], "\\|")[[1]]
      toks <- t2[!grepl("^n=\\d+$", t2, ignore.case = TRUE)]
    }

    # Remove annotation part (after :, e.g., PMID string), standardize
    toks <- std_id(toks)

    # Keep only gene-symbol_site format (GENE_[STY]pos)
    keep <- grepl("^[A-Z0-9._-]+_[STY][0-9]+$", toks)
    toks <- unique(toks[keep & nzchar(toks)])

    res[[i]] <- toks
  }

  names(res) <- nm
  res[lengths(res) == 0] <- NULL
  # If all empty, usually means this GMT is **UniProt;site** (e.g., O14974;T696), please use .xlsx or gene version GMT
  if (!length(res)) {
    stop(
      "[PTMsigDB] Parsing result is empty. This GMT is likely UniProt version (e.g., O14974;T696).",
      " Please use source with gene-symbol (.xlsx or gene version GMT)."
    )
  }
  res
}

## ===== PTMsigDB Configuration =====
## PTMsigDB v2.0.0 (Krug et al., 2019) - required for DPS heatmap annotation
## Download from: https://proteogenomics.shinyapps.io/ptmsigdb/ptmsigdb.Rmd
## File: data_PTMsigDB_all_sites_v2.0.0.xlsx (place in the same directory as this script)
PTMSIGDB_GMT_FP <- file.path(getwd(), "data_PTMsigDB_all_sites_v2.0.0.xlsx")


genesets_by_group_ptm <- list()
if (nzchar(PTMSIGDB_GMT_FP) && file.exists(PTMSIGDB_GMT_FP)) {
  genesets_by_group_ptm[["PTMsigDB"]] <- read_ptmsigdb_gmt(PTMSIGDB_GMT_FP)
  log_msg("Loaded PTMsigDB site collections: %d sets", length(genesets_by_group_ptm[["PTMsigDB"]]))
} else {
  log_msg("(Note) PTMSIGDB_GMT not set; site-level PTM-SEA will not be available.")
}
## ---- [END NEW] ----


## ===== File reading and shared utilities =====
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
  # If directory has case_list, try to match
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

# Safe version CSN SCORE (try original method first; fallback on failure)
build_csn_score_safe <- function(mat0, subunits, combine_7AB = TRUE,
                                 min_members = 5L, pca_min_samples = 10L) {
  sub <- intersect(subunits, rownames(mat0))
  out_na <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
  if (length(sub) < min_members) {
    return(out_na)
  }

  cs_try <- try(
    {
      build_csn_score(prot0, subunits = sub, combine_7AB = combine_7AB, min_members = min_members)
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
      "[CSN_SCORE-safe] fallback：genes=%d；PC1%%=%.1f；nonNA=%d/%d",
      nrow(X), 100 * varpc1, sum(is.finite(out)), length(out)
    )
  }
  out
}

# Safe version audit: try original audit first; on failure use clean matrix for prcomp (record only, don't stop flow)
audit_csn_score_feasibility_safe <- function(ds_id, stratum, mat0, prot0, present_sub,
                                             min_members = 5L, pca_min_samples = 10L,
                                             min_per_group = get0("min_per_group", ifnotfound = 8L),
                                             out_dir = file.path("run_info", "csn_score_audit")) {
  ok <- try(
    {
      audit_csn_score_feasibility(
        ds_id = ds_id,
        stratum = stratum,
        mat0 = mat0,
        prot0 = prot0,
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

  X <- .clean_for_pca(prot0[intersect(present_sub, rownames(prot0)), , drop = FALSE],
    min_samples = pca_min_samples, min_genes = min_members
  )
  if (is.null(X)) {
    if (exists("log_msg", mode = "function")) {
      log_msg("[CSN-audit-safe] %s | %s: Insufficient subunits or samples, skipping audit (not stopping)", ds_id, stratum)
    }
    return(invisible(FALSE))
  }
  pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) {
    if (exists("log_msg", mode = "function")) {
      log_msg("[CSN-audit-safe] %s | %s: fallback prcomp still failed, skipping (not stopping)", ds_id, stratum)
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

# --- helper: ensure limma t-stat has correct gene names ---
._ensure_stats_names <- function(stats, gene_names, label = NULL) {
  # stats: numeric vector (limma t, etc.)
  # gene_names: corresponding matrix rownames (gene IDs)
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
  # Don't filter yet, leave for next step cleanup (will log)
  v
}



.align_limma_inputs <- function(M, design_df, rhs_terms, label = "limma") {
  stopifnot(!is.null(M), is.matrix(M), ncol(M) > 0, !is.null(colnames(M)))
  # First convert design table to data.frame, ensure rownames are sample IDs
  design_df <- as.data.frame(design_df, check.names = FALSE)
  stopifnot(!is.null(rownames(design_df)))
  samp_M <- as.character(colnames(M))
  samp_des <- as.character(rownames(design_df))

  # First remove NA samples from design table (predictor / numeric covariates / factor NA all considered unavailable)
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


## ===== impute_and_filter: fix seed before imputation (5) =====
impute_and_filter <- function(mat, min_frac = 0.75) {
  keep <- rowMeans(!is.na(mat)) >= min_frac
  m <- mat[keep, , drop = FALSE]
  if (any(is.na(m))) {
    set.seed(1234)
    m <- imputeLCMD::impute.MinProb(m, q = 0.01)
  }
  m
}

## ===== Coverage and covariate audit: audit_covars_coverage() =====
## Recommended to place before run_one_stratum()
## Will append summary to run_info/covars_audit/audit_rows.csv and print a log line
if (!exists(".ensure_mat_or_null", mode = "function")) {
  .ensure_mat_or_null <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    m <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (is.null(m)) {
      return(NULL)
    }
    if (is.null(dim(m))) { # vector -> single column matrix
      m <- matrix(m, ncol = 1)
    }
    if (ncol(m) == 0) {
      return(NULL)
    }
    m
  }
}

## ==== PATCH: helpers ====


## === NEW: CSN complex score (PC1) & residualization helpers ==================

# Build CSN complex score (PC1) using cross-sample z-values of subunits; direction corrected to match average z sign
build_csn_score <- function(mat0,
                            subunits = csn_subunits,
                            combine_7AB = TRUE,
                            min_members = 5L) {
  present <- intersect(subunits, rownames(mat0))
  # Pre-build return skeleton (always has names)
  s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
  if (!length(present)) {
    return(s)
  }

  # z-score (keep sample names)
  get_z <- function(v) {
    nm <- names(v) # Store sample names first
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
    # Direction correction: same sign as average subunit z
    mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)
    if (suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs")) < 0) sc <- -sc
    s[keep_sam] <- sc # Here alignment by sample name will definitely succeed
  }

  s
}



audit_covars_coverage <- function(tag, ds_id, stratum, su,
                                  sample_ids,
                                  batch = NULL,
                                  covars = NULL) {
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
      mean(!is.na(v)) * 100 # Factor/string -> use non-NA proportion
    }
  }
  cov_purity <- get_cov("purity")
  cov_sex <- get_cov("sex")
  cov_age <- get_cov("age")

  line <- sprintf(
    "  [audit:%s] covars={%s} | coverage: purity=%.1f%%, sex=%.1f%%, age=%.1f%% | batch_levels=%s",
    tag,
    if (length(cov_nms)) paste(cov_nms, collapse = ",") else "NULL",
    cov_purity %||% NaN, cov_sex %||% NaN, cov_age %||% NaN,
    if (!is.null(batch)) nlevels(batch) else 0L
  )
  log_msg(line)

  # Append one row to summary table
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


## ---- Auto-detect batch column & read (expanded candidate columns) ----
## Small utilities (only used in this file)
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
  # If provided is named vector/factor/data.frame, align by colnames(mat); otherwise assume already aligned
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
  # Optional: merge very small batches (default single instance merged, set 0 to disable)
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

## ===== Batch cleanup settings (you can adjust as needed) =====
BATCH_PIPE_POLICY <- "NA" # "NA" or "b_small": values with '|' should be set to NA or merged to b_small
BATCH_MIN_PER_LEVEL <- 2 # Levels below this threshold will be merged to b_small (avoid non-estimable coefficients)

## ---- Batch value cleanup: handle '|', legalize names, merge sparse levels ----
sanitize_batch_levels <- function(x,
                                  pipe_policy = BATCH_PIPE_POLICY,
                                  min_per_level = BATCH_MIN_PER_LEVEL) {
  x0 <- as.character(x)
  has_pipe <- grepl("\\|", x0 %||% "")
  if (any(has_pipe)) {
    log_msg("  [batch] Detected {sum(has_pipe)} values with '|' -> processing per policy {pipe_policy}")
    x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
  }
  # Legalize names (avoid eBayes/design matrix column name issues)
  fac <- factor(make.names(x0))
  fac <- droplevels(fac)

  # Merge sparse levels (e.g., only 1 sample)
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

## ---- Auto-detect batch column (with cleanup pipeline) ----
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
    # Need at least 2 valid levels and valid sample count >= 3
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

  # First align by sample_ids to ensure detection and returned vector have same length
  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
  rownames(meta) <- sample_ids

  det <- detect_batch_column(meta,
    pipe_policy   = pipe_policy,
    min_per_level = min_per_level
  )
  if (is.null(det)) {
    ## === [NEW] Fallback: infer TMT-plex from TMT_protein.csv as batch ===
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

        ## Build sample_id -> plex mapping (each row is a plex; all samples in tmt_* columns belong to this plex)
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
          ## Success condition: at least two valid levels and valid sample count >= 3
          if (nlevels(fac2) >= 2 && sum(!is.na(fac2)) >= 3) {
            names(fac2) <- sample_ids
            return(list(name = "TMT_protein.csv:RunMetadataID", fac = fac2))
          }
        }
      }
    }
    ## Fallback also failed -> return NULL (maintain original semantics)
    return(NULL)
  }

  fac <- det$fac # Already aligned with sample_ids
  names(fac) <- sample_ids
  list(name = det$name, fac = fac)
}


# --- [NEW] phospho version batch get: clinical empty -> fallback to TMT_phos.csv ---------------
get_batch_factor_phospho <- function(ds_dir, sample_ids,
                                     pipe_policy = c("strict", "lenient", "NA"),
                                     min_per_level = 3L) {
  pipe_policy <- match.arg(pipe_policy)
  # First use clinical detection (same as protein version)
  meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  if (file.exists(meta_fp)) {
    meta <- suppressMessages(
      readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")
    ) |> as.data.frame()
    id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"), names(meta))
    if (length(id_cols)) {
      meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
      meta$SAMPLE_ID <- as.character(meta$SAMPLE_ID)
      ## Same as protein version: use match to align and preserve sample_ids order and length
      meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
      rownames(meta) <- sample_ids
      det <- detect_batch_column(meta,
        pipe_policy   = pipe_policy,
        min_per_level = min_per_level
      )
      if (!is.null(det)) {
        fac <- det$fac
        names(fac) <- sample_ids
        return(list(name = det$name, fac = fac))
      }
    }
  }

  # FALLBACK: Read TMT_phos.csv
  tmt_fp <- file.path(ds_dir, "TMT_phos.csv")
  if (!file.exists(tmt_fp)) {
    return(NULL)
  }

  tmt <- suppressMessages(readr::read_csv(tmt_fp, show_col_types = FALSE)) |> as.data.frame()
  cn <- names(tmt)
  cn_trim <- trimws(cn)
  names(tmt) <- cn_trim
  run_hits <- grep("^Run\\s*Metadata\\s*ID$", cn_trim, ignore.case = TRUE, value = TRUE)
  tmt_cols <- grep("^tmt_", cn_trim, ignore.case = TRUE, value = TRUE)
  if (!length(run_hits) || !length(tmt_cols)) {
    return(NULL)
  }

  run_col <- run_hits[1]
  plex_by_sample <- list()
  for (i in seq_len(nrow(tmt))) {
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

  if (!length(plex_by_sample)) {
    return(NULL)
  }
  v <- rep(NA_character_, length(sample_ids))
  names(v) <- sample_ids
  mm <- intersect(names(plex_by_sample), sample_ids)
  if (length(mm)) v[mm] <- unlist(plex_by_sample[mm], use.names = FALSE)

  fac2 <- sanitize_batch_levels(v,
    pipe_policy   = pipe_policy,
    min_per_level = min_per_level
  )
  if (nlevels(fac2) >= 2 && sum(!is.na(fac2)) >= 3) {
    names(fac2) <- sample_ids
    return(list(name = "TMT_phos.csv:RunMetadataID", fac = fac2))
  }
  return(NULL)
}


## ---- Batch need screening (using your logic) ----
screen_batch_need <- function(ds_dir, min_frac_complete = 0.75) {
  log_msg("== Batch check: %s ==", basename(ds_dir))
  mat0 <- load_phospho_matrix_from_dataset_dir(ds_dir)
  mx <- suppressWarnings(max(mat0, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)
  m <- impute_and_filter(mat0, min_frac = min_frac_complete)

  bi <- get_batch_factor_phospho(ds_dir, colnames(m))
  if (is.null(bi)) {
    log_msg("  [batch] Cannot find clear batch column -> no correction for now")
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

  ## Safe PCA: first clean NA/Inf, keep at least 10 valid samples and >=5 genes; otherwise skip PCA
  X <- .clean_for_pca(m, min_samples = 10L, min_genes = 5L)
  if (is.null(X)) {
    r2 <- pv <- numeric(0)
    log_msg("  PCA (by batch) skipped: insufficient genes/samples or still contains non-finite values")
  } else {
    Xs <- t(scale(t(X), center = TRUE, scale = TRUE))
    pc <- stats::prcomp(t(Xs), scale. = FALSE)
    K <- min(5L, ncol(pc$x))
    r2 <- vapply(seq_len(K), function(i) summary(lm(pc$x[, i] ~ batch))$r.squared, numeric(1))
    pv <- vapply(seq_len(K), function(i) {
      a <- anova(lm(pc$x[, i] ~ batch))
      as.numeric(a$`Pr(>F)`[1])
    }, numeric(1))
    log_msg(
      "  PCA (by batch) R2: %s ; p: %s", paste(round(r2, 3), collapse = ", "),
      paste(signif(pv, 3), collapse = ", ")
    )
  }

  design <- model.matrix(~batch)
  fit <- limma::lmFit(m, design)
  fit <- limma::eBayes(fit)
  padj <- p.adjust(fit$F.p.value, "BH")
  prop05 <- mean(padj < 0.05, na.rm = TRUE)
  prop25 <- mean(padj < 0.25, na.rm = TRUE)
  log_msg(
    "  Gene-level F test: FDR<0.05 proportion = %.1f%%; FDR<0.25 = %.1f%%",
    100 * prop05, 100 * prop25
  )

  recommend <- (any(r2 >= 0.10 & pv[seq_along(r2)] < 0.01) || prop05 >= 0.05)
  if (recommend) {
    log_msg("  **Recommend correction**: R2 or gene proportion reached threshold (>=10%% R2 and p<0.01, or FDR<0.05 genes >=5%%)")
  } else {
    log_msg("  **No correction needed for now**: No obvious batch effect detected (record kept)")
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

dataset_ids <- c(
  "brca_cptac_2020",
  "coad_cptac_2019",
  "gbm_cptac_2021",
  "luad_cptac_2020",
  "lusc_cptac_2021",
  "paad_cptac_2021",
  "ucec_cptac_2020"
)

## PTM-SEA output prefix (hardcoded for consistent output paths)
COMBO_PREFIX <- "phospho_combo_2"

## Pipeline settings (limma_t only)
PIPELINES_TO_RUN <- c("limma_t")
.RUN_LIMMA <- TRUE



dataset_dirs <- setNames(file.path(datasets_root, dataset_ids), dataset_ids)
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")
message("datasets_root = ", datasets_root)

## Old version was stopifnot(all(dir.exists(dataset_dirs))), which would stop entire batch if some directories missing.
## Changed to: list missing ones, use dataset_dirs_run to filter "actually runnable" set for this run.
missing_dirs <- names(dataset_dirs)[!dir.exists(dataset_dirs)]
if (length(missing_dirs)) {
  log_msg("Detected %d missing folders, will skip: %s", length(missing_dirs), paste(missing_dirs, collapse = ", "))
}


## List of datasets to actually run and available this run (folder exists and contains data_phosphoprotein_quantification.txt)
dataset_dirs_run <- dataset_dirs[
  dir.exists(dataset_dirs) &
    file.exists(file.path(dataset_dirs, "data_phosphoprotein_quantification.txt"))
]

if (!length(dataset_dirs_run)) stop("dataset_dirs_run is empty, please verify folders and data_phosphoprotein_quantification.txt exist")

log_msg("Available datasets this run (%d): %s", length(dataset_dirs_run), paste(names(dataset_dirs_run), collapse = ", "))


# --- [NEW] phospho version matrix loading (aggregated by gene symbol) -----------------------------
load_phospho_matrix_from_dataset_dir <- function(dir,
                                                 site_collapse = c("median", "max"),
                                                 protein_adjust = FALSE,
                                                 protein_fp = file.path(dir, "data_protein_quantification.txt")) {
  site_collapse <- match.arg(site_collapse)
  fp <- file.path(dir, "data_phosphoprotein_quantification.txt")
  if (!file.exists(fp)) stop("[phospho] File does not exist: ", fp)

  suppressMessages({
    df <- readr::read_tsv(fp, show_col_types = FALSE, progress = FALSE)
  })

  # Compatible columns: ENTITY_STABLE_ID, NAME, DESCRIPTION, GENE_SYMBOL, PHOSPHOSITES, <samples...>
  meta_cols <- intersect(
    colnames(df),
    c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "GENE_SYMBOL", "PHOSPHOSITES")
  )
  sample_cols <- setdiff(colnames(df), meta_cols)
  if (!length(sample_cols)) stop("[phospho] Cannot find sample columns")

  # --- site x sample matrix (using GENE_SYMBOL as group key), values are already log2 ratio (CPTAC standard)
  M_site <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(M_site) <- "numeric"
  rownames(M_site) <- df$GENE_SYMBOL

  # --- Aggregate to gene x sample (consistent interface with protein version)
  split_idx <- split(seq_len(nrow(M_site)), rownames(M_site))
  agg_fun <- if (site_collapse == "median") {
    function(x) stats::median(x, na.rm = TRUE)
  } else {
    function(x) suppressWarnings(max(x, na.rm = TRUE))
  }
  M_gene <- do.call(rbind, lapply(split_idx, function(ii) {
    apply(M_site[ii, , drop = FALSE], 2, agg_fun)
  }))
  M_gene <- M_gene[order(rownames(M_gene)), , drop = FALSE]

  # --- Optional: adjust by protein abundance (remove same-gene protein effect; see "Method description" section)
  if (isTRUE(protein_adjust) && file.exists(protein_fp)) {
    suppressMessages({
      prot_df <- readr::read_tsv(protein_fp, show_col_types = FALSE, progress = FALSE)
    })
    # Protein file sample columns: exclude common annotation columns
    prot_meta <- intersect(
      colnames(prot_df),
      c(
        "Hugo_Symbol", "GENE_SYMBOL", "GeneSymbol", "Description",
        "GENE", "GENE_ID", "GENE.STABLE.ID", "Composite.Element.REF"
      )
    )
    prot_sample <- setdiff(colnames(prot_df), prot_meta)
    # Try to find gene symbol column (including Composite.Element.REF; case insensitive)
    cand <- c("Hugo_Symbol", "GENE_SYMBOL", "GeneSymbol", "GENE", "Composite.Element.REF")
    cn <- colnames(prot_df)
    lc <- tolower(cn)
    hit <- match(tolower(cand), lc, nomatch = NA_integer_)
    hit <- hit[!is.na(hit)]
    gcol <- if (length(hit)) cn[hit[1]] else NA_character_
    if (is.na(gcol)) stop("[phospho] protein_adjust=TRUE but protein file missing gene column")
    # Ensure sample columns don't contain gene column
    prot_sample <- setdiff(prot_sample, gcol)
    M_prot <- as.matrix(prot_df[, prot_sample, drop = FALSE])
    storage.mode(M_prot) <- "numeric"
    rownames(M_prot) <- prot_df[[gcol]]
    rownames(M_prot) <- sub("\\|.*$", "", rownames(M_prot))

    # Gene and sample intersection
    g_common <- intersect(rownames(M_gene), rownames(M_prot))
    s_common <- intersect(colnames(M_gene), colnames(M_prot))
    if (length(g_common) >= 1 && length(s_common) >= 3) {
      Y <- M_gene[g_common, s_common, drop = FALSE]
      X <- M_prot[g_common, s_common, drop = FALSE]
      Y_adj <- Y
      for (g in g_common) {
        y <- as.numeric(Y[g, ])
        x <- as.numeric(X[g, ])
        ok <- is.finite(y) & is.finite(x)
        if (sum(ok) >= 3) {
          fit <- stats::lm(y[ok] ~ x[ok])
          res <- rep(NA_real_, length(y))
          res[ok] <- stats::residuals(fit)
          Y_adj[g, ] <- res
        } else {
          Y_adj[g, ] <- y
        }
      }
      M_gene[g_common, s_common] <- Y_adj
      attr(M_gene, "protein_adjusted") <- TRUE
    }
  }

  # Sample name minor cleanup (consistent with protein version)
  colnames(M_gene) <- gsub("\\s+", "", colnames(M_gene))
  if (exists("log_msg")) log_msg("[phospho] Matrix dimensions: %d genes x %d samples", nrow(M_gene), ncol(M_gene))
  return(M_gene)
}


.build_design_for_version <- function(version, pred_centered, sample_order,
                                      batch_all, purity_all, sa_all,
                                      USE_AGE_MISSING_INDICATOR = FALSE,
                                      ds_id = NULL, mat_for_sv = NULL) {
  # Consistent with run_predictor_analyses design: base covariates
  DF <- data.frame(
    pred_centered = as.numeric(pred_centered[sample_order]),
    row.names = sample_order, check.names = FALSE
  )
  if (!is.null(sa_all)) {
    if ("sex" %in% colnames(sa_all)) DF$sex <- factor(sa_all[sample_order, "sex"])
    if ("age" %in% colnames(sa_all)) DF$age <- suppressWarnings(as.numeric(sa_all[sample_order, "age"]))
    if (USE_AGE_MISSING_INDICATOR) {
      if ("age_missing" %in% colnames(sa_all)) DF$age_missing <- suppressWarnings(as.numeric(sa_all[sample_order, "age_missing"]))
      if ("age_z_imputed" %in% colnames(sa_all)) DF$age_z_imputed <- suppressWarnings(as.numeric(sa_all[sample_order, "age_z_imputed"]))
    }
  }
  if (identical(version, "BatchAdj")) {
    if (!is.null(batch_all)) DF$batch <- droplevels(batch_all[sample_order])
  }
  # Clean singular or empty columns
  DF <- coerce_covariates_safely(DF)
  DF
}

## ---- Small utility: align any vector/data.frame to specified sample order ----
.align_to_samples <- function(x, sam, what = "covariate") {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.vector(x) || is.factor(x)) {
    # Allow length exactly equal to sam without names; else prefer align by name
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

# Any object -> matrix; force set rownames; return NULL if 0 columns or row count doesn't match specified rows
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
# More general safe matrix conversion (allows vector, 1 column, or zero length; doesn't specify rows)
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


## Safely convert vector to factor, and make NA explicit as "NA" (then add prefix later)
.factorize_with_explicit_NA <- function(x) {
  x_chr <- as.character(x)
  x_chr[!nzchar(x_chr) | is.na(x_chr)] <- "NA" # Empty string/NA both treated as "NA"
  factor(x_chr)
}

## ---- limma t-rank (can merge batch / other covariates; force alignment) ----
## ==== PATCH: limma t with covariates (guard p~n) ====
limma_t_with_covars <- function(mat, grp2, batch = NULL, covars = NULL) {
  stopifnot(!is.null(colnames(mat)))
  sam <- colnames(mat)

  if (is.null(names(grp2))) names(grp2) <- sam
  grp2 <- factor(as.character(grp2[sam]), levels = c("Low", "High"))
  if (anyNA(grp2)) stop("grp2 still has NA, please handle group samples first")

  des_grp <- model.matrix(~ 0 + grp2)
  colnames(des_grp) <- c("Low", "High")

  Xb <- NULL
  if (!is.null(batch)) {
    b <- .align_to_samples(batch, sam, what = "batch")
    b <- .factorize_with_explicit_NA(b)
    levels(b) <- paste0("b", make.names(levels(b)))
    if (nlevels(b) >= 1) {
      Xb <- model.matrix(~ 0 + b)
      colnames(Xb) <- levels(b)
    }
  }

  Xc <- NULL
  if (!is.null(covars)) {
    cv <- .align_to_samples(covars, sam, what = "covars")
    cv <- as.matrix(cv)
    if (!is.null(cv) && nrow(cv) == length(sam) && ncol(cv) > 0) Xc <- cv
  }
  if (!is.null(Xc)) {
    Xc[!is.finite(Xc)] <- 0
  }

  parts <- list(des_grp)
  if (!is.null(Xb) && nrow(Xb) == length(sam) && ncol(Xb) > 0) parts <- c(parts, list(Xb))
  if (!is.null(Xc) && nrow(Xc) == length(sam) && ncol(Xc) > 0) parts <- c(parts, list(Xc))
  design <- do.call(cbind, parts)
  storage.mode(design) <- "double"

  log_msg("  [design] Initial columns: %s", paste(colnames(design), collapse = ", "))
  log_msg("  [design] dim=%dx%d", nrow(design), ncol(design))

  # Allow NA; at least 3 finite values and variance>0
  var_ok <- vapply(seq_len(ncol(design)), function(j) {
    z <- design[, j]
    fin <- is.finite(z)
    if (sum(fin) < 3) {
      return(FALSE)
    }
    v <- var(z[fin])
    is.finite(v) && v > 0
  }, logical(1))
  if (!all(var_ok)) {
    drop <- setdiff(colnames(design)[!var_ok], c("Low", "High"))
    if (length(drop)) log_msg("  [design] Zero/near-zero variance -> delete: %s", paste(drop, collapse = ", "))
    design <- design[, var_ok, drop = FALSE]
  }

  if (!all(c("Low", "High") %in% colnames(design))) {
    log_msg("  [design] Low/High incomplete -> skip this limma")
    return(NULL)
  }

  # Remove collinearity
  q <- qr(design)
  keep_idx <- q$pivot[seq_len(q$rank)]
  keep_nms <- colnames(design)[keep_idx]
  removed <- setdiff(colnames(design), keep_nms)
  removed <- setdiff(removed, c("Low", "High"))
  if (length(removed)) log_msg("  [design] QR remove collinearity: remove %s", paste(removed, collapse = ", "))
  design <- design[, keep_nms, drop = FALSE]

  # p~n safeguard
  if (ncol(design) > (ncol(mat) - 2)) {
    keep <- intersect(
      colnames(design),
      c(
        "Low", "High",
        grep("^b", colnames(design), value = TRUE),
        grep("^SV", colnames(design), value = TRUE), # <- Keep SV
        "purity", "sex", "age", "PC1", "PC2"
      )
    )
    if (length(keep) >= 2) design <- design[, keep, drop = FALSE]
  }

  fit <- limma::lmFit(mat, design)
  fit2 <- limma::contrasts.fit(fit, contrasts = limma::makeContrasts(High - Low, levels = design))
  fit2 <- limma::eBayes(fit2)
  t <- fit2$t[, 1]
  t <- sort(t[is.finite(t)], decreasing = TRUE)
  attr(t, "design_cols") <- colnames(design)
  t
}

## ---- NEW: limma continuous predictor with covariates ----
limma_cont_with_covars <- function(mat, x, batch = NULL, covars = NULL, min_per_group = 8L) {
  stopifnot(!is.null(colnames(mat)))
  sam <- colnames(mat)

  if (is.null(names(x))) names(x) <- sam
  x <- as.numeric(x[sam])
  names(x) <- sam
  if (sum(is.finite(x)) < 4) stop("Continuous predictor has too few valid samples")

  # 1) Compose single DF -> model.matrix
  DF <- data.frame(x = x, row.names = sam, check.names = FALSE)
  if (!is.null(batch)) {
    b <- .align_to_samples(batch, sam, what = "batch")
    DF[["batch"]] <- droplevels(b)
  }
  if (!is.null(covars)) {
    cv <- .align_to_samples(covars, sam, what = "covars")
    cv <- as.data.frame(cv, check.names = FALSE)
    cv <- coerce_covariates_safely(cv) # <- Added
    for (cn in colnames(cv)) DF[[cn]] <- cv[[cn]] # <- Keep original fill into DF
  }
  ok <- stats::complete.cases(DF)
  if (sum(ok) < (2L * min_per_group)) stop(sprintf("Insufficient usable samples (%d < %d)", sum(ok), 2L * min_per_group))

  DF <- DF[ok, , drop = FALSE]
  M <- as.matrix(mat[, rownames(DF), drop = FALSE])
  des <- stats::model.matrix(~ 1 + ., data = DF, na.action = stats::na.fail)
  stopifnot(nrow(des) == ncol(M))

  fit <- limma::lmFit(M, des)
  # Coefficient name "x" always exists
  coef_idx <- which(colnames(des) == "x")
  if (length(coef_idx) != 1) stop("Cannot find x coefficient")
  # contrasts.fit requires contrast matrix; build a unit vector to extract x's t value
  C <- matrix(0, nrow = ncol(des), ncol = 1, dimnames = list(colnames(des), "x"))
  C["x", 1] <- 1
  fit2 <- limma::contrasts.fit(fit, C)
  fit2 <- limma::eBayes(fit2)

  t <- fit2$t[, 1]
  t <- sort(t[is.finite(t)], decreasing = TRUE)
  attr(t, "design_cols") <- colnames(des)
  t
}

## ===== Get TP53 mutation status =====
## ===== TP53 status (protein-altering baseline) =====
## Only count protein-altering variants as TP53-mutant; others (Silent, UTR, Intron, IGR, RNA, lincRNA, Flank...) all treated as wild type

TP53_KEEP_CLASSES <- c(
  "MISSENSE_MUTATION", "NONSENSE_MUTATION",
  "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
  "IN_FRAME_DEL", "IN_FRAME_INS",
  "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
)

normalize_vc <- function(x) {
  # Normalize Variant_Classification: uppercase, unify various separators to underscore
  x <- toupper(trimws(as.character(x)))
  gsub("[^A-Z0-9]+", "_", x)
}

get_tp53_status <- function(ds_dir, sample_ids) {
  # Default all wild-type
  status <- setNames(rep("TP53_wild_type", length(sample_ids)), sample_ids)

  mut_fp <- file.path(ds_dir, "data_mutations.txt")
  if (!file.exists(mut_fp)) {
    log_msg("  [TP53] Cannot find data_mutations.txt, treating entire batch as wild-type/ALL available")
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
##      (original Variant_Classification -> "sample count"; also contains Any_TP53_mutation and protein_altering)
##   - run_info/tp53_status/tp53_class_sample_counts_long.csv (aggregated long format)
##   - run_info/tp53_status/tp53_binary_counts_by_dataset.csv (wild type / mutant summary table)

dir.create(file.path("run_info", "tp53_status"), recursive = TRUE, showWarnings = FALSE)

# Safety: if not defined earlier, provide one here
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

# Single dataset TP53 sample count summary
summarize_tp53_counts_for_dataset <- function(ds_dir) {
  ds_id <- basename(ds_dir)

  # Use protein matrix samples as "population"
  M <- try(load_phospho_matrix_from_dataset_dir(ds_dir), silent = TRUE)
  if (inherits(M, "try-error")) {
    log_msg("[TP53-audit] %s: Cannot read matrix, skipping", ds_id)
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

  # Read MAF (original Variant_Classification distribution -> "sample count")
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

        # Each classification -> how many "unique samples" hit
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

  # Write individual files
  if (!is.null(class_df)) {
    data.table::fwrite(
      class_df,
      file.path("run_info", "tp53_status", paste0(ds_id, "_tp53_class_sample_counts.csv"))
    )
  } else {
    # If no MAF or no TP53 records, output empty shell
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
    log_msg("[TP53-audit] Starting: %s", ds)
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
    log_msg("[TP53-audit] Written: tp53_binary_counts_by_dataset.csv")
  }
  if (length(all_class)) {
    class_df <- dplyr::bind_rows(all_class)
    data.table::fwrite(class_df, file.path("run_info", "tp53_status", "tp53_class_sample_counts_long.csv"))
    log_msg("[TP53-audit] Written: tp53_class_sample_counts_long.csv")
  }
}


## ===== CSN subunits coverage rate & CSN_SCORE (PC1) feasibility audit (robust version) =====
audit_csn_score_feasibility <- function(ds_id, stratum, mat0, prot0, present_sub,
                                        min_members = 5L,
                                        pca_min_samples = 10L, # Same threshold as build_csn_score
                                        min_per_group = 8L, # Your limma threshold
                                        out_dir = file.path("run_info", "csn_score_audit")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  n_samples <- ncol(prot0)
  if (!is.matrix(mat0) || is.null(rownames(mat0)) || is.null(colnames(mat0))) {
    log_msg("[CSN-audit] %s | %s: mat0 structure incomplete, skip", ds_id, stratum)
    return(invisible(NULL))
  }
  present_sub <- intersect(present_sub, rownames(prot0))
  if (!length(present_sub) || n_samples == 0) {
    log_msg("[CSN-audit] %s | %s: No available CSN subunit or sample count is 0, skip", ds_id, stratum)
    return(invisible(NULL))
  }

  ## 1) Coverage rate for each subunit
  sub_cov <- vapply(present_sub, function(g) mean(is.finite(prot0[g, ])) * 100, numeric(1))
  sub_tbl <- data.frame(
    dataset = ds_id,
    stratum = stratum,
    subunit = present_sub,
    nonNA_pct = round(sub_cov, 1),
    nonNA_n = vapply(present_sub, function(g) sum(is.finite(prot0[g, ])), integer(1)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  cov_min <- if (length(sub_cov)) round(min(sub_cov, na.rm = TRUE), 1) else NA_real_
  cov_med <- if (length(sub_cov)) round(stats::median(sub_cov, na.rm = TRUE), 1) else NA_real_
  cov_max <- if (length(sub_cov)) round(max(sub_cov, na.rm = TRUE), 1) else NA_real_

  ## 2) How many subunits are non-NA for each sample; count of samples satisfying >= min_members
  sample_counts <- colSums(is.finite(prot0[present_sub, , drop = FALSE]))
  enough <- sample_counts >= min_members
  n_enough <- sum(enough)

  ## 3) Try to compute CSN_SCORE (PC1), and record feasibility
  csn_score <- build_csn_score(prot0,
    subunits = present_sub,
    combine_7AB = TRUE, min_members = min_members
  )
  csn_nonNA <- sum(is.finite(csn_score))
  csn_can_pca <- (n_enough >= pca_min_samples) && (csn_nonNA >= pca_min_samples)

  ## Additional: PC1 explained variance (when feasible)
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
    Z <- do.call(rbind, lapply(present_sub, function(g) get_z(prot0[g, ok_sam, drop = FALSE])))
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

  ## 4) Write out files
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

  # 4b) Subunit coverage rate (one file per stratum)
  tag <- paste(ds_id, stratum, sep = "_")
  fp_sub <- file.path(out_dir, sprintf("%s_subunit_coverage.csv", tag))
  data.table::fwrite(sub_tbl, fp_sub)

  # 4c) Non-NA subunit count per sample (one file per stratum)
  sample_tbl <- data.frame(
    dataset = ds_id, stratum = stratum,
    sample_id = colnames(prot0),
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


## Inspect current PLEX values (raw vs cleaned)
inspect_plex <- function(ds_dir, col = "TMT_PLEX",
                         pipe_policy = BATCH_PIPE_POLICY,
                         min_per_level = BATCH_MIN_PER_LEVEL) {
  # Read samples actually used (based on protein matrix)
  mat0 <- load_phospho_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(mat0)

  # Read clinical table and align samples
  meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  stopifnot(file.exists(meta_fp))
  meta <- suppressMessages(readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
  id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"), names(meta))
  stopifnot(length(id_cols) > 0)
  meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]

  # Get raw and cleaned PLEX values
  raw <- as.character(meta[[col]])
  clean <- sanitize_batch_levels(raw, pipe_policy = pipe_policy, min_per_level = min_per_level)

  # Summary
  cat("\n==== ", basename(ds_dir), " | Column: ", col, " ====\n", sep = "")
  cat("[Raw PLEX levels]:\n")
  print(sort(table(raw), decreasing = TRUE))
  cat("\n[Raw values containing '|' (first few)]:\n")
  print(utils::head(unique(raw[grepl("\\|", raw %||% "")]), 10))

  cat("\n[Cleaned PLEX levels] (pipe_policy = ", pipe_policy,
    ", min_per_level = ", min_per_level, "):\n",
    sep = ""
  )
  print(sort(table(clean, useNA = "ifany"), decreasing = TRUE))

  # Output mapping table (first 10 rows as demo)
  df_map <- data.frame(SAMPLE_ID = sample_ids, raw_plex = raw, clean_plex = as.character(clean))
  cat("\n[Sample mapping table (first 10 rows)]\n")
  print(utils::head(df_map, 10))

  invisible(list(
    raw_counts = sort(table(raw), decreasing = TRUE),
    clean_counts = sort(table(clean, useNA = "ifany"), decreasing = TRUE),
    map = df_map
  ))
}


## ===== Missingness sensitivity（AGE）=====
USE_AGE_MISSING_INDICATOR <- FALSE # Main analysis default: no missing-indicator, keep NA only

## Safe z-score: estimate mean/sd only from finite values, keep NA, no mean imputation
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

## Utility: column name normalization, z-score, 0~1 scaling
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

## For PAAD: get median from semicolon-separated values
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
## Get sex/age (sample aligned; sex: 0/1; age: z-score)
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

## For main workflow: strictly generate covariates from patient file SEX and AGE (sex: 0/1; age: z-score)
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
  # If .zscore (mean-imputed z) not defined externally, provide fallback version for consistent logic
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
    stop("[get_sex_age_covariates] Cannot establish sample-patient mapping (missing SAMPLE_ID/PATIENT_ID)")
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

  ## AGE: Main analysis no imputation; sensitivity can use missing-indicator + mean-imputed z
  age_col <- intersect(c("AGE", "AGE_AT_DIAGNOSIS", "AGE_AT_INDEX", "AGE_YEARS"), names(pat))[1]
  age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  age_missing <- age_z_imputed <- NULL
  if (!is.na(age_col)) {
    v <- suppressWarnings(as.numeric(pat[[age_col]]))
    names(v) <- as.character(pat[[pid_pat]])
    age_raw <- unname(v[map_pt[sample_ids]])
    age[] <- .z_no_impute(age_raw) # Keep NA (main analysis)
    if (!exists("USE_AGE_MISSING_INDICATOR", inherits = FALSE)) USE_AGE_MISSING_INDICATOR <- FALSE
    if (isTRUE(USE_AGE_MISSING_INDICATOR)) {
      age_missing <- as.numeric(is.na(age_raw))
      age_z_imputed <- .zscore(age_raw) # Mean-imputed z (sensitivity)
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
  ## NEW: Per-column drop reason collector
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
    logf(sprintf("  [covars-%s] policy: drop TP53 columns from covariates → %s", label, paste(tp_cols, collapse = ",")))
  }

  logf("  [covars-%s] after  align: C_dim=%d x %d", label, NROW(df), NCOL(df))


  # 2) Always use data.frame; preserve factors
  df <- as.data.frame(df, stringsAsFactors = TRUE, check.names = FALSE)

  # 3) Coverage: unnamed columns use global minimum threshold (safe indexing, not using [[cn]])
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  global_min <- if (length(min_cov_named)) min(min_cov_named, na.rm = TRUE) else 0
  keep_cov <- rep(TRUE, ncol(df))
  names(keep_cov) <- colnames(df)

  for (cn in colnames(df)) {
    v <- df[[cn]]
    # Unified coverage definition: numeric uses is.finite; non-numeric uses !is.na
    cover <- if (is.numeric(v)) mean(is.finite(v)) else mean(!is.na(v))

    thr <- if (!is.null(names(min_cov_named)) && (cn %in% names(min_cov_named))) {
      min_cov_named[[cn]] # Safe: only use [[cn]] when name exists
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

  ## NEW: Re-initialize keep_cov to match length and column names of the shrunk df
  keep_cov <- rep(TRUE, ncol(df))
  names(keep_cov) <- colnames(df)
  ## /NEW

  # 4) Correlation filtering with biology (only for numeric columns)
  if (!is.null(y)) {
    y <- suppressWarnings(as.numeric(y))
    ## Do not apply rho gate to these covariates
    skip_rho_gate <- c("purity", "age", "sex")

    for (cn in colnames(df)) {
      v <- df[[cn]]
      ## Only apply |rho| filtering for numeric columns not in skip list
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


  # 5) Return: row order = sample_order; for model.matrix expansion
  stopifnot(nrow(df) == length(so), identical(rownames(df), so))
  ## NEW: Attach drop reason for each column as attribute (caller decides whether to write CSV)
  attr(df, "drop_log") <- if (length(drop_log)) do.call(rbind, drop_log) else NULL
  df
}


## ==============================
## Mini audit: sex/age + purity + batch overall check
## ==============================

# Safe: utility functions (if not defined yet)
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

# --- NEW: Audit purity source and coverage by dataset rules ---
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
      # Try to estimate from immune/stroma (if columns exist)
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
## Mini audit: sex/age + purity + batch overall check (improved version)
## -- Added: policy/batch_min/limma_min to log and output columns
## ==============================
audit_one_dataset_sa_batch <- function(ds_dir, pipe_policy = "NA", min_per_level = 2) {
  ds_id <- basename(ds_dir)
  # 1) Get samples from protein matrix
  m <- load_phospho_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(m)

  # 2) Get SEX/AGE only from patient file
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

      # AGE (raw value statistics; z-score conversion in main workflow)
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

  # 3) Purity audit
  pur <- .audit_purity_for_dataset(ds_id, ds_dir, sample_ids)

  # 4) Batch check (using existing get_batch_factor_phospho)
  bi <- get_batch_factor_phospho(ds_dir, sample_ids,
    pipe_policy   = pipe_policy,
    min_per_level = min_per_level
  )
  batch_col <- if (!is.null(bi)) bi$name else "NONE"
  batch_levels <- if (!is.null(bi)) nlevels(bi$fac) else 0L
  batch_nonNA <- if (!is.null(bi)) mean(!is.na(bi$fac)) * 100 else NA_real_
  batch_sizes <- if (!is.null(bi)) .format_batch_sizes(bi$fac) else NA_character_

  # Print one-line summary (with policy/batch_min/limma_min)
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


audit_all_datasets_sa_batch <- function(dataset_dirs, pipe_policy = BATCH_PIPE_POLICY, min_per_level = 2) {
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


## ===== DO_AUDIT switch for QC/exploration (4) =====
DO_AUDIT <- FALSE # <-- Default off; change to TRUE when needed

if (DO_AUDIT) {
  results <- list()
  for (ds in names(dataset_dirs)) {
    ds_dir <- dataset_dirs[[ds]]
    if (!dir.exists(ds_dir)) {
      log_msg("Skipping: folder not found %s", ds_dir)
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
  strata = strata
), file.path("run_info", "run_manifest.yml"))


## Start re-running from this dataset (can be changed)
start_from <- "brca_cptac_2020"
ord <- dataset_ids
ix <- match(start_from, ord)
if (is.na(ix)) stop(sprintf("'%s' not in dataset_ids", start_from))
dataset_dirs_run <- dataset_dirs[ord[ix:length(ord)]]

exists("coerce_covariates_safely")
getAnywhere("coerce_covariates_safely")
exists("opt") # Should return TRUE
getAnywhere("opt") # Should show in .GlobalEnv


## ---- Optional: for testing, run only specific strata ----
# only_strata <- c("TP53_mutant")   # To run multiple strata, use c("ALL","TP53_mutant") etc.

## Utility: if only_strata not defined, run all by default
.should_run <- function(tag) {
  if (!exists("only_strata") || is.null(only_strata) || !length(only_strata)) {
    return(TRUE)
  }
  tag %in% only_strata
}



# --- [NEW | SITE-LEVEL] phospho site matrix loading (retain site; can do protein adjustment) ---
load_phosphosite_matrix_from_dataset_dir <- function(ds_dir,
                                                     protein_adjust = TRUE,
                                                     min_nonNA_row_frac = 0) {
  fp <- file.path(ds_dir, "data_phosphoprotein_quantification.txt")
  if (!file.exists(fp)) stop("Cannot find phosphoprotein file: ", fp)
  log_msg("Reading phospho (site-level) matrix: data_phosphoprotein_quantification.txt")

  df <- suppressMessages(vroom::vroom(fp, delim = "\t", col_types = vroom::cols(
    .default = "d",
    ENTITY_STABLE_ID = "c",
    NAME = "c",
    DESCRIPTION = "c",
    GENE_SYMBOL = "c",
    PHOSPHOSITES = "c", PHOSPHOSITE = "c"
  )))
  df <- as.data.frame(df, check.names = FALSE)

  # Detect annotation columns + sample columns
  anno_cols <- intersect(c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "GENE_SYMBOL", "PHOSPHOSITES", "PHOSPHOSITE"), colnames(df))
  samp_cols <- setdiff(colnames(df), anno_cols)
  stopifnot(length(samp_cols) > 0)

  # Generate site_id and normalize
  site_id <- .infer_site_id(df)
  keep <- !is.na(site_id) & nzchar(site_id)
  df <- df[keep, , drop = FALSE]
  site_id <- site_id[keep]

  M_site <- as.matrix(df[, samp_cols, drop = FALSE])
  storage.mode(M_site) <- "numeric"
  rownames(M_site) <- site_id
  colnames(M_site) <- samp_cols

  # (Optional) Filter by minimum non-missing ratio
  if (min_nonNA_row_frac > 0) {
    ok <- rowMeans(is.finite(M_site)) >= min_nonNA_row_frac
    M_site <- M_site[ok, , drop = FALSE]
  }

  # Protein adjustment: subtract corresponding gene's protein abundance (site residual-like)
  if (isTRUE(protein_adjust)) {
    prot_fp <- file.path(ds_dir, "data_protein_quantification.txt")
    if (!file.exists(prot_fp)) stop("[phospho-site] protein_adjust=TRUE but protein file not found")
    log_msg("Site-level protein adjustment: reading protein matrix: data_protein_quantification.txt")

    prot_df <- suppressMessages(vroom::vroom(prot_fp, delim = "\t"))
    prot_df <- as.data.frame(prot_df, check.names = FALSE)

    # Try to find protein gene column (consistent with your previous approach)
    gcol <- which(tolower(colnames(prot_df)) %in% tolower(c("Composite.Element.REF", "Gene", "GENE", "GENE_SYMBOL", "GENE_NAME")))[1]
    if (is.na(gcol)) stop("[phospho-site] protein_adjust=TRUE but protein file missing gene column")
    prot_genes <- as.character(prot_df[[gcol]])
    prot_mat <- as.matrix(prot_df[, setdiff(colnames(prot_df), colnames(prot_df)[gcol]), drop = FALSE])
    rownames(prot_mat) <- prot_genes
    storage.mode(prot_mat) <- "numeric"

    # Map site back to gene: directly use .make_site_id gene input (df$GENE_SYMBOL is aligned)
    # Recalculate here, only taking gene
    gene_of_site <- as.character(df$GENE_SYMBOL)
    gene_of_site <- gene_of_site[keep]

    # For each site, subtract corresponding gene's protein vector (based on common samples)
    common <- intersect(colnames(M_site), colnames(prot_mat))
    if (length(common) >= 2) {
      M_site <- M_site[, common, drop = FALSE]
      prot_mat <- prot_mat[, common, drop = FALSE]

      # Vectorized subtract: use split indexing to avoid row-by-row loop
      idx <- match(gene_of_site, rownames(prot_mat))
      hasp <- !is.na(idx)
      if (any(hasp)) {
        M_site[hasp, ] <- M_site[hasp, , drop = FALSE] - prot_mat[idx[hasp], , drop = FALSE]
      }
    } else {
      log_msg("(Warning) phospho and protein sample intersection < 2; skipping protein adjustment")
    }
  }

  # Return "site-level" matrix (downstream imputation / batch adjustment handled consistently with gene-level)
  M_site
}
# --- [END NEW] ---


## =========================================================
## [NEW | DPS] Differential Phosphorylation Site analysis (limma)
## - predictors: CSN subunits, CSN_SCORE, RESIDUAL_<SU> (consistent with template)
## - strata: ALL / TP53_mutant / TP53_wild_type (consistent with template)
## - version: BatchAdj (controlled by options(csn.run_passes))
## - No imputation; limma handles NA per site
## - Output:
##   With COMBO_MODE -> ./phospho_<combo>_ DPS/<dataset>/<STRATUM>/DPS/<version>/
##   Without COMBO_MODE -> <dataset>/phospho_DPS_results_TP53/<STRATUM>/DPS/<version>/
## =========================================================

.run_dps_fit <- function(M, pred, version, ds_id, out_dir,
                         batch_all, purity_all, sa_all, mat_for_sv = NULL,
                         extra_covars = NULL) { # <--- Added extra_covars argument
  ## 1) Align samples and predictors
  sample_order <- intersect(colnames(M), names(pred))
  pred <- pred[sample_order]
  M <- as.matrix(M[, sample_order, drop = FALSE])

  ## 2) Assemble design table from template (using existing .build_design_for_version)
  DF <- .build_design_for_version(
    version = version,
    pred_centered = pred, # Template already centered: use pred_centered to indicate
    sample_order = sample_order,
    batch_all = batch_all,
    purity_all = purity_all,
    sa_all = sa_all,
    USE_AGE_MISSING_INDICATOR = isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE)),
    ds_id = ds_id,
    mat_for_sv = mat_for_sv
  )

  ## [NEW] Merge extra_covars if provided (e.g. CSN_SCORE for one-step analysis)
  if (!is.null(extra_covars)) {
    extra_df <- .align_to_samples(extra_covars, sample_order, what = "extra_covars")
    extra_df <- as.data.frame(extra_df, check.names = FALSE)
    # Combine (be careful with name collisions, though usually distinctive)
    DF <- cbind(DF, extra_df)
  }

  ## 3) Select only template-allowed RHS variables; rename pred_centered to predictor as main coefficient
  # [Updated] Include columns from extra_covars in RHS
  rhs_candidates <- c("sex", "age", "age_missing", "age_z_imputed", "batch", paste0("SV", 1:5))
  if (!is.null(extra_covars)) rhs_candidates <- c(rhs_candidates, colnames(extra_df))

  rhs <- intersect(rhs_candidates, colnames(DF))

  DF2 <- data.frame(
    predictor = DF$pred_centered,
    DF[, rhs, drop = FALSE],
    row.names = rownames(DF),
    check.names = FALSE
  )

  ## 4) Drop constant/single-level/zero-variance columns (always keep predictor); also drop NA rows and align M
  keep_cols <- rep(TRUE, ncol(DF2))
  names(keep_cols) <- colnames(DF2)
  if (ncol(DF2) > 1L) {
    keep_cols <- vapply(DF2, function(z) {
      if (is.factor(z)) {
        nz <- droplevels(z)
        (nlevels(nz) >= 2) && (sum(!is.na(nz)) >= 2)
      } else {
        vz <- suppressWarnings(stats::var(as.numeric(z), na.rm = TRUE))
        is.finite(vz) && vz > 0 && sum(is.finite(z)) >= 2
      }
    }, logical(1))
    if ("predictor" %in% names(keep_cols)) keep_cols["predictor"] <- TRUE
    # [Optional] Force keep extra_covars if they are meaningful?
    # Usually better to let variance check decide. If CSN_SCORE is constant (impossible), it should be dropped.
  }
  DF2 <- DF2[, keep_cols, drop = FALSE]

  row_keep <- stats::complete.cases(DF2)
  if (sum(row_keep) < 3L) {
    log_msg("[DPS|%s|%s] Design has insufficient usable samples (<3) -> skip", ds_id, version)
    return(invisible(NULL))
  }
  DF2 <- DF2[row_keep, , drop = FALSE]
  M <- M[, rownames(DF2), drop = FALSE]
  for (cn in colnames(DF2)) if (is.factor(DF2[[cn]])) DF2[[cn]] <- droplevels(DF2[[cn]])

  stopifnot(identical(colnames(M), rownames(DF2)))

  ## 5) Build design matrix (with intercept): ~ predictor + RHS
  make_X <- function(df, rhs_vars) {
    if (length(rhs_vars)) {
      stats::model.matrix(stats::as.formula(paste("~ predictor +", paste(rhs_vars, collapse = " + "))), data = df)
    } else {
      stats::model.matrix(~predictor, data = df)
    }
  }
  rhs_use <- setdiff(colnames(DF2), "predictor")
  des <- make_X(DF2, rhs_use)

  ## 6) Saturation protection: if residual df <= 0, sequentially remove covariates until df_res > 0
  rank_d <- qr(des)$rank
  res_df <- nrow(des) - rank_d
  if (res_df <= 0L) {
    # [Updated] Drop order: SVs first, then batch/clinical. Keep extra_covars (biology) as long as possible.
    drop_candidates <- intersect(
      c(paste0("SV", 5:1), "batch", "sex", "age_missing", "age_z_imputed", "age"),
      rhs_use
    )
    # If we need to drop extra_covars (e.g. CSN_SCORE) to run at all, that might be needed too, but put them last.
    if (!is.null(extra_covars)) drop_candidates <- c(drop_candidates, colnames(extra_df))

    for (cv in drop_candidates) {
      rhs_use <- setdiff(rhs_use, cv)
      des <- make_X(DF2, rhs_use)
      rank_d <- qr(des)$rank
      res_df <- nrow(des) - rank_d
      log_msg("[DPS|%s|%s] Saturation protection: remove covariate=%s -> df_res=%d", ds_id, version, cv, res_df)
      if (res_df > 0L) break
    }
  }
  if (res_df <= 0L) {
    log_msg("[DPS|%s|%s] No residual df after removing all optional covariates -> skip", ds_id, version)
    return(invisible(NULL))
  }
  if (nrow(des) != ncol(M)) {
    stop(sprintf(
      "[DPS|%s|%s] Design rows (%d) != matrix columns (%d): sample alignment failed",
      ds_id, version, nrow(des), ncol(M)
    ))
  }

  ## 7) Fit and output (using predictor coefficient as DPS statistic)
  fit <- limma::eBayes(limma::lmFit(M, des))
  j <- match("predictor", colnames(des))
  if (is.na(j)) stop("[DPS] Design matrix coefficient not found: predictor")
  tt <- limma::topTable(fit, coef = j, number = Inf, sort.by = "none")

  # Standardize output: site ID and observation count
  if ("ID" %in% colnames(tt)) {
    # If limma has the actual site ID (e.g. "AHNAK_S18") in tt$ID
    tt$site_id <- tt$ID
  } else {
    # Otherwise use rownames as site ID
    tt$site_id <- rownames(tt)
  }

  tt$n_obs <- rowSums(is.finite(M))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(tt, file.path(out_dir, "DPS_results.csv"))
  invisible(tt)
}


run_dps_stratum <- function(ds_id, ds_dir, M_site_full, sample_keep, out_root, prot0_full) {
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  sam <- intersect(colnames(M_site_full), sample_keep)
  if (length(sam) < (2L * min_per_group)) {
    log_msg("  [DPS] %s insufficient samples -> skip", basename(out_root))
    return(invisible(NULL))
  }
  M0 <- M_site_full[, sam, drop = FALSE]
  prot0 <- prot0_full[, sam, drop = FALSE]

  present_sub <- intersect(csn_subunits, rownames(prot0))
  if (!length(present_sub)) {
    log_msg("  [DPS] No CSN subunit -> skip")
    return(invisible(NULL))
  }

  csn_score <- build_csn_score_safe(
    prot0,
    subunits = present_sub, combine_7AB = TRUE,
    min_members = 5L, pca_min_samples = 10L
  )

  # Covariates / batch consistent with template
  batch_all <- get_batch_factor_phospho(ds_dir, sam)
  batch_all <- if (!is.null(batch_all)) droplevels(batch_all$fac[sam]) else NULL
  purity_all <- get_purity_covariate(ds_id, ds_dir, sam)
  sa_all <- get_sex_age_covariates(ds_dir, sam)

  run_passes <- getOption("csn.run_passes", c("BatchAdj"))
  for (version in run_passes) {
    out_ver <- file.path(out_root, version)
    dir.create(out_ver, recursive = TRUE, showWarnings = FALSE)

    ## (1) subunits
    for (su in present_sub) {
      log_msg("  [DPS|%s|%s] predictor=%s", version, basename(out_root), su)
      .run_dps_fit(
        M = M0, pred = prot0[su, ], version = version, ds_id = ds_id,
        out_dir = file.path(out_ver, paste0("SU__", safe_fs_name(su))),
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
        mat_for_sv = M0
      )
    }

    ## (2) CSN_SCORE
    if (sum(is.finite(csn_score)) >= (2L * min_per_group)) {
      log_msg("  [DPS|%s|%s] predictor=CSN_SCORE", version, basename(out_root))
      .run_dps_fit(
        M = M0, pred = csn_score, version = version, ds_id = ds_id,
        out_dir = file.path(out_ver, "SU__CSN_SCORE"),
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
        mat_for_sv = M0
      )
    }

    ## (3) RESIDUAL_<SU> (Refactored to One-Step Analysis)
    ## Old method: Calculate residuals -> run limma (Two-Step, "Generated Regressor Problem")
    ## New method: Run limma with ~ Subunit + CSN_SCORE (One-Step Multivariable Regression)
    ## This estimates the effect of Subunit *adjusting for* CSN_SCORE directly.
    for (su in present_sub) {
      if (!sum(is.finite(csn_score))) next

      # Use original subunit expression as predictor
      # But Pass CSN_SCORE as an EXTRA COVARIATE to .run_dps_fit
      log_msg("  [DPS|%s|%s] predictor=RESIDUAL_%s (One-Step: ~ %s + CSN_SCORE)", version, basename(out_root), su, su)

      .run_dps_fit(
        M = M0, pred = prot0[su, ], version = version, ds_id = ds_id,
        out_dir = file.path(out_ver, paste0("SU__RESIDUAL_", safe_fs_name(su))),
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
        mat_for_sv = M0,
        extra_covars = data.frame(row.names = colnames(prot0), CSN_SCORE = csn_score) # <--- Key change: Add CSN_SCORE here
      )
    }

    ## (4) [NEW | DPS summary_wide] Summarize all predictor results within this stratum x version into a wide table
    {
      ## out_ver is currently the output folder for this version, e.g.:
      ##   <dataset>/<STRATUM>/DPS/<version>/
      ## It contains many subdirectories: SU__GPS1, SU__COPS2, SU__CSN_SCORE, SU__RESIDUAL_GPS1, ...
      ## We will extract logFC / adj.P.Val for each predictor, wide-merge into one table
      ## And output to out_ver/DPS_summary_wide.csv

      pred_dirs <- list.dirs(out_ver, full.names = TRUE, recursive = FALSE)
      pred_dirs <- pred_dirs[grepl("^SU__", basename(pred_dirs))]

      wide_list <- list()

      for (pd in pred_dirs) {
        pred_label <- sub("^SU__", "", basename(pd)) # e.g. "GPS1", "CSN_SCORE", "RESIDUAL_GPS1"
        f_csv <- file.path(pd, "DPS_results.csv")
        if (!file.exists(f_csv)) next

        tt_local <- try(data.table::fread(f_csv), silent = TRUE)
        if (inherits(tt_local, "try-error") || !nrow(tt_local)) next

        ## Get site name. v164 .run_dps_fit guarantees site_id presence.
        ## But for compatibility with old output, we still keep the fallback.
        if ("site_id" %in% colnames(tt_local)) {
          site_vec <- tt_local$site_id
        } else if ("ID" %in% colnames(tt_local)) {
          site_vec <- tt_local$ID
        } else {
          site_vec <- rownames(tt_local)
        }

        tmp <- data.frame(
          site_id = site_vec,
          logFC = tt_local$logFC,
          adj.P.Val = tt_local$adj.P.Val,
          stringsAsFactors = FALSE
        )

        ## If same site_id appears multiple times for this predictor (e.g. isoform/duplicate rownames), keep only first row to avoid Cartesian explosion during merge
        tmp <- tmp[!duplicated(tmp$site_id), , drop = FALSE]

        ## Rename to <predictor>.logFC / <predictor>.adj.P.Val
        colnames(tmp)[colnames(tmp) == "logFC"] <- paste0(pred_label, ".logFC")
        colnames(tmp)[colnames(tmp) == "adj.P.Val"] <- paste0(pred_label, ".adj.P.Val")

        wide_list[[pred_label]] <- tmp
      }

      if (length(wide_list)) {
        summary_wide <- Reduce(
          function(x, y) merge(x, y, by = "site_id", all = TRUE),
          wide_list
        )
        data.table::fwrite(summary_wide, file.path(out_ver, "DPS_summary_wide.csv"))
      }
    }
  }
  invisible(NULL)
}

## ---- [NEW | DPS main program] Process datasets sequentially (run DPS only, no PTM-SEA) ----
if (isTRUE(RUN_DPS_ONLY)) {
  # [FIX] Directly use dataset_dirs_run already built in main workflow (including combo and start_from rules)
  dataset_dirs_run_dps <- get("dataset_dirs_run", inherits = TRUE)
  if (is.null(dataset_dirs_run_dps) || !length(dataset_dirs_run_dps)) {
    stop("dataset_dirs_run is not defined or empty; please build dataset and combo selection in main workflow first.")
  }
  log_msg("[DPS] Will process datasets: %s", paste(names(dataset_dirs_run_dps), collapse = ","))

  for (ds in names(dataset_dirs_run_dps)) {
    ds_dir <- dataset_dirs_run_dps[[ds]]
    ds_id <- ds

    # Read site-level matrix (same as template)
    M_site_full <- try(load_phosphosite_matrix_from_dataset_dir(ds_dir), silent = TRUE)
    if (inherits(M_site_full, "try-error") || is.null(dim(M_site_full))) {
      log_msg("[DPS] %s unable to read site matrix -> skip", ds)
      next
    }

    # Read protein matrix (for building predictors)
    prot0_full <- try(load_matrix_from_dataset_dir(ds_dir), silent = TRUE)
    if (inherits(prot0_full, "try-error") || is.null(dim(prot0_full))) {
      log_msg("[DPS] %s no protein matrix -> skip", ds)
      next
    }

    # Consistent with template: LUSC special case (incompatible column structure) skip directly
    if (identical(ds, "lusc_cptac_2021")) {
      log_msg("[DPS] LUSC special case: skip site-level")
      next
    }

    # Stratified samples
    all_ids <- colnames(M_site_full)
    tp53_all <- get_tp53_status(ds_dir, all_ids)
    keep_all <- all_ids
    keep_mt <- names(tp53_all)[tp53_all == "TP53_mutant"]
    keep_wt <- names(tp53_all)[tp53_all == "TP53_wild_type"]

    # Determine output root directory (exactly following your naming: phospho_<combo>_ DPS)
    base_root <- if (!is.null(COMBO_PREFIX_DPS)) {
      file.path(COMBO_PREFIX_DPS, ds)
    } else {
      file.path(ds_dir, "phospho_DPS_results_TP53")
    }

    # Execute each stratum (version controlled by options(csn.run_passes), combo_2 defaults to BatchAdj only)
    run_dps_stratum(ds_id, ds_dir, M_site_full, keep_all, file.path(base_root, "ALL", "DPS"), prot0_full)
    if (length(keep_mt) >= (2L * min_per_group)) {
      run_dps_stratum(ds_id, ds_dir, M_site_full, keep_mt, file.path(base_root, "TP53_mutant", "DPS"), prot0_full)
    }
    if (length(keep_wt) >= (2L * min_per_group)) {
      run_dps_stratum(ds_id, ds_dir, M_site_full, keep_wt, file.path(base_root, "TP53_wild_type", "DPS"), prot0_full)
    }

    log_msg("[DPS] Completed dataset: %s -> %s", ds, base_root)
  }
}


## =========================================================
## =========================================================
## [DPS-meta] Generate cross-dataset Stouffer combined results (DPS)
## Condition: the above for-loop (DPS) has completed and
##       all dataset / stratum / version / SU__*/DPS_results.csv
##       have been written to disk
##
## Output:
##   <OUT_ROOT>/<STRATUM>/<VERSION>/<SUBUNIT>/DPS_meta_stouffer.csv
##   Where:
##     - site_id              (equivalent to pathway in PTM-SEA)
##     - logFC_mean           (equivalent to NES mean direction/magnitude in PTM-SEA)
##     - Z                    (Stouffer combined Z, direction uses logFC sign)
##     - p_meta, padj_meta    (BH-adjusted meta FDR)
##     - n_ds                 (how many datasets contribute to this site)
## =========================================================
dps_meta_stouffer <- function(
  dataset_dirs,
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  versions = c("BatchAdj"),
  out_root = if (is.null(COMBO_PREFIX_DPS)) {
    "phospho_DPS_meta"
  } else {
    file.path(COMBO_PREFIX_DPS, "phospho_DPS_meta")
  }
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install data.table first")
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("stats package required")
  }

  .safe_read_dps <- function(fp) {
    if (!file.exists(fp)) {
      return(NULL)
    }
    dt <- tryCatch(
      data.table::fread(fp, na.strings = c("NA", "NaN", "")),
      error = function(e) NULL
    )
    if (is.null(dt) || !nrow(dt)) {
      return(NULL)
    }

    ## Old version might put site ID in ID instead of site_id
    if (!"site_id" %in% names(dt) && "ID" %in% names(dt)) {
      dt$site_id <- dt$ID
    }
    dt
  }

  .combine_one_predictor <- function(tbl_list) {
    ## tbl_list: named list(dataset_id -> data.table with site_id, logFC, P.Value, adj.P.Val)
    keep <- tbl_list[!vapply(tbl_list, is.null, logical(1))]
    if (!length(keep)) {
      return(NULL)
    }

    long <- data.table::rbindlist(
      lapply(names(keep), function(ds_nm) {
        x <- keep[[ds_nm]]
        if (!all(c("site_id", "logFC", "P.Value", "adj.P.Val") %in% names(x))) {
          return(NULL)
        }
        data.table::data.table(
          dataset = ds_nm,
          site_id = x$site_id,
          logFC   = as.numeric(x$logFC),
          pval    = as.numeric(x$P.Value),
          padj    = as.numeric(x$adj.P.Val)
        )
      }),
      fill = TRUE
    )
    if (is.null(long) || !nrow(long)) {
      return(NULL)
    }

    ## Direction = sign of logFC, magnitude = two-tailed p-value converted to z
    long[, z_i := {
      # [Safeguard] Set p-value lower limit to avoid p=0 -> Z=Inf -> Meta=NaN
      # Use .Machine$double.xmin (approx 2.2e-308) to preserve ranking of extremely significant hits
      safe_p <- pmax(pval, .Machine$double.xmin)

      dirn <- sign(logFC)
      zz <- stats::qnorm(1 - safe_p / 2)
      dirn * zz
    }]

    ## Aggregate by site_id to get meta
    summ <- long[
      , .(
        n_ds       = sum(is.finite(z_i)),
        Z          = if (sum(is.finite(z_i)) > 0) sum(z_i, na.rm = TRUE) / sqrt(sum(is.finite(z_i))) else NA_real_,
        logFC_mean = mean(logFC, na.rm = TRUE),
        pval_min   = suppressWarnings(min(pval, na.rm = TRUE)),
        padj_min   = suppressWarnings(min(padj, na.rm = TRUE))
      ),
      by = .(site_id)
    ]

    ## Z -> two-tailed p, then BH to get meta FDR
    summ[, p_meta := 2 * stats::pnorm(-abs(Z))]
    summ[, padj_meta := stats::p.adjust(p_meta, method = "BH")]

    ## Sorting criteria:
    ## 1. padj_meta ascending (more significant first)
    ## 2. p_meta ascending
    ## 3. |Z| descending (direction strength)
    ## 4. site_id alphabetical
    summ[, absZ := abs(Z)]
    data.table::setorder(summ, padj_meta, p_meta, -absZ, site_id)
    summ[, absZ := NULL]

    summ[]
  }

  for (st in strata) {
    for (ver in versions) {
      ## ------------------------------------------------------------------
      ## Collect subunit/predictor directories (SU__*) from all datasets, take union
      ## This way we won't only look at first dataset and miss like COPS3 / RESIDUAL_COPS3
      ## ------------------------------------------------------------------
      su_tags_all <- character(0)

      for (ds_nm in names(dataset_dirs)) {
        probe_dir <- file.path(dataset_dirs[[ds_nm]], st, "DPS", ver)
        if (!dir.exists(probe_dir)) next

        # Here using full.names = FALSE, directly get directory names, e.g. "SU__COPS1"
        this_su_dirs <- list.dirs(probe_dir, full.names = FALSE, recursive = FALSE)
        this_su_tags <- this_su_dirs[grepl("^SU__", this_su_dirs)]
        if (length(this_su_tags)) {
          su_tags_all <- union(su_tags_all, this_su_tags)
        }
      }

      ## If this stratum/version has no SU__* in any dataset, skip
      if (!length(su_tags_all)) next

      ## ------------------------------------------------------------------
      ## For each predictor/subunit (e.g. SU__COPS1, SU__COPS3, SU__RESIDUAL_COPS3, ...)
      ## Collect results from all datasets with same name and do Stouffer combination
      ## ------------------------------------------------------------------
      for (su_tag in su_tags_all) {
        su_name <- sub("^SU__", "", su_tag)

        ## Collect DPS_results.csv from all datasets for the same predictor
        per_ds <- lapply(names(dataset_dirs), function(ds_nm) {
          fp <- file.path(
            dataset_dirs[[ds_nm]],
            st, "DPS", ver,
            su_tag,
            "DPS_results.csv"
          )
          .safe_read_dps(fp)
        })
        names(per_ds) <- names(dataset_dirs)

        ## Perform Stouffer combination (including Z, p_meta, padj_meta, logFC_mean, n_ds ...)
        comb <- .combine_one_predictor(per_ds)
        if (is.null(comb) || !nrow(comb)) next

        ## Write out meta results:
        ## <out_root>/<STRATUM>/<VERSION>/<SUBUNIT>/DPS_meta_stouffer.csv
        out_dir <- file.path(out_root, st, ver, su_name)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        out_fp <- file.path(out_dir, "DPS_meta_stouffer.csv")
        data.table::fwrite(comb, out_fp)

        if (exists("log_msg", mode = "function")) {
          try(
            log_msg("  [DPS meta] stratum={st} ver={ver} su={su_name} -> {out_fp}"),
            silent = TRUE
          )
        }
      } # end for each su_tag
    } # end for ver
  } # end for st

  invisible(TRUE)
}

## =========================================================
## [DPS-meta] Execute cross-dataset Stouffer combination
## =========================================================

## Build mapping table of actual DPS result root directories:
##   - combo mode (e.g. COMBO_MODE == "combo_2"): phospho_<combo>_ DPS/<dataset>
##   - non-combo mode: <dataset>/phospho_DPS_results_TP53
dataset_dirs_run_dps_out <- setNames(
  if (!is.null(COMBO_PREFIX_DPS)) {
    file.path(COMBO_PREFIX_DPS, names(dataset_dirs_run))
  } else {
    file.path(dataset_dirs_run, "phospho_DPS_results_TP53")
  },
  names(dataset_dirs_run)
)

dps_meta_stouffer(
  dataset_dirs = dataset_dirs_run_dps_out, # This is each dataset's DPS base_root
  strata = strata, # c("ALL","TP53_mutant","TP53_wild_type")
  versions = c("BatchAdj"),
  out_root = if (is.null(COMBO_PREFIX_DPS)) {
    "phospho_DPS_meta"
  } else {
    file.path(COMBO_PREFIX_DPS, "phospho_DPS_meta")
  }
)


## =========================================================
## =========================================================
## [DPS-meta summary] Summarize DPS Stouffer meta across subunits
##
##   - Input:
##       <meta_root>/<STRATUM>/BatchAdj/<SUBUNIT>/DPS_meta_stouffer.csv
##     Where each file comes from dps_meta_stouffer(), already containing
##     site-level (site_id) cross-dataset Stouffer combined results
##
##   - Output:
##       <meta_root>/summary/<STRATUM>/BatchAdj/DPS/DPS_meta_stouffer/
##          ├─ Summary_DPS_DPS_meta_stouffer_ALL.csv
##          ├─ Summary_DPS_DPS_meta_stouffer_padjLT0.05.csv
##          ├─ Summary_DPS_DPS_meta_stouffer_padjLT0.25.csv
##          └─ Summary_DPS_DPS_meta_stouffer.xlsx
##
##   - Consolidation method:
##       1. Read (Z, padj_meta) for each subunit/predictor
##          -> Convert to columns Z_<SUBUNIT>, padj_meta_<SUBUNIT>
##       2. full_join by site_id (site_id is equivalent to pathway in PTM-SEA)
##       3. Calculate how many subunits are significant for each site, how many positive/negative
##          (padj_meta < 0.05 / 0.25, direction based on Z_*)
## =========================================================

.read_dps_meta_table <- function(meta_root, stratum, version, subunit,
                                 stat_tag = "DPS_meta_stouffer") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")

  fp <- file.path(meta_root, stratum, version, subunit, "DPS_meta_stouffer.csv")
  if (!file.exists(fp)) {
    return(NULL)
  }

  dt <- tryCatch(
    data.table::fread(fp, na.strings = c("NA", "NaN", "")),
    error = function(e) NULL
  )
  if (is.null(dt)) {
    return(NULL)
  }

  need <- c("site_id", "Z", "padj_meta")
  ## Fault-tolerant column name mapping
  nms <- tolower(names(dt))
  if (!"site_id" %in% names(dt) && "site_id" %in% nms) names(dt)[match("site_id", nms)] <- "site_id"
  if (!"Z" %in% names(dt) && "z" %in% nms) names(dt)[match("z", nms)] <- "Z"
  if (!"padj_meta" %in% names(dt) && "padj_meta" %in% nms) names(dt)[match("padj_meta", nms)] <- "padj_meta"
  if (!all(need %in% names(dt))) {
    return(NULL)
  }

  out <- as.data.frame(dt[, ..need])
  names(out)[names(out) == "Z"] <- paste0("Z_", subunit)
  names(out)[names(out) == "padj_meta"] <- paste0("padj_meta_", subunit)
  out
}

.merge_subunit_tables_dps <- function(tbl_list) {
  keep <- tbl_list[!vapply(tbl_list, is.null, logical(1))]
  if (!length(keep)) {
    return(NULL)
  }
  out <- Reduce(function(x, y) dplyr::full_join(x, y, by = "site_id"), keep)
  num_cols <- setdiff(names(out), "site_id")
  out[num_cols] <- lapply(out[num_cols], function(z) suppressWarnings(as.numeric(z)))
  out
}

.add_sig_counts_meta_dps <- function(df, alphas = c(0.05, 0.25)) {
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

    ## How many subunits reached padj_meta < a for this site (regardless of direction)
    sig_mat <- sapply(padj_cols, function(p) {
      pv <- df[[p]]
      as.integer(is.finite(pv) & pv < a)
    })
    if (is.null(dim(sig_mat))) sig_mat <- matrix(sig_mat, ncol = 1)
    df[[sprintf("sig_n_padj_meta_%s", a_tag)]] <- rowSums(sig_mat, na.rm = TRUE)

    ## Positive significant: padj_meta < a and Z > 0
    pos_mat <- mapply(function(p, zc) {
      pv <- df[[p]]
      zv <- df[[zc]]
      as.integer(is.finite(pv) & pv < a & is.finite(zv) & zv > 0)
    }, padj_cols, z_cols, SIMPLIFY = TRUE)
    if (is.null(dim(pos_mat))) pos_mat <- matrix(pos_mat, ncol = 1)
    df[[sprintf("pos_n_padj_meta_%s", a_tag)]] <- rowSums(pos_mat, na.rm = TRUE)

    ## Negative significant: padj_meta < a and Z < 0
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

.write_summary_outputs_meta_csv_dps <- function(df, out_dir,
                                                group_name = "DPS",
                                                stat_tag = "DPS_meta_stouffer") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install openxlsx first")

  if (is.null(df) || !nrow(df)) {
    if (exists("log_msg", mode = "function")) {
      try(log_msg("  [Skip output] {group_name} | {stat_tag} no available results"), silent = TRUE)
    }
    return(invisible(NULL))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(out_dir, paste0("Summary_", safe_fs_name(group_name), "_", stat_tag))

  ## ---- Sorting: prioritize by number of subunits with padj_meta<0.05, then positive/negative, then site_id
  ord_keys <- c("sig_n_padj_meta_0_05", "pos_n_padj_meta_0_05", "neg_n_padj_meta_0_05")
  ord_keys <- intersect(ord_keys, names(df))
  if (length(ord_keys)) {
    df <- df |>
      dplyr::arrange(
        dplyr::desc(.data[[ord_keys[1]]]),
        dplyr::desc(ifelse(length(ord_keys) > 1, .data[[ord_keys[2]]], 0)),
        dplyr::desc(ifelse(length(ord_keys) > 2, .data[[ord_keys[3]]], 0)),
        .data[["site_id"]]
      )
  } else {
    df <- df |> dplyr::arrange(.data[["site_id"]])
  }

  ## ---- Output CSV (ALL / padjLT0.05 / padjLT0.25)
  data.table::fwrite(df, paste0(base, "_ALL.csv"))

  df005 <- NULL
  df025 <- NULL
  if ("sig_n_padj_meta_0_05" %in% names(df)) {
    df005 <- df |> dplyr::filter(.data[["sig_n_padj_meta_0_05"]] > 0)
    dir.create(dirname(paste0(base, "_padjLT0.05.csv")), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(df005, paste0(base, "_padjLT0.05.csv"))
  }
  if ("sig_n_padj_meta_0_25" %in% names(df)) {
    df025 <- df |> dplyr::filter(.data[["sig_n_padj_meta_0_25"]] > 0)
    dir.create(dirname(paste0(base, "_padjLT0.25.csv")), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(df025, paste0(base, "_padjLT0.25.csv"))
  }

  ## ---- XLSX (three sheets: ALL / padjLT0.05 / padjLT0.25)
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "ALL")
  openxlsx::writeData(wb, "ALL", df)
  openxlsx::freezePane(wb, "ALL", firstRow = TRUE)
  openxlsx::setColWidths(wb, "ALL", cols = 1:ncol(df), widths = "auto")

  openxlsx::addWorksheet(wb, "padjLT0.05")
  if (!is.null(df005)) {
    openxlsx::writeData(wb, "padjLT0.05", df005)
    openxlsx::freezePane(wb, "padjLT0.05", firstRow = TRUE)
    openxlsx::setColWidths(wb, "padjLT0.05", cols = 1:ncol(df), widths = "auto")
  }

  openxlsx::addWorksheet(wb, "padjLT0.25")
  if (!is.null(df025)) {
    openxlsx::writeData(wb, "padjLT0.25", df025)
    openxlsx::freezePane(wb, "padjLT0.25", firstRow = TRUE)
    openxlsx::setColWidths(wb, "padjLT0.25", cols = 1:ncol(df), widths = "auto")
  }

  openxlsx::saveWorkbook(wb, paste0(base, ".xlsx"), overwrite = TRUE)

  if (exists("log_msg", mode = "function")) {
    try(log_msg("  [Output complete] {group_name} | {stat_tag} -> {dirname(base)} (.csv + .xlsx)"),
      silent = TRUE
    )
  }

  invisible(NULL)
}

summarize_meta_fdr_across_subunits <- function(
  meta_root = if (is.null(COMBO_PREFIX_DPS)) {
    "phospho_DPS_meta"
  } else {
    file.path(COMBO_PREFIX_DPS, "phospho_DPS_meta")
  },
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  versions = c("BatchAdj"),
  stat_tag = "DPS_meta_stouffer"
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")

  for (st in strata) {
    for (ver in versions) {
      base_dir <- file.path(meta_root, st, ver)
      if (!dir.exists(base_dir)) next

      ## Auto-detect subunit/predictor names (directory names)
      subs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
      subs <- subs[nzchar(subs)]
      if (!length(subs)) next

      if (exists("log_msg", mode = "function")) {
        try(log_msg("== [DPS-meta-summary] stratum=%s | version=%s ==", st, ver), silent = TRUE)
      }

      ## Collect meta results (Z, padj_meta) for each subunit
      lst <- setNames(vector("list", length(subs)), subs)
      for (su in subs) {
        lst[[su]] <- .read_dps_meta_table(meta_root, st, ver, su, stat_tag)
      }

      ## wide: site_id as rows, Z_<SU> / padj_meta_<SU> as columns
      wide <- .merge_subunit_tables_dps(lst)

      ## Add significance statistics (padj_meta<0.05 / 0.25, positive/negative)
      wide <- .add_sig_counts_meta_dps(wide, alphas = c(0.05, 0.25))

      ## Write out summary
      out_dir <- file.path(meta_root, "summary", st, ver, "DPS", stat_tag)
      .write_summary_outputs_meta_csv_dps(
        wide,
        out_dir,
        group_name = "DPS",
        stat_tag   = stat_tag
      )
    }
  }

  invisible(TRUE)
}

## ---- One-click execution (no re-run DPS; only read existing DPS_meta_stouffer.csv then summarize) ----
posthoc_summary_meta_fdr <- function() {
  summarize_meta_fdr_across_subunits(
    meta_root = if (is.null(COMBO_PREFIX_DPS)) {
      "phospho_DPS_meta"
    } else {
      file.path(COMBO_PREFIX_DPS, "phospho_DPS_meta")
    },
    strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
    versions = c("BatchAdj"),
    stat_tag = "DPS_meta_stouffer"
  )
  invisible(TRUE)
}

## Run directly once (can be executed repeatedly; will overwrite same-name csv / xlsx)
posthoc_summary_meta_fdr()




## ---- [NEW | DPS meta-FDR heatmap (ggplot2)] ----
## Purpose:
## 1. Read your DPS pipeline's meta Stouffer FDR summary file
##    (phosphoproteomic_DPS/phospho_DPS_meta/summary/<STRATUM>/BatchAdj/DPS/DPS_meta_stouffer/
##     Summary_DPS_DPS_meta_stouffer_ALL.csv)
##    Three strata: ALL / TP53_mutant / TP53_wild_type
##
## 2. Select one or more phosphosite sets (pathways) from PTMsigDB, extract all site_ids from them
##
## 3. Use ALL stratum's Z_CSN_SCORE to determine:
##    - y-axis order (larger Z_CSN_SCORE on top)
##    - Which sites to plot (default takes top top_n_pos with largest Z_CSN_SCORE + top_n_neg with smallest)
##    Then TP53_mutant, TP53_wild_type use the same site_ids and same order for plotting,
##    for easier comparison.
##
## 4. Plot heatmap (ggplot2::geom_tile):
##    x-axis = predictor (COPS subunits, CSN_SCORE, GPS1, RESIDUAL_*...)
##    y-axis = phosphosite (site_id)
##    fill = Z_<predictor>
##    Overlay significance marker: padj_meta_<predictor> < 0.05
##    Color scale positive/negative colors can be specified (col_up for Z>0, col_down for Z<0)
##
## 5. Output simultaneously:
##    - Image file (default only .tiff; if out_formats includes "png"/"pdf"/"jpg" will also output)
##    - Corresponding data file .csv (each row is a site_id, columns are all predictor's Z_* and padj_meta_*;
##      exported is the full set of target pathway(s), not just top/bottom shown in heatmap)
##
## Assumptions:
## - You already have global genesets_by_group_ptm, and
##   genesets_by_group_ptm[["PTMsigDB"]] is a list:
##   names = pathway names, value = phosphosite vectors for this pathway (like "SORBS1_S1")
##
## - COMBO_PREFIX_DPS (if combo mode, e.g. "phosphoproteomic_DPS") already exists,
##   path logic same as summarize_meta_fdr_across_subunits() / posthoc_summary_meta_fdr()
##
## - Currently meta summary is stored under BatchAdj version, filename is fixed as Summary_DPS_DPS_meta_stouffer_ALL.csv
##   Same filename in all three stratum (ALL / TP53_mutant / TP53_wild_type) paths
##
## Dependencies: data.table, ggplot2, dplyr, scales
## ------------------------------------------------------------

.read_dps_summary_csv <- function(
  stratum,
  base_root = if (is.null(COMBO_PREFIX_DPS)) {
    "phospho_DPS_meta"
  } else {
    file.path(COMBO_PREFIX_DPS, "phospho_DPS_meta")
  },
  version = "BatchAdj"
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
  fp <- file.path(
    base_root,
    "summary",
    stratum,
    version,
    "DPS",
    "DPS_meta_stouffer",
    "Summary_DPS_DPS_meta_stouffer_ALL.csv"
  )
  if (!file.exists(fp)) {
    stop("Cannot find DPS summary file: ", fp)
  }
  data.table::fread(fp, na.strings = c("NA", "NaN", "Inf", "-Inf", ""))
}

.get_sites_from_ptmsigdb <- function(
  target_pathways,
  ptm_sets = genesets_by_group_ptm[["PTMsigDB"]]
) {
  ## Map user-specified PTMsigDB pathway(s) to all phosphosites
  ## If a pathway name doesn't exist, it will be automatically ignored
  if (is.null(ptm_sets)) {
    stop("genesets_by_group_ptm[['PTMsigDB']] does not exist, please first build PTMsigDB phosphosite sets")
  }
  keep_list <- ptm_sets[target_pathways]
  keep_list <- keep_list[!vapply(keep_list, is.null, logical(1))]
  unique(unlist(keep_list, use.names = FALSE))
}

.build_dps_long_df_for_plot <- function(df_wide, keep_sites_ordered) {
  ## Convert wide table (site_id, Z_xxx, padj_meta_xxx, ...) to long format
  ## And build:
  ##   - y-axis order:
  ##       * Based on keep_sites_ordered sorted by ALL's Z_CSN_SCORE
  ##       * But we want largest Z_CSN_SCORE at top of heatmap
  ##         -> Because ggplot factor levels are drawn bottom to top, we need to reverse levels
  ##
  ##   - x-axis order & group spacing:
  ##       * Sorting order follows your site level PTM-SEA meta FDR heatmap:
  ##           Z_CSN_SCORE,
  ##           Z_GPS1,
  ##           Z_COPS2..Z_COPS9,
  ##           Z_RESIDUAL_GPS1,
  ##           Z_RESIDUAL_COPS2..Z_RESIDUAL_COPS9
  ##       * Only insert a gap between
  ##           (Z_CSN_SCORE <-> Z_GPS1)
  ##         and
  ##           (Z_COPS9 <-> Z_RESIDUAL_GPS1)
  ##         (this is the blank position in your site level PTM-SEA heatmap)
  ##
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr first")

  df_sub <- df_wide[df_wide$site_id %in% keep_sites_ordered, , drop = FALSE]

  ## Find all predictor names (Z_*)
  z_cols <- grep("^Z_", names(df_sub), value = TRUE)
  pred_names <- sub("^Z_", "", z_cols)

  ## long format
  long_lst <- lapply(pred_names, function(pr) {
    data.frame(
      site_id = df_sub$site_id,
      predictor = paste0("Z_", pr),
      Z = df_sub[[paste0("Z_", pr)]],
      padj_meta = df_sub[[paste0("padj_meta_", pr)]],
      stringsAsFactors = FALSE
    )
  })
  df_long <- do.call(rbind, long_lst)

  ## y-axis order (based on keep_sites_ordered; this vector is already sorted from small to large)
  df_long$site_id <- factor(
    df_long$site_id,
    levels  = keep_sites_ordered,
    ordered = TRUE
  )

  ## ---- x-axis order & group spacing (follows PTM-SEA meta FDR heatmap) ----
  canonical_order <- c(
    "Z_CSN_SCORE",
    "Z_GPS1",
    "Z_COPS2", "Z_COPS3", "Z_COPS4", "Z_COPS5", "Z_COPS6", "Z_COPS7A", "Z_COPS7B", "Z_COPS8", "Z_COPS9",
    "Z_RESIDUAL_GPS1",
    "Z_RESIDUAL_COPS2", "Z_RESIDUAL_COPS3", "Z_RESIDUAL_COPS4", "Z_RESIDUAL_COPS5", "Z_RESIDUAL_COPS6",
    "Z_RESIDUAL_COPS7A", "Z_RESIDUAL_COPS7B", "Z_RESIDUAL_COPS8", "Z_RESIDUAL_COPS9"
  )

  present_preds <- intersect(canonical_order, unique(df_long$predictor))
  leftover <- setdiff(unique(df_long$predictor), present_preds)
  ## leftover (theoretically almost never used) goes to the end as a fallback
  present_preds <- c(present_preds, sort(leftover))

  ## Gap strategy:
  ##   gap1: between Z_CSN_SCORE <-> Z_GPS1
  ##   gap2: between Z_COPS9 <-> Z_RESIDUAL_GPS1
  ##
  gap_after_preds <- c("Z_CSN_SCORE", "Z_COPS9")
  gap_size <- 0.4

  pos_map <- list()
  cur_pos <- 0
  for (p in present_preds) {
    cur_pos <- cur_pos + 1
    pos_map[[p]] <- cur_pos
    if (p %in% gap_after_preds) {
      cur_pos <- cur_pos + gap_size
    }
  }

  ## Write back corresponding continuous x coordinates
  df_long$xpos <- unname(
    vapply(df_long$predictor, function(p) pos_map[[p]], numeric(1))
  )

  list(
    df_long       = df_long,
    present_preds = present_preds,
    xpos_map      = pos_map
  )
}


.save_dps_heatmap_outputs <- function(
  p,
  df_export,
  out_dir,
  out_stub,
  width_cm,
  height_cm,
  dpi = 600,
  formats = c("tiff")
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2 first")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ## 1. Save image files (multiple formats)
  for (fmt in formats) {
    fp <- file.path(out_dir, paste0(out_stub, ".", fmt))
    if (fmt %in% c("tiff", "tif")) {
      ggplot2::ggsave(
        filename = fp, plot = p,
        width = width_cm / 2.54, height = height_cm / 2.54,
        dpi = dpi, compression = "lzw", units = "in"
      )
    } else {
      ggplot2::ggsave(
        filename = fp, plot = p,
        width = width_cm / 2.54, height = height_cm / 2.54,
        dpi = dpi, units = "in"
      )
    }
  }

  ## 2. Save heatmap data (rows = site_id, full set, not just top/bottom; columns = Z_* and padj_meta_*)
  data.table::fwrite(
    df_export,
    file.path(out_dir, paste0(out_stub, "_data.csv"))
  )

  invisible(TRUE)
}

plot_dps_meta_heatmap_for_pathways <- function(
  target_pathways,
  top_n_pos = 25, ## Take top N with largest Z_CSN_SCORE
  top_n_neg = 25, ## Take top N with smallest Z_CSN_SCORE
  version = "BatchAdj",
  base_root = if (is.null(COMBO_PREFIX_DPS)) {
    "phospho_DPS_meta"
  } else {
    file.path(COMBO_PREFIX_DPS, "phospho_DPS_meta")
  },
  col_up = "#67001F", ## Color for Z > 0 (upregulated)
  col_down = "#053061", ## Color for Z < 0 (downregulated)
  mid_col = "#FFFFFF",
  alpha_sig = 0.05, ## Threshold for padj_meta significance
  sig_dot_size = 0.6, ## [NEW] Significance dot size (applied to all three strata)
  out_formats = c("tiff"), ## Can also specify c("tiff","png","pdf","jpg")
  fig_width_in = NULL, ## [NEW] Specify output width (inches); NULL = auto-calculate based on data
  fig_height_in = NULL ## [NEW] Specify output height (inches); NULL = auto-calculate based on data
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2 first")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr first")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Please install scales first")

  ## --------------------------------------------------------
  ## (1) Read DPS meta summary (wide table) for three strata
  ## --------------------------------------------------------
  wide_ALL <- .read_dps_summary_csv("ALL", base_root, version)
  wide_MT <- .read_dps_summary_csv("TP53_mutant", base_root, version)
  wide_WT <- .read_dps_summary_csv("TP53_wild_type", base_root, version)

  ## --------------------------------------------------------
  ## (2) Get all phosphosites from user-specified pathway(s) in PTMsigDB
  ## --------------------------------------------------------
  all_sites_in_sets <- .get_sites_from_ptmsigdb(
    target_pathways,
    ptm_sets = genesets_by_group_ptm[["PTMsigDB"]]
  )

  ## --------------------------------------------------------
  ## (3) Use ALL stratum to determine:
  ##     - y-axis sorting (by Z_CSN_SCORE)
  ##     - Which sites to plot (top_n_pos + top_n_neg)
  ## --------------------------------------------------------
  wide_ALL_keep <- wide_ALL[wide_ALL$site_id %in% all_sites_in_sets, , drop = FALSE]

  if (!"Z_CSN_SCORE" %in% names(wide_ALL_keep)) {
    stop("ALL stratum summary missing Z_CSN_SCORE column, cannot sort by Z_CSN_SCORE")
  }

  ## Sort by Z_CSN_SCORE descending / ascending
  ord_pos <- wide_ALL_keep$site_id[order(-wide_ALL_keep$Z_CSN_SCORE)]
  ord_pos <- unique(ord_pos)

  ord_neg <- wide_ALL_keep$site_id[order(wide_ALL_keep$Z_CSN_SCORE)]
  ord_neg <- unique(ord_neg)

  sel_pos <- head(ord_pos, top_n_pos)
  sel_neg <- head(ord_neg, top_n_neg)

  ## y-axis target order:
  ##   - We want "larger Z_CSN_SCORE on top"
  ##   - ggplot draws y=factor with levels from bottom to top
  ##   - So we put most negative first, most positive last
  ##   - sel_neg starts from most negative, sel_pos is rev() of largest positive
  ##   -> Put sel_neg first, then rev(sel_pos) after
  site_order_raw <- c(sel_neg, rev(sel_pos))

  ## Remove duplicates (some sites may appear in both sel_neg and sel_pos)
  site_order <- site_order_raw[!duplicated(site_order_raw)]

  ## Only plot this set in heatmap
  heatmap_sites <- site_order

  ## --------------------------------------------------------
  ## (4) Prepare long format (ALL / TP53_mutant / TP53_wild_type)
  ##     But y-axis levels all use site_order
  ## --------------------------------------------------------
  prep_ALL <- .build_dps_long_df_for_plot(wide_ALL, site_order)
  prep_MT <- .build_dps_long_df_for_plot(wide_MT, site_order)
  prep_WT <- .build_dps_long_df_for_plot(wide_WT, site_order)

  ## --------------------------------------------------------
  ## (5) Plot single stratum heatmap (y uses site_order, only plot heatmap_sites)
  ## --------------------------------------------------------
  .make_stratum_plot <- function(prep_obj, stratum_label) {
    df_long <- prep_obj$df_long
    present_preds <- prep_obj$present_preds
    pos_map <- prep_obj$xpos_map

    ## Only plot heatmap_sites determined by top_n_pos / top_n_neg
    df_long <- df_long[df_long$site_id %in% heatmap_sites, , drop = FALSE]

    ## Symmetric color scale: positive red, negative blue, middle white
    lim <- max(abs(df_long$Z), na.rm = TRUE)
    if (!is.finite(lim) || lim <= 0) lim <- 1

    ## x-axis breaks/labels follow present_preds order, using gap-adjusted coordinates from pos_map
    x_breaks <- unname(vapply(present_preds, function(p) pos_map[[p]], numeric(1)))
    x_labs <- present_preds

    ggplot2::ggplot(
      df_long,
      ggplot2::aes(x = xpos, y = site_id, fill = Z)
    ) +
      ## Tile itself: no border, like your site level PTM-SEA heatmap
      ggplot2::geom_tile(
        width  = 0.9,
        height = 0.9,
        color  = NA
      ) +
      ## Significance dots: padj_meta < alpha_sig -> draw small black filled dots
      ## Dots should be small (consistent with your site level PTM-SEA heatmap), avoid covering main color blocks
      ggplot2::geom_point(
        data = df_long[df_long$padj_meta < alpha_sig &
          is.finite(df_long$padj_meta), , drop = FALSE],
        ggplot2::aes(x = xpos, y = site_id),
        inherit.aes = FALSE,
        shape = 16, # Filled circle
        size = sig_dot_size, # [NEW] Adjustable dot size (default 0.6)
        colour = "black"
      ) +
      ggplot2::scale_fill_gradient2(
        low       = col_down,
        mid       = mid_col,
        high      = col_up,
        midpoint  = 0,
        limits    = c(-lim, lim),
        oob       = scales::squish,
        name      = "Z"
      ) +
      ## y-axis keeps the levels we set in build_dps_long_df_for_plot()
      ggplot2::scale_y_discrete(drop = FALSE) +
      ## x-axis uses continuous coordinates + custom breaks (with gaps)
      ggplot2::scale_x_continuous(
        breaks   = x_breaks,
        labels   = x_labs,
        expand   = ggplot2::expansion(mult = c(0.02, 0.02)),
        position = "top"
      ) +
      ggplot2::labs(
        x = paste0("predictor (", stratum_label, ")"),
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        ## Reduce y-axis text size to avoid "overlapping" as you mentioned
        axis.text.y = ggplot2::element_text(size = 5),
        ## x-axis label font also reduced, maintaining vertical placement and left alignment
        axis.text.x = ggplot2::element_text(
          angle = 90,
          hjust = 0,
          vjust = 0.5,
          size  = 5
        ),
        ## margin condensed similar to your site level PTM-SEA layout, not too wide nor squished
        plot.margin = ggplot2::margin(5.5, 20, 5.5, 5.5, "pt"),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 7),
        legend.text = ggplot2::element_text(size = 6)
      )
  }

  p_ALL <- .make_stratum_plot(prep_ALL, "ALL")
  p_MT <- .make_stratum_plot(prep_MT, "TP53_mutant")
  p_WT <- .make_stratum_plot(prep_WT, "TP53_wild_type")

  ## --------------------------------------------------------
  ## (6) Export
  ##     - Image files (tiff / png / pdf / jpg)
  ##     - Data file (all site_in_sets, not just top/bottom)
  ## --------------------------------------------------------
  path_stub <- safe_fs_name(paste(target_pathways, collapse = "_"))

  out_stub_base <- paste0(
    "DPS_meta_heatmap_",
    path_stub,
    "_top", top_n_pos, "pos_bot", top_n_neg, "neg_",
    version
  )

  ## ALL
  out_dir_ALL <- file.path(
    base_root, "summary", "ALL", "BatchAdj", "DPS", "DPS_meta_stouffer"
  )
  .save_dps_heatmap_outputs(
    p           = p_ALL,
    df_export   = wide_ALL[wide_ALL$site_id %in% all_sites_in_sets, , drop = FALSE],
    out_dir     = out_dir_ALL,
    out_stub    = paste0(out_stub_base, "_ALL"),
    width_cm    = if (is.null(fig_width_in)) max(6, length(prep_ALL$present_preds) * 0.6) else fig_width_in * 2.54,
    height_cm   = if (is.null(fig_height_in)) max(6, length(site_order) * 0.25) else fig_height_in * 2.54,
    formats     = out_formats
  )

  ## TP53_mutant
  out_dir_MT <- file.path(
    base_root, "summary", "TP53_mutant", "BatchAdj", "DPS", "DPS_meta_stouffer"
  )
  .save_dps_heatmap_outputs(
    p           = p_MT,
    df_export   = wide_MT[wide_MT$site_id %in% all_sites_in_sets, , drop = FALSE],
    out_dir     = out_dir_MT,
    out_stub    = paste0(out_stub_base, "_TP53_mutant"),
    width_cm    = if (is.null(fig_width_in)) max(6, length(prep_MT$present_preds) * 0.6) else fig_width_in * 2.54,
    height_cm   = if (is.null(fig_height_in)) max(6, length(site_order) * 0.25) else fig_height_in * 2.54,
    formats     = out_formats
  )

  ## TP53_wild_type
  out_dir_WT <- file.path(
    base_root, "summary", "TP53_wild_type", "BatchAdj", "DPS", "DPS_meta_stouffer"
  )
  .save_dps_heatmap_outputs(
    p           = p_WT,
    df_export   = wide_WT[wide_WT$site_id %in% all_sites_in_sets, , drop = FALSE],
    out_dir     = out_dir_WT,
    out_stub    = paste0(out_stub_base, "_TP53_wild_type"),
    width_cm    = if (is.null(fig_width_in)) max(6, length(prep_WT$present_preds) * 0.6) else fig_width_in * 2.54,
    height_cm   = if (is.null(fig_height_in)) max(6, length(site_order) * 0.25) else fig_height_in * 2.54,
    formats     = out_formats
  )

  ## Invisibly return plots + order (convenient for interactive checking, doesn't affect file writing)
  invisible(list(
    ALL            = list(plot = p_ALL, order = site_order),
    TP53_mutant    = list(plot = p_MT, order = site_order),
    TP53_wild_type = list(plot = p_WT, order = site_order)
  ))
}
## ---- [END NEW | DPS meta-FDR heatmap (ggplot2)] ----


plot_dps_meta_heatmap_for_pathways(
  target_pathways = c(  #"PATH-NP_NOTCH_PATHWAY"
    # "KINASE-PSP_CDK1"
    # "KINASE-PSP_CDK2"
    #"KINASE-PSP_Akt1/AKT1"
    # "KINASE-PSP_PKACA/PRKACA"
    # "PATH-NP_EGFR1_PATHWAY"
     "PATH-BI_ISCHEMIA"
    # "PATH-NP_NOTCH_PATHWAY"
  ),
  top_n_pos = 10, top_n_neg = 10,
  version = "BatchAdj",
  base_root = if (is.null(COMBO_PREFIX_DPS)) "phospho_DPS_meta" else file.path(COMBO_PREFIX_DPS, "phospho_DPS_meta"),
  col_up = "#67001F", col_down = "#053061",
  out_formats = c("tiff"),
  fig_width_in = 12,
  fig_height_in = 9,
  sig_dot_size = 3 ## <- Add this line to adjust dot size
)
