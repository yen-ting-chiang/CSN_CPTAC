## =========================================================
## CSN subunits → proteomic GSEA (CPTAC/TCGA, add TP53 layering)
## - Output to strata = c("ALL","TP53_mutant","TP53_wild_type")
## <dataset>/csn_gsea_results_TP53/<STRATUM>/...
## - Also perform pan-cancer summary
## =========================================================

# --- Working directory (portable ) ---
wd_candidates <- c(
  Sys.getenv("CSN_CPTAC_ROOT", unset = NA_character_),
  "C:/Users/danny/Documents/R_project/CSN_CPTAC",
  "C:/Users/cmuh/Documents/YenTing_document/CSN_CPTAC"
)
wd_candidates <- wd_candidates[!is.na(wd_candidates) & wd_candidates != ""]
hit <- wd_candidates[dir.exists(wd_candidates)][1]
if (!is.na(hit)) setwd(hit)
getwd()


## ===== Package=====
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(janitor)
  library(glue)
  library(limma)
  library(imputeLCMD)
  library(msigdbr)
  library(fgsea)
  library(openxlsx)
  library(cowplot)
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
        if (length(lv) <= 1) {
          keep[cn] <- FALSE
          if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
        } else {
          df[[cn]] <- v
        }
      } else {
        df[[cn]] <- suppressWarnings(as.numeric(v))
      }
    }
    df <- df[, keep, drop = FALSE]
    df
  }
}


## ====Force single-threaded execution, disable all parallel backends ====
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

  options(.fgsea_nproc = 1L)
}
.force_serial_execution()

## Create the run_info directory
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)


## Explanation of Multiple-testing policy
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


## ===== USER CONFIG | parameter =====
set.seed(1234)
csn_subunits <- c("GPS1", "COPS2", "COPS3", "COPS4", "COPS5", "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9")
min_frac_complete <- 0.75 # Each gene must contain at least 75% non-NA.
minSize <- 5
maxSize <- 1000
fgsea_eps <- 1e-10
min_per_group <- 8 #  Minimum number of samples for High/Low (skip limma if insufficient)
MAKE_PLOTS <- FALSE
RUN_PASSES <- "BatchAdj"
options(csn.run_passes = RUN_PASSES)
## ---- [END USER CONFIG | parameter] ----

## ---- [USER CONFIG | heatmap top/bottom limits] ----
## Non-H Collection Heatmap Filtering Thresholds (NES for single datasets, Z for PANs; both based on CSN_SCORE)
## - Single dataset: Only retain the top N and top M smallest NES values with padj < 0.05 (based on CSN_SCORE)
## - PAN (meta-FDR): Only retain the top N and top M smallest Z values with padj_meta < 0.05 (based on CSN_SCORE)
DATASET_HEATMAP_TOP_N <- 25 # Default top 25 (NES)
DATASET_HEATMAP_BOTTOM_N <- 25 # Default bottom 25 (NES)
PAN_HEATMAP_TOP_N <- 25 # Default top 25（Z）
PAN_HEATMAP_BOTTOM_N <- 25 # Default bottom 25 (NES)
## ---- [END USER CONFIG | heatmap top/bottom limits] ----


## ---- [USER CONFIG | heatmap collections] ----
## Which collections' heatmaps to plot (string vectors; e.g., c("H", "C6", "C2:CP:BIOCARTA"))
## - Single dataset heatmap: Leave blank or NULL to plot all available collections.
## - PAN heatmap: Leave blank or NULL to plot all available collections.
PLOT_DATASET_COLLECTIONS <- NULL
PLOT_PAN_COLLECTIONS <- NULL
## ---- [END USER CONFIG | heatmap collections] ----


## ---- Pipeline----
.RUN_LIMMA <- TRUE

log_msg <- function(text, ..., .envir = parent.frame()) {
  ts <- format(Sys.time(), "%H:%M:%S")
  msg <- tryCatch(
    {
      if (grepl("\\{[^}]+\\}", text)) {
        glue::glue(text, ..., .envir = .envir)
      } else if (grepl("%", text)) {
        do.call(sprintf, c(list(fmt = text), list(...)))
      } else {
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

## =====Gene set selection (supports all MSigDB collections/subcollections)=====
GENESET_GROUPS_TO_RUN <- c("H")


## ===== Gene sets=====
log_msg(
  "Prepare MSigDB gene sets; build them according to GENESET_GROUPS_TO_RUN：%s",
  paste(GENESET_GROUPS_TO_RUN, collapse = ", ")
)

df0 <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens"), error = function(e) NULL)

if (!is.null(df0) && NROW(df0) > 0) {
  df0 <- as.data.frame(df0)

  ## 1) Detect the actual field name (gs_collection / gs_subcollection in 2025.1)
  cat_candidates <- c("gs_collection", "gs_cat", "gs_category", "category", "collection")
  sub_candidates <- c("gs_subcollection", "gs_subcat", "gs_subcategory", "subcategory", "sub_category")
  name_candidates <- c("gs_name", "geneset_name", "set_name")

  cat_hits <- intersect(cat_candidates, names(df0))
  sub_hits <- intersect(sub_candidates, names(df0))
  name_hits <- intersect(name_candidates, names(df0))

  if (length(cat_hits) == 0 || length(name_hits) == 0) {
    log_msg("(Omitted) The required field (collection or gs_name) could not be found; no collections list was output.")
  } else {
    catcol <- cat_hits[1]
    subcol <- if (length(sub_hits) >= 1) sub_hits[1] else NULL
    namecol <- name_hits[1]

    ## 2) Group the data frames and count the number of unique genesets (gs_name).
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

    ## 3) Sorting (treat NA subcategories as empty strings to avoid NA in the order)
    ord <- order(counts$gs_cat, ifelse(is.na(counts$gs_subcat), "", counts$gs_subcat))
    counts <- counts[ord, , drop = FALSE]

    ## 4) Output the manifest to run_info
    dir.create("run_info", showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(counts, file.path("run_info", "msigdb_available_collections.csv"))
    log_msg("The available collections for msigdb (including geneset counts) have been output to run_info/msigdb_available_collections.csv")
  }
} else {
  log_msg("(Omitted) msigdbr() returns null or the call failed; no collections list is output.")
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
    log_msg("  Note: msigdbr cannot provide category %s；Skip %s", cat, pp$raw)
    next
  }

  if (!is.null(pp$sub)) {
    subcol <- .get_subcat_col(df)
    subtxt <- if (!is.null(subcol)) as.character(df[[subcol]]) else ""
    # Match only at the "beginning" of the subclass string (e.g., MIR, CP:REACTOME, GO:BP)
    # and allows it to be followed by ":" or the end.
    esc <- function(s) gsub("([\\^\\$\\.|\\(\\)\\[\\]\\{\\}\\*\\+\\?\\\\])", "\\\\\\1", s)
    pattern <- paste0("^", esc(pp$sub), "($|:)")
    keep <- grepl(pattern, subtxt, ignore.case = TRUE)
    df <- df[keep, , drop = FALSE]
    if (!nrow(df)) {
      log_msg("  Note: %s does not find any subset of %s; skip this step.", cat, pp$sub)
      next
    }
    grp <- .make_group_label(cat, df, pp$sub)
  } else {
    grp <- cat
  }

  genesets_by_group[[grp]] <- lapply(split(df$gene_symbol, df$gs_name), unique)
  log_msg("  gey %s：%d  sets", grp, length(genesets_by_group[[grp]]))
}

if (!length(genesets_by_group)) {
  stop("No gene-set is available. Please check the GENESET_GROUPS_TO_RUN setting.")
}
log_msg("final gene-set groups：%s", paste(names(genesets_by_group), collapse = ", "))


## ===== File reading and sharing tools=====
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
  if (!file.exists(fp)) stop(glue::glue("File not found：{fp}"))
  log_msg("Reading the protein matrix：{basename(fp)}")
  dat <- suppressMessages(readr::read_tsv(fp, guess_max = 200000, show_col_types = FALSE))
  gene_cols <- c("Hugo_Symbol", "hugo_symbol", "Gene", "Gene_Symbol", "HugoSymbol", "GENE_SYMBOL", "gene", "gene_symbol")
  gcol <- intersect(gene_cols, names(dat))
  if (!length(gcol)) gcol <- names(dat)[1]
  dat <- dplyr::rename(dat, Gene = !!gcol[1])
  dat$Gene <- sub("\\|.*$", "", dat$Gene)
  not_sample <- c("Gene", "Entrez_Gene_Id", "Entrez_Gene_Id.", "ENTREZ_GENE_ID", "Description", "Gene_Name", "GeneName", "Gene_Symbol")
  sample_cols_all <- setdiff(names(dat), not_sample)

  case_file <- file.path(dir, "case_lists", "cases_protein_quantification.txt")
  keep_ids <- read_case_list(case_file)
  if (length(keep_ids)) {
    inter <- intersect(sample_cols_all, keep_ids)
    sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
    if (length(inter) < 10) log_msg("hint：case_list intersection too small（{length(inter)}）Use all sample fields instead")
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
  if (ncol(m) == 0) stop("0 sample fields read")
  if (anyDuplicated(rownames(m))) {
    log_msg("Duplicate genes are detected, the average of the duplicate rows is taken.")
    m <- rowsum(m, group = rownames(m), reorder = FALSE) / as.vector(table(rownames(m)))
  }
  log_msg("Matrix Dimension：{nrow(m)} genes × {ncol(m)} samples")
  m
}

write_geneset_manifest <- function(genesets_by_group, out_csv = file.path("run_info", "geneset_manifest.csv")) {
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  ver <- tryCatch(as.character(utils::packageVersion("msigdbr")), error = function(e) NA_character_)
  rows <- lapply(names(genesets_by_group), function(g) {
    gs <- genesets_by_group[[g]]
    if (is.null(gs) || !length(gs)) {
      return(NULL)
    }
    data.frame(
      group = g,
      pathway = names(gs),
      genes_n = vapply(gs, function(v) length(unique(v)), integer(1)),
      stringsAsFactors = FALSE
    )
  })
  df <- dplyr::bind_rows(rows)
  if (!is.null(df)) data.table::fwrite(df, out_csv)
  return(ver)
}


# Cleanup before PCA: Remove columns that are Inf→NA, retain at least min_samples finite values; interpolate inline median to NA.
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

# safe version CSN SCORE
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
    if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] Insufficient available subcells or samples → All NA")
    return(out_na)
  }
  pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) {
    if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] prcomp failed → all NA")
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

# safe version audit
audit_csn_score_feasibility_safe <- function(ds_id, stratum, mat0, present_sub,
                                             min_members = 5L, pca_min_samples = 10L,
                                             min_per_group = get0("min_per_group", ifnotfound = 8L),
                                             out_dir = file.path("run_info", "csn_score_audit")) {
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
      log_msg("[CSN-audit-safe] %s | %s：insufficient available sub-units or samples, audit skipped (without termination).", ds_id, stratum)
    }
    return(invisible(FALSE))
  }
  pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) {
    if (exists("log_msg", mode = "function")) {
      log_msg("[CSN-audit-safe] %s | %s：fallback prcomp still failed，skip", ds_id, stratum)
    }
    return(invisible(FALSE))
  }
  varpc1 <- (pc$sdev[1]^2) / sum(pc$sdev^2)
  if (exists("log_msg", mode = "function")) {
    log_msg(
      "[CSN-audit-safe] %s | %s：fallback OK；genes=%d；PC1%%=%.1f",
      ds_id, stratum, nrow(X), 100 * varpc1
    )
  }
  invisible(TRUE)
}

# --- helper: Ensure limma t-stat has correct gene names ---
._ensure_stats_names <- function(stats, gene_names, label = NULL) {
  if (is.null(stats)) {
    return(NULL)
  }
  v <- suppressWarnings(as.numeric(stats))
  if (length(v) != length(gene_names)) {
    stop(sprintf(
      "[ensure-names%s] stats(%d) 與 gene_names(%d) Length mismatch",
      if (!is.null(label)) paste0("-", label) else "",
      length(v), length(gene_names)
    ))
  }
  names(v) <- as.character(gene_names)

  v
}


# --- Tools: Clean and sort GSEA stats (keep limited values, add names, remove duplicates, sort) ---
._finite_rank_stats <- function(stats, gene_names = NULL, label = "", min_n = 20L) {
  nm <- names(stats)
  stats <- suppressWarnings(as.numeric(stats))
  if (is.null(nm) && !is.null(gene_names) && length(stats) == length(gene_names)) {
    nm <- gene_names
  }
  names(stats) <- nm

  ok <- is.finite(stats) & !is.na(stats)
  if (sum(!ok) > 0L) {
    if (!is.null(label) && nzchar(label)) {
      log_msg("[gsea-%s] drop non-finite stats: %d", label, sum(!ok))
    }
  }
  stats <- stats[ok]
  nm <- nm[ok]
  names(stats) <- nm
  # The name must exist
  if (is.null(nm)) {
    log_msg("[gsea-%s] stats has no names after filtering → skip", label)
    return(NULL)
  }
  # Remove duplicate names (keep first name)
  dup <- duplicated(nm)
  if (any(dup)) {
    if (!is.null(label) && nzchar(label)) {
      log_msg("[gsea-%s] drop duplicated gene names in stats: %d", label, sum(dup))
    }
    stats <- stats[!dup]
    nm <- nm[!dup]
  }
  names(stats) <- nm
  if (length(stats) < min_n) {
    log_msg("[gsea-%s] too few finite stats after filtering: %d < %d → skip", label, length(stats), min_n)
    return(NULL)
  }
  stats[order(stats, decreasing = TRUE)]
}

# --- Tool: Map pathways to the universe of stats and filter by size ---
._intersect_and_filter_pathways <- function(pathways, universe, minSize, maxSize, label = "") {
  pw <- lapply(pathways, function(gs) unique(intersect(gs, universe)))
  lens <- vapply(pw, length, integer(1))
  keep <- which(lens >= minSize & lens <= maxSize)
  if (length(keep) == 0L) {
    log_msg("[gsea-%s] no pathways within size bounds after intersect (min=%d, max=%d)", label, minSize, maxSize)
    return(NULL)
  }
  pw[keep]
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
  head_show <- utils::head(df[, c("pos", "M_sample", "interest_rn", "nuisance_rn", "interest_ok", "nuisance_ok")], 6)
  log_msg("[align:%s] head:\n%s", tag, utils::capture.output(print(head_show)) |> paste(collapse = "\n"))
  invisible(df)
}


## ===== impute_and_filter: fix seed(5) before imputation =====
impute_and_filter <- function(mat, min_frac = 0.75) {
  keep <- rowMeans(!is.na(mat)) >= min_frac
  m <- mat[keep, , drop = FALSE]
  if (any(is.na(m))) {
    set.seed(1234)
    m <- imputeLCMD::impute.MinProb(m, q = 0.01)
  }
  m
}


read_gsea_table <- function(out_root, subunit, group_name, stat_tag) {
  grp <- safe_fs_name(group_name)
  su <- subunit
  f <- paste0(stat_tag, ".csv")

  candidates <- c(
    file.path(out_root, su, f), #  collection-first: ./<subunit>/STAT.csv
    file.path(out_root, grp, su, "H", f), # ./<group>/<subunit>/H/STAT.csv
    file.path(out_root, grp, su, f), # no H file
    file.path(out_root, su, grp, "H", f), # ./<subunit>/<group>/H/STAT.csv
    file.path(out_root, su, grp, f) # no H file
  )

  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    if (exists("log_msg", mode = "function")) try(log_msg(" file not exist：%s", paste(candidates, collapse = " | ")), silent = TRUE)
    return(NULL)
  }
  fp <- hit[1L]

  dt <- tryCatch(data.table::fread(fp, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
  if (is.null(dt)) {
    return(NULL)
  }

  need_cols <- c("pathway", "NES", "padj")
  if (!all(need_cols %in% names(dt))) {
    nms <- tolower(names(dt))
    if (!"pathway" %in% names(dt) && "pathway" %in% nms) names(dt)[match("pathway", nms)] <- "pathway"
    if (!"NES" %in% names(dt) && "nes" %in% nms) names(dt)[match("nes", nms)] <- "NES"
    if (!"padj" %in% names(dt)) dt$padj <- NA_real_
  }
  dt <- as.data.frame(dt[, c("pathway", "NES", "padj")])
  names(dt)[names(dt) == "NES"] <- paste0("NES_", subunit)
  names(dt)[names(dt) == "padj"] <- paste0("padj_", subunit)
  dt
}

merge_subunit_tables <- function(tbl_list) {
  keep <- tbl_list[!vapply(tbl_list, is.null, logical(1))]
  if (!length(keep)) {
    return(NULL)
  }
  out <- Reduce(function(x, y) dplyr::full_join(x, y, by = "pathway"), keep)
  num_cols <- setdiff(names(out), "pathway")
  out[num_cols] <- lapply(out[num_cols], function(z) suppressWarnings(as.numeric(z)))
  out
}

add_sig_counts <- function(df, alphas = c(0.05, 0.25)) {
  if (is.null(df) || !nrow(df)) {
    return(df)
  }
  padj_cols <- grep("^padj_", names(df), value = TRUE)
  subunits <- sub("^padj_", "", padj_cols)
  nes_cols <- paste0("NES_", subunits)
  keep_idx <- nes_cols %in% names(df)
  padj_cols <- padj_cols[keep_idx]
  nes_cols <- nes_cols[keep_idx]
  if (!length(padj_cols)) {
    return(df)
  }
  for (a in alphas) {
    a_tag <- gsub("\\.", "_", as.character(a))
    sig_mat <- mapply(function(p, n) {
      pv <- df[[p]]
      as.integer(!is.na(pv) & pv < a)
    }, padj_cols, nes_cols, SIMPLIFY = TRUE)
    if (is.null(dim(sig_mat))) sig_mat <- matrix(sig_mat, ncol = 1)
    df[[sprintf("sig_n_padj_%s", a_tag)]] <- rowSums(sig_mat, na.rm = TRUE)
    pos_mat <- mapply(function(p, n) {
      pv <- df[[p]]
      nv <- df[[n]]
      as.integer(!is.na(pv) & pv < a & !is.na(nv) & nv > 0)
    }, padj_cols, nes_cols, SIMPLIFY = TRUE)
    if (is.null(dim(pos_mat))) pos_mat <- matrix(pos_mat, ncol = 1)
    df[[sprintf("pos_n_padj_%s", a_tag)]] <- rowSums(pos_mat, na.rm = TRUE)
    neg_mat <- mapply(function(p, n) {
      pv <- df[[p]]
      nv <- df[[n]]
      as.integer(!is.na(pv) & pv < a & !is.na(nv) & nv < 0)
    }, padj_cols, nes_cols, SIMPLIFY = TRUE)
    if (is.null(dim(neg_mat))) neg_mat <- matrix(neg_mat, ncol = 1)
    df[[sprintf("neg_n_padj_%s", a_tag)]] <- rowSums(neg_mat, na.rm = TRUE)
  }
  df
}

write_summary_outputs <- function(df, out_dir, group_name, stat_tag) {
  if (is.null(df) || !nrow(df)) {
    log_msg("  [Skip the output] {group_name} | {stat_tag} No results available")
    return(invisible(NULL))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(out_dir, paste0("Summary_", safe_fs_name(group_name), "_", stat_tag))
  ord_keys <- c("sig_n_padj_0_05", "pos_n_padj_0_05", "neg_n_padj_0_05")
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
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "ALL")
  openxlsx::writeData(wb, "ALL", df)
  if ("sig_n_padj_0_05" %in% names(df)) {
    df005 <- df |> dplyr::filter(.data[["sig_n_padj_0_05"]] > 0)
    openxlsx::addWorksheet(wb, "padj_lt_0.05")
    openxlsx::writeData(wb, "padj_lt_0.05", df005)
    data.table::fwrite(df005, paste0(base, "_padjLT0.05.csv"))
  }
  if ("sig_n_padj_0_25" %in% names(df)) {
    df025 <- df |> dplyr::filter(.data[["sig_n_padj_0_25"]] > 0)
    openxlsx::addWorksheet(wb, "padj_lt_0.25")
    openxlsx::writeData(wb, "padj_lt_0.25", df025)
    data.table::fwrite(df025, paste0(base, "_padjLT0.25.csv"))
  }
  openxlsx::saveWorkbook(wb, paste0(base, ".xlsx"), overwrite = TRUE)
  log_msg("  [Complete output] {group_name} | {stat_tag} -> {dirname(base)}")
}

summarize_all_groups <- function(out_root, csn_subunits, genesets_by_group,
                                 stat_tags = c("GSEA_limma_t")) {
  sum_root <- file.path(out_root, "summary")
  dir.create(sum_root, recursive = TRUE, showWarnings = FALSE)
  for (grp_name in names(genesets_by_group)) {
    grp_safe <- safe_fs_name(grp_name)
    for (stat_tag in stat_tags) {
      log_msg("== Summary：group={grp_name} | stat={stat_tag} ==")
      tbl_list <- setNames(vector("list", length(csn_subunits)), csn_subunits)
      for (su in csn_subunits) tbl_list[[su]] <- read_gsea_table(out_root, su, grp_name, stat_tag)
      wide <- merge_subunit_tables(tbl_list)
      wide <- add_sig_counts(wide, alphas = c(0.05, 0.25))
      out_dir <- file.path(sum_root, stat_tag)
      write_summary_outputs(wide, out_dir, grp_name, stat_tag)
    }
  }
}

## ===== Meta-FDR via Stouffer (per stratum × predictor × group × stat) =====

meta_fdr_stouffer <- function(dataset_dirs,
                              strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
                              stat_tags = c("GSEA_limma_t_cont"),
                              groups = names(genesets_by_group),
                              out_root = if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr")) {
  stopifnot(length(dataset_dirs) > 0)
  ## [EXCLUDE-FROM-META] exclude lusc_cptac_2021
  if (!is.null(names(dataset_dirs))) {
    keep_idx <- names(dataset_dirs) != "lusc_cptac_2021"
  } else {
    keep_idx <- basename(dataset_dirs) != "lusc_cptac_2021"
  }
  dataset_dirs <- dataset_dirs[keep_idx]
  if (!length(dataset_dirs)) {
    message("[meta] All datasets were excluded or no datasets were available (lusc_cptac_2021 was excluded); skipped.")
    return(invisible(TRUE))
  }
  if (!requireNamespace("data.table", quietly = TRUE)) stop("install data.table first.")
  if (!requireNamespace("stats", quietly = TRUE)) stop("missing stats package")
  data.table::setDTthreads(1L)

  # ---- helpers --------------------------------------------------------------
  normalize_cols <- function(DT) {
    # Map common column names to generic names: pathway, NES, pval, padj, size
    nm <- names(DT)

    low <- tolower(nm)
    map <- c(
      "pathway" = "pathway",
      "term" = "pathway",
      "pathway_name" = "pathway",
      "nes" = "NES",
      "enrichment_score" = "NES",
      "pval" = "pval",
      "p.value" = "pval",
      "pvalue" = "pval",
      "p" = "pval",
      "padj" = "padj",
      "fdr" = "padj",
      "qval" = "padj",
      "size" = "size",
      "setsize" = "size",
      "n" = "size"
    )

    for (i in seq_along(nm)) {
      key <- low[i]
      if (key %in% names(map)) data.table::setnames(DT, nm[i], map[[key]])
    }
    DT
  }

  pick_cols <- function(DT, need) {
    need <- unique(need)
    if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)

    miss <- setdiff(need, names(DT))
    for (m in miss) DT[, (m) := NA]

    DT[, ..need]
  }

  stouffer_df <- function(D) {
    D <- D[is.finite(pval) & is.finite(NES)]
    if (nrow(D) == 0) {
      return(NULL)
    }
    # Construct a signed Z using "one-tailed p-values + NES direction"
    # Clamp pval to avoid Inf (when p=0) or -Inf (when p=1)
    D[, pval := pmin(pmax(pval, 1e-10), 1 - 1e-10)]
    D[, z := sign(NES) * stats::qnorm(p = pval, lower.tail = FALSE)]
    # Unweighted Stouffer
    Out <- D[, .(k = .N, Z = sum(z) / sqrt(.N)), by = .(pathway)]
    Out[, p_meta := 2 * stats::pnorm(abs(Z), lower.tail = FALSE)]
    Out[, padj_meta := p.adjust(p_meta, method = "BH")]
    data.table::setorder(Out, p_meta)
    Out[]
  }

  # ---- main loop ------------------------------------------------------------
  passes <- c("BatchAdj")
  base_out <- out_root
  dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

  for (stratum in strata) {
    for (pass_label in passes) {
      for (stat_tag in stat_tags) {
        for (grp in groups) {
          # 1) Find all subunits (by unioning the folders of each dataset)
          subunits <- unique(unlist(lapply(dataset_dirs, function(dsdir) {
            base <- file.path(
              GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection",
              safe_fs_name(grp), pass_label, basename(dsdir), stratum
            )
            if (!dir.exists(base)) {
              return(character(0))
            }
            # List only the top-level subdirectory names (each subunit)
            subs <- list.dirs(base, full.names = FALSE, recursive = FALSE)
            subs[subs != ""]
          })))

          if (length(subunits) == 0) {
            message(sprintf(
              "[meta] %s/%s/%s/%s no available subunit，skip",
              stratum, pass_label, stat_tag, grp
            ))
            next
          }

          # 2) For each subunit, collect the CSV files from each dataset, perform Stouffer merging, and archive them.
          for (su in subunits) {
            parts <- list()
            for (dsdir in dataset_dirs) {
              f <- file.path(
                GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection",
                safe_fs_name(grp), pass_label, basename(dsdir), stratum, su,
                sprintf("%s.csv", stat_tag)
              )
              if (!file.exists(f)) next
              dt <- tryCatch(data.table::fread(f), error = function(e) NULL)
              if (is.null(dt) || nrow(dt) == 0) next
              dt <- normalize_cols(dt)
              dt <- pick_cols(dt, need = c("pathway", "NES", "pval", "padj", "size"))
              dt[, dataset := basename(dsdir)]
              parts[[length(parts) + 1L]] <- dt
            }

            if (length(parts) == 0) {
              next
            }

            D <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)

            D <- D[!is.na(pathway)]
            res <- stouffer_df(D)
            if (is.null(res)) next

            # 3) write file
            grp_safe <- safe_fs_name(grp)
            out_dir <- file.path(base_out, stratum, pass_label, su, grp_safe)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            out_file <- file.path(out_dir, paste0(stat_tag, "_meta_fdr.csv"))
            data.table::fwrite(res, out_file)
            message(sprintf(
              "[meta] %s/%s/%s/%s -> %s (n=%d pathways, k>=1 datasets)",
              stratum, pass_label, su, stat_tag, out_file, nrow(res)
            ))
          } # su
        } # grp
      } # stat_tag
    } # pass
  } # stratum

  invisible(TRUE)
}


# Perform PCA using the "cross-sample z-values" of subunits, and take PC1 as the CSN score; adjust the direction to have the same sign as the mean z.
build_csn_score <- function(mat0,
                            subunits = csn_subunits,
                            combine_7AB = TRUE,
                            min_members = 5L) {
  present <- intersect(subunits, rownames(mat0))
  # Pre-create the return skeleton
  s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
  if (!length(present)) {
    return(s)
  }

  # z-score (preserve sample names)
  get_z <- function(v) {
    nm <- names(v) # Save the sample name first
    v <- as.numeric(v)
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    out <- (v - mu) / sdv
    names(out) <- nm # Replace sample name
    out
  }

  X <- do.call(rbind, lapply(present, function(g) get_z(mat0[g, ])))
  rownames(X) <- present
  colnames(X) <- colnames(mat0) # Ensure sample names are available

  # combine COPS7A/7B
  if (combine_7AB && all(c("COPS7A", "COPS7B") %in% rownames(X))) {
    Z7 <- colMeans(X[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
    X <- rbind(X[setdiff(rownames(X), c("COPS7A", "COPS7B")), , drop = FALSE],
      "COPS7*" = Z7
    )
  }


  enough <- colSums(is.finite(mat0[present, , drop = FALSE])) >= min_members
  keep_sam <- names(s)[enough]

  if (length(keep_sam) >= 10) {
    pc <- stats::prcomp(t(X[, keep_sam, drop = FALSE]), center = TRUE, scale. = FALSE)
    sc <- pc$x[, 1]
    # Direction correction: Same sign as subunit average z
    mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)
    if (suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs")) < 0) sc <- -sc
    s[keep_sam] <- sc
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
      mean(!is.na(v)) * 100
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


## ===== Batch Cleanup Settings  =====
BATCH_PIPE_POLICY <- "NA" # "NA" or "b_small": Should values containing '|' be set to NA or merged into b_small?
BATCH_MIN_PER_LEVEL <- 2 # Levels below this threshold will be merged into b_small.

## ---- Batch value cleanup: Handling values containing '|', validating names, merging sparse levels ----
sanitize_batch_levels <- function(x,
                                  pipe_policy = BATCH_PIPE_POLICY,
                                  min_per_level = BATCH_MIN_PER_LEVEL) {
  x0 <- as.character(x)
  has_pipe <- grepl("\\|", x0 %||% "")
  if (any(has_pipe)) {
    log_msg("  [batch] Detects {sum(has_pipe)} values containing '|' → Processes according to policy {pipe_policy}")
    x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
  }

  fac <- factor(make.names(x0))
  fac <- droplevels(fac)

  # Merge sparse levels
  if (!is.null(min_per_level) && min_per_level > 1) {
    tab <- table(fac, useNA = "no")
    small <- names(tab)[tab < min_per_level]
    if (length(small)) {
      log_msg(
        "  [batch] Merge sparse level to 'b_small'：%s",
        paste(sprintf("%s(n=%d)", small, as.integer(tab[small])), collapse = ", ")
      )
      fac_chr <- as.character(fac)
      fac_chr[fac_chr %in% small] <- "b_small"
      fac <- droplevels(factor(fac_chr))
    }
  }
  fac
}

## ---- Automatically detect batch field  ----
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
    # There must be at least 2 valid levels and the number of valid samples must be >= 3
    if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
      return(list(name = cn, fac = fac))
    }
  }
  NULL
}

## ---- Retrieve batch factor in order of sample ----
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

  # align according to sample_ids to ensure that the detected and returned vector lengths are consistent.
  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
  rownames(meta) <- sample_ids

  det <- detect_batch_column(meta,
    pipe_policy   = pipe_policy,
    min_per_level = min_per_level
  )
  if (is.null(det)) {
    ## === Fallback: Pushing back TMT-plex from TMT_protein.csv when batch ===
    tmt_fp <- file.path(ds_dir, "TMT_protein.csv")
    if (file.exists(tmt_fp)) {
      tmt <- suppressMessages(readr::read_csv(tmt_fp, show_col_types = FALSE)) |> as.data.frame()

      ## Normalize field names (remove leading and trailing whitespace; handle "Run Metadata ID")
      cn <- names(tmt)
      cn_trim <- trimws(cn)
      names(tmt) <- cn_trim
      run_hits <- grep("^Run\\s*Metadata\\s*ID$", cn_trim, ignore.case = TRUE, value = TRUE)
      tmt_cols <- grep("^tmt_", cn_trim, ignore.case = TRUE, value = TRUE)

      if (length(run_hits) >= 1 && length(tmt_cols) >= 1) {
        run_col <- run_hits[1]

        ## Create a sample_id → plex mapping
        plex_by_sample <- list()
        nR <- nrow(tmt)
        for (i in seq_len(nR)) {
          run_id <- as.character(tmt[[run_col]][i])
          if (!nzchar(run_id) || is.na(run_id)) next
          plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id) # Take the first two codes
          if (!nzchar(plex2) || is.na(plex2)) next

          for (tc in tmt_cols) {
            cell <- tmt[[tc]][i]
            if (is.na(cell)) next
            cell <- as.character(cell)
            if (!nzchar(cell)) next
            sid <- sub("\\r?\\n.*$", "", cell)
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
          ## Success Conditions: At least two valid levels and ≥ 3 valid samples
          if (nlevels(fac2) >= 2 && sum(!is.na(fac2)) >= 3) {
            names(fac2) <- sample_ids
            return(list(name = "TMT_protein.csv:RunMetadataID", fac = fac2))
          }
        }
      }
    }
    ## If fallback also fails to retrieve the value, return NULL.
    return(NULL)
  }

  fac <- det$fac
  names(fac) <- sample_ids
  list(name = det$name, fac = fac)
}

## ---- Batch Demand Quick Screening ----
screen_batch_need <- function(ds_dir, min_frac_complete = 0.75) {
  log_msg("== Batch check：%s ==", basename(ds_dir))
  mat0 <- load_matrix_from_dataset_dir(ds_dir)
  mx <- suppressWarnings(max(mat0, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)
  m <- impute_and_filter(mat0, min_frac = min_frac_complete)

  bi <- get_batch_factor(ds_dir, colnames(m))
  if (is.null(bi)) {
    log_msg("  [batch] No specific batch field found → Do not correct")
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
    "  [batch] column：%s；levels=%d；each level n：%s",
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
    "  PCA (by batch) R²：%s ; p：%s", paste(round(r2, 3), collapse = ", "),
    paste(signif(pv, 3), collapse = ", ")
  )

  design <- model.matrix(~batch)
  fit <- limma::lmFit(m, design)
  fit <- limma::eBayes(fit)
  padj <- p.adjust(fit$F.p.value, "BH")
  prop05 <- mean(padj < 0.05, na.rm = TRUE)
  prop25 <- mean(padj < 0.25, na.rm = TRUE)
  log_msg(
    "  gene level F test：FDR<0.05 ratio = %.1f%%；FDR<0.25 = %.1f%%",
    100 * prop05, 100 * prop25
  )

  recommend <- (any(r2 >= 0.10 & pv[seq_along(r2)] < 0.01) || prop05 >= 0.05)
  if (recommend) {
    log_msg("  **Recommended Correction**: R² or gene ratio reaches the threshold.（≥10%% R² and p<0.01，or FDR<0.05 gene ≥5%%）")
  } else {
    log_msg("  **No correction needed for now:** No significant batch impact observed (record retained)")
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


## ---- [CLEAN | FIXED DATASETS & PREFIX] ----
## 7 CPTAC datasets actually used
dataset_ids <- c(
  "brca_cptac_2020",
  "luad_cptac_2020",
  "lusc_cptac_2021",
  "ucec_cptac_2020",
  "coad_cptac_2019",
  "gbm_cptac_2021",
  "paad_cptac_2021"
)

## Unified output root directory prefix:
GSEA_PREFIX <- "proteomic_GSEA"
## ---- [END CLEAN | FIXED DATASETS & PREFIX] ----


dataset_dirs <- setNames(file.path(datasets_root, dataset_ids), dataset_ids)
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")
message("datasets_root = ", datasets_root)


missing_dirs <- names(dataset_dirs)[!dir.exists(dataset_dirs)]
if (length(missing_dirs)) {
  log_msg("detected %d files，will skip：%s", length(missing_dirs), paste(missing_dirs, collapse = ", "))
}


## List of datasets to be run and available in this round (data folders containing data_protein_quantification.txt)
dataset_dirs_run <- dataset_dirs[
  dir.exists(dataset_dirs) &
    file.exists(file.path(dataset_dirs, "data_protein_quantification.txt"))
]

if (!length(dataset_dirs_run)) stop("`dataset_dirs_run` is empty,  check if the folder and `data_protein_quantification.txt` exist.")

log_msg("dataset available this round（%d）：%s", length(dataset_dirs_run), paste(names(dataset_dirs_run), collapse = ", "))


## ===== Obtain TP53 mutation status =====
## ===== TP53 status (protein-altering baseline) =====
## Only mutations that alter the protein are considered TP53-mutants; all others (Silent, UTR, Intron, IGR, RNA, lincRNA, Flank, etc.) are considered wild-type.

TP53_KEEP_CLASSES <- c(
  "MISSENSE_MUTATION", "NONSENSE_MUTATION",
  "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
  "IN_FRAME_DEL", "IN_FRAME_INS",
  "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
)

normalize_vc <- function(x) {
  x <- toupper(trimws(as.character(x)))
  gsub("[^A-Z0-9]+", "_", x)
}

get_tp53_status <- function(ds_dir, sample_ids) {
  status <- setNames(rep("TP53_wild_type", length(sample_ids)), sample_ids)

  mut_fp <- file.path(ds_dir, "data_mutations.txt")
  if (!file.exists(mut_fp)) {
    log_msg("  [TP53] data_mutations.txt not found; the entire batch is considered wild-type/ALL and is available.")
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
    log_msg("  [TP53] data_mutations.txt missing column：%s → treated as wild-type", paste(miss, collapse = ", "))
    return(status)
  }

  tp53_df <- subset(mutation_df, Hugo_Symbol == "TP53",
    select = c("Variant_Classification", "Tumor_Sample_Barcode")
  )
  if (!nrow(tp53_df)) {
    return(status)
  }

  vc_norm <- normalize_vc(tp53_df$Variant_Classification)

  # Strictly adheres to the protein-altering category; while tolerating prefixes such as FRAME_SHIFT_* / IN_FRAME_*.
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
## output：
##   - run_info/tp53_status/<dataset>_tp53_class_sample_counts.csv
##   - run_info/tp53_status/tp53_class_sample_counts_long.csv
##   - run_info/tp53_status/tp53_binary_counts_by_dataset.csv

dir.create(file.path("run_info", "tp53_status"), recursive = TRUE, showWarnings = FALSE)


# Summary of TP53 sample counts for a single dataset
summarize_tp53_counts_for_dataset <- function(ds_dir) {
  ds_id <- basename(ds_dir)


  M <- try(load_matrix_from_dataset_dir(ds_dir), silent = TRUE)
  if (inherits(M, "try-error")) {
    log_msg("[TP53-audit] %s: cannot read matrix, skip", ds_id)
    return(NULL)
  }
  sample_ids <- colnames(M)
  sid_up <- toupper(sample_ids)
  n_all <- length(sample_ids)

  # obtain mutant/wild type
  status <- get_tp53_status(ds_dir, sample_ids)
  tb_bin <- table(factor(status, levels = c("TP53_wild_type", "TP53_mutant")))
  bin_row <- data.frame(
    dataset = ds_id,
    in_matrix_n = n_all,
    WT_n = as.integer(tb_bin["TP53_wild_type"]),
    MUT_n = as.integer(tb_bin["TP53_mutant"]),
    stringsAsFactors = FALSE
  )

  # read MAF
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

        # Only samples "in the protein matrix" are counted.
        tp53 <- tp53[tp53$Tumor_Sample_Barcode %in% sid_up, , drop = FALSE]

        # Each classification → How many "independent samples" are matched?
        class_counts <- tp53 |>
          dplyr::group_by(Variant_Classification) |>
          dplyr::summarise(sample_n = dplyr::n_distinct(Tumor_Sample_Barcode), .groups = "drop") |>
          dplyr::arrange(dplyr::desc(sample_n), Variant_Classification)

        # Add summary columns for "Any_TP53_mutation" and "protein_altering" (sample level)
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

  # write file
  if (!is.null(class_df)) {
    data.table::fwrite(
      class_df,
      file.path("run_info", "tp53_status", paste0(ds_id, "_tp53_class_sample_counts.csv"))
    )
  } else {
    data.table::fwrite(
      data.frame(dataset = ds_id, Variant_Classification = NA_character_, sample_n = NA_integer_),
      file.path("run_info", "tp53_status", paste0(ds_id, "_tp53_class_sample_counts.csv"))
    )
  }

  list(binary = bin_row, class_long = class_df)
}

# Run all datasets and compile a summary table
summarize_tp53_counts_all_datasets <- function(dataset_dirs) {
  all_bin <- list()
  all_class <- list()
  k <- 1L
  j <- 1L
  for (ds in names(dataset_dirs)) {
    ds_dir <- dataset_dirs[[ds]]
    if (!dir.exists(ds_dir)) next
    log_msg("[TP53-audit] starts：%s", ds)
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
    log_msg("[TP53-audit] write：tp53_binary_counts_by_dataset.csv")
  }
  if (length(all_class)) {
    class_df <- dplyr::bind_rows(all_class)
    data.table::fwrite(class_df, file.path("run_info", "tp53_status", "tp53_class_sample_counts_long.csv"))
    log_msg("[TP53-audit] write：tp53_class_sample_counts_long.csv")
  }
}

## === call ===
## summarize_tp53_counts_all_datasets(dataset_dirs)

## === robust run_predictor_analyses (center predictors for limma models) ===
run_predictor_analyses <- function(
  predictor_name,
  predictor_vec, # named numeric by sample
  exclude_genes = NULL, # Genes to be excluded from the ranking (e.g., autologous/entire complex)
  ds_id, ds_dir,
  mat0, # limma continuous method
  mat, # after impute+filter
  out_root,
  genesets_by_group, # list(group -> pathways)
  batch_all = NULL, # factor by sample or NULL
  purity_all = NULL, # named numeric by sample or NULL
  sa_all = NULL, # data.frame(sex, age[, age_missing, age_z_imputed]); rownames=sample
  tp53_num_all = NULL, # named numeric (0/1) or NULL
  is_ALL = FALSE, # TP53_mutant only be added when TRUE is active.
  extra_covars_df = NULL # Optional: Additional covariates (e.g. CSN_SCORE)
) {
  ## -------- parameter --------
  opt <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
  min_per_group <- opt("min_per_group", 8L)
  minSize <- opt("minSize", 15L)
  maxSize <- opt("maxSize", 500L)
  fgsea_eps <- opt("fgsea_eps", 0)
  USE_AGE_MISSING_INDICATOR <- isTRUE(opt("USE_AGE_MISSING_INDICATOR", FALSE))

  logf <- function(...) if (exists("log_msg", mode = "function")) try(log_msg(...), silent = TRUE)


  stopifnot(!is.null(genesets_by_group), length(genesets_by_group) > 0)

  ## -------- Aligning Samples ( based on the predictor) --------
  if (is.null(names(predictor_vec))) names(predictor_vec) <- colnames(mat)
  keep <- intersect(colnames(mat), names(predictor_vec))
  pred_all <- suppressWarnings(as.numeric(predictor_vec[keep]))
  names(pred_all) <- keep
  fin <- is.finite(pred_all)
  if (sum(fin) < (2L * min_per_group)) {
    logf("  [%s] predictor Insufficient non-NA samples（%d < %d）→ skip", predictor_name, sum(fin), 2L * min_per_group)
    return(invisible(NULL))
  }
  sample_order <- keep[fin]
  pred <- pred_all[fin]
  stopifnot(identical(names(pred), sample_order))

  ## -------- Construct common variable data frames --------
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
    # Merge extra covariates if provided
    if (!is.null(extra_covars_df)) {
      # Ensure extra_covars_df has rownames
      if (!is.null(rownames(extra_covars_df))) {
        common_ex <- intersect(so, rownames(extra_covars_df))
        if (length(common_ex) > 0) {
          ex_sub <- extra_covars_df[so, , drop = FALSE]
          df <- cbind(df, ex_sub)
        }
      }
    }

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
        v <- factor(v)

        lv <- levels(droplevels(v[!is.na(v)]))
        if (length(lv) <= 1) {
          keep[cn] <- FALSE
          if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
        } else {
          df[[cn]] <- v
        }
      } else {
        df[[cn]] <- suppressWarnings(as.numeric(v))
      }
    }

    df <- df[, keep, drop = FALSE]
    df
  }


  ## -------- coverage gate --------
  base_thr <- c(purity = 0.60, sex = 0.80, age = 0.80)
  present <- colnames(df_covars0)
  extra <- setdiff(present, names(base_thr))
  if (length(extra)) base_thr <- c(base_thr, stats::setNames(rep(min(base_thr), length(extra)), extra))

  ## -------- Select common variables-------
  pick_covars_df <- function(label) {
    ## Determine the source of common variables (if covars_all is not in this scope, fall back to df_covars0)
    cov_src <- tryCatch(get("covars_all", inherits = TRUE), error = function(e) NULL)
    if (is.null(cov_src)) {
      cov_src <- tryCatch(get("df_covars0", inherits = TRUE), error = function(e) NULL)
    }
    if (is.null(cov_src)) {
      logf("  [covars-%s] NO covariate table in scope → skip", label)
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

    dl_sel <- tryCatch(attr(sel, "drop_log"), error = function(e) NULL)

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
          logf("  [covars-%s] Return shape/alignment unclear, discard.", label)
          return(NULL)
        }
      }
      X <- X[sample_order, , drop = FALSE]
    }
    if (is.null(X) || !ncol(X)) {
      return(NULL)
    }
    X <- coerce_covariates_safely(X)

    good <- vapply(X, function(v) sum(!is.na(v)) >= 3, logical(1))
    if (any(!good)) {
      logf("  [covars-%s] drop low-coverage columns: %s", label, paste(colnames(X)[!good], collapse = ","))
      X <- X[, good, drop = FALSE]
    }
    if (!ncol(X)) {
      return(NULL)
    }
    logf("  [covars-picked:%s] kept = {%s}", label, if (!is.null(X)) paste(colnames(X), collapse = ",") else "NULL")
    ## === Mapping the discarded covariates and their causes to a CSV file ===
    if (!is.null(dl_sel) && nrow(dl_sel)) {
      dl_sel$dataset <- ds_id
      dl_sel$stratum <- basename(out_root)
      dl_sel$subunit <- predictor_name
      dl_sel$pass <- label # dl_sel$pass <- label # "limma-cont:base"
      dl_sel$samples <- length(sample_order)

      out_dir <- file.path("run_info", "covars_audit")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      fp <- file.path(out_dir, "covariate_drop_reasons.csv")

      data.table::fwrite(dl_sel, fp, append = file.exists(fp))
    }

    X
  }

  X_ba_cov <- pick_covars_df("limma-cont:base")
  if (!is.null(X_ba_cov)) X_ba_cov <- coerce_covariates_safely(X_ba_cov)

  ## -------- BatchAdj: Same as above; batch is included. --------
  DF_ba <- data.frame(predictor = as.numeric(pred), row.names = sample_order, check.names = FALSE)
  if (!is.null(X_ba_cov)) for (nm in colnames(X_ba_cov)) DF_ba[[nm]] <- X_ba_cov[[nm]]
  if (!is.null(batch_all)) DF_ba$batch <- droplevels(as.factor(batch_all[sample_order]))

  ## delete all NA column
  all_na_col <- vapply(DF_ba, function(v) all(is.na(v)), logical(1))
  if (any(all_na_col)) {
    logf(
      "[BatchAdj] drop all-NA columns: %s",
      paste(names(all_na_col)[all_na_col], collapse = ",")
    )
    DF_ba <- DF_ba[, !all_na_col, drop = FALSE]
  }

  ## Retrieve the columns that will actually be used in this round (predictor + selected covars)
  use_cols <- unique(c(
    "predictor",
    intersect(c("purity", "sex", "age", "batch", "TP53_mutant"), colnames(DF_ba))
  ))
  use_cols <- use_cols[use_cols %in% colnames(DF_ba)] # 防守

  ok_rows <- complete.cases(DF_ba[, use_cols, drop = FALSE])
  n_ok <- sum(ok_rows)

  if (n_ok < 16L) {
    logf("  [%s|BatchAdj] available samples not enough（%d < 16）→ skip", predictor_name, n_ok)
    return(invisible(NULL))
  }

  DF_ba <- DF_ba[ok_rows, , drop = FALSE]

  ## Centralization
  DF_ba$predictor <- DF_ba$predictor - mean(DF_ba$predictor, na.rm = TRUE)

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
    covars     = DF_ba[, intersect(c("purity", "sex", "age", "batch"), colnames(DF_ba)), drop = FALSE]
  )


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


  ## -------- tools --------
  coef_name <- function(des) if ("predictor" %in% colnames(des)) "predictor" else colnames(des)[ncol(des)]

  gsea_from <- function(pw, stats) {
    if (is.null(pw) || !length(pw) || is.null(stats) || !length(stats)) {
      return(NULL)
    }

    minSize <- opt("minSize", 15L)
    maxSize <- opt("maxSize", 500L)

    orig_n <- length(pw)
    pw_use <- ._intersect_and_filter_pathways(pw, names(stats), minSize = minSize, maxSize = maxSize)
    message(sprintf("[gsea] (gsea_from) |H| total=%d → use=%d", orig_n, length(pw_use)))
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
        message(sprintf("[gsea] Multilevel failed：%s → change to fgseaSimple", conditionMessage(e)))
        suppressWarnings(fgsea::fgseaSimple(
          pathways = pw_use, stats = stats,
          nperm = 10000, minSize = minSize, maxSize = maxSize
        ))
      }
    )
    res
  }

  ## -------- run each group：BatchAdj ----------
  if (.RUN_LIMMA) {
    for (grp_name in names(genesets_by_group)) {
      pw <- genesets_by_group[[grp_name]]


      out_root_coll <- file.path(GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), ds_id, basename(out_root))

      ## A2
      if ("BatchAdj" %in% getOption("csn.run_passes", c("BatchAdj"))) {
        out_dir_A2 <- file.path(GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root), predictor_name)
        dir.create(out_dir_A2, recursive = TRUE, showWarnings = FALSE)
        fit2 <- limma::eBayes(limma::lmFit(M_ba, des_ba))
        t2 <- fit2$t[, coef_name(des_ba)]
        t2 <- ._ensure_stats_names(t2, rownames(M_ba))
        if (exists("exclude_genes") && length(exclude_genes)) {
          t2 <- t2[setdiff(names(t2), exclude_genes)]
        }
        t2 <- ._finite_rank_stats(t2, label = paste0("A2-", predictor_name))
        if (!is.null(t2)) {
          res2 <- gsea_from(pw, t2)
          if (!is.null(res2) && nrow(res2)) data.table::fwrite(res2, file.path(out_dir_A2, "GSEA_limma_t_cont.csv"))
          ## --- interaction (BatchAdj; A2) ---
          if (isTRUE(is_ALL) && !is.null(tp53_num_all)) {
            tp53_ba_num <- suppressWarnings(as.numeric(tp53_num_all[rownames(DF_ba)]))
            tp53_ba_fac <- factor(ifelse(tp53_ba_num == 1, "MT", "WT"), levels = c("WT", "MT"))

            ## If TP53 in the ALL layer has only a single level, skip it directly.
            if (nlevels(droplevels(tp53_ba_fac)) >= 2) {
              covars_only_ba <- DF_ba[, setdiff(
                colnames(DF_ba),
                c("predictor", "TP53_mutant", "TP53_status", "TP53")
              ),
              drop = FALSE
              ]
              df_int_ba <- data.frame(
                pred = DF_ba$predictor,
                tp53 = tp53_ba_fac,
                covars_only_ba,
                row.names = rownames(DF_ba),
                check.names = FALSE
              )
              des_int_ba <- stats::model.matrix(~ pred * tp53 + ., data = df_int_ba)

              if (.RUN_LIMMA) {
                fit_int_ba <- limma::eBayes(limma::lmFit(M_ba, des_int_ba))
                coef_int <- "pred:tp53MT"
                if (coef_int %in% colnames(coef(fit_int_ba))) {
                  tt <- limma::topTable(fit_int_ba, coef = coef_int, number = nrow(M_ba), sort.by = "none")
                  tvec <- if ("t" %in% names(tt)) tt$t else tt$logFC
                  names(tvec) <- rownames(tt)

                  ## Ranking Processing Consistent with the Main pipeline
                  ranks <- ._ensure_stats_names(tvec, rownames(M_ba))
                  ranks <- ._finite_rank_stats(ranks, label = paste0("A2-", predictor_name, "-TP53_interaction"))

                  ## Perform GSEA on each collection and write to GSEA_limma_interaction.csv
                  res_fg <- gsea_from(pw, ranks)
                  if (!is.null(res_fg) && nrow(res_fg)) {
                    data.table::fwrite(
                      as.data.frame(res_fg),
                      file.path(out_dir_A2, "GSEA_limma_interaction.csv")
                    )
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


## ===== CSN subunits coverage & CSN_SCORE feasibility  =====
audit_csn_score_feasibility <- function(ds_id, stratum, mat0, present_sub,
                                        min_members = 5L,
                                        pca_min_samples = 10L,
                                        min_per_group = 8L,
                                        out_dir = file.path("run_info", "csn_score_audit")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  n_samples <- ncol(mat0)
  if (!is.matrix(mat0) || is.null(rownames(mat0)) || is.null(colnames(mat0))) {
    log_msg("[CSN-audit] %s | %s：mat0 imcomplete，skip", ds_id, stratum)
    return(invisible(NULL))
  }
  present_sub <- intersect(present_sub, rownames(mat0))
  if (!length(present_sub) || n_samples == 0) {
    log_msg("[CSN-audit] %s | %s：no available CSN subunit or samples= 0，skip", ds_id, stratum)
    return(invisible(NULL))
  }

  ## 1) Coverage per subunit
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

  ## 2) How many non-NA subunits are there for each sample; the number of samples that satisfy >= min_members.
  sample_counts <- colSums(is.finite(mat0[present_sub, , drop = FALSE]))
  enough <- sample_counts >= min_members
  n_enough <- sum(enough)

  ## 3) Attempt to calculate CSN_SCORE(PC1) and record the feasibility.
  csn_score <- build_csn_score(mat0,
    subunits = present_sub,
    combine_7AB = TRUE, min_members = min_members
  )
  csn_nonNA <- sum(is.finite(csn_score))
  csn_can_pca <- (n_enough >= pca_min_samples) && (csn_nonNA >= pca_min_samples)

  ## Additional: PC1 explained variance
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

  ## 4) write files
  # 4a) summary
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

  # 4b) subunit coverage
  tag <- paste(ds_id, stratum, sep = "_")
  fp_sub <- file.path(out_dir, sprintf("%s_subunit_coverage.csv", tag))
  data.table::fwrite(sub_tbl, fp_sub)

  # 4c) Number of non-NA subunits per sample（每分層一個檔）
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
    "[CSN-audit] %s | %s：subunits=%d；min/median/max coverage=%.1f/%.1f/%.1f%%；eligible_samples=%d；CSN_SCORE nonNA=%d；PC1 available=%s（PC1%%=%.1f）",
    ds_id, stratum, length(present_sub), cov_min %||% NaN, cov_med %||% NaN, cov_max %||% NaN,
    n_enough, csn_nonNA, as.character(csn_can_pca), pc1_var_pct %||% NaN
  )

  invisible(list(summary = sum_row, per_subunit = sub_tbl, per_sample = sample_tbl))
}


## ===== Execute a complete GSEA function once for the "specified sample set" =====
## ===== run_one_stratum: Standardize filenames, write sample_counts, and  batch comments =====

run_one_stratum <- function(ds_id, ds_dir, mat0_full, sample_keep, out_root, genesets_by_group) {
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  log_msg("  -- stratum：%s | N(sample_keep)=%d", basename(out_root), length(sample_keep))

  ## 1) Subset samples + size check + log scaling
  keep <- intersect(colnames(mat0_full), sample_keep)
  if (length(keep) < 4) {
    log_msg("  [skip] Too few samples：%d", length(keep))
    return(invisible(NULL))
  }
  mat0 <- mat0_full[, keep, drop = FALSE]
  mx <- suppressWarnings(max(mat0, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)

  ## 2) Feasibility audit of limma using imputation + filtering matrix + CSN_SCORE
  mat <- impute_and_filter(mat0, min_frac = min_frac_complete)
  present_sub <- intersect(csn_subunits, rownames(mat0))
  if (!length(present_sub)) {
    log_msg("  [skip] no CSN subunit in this stratum ")
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

  ## === In stratum ALL Analysis: Preparing TP53 State  ===
  is_ALL <- identical(basename(out_root), "ALL")
  tp53_num_all <- NULL
  if (is_ALL) {
    tp53_status_all <- get_tp53_status(ds_dir, colnames(mat0))
    tp53_num_all <- as.numeric(tp53_status_all == "TP53_mutant")
    names(tp53_num_all) <- colnames(mat0)
  }

  ## 3) Covariates / Batch (Aligned Samples) — Establish limma using covariates
  bi_all <- get_batch_factor(ds_dir, colnames(mat0))
  batch_all <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL
  purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
  sa_all <- get_sex_age_covariates(ds_dir, colnames(mat0)) # data.frame(sex, age [, age_missing, age_z_imputed])
  sa_all_limma <- sa_all
  if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing", "age_z_imputed") %in% colnames(sa_all))) {
    keep_cols <- intersect(c("sex", "age", "age_missing", "age_z_imputed"), colnames(sa_all))
    sa_all_limma <- sa_all[, keep_cols, drop = FALSE]
  }

  ## 4) Perform analysis on a per-prediction vector basis: original value of subunit, CSN_SCORE, RESIDUAL_<SU>
  {
    present_sub <- intersect(csn_subunits, rownames(mat0))
    if (!length(present_sub)) {
      log_msg("  [skip] no CSN subunit in this stratum")
      return(invisible(NULL))
    }


    is_ALL <- identical(basename(out_root), "ALL")
    tp53_num_all <- if (is_ALL) {
      status_all <- get_tp53_status(ds_dir, colnames(mat0))
      setNames(as.numeric(status_all == "TP53_mutant"), colnames(mat0))
    } else {
      NULL
    }

    bi_all <- get_batch_factor(ds_dir, colnames(mat0))
    batch_all <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL
    purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
    sa_all <- get_sex_age_covariates(ds_dir, colnames(mat0))
    sa_all_limma <- sa_all
    if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing", "age_z_imputed") %in% colnames(sa_all))) {
      keep_cols <- intersect(c("sex", "age", "age_missing", "age_z_imputed"), colnames(sa_all))
      sa_all_limma <- sa_all[, keep_cols, drop = FALSE]
    }

    ## 4a) each subunits（exclude itself ）
    for (su in present_sub) {
      run_predictor_analyses(
        predictor_name = su,
        predictor_vec = mat0[su, ],
        exclude_genes = su,
        ds_id = ds_id, ds_dir = ds_dir,
        mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
        out_root = out_root,
        genesets_by_group = genesets_by_group,
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
        tp53_num_all = tp53_num_all, is_ALL = is_ALL
      )
    }

    ## 4b) CSN complex score ; all CSN members are excluded from ranking.
    csn_score <- build_csn_score_safe(
      mat0,
      subunits = present_sub, combine_7AB = TRUE,
      min_members = 5L, pca_min_samples = 10L
    )

    if (sum(is.finite(csn_score)) >= (2 * min_per_group)) {
      run_predictor_analyses(
        predictor_name = "CSN_SCORE",
        predictor_vec = csn_score,
        exclude_genes = present_sub,
        ds_id = ds_id, ds_dir = ds_dir,
        mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
        out_root = out_root,
        genesets_by_group = genesets_by_group,
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
        tp53_num_all = tp53_num_all, is_ALL = is_ALL
      )
    } else {
      log_msg("  [CSN_SCORE] non NA samples not enough，skip")
    }

    ## 4c) RESIDUAL_<SU> (Refactored to One-Step Analysis: Gene ~ Subunit + CSN_SCORE + Covariates)
    min_n_resid <- min_per_group
    csn_nonNA <- sum(is.finite(csn_score))
    if (csn_nonNA < min_n_resid) {
      log_msg("  [RESIDUAL] CSN score non NA samples not enough（%d < %d），skip residual_*", csn_nonNA, min_n_resid)
    } else {
      # Prepare CSN_SCORE as an extra covariate
      ex_df <- data.frame(CSN_SCORE = csn_score, row.names = names(csn_score))

      for (su in present_sub) {
        # Check if subunit has enough data
        if (sum(is.finite(mat0[su, ])) < (2 * min_per_group)) {
          log_msg("  [RESIDUAL_%s] non NA samples not enough，skip", su)
          next
        }

        # One-Step Analysis:
        # Predictor = Raw Subunit Expression (mat0[su, ])
        # Covariate = CSN_SCORE (passed via extra_covars_df)
        # Limma Model: Gene ~ Subunit + CSN_SCORE + Batch + Purity + Sex + Age ...
        run_predictor_analyses(
          predictor_name = paste0("RESIDUAL_", su), # Keep naming for compatibility
          predictor_vec = mat0[su, ], # Use RAW subunit expression
          exclude_genes = su,
          ds_id = ds_id, ds_dir = ds_dir,
          mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
          out_root = out_root,
          genesets_by_group = genesets_by_group,
          batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
          tp53_num_all = tp53_num_all, is_ALL = is_ALL,
          extra_covars_df = ex_df # Pass CSN_SCORE as covariate
        )
      }
    }
  }

  ## 5) Each version is compiled separately
  present_sub <- intersect(csn_subunits, rownames(mat0))
  sum_units <- c(present_sub, "CSN_SCORE", paste0("RESIDUAL_", present_sub))


  # BatchAdj：limma_t_cont  + interaction
  for (grp_name in names(genesets_by_group)) {
    out_root_coll <- file.path(GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), ds_id, basename(out_root))
    ver_root <- file.path(GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root))
    summarize_all_groups(
      out_root = ver_root,
      csn_subunits = sum_units,
      genesets_by_group = setNames(list(genesets_by_group[[grp_name]]), grp_name),
      stat_tags = {
        st <- c("GSEA_limma_t_cont", "GSEA_limma_interaction")
        st
      }
    )
  }

  invisible(NULL)
}


## Examine which PLEX have been retrieved so far (raw vs. cleaned).
inspect_plex <- function(ds_dir, col = "TMT_PLEX",
                         pipe_policy = BATCH_PIPE_POLICY,
                         min_per_level = BATCH_MIN_PER_LEVEL) {
  # Read the actual samples used (based on the protein matrix).
  mat0 <- load_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(mat0)

  # Read the clinical table and align the samples
  meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  stopifnot(file.exists(meta_fp))
  meta <- suppressMessages(readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
  id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"), names(meta))
  stopifnot(length(id_cols) > 0)
  meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]

  # Get the original and cleaned PLEX
  raw <- as.character(meta[[col]])
  clean <- sanitize_batch_levels(raw, pipe_policy = pipe_policy, min_per_level = min_per_level)


  cat("\n==== ", basename(ds_dir), " | column：", col, " ====\n", sep = "")
  cat("[original PLEX levels]：\n")
  print(sort(table(raw), decreasing = TRUE))
  cat("\n[含 '|' 的Original value]：\n")
  print(utils::head(unique(raw[grepl("\\|", raw %||% "")]), 10))

  cat("\n[cleaned PLEX levels]（pipe_policy = ", pipe_policy,
    ", min_per_level = ", min_per_level, ")：\n",
    sep = ""
  )
  print(sort(table(clean, useNA = "ifany"), decreasing = TRUE))

  # Output lookup table (first 10 columns as an example)
  df_map <- data.frame(SAMPLE_ID = sample_ids, raw_plex = raw, clean_plex = as.character(clean))
  cat("\n[lookup table (first 10 columns)]\n")
  print(utils::head(df_map, 10))

  invisible(list(
    raw_counts = sort(table(raw), decreasing = TRUE),
    clean_counts = sort(table(clean, useNA = "ifany"), decreasing = TRUE),
    map = df_map
  ))
}


## ===== Missingness sensitivity（AGE）=====
USE_AGE_MISSING_INDICATOR <- FALSE

## Safe z-score: Estimates mean/sd using only finite values, retains NA.
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

## PAAD: Used to find the median of values cut off by a semicolon.
.median_from_semicolon <- function(x_chr) {
  vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
  vv <- vv[is.finite(vv)]
  if (!length(vv)) {
    return(NA_real_)
  }
  stats::median(vv)
}


## 3) get purity（named numeric，names=sample_ids，0~1）
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
## get sex/age（sample alignment；sex: 0/1；age: z-score）
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

## Main workflow: strictly using patient file's SEX and AGE to generate covariates (sex: 0/1; age: z-score)
get_sex_age_covariates <- function(ds_dir, sample_ids) {
  samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")

  if (!file.exists(samp_fp) || !file.exists(pat_fp)) {
    log_msg("[covars] The sample/patient file is missing; sex/age is replaced with NA.")
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
    stop("[get_sex_age_covariates] Unable to establish a sample-patient correspondence (missing SAMPLE_ID/PATIENT_ID)")
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

  ## AGE: Main analysis without imputation
  age_col <- intersect(c("AGE", "AGE_AT_DIAGNOSIS", "AGE_AT_INDEX", "AGE_YEARS"), names(pat))[1]
  age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  age_missing <- age_z_imputed <- NULL
  if (!is.na(age_col)) {
    v <- suppressWarnings(as.numeric(pat[[age_col]]))
    names(v) <- as.character(pat[[pid_pat]])
    age_raw <- unname(v[map_pt[sample_ids]])
    age[] <- .z_no_impute(age_raw) # NA keep
    if (!exists("USE_AGE_MISSING_INDICATOR", inherits = FALSE)) USE_AGE_MISSING_INDICATOR <- FALSE
    if (isTRUE(USE_AGE_MISSING_INDICATOR)) {
      age_missing <- as.numeric(is.na(age_raw))
      age_z_imputed <- .zscore(age_raw)
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
  df, # Candidate covariates；rownames=sample_order
  sample_order,
  label = "covars",
  y = NULL, # contim=nuous predictor；can be NULL
  min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
  max_abs_cor = 0.30,
  min_pairs = 20L
) {
  logf <- function(...) if (exists("log_msg", mode = "function")) try(log_msg(...), silent = TRUE)
  ## Drop Cause Collector
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


  # 1) Align the sample order
  if (is.null(rownames(df))) {
    logf("  [covars-%s] df no rownames → skip", label)
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


  # 2)Use data.frame ; retain factors
  df <- as.data.frame(df, stringsAsFactors = TRUE, check.names = FALSE)

  # 3) coverage: Unnamed fields use the global minimum threshold.
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  global_min <- if (length(min_cov_named)) min(min_cov_named, na.rm = TRUE) else 0
  keep_cov <- rep(TRUE, ncol(df))
  names(keep_cov) <- colnames(df)

  for (cn in colnames(df)) {
    v <- df[[cn]]

    cover <- if (is.numeric(v)) mean(is.finite(v)) else mean(!is.na(v))

    thr <- if (!is.null(names(min_cov_named)) && (cn %in% names(min_cov_named))) {
      min_cov_named[[cn]]
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

  ## Reinitialize keep_cov.
  keep_cov <- rep(TRUE, ncol(df))
  names(keep_cov) <- colnames(df)


  # 4) Biology-related correlation filtering
  if (!is.null(y)) {
    y <- suppressWarnings(as.numeric(y))
    ## Do not perform rho gate on these covariates
    skip_rho_gate <- c("purity", "age", "sex", "CSN_SCORE")

    for (cn in colnames(df)) {
      v <- df[[cn]]

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


  # 5) Return: row order = sample_order; used for model.matrix expansion
  stopifnot(nrow(df) == length(so), identical(rownames(df), so))

  attr(df, "drop_log") <- if (length(drop_log)) do.call(rbind, drop_log) else NULL
  df
}


# Utility Functions
if (!exists(".norm_names", inherits = FALSE)) {
  .norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))
}
if (!exists("%||%", inherits = FALSE)) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
if (!exists(".to01", inherits = FALSE)) {
  .to01 <- function(v) {
    v <- suppressWarnings(as.numeric(v))

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

#  tool: Rearrange the number of samples in each layer of a batch into a string.
.format_batch_sizes <- function(fac) {
  if (is.null(fac)) {
    return(NA_character_)
  }
  tb <- sort(table(fac), decreasing = TRUE)
  paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
}


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

  # Create a corresponding sample→patient profile
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
## —— Add: policy/batch_min/limma_min to the log and output fields
## ==============================
audit_one_dataset_sa_batch <- function(ds_dir, pipe_policy = "NA", min_per_level = 2) {
  ds_id <- basename(ds_dir)
  # 1) get protein matrix
  m <- load_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(m)

  # 2) get SEX/AGE from patient files
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

      # AGE
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

  # 4) batch check
  bi <- get_batch_factor(ds_dir, sample_ids,
    pipe_policy   = pipe_policy,
    min_per_level = min_per_level
  )
  batch_col <- if (!is.null(bi)) bi$name else "NONE"
  batch_levels <- if (!is.null(bi)) nlevels(bi$fac) else 0L
  batch_nonNA <- if (!is.null(bi)) mean(!is.na(bi$fac)) * 100 else NA_real_
  batch_sizes <- if (!is.null(bi)) .format_batch_sizes(bi$fac) else NA_character_

  # Print a summary line
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


## ===== DO_AUDIT switch=====
DO_AUDIT <- TRUE

if (DO_AUDIT) {
  results <- list()
  for (ds in names(dataset_dirs)) {
    ds_dir <- dataset_dirs[[ds]]
    if (!dir.exists(ds_dir)) {
      log_msg("skip：can not find the folder %s", ds_dir)
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

msigdbr_version <- write_geneset_manifest(genesets_by_group)
yaml::write_yaml(list(
  seed = 1234,
  min_frac_complete = min_frac_complete,
  minSize = minSize, maxSize = maxSize, fgsea_eps = fgsea_eps,
  min_per_group = min_per_group,
  MAKE_PLOTS = MAKE_PLOTS,
  datasets_root = datasets_root,
  dataset_ids = dataset_ids,
  strata = strata,
  geneset_groups_selected = GENESET_GROUPS_TO_RUN,
  msigdbr_version = msigdbr_version
), file.path("run_info", "run_manifest.yml"))


## Which one to start from? (Can be changed)
start_from <- "brca_cptac_2020"
ord <- dataset_ids
ix <- match(start_from, ord)
if (is.na(ix)) stop(sprintf("'%s' not in dataset_ids ", start_from))
dataset_dirs_run <- dataset_dirs[ord[ix:length(ord)]]

exists("coerce_covariates_safely")
getAnywhere("coerce_covariates_safely")
exists("opt")
getAnywhere("opt")


## ---- Optional: only run specific strata ----
# only_strata <- c("TP53_mutant")

## Tool: If only_strata is not set, it will run all stratum by default.
.should_run <- function(tag) {
  if (!exists("only_strata") || is.null(only_strata) || !length(only_strata)) {
    return(TRUE)
  }
  tag %in% only_strata
}


## =========================================================
## Process the datasets to be run in this round in sequence (complete the main loop + define the sample sets for each layer)
## =========================================================
for (ds in names(dataset_dirs_run)) {
  ds_dir <- dataset_dirs_run[[ds]]
  log_msg("== Start processing Dataset：%s ==", ds)

  ## ---- [RUN_PASSES] ----
  options(csn.run_passes = c("BatchAdj"))
  ## ---- [END RUN_PASSES]] ----


  ## Protein Matrix & TP53 Status
  mat0_full <- load_matrix_from_dataset_dir(ds_dir)
  tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full))

  ## three stratum
  samples_ALL <- colnames(mat0_full)
  samples_MUT <- names(tp53_status)[tp53_status == "TP53_mutant"]
  samples_WT <- names(tp53_status)[tp53_status == "TP53_wild_type"]

  ## The output root directory for each stratum
  base_tp53_root <- if (is.null(GSEA_PREFIX)) file.path(ds_dir, "csn_gsea_results_TP53") else file.path(GSEA_PREFIX, ds, "csn_gsea_results_TP53")


  ## Run three strata in sequence
  if (.should_run("ALL")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full,
      sample_keep = samples_ALL,
      out_root = file.path(base_tp53_root, "ALL"),
      genesets_by_group = genesets_by_group
    )
  }

  if (.should_run("TP53_mutant")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full,
      sample_keep = samples_MUT,
      out_root = file.path(base_tp53_root, "TP53_mutant"),
      genesets_by_group = genesets_by_group
    )
  }

  if (.should_run("TP53_wild_type")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full,
      sample_keep = samples_WT,
      out_root = file.path(base_tp53_root, "TP53_wild_type"),
      genesets_by_group = genesets_by_group
    )
  }

  log_msg("== completed dataset：%s（TP53 stratum → %s）==", ds, base_tp53_root)
}


## ----------  WT vs MT ΔNES aggregation ----------
tp53_delta_nes_aggregate <- function(datasets_root,
                                     dataset_ids = NULL,
                                     versions = c("BatchAdj"),
                                     groups = names(genesets_by_group),
                                     stat_tags = c("GSEA_limma_t_cont")) {
  if (is.null(dataset_ids)) {
    dataset_ids <- list.dirs(datasets_root, full.names = FALSE, recursive = FALSE)
  }

  for (ds in dataset_ids) {
    base_dir <- file.path(datasets_root, ds, "csn_gsea_results_TP53")
    for (ver in versions) {
      # Automatically detect subunits: Read subdirectory names under the WT/MT folder
      subunits <- unique(unlist(lapply(groups, function(g) {
        b1 <- file.path(GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(g), ver, ds, "TP53_wild_type")
        b2 <- file.path(GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(g), ver, ds, "TP53_mutant")
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
            fp_wt <- file.path(GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection", grp_safe, ver, ds, "TP53_wild_type", su, paste0(st, ".csv"))
            fp_mt <- file.path(GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection", grp_safe, ver, ds, "TP53_mutant", su, paste0(st, ".csv"))
            if (!file.exists(fp_wt) || !file.exists(fp_mt)) next
            wt <- tryCatch(data.table::fread(fp_wt, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
            mt <- tryCatch(data.table::fread(fp_mt, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
            if (is.null(wt) || is.null(mt) || !nrow(wt) || !nrow(mt)) next
            wt <- setNames(as.data.frame(wt), tolower(names(wt)))
            mt <- setNames(as.data.frame(mt), tolower(names(mt)))
            # Standardize column name
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
            out_dir <- file.path(GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection", grp_safe, ver, ds, "DeltaWT_MT", su)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            out_csv <- file.path(out_dir, paste0(st, "_deltaWTMT.csv"))
            data.table::fwrite(df, out_csv)
          }
        }
      }
      message(sprintf("[ΔNES] %s/%s completed", ds, ver))
    }
  }
  invisible(TRUE)
}

## =========================================================
## Generating meta-FDR across datasets (Stouffer → BH)
## Conditions: The above for loop has been completed and all datasets/stratifications/statistics have been implemented.
## =========================================================
meta_fdr_stouffer(
  dataset_dirs = dataset_dirs_run,
  strata = strata,
  stat_tags = c(
    "GSEA_limma_t_cont",
    "GSEA_limma_interaction"
  ),
  groups = names(genesets_by_group),
  out_root = if (is.null(GSEA_PREFIX)) {
    "csn_gsea_pan_summary_TP53/meta_fdr"
  } else {
    file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr")
  }
)


## [1b] WT vs MT ΔNES aggregation
try(tp53_delta_nes_aggregate(
  datasets_root = datasets_root,
  dataset_ids   = names(dataset_dirs_run),
  versions      = c("BatchAdj"),
  groups        = names(genesets_by_group),
  stat_tags     = c("GSEA_limma_t_cont")
), silent = TRUE)


## =========================================================
##  CSN subunits pairwise correlations per dataset
## =========================================================

if (!exists("opt", mode = "function")) {
  opt <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
}

.safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

# residualize
residualize_to_covars <- function(y, batch = NULL, covars = NULL, min_n = 8L) {
  stopifnot(!is.null(names(y)))
  sam <- names(y)
  DF <- data.frame(row.names = sam, check.names = FALSE)
  if (!is.null(batch)) DF[["batch"]] <- droplevels(batch[sam])
  if (!is.null(covars)) {
    C <- as.data.frame(covars, check.names = FALSE)
    rn <- rownames(C)
    if (is.null(rn)) stop("[residualize_to_covars] covars need rownames=sample")
    C <- C[sam, , drop = FALSE]

    if (exists("coerce_covariates_safely", mode = "function")) C <- coerce_covariates_safely(C)
    for (cn in colnames(C)) DF[[cn]] <- C[[cn]]
  }
  yv <- suppressWarnings(as.numeric(y[sam]))
  # Clear all NA/Constant/Single Level columns
  if (ncol(DF) > 0) {
    all_na <- vapply(DF, function(v) all(is.na(v)), logical(1))
    if (any(all_na)) DF <- DF[, !all_na, drop = FALSE]
    is_const <- vapply(DF, function(v) {
      vv <- v[!is.na(v)]
      if (!length(vv)) TRUE else if (is.factor(v)) nlevels(droplevels(vv)) <= 1 else stats::var(as.numeric(vv)) == 0
    }, logical(1))
    if (any(is_const)) DF <- DF[, !is_const, drop = FALSE]
  }
  ok <- is.finite(yv) & (if (ncol(DF) == 0) TRUE else stats::complete.cases(DF))
  if (sum(ok) < min_n) {
    return(setNames(rep(NA_real_, length(y)), names(y)))
  }
  des <- if (ncol(DF) == 0) model.matrix(~1) else model.matrix(~ 1 + ., data = DF[ok, , drop = FALSE])
  fit <- lm.fit(x = des, y = yv[ok])
  out <- setNames(rep(NA_real_, length(y)), names(y))
  out[ok] <- fit$residuals
  out
}

# single dataset： pairwise correlation
csn_pairwise_correlation_one_ds <- function(
  ds_id, ds_dir,
  out_root = "CSN_subunits_correlation_coefficient",
  subunits = csn_subunits,
  min_pairs = 10L,
  dpi = 600,
  width = 4.5, height = 4.5, # inches
  point_size = 1.6
) {
  dir.create(out_root, showWarnings = FALSE)
  ds_out <- file.path(out_root, ds_id)
  dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)

  # read matrix
  mat0 <- load_matrix_from_dataset_dir(ds_dir)
  present <- intersect(subunits, rownames(mat0))
  if (length(present) < 2) {
    log_msg("[pairwise] %s: available CSN subunits < 2，skip", ds_id)
    return(invisible(NULL))
  }
  sam_all <- colnames(mat0)

  # covariate candidte（purity/sex/age）
  build_covars_df <- function(ds_id, ds_dir, sample_ids) {
    pur <- get_purity_covariate(ds_id, ds_dir, sample_ids)
    sa <- get_sex_age_covariates(ds_dir, sample_ids)
    df <- data.frame(
      purity = as.numeric(pur),
      sex = as.numeric(sa[, "sex"]),
      age = as.numeric(sa[, "age"]),
      row.names = sample_ids,
      check.names = FALSE
    )
    if (exists("coerce_covariates_safely", mode = "function")) {
      df <- coerce_covariates_safely(df)
    }
    df
  }
  cov0 <- build_covars_df(ds_id, ds_dir, sam_all)
  cov_raw <- cov0


  bi <- get_batch_factor(ds_dir, sam_all)
  batch <- NULL
  if (!is.null(bi) && !is.null(bi$fac)) {
    b <- bi$fac
    if (is.null(names(b))) {
      if (length(b) == length(sam_all)) {
        b <- stats::setNames(b, sam_all)
      } else {
        b <- NULL
      }
    }
    if (!is.null(b)) batch <- droplevels(b[sam_all])
  }


  keep_rows <- rowMeans(is.finite(mat0)) >= opt("min_frac_complete", 0.75)
  Mimp <- mat0[keep_rows, , drop = FALSE]
  set.seed(1234)
  Mimp <- imputeLCMD::impute.MinProb(Mimp, q = 0.01)


  pairs <- utils::combn(present, 2, simplify = FALSE)


  all_rows <- list()

  # Helper for drawing output (not executed by default; controlled by MAKE_PLOTS)
  .save_scatter <- function(df_xy, title, out_base) {
    gp <- ggplot2::ggplot(df_xy, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(size = point_size, alpha = 0.8) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
      ggplot2::labs(title = title, x = df_xy$gx[1], y = df_xy$gy[1]) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA)
      )
    if (isTRUE(opt("MAKE_PLOTS", FALSE))) {
      ggplot2::ggsave(paste0(out_base, ".tiff"), gp, width = width, height = height, dpi = dpi, compression = "lzw")
      ggplot2::ggsave(paste0(out_base, ".png"), gp, width = width, height = height, dpi = dpi)
      ggplot2::ggsave(paste0(out_base, ".jpg"), gp, width = width, height = height, dpi = dpi)
    }
  }

  # Main loop: Each pair of subunits
  for (p in pairs) {
    gx <- p[1]
    gy <- p[2]
    x <- as.numeric(mat0[gx, sam_all])
    names(x) <- sam_all
    y <- as.numeric(mat0[gy, sam_all])
    names(y) <- sam_all
    ok <- stats::complete.cases(x, y)
    if (sum(ok) < min_pairs) next

    ## ---- version 1: NoCovariate ----
    rr <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "pearson"))
    fit <- stats::lm(y ~ x, data = data.frame(x = x[ok], y = y[ok]))
    all_rows[[length(all_rows) + 1L]] <- data.frame(
      dataset = ds_id, version = "NoCovariate",
      gene_x = gx, gene_y = gy, n = sum(ok),
      pearson_r = as.numeric(rr$estimate), pearson_p = rr$p.value,
      R2 = summary(fit)$r.squared,
      slope = unname(stats::coef(fit)[["x"]]),
      intercept = unname(stats::coef(fit)[["(Intercept)"]]),
      stringsAsFactors = FALSE, check.names = FALSE
    )


    .save_scatter(
      df_xy = data.frame(x = x[ok], y = y[ok], gx = gx, gy = gy),
      title = sprintf("%s | %s vs %s (NoCovariate)", ds_id, gx, gy),
      out_base = file.path(ds_out, sprintf("%s_vs_%s_NoCovariate", .safe_fs(gx), .safe_fs(gy)))
    )

    ## ---- version 3: BatchAdj_covars----
    cov_ba <- cov_raw
    batch_for_res <- NULL

    if (!is.null(batch)) {
      # 有 batch
      batch_for_res <- batch
    }

    xb <- residualize_to_covars(x, batch = batch_for_res, covars = cov_ba)
    yb <- residualize_to_covars(y, batch = batch_for_res, covars = cov_ba)
    ok3 <- stats::complete.cases(xb, yb)

    if (sum(ok3) >= min_pairs) {
      rr3 <- suppressWarnings(stats::cor.test(xb[ok3], yb[ok3], method = "pearson"))

      fit3 <- stats::lm(yb ~ xb, data = data.frame(xb = xb[ok3], yb = yb[ok3]))
      all_rows[[length(all_rows) + 1L]] <- data.frame(
        dataset = ds_id, version = "BatchAdj_covars",
        gene_x = gx, gene_y = gy, n = sum(ok3),
        pearson_r = as.numeric(rr3$estimate), pearson_p = rr3$p.value,
        R2 = summary(fit3)$r.squared,
        slope = unname(stats::coef(fit3)[["xb"]]),
        intercept = unname(stats::coef(fit3)[["(Intercept)"]]),
        stringsAsFactors = FALSE, check.names = FALSE
      )
      .save_scatter(
        df_xy = data.frame(x = xb[ok3], y = yb[ok3], gx = gx, gy = gy),
        title = sprintf("%s | %s vs %s (BatchAdj_covars)", ds_id, gx, gy),
        out_base = file.path(ds_out, sprintf("%s_vs_%s_BatchAdj_covars", .safe_fs(gx), .safe_fs(gy)))
      )
    }
  } # end for pairs

  # Summary and Output
  if (!length(all_rows)) {
    log_msg("[pairwise] %s: No available pairs (insufficient samples or too many NAs)", ds_id)
    return(invisible(NULL))
  }
  RES <- do.call(rbind, all_rows)
  # BH correction based on dataset + version
  RES$pearson_padj <- ave(RES$pearson_p, interaction(RES$dataset, RES$version, drop = TRUE),
    FUN = function(p) stats::p.adjust(p, method = "BH")
  )


  out_csv <- file.path(ds_out, "pairwise_correlations_all_versions.csv")
  data.table::fwrite(RES, out_csv)

  out_xlsx <- file.path(ds_out, "pairwise_correlations_all_versions.xlsx")
  wb <- openxlsx::createWorkbook()
  for (ver in unique(RES$version)) {
    openxlsx::addWorksheet(wb, ver)
    openxlsx::writeData(wb, ver, RES[RES$version == ver, ], withFilter = TRUE)
  }
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

  log_msg("[pairwise] %s: completed. CSV: %s；XLSX: %s", ds_id, basename(out_csv), basename(out_xlsx))
  invisible(RES)
}

# Run all datasets
run_csn_subunits_pairwise_correlations <- function(
  dataset_dirs_map = NULL,
  out_root = "CSN_subunits_correlation_coefficient"
) {
  if (is.null(dataset_dirs_map)) {
    if (exists("dataset_dirs_run")) {
      dataset_dirs_map <- get("dataset_dirs_run", inherits = TRUE)
    } else if (exists("dataset_dirs")) {
      dd <- get("dataset_dirs", inherits = TRUE)
      dataset_dirs_map <- dd[dir.exists(dd) & file.exists(file.path(dd, "data_protein_quantification.txt"))]
    } else {
      stop("cannot find dataset_dirs_run or dataset_dirs")
    }
  }
  dir.create(out_root, showWarnings = FALSE)
  for (ds in names(dataset_dirs_map)) {
    try(csn_pairwise_correlation_one_ds(ds_id = ds, ds_dir = dataset_dirs_map[[ds]], out_root = out_root), silent = FALSE)
  }
  invisible(TRUE)
}

## =========================================================
## Direct execution: Pairwise correlation of CSN subunits
## - Scatter plots are not output by default; to output .tiff/.png/.jpg simultaneously, change MAKE_PLOTS <- TRUE
## =========================================================
MAKE_PLOTS <- FALSE

run_csn_subunits_pairwise_correlations(
  dataset_dirs_map = if (exists("dataset_dirs_run")) dataset_dirs_run else NULL,
  out_root = "CSN_subunits_correlation_coefficient"
)


## =========================================================
##  Correlation coefficient heatmaps from pairwise CSVs
##   - For each dataset under CSN_subunits_correlation_coefficient/<ds>/
##   - Read: pairwise_correlations_all_versions.csv
## - Make versions: NoCovariate / BatchAdj_covars
## =========================================================
suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr", repos = "https://cloud.r-project.org")
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

CSN_SUB_ORDER <- c("GPS1", "COPS2", "COPS3", "COPS4", "COPS5", "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9")
CELL_BLUE <- "#3B4CC0"
CELL_WHITE <- "#F7F7F7"
CELL_RED <- "#B40426"


.safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))


.build_corr_grid <- function(csv_path, version, sub_order = CSN_SUB_ORDER,
                             triangle = getOption("csn_heatmap_triangle", "full"),
                             method = getOption("csn_heatmap_method", "pearson")) {
  stopifnot(file.exists(csv_path))
  df <- suppressMessages(readr::read_csv(csv_path, show_col_types = FALSE))
  req <- c("version", "gene_x", "gene_y", "pearson_r", "pearson_padj")
  if (!all(req %in% names(df))) stop("CSV missing columns：", paste(setdiff(req, names(df)), collapse = ", "))


  rcol <- "pearson_r"
  pcol <- "pearson_padj"


  dfv <- df %>%
    dplyr::filter(.data$version == !!version) %>%
    dplyr::mutate(
      gene_x = stringr::str_trim(gene_x),
      gene_y = stringr::str_trim(gene_y)
    )
  tmp_r <- dfv[[rcol]]
  tmp_p <- dfv[[pcol]]
  dfv <- dfv %>% dplyr::transmute(gene_x, gene_y, r = tmp_r, padj = tmp_p)

  present <- intersect(sub_order, unique(c(dfv$gene_x, dfv$gene_y)))
  if (length(present) == 0L) stop("This version has no CSN subunits in this dataset.")

  grid <- tidyr::expand_grid(gene_x = present, gene_y = present) %>%
    dplyr::mutate(
      key     = paste(gene_x, gene_y, sep = "|"),
      key_rev = paste(gene_y, gene_x, sep = "|")
    )

  dfv_key <- dfv %>%
    dplyr::mutate(key = paste(gene_x, gene_y, sep = "|")) %>%
    dplyr::select(key, r, padj)

  grid2 <- grid %>%
    dplyr::left_join(dfv_key, by = "key") %>%
    dplyr::left_join(dplyr::rename(dfv_key, key_rev = key, r2 = r, padj2 = padj), by = "key_rev") %>%
    dplyr::mutate(
      r    = dplyr::coalesce(r, r2, ifelse(gene_x == gene_y, 1, NA_real_)),
      padj = dplyr::coalesce(padj, padj2, ifelse(gene_x == gene_y, NA_real_, NA_real_))
    ) %>%
    dplyr::select(gene_x, gene_y, r, padj) %>%
    dplyr::mutate(
      gene_x = factor(gene_x, levels = present),
      gene_y = factor(gene_y, levels = present),
      signif = !is.na(padj) & padj < 0.05 & abs(r) < 0.999999
    )

  tri <- match.arg(tolower(triangle), c("full", "lower", "upper"))
  if (tri != "full") {
    ix <- as.integer(grid2$gene_x)
    iy <- as.integer(grid2$gene_y)
    if (tri == "lower") grid2 <- grid2[iy >= ix, , drop = FALSE]
    if (tri == "upper") grid2 <- grid2[iy <= ix, , drop = FALSE]
  }
  grid2
}


.plot_corr_heatmap_save <- function(df_grid, title, out_base,
                                    width = 6.5, height = 6.5, dpi = 600) {
  # Using blue-white-red diverging colors, limited to -1~1; NA pale gray
  p <- ggplot(df_grid, aes(x = gene_x, y = gene_y, fill = r)) +
    geom_tile(color = "white", size = 0.3, na.rm = FALSE) +
    #  black spots
    geom_point(
      data = subset(df_grid, signif), aes(x = gene_x, y = gene_y),
      inherit.aes = FALSE, shape = 16, size = 2.0, color = "black", alpha = 0.9
    ) +
    scale_fill_gradient2(
      low = CELL_BLUE, mid = CELL_WHITE, high = CELL_RED,
      midpoint = 0, na.value = "grey90", limits = c(-1, 1), oob = scales::squish,
      name = "Pearson r"
    ) +
    scale_x_discrete(position = "top") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.text.x = element_text(angle = 45, hjust = 0),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(),
      legend.position = "right"
    )


  ggsave(paste0(out_base, ".tiff"), p, width = width, height = height, dpi = dpi, compression = "lzw")
  ## ggsave(paste0(out_base, ".png"),  p, width = width, height = height, dpi = dpi)
  ## ggsave(paste0(out_base, ".jpg"),  p, width = width, height = height, dpi = dpi)
  ## ggsave(paste0(out_base, ".pdf"),  p, width = width, height = height)  # 向量化
}

# Find the CSV files for all datasets and plot them.
run_csn_corrcoef_plots <- function(out_root = "CSN_subunits_correlation_coefficient") {
  if (!dir.exists(out_root)) {
    message("[corrplot] Root directory not found：", out_root)
    return(invisible(FALSE))
  }
  csvs <- list.files(out_root, pattern = "^pairwise_correlations_all_versions\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(csvs)) {
    message("[corrplot] cannot find any pairwise_correlations_all_versions.csv")
    return(invisible(FALSE))
  }

  versions <- c("NoCovariate", "BatchAdj_covars")
  for (csv in csvs) {
    ds_dir <- dirname(csv)
    ds_id <- basename(ds_dir)
    plots_dir <- file.path(ds_dir)


    for (ver in versions) {
      grid_p <- try(.build_corr_grid(csv, ver, CSN_SUB_ORDER,
        triangle = getOption("csn_heatmap_triangle", "full"),
        method   = "pearson"
      ), silent = TRUE)
      if (!inherits(grid_p, "try-error")) {
        title_p <- sprintf("%s — %s (Pearson r)", ds_id, ver)
        out_base_p <- file.path(plots_dir, sprintf("%s_corrcoef_%s", .safe_fs(ds_id), .safe_fs(ver)))
        try(.plot_corr_heatmap_save(grid_p, title_p, out_base_p), silent = FALSE)
      }
    }
    message("[corrplot] compleyed：", ds_id)
  }
  invisible(TRUE)
}

## ----  Execute drawing ----
run_csn_corrcoef_plots(out_root = "CSN_subunits_correlation_coefficient")


## =========================================================
## Batch plotting GSEA Summary CSV of a single dataset/stratum into a heatmap
## =========================================================

# --- Working directory (portable) ---
wd_candidates <- c(
  Sys.getenv("CSN_CPTAC_ROOT", unset = NA_character_),
  "C:/Users/danny/Documents/R_project/CSN_CPTAC",
  "C:/Users/cmuh/Documents/YenTing_document/CSN_CPTAC"
)
wd_candidates <- wd_candidates[!is.na(wd_candidates) & wd_candidates != ""]
hit <- wd_candidates[dir.exists(wd_candidates)][1]
if (!is.na(hit)) setwd(hit)
getwd()


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})


CELL_OG_PALETTES <- list(
  cell_og1 = c(neg = "#1B9E77", mid = "#FFFFFF", pos = "#D95F02"),
  cell_og2 = c(neg = "#2E7D32", mid = "#FFFFFF", pos = "#FB8C00"),
  cell_og3 = c(neg = "#009E73", mid = "#FFFFFF", pos = "#E69F00"),
  cell_og4 = c(neg = "#006D2C", mid = "#FFFFFF", pos = "#E6550D"),
  cell_og5 = c(neg = "#238B45", mid = "#FFFFFF", pos = "#F16913")
)


.get_interaction_palette <- function() {
  key <- getOption("csn_interaction_palette", "cell_og2")
  pal <- CELL_OG_PALETTES[[key]]
  if (is.null(pal)) pal <- CELL_OG_PALETTES[[1]]
  pal
}


.make_heatmap_plot_with_yorder_colored <- function(csv_file, y_order, palette = .get_interaction_palette()) {
  df_raw <- suppressMessages(.read_csv_safe(csv_file))
  path_candidates <- c("pathway", "Pathway", "term", "Term", "gs_name", "NAME", "set", "Set")
  path_col <- intersect(path_candidates, names(df_raw))[1]
  if (is.na(path_col)) stop("cannot find pathway ：", paste(path_candidates, collapse = ", "))

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
    if (!length(nes_cols)) stop("wide list：cannot find NES_* ")
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


  df_long <- df_long %>% dplyr::filter(.data$pathway %in% y_order)
  y_levels <- y_order[y_order %in% df_long$pathway]
  if (!length(y_levels)) stop("ordered version：The pathway order does not overlap with ALL.")
  df_long <- df_long %>% dplyr::mutate(pathway = factor(.data$pathway, levels = rev(y_levels)))


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


  if (!file.exists(p) && .Platform$OS.type == "windows") {
    p2 <- gsub("/", "\\\\", p, fixed = TRUE)
    if (file.exists(p2)) p <- p2
  }
  if (!file.exists(p)) stop("file not exiests：", p)

  fi <- suppressWarnings(file.info(p))
  if (isTRUE(fi$isdir)) stop("The target is not the file.：", p)
  if (!is.na(fi$size) && fi$size == 0) stop("empty file（size = 0）：", p)


  out <- try(suppressMessages(readr::read_csv(p, show_col_types = FALSE, progress = FALSE)), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(out)
  }


  out <- try(suppressMessages(readr::read_csv(
    p,
    locale = readr::locale(encoding = "UTF-8-BOM"),
    show_col_types = FALSE, progress = FALSE
  )), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(out)
  }


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


  out <- try(utils::read.csv(p, check.names = FALSE), silent = TRUE)
  if (inherits(out, "try-error")) stop(out)
  tibble::as_tibble(out)
}


.safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

# X-axis order
.pred_order_all <- c(
  "NES_CSN_SCORE", "NES_GPS1",
  "NES_COPS2", "NES_COPS3", "NES_COPS4", "NES_COPS5", "NES_COPS6",
  "NES_COPS7A", "NES_COPS7B", "NES_COPS8", "NES_COPS9",
  "NES_RESIDUAL_GPS1",
  "NES_RESIDUAL_COPS2", "NES_RESIDUAL_COPS3", "NES_RESIDUAL_COPS4", "NES_RESIDUAL_COPS5", "NES_RESIDUAL_COPS6",
  "NES_RESIDUAL_COPS7A", "NES_RESIDUAL_COPS7B", "NES_RESIDUAL_COPS8", "NES_RESIDUAL_COPS9"
)


.make_heatmap_plot <- function(csv_file) {
  df_raw <- suppressMessages(.read_csv_safe(csv_file))


  path_candidates <- c("pathway", "Pathway", "term", "Term", "gs_name", "NAME", "set", "Set")
  path_col <- intersect(path_candidates, names(df_raw))[1]
  if (is.na(path_col)) stop("找不到 pathway 欄位：", paste(path_candidates, collapse = ", "))


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
    if (!length(nes_cols)) stop("Wide table format：cannot find NES_* 。")

    df_wide <- df_raw %>% rename(pathway = all_of(path_col))
    nes_map <- tibble(nes_col = nes_cols, key = sub("^NES_", "", nes_cols))
    padj_map <- tibble(padj_col = padj_cols, key = sub("^padj_", "", padj_cols))
    pair_map <- nes_map %>% left_join(padj_map, by = "key")

    df_list <- lapply(seq_len(nrow(pair_map)), function(i) {
      nes_c <- pair_map$nes_col[i]
      padj_c <- pair_map$padj_col[i]
      tibble(
        pathway   = df_wide$pathway,
        predictor = paste0("NES_", pair_map$key[i]),
        NES       = suppressWarnings(as.numeric(df_wide[[nes_c]])),
        padj      = if (!is.na(padj_c)) suppressWarnings(as.numeric(df_wide[[padj_c]])) else NA_real_
      )
    })
    df_long <- bind_rows(df_list)
  }


  present_preds <- intersect(.pred_order_all, unique(df_long$predictor))
  if (!length(present_preds)) stop("The data does not specify any predictors.")
  df_long <- df_long %>% filter(predictor %in% present_preds)

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

  # Y-axis sorting: Prioritize using NES_CSN_SCORE from largest to smallest; otherwise, use the row average.
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
  df_long <- df_long %>% mutate(xpos = unname(pos_map[predictor]))


  L <- max(abs(df_long$NES), na.rm = TRUE)
  if (!is.finite(L) || L == 0) L <- 1
  col_low <- "#053061"
  col_mid <- "#FFFFFF"
  col_high <- "#67001F"

  # Plot the graph (black dots: padj < 0.05)
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


  n_path <- nlevels(df_long$pathway)
  W <- max(8, length(present_preds) * 0.45) # inches
  H <- max(6, n_path * 0.22) # inches
  list(plot = p, width = W, height = H)
}

#
.compute_pathway_order_from_csv <- function(csv_file) {
  df_raw <- suppressMessages(.read_csv_safe(csv_file))
  path_candidates <- c("pathway", "Pathway", "term", "Term", "gs_name", "NAME", "set", "Set")
  path_col <- intersect(path_candidates, names(df_raw))[1]
  if (is.na(path_col)) stop("cannot find pathway ：", paste(path_candidates, collapse = ", "))

  has_long <- all(c("predictor", "NES", "padj") %in% names(df_raw))
  if (has_long) {
    df_long <- df_raw %>%
      dplyr::rename(pathway = dplyr::all_of(path_col)) %>%
      dplyr::mutate(
        predictor = as.character(.data$predictor),
        NES = as.numeric(.data$NES)
      )
  } else {
    nes_cols <- grep("^NES_", names(df_raw), value = TRUE)
    if (!length(nes_cols)) stop("wide table：cannot find NES_* ")
    df_wide <- df_raw %>% dplyr::rename(pathway = dplyr::all_of(path_col))
    nes_map <- tibble::tibble(nes_col = nes_cols, key = sub("^NES_", "", nes_cols))
    df_list <- lapply(seq_len(nrow(nes_map)), function(i) {
      nes_c <- nes_map$nes_col[i]
      tibble::tibble(
        pathway   = df_wide$pathway,
        predictor = paste0("NES_", nes_map$key[i]),
        NES       = suppressWarnings(as.numeric(df_wide[[nes_c]]))
      )
    })
    df_long <- dplyr::bind_rows(df_list)
  }

  present_preds <- intersect(.pred_order_all, unique(df_long$predictor))
  if (!length(present_preds)) stop("The data does not specify any predictors.")
  df_long <- df_long %>% dplyr::filter(.data$predictor %in% present_preds)


  .coll_tok <- sub("^Summary_(.+?)_GSEA.*$", "\\1", basename(csv_file))
  if ("NES_CSN_SCORE" %in% present_preds) {
    csn_tbl <- df_long %>%
      dplyr::filter(.data$predictor == "NES_CSN_SCORE") %>%
      dplyr::select(.data$pathway, .data$NES) %>%
      dplyr::distinct()
    .top_n <- get0("DATASET_HEATMAP_TOP_N", ifnotfound = 25L)
    .bot_n <- get0("DATASET_HEATMAP_BOTTOM_N", ifnotfound = 25L)

    if (!identical(toupper(.coll_tok), "H")) {
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
        csn_tbl <- csn_tbl %>% dplyr::filter(.data$pathway %in% keep)
      }
    } else {
      keep <- csn_tbl$pathway
      df_long <- df_long %>% dplyr::filter(.data$pathway %in% keep)
      csn_tbl <- csn_tbl %>% dplyr::filter(.data$pathway %in% keep)
    }
  }


  if ("NES_CSN_SCORE" %in% present_preds) {
    nes_csn <- df_long %>%
      dplyr::filter(.data$predictor == "NES_CSN_SCORE") %>%
      dplyr::select(.data$pathway, .data$NES) %>%
      dplyr::distinct()
    path_order <- nes_csn %>%
      dplyr::arrange(dplyr::desc(.data$NES)) %>%
      dplyr::pull(.data$pathway)
  } else {
    path_order <- df_long %>%
      dplyr::group_by(.data$pathway) %>%
      dplyr::summarise(m = mean(.data$NES, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(.data$m)) %>%
      dplyr::pull(.data$pathway)
  }
  path_order
}


# # Plotting based on a given y_order (pathway levels)
.make_heatmap_plot_with_yorder <- function(csv_file, y_order) {
  df_raw <- suppressMessages(.read_csv_safe(csv_file))
  path_candidates <- c("pathway", "Pathway", "term", "Term", "gs_name", "NAME", "set", "Set")
  path_col <- intersect(path_candidates, names(df_raw))[1]
  if (is.na(path_col)) stop("cannot find pathway ：", paste(path_candidates, collapse = ", "))

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
    if (!length(nes_cols)) stop("wide table：cannot find NES_* ")
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


  df_long <- df_long %>% dplyr::filter(.data$pathway %in% y_order)
  y_levels <- y_order[y_order %in% df_long$pathway]
  if (!length(y_levels)) stop("The ordered version has no overlap with the pathway order of ALL.")
  df_long <- df_long %>% dplyr::mutate(pathway = factor(.data$pathway, levels = rev(y_levels)))


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

  L <- max(abs(df_long$NES), na.rm = TRUE)
  if (!is.finite(L) || L == 0) L <- 1
  col_low <- "#053061"
  col_mid <- "#FFFFFF"
  col_high <- "#67001F"

  p <- ggplot(df_long, aes(x = xpos, y = pathway, fill = NES)) +
    geom_tile(width = 1, height = 0.9, color = NA) +
    geom_point(
      data = df_long %>% dplyr::filter(is.finite(.data$padj), .data$padj < 0.05),
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

  n_path <- nlevels(df_long$pathway)
  W <- max(8, length(present_preds) * 0.45)
  H <- max(6, n_path * 0.22)
  list(plot = p, width = W, height = H)
}


.save_near_only <- function(p, width, height, csv_file, meta, suffix = "_ordered", dpi = 600) {
  bn <- basename(csv_file)
  dir_near <- dirname(csv_file)
  near_prefix <- paste0(
    "heatmap_", .safe_fs(meta$dataset), "_", .safe_fs(meta$stratum),
    "_", .safe_fs(meta$variant), "_"
  )
  near_base <- file.path(dir_near, paste0(near_prefix, bn, suffix))
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
}


.parse_meta_from_path <- function(csv_file) {
  parts <- strsplit(normalizePath(csv_file, winslash = "/"), "/")[[1]]
  parts_lc <- tolower(parts)
  bn <- basename(csv_file)


  vidx <- which(parts_lc %in% c("batchadj"))
  variant <- if (length(vidx)) parts[vidx[1]] else NA_character_


  stratum <- NA_character_
  if (length(vidx) && vidx[1] - 1 >= 1) {
    cand <- parts[vidx[1] - 1]
    if (tolower(cand) %in% c("all", "tp53_mutant", "tp53_wild_type")) {
      stratum <- cand
    }
  }

  if (is.na(stratum) || !nzchar(stratum)) {
    hit <- which(parts_lc %in% c("all", "tp53_mutant", "tp53_wild_type"))
    if (length(hit)) stratum <- parts[hit[1]]
  }


  dataset <- NA_character_
  hit2 <- which(parts_lc == "csn_gsea_results_tp53")
  if (length(hit2) && hit2[1] > 1) dataset <- parts[hit2[1] - 1]

  list(dataset = dataset, variant = variant, stratum = stratum)
}


# output
.save_both_places <- function(p, width, height, csv_file, meta,
                              collect_root = "single_dataset_GSEA_heatmap", dpi = 600) {
  bn <- basename(csv_file) # e.g. Summary_H_GSEA_limma_t_cont_ALL.csv
  dir_near <- dirname(csv_file)


  near_prefix <- paste0(
    "heatmap_", .safe_fs(meta$dataset), "_", .safe_fs(meta$stratum),
    "_", .safe_fs(meta$variant), "_"
  )
  near_base <- file.path(dir_near, paste0(near_prefix, bn))


  out_dir2 <- file.path(collect_root, .safe_fs(meta$dataset))
  dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)
  collect_base <- file.path(out_dir2, paste0("heatmap_", bn))


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


.find_target_csvs_for_dataset <- function(ds_dir) {
  files <- list.files(
    ds_dir,
    pattern = "^Summary_.+_GSEA_(limma_t(_cont)?|limma_interaction)_(ALL|TP53_mutant|TP53_wild_type)\\.csv$",
    recursive = TRUE, full.names = TRUE
  )
  if (!length(files)) {
    return(files)
  }


  sep_pat <- "[/\\\\]"
  keep <- grepl(paste0(sep_pat, "batchadj", sep_pat), files, ignore.case = TRUE)

  files <- files[keep]


  if (length(files)) {
    finfo <- suppressWarnings(file.info(files))
    files <- files[is.finite(finfo$size) & finfo$size > 0 & !finfo$isdir]
  }
  ## ====================================================

  files
}


#  Main program: Run all datasets
run_gsea_heatmaps_for_all_datasets <- function(dataset_dirs_map = NULL,
                                               collect_root = "single_dataset_GSEA_heatmap") {
  if (is.null(dataset_dirs_map)) {
    if (exists("dataset_dirs_run")) {
      dataset_dirs_map <- get("dataset_dirs_run", inherits = TRUE)
    } else if (exists("dataset_dirs")) {
      dataset_dirs_map <- get("dataset_dirs", inherits = TRUE)
    } else {
      stop("Please provide dataset_dirs_map or define dataset_dirs_run / dataset_dirs in the environment.")
    }
  }
  for (ds in names(dataset_dirs_map)) {
    ds_dir <- dataset_dirs_map[[ds]]
    scan_dir <- if (is.null(GSEA_PREFIX)) ds_dir else GSEA_PREFIX
    csvs <- .find_target_csvs_for_dataset(scan_dir)

    if (!is.null(GSEA_PREFIX)) {
      ds_name <- basename(ds_dir)
      sep_pat <- "[/\\\\]"
      pat_ds <- paste0(sep_pat, ds_name, sep_pat)
      csvs <- csvs[grepl(pat_ds, csvs, ignore.case = TRUE)]
      csvs <- csvs[grepl(paste0(sep_pat, "summary", sep_pat), tolower(csvs))]

      .grp_labels <- names(get0("genesets_by_group", ifnotfound = list()))
      if (length(.grp_labels)) {
        .gs <- unique(vapply(.grp_labels, safe_fs_name, FUN.VALUE = character(1)))
        .pat <- paste(paste0(sep_pat, tolower(.gs), sep_pat), collapse = "|")
        csvs <- csvs[grepl(.pat, tolower(csvs))]
      }
    }

    .cols_ds <- get0("PLOT_DATASET_COLLECTIONS", ifnotfound = NULL)
    if (!is.null(.cols_ds) && length(.cols_ds)) {
      .al <- unique(vapply(.cols_ds, safe_fs_name, FUN.VALUE = character(1)))
      .bas <- basename(csvs)
      .pref <- paste0("Summary_", .al, "_")
      .keep <- Reduce(`|`, lapply(.pref, function(p) startsWith(.bas, p)))
      csvs <- csvs[.keep]
    }
    if (!length(csvs)) {
      message("[heatmap] ", ds, ": No matching Summary CSV found, skipped.")
      next
    }
    for (csv in csvs) {
      meta <- .parse_meta_from_path(csv)

      if (!length(meta$dataset) || is.na(meta$dataset) || meta$dataset == "") {
        meta$dataset <- ds
      }


      if (grepl("GSEA_limma_interaction", basename(csv), fixed = TRUE) &&
        !identical(meta$stratum, "ALL")) {
        next
      }


      bn_lc <- tolower(basename(csv))
      is_interaction_now <- grepl("gsea_limma_interaction", bn_lc, fixed = TRUE)

      h <- try(.make_heatmap_plot(csv), silent = TRUE)
      if (inherits(h, "try-error")) {
        message("[heatmap] 失敗：", csv, " | ", as.character(h))
        next
      }
      .save_both_places(h$plot, h$width, h$height, csv, meta, collect_root = collect_root)
      message("[heatmap] 完成：", csv)


      bn_lc <- tolower(basename(csv))
      if (grepl("gsea_limma_interaction", bn_lc, fixed = TRUE) &&
        identical(meta$stratum, "ALL")) {
        dir_now <- normalizePath(dirname(csv), winslash = "/", mustWork = FALSE)
        dir_cont <- sub("/GSEA_limma_interaction/?$", "/GSEA_limma_t_cont", dir_now, ignore.case = TRUE)
        bn_cont <- gsub("GSEA_limma_interaction", "GSEA_limma_t_cont", basename(csv), ignore.case = TRUE)
        cont_csv <- file.path(dir_cont, bn_cont)

        if (file.exists(cont_csv)) {
          y_order <- try(.compute_pathway_order_from_csv(cont_csv), silent = TRUE)
          if (!inherits(y_order, "try-error")) {
            h2 <- try(.make_heatmap_plot_with_yorder_colored(csv, y_order), silent = TRUE)
            if (!inherits(h2, "try-error")) {
              .save_near_only(h2$plot, h2$width, h2$height, csv, meta, suffix = "_ordered")
              message("[heatmap][interaction-ordered] completed：", csv)
            } else {
              message("[heatmap][interaction-ordered] plotting failed, skipped.：", csv, " | ", as.character(h2))
            }
          } else {
            message("[heatmap][interaction-ordered] Failed to retrieve limma_t_cont ALL in the correct order, skipped:", cont_csv, " | ", as.character(y_order))
          }
        } else {
          message("[heatmap][interaction-ordered] No corresponding limma_t_cont ALL found, skipped:", cont_csv)
        }
      }

      ## === : Ordered versions of TP53_mutant / TP53_wild_type (in y-axis order of ALL) ===

      if (tolower(meta$stratum) %in% c("tp53_mutant", "tp53_wild_type")) {
        bn_lc <- tolower(basename(csv))
        is_cont <- grepl("gsea_limma_t_cont", bn_lc, fixed = TRUE)

        if (is_cont) {
          parts <- strsplit(normalizePath(csv, winslash = "/"), "/")[[1]]
          parts_lc <- tolower(parts)
          sidx <- which(parts_lc %in% c("all", "tp53_mutant", "tp53_wild_type"))
          if (length(sidx)) {
            parts2 <- parts
            parts2[sidx[length(sidx)]] <- "ALL"
            all_csv <- paste(parts2, collapse = "/")
            if (file.exists(all_csv)) {
              y_order <- try(.compute_pathway_order_from_csv(all_csv), silent = TRUE)
              if (!inherits(y_order, "try-error")) {
                h2 <- try(.make_heatmap_plot_with_yorder(csv, y_order), silent = TRUE)
                if (!inherits(h2, "try-error")) {
                  .save_near_only(h2$plot, h2$width, h2$height, csv, meta, suffix = "_ordered")
                  message("[heatmap][ordered] completed：", csv)
                } else {
                  message("[heatmap][ordered] Plotting failed, skipped:", csv, " | ", as.character(h2))
                }
              } else {
                message("[heatmap][ordered] Failed to retrieve ALL in sequence, skipped:", all_csv, " | ", as.character(y_order))
              }
            } else {
              message("[heatmap][ordered] No corresponding ALL file found, skipped:", all_csv)
            }
          }
        }
      }

      try(closeAllConnections(), silent = TRUE)
      invisible(gc(FALSE))
    }
  }
  invisible(TRUE)
}

options(csn_interaction_palette = "cell_og2")

## ---- execute----
run_gsea_heatmaps_for_all_datasets(
  dataset_dirs_map = if (exists("dataset_dirs_run")) dataset_dirs_run else NULL,
  collect_root = if (is.null(GSEA_PREFIX)) "single_dataset_GSEA_heatmap" else file.path(GSEA_PREFIX, "single_dataset_GSEA_heatmap")
)


## =========================================================
##  Summarize meta-FDR (Stouffer → BH) across subunits
## - Input: csn_gsea_pan_summary_TP53/meta_fdr/<STRATUM>/BatchAdj/<SUBUNIT>/<GROUP>/GSEA_limma_t_cont_meta_fdr.csv
## =========================================================

.read_meta_fdr_table <- function(meta_root, stratum, version, subunit,
                                 group_name, stat_tag = "GSEA_limma_t_cont_meta_fdr") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("please install data.table")
  grp <- safe_fs_name(group_name)
  fp <- file.path(meta_root, stratum, version, subunit, grp, paste0(stat_tag, ".csv"))
  if (!file.exists(fp)) {
    return(NULL)
  }
  dt <- tryCatch(data.table::fread(fp, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
  if (is.null(dt)) {
    return(NULL)
  }
  need <- c("pathway", "Z", "padj_meta")

  nms <- tolower(names(dt))
  if (!"pathway" %in% names(dt) && "pathway" %in% nms) names(dt)[match("pathway", nms)] <- "pathway"
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

.merge_subunit_tables_meta <- function(tbl_list) {
  keep <- tbl_list[!vapply(tbl_list, is.null, logical(1))]
  if (!length(keep)) {
    return(NULL)
  }
  out <- Reduce(function(x, y) dplyr::full_join(x, y, by = "pathway"), keep)
  num_cols <- setdiff(names(out), "pathway")
  out[num_cols] <- lapply(out[num_cols], function(z) suppressWarnings(as.numeric(z)))
  out
}

.add_sig_counts_meta <- function(df, alphas = c(0.05, 0.25)) {
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

    pos_mat <- mapply(
      function(p, zc) {
        pv <- df[[p]]
        zv <- df[[zc]]
        as.integer(is.finite(pv) & pv < a & is.finite(zv) & zv > 0)
      },
      padj_cols, z_cols,
      SIMPLIFY = TRUE
    )
    if (is.null(dim(pos_mat))) pos_mat <- matrix(pos_mat, ncol = 1)
    df[[sprintf("pos_n_padj_meta_%s", a_tag)]] <- rowSums(pos_mat, na.rm = TRUE)

    neg_mat <- mapply(
      function(p, zc) {
        pv <- df[[p]]
        zv <- df[[zc]]
        as.integer(is.finite(pv) & pv < a & is.finite(zv) & zv < 0)
      },
      padj_cols, z_cols,
      SIMPLIFY = TRUE
    )
    if (is.null(dim(neg_mat))) neg_mat <- matrix(neg_mat, ncol = 1)
    df[[sprintf("neg_n_padj_meta_%s", a_tag)]] <- rowSums(neg_mat, na.rm = TRUE)
  }
  df
}

.write_summary_outputs_meta_csv <- function(df, out_dir, group_name,
                                            stat_tag = "GSEA_limma_t_cont_meta_fdr") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("please install data.table")
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("please install openxlsx")
  if (is.null(df) || !nrow(df)) {
    if (exists("log_msg", mode = "function")) try(log_msg("  [skip output] {group_name} | {stat_tag} No results available"), silent = TRUE)
    return(invisible(NULL))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(out_dir, paste0("Summary_", safe_fs_name(group_name), "_", stat_tag))


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

  ## ---- XLSX output----
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
    try(log_msg("  [output completed] {group_name} | {stat_tag} -> {dirname(base)} (.csv + .xlsx)"), silent = TRUE)
  }
  invisible(NULL)
}


summarize_meta_fdr_across_subunits <- function(
  meta_root = if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  versions = c("BatchAdj"),
  genesets_by_group = genesets_by_group,
  stat_tag = "GSEA_limma_t_cont_meta_fdr"
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("請先安裝 data.table")
  for (st in strata) {
    for (ver in versions) {
      base_dir <- file.path(meta_root, st, ver)
      if (!dir.exists(base_dir)) next

      subs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
      subs <- subs[nzchar(subs)]
      if (!length(subs)) next

      for (grp in names(genesets_by_group)) {
        if (exists("log_msg", mode = "function")) try(log_msg("== [meta-summary] stratum=%s | version=%s | group=%s ==", st, ver, grp), silent = TRUE)
        lst <- setNames(vector("list", length(subs)), subs)
        for (su in subs) lst[[su]] <- .read_meta_fdr_table(meta_root, st, ver, su, grp, stat_tag)
        wide <- .merge_subunit_tables_meta(lst)
        wide <- .add_sig_counts_meta(wide, alphas = c(0.05, 0.25))
        out_dir <- file.path(meta_root, "summary", st, ver, safe_fs_name(grp), stat_tag)
        .write_summary_outputs_meta_csv(wide, out_dir, grp, stat_tag)
      }
    }
  }
  invisible(TRUE)
}

## ----execute----
posthoc_summary_meta_fdr <- function() {
  summarize_meta_fdr_across_subunits(
    meta_root = if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
    strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
    versions = c("BatchAdj"),
    genesets_by_group = genesets_by_group,
    stat_tag = "GSEA_limma_t_cont_meta_fdr"
  )
  invisible(TRUE)
}

posthoc_summary_meta_fdr()


## =========================================================
##  Summarize meta-FDR for limma interaction across subunits
## =========================================================

posthoc_summary_meta_fdr_interaction <- function() {
  summarize_meta_fdr_across_subunits(
    meta_root = if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
    strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
    versions = c("BatchAdj"),
    genesets_by_group = genesets_by_group,
    stat_tag = "GSEA_limma_interaction_meta_fdr"
  )
  invisible(TRUE)
}


posthoc_summary_meta_fdr_interaction()


## =========================================================
## Heatmaps for meta-FDR summary CSVs (Z + padj_meta)
## =========================================================


if (!exists(".read_csv_safe", mode = "function")) {
  .read_csv_safe <- function(path) {
    p <- path

    if (!file.exists(p) && .Platform$OS.type == "windows") {
      p2 <- gsub("/", "\\\\", p, fixed = TRUE)
      if (file.exists(p2)) p <- p2
    }
    if (!file.exists(p)) stop("The file does not exist:", p)

    fi <- suppressWarnings(file.info(p))
    if (isTRUE(fi$isdir)) stop("The target is not the file:", p)
    if (!is.na(fi$size) && fi$size == 0) stop("empty file（size = 0）：", p)


    out <- try(suppressMessages(readr::read_csv(p, show_col_types = FALSE, progress = FALSE)), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }


    out <- try(suppressMessages(readr::read_csv(
      p,
      locale = readr::locale(encoding = "UTF-8-BOM"),
      show_col_types = FALSE, progress = FALSE
    )), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }


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


    out <- try(utils::read.csv(p, check.names = FALSE), silent = TRUE)
    if (inherits(out, "try-error")) stop(out)
    tibble::as_tibble(out)
  }
}

if (!exists(".safe_fs", mode = "function")) {
  .safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

if (!exists(".pred_order_all", mode = "any")) {
  .pred_order_all <- c(
    "NES_CSN_SCORE", "NES_GPS1",
    "NES_COPS2", "NES_COPS3", "NES_COPS4", "NES_COPS5", "NES_COPS6",
    "NES_COPS7A", "NES_COPS7B", "NES_COPS8", "NES_COPS9",
    "NES_RESIDUAL_GPS1",
    "NES_RESIDUAL_COPS2", "NES_RESIDUAL_COPS3", "NES_RESIDUAL_COPS4", "NES_RESIDUAL_COPS5", "NES_RESIDUAL_COPS6",
    "NES_RESIDUAL_COPS7A", "NES_RESIDUAL_COPS7B", "NES_RESIDUAL_COPS8", "NES_RESIDUAL_COPS9"
  )
}


.pred_order_meta <- gsub("^NES_", "Z_", .pred_order_all)


.meta_read_summary_long <- function(csv_file) {
  df <- suppressMessages(.read_csv_safe(csv_file))

  path_candidates <- c("pathway", "Pathway", "term", "Term", "gs_name", "NAME", "set", "Set")
  path_col <- intersect(path_candidates, names(df))[1]
  if (is.na(path_col)) stop("cannot find pathway in: ", csv_file)

  z_cols <- grep("^Z_", names(df), value = TRUE)
  padj_cols <- grep("^padj_meta_", names(df), value = TRUE)
  if (!length(z_cols)) stop("cannot find Z_* in: ", csv_file)

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


.meta_compute_y_order <- function(csv_file) {
  df_long <- .meta_read_summary_long(csv_file)
  present <- intersect(.pred_order_meta, unique(df_long$predictor))
  if (!length(present)) stop("無可用 predictors 於: ", csv_file)
  df_long <- dplyr::filter(df_long, .data$predictor %in% present)


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


.meta_make_heatmap_plot <- function(csv_file, y_order = NULL,
                                    palette = c(low = "#053061", mid = "#FFFFFF", high = "#67001F")) {
  df_long <- .meta_read_summary_long(csv_file)

  present <- intersect(.pred_order_meta, unique(df_long$predictor))
  if (!length(present)) stop("The data does not specify any predictors:", csv_file)
  df_long <- dplyr::filter(df_long, .data$predictor %in% present)

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


  if (is.null(y_order)) {
    y_order <- .meta_compute_y_order(csv_file)
  }

  df_long <- df_long %>% dplyr::filter(.data$pathway %in% y_order)
  y_levels <- y_order[y_order %in% df_long$pathway]
  if (!length(y_levels)) stop("y-order has no overlap with the data:", csv_file)
  df_long <- dplyr::mutate(df_long, pathway = factor(.data$pathway, levels = rev(y_levels)))


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

  z_fin <- df_long$Z[is.finite(df_long$Z)]
  L <- if (length(z_fin)) max(abs(z_fin), na.rm = TRUE) else 0
  if (L == 0) L <- 1

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


.meta_parse_summary_path <- function(csv_file) {
  parts <- strsplit(normalizePath(csv_file, winslash = "/"), "/")[[1]]
  parts_lc <- tolower(parts)
  ver_idx <- which(parts_lc %in% c("BatchAdj"))
  version <- if (length(ver_idx)) parts[ver_idx[1]] else NA_character_
  st_idx <- which(parts_lc %in% c("all", "tp53_mutant", "tp53_wild_type"))
  stratum <- if (length(st_idx)) parts[st_idx[length(st_idx)]] else NA_character_
  list(version = version, stratum = stratum)
}


.meta_save_near_v3_shortname <- function(p, width, height, csv_file, meta, suffix = "", dpi = 600) {
  out_dir <- dirname(csv_file)
  bn <- basename(csv_file)

  grp_tok <- sub("^Summary_(.+?)_GSEA.*$", "\\1", bn, perl = TRUE)
  grp_tok <- .safe_fs(grp_tok)
  method_tok <- if (grepl("limma_interaction", bn, ignore.case = TRUE)) {
    "interaction"
  } else if (grepl("limma_t_cont", bn, ignore.case = TRUE)) {
    "tcont"
  } else {
    "method"
  }
  pre <- paste0("heatmap_", .safe_fs(meta$stratum), "_", .safe_fs(meta$version), "_")
  base_name <- paste0(pre, grp_tok, "_", method_tok, suffix)


  target <- file.path(out_dir, paste0(base_name, ".tiff"))
  if (.Platform$OS.type == "windows") {
    dir_len <- nchar(normalizePath(out_dir, winslash = "/", mustWork = FALSE))
    max_full <- 240L
    max_bn <- max_full - dir_len - 1L - nchar(".tiff")
    if (max_bn < nchar(base_name)) {
      base_name <- substr(base_name, 1L, max(1L, max_bn))
      target <- file.path(out_dir, paste0(base_name, ".tiff"))
    }
  }

  ggplot2::ggsave(target, p,
    width = width, height = height, units = "in",
    dpi = dpi, bg = "white", compression = "lzw"
  )
}


.meta_find_target_csvs <- function(root = "csn_gsea_pan_summary_TP53/meta_fdr/summary") {
  if (!dir.exists(root)) {
    return(character(0))
  }
  pat_all_t <- "^Summary_.+_GSEA_limma_t_cont_meta_fdr_ALL\\.csv$"
  pat_all_i <- "^Summary_.+_GSEA_limma_interaction_meta_fdr_ALL\\.csv$"
  pat_mut_t <- "^Summary_.+_GSEA_limma_t_cont_meta_fdr_ALL\\.csv$"
  pat_wt_t <- "^Summary_.+_GSEA_limma_t_cont_meta_fdr_ALL\\.csv$"


  files <- list.files(root, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  # 只接受 BatchAdj
  keep <- grepl("(/|\\\\)(BatchAdj)(/|\\\\)", files, ignore.case = TRUE)

  files <- files[keep]
  if (!length(files)) {
    return(files)
  }

  bas <- basename(files)
  is_all <- grepl("(/|\\\\)ALL(/|\\\\)", files, ignore.case = TRUE)
  is_mut <- grepl("(/|\\\\)TP53_mutant(/|\\\\)", files, ignore.case = TRUE)
  is_wt <- grepl("(/|\\\\)TP53_wild_type(/|\\\\)", files, ignore.case = TRUE)

  tgt <- c(
    files[is_all & grepl(pat_all_t, bas)],
    files[is_all & grepl(pat_all_i, bas)],
    files[is_mut & grepl(pat_mut_t, bas)],
    files[is_wt & grepl(pat_wt_t, bas)]
  )


  if (length(tgt)) {
    finfo <- suppressWarnings(file.info(tgt))
    tgt <- tgt[is.finite(finfo$size) & finfo$size > 0 & !finfo$isdir]
  }
  unique(tgt)
}


if (!exists(".find_target_csvs_for_dataset", mode = "function")) {
  .find_target_csvs_for_dataset <- function(ds_dir) {
    # 搜尋名稱匹配的 Summary 檔
    files <- list.files(
      ds_dir,
      pattern = "^Summary_.+_GSEA_(limma_t(_cont)?|limma_interaction)_(ALL|TP53_mutant|TP53_wild_type)\\.csv$",
      recursive = TRUE, full.names = TRUE
    )
    if (!length(files)) {
      return(files)
    }


    sep_pat <- "[/\\\\]"
    keep <-
      grepl(paste0(sep_pat, "batchadj", sep_pat), files, ignore.case = TRUE)
    files <- files[keep]


    if (length(files)) {
      finfo <- suppressWarnings(file.info(files))
      files <- files[is.finite(finfo$size) & finfo$size > 0 & !finfo$isdir]
    }
    files
  }
}


if (!exists(".parse_meta_from_path", mode = "function")) {
  .parse_meta_from_path <- function(csv_file) {
    parts <- strsplit(normalizePath(csv_file, winslash = "/"), "/")[[1]]
    parts_lc <- tolower(parts)
    bn <- basename(csv_file)


    vidx <- which(parts_lc %in% c("batchadj"))
    variant <- if (length(vidx)) parts[vidx[1]] else NA_character_


    stratum <- NA_character_
    if (length(vidx) && vidx[1] - 1 >= 1) {
      cand <- parts[vidx[1] - 1]
      if (tolower(cand) %in% c("all", "tp53_mutant", "tp53_wild_type")) {
        stratum <- cand
      }
    }

    if (is.na(stratum) || !nzchar(stratum)) {
      hit <- which(parts_lc %in% c("all", "tp53_mutant", "tp53_wild_type"))
      if (length(hit)) stratum <- parts[hit[1]]
    }


    dataset <- NA_character_
    hit2 <- which(parts_lc == "csn_gsea_results_tp53")
    if (length(hit2) && hit2[1] > 1) dataset <- parts[hit2[1] - 1]

    list(dataset = dataset, variant = variant, stratum = stratum)
  }
}


# Main program: Plots the found meta_fdr Summary file and generates an ordered version (depending on ALL/t_cont).
run_meta_fdr_heatmaps <- function(root = if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr/summary" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr/summary")) {
  csvs <- .meta_find_target_csvs(root)

  .cols_pan <- get0("PLOT_PAN_COLLECTIONS", ifnotfound = NULL)
  if (!is.null(.cols_pan) && length(.cols_pan)) {
    .al <- unique(vapply(.cols_pan, safe_fs_name, FUN.VALUE = character(1)))
    .bas <- basename(csvs)
    .pref <- paste0("Summary_", .al, "_")
    .keep <- Reduce(`|`, lapply(.pref, function(p) startsWith(.bas, p)))
    csvs <- csvs[.keep]
  }
  if (!length(csvs)) {
    message("[meta-heatmap] The target file Summary_*_meta_fdr_ALL.csv could not be found, so it will be skipped.")
    return(invisible(TRUE))
  }
  for (csv in csvs) {
    meta <- .meta_parse_summary_path(csv)
    bn_lc <- tolower(basename(csv))
    is_interaction <- grepl("gsea_limma_interaction_meta_fdr", bn_lc, fixed = TRUE)


    pal_cell <- NULL
    if (is_interaction) {
      pal0 <- .get_interaction_palette() # c(neg=..., mid=..., pos=...)
      pal_cell <- c(low = pal0[["neg"]], mid = pal0[["mid"]], high = pal0[["pos"]])
    }


    h <- try(
      {
        if (is_interaction) {
          .meta_make_heatmap_plot(csv, y_order = NULL, palette = pal_cell)
        } else {
          .meta_make_heatmap_plot(csv)
        }
      },
      silent = TRUE
    )

    if (inherits(h, "try-error")) {
      message("[meta-heatmap] plotting failed：", csv, " | ", as.character(h))
      next
    }
    .meta_save_near_v3_shortname(h$plot, h$width, h$height, csv, meta, suffix = "")


    is_interaction_all <- grepl("gsea_limma_interaction_meta_fdr_all\\.csv$", bn_lc)
    if (is_interaction_all && identical(tolower(meta$stratum), "all")) {
      dir_now <- normalizePath(dirname(csv), winslash = "/", mustWork = FALSE)
      parent_dir <- dirname(dir_now)
      cont_dir <- file.path(parent_dir, "GSEA_limma_t_cont_meta_fdr")
      cont_csvs <- list.files(cont_dir, pattern = "^Summary_.+_GSEA_limma_t_cont_meta_fdr_ALL\\.csv$", full.names = TRUE)
      cont_csv <- if (length(cont_csvs)) cont_csvs[1L] else ""
      if (file.exists(cont_csv)) {
        yo_full <- try(.meta_compute_y_order(cont_csv), silent = TRUE)
        if (!inherits(yo_full, "try-error")) {
          .coll_tok2 <- sub("^Summary_(.+?)_GSEA.*$", "\\1", basename(cont_csv))
          if (!identical(toupper(.coll_tok2), "H")) {
            .top_n <- get0("PAN_HEATMAP_TOP_N", ifnotfound = 25)
            .bot_n <- get0("PAN_HEATMAP_BOTTOM_N", ifnotfound = 25)
            d50 <- .meta_read_summary_long(cont_csv) %>%
              dplyr::filter(.data$predictor == "Z_CSN_SCORE")
            keep_up <- d50 %>%
              dplyr::arrange(dplyr::desc(.data$Z)) %>%
              dplyr::slice_head(n = .top_n) %>%
              dplyr::pull(.data$pathway)
            keep_dn <- d50 %>%
              dplyr::arrange(.data$Z) %>%
              dplyr::slice_head(n = .bot_n) %>%
              dplyr::pull(.data$pathway)
            keep50 <- unique(c(keep_up, keep_dn))
          } else {
            keep50 <- yo_full
          }

          y_levels <- yo_full[yo_full %in% keep50]
          h2 <- try(.meta_make_heatmap_plot(csv, y_order = y_levels, palette = pal_cell), silent = TRUE)
          if (!inherits(h2, "try-error")) {
            .meta_save_near_v3_shortname(h2$plot, h2$width, h2$height, csv, meta, suffix = "_ordered")
          } else {
            message("[meta-heatmap][ordered] plotting failed:", csv, " | ", as.character(h2))
          }
        } else {
          message("[meta-heatmap][ordered] failed to get y-order（ALL/t_cont）：", cont_csv, " | ", as.character(yo_full))
        }
      } else {
        message("[meta-heatmap][ordered] cannot find ALL/t_cont file：", cont_csv)
      }
    }


    is_tcont_all <- grepl("gsea_limma_t_cont_meta_fdr_all\\.csv$", bn_lc)
    if (is_tcont_all && tolower(meta$stratum) %in% c("tp53_mutant", "tp53_wild_type")) {
      parts <- strsplit(normalizePath(csv, winslash = "/"), "/")[[1]]
      parts_lc <- tolower(parts)
      st_idx <- which(parts_lc %in% c("all", "tp53_mutant", "tp53_wild_type"))
      if (length(st_idx)) {
        parts2 <- parts
        parts2[st_idx[length(st_idx)]] <- "ALL"
        all_csv <- paste(parts2, collapse = "/")
        if (file.exists(all_csv)) {
          yo_full <- try(.meta_compute_y_order(all_csv), silent = TRUE)
          if (!inherits(yo_full, "try-error")) {
            #
            .coll_tok2 <- sub("^Summary_(.+?)_GSEA.*$", "\\1", basename(all_csv))
            if (!identical(toupper(.coll_tok2), "H")) {
              .top_n <- get0("PAN_HEATMAP_TOP_N", ifnotfound = 25)
              .bot_n <- get0("PAN_HEATMAP_BOTTOM_N", ifnotfound = 25)
              d50 <- .meta_read_summary_long(all_csv) %>%
                dplyr::filter(.data$predictor == "Z_CSN_SCORE")
              keep_up <- d50 %>%
                dplyr::arrange(dplyr::desc(.data$Z)) %>%
                dplyr::slice_head(n = .top_n) %>%
                dplyr::pull(.data$pathway)
              keep_dn <- d50 %>%
                dplyr::arrange(.data$Z) %>%
                dplyr::slice_head(n = .bot_n) %>%
                dplyr::pull(.data$pathway)
              keep50 <- unique(c(keep_up, keep_dn))
            } else {
              keep50 <- yo_full
            }

            y_levels <- yo_full[yo_full %in% keep50]
            h3 <- try(.meta_make_heatmap_plot(csv, y_order = y_levels), silent = TRUE)
            if (!inherits(h3, "try-error")) {
              .meta_save_near_v3_shortname(h3$plot, h3$width, h3$height, csv, meta, suffix = "_ordered")
            } else {
              message("[meta-heatmap][ordered] plotting failed", csv, " | ", as.character(h3))
            }
          } else {
            message("[meta-heatmap][ordered] Retrieving the y-order of ALL/t_cont failed:", all_csv, " | ", as.character(yo_full))
          }
        } else {
          message("[meta-heatmap][ordered] No corresponding ALL/t_cont file found:", all_csv)
        }
      }
    }

    try(closeAllConnections(), silent = TRUE)
    invisible(gc(FALSE))
    message("[meta-heatmap] completed：", csv)
  }
  invisible(TRUE)
}


## ---- execute----
run_meta_fdr_heatmaps(
  root = if (is.null(GSEA_PREFIX)) {
    "csn_gsea_pan_summary_TP53/meta_fdr/summary"
  } else {
    file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr/summary")
  }
)
