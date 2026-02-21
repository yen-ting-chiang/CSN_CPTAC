setwd("C:/Users/danny/Documents/R_project/CSN_CPTAC") ## YTC laptop
## setwd("C:/Users/cmuh/Documents/YenTing_document/CSN_CPTAC") ## lab computer
getwd()

## ----------------------------------------------------------


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


dir.create("run_info", recursive = TRUE, showWarnings = FALSE)



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



set.seed(1234)

csn_subunits <- c("GPS1", "COPS2", "COPS3", "COPS4", "COPS5", "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9")

min_frac_complete <- 0.75
hi_lo_quantile <- 0.25
minSize <- 5
maxSize <- 1000
fgsea_eps <- 1e-10

min_per_group <- 8

MAKE_PLOTS <- FALSE


RUN_HILO_SENSITIVITY <- FALSE


RUN_PASSES <- c("BatchAdj")


options(csn.run_passes = RUN_PASSES)
## ---- [END USER CONFIG] ----

## ---- [USER CONFIG | heatmap top/bottom limits] ----

DATASET_HEATMAP_TOP_N <- 25
DATASET_HEATMAP_BOTTOM_N <- 25
PAN_HEATMAP_TOP_N <- 40
PAN_HEATMAP_BOTTOM_N <- 25
## ---- [END USER CONFIG | heatmap top/bottom limits] ----


## ---- [USER CONFIG | heatmap collections] ----
PLOT_DATASET_COLLECTIONS <- NULL
PLOT_PAN_COLLECTIONS <- NULL
## ---- [END USER CONFIG | heatmap collections] ----

## ---- [USER CONFIG | heatmap pathway prefix filter] ----
options(HEATMAP_PATHWAY_PREFIX = get0("HEATMAP_PATHWAY_PREFIX", ifnotfound = NULL))
## ---- [END USER CONFIG | heatmap pathway prefix filter] ----
## ---- [USER CONFIG | heatmap dot sizes] ----

options(csn.dot_size_dataset = get0("csn.dot_size_dataset", ifnotfound = 1.6))
options(csn.dot_size_meta = get0("csn.dot_size_meta", ifnotfound = 1.6))
## ---- [END USER CONFIG | heatmap dot sizes] ----



PIPELINES_TO_RUN <- c("limma_t")
PIPELINES_TO_RUN <- get0("PIPELINES_TO_RUN", ifnotfound = c("limma_t"))
.RUN_LIMMA <- any(tolower(PIPELINES_TO_RUN) %in% c("limma_t", "limma", "both", "all"))


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


GENESET_GROUPS_TO_RUN <- c("H")




log_msg(
  "Prepare MSigDB gene sets; build according to GENESET_GROUPS_TO_RUN: %s",
  paste(GENESET_GROUPS_TO_RUN, collapse = ", ")
)

df0 <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens"), error = function(e) NULL)

if (!is.null(df0) && NROW(df0) > 0) {
  df0 <- as.data.frame(df0)


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


    ord <- order(counts$gs_cat, ifelse(is.na(counts$gs_subcat), "", counts$gs_subcat))
    counts <- counts[ord, , drop = FALSE]


    dir.create("run_info", showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(counts, file.path("run_info", "msigdb_available_collections.csv"))
    log_msg("The msigdb available collections (including geneset counts) have been output to run_info/msigdb_available_collections.csv")
  }
} else {
  log_msg("(Omitted) msigdbr() returns null or the call failed; no collections list is output.")
}

genesets_by_group <- list()

Sys.setenv(PTMSIGDB_GMT = "C:/Users/danny/Documents/R_project/CSN_CPTAC/ptm.sig.db.all.uniprot.human.v2.0.0.gmt")

## ---- [NEW | PTMsigDB helper for site-level PTM-SEA] ----
.std_site_id <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("[:\\s]", "_", x)
  x <- gsub("_([STY])(\\d+)[A-Z]*$", "_\\1\\2", x)
  x
}


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

  # helper
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


  cand <- id_from_name
  need <- is.na(cand) | !nzchar(cand)
  cand[need] <- id_from_gp[need]
  need <- is.na(cand) | !nzchar(cand)
  cand[need] <- id_from_ent[need]

  .std_site_id(cand)
}


read_ptmsigdb_gmt <- function(fp) {
  if (missing(fp) || !nzchar(fp) || !file.exists(fp)) {
    stop("[PTMsigDB] File does not exist:", fp %||% "<empty>")
  }

  std_id <- function(x) {
    x <- toupper(trimws(x))
    x <- sub(":.*$", "", x)
    x <- gsub("[^A-Z0-9_]", "", x)
    x
  }


  if (grepl("\\.xlsx$", fp, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("[PTMsigDB] requires the readxl package to read .xlsx files.")
    }
    df <- readxl::read_xlsx(fp)
    req <- c("signature", "site.annotation")
    if (!all(req %in% names(df))) {
      stop("[PTMsigDB].xlsx is missing necessary fields:", paste(setdiff(req, names(df)), collapse = ", "))
    }

    df$gene_site <- std_id(df$site.annotation)
    by_sig <- split(df$gene_site, df$signature)
    out <- lapply(by_sig, function(v) unique(v[nzchar(v)]))
    out[lengths(out) == 0] <- NULL
    return(out)
  }


  lines <- readr::read_lines(fp)
  res <- vector("list", length(lines))
  nm <- character(length(lines))

  for (i in seq_along(lines)) {
    fields <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    if (length(fields) < 2) next
    sig <- fields[1]
    nm[i] <- sig


    toks <- if (length(fields) >= 3) fields[-(1:2)] else character(0)
    if (length(toks) == 0 && grepl("\\|", fields[2], fixed = TRUE)) {
      t2 <- strsplit(fields[2], "\\|")[[1]]
      toks <- t2[!grepl("^n=\\d+$", t2, ignore.case = TRUE)]
    }


    toks <- std_id(toks)


    keep <- grepl("^[A-Z0-9._-]+_[STY][0-9]+$", toks)
    toks <- unique(toks[keep & nzchar(toks)])

    res[[i]] <- toks
  }

  names(res) <- nm
  res[lengths(res) == 0] <- NULL

  if (!length(res)) {
    stop(
      "The [PTMsigDB] parsing result is empty. This GMT is likely a UniProt version (e.g., O14974; T696).",
      " Please use a source with gene-symbol instead (.xlsx or gene version of GMT)."
    )
  }
  res
}



Sys.setenv(PTMSIGDB_GMT = "C:/Users/danny/Documents/R_project/CSN_CPTAC/data_PTMsigDB_all_sites_v2.0.0.xlsx")

PTMSIGDB_GMT_FP <- Sys.getenv("PTMSIGDB_GMT", unset = "")
genesets_by_group_ptm <- list()
if (nzchar(PTMSIGDB_GMT_FP) && file.exists(PTMSIGDB_GMT_FP)) {
  genesets_by_group_ptm[["PTMsigDB"]] <- read_ptmsigdb_gmt(PTMSIGDB_GMT_FP)
  log_msg("Loaded PTMsigDB site collections：%d sets", length(genesets_by_group_ptm[["PTMsigDB"]]))
} else {
  log_msg("(Reminder) PTMSIGDB_GMT has not been set; site-level PTM-SEA will not be able to execute.")
}
## ---- [END NEW] ----


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
    log_msg("Note: msigdbr cannot provide categories %s；skip %s", cat, pp$raw)
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
      log_msg("Notice: %s does not contain any subset of %s; skip this step.", cat, pp$sub)
      next
    }
    grp <- .make_group_label(cat, df, pp$sub)
  } else {
    grp <- cat
  }

  genesets_by_group[[grp]] <- lapply(split(df$gene_symbol, df$gs_name), unique)
  log_msg("Get %s:%d sets", grp, length(genesets_by_group[[grp]]))
}

if (!length(genesets_by_group)) {
  stop("No gene-set is available. Please check the GENESET_GROUPS_TO_RUN setting.")
}
log_msg("Final gene-set groups: %s", paste(names(genesets_by_group), collapse = ", "))



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
  if (!file.exists(fp)) stop(glue::glue("cannot find files：{fp}"))
  log_msg("read protein matrix：{basename(fp)}")
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
    if (length(inter) < 10) log_msg("Note: The intersection of case_list is too small ({length(inter)}), use all sample fields instead.")
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
  if (ncol(m) == 0) stop("0 sample columns read")
  if (anyDuplicated(rownames(m))) {
    log_msg("When duplicate genes are detected, the average of the duplicate rows is taken.")
    m <- rowsum(m, group = rownames(m), reorder = FALSE) / as.vector(table(rownames(m)))
  }
  log_msg("matrix: {nrow(m)} genes × {ncol(m)} samples")
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
    if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] Insufficient available subcells or samples → All NA")
    return(out_na)
  }
  pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) {
    if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] prcomp failed → All NA")
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
      log_msg("[CSN-audit-safe] %s | %s: Insufficient available sub-cells or samples, audit skipped (without aborting).", ds_id, stratum)
    }
    return(invisible(FALSE))
  }
  pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) {
    if (exists("log_msg", mode = "function")) {
      log_msg("[CSN-audit-safe] %s | %s: fallback prcomp still failed, skipped (without aborting).", ds_id, stratum)
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


._ensure_stats_names <- function(stats, gene_names, label = NULL) {
  if (is.null(stats)) {
    return(NULL)
  }
  v <- suppressWarnings(as.numeric(stats))
  if (length(v) != length(gene_names)) {
    stop(sprintf(
      "The lengths of `[ensure-names%s] stats(%d)` and `gene_names(%d)` do not match.",
      if (!is.null(label)) paste0("-", label) else "",
      length(v), length(gene_names)
    ))
  }
  names(v) <- as.character(gene_names)

  v
}


.gsea_from_ranks <- function(pathways, stats, minSize, maxSize, gsea_eps = 0, label = NULL) {
  if (is.null(pathways) || !length(pathways) || is.null(stats) || !length(stats)) {
    return(NULL)
  }
  orig_n <- length(pathways)
  pw_use <- ._intersect_and_filter_pathways(pathways, names(stats), minSize = minSize, maxSize = maxSize)
  message(sprintf(
    "[gsea:%s] original=%d → use=%d",
    ifelse(is.null(label), "NA", label),
    orig_n, length(pw_use)
  ))
  if (!length(pw_use)) {
    return(NULL)
  }

  nproc <- getOption(".fgsea_nproc", 1L)
  set.seed(1L)
  res <- tryCatch(
    {
      suppressWarnings(fgsea::fgseaMultilevel(
        pathways = pw_use, stats = stats,
        minSize = minSize, maxSize = maxSize,
        eps = gsea_eps, nproc = nproc
      ))
    },
    error = function(e) {
      message(sprintf("[gsea_from_ranks] fgseaMultilevel failed: %s; fallback to fgseaSimple", conditionMessage(e)))
      suppressWarnings(fgsea::fgseaSimple(
        pathways = pw_use, stats = stats,
        nperm = 10000, minSize = minSize, maxSize = maxSize
      ))
    }
  )
  if (!is.null(label)) attr(res, "label") <- label
  res
}


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

  if (is.null(nm)) {
    log_msg("[gsea-%s] stats has no names after filtering → skip", label)
    return(NULL)
  }

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

.align_limma_inputs <- function(M, design_df, rhs_terms, label = "limma") {
  stopifnot(!is.null(M), is.matrix(M), ncol(M) > 0, !is.null(colnames(M)))

  design_df <- as.data.frame(design_df, check.names = FALSE)
  stopifnot(!is.null(rownames(design_df)))
  samp_M <- as.character(colnames(M))
  samp_des <- as.character(rownames(design_df))


  na_mask <- rep(FALSE, nrow(design_df))
  rn <- rownames(design_df)


  if (!("predictor" %in% colnames(design_df))) {
    stop(sprintf("[%s] design_df lack predictor colums", label))
  }


  v <- design_df$predictor
  na_mask <- na_mask | is.na(v) | !is.finite(v)


  num_cols <- setdiff(colnames(design_df), c("predictor"))
  for (cn in num_cols) {
    v <- design_df[[cn]]
    if (is.numeric(v)) {
      na_mask <- na_mask | is.na(v) | !is.finite(v)
    } else {
      na_mask <- na_mask | is.na(v)
    }
  }
  design_df <- design_df[!na_mask, , drop = FALSE]


  common <- intersect(as.character(colnames(M)), as.character(rownames(design_df)))
  if (!length(common)) {
    return(NULL)
  }
  common <- sort(common)

  M2 <- M[, common, drop = FALSE]
  design_df2 <- design_df[common, , drop = FALSE]


  rhs <- unique(rhs_terms)
  des2 <- stats::model.matrix(stats::reformulate(rhs), data = design_df2)


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

  head_show <- utils::head(df[, c("pos", "M_sample", "interest_rn", "nuisance_rn", "interest_ok", "nuisance_ok")], 6)
  log_msg("[align:%s] head:\n%s", tag, utils::capture.output(print(head_show)) |> paste(collapse = "\n"))
  invisible(df)
}


impute_and_filter <- function(mat, min_frac = 0.75) {
  keep <- rowMeans(!is.na(mat)) >= min_frac
  m <- mat[keep, , drop = FALSE]
  if (any(is.na(m))) {
    set.seed(1234)
    m <- imputeLCMD::impute.MinProb(m, q = 0.01)
  }
  m
}





if (!exists(".ensure_mat_or_null", mode = "function")) {
  .ensure_mat_or_null <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    m <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (is.null(m)) {
      return(NULL)
    }
    if (is.null(dim(m))) { # 向量 → 單欄矩陣
      m <- matrix(m, ncol = 1)
    }
    if (ncol(m) == 0) {
      return(NULL)
    }
    m
  }
}



gsea_from_diffcorr <- function(pathways, z_stats, label = "diff-corr", out_prefix = NULL) {
  stopifnot(is.numeric(z_stats), !is.null(names(z_stats)))
  minSize <- opt("minSize", 15L)
  maxSize <- opt("maxSize", 500L)
  orig_n <- length(pathways)
  pw_use <- ._intersect_and_filter_pathways(pathways, names(z_stats), minSize = minSize, maxSize = maxSize)
  if (!is.null(out_prefix)) {
    message(sprintf("[gsea:%s] %s original=%d → use=%d", label, out_prefix, orig_n, length(pw_use)))
  } else {
    message(sprintf("[gsea:%s] original=%d → use=%d", label, orig_n, length(pw_use)))
  }
  if (!length(pw_use)) {
    return(invisible(NULL))
  }

  set.seed(1L)
  res <- tryCatch(
    {
      suppressWarnings(
        fgsea::fgseaMultilevel(
          pathways = pw_use, stats = z_stats,
          minSize = minSize, maxSize = maxSize,
          eps = 1e-10, nproc = getOption(".fgsea_nproc", 1L)
        )
      )
    },
    error = function(e) {
      message(sprintf("[gsea_from_diffcorr] Multilevel failed: %s → Use fgseaSimple instead.", conditionMessage(e)))
      suppressWarnings(
        fgsea::fgseaSimple(
          pathways = pw_use, stats = z_stats,
          nperm = 10000, minSize = minSize, maxSize = maxSize
        )
      )
    }
  )
  invisible(as.data.frame(res))
}

run_fgsea_save <- function(stats, genesets, out_prefix, top_plot_n = 0, plot_title = "") {
  set.seed(1234)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

  minSize <- opt("minSize", 15L)
  maxSize <- opt("maxSize", 500L)

  orig_n <- length(genesets)
  pw_use <- ._intersect_and_filter_pathways(genesets, names(stats), minSize = minSize, maxSize = maxSize)
  message(sprintf("[gsea] %s original=%d → use=%d", out_prefix, orig_n, length(pw_use)))
  if (!length(pw_use)) {
    data.table::fwrite(data.frame(), paste0(out_prefix, ".csv"))
    return(invisible(NULL))
  }

  bp <- tryCatch(
    {
      if (requireNamespace("BiocParallel", quietly = TRUE)) BiocParallel::SerialParam() else NULL
    },
    error = function(e) NULL
  )

  res <- tryCatch(
    {
      suppressWarnings(
        fgsea::fgseaMultilevel(
          pathways = pw_use,
          stats    = stats,
          minSize  = minSize,
          maxSize  = maxSize,
          eps      = 1e-10,
          BPPARAM  = bp
        )
      )
    },
    error = function(e) {
      log_msg("    [fgseaMultilevel] failed: %s → use fgsea()", e$message)
      suppressWarnings(
        fgsea::fgsea(
          pathways = pw_use,
          stats    = stats,
          minSize  = minSize,
          maxSize  = maxSize,
          nperm    = 1000,
          BPPARAM  = bp
        )
      )
    }
  )
  res <- as.data.frame(res)
  data.table::fwrite(res, paste0(out_prefix, ".csv"))
  if (top_plot_n > 0) {}
  invisible(res)
}


read_gsea_table <- function(out_root, subunit, group_name, stat_tag) {
  grp <- safe_fs_name(group_name)
  su <- subunit
  f <- paste0(stat_tag, ".csv")


  candidates <- c(
    file.path(out_root, su, f),
    file.path(out_root, grp, su, "H", f),
    file.path(out_root, grp, su, f),
    file.path(out_root, su, grp, "H", f),
    file.path(out_root, su, grp, f)
  )

  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    if (exists("log_msg", mode = "function")) try(log_msg("file not found: %s", paste(candidates, collapse = " | ")), silent = TRUE)
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
    log_msg("[Output omitted] {group_name} | {stat_tag} No results available")
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
  log_msg("  [Output complete] {group_name} | {stat_tag} -> {dirname(base)}")
}

summarize_all_groups <- function(out_root, csn_subunits, genesets_by_group,
                                 stat_tags = c("GSEA_limma_t")) {
  sum_root <- file.path(out_root, "summary")
  dir.create(sum_root, recursive = TRUE, showWarnings = FALSE)
  for (grp_name in names(genesets_by_group)) {
    grp_safe <- safe_fs_name(grp_name)
    for (stat_tag in stat_tags) {
      log_msg("==Summary:group={grp_name} | stat={stat_tag} ==")
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
safe_read_gsea <- function(fp) {
  if (!file.exists(fp)) {
    return(NULL)
  }
  dt <- tryCatch(data.table::fread(fp, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) {
    return(NULL)
  }
  need <- c("pathway", "NES", "pval")
  if (!all(need %in% names(dt))) {
    return(NULL)
  }
  as.data.frame(dt[, need])
}


meta_fdr_stouffer <- function(dataset_dirs,
                              strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
                              stat_tags = c("GSEA_limma_t_cont"),
                              groups = names(genesets_by_group),
                              out_root = if (is.null(COMBO_PREFIX)) "phospho_csn_gsea_pan_summary_TP53/meta_fdr" else file.path(COMBO_PREFIX, "phospho_csn_gsea_pan_summary_TP53/meta_fdr")) {
  stopifnot(length(dataset_dirs) > 0)

  if (!is.null(names(dataset_dirs))) {
    keep_idx <- names(dataset_dirs) != "lusc_cptac_2021"
  } else {
    keep_idx <- basename(dataset_dirs) != "lusc_cptac_2021"
  }
  dataset_dirs <- dataset_dirs[keep_idx]
  if (!length(dataset_dirs)) {
    message("[meta] All datasets were excluded or no datasets were available (lusc_cptac_2021 was excluded); omitted.")
    return(invisible(TRUE))
  }
  if (!requireNamespace("data.table", quietly = TRUE)) stop("install data.table first.")
  if (!requireNamespace("stats", quietly = TRUE)) stop("Missing stats suite")
  data.table::setDTthreads(1L)

  # ---- helpers --------------------------------------------------------------
  normalize_cols <- function(DT) {
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

    D[, z := sign(NES) * stats::qnorm(p = pmax(pmin(pval, 1 - 1e-10), 1e-10), lower.tail = FALSE)]

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
          subunits <- unique(unlist(lapply(dataset_dirs, function(dsdir) {
            base <- file.path(
              COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection,
              safe_fs_name(grp), pass_label, basename(dsdir), stratum
            )
            if (!dir.exists(base)) {
              return(character(0))
            }

            subs <- list.dirs(base, full.names = FALSE, recursive = FALSE)
            subs[subs != ""]
          })))

          if (length(subunits) == 0) {
            message(sprintf(
              "[meta] %s/%s/%s/%sNo subunit available, skip this step.",
              stratum, pass_label, stat_tag, grp
            ))
            next
          }

          # 2) 對每個 subunit，收集各 dataset 的 CSV、做 Stouffer 合併並落檔
          for (su in subunits) {
            parts <- list()
            for (dsdir in dataset_dirs) {
              f <- file.path(
                COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection,
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

## ==== PATCH: helpers ====

## === NEW: CSN complex score (PC1) & residualization helpers ==================


build_csn_score <- function(mat0,
                            subunits = csn_subunits,
                            combine_7AB = TRUE,
                            min_members = 5L) {
  present <- intersect(subunits, rownames(mat0))

  s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
  if (!length(present)) {
    return(s)
  }


  get_z <- function(v) {
    nm <- names(v)
    v <- as.numeric(v)
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    out <- (v - mu) / sdv
    names(out) <- nm
    out
  }

  X <- do.call(rbind, lapply(present, function(g) get_z(mat0[g, ])))
  rownames(X) <- present
  colnames(X) <- colnames(mat0)

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

    mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)
    if (suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs")) < 0) sc <- -sc
    s[keep_sam] <- sc
  }

  s
}


residualize_vector <- function(y, csn_score, batch = NULL, covars = NULL, min_n = 8L) {
  if (is.null(names(y))) stop("[residualize_vector]y must be a numeric vector with sample names.")
  sam <- names(y)
  common <- sam
  if (!is.null(csn_score)) common <- intersect(common, names(csn_score))
  if (!is.null(batch)) common <- intersect(common, names(batch))
  if (!is.null(covars)) {
    rn <- rownames(as.data.frame(covars, check.names = FALSE))
    if (is.null(rn)) stop("[residualize_vector]covars must have rownames=sample")
    common <- intersect(common, rn)
  }
  if (length(common) < min_n) {
    return(setNames(rep(NA_real_, length(y)), names(y)))
  }


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
    C <- coerce_covariates_safely(C)
    for (cn in colnames(C)) DF[[cn]] <- C[[cn]]
  }
  DF_y <- suppressWarnings(as.numeric(y[common]))


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


  des <- stats::model.matrix(~ 1 + ., data = DF, na.action = stats::na.fail)
  stopifnot(nrow(des) == length(y_ok))


  fit <- lm.fit(x = des, y = y_ok)
  res <- rep(NA_real_, length(DF_y))
  res[ok] <- fit$residuals


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
      mean(!is.na(v)) * 100
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


  X <- cbind("(Intercept)" = rep(1, nrow(X)), X)


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

gate_tech_pcs <- function(cv, v, thr = 0.3, alpha = 0.01) {
  if (is.null(cv) || ncol(cv) == 0) {
    return(cv)
  }
  pc_cols <- grep("^PC[0-9]+$", colnames(cv), value = TRUE)
  if (!length(pc_cols)) {
    return(cv)
  }
  keep <- sapply(pc_cols, function(cn) {
    r <- suppressWarnings(stats::cor(cv[, cn], v, use = "pairwise.complete.obs"))
    if (!is.finite(r)) {
      return(TRUE)
    }
    n <- sum(stats::complete.cases(cv[, cn], v))
    if (n < 15) {
      return(TRUE)
    }
    t <- r * sqrt((n - 2) / (1 - r^2))
    p <- 2 * stats::pt(-abs(t), df = n - 2)
    !(abs(r) >= thr && p < alpha)
  })
  cbind(
    cv[, setdiff(colnames(cv), pc_cols), drop = FALSE],
    cv[, pc_cols[keep], drop = FALSE]
  )
}

.take_first <- function(x) sub("\\|.*$", "", as.character(x))
.sanitize_levels <- function(f) {
  if (is.null(f)) {
    return(NULL)
  }
  f <- droplevels(factor(f))
  lv <- levels(f)
  lv2 <- make.names(lv)
  lv2 <- paste0("b_", lv2)
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


BATCH_PIPE_POLICY <- "NA"
BATCH_MIN_PER_LEVEL <- 2


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


  if (!is.null(min_per_level) && min_per_level > 1) {
    tab <- table(fac, useNA = "no")
    small <- names(tab)[tab < min_per_level]
    if (length(small)) {
      log_msg(
        "  [batch] Merge sparse levels to 'b_small': %s",
        paste(sprintf("%s(n=%d)", small, as.integer(tab[small])), collapse = ", ")
      )
      fac_chr <- as.character(fac)
      fac_chr[fac_chr %in% small] <- "b_small"
      fac <- droplevels(factor(fac_chr))
    }
  }
  fac
}


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

    if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
      return(list(name = cn, fac = fac))
    }
  }
  NULL
}


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


  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
  rownames(meta) <- sample_ids

  det <- detect_batch_column(meta,
    pipe_policy   = pipe_policy,
    min_per_level = min_per_level
  )
  if (is.null(det)) {
    tmt_fp <- file.path(ds_dir, "TMT_protein.csv")
    if (file.exists(tmt_fp)) {
      tmt <- suppressMessages(readr::read_csv(tmt_fp, show_col_types = FALSE)) |> as.data.frame()


      cn <- names(tmt)
      cn_trim <- trimws(cn)
      names(tmt) <- cn_trim
      run_hits <- grep("^Run\\s*Metadata\\s*ID$", cn_trim, ignore.case = TRUE, value = TRUE)
      tmt_cols <- grep("^tmt_", cn_trim, ignore.case = TRUE, value = TRUE)

      if (length(run_hits) >= 1 && length(tmt_cols) >= 1) {
        run_col <- run_hits[1]


        plex_by_sample <- list()
        nR <- nrow(tmt)
        for (i in seq_len(nR)) {
          run_id <- as.character(tmt[[run_col]][i])
          if (!nzchar(run_id) || is.na(run_id)) next
          plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id)
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

          if (nlevels(fac2) >= 2 && sum(!is.na(fac2)) >= 3) {
            names(fac2) <- sample_ids
            return(list(name = "TMT_protein.csv:RunMetadataID", fac = fac2))
          }
        }
      }
    }

    return(NULL)
  }

  fac <- det$fac
  names(fac) <- sample_ids
  list(name = det$name, fac = fac)
}


get_batch_factor_phospho <- function(ds_dir, sample_ids,
                                     pipe_policy = c("strict", "lenient", "NA"),
                                     min_per_level = 3L) {
  pipe_policy <- match.arg(pipe_policy)

  meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  if (file.exists(meta_fp)) {
    meta <- suppressMessages(
      readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")
    ) |> as.data.frame()
    id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"), names(meta))
    if (length(id_cols)) {
      meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
      meta$SAMPLE_ID <- as.character(meta$SAMPLE_ID)

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
    plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id)
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


screen_batch_need <- function(ds_dir, min_frac_complete = 0.75) {
  log_msg("== Batch check：%s ==", basename(ds_dir))
  mat0 <- load_phospho_matrix_from_dataset_dir(ds_dir)
  mx <- suppressWarnings(max(mat0, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)
  m <- impute_and_filter(mat0, min_frac = min_frac_complete)

  bi <- get_batch_factor_phospho(ds_dir, colnames(m))
  if (is.null(bi)) {
    log_msg("  [batch] No clear batch field found → Do not correct yet")
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
    "  [batch] Fields: %s; levels=%d; n per level: %s",
    bi$name, nlevels(batch),
    paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")
  )


  X <- .clean_for_pca(m, min_samples = 10L, min_genes = 5L)
  if (is.null(X)) {
    r2 <- pv <- numeric(0)
    log_msg("  PCA (by batch) omitted: Insufficient available genes/samples or still contains non-finite values.")
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
      "  PCA (by batch) R²：%s ; p：%s", paste(round(r2, 3), collapse = ", "),
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
    "  Genetic level F test: FDR < 0.05 = %.1f%%; FDR < 0.25 = %.1f%%",
    100 * prop05, 100 * prop25
  )

  recommend <- (any(r2 >= 0.10 & pv[seq_along(r2)] < 0.01) || prop05 >= 0.05)
  if (recommend) {
    log_msg("  **Recommended Correction**: R² or gene percentage reaches the threshold (≥10% R² and p<0.01, or FDR<0.05 and gene percentage ≥5%).")
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

dataset_ids <- c(
  "brca_cptac_2020",
  "brca_tcga_pan_can_atlas_2018",
  "coad_cptac_2019",
  "coadread_tcga_pan_can_atlas_2018",
  "gbm_cptac_2021",
  "luad_cptac_2020",
  "lusc_cptac_2021",
  "ov_tcga_pan_can_atlas_2018",
  "paad_cptac_2021",
  "ucec_cptac_2020"
)


## ---- [NEW | COMBO PRESETS] ----

COMBO_PREFIX <- "phosphoproteomic_site_level_GSEA"

.combo_ds_1 <- c("brca_cptac_2020", "luad_cptac_2020", "lusc_cptac_2021", "ucec_cptac_2020")
.combo_ds_7 <- c(
  "brca_cptac_2020", "luad_cptac_2020", "lusc_cptac_2021", "ucec_cptac_2020",
  "coad_cptac_2019", "gbm_cptac_2021", "paad_cptac_2021"
)


dataset_ids <- .combo_ds_7
RUN_PASSES <<- c("BatchAdj")
PIPELINES_TO_RUN <<- c("limma_t")
.RUN_LIMMA <<- TRUE
## ---- [END NEW | COMBO PRESETS] ----




dataset_dirs <- setNames(file.path(datasets_root, dataset_ids), dataset_ids)
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")
message("datasets_root = ", datasets_root)


missing_dirs <- names(dataset_dirs)[!dir.exists(dataset_dirs)]
if (length(missing_dirs)) {
  log_msg("%d missing folders were detected; %s will be skipped.", length(missing_dirs), paste(missing_dirs, collapse = ", "))
}



dataset_dirs_run <- dataset_dirs[
  dir.exists(dataset_dirs) &
    file.exists(file.path(dataset_dirs, "data_phosphoprotein_quantification.txt"))
]

if (!length(dataset_dirs_run)) stop(" `dataset_dirs_run` is empty, please check if the folder and `data_phosphoprotein_quantification.txt` exist.")

log_msg("Available data set for this round（%d）：%s", length(dataset_dirs_run), paste(names(dataset_dirs_run), collapse = ", "))


build_tech_covars <- function(M) {
  med <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colMedians(M, na.rm = TRUE)
  } else {
    apply(M, 2, stats::median, na.rm = TRUE)
  }
  miss <- colMeans(!is.finite(M))


  Z <- t(scale(t(M), center = TRUE, scale = TRUE))
  Y <- t(Z)
  sds <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colSds(Y, na.rm = TRUE)
  } else {
    apply(Y, 2, stats::sd, na.rm = TRUE)
  }
  keep <- which(is.finite(sds) & sds > 0)

  pc1 <- rep(0, nrow(Y))
  pc2 <- rep(0, nrow(Y))
  if (length(keep) >= 2) {
    Yk <- Y[, keep, drop = FALSE]
    pc <- tryCatch(
      {
        if (requireNamespace("irlba", quietly = TRUE)) {
          irlba::prcomp_irlba(Yk, n = 2, center = TRUE, scale. = FALSE)
        } else {
          stats::prcomp(Yk, rank. = 2, center = TRUE, scale. = FALSE)
        }
      },
      error = function(e) NULL
    )
    if (!is.null(pc)) {
      pc1 <- pc$x[, 1]
      if (ncol(pc$x) >= 2) pc2 <- pc$x[, 2]
    }
  }
  out <- cbind(median = med, missrate = miss, PC1 = pc1, PC2 = pc2)
  rownames(out) <- colnames(M)
  as.matrix(out)
}


load_phospho_matrix_from_dataset_dir <- function(dir,
                                                 site_collapse = c("median", "max"),
                                                 protein_adjust = FALSE,
                                                 protein_fp = file.path(dir, "data_protein_quantification.txt")) {
  site_collapse <- match.arg(site_collapse)
  fp <- file.path(dir, "data_phosphoprotein_quantification.txt")
  if (!file.exists(fp)) stop("[phospho] The file does not exist.：", fp)

  suppressMessages({
    df <- readr::read_tsv(fp, show_col_types = FALSE, progress = FALSE)
  })


  meta_cols <- intersect(
    colnames(df),
    c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "GENE_SYMBOL", "PHOSPHOSITES")
  )
  sample_cols <- setdiff(colnames(df), meta_cols)
  if (!length(sample_cols)) stop("[phospho] can not find sample column")


  M_site <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(M_site) <- "numeric"
  rownames(M_site) <- df$GENE_SYMBOL


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


  if (isTRUE(protein_adjust) && file.exists(protein_fp)) {
    suppressMessages({
      prot_df <- readr::read_tsv(protein_fp, show_col_types = FALSE, progress = FALSE)
    })

    prot_meta <- intersect(
      colnames(prot_df),
      c(
        "Hugo_Symbol", "GENE_SYMBOL", "GeneSymbol", "Description",
        "GENE", "GENE_ID", "GENE.STABLE.ID", "Composite.Element.REF"
      )
    )
    prot_sample <- setdiff(colnames(prot_df), prot_meta)

    cand <- c("Hugo_Symbol", "GENE_SYMBOL", "GeneSymbol", "GENE", "Composite.Element.REF")
    cn <- colnames(prot_df)
    lc <- tolower(cn)
    hit <- match(tolower(cand), lc, nomatch = NA_integer_)
    hit <- hit[!is.na(hit)]
    gcol <- if (length(hit)) cn[hit[1]] else NA_character_
    if (is.na(gcol)) stop("[phospho] protein_adjust=TRUE but protein file lack gene column")

    prot_sample <- setdiff(prot_sample, gcol)
    M_prot <- as.matrix(prot_df[, prot_sample, drop = FALSE])
    storage.mode(M_prot) <- "numeric"
    rownames(M_prot) <- prot_df[[gcol]]
    rownames(M_prot) <- sub("\\|.*$", "", rownames(M_prot))


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


  colnames(M_gene) <- gsub("\\s+", "", colnames(M_gene))
  if (exists("log_msg")) log_msg("【phospho】matrix dimension：%d genes × %d samples", nrow(M_gene), ncol(M_gene))
  return(M_gene)
}


## ===== Imputation sensitivity =====





.align_to_samples <- function(x, sam, what = "covariate") {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.vector(x) || is.factor(x)) {
    if (is.null(names(x))) {
      if (length(x) != length(sam)) stop(sprintf("[%s] Length mismatch: %d vs %d", what, length(x), length(sam)))
      names(x) <- sam
    }
    if (!all(sam %in% names(x))) stop(sprintf("[%s] Missing sample: %s", what, paste(setdiff(sam, names(x)), collapse = ", ")))
    out <- x[sam]
    return(out)
  } else {
    x <- as.data.frame(x)

    if (is.null(rownames(x))) {
      if (nrow(x) != length(sam)) stop(sprintf("[%s] row number mismatch: %d vs %d and no rownames to align.", what, nrow(x), length(sam)))
      rownames(x) <- sam
    }
    if (!all(sam %in% rownames(x))) stop(sprintf("[%s] lack sample: %s", what, paste(setdiff(sam, rownames(x)), collapse = ", ")))
    out <- x[sam, , drop = FALSE]
    return(out)
  }
}


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


.factorize_with_explicit_NA <- function(x) {
  x_chr <- as.character(x)
  x_chr[!nzchar(x_chr) | is.na(x_chr)] <- "NA"
  factor(x_chr)
}


## ==== PATCH: limma t with covariates (guard p≈n) ====
limma_t_with_covars <- function(mat, grp2, batch = NULL, covars = NULL) {
  stopifnot(!is.null(colnames(mat)))
  sam <- colnames(mat)

  if (is.null(names(grp2))) names(grp2) <- sam
  grp2 <- factor(as.character(grp2[sam]), levels = c("Low", "High"))
  if (anyNA(grp2)) stop("grp2 still has NAs; please process the grouped samples first.")

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
    if (length(drop)) log_msg("  [design] Zero/near-zero mutation → Delete: %s", paste(drop, collapse = ", "))
    design <- design[, var_ok, drop = FALSE]
  }

  if (!all(c("Low", "High") %in% colnames(design))) {
    log_msg("  [Design] Low/High incomplete, Skip this time limma")
    return(NULL)
  }


  q <- qr(design)
  keep_idx <- q$pivot[seq_len(q$rank)]
  keep_nms <- colnames(design)[keep_idx]
  removed <- setdiff(colnames(design), keep_nms)
  removed <- setdiff(removed, c("Low", "High"))
  if (length(removed)) log_msg("[Design] QR code de-collinearity: Remove %s", paste(removed, collapse = ", "))
  design <- design[, keep_nms, drop = FALSE]


  if (ncol(design) > (ncol(mat) - 2)) {
    keep <- intersect(
      colnames(design),
      c(
        "Low", "High",
        grep("^b", colnames(design), value = TRUE),
        grep("^SV", colnames(design), value = TRUE),
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
  if (sum(is.finite(x)) < 4) stop("The continuous predictor has too few valid samples.")


  DF <- data.frame(x = x, row.names = sam, check.names = FALSE)
  if (!is.null(batch)) {
    b <- .align_to_samples(batch, sam, what = "batch")
    DF[["batch"]] <- droplevels(b)
  }
  if (!is.null(covars)) {
    cv <- .align_to_samples(covars, sam, what = "covars")
    cv <- as.data.frame(cv, check.names = FALSE)
    cv <- coerce_covariates_safely(cv)
    for (cn in colnames(cv)) DF[[cn]] <- cv[[cn]]
  }
  ok <- stats::complete.cases(DF)
  if (sum(ok) < (2L * min_per_group)) stop(sprintf("Insufficient available samples (%d < %d)", sum(ok), 2L * min_per_group))

  DF <- DF[ok, , drop = FALSE]
  M <- as.matrix(mat[, rownames(DF), drop = FALSE])
  des <- stats::model.matrix(~ 1 + ., data = DF, na.action = stats::na.fail)
  stopifnot(nrow(des) == ncol(M))

  fit <- limma::lmFit(M, des)

  coef_idx <- which(colnames(des) == "x")
  if (length(coef_idx) != 1) stop("x cannot be found.")

  C <- matrix(0, nrow = ncol(des), ncol = 1, dimnames = list(colnames(des), "x"))
  C["x", 1] <- 1
  fit2 <- limma::contrasts.fit(fit, C)
  fit2 <- limma::eBayes(fit2)

  t <- fit2$t[, 1]
  t <- sort(t[is.finite(t)], decreasing = TRUE)
  attr(t, "design_cols") <- colnames(des)
  t
}




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
    log_msg("  [TP53] data_mutations.txt not found; the entire batch is considered wild-type/ALL and available.")
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
    log_msg("[TP53] data_mutations.txt Missing field: %s, Treat as wild-type", paste(miss, collapse = ", "))
    return(status)
  }

  tp53_df <- subset(mutation_df, Hugo_Symbol == "TP53",
    select = c("Variant_Classification", "Tumor_Sample_Barcode")
  )
  if (!nrow(tp53_df)) {
    return(status)
  }

  vc_norm <- normalize_vc(tp53_df$Variant_Classification)


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


dir.create(file.path("run_info", "tp53_status"), recursive = TRUE, showWarnings = FALSE)


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


summarize_tp53_counts_for_dataset <- function(ds_dir) {
  ds_id <- basename(ds_dir)


  M <- try(load_phospho_matrix_from_dataset_dir(ds_dir), silent = TRUE)
  if (inherits(M, "try-error")) {
    log_msg("[TP53-audit] %s: Unable to read matrix, skipped.", ds_id)
    return(NULL)
  }
  sample_ids <- colnames(M)
  sid_up <- toupper(sample_ids)
  n_all <- length(sample_ids)


  status <- get_tp53_status(ds_dir, sample_ids)
  tb_bin <- table(factor(status, levels = c("TP53_wild_type", "TP53_mutant")))
  bin_row <- data.frame(
    dataset = ds_id,
    in_matrix_n = n_all,
    WT_n = as.integer(tb_bin["TP53_wild_type"]),
    MUT_n = as.integer(tb_bin["TP53_mutant"]),
    stringsAsFactors = FALSE
  )


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


        tp53 <- tp53[tp53$Tumor_Sample_Barcode %in% sid_up, , drop = FALSE]


        class_counts <- tp53 |>
          dplyr::group_by(Variant_Classification) |>
          dplyr::summarise(sample_n = dplyr::n_distinct(Tumor_Sample_Barcode), .groups = "drop") |>
          dplyr::arrange(dplyr::desc(sample_n), Variant_Classification)


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


summarize_tp53_counts_all_datasets <- function(dataset_dirs) {
  all_bin <- list()
  all_class <- list()
  k <- 1L
  j <- 1L
  for (ds in names(dataset_dirs)) {
    ds_dir <- dataset_dirs[[ds]]
    if (!dir.exists(ds_dir)) next
    log_msg("[TP53-audit] start: %s", ds)
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
    log_msg("[TP53-audit] Write out: tp53_binary_counts_by_dataset.csv")
  }
  if (length(all_class)) {
    class_df <- dplyr::bind_rows(all_class)
    data.table::fwrite(class_df, file.path("run_info", "tp53_status", "tp53_class_sample_counts_long.csv"))
    log_msg("[TP53-audit] Write out: tp53_class_sample_counts_long.csv")
  }
}


## summarize_tp53_counts_all_datasets(dataset_dirs)

## === Predictor analysis function ===
run_predictor_analyses <- function(
  predictor_name,
  predictor_vec,
  exclude_genes = NULL,
  ds_id, ds_dir,
  mat0,
  mat,
  out_root,
  genesets_by_group,
  batch_all = NULL,
  purity_all = NULL,
  sa_all = NULL,
  tp53_num_all = NULL,
  is_ALL = FALSE,
  additional_covars = NULL
) {
  opt <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
  min_per_group <- opt("min_per_group", 8L)
  minSize <- opt("minSize", 15L)
  maxSize <- opt("maxSize", 500L)
  fgsea_eps <- opt("fgsea_eps", 0)
  USE_AGE_MISSING_INDICATOR <- isTRUE(opt("USE_AGE_MISSING_INDICATOR", FALSE))

  logf <- function(...) if (exists("log_msg", mode = "function")) try(log_msg(...), silent = TRUE)
  sfn <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

  stopifnot(!is.null(genesets_by_group), length(genesets_by_group) > 0)


  if (is.null(names(predictor_vec))) names(predictor_vec) <- colnames(mat)
  keep <- intersect(colnames(mat), names(predictor_vec))
  pred_all <- suppressWarnings(as.numeric(predictor_vec[keep]))
  names(pred_all) <- keep
  fin <- is.finite(pred_all)
  if (sum(fin) < (2L * min_per_group)) {
    logf("[%s] predictor - insufficient non-NA samples (%d < %d), skip", predictor_name, sum(fin), 2L * min_per_group)
    return(invisible(NULL))
  }
  sample_order <- keep[fin]
  pred <- pred_all[fin]
  stopifnot(identical(names(pred), sample_order))


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
    if (!is.null(additional_covars)) {
      ac <- as.data.frame(additional_covars, check.names = FALSE)
      common_ac <- intersect(rownames(ac), so)
      if (length(common_ac) > 0) {
        ac <- ac[common_ac, , drop = FALSE]
        for (cn in colnames(ac)) {
          df[common_ac, cn] <- ac[common_ac, cn]
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


  ## --------coverage--------
  base_thr <- c(purity = 0.60, sex = 0.80, age = 0.80)
  present <- colnames(df_covars0)
  extra <- setdiff(present, names(base_thr))
  if (length(extra)) base_thr <- c(base_thr, stats::setNames(rep(min(base_thr), length(extra)), extra))


  pick_covars_df <- function(label) {
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
          logf("[covars-%s] Return shape/alignment unclear, discard.", label)
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

    if (!is.null(dl_sel) && nrow(dl_sel)) {
      dl_sel$dataset <- ds_id
      dl_sel$stratum <- basename(out_root)
      dl_sel$subunit <- predictor_name
      dl_sel$pass <- label
      dl_sel$samples <- length(sample_order)

      out_dir <- file.path("run_info", "covars_audit")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      fp <- file.path(out_dir, "covariate_drop_reasons.csv")

      data.table::fwrite(dl_sel, fp, append = file.exists(fp))
    }

    if (!is.null(X) && !is.null(additional_covars)) {
      ac_names <- colnames(as.data.frame(additional_covars))
      # Ensure these columns are present in X if they are in df_covars0 (which they should be)
      for (acn in ac_names) {
        if (acn %in% colnames(df_covars0) && !(acn %in% colnames(X))) {
          X[[acn]] <- df_covars0[rownames(X), acn]
        }
      }
    }

    X
  }
  X_raw_cov <- pick_covars_df("limma-cont:RAW")
  if (!is.null(X_raw_cov)) X_raw_cov <- coerce_covariates_safely(X_raw_cov)
  X_ba_cov <- pick_covars_df("limma-cont:base")
  if (!is.null(X_ba_cov)) X_ba_cov <- coerce_covariates_safely(X_ba_cov)

  DF_raw <- data.frame(predictor = as.numeric(pred), row.names = sample_order, check.names = FALSE)
  if (!is.null(X_raw_cov)) for (nm in colnames(X_raw_cov)) DF_raw[[nm]] <- X_raw_cov[[nm]]
  ok_raw <- stats::complete.cases(DF_raw)
  if (sum(ok_raw) < (2L * min_per_group)) {
    logf("[%s|RAW] Insufficient available samples (%d < %d), skip", predictor_name, sum(ok_raw), 2L * min_per_group)
    return(invisible(NULL))
  }
  DF_raw <- DF_raw[ok_raw, , drop = FALSE]

  DF_raw$predictor <- DF_raw$predictor - mean(DF_raw$predictor, na.rm = TRUE)
  M_raw <- as.matrix(mat[, rownames(DF_raw), drop = FALSE])
  des_raw <- stats::model.matrix(~ 1 + ., data = DF_raw, na.action = stats::na.fail)
  stopifnot(nrow(des_raw) == ncol(M_raw))

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
  )

  audit_design_alignment(
    tag = sprintf("%s_%s_%s|RAW", ds_id, basename(out_root), predictor_name),
    samples = colnames(M_raw),
    mod_interest = des_raw,
    mod_nuisance = NULL,
    out_dir = file.path("run_info", "covars_audit")
  )


  DF_ba <- data.frame(predictor = as.numeric(pred), row.names = sample_order, check.names = FALSE)
  if (!is.null(X_ba_cov)) for (nm in colnames(X_ba_cov)) DF_ba[[nm]] <- X_ba_cov[[nm]]
  if (!is.null(batch_all)) DF_ba$batch <- droplevels(as.factor(batch_all[sample_order]))

  use_cols <- unique(c(
    "predictor",
    intersect(c("purity", "sex", "age", "batch", "TP53_mutant"), colnames(DF_ba)),
    grep("^SV\\d+$|^PC\\d+$", colnames(DF_ba), value = TRUE)
  ))
  use_cols <- use_cols[use_cols %in% colnames(DF_ba)]


  ok_rows <- complete.cases(DF_ba[, use_cols, drop = FALSE])
  n_ok <- sum(ok_rows)

  if (n_ok < 16L) {
    logf("[%s|BatchAdj] Insufficient available samples (%d < 16) → skip", predictor_name, n_ok)
    return(invisible(NULL))
  }


  DF_ba <- DF_ba[ok_rows, , drop = FALSE]


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
    covars     = DF_ba[, intersect(c("purity", "sex", "age", "batch"), colnames(DF_ba)), drop = FALSE],
    sv         = DF_ba[, grep("^SV\\d+$", colnames(DF_ba), perl = TRUE), drop = FALSE],
    tech       = DF_ba[, grep("^PC\\d+$", colnames(DF_ba), perl = TRUE), drop = FALSE]
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


  if (isTRUE(is_ALL) && !is.null(tp53_num_all) &&
    length(genesets_by_group) &&
    any(vapply(genesets_by_group, function(x) length(x) > 0L, logical(1)))) {
    ## ---------- BatchAdj branch（interaction + diffcorr）----------
    try(
      {
        tp53_ba_num <- suppressWarnings(as.numeric(tp53_num_all[rownames(DF_ba)]))
        tp53_ba_fac <- factor(ifelse(tp53_ba_num == 1, "MT", "WT"), levels = c("WT", "MT"))
        covars_only_ba <- DF_ba[, setdiff(
          colnames(DF_ba),
          c("predictor", "TP53_mutant", "TP53_status", "TP53")
        ), drop = FALSE]

        df_int_ba <- data.frame(
          pred = DF_ba$predictor, tp53 = tp53_ba_fac,
          covars_only_ba, row.names = rownames(DF_ba), check.names = FALSE
        )
        des_int_ba <- stats::model.matrix(~ pred * tp53 + ., data = df_int_ba)
        if (.RUN_LIMMA) {
          fit_int_ba <- limma::eBayes(limma::lmFit(M_ba, des_int_ba))
          coef_int <- "pred:tp53MT"
          if (coef_int %in% colnames(coef(fit_int_ba))) {
            tt <- limma::topTable(fit_int_ba, coef = coef_int, number = nrow(M_ba), sort.by = "none")
            tvec <- if ("t" %in% names(tt)) tt$t else tt$logFC
            names(tvec) <- rownames(tt)

            ranks <- ._ensure_stats_names(tvec, rownames(M_ba))
            ranks <- ._finite_rank_stats(ranks, label = paste0("A2-", predictor_name, "-TP53_interaction"))
            for (grp in names(genesets_by_group)) {
              pw <- genesets_by_group[[grp]]
              if (is.null(pw) || !length(pw)) next
              out_dir_int_ba <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, safe_fs_name(grp), "BatchAdj", ds_id, basename(out_root), predictor_name)
              dir.create(out_dir_int_ba, recursive = TRUE, showWarnings = FALSE)

              # ---- fgsea ----
              res_fg <- .gsea_from_ranks(
                pathways = pw,
                stats    = ranks,
                minSize  = opt("minSize", 15L),
                maxSize  = opt("maxSize", 500L),
                gsea_eps = 1e-10,
                label    = paste0("A2-", predictor_name, "-TP53_interaction")
              )
              data.table::fwrite(as.data.frame(res_fg), file.path(out_dir_int_ba, "GSEA_limma_interaction.csv"))
            }
          }
        }
        sWT <- which(tp53_ba_fac == "WT")
        sMT <- which(tp53_ba_fac == "MT")
        min_pairs <- 10L
        if (length(sWT) >= min_pairs && length(sMT) >= min_pairs) {
          orig_pred_ba <- pred[rownames(DF_ba)]
          rho_wt <- apply(M_ba, 1, function(g) suppressWarnings(cor(g[sWT], orig_pred_ba[sWT], method = "spearman", use = "pairwise.complete.obs")))
          rho_mt <- apply(M_ba, 1, function(g) suppressWarnings(cor(g[sMT], orig_pred_ba[sMT], method = "spearman", use = "pairwise.complete.obs")))
          z_wt <- atanh(pmax(pmin(rho_wt, 0.999999), -0.999999))
          z_mt <- atanh(pmax(pmin(rho_mt, 0.999999), -0.999999))
          n_wt <- sum(is.finite(orig_pred_ba[sWT]))
          n_mt <- sum(is.finite(orig_pred_ba[sMT]))
          if (n_wt >= 4 && n_mt >= 4) {}
        }
      },
      silent = TRUE
    )
  }


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


    minSize <- opt("minSize", 15L)
    maxSize <- opt("maxSize", 500L)

    orig_n <- length(pw)
    pw_use <- ._intersect_and_filter_pathways(pw, names(stats), minSize = minSize, maxSize = maxSize)
    message(sprintf("[gsea] (gsea_from) |H| original=%d → use=%d", orig_n, length(pw_use)))
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
        message(sprintf("[gsea] Multilevel failed: %s → Use fgseaSimple instead.", conditionMessage(e)))
        suppressWarnings(fgsea::fgseaSimple(
          pathways = pw_use, stats = stats,
          nperm = 10000, minSize = minSize, maxSize = maxSize
        ))
      }
    )
    res
  }
  ## ---------- Continuous Analysis (BatchAdj only) ----------
  if (.RUN_LIMMA) {
    for (grp_name in names(genesets_by_group)) {
      pw <- genesets_by_group[[grp_name]]

      out_root_coll <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, safe_fs_name(grp_name), ds_id, basename(out_root))

      out_dir_A2 <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root), predictor_name)
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
      }
    }
  }
}





audit_csn_score_feasibility <- function(ds_id, stratum, mat0, prot0, present_sub,
                                        min_members = 5L,
                                        pca_min_samples = 10L,
                                        min_per_group = 8L,
                                        out_dir = file.path("run_info", "csn_score_audit")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  n_samples <- ncol(prot0)
  if (!is.matrix(mat0) || is.null(rownames(mat0)) || is.null(colnames(mat0))) {
    log_msg("[CSN-audit] %s | %s: mat0 The structure is incomplete, so it's skipped.", ds_id, stratum)
    return(invisible(NULL))
  }
  present_sub <- intersect(present_sub, rownames(prot0))
  if (!length(present_sub) || n_samples == 0) {
    log_msg("[CSN-audit] %s | %s: No CSN subunit available or 0 samples, skipped.", ds_id, stratum)
    return(invisible(NULL))
  }


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


  sample_counts <- colSums(is.finite(prot0[present_sub, , drop = FALSE]))
  enough <- sample_counts >= min_members
  n_enough <- sum(enough)


  csn_score <- build_csn_score(prot0,
    subunits = present_sub,
    combine_7AB = TRUE, min_members = min_members
  )
  csn_nonNA <- sum(is.finite(csn_score))
  csn_can_pca <- (n_enough >= pca_min_samples) && (csn_nonNA >= pca_min_samples)


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


  tag <- paste(ds_id, stratum, sep = "_")
  fp_sub <- file.path(out_dir, sprintf("%s_subunit_coverage.csv", tag))
  data.table::fwrite(sub_tbl, fp_sub)


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


run_one_stratum <- function(ds_id, ds_dir, mat0_full, sample_keep, out_root, genesets_by_group) {
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  log_msg("  -- stratum：%s | N(sample_keep)=%d", basename(out_root), length(sample_keep))


  keep <- intersect(colnames(mat0_full), sample_keep)
  if (length(keep) < 4) {
    log_msg("[Skip] Insufficient sample size: %d", length(keep))
    return(invisible(NULL))
  }
  mat0 <- mat0_full[, keep, drop = FALSE]
  mx <- suppressWarnings(max(mat0, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)

  prot0_full <- load_matrix_from_dataset_dir(ds_dir)

  sam_keep <- intersect(colnames(mat0), colnames(prot0_full))
  prot0 <- prot0_full[, sam_keep, drop = FALSE]
  mat0 <- mat0[, sam_keep, drop = FALSE]

  mxp <- suppressWarnings(max(prot0, na.rm = TRUE))
  if (is.finite(mxp) && mxp > 100) prot0 <- log2(prot0 + 1)


  mat <- impute_and_filter(mat0, min_frac = min_frac_complete)
  present_sub <- intersect(csn_subunits, rownames(prot0))
  if (!length(present_sub)) {
    log_msg("[Skip] This layer has no CSN subunits.")
    return(invisible(NULL))
  }
  audit_csn_score_feasibility_safe(
    ds_id = ds_id,
    stratum = basename(out_root),
    mat0 = mat0,
    prot0 = prot0,
    present_sub = present_sub,
    min_members = 5L,
    pca_min_samples = 10L,
    min_per_group = min_per_group,
    out_dir = file.path("run_info", "csn_score_audit")
  )


  is_ALL <- identical(basename(out_root), "ALL")
  tp53_num_all <- NULL
  if (is_ALL) {
    tp53_status_all <- get_tp53_status(ds_dir, colnames(mat0))
    tp53_num_all <- as.numeric(tp53_status_all == "TP53_mutant")
    names(tp53_num_all) <- colnames(mat0)
  }


  bi_all <- get_batch_factor_phospho(ds_dir, colnames(mat0))
  batch_all <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL
  purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
  sa_all <- get_sex_age_covariates(ds_dir, colnames(mat0))
  sa_all_limma <- sa_all
  if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing", "age_z_imputed") %in% colnames(sa_all))) {
    keep_cols <- intersect(c("sex", "age", "age_missing", "age_z_imputed"), colnames(sa_all))
    sa_all_limma <- sa_all[, keep_cols, drop = FALSE]
  }


  {
    present_sub <- intersect(csn_subunits, rownames(prot0))
    if (!length(present_sub)) {
      log_msg("[Skip] This layer has no CSN subunits.")
      return(invisible(NULL))
    }


    is_ALL <- identical(basename(out_root), "ALL")
    tp53_num_all <- if (is_ALL) {
      status_all <- get_tp53_status(ds_dir, colnames(mat0))
      setNames(as.numeric(status_all == "TP53_mutant"), colnames(mat0))
    } else {
      NULL
    }

    bi_all <- get_batch_factor_phospho(ds_dir, colnames(mat0))
    batch_all <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL
    purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
    sa_all <- get_sex_age_covariates(ds_dir, colnames(mat0))
    sa_all_limma <- sa_all
    if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing", "age_z_imputed") %in% colnames(sa_all))) {
      keep_cols <- intersect(c("sex", "age", "age_missing", "age_z_imputed"), colnames(sa_all))
      sa_all_limma <- sa_all[, keep_cols, drop = FALSE]
    }


    for (su in present_sub) {
      run_predictor_analyses(
        predictor_name = su,
        predictor_vec = prot0[su, ],
        exclude_genes = su,
        ds_id = ds_id, ds_dir = ds_dir,
        mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
        out_root = out_root,
        genesets_by_group = genesets_by_group,
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
        tp53_num_all = tp53_num_all, is_ALL = is_ALL
      )
    }


    csn_score <- build_csn_score_safe(
      prot0,
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
      log_msg("[CSN_SCORE] Insufficient non-NA samples, skipped.")
    }


    min_n_resid <- min_per_group
    csn_nonNA <- sum(is.finite(csn_score))
    if (csn_nonNA < min_n_resid) {
      log_msg("[RESIDUAL] CSN score: Insufficient non-NA samples (%d < %d), skip residual_* in the entire batch.", csn_nonNA, min_n_resid)
    } else {
      base_covars_all <- data.frame(
        purity = as.numeric(purity_all[colnames(mat0)]),
        sex = as.numeric(sa_all[colnames(mat0), "sex"]),
        age = as.numeric(sa_all[colnames(mat0), "age"]),
        row.names = colnames(mat0), check.names = FALSE
      )

      if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing", "age_z_imputed") %in% colnames(sa_all))) {
        base_covars_all$age_missing <- as.numeric(sa_all[colnames(mat0), "age_missing"])
        base_covars_all$age_z_imputed <- as.numeric(sa_all[colnames(mat0), "age_z_imputed"])
      }
      if (is_ALL && !is.null(tp53_num_all)) {
        base_covars_all$TP53_mutant <- as.numeric(tp53_num_all[colnames(mat0)])
      }

      for (su in present_sub) {
        # [MODIFIED] One-step analysis: Do not manually residualize.
        # Instead, pass CSN_SCORE (and other covars) into the design matrix.

        # We still check for sufficient samples based on CSN_SCORE availability
        if (sum(is.finite(csn_score)) < min_n_resid) {
          log_msg("[RESIDUAL_%s] Insufficient non-NA samples (CSN_SCORE), skipped.", su)
          next
        }

        # Prepare additional covariates (CSN_SCORE)
        # Note: batch, purity, sex, age, tp53 are already handled by run_predictor_analyses via its arguments
        # We only need to pass CSN_SCORE here.
        add_cov <- data.frame(CSN_SCORE = csn_score, row.names = names(csn_score), check.names = FALSE)

        # Use the RAW subunit expression as the predictor vector
        # But name the predictor "RESIDUAL_su" to keep output filenames unchanged
        run_predictor_analyses(
          predictor_name = paste0("RESIDUAL_", su),
          predictor_vec = prot0[su, ],
          exclude_genes = su,
          ds_id = ds_id, ds_dir = ds_dir,
          mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
          out_root = out_root,
          genesets_by_group = genesets_by_group,
          batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
          tp53_num_all = tp53_num_all, is_ALL = is_ALL,
          additional_covars = add_cov
        )
      }
    }
  }


  present_sub <- intersect(csn_subunits, rownames(prot0))
  sum_units <- c(present_sub, "CSN_SCORE", paste0("RESIDUAL_", present_sub))


  for (grp_name in names(genesets_by_group)) {
    out_root_coll <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, safe_fs_name(grp_name), ds_id, basename(out_root))
    ver_root <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, safe_fs_name(grp_name), "RAW", ds_id, basename(out_root))
    summarize_all_groups(
      out_root = ver_root,
      csn_subunits = sum_units,
      genesets_by_group = setNames(list(genesets_by_group[[grp_name]]), grp_name),
      stat_tags = c("GSEA_limma_t_cont", "GSEA_limma_interaction")
    )
  }


  for (grp_name in names(genesets_by_group)) {
    out_root_coll <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, safe_fs_name(grp_name), ds_id, basename(out_root))
    ver_root <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root))
    summarize_all_groups(
      out_root = ver_root,
      csn_subunits = sum_units,
      genesets_by_group = setNames(list(genesets_by_group[[grp_name]]), grp_name),
      stat_tags = {
        st <- c("GSEA_limma_t_cont", "GSEA_limma_interaction")
        if (isTRUE(RUN_HILO_SENSITIVITY)) st <- c(st, "GSEA_limma_t_hilo")
        st
      }
    )
  }

  invisible(NULL)
}


inspect_plex <- function(ds_dir, col = "TMT_PLEX",
                         pipe_policy = BATCH_PIPE_POLICY,
                         min_per_level = BATCH_MIN_PER_LEVEL) {
  mat0 <- load_phospho_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(mat0)


  meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  stopifnot(file.exists(meta_fp))
  meta <- suppressMessages(readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
  id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"), names(meta))
  stopifnot(length(id_cols) > 0)
  meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]


  raw <- as.character(meta[[col]])
  clean <- sanitize_batch_levels(raw, pipe_policy = pipe_policy, min_per_level = min_per_level)


  cat("\n==== ", basename(ds_dir), " | column: ", col, " ====\n", sep = "")
  cat("[Original PLEX levels]:\n")
  print(sort(table(raw), decreasing = TRUE))
  cat("[The first few primitive values containing '|']:")
  print(utils::head(unique(raw[grepl("\\|", raw %||% "")]), 10))

  cat("[Cleaned PLEX levels] (pipe_policy =", pipe_policy,
    ", min_per_level = ", min_per_level, ")：\n",
    sep = ""
  )
  print(sort(table(clean, useNA = "ifany"), decreasing = TRUE))


  df_map <- data.frame(SAMPLE_ID = sample_ids, raw_plex = raw, clean_plex = as.character(clean))
  cat("[Sample Comparison Table (First 10)]")
  print(utils::head(df_map, 10))

  invisible(list(
    raw_counts = sort(table(raw), decreasing = TRUE),
    clean_counts = sort(table(clean, useNA = "ifany"), decreasing = TRUE),
    map = df_map
  ))
}


ALWAYS_ADD_TECH_COVARS <- FALSE


USE_AGE_MISSING_INDICATOR <- FALSE


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

## PAAD
.median_from_semicolon <- function(x_chr) {
  vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
  vv <- vv[is.finite(vv)]
  if (!length(vv)) {
    return(NA_real_)
  }
  stats::median(vv)
}


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
    stop("[get_sex_age_covariates] Unable to establish a sample-patient correspondence (missing SAMPLE_ID/PATIENT_ID).")
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

  ## AGE
  age_col <- intersect(c("AGE", "AGE_AT_DIAGNOSIS", "AGE_AT_INDEX", "AGE_YEARS"), names(pat))[1]
  age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  age_missing <- age_z_imputed <- NULL
  if (!is.na(age_col)) {
    v <- suppressWarnings(as.numeric(pat[[age_col]]))
    names(v) <- as.character(pat[[pid_pat]])
    age_raw <- unname(v[map_pt[sample_ids]])
    age[] <- .z_no_impute(age_raw)
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
  df,
  sample_order,
  label = "covars",
  y = NULL,
  min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
  max_abs_cor = 0.30,
  min_pairs = 20L
) {
  logf <- function(...) if (exists("log_msg", mode = "function")) try(log_msg(...), silent = TRUE)

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


  # 1)
  if (is.null(rownames(df))) {
    logf("[covars-%s] df none rownames, skip", label)
    return(NULL)
  }
  so <- as.character(sample_order)
  logf("  [covars-%s] before align: C_dim=%d x %d; has_rownames=%s", label, NROW(df), NCOL(df), !is.null(rownames(df)))
  df <- df[so, , drop = FALSE]


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


  # 2)
  df <- as.data.frame(df, stringsAsFactors = TRUE, check.names = FALSE)

  # 3) coverage
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


  keep_cov <- rep(TRUE, ncol(df))
  names(keep_cov) <- colnames(df)


  # 4)
  if (!is.null(y)) {
    y <- suppressWarnings(as.numeric(y))

    skip_rho_gate <- c("purity", "age", "sex")

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


  # 5)
  stopifnot(nrow(df) == length(so), identical(rownames(df), so))

  attr(df, "drop_log") <- if (length(drop_log)) do.call(rbind, drop_log) else NULL
  df
}


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


audit_one_dataset_sa_batch <- function(ds_dir, pipe_policy = "NA", min_per_level = 2) {
  ds_id <- basename(ds_dir)
  # 1)
  m <- load_phospho_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(m)

  # 2)
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

  # 3) purity
  pur <- .audit_purity_for_dataset(ds_id, ds_dir, sample_ids)

  # 4) batch
  bi <- get_batch_factor_phospho(ds_dir, sample_ids,
    pipe_policy   = pipe_policy,
    min_per_level = min_per_level
  )
  batch_col <- if (!is.null(bi)) bi$name else "NONE"
  batch_levels <- if (!is.null(bi)) nlevels(bi$fac) else 0L
  batch_nonNA <- if (!is.null(bi)) mean(!is.na(bi$fac)) * 100 else NA_real_
  batch_sizes <- if (!is.null(bi)) .format_batch_sizes(bi$fac) else NA_character_


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


DO_AUDIT <- FALSE

if (DO_AUDIT) {
  results <- list()
  for (ds in names(dataset_dirs)) {
    ds_dir <- dataset_dirs[[ds]]
    if (!dir.exists(ds_dir)) {
      log_msg("Skip: Folder not found %s", ds_dir)
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
  hi_lo_quantile = hi_lo_quantile,
  minSize = minSize, maxSize = maxSize, fgsea_eps = fgsea_eps,
  min_per_group = min_per_group,
  MAKE_PLOTS = MAKE_PLOTS,
  datasets_root = datasets_root,
  dataset_ids = dataset_ids,
  strata = strata,
  geneset_groups_selected = GENESET_GROUPS_TO_RUN,
  msigdbr_version = msigdbr_version
), file.path("run_info", "run_manifest.yml"))



start_from <- "brca_cptac_2020"
ord <- dataset_ids
ix <- match(start_from, ord)
if (is.na(ix)) stop(sprintf("'%s' is not in dataset_ids", start_from))
dataset_dirs_run <- dataset_dirs[ord[ix:length(ord)]]

exists("coerce_covariates_safely")
getAnywhere("coerce_covariates_safely")
exists("opt")
getAnywhere("opt")




.should_run <- function(tag) {
  if (!exists("only_strata") || is.null(only_strata) || !length(only_strata)) {
    return(TRUE)
  }
  tag %in% only_strata
}


# --- [NEW | SITE-LEVEL] ---
load_phosphosite_matrix_from_dataset_dir <- function(ds_dir,
                                                     protein_adjust = TRUE,
                                                     min_nonNA_row_frac = 0) {
  fp <- file.path(ds_dir, "data_phosphoprotein_quantification.txt")
  if (!file.exists(fp)) stop("The phosphoprotein file could not be found: ", fp)
  log_msg("Read the phospho (site-level) matrix: data_phosphoprotein_quantification.txt")

  df <- suppressMessages(vroom::vroom(fp, delim = "\t", col_types = vroom::cols(
    .default = "d",
    ENTITY_STABLE_ID = "c",
    NAME = "c",
    DESCRIPTION = "c",
    GENE_SYMBOL = "c",
    PHOSPHOSITES = "c", PHOSPHOSITE = "c"
  )))
  df <- as.data.frame(df, check.names = FALSE)


  anno_cols <- intersect(c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "GENE_SYMBOL", "PHOSPHOSITES", "PHOSPHOSITE"), colnames(df))
  samp_cols <- setdiff(colnames(df), anno_cols)
  stopifnot(length(samp_cols) > 0)


  site_id <- .infer_site_id(df)
  keep <- !is.na(site_id) & nzchar(site_id)
  df <- df[keep, , drop = FALSE]
  site_id <- site_id[keep]

  M_site <- as.matrix(df[, samp_cols, drop = FALSE])
  storage.mode(M_site) <- "numeric"
  rownames(M_site) <- site_id
  colnames(M_site) <- samp_cols


  if (min_nonNA_row_frac > 0) {
    ok <- rowMeans(is.finite(M_site)) >= min_nonNA_row_frac
    M_site <- M_site[ok, , drop = FALSE]
  }


  if (isTRUE(protein_adjust)) {
    prot_fp <- file.path(ds_dir, "data_protein_quantification.txt")
    if (!file.exists(prot_fp)) stop("[phospho-site] protein_adjust=TRUE but the protein file cannot be found.")
    log_msg("Site-level protein quantification: Reading the protein matrix: data_protein_quantification.txt")

    prot_df <- suppressMessages(vroom::vroom(prot_fp, delim = "\t"))
    prot_df <- as.data.frame(prot_df, check.names = FALSE)


    gcol <- which(tolower(colnames(prot_df)) %in% tolower(c("Composite.Element.REF", "Gene", "GENE", "GENE_SYMBOL", "GENE_NAME")))[1]
    if (is.na(gcol)) stop("[phospho-site] protein_adjust=TRUE, but the protein file is missing the gene field.")
    prot_genes <- as.character(prot_df[[gcol]])
    prot_mat <- as.matrix(prot_df[, setdiff(colnames(prot_df), colnames(prot_df)[gcol]), drop = FALSE])
    rownames(prot_mat) <- prot_genes
    storage.mode(prot_mat) <- "numeric"


    gene_of_site <- as.character(df$GENE_SYMBOL)
    gene_of_site <- gene_of_site[keep]


    common <- intersect(colnames(M_site), colnames(prot_mat))
    if (length(common) >= 2) {
      M_site <- M_site[, common, drop = FALSE]
      prot_mat <- prot_mat[, common, drop = FALSE]


      idx <- match(gene_of_site, rownames(prot_mat))
      hasp <- !is.na(idx)
      if (any(hasp)) {
        M_site[hasp, ] <- M_site[hasp, , drop = FALSE] - prot_mat[idx[hasp], , drop = FALSE]
      }
    } else {
      log_msg("(WARNING) The intersection of phospho and protein samples is < 2; protein correction is skipped.")
    }
  }


  M_site
}
# --- [END NEW] ---


## =========================================================
## [COMMENTED OUT] Process the datasets to be run in this round in sequence (complete the main loop + define the sample sets for each layer).
## This section performs phosphoprotein-level GSEA using MsigDB (genesets_by_group).
## It has been disabled to focus on phosphosite-level PTM-SEA using PTMsigDB (genesets_by_group_ptm).
## =========================================================
# for (ds in names(dataset_dirs_run)) {
#   ds_dir <- dataset_dirs_run[[ds]]
#   log_msg("== Start Data Set: %s ==", ds)
#

#
#
#   mat0_full <- load_phospho_matrix_from_dataset_dir(ds_dir, protein_adjust = TRUE)
#   tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full))
#
#
#   samples_ALL <- colnames(mat0_full)
#   samples_MUT <- names(tp53_status)[tp53_status == "TP53_mutant"]
#   samples_WT <- names(tp53_status)[tp53_status == "TP53_wild_type"]
#
#
#   base_tp53_root <- if (is.null(COMBO_PREFIX)) file.path(ds_dir, "phospho_csn_gsea_results_TP53") else file.path(COMBO_PREFIX, ds, "phospho_csn_gsea_results_TP53")
#
#
#   if (.should_run("ALL")) {
#     run_one_stratum(
#       ds_id = ds, ds_dir = ds_dir,
#       mat0_full = mat0_full,
#       sample_keep = samples_ALL,
#       out_root = file.path(base_tp53_root, "ALL"),
#       genesets_by_group = genesets_by_group
#     )
#   }
#
#   if (.should_run("TP53_mutant")) {
#     run_one_stratum(
#       ds_id = ds, ds_dir = ds_dir,
#       mat0_full = mat0_full,
#       sample_keep = samples_MUT,
#       out_root = file.path(base_tp53_root, "TP53_mutant"),
#       genesets_by_group = genesets_by_group
#     )
#   }
#
#   if (.should_run("TP53_wild_type")) {
#     run_one_stratum(
#       ds_id = ds, ds_dir = ds_dir,
#       mat0_full = mat0_full,
#       sample_keep = samples_WT,
#       out_root = file.path(base_tp53_root, "TP53_wild_type"),
#       genesets_by_group = genesets_by_group
#     )
#   }
#
#   log_msg("== Completed dataset: %s (TP53 layered output → %s) ==", ds, base_tp53_root)
# }


## =========================================================
## [NEW | SITE-LEVEL PTM-SEA] Sequential processing of datasets (site level; PTMsigDB)
## =========================================================
for (ds in names(dataset_dirs_run)) {
  ds_dir <- dataset_dirs_run[[ds]]
  log_msg("== [PTM-SEA | site-level] Start Dataset: %s ==", ds)
  ## [SITE] Special Case: LUSC's Phospho column structure is incompatible; site-level entries will be skipped.
  if (identical(ds, "lusc_cptac_2021")) {
    log_msg("[PTM-SEA] skip: %s (LUSC phospho column structure incompatible; site-level not executed)", ds)
    next
  }
  ## ---- RUN_PASSES ----
  options(csn.run_passes = c("BatchAdj"))
  ## ---- [END RUN_PASSES] ----


  .old_min_frac_complete <- min_frac_complete
  min_frac_complete <- 0.30
  mat0_full_site <- load_phosphosite_matrix_from_dataset_dir(ds_dir, protein_adjust = TRUE)


  tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full_site))
  samples_ALL <- colnames(mat0_full_site)
  samples_MUT <- names(tp53_status)[tp53_status == "TP53_mutant"]
  samples_WT <- names(tp53_status)[tp53_status == "TP53_wild_type"]

  base_tp53_root_site <- file.path("phospho_site_combo_2", ds, "phospho_site_csn_gsea_results_TP53")


  if (length(genesets_by_group_ptm) == 0) {
    log_msg("[PTM-SEA] PTMsigDB not loaded, skipping dataset %s", ds)
    next
  }

  if (.should_run("ALL")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full_site,
      sample_keep = samples_ALL,
      out_root = file.path(base_tp53_root_site, "ALL"),
      genesets_by_group = genesets_by_group_ptm
    )
  }
  if (.should_run("TP53_mutant")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full_site,
      sample_keep = samples_MUT,
      out_root = file.path(base_tp53_root_site, "TP53_mutant"),
      genesets_by_group = genesets_by_group_ptm
    )
  }
  if (.should_run("TP53_wild_type")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full_site,
      sample_keep = samples_WT,
      out_root = file.path(base_tp53_root_site, "TP53_wild_type"),
      genesets_by_group = genesets_by_group_ptm
    )
  }


  min_frac_complete <- .old_min_frac_complete
  log_msg("== [PTM-SEA | site-level] Complete dataset: %s (TP53 layered output → %s) ==", ds, base_tp53_root_site)
}


## ---------- [NEW] WT vs MT ΔNES Aggregation (Summarizing from Existing Outputs) ----------
tp53_delta_nes_aggregate <- function(datasets_root,
                                     dataset_ids = NULL,
                                     versions = c("BatchAdj"),
                                     groups = names(genesets_by_group),
                                     stat_tags = c("GSEA_limma_t_cont")) {
  if (is.null(dataset_ids)) {
    dataset_ids <- list.dirs(datasets_root, full.names = FALSE, recursive = FALSE)
  }
  sfn <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  for (ds in dataset_ids) {
    base_dir <- file.path(datasets_root, ds, "phospho_csn_gsea_results_TP53")
    for (ver in versions) {
      subunits <- unique(unlist(lapply(groups, function(g) {
        b1 <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, safe_fs_name(g), ver, ds, "TP53_wild_type")
        b2 <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, safe_fs_name(g), ver, ds, "TP53_mutant")
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
            fp_wt <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, grp_safe, ver, ds, "TP53_wild_type", su, paste0(st, ".csv"))
            fp_mt <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, grp_safe, ver, ds, "TP53_mutant", su, paste0(st, ".csv"))
            if (!file.exists(fp_wt) || !file.exists(fp_mt)) next
            wt <- tryCatch(data.table::fread(fp_wt, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
            mt <- tryCatch(data.table::fread(fp_mt, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
            if (is.null(wt) || is.null(mt) || !nrow(wt) || !nrow(mt)) next
            wt <- setNames(as.data.frame(wt), tolower(names(wt)))
            mt <- setNames(as.data.frame(mt), tolower(names(mt)))
            # 標準化欄名
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
            out_dir <- file.path(COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection, grp_safe, ver, ds, "DeltaWT_MT", su)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            out_csv <- file.path(out_dir, paste0(st, "_deltaWTMT.csv"))
            data.table::fwrite(df, out_csv)
          }
        }
      }
      message(sprintf("[ΔNES] %s/%s Complete", ds, ver))
    }
  }
  invisible(TRUE)
}

## =========================================================
## [NEW-1] Producing meta-FDR across datasets (Stouffer → BH)

## Conditions: The above for loop has been completed.
## =========================================================
meta_fdr_stouffer(
  dataset_dirs = dataset_dirs_run,
  strata = strata,
  stat_tags = c(
    "GSEA_limma_t_cont",
    "GSEA_limma_interaction"
  ),
  groups = names(genesets_by_group_ptm),
  out_root = if (is.null(COMBO_PREFIX)) {
    file.path("PTMsigDB", "meta_fdr")
  } else {
    file.path(COMBO_PREFIX, "PTMsigDB", "meta_fdr")
  }
)


## [NEW-1b] WT vs MT ΔNES Aggregation (Compilation of Results from Landed Results)
try(tp53_delta_nes_aggregate(
  datasets_root = datasets_root,
  dataset_ids   = names(dataset_dirs_run),
  versions      = c("BatchAdj"),
  groups        = names(genesets_by_group_ptm),
  stat_tags     = c("GSEA_limma_t_cont")
), silent = TRUE)


## =========================================================
## Batch plotting GSEA Summary CSV of a single dataset/stratum into a heatmap
## =========================================================

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



options(
  csn_ylabel_wrap = 40,
  csn_left_margin_pts = 160,
  csn_ylabel_extra_width = 7.5,
  csn_ds_heatmap_width_in = NA_real_,
  csn_ds_heatmap_height_in = NA_real_,
  csn_meta_heatmap_width_in = NA_real_,
  csn_meta_heatmap_height_in = NA_real_,
  csn_heatmap_dpi = 600,
  csn_meta_heatmap_dpi = 600
)


.choose_plot_size <- function(n_path, n_preds, context = c("dataset", "meta")) {
  context <- match.arg(context)
  if (context == "dataset") {
    W_fix <- getOption("csn_ds_heatmap_width_in", NA_real_)
    H_fix <- getOption("csn_ds_heatmap_height_in", NA_real_)
    if (is.finite(W_fix) || is.finite(H_fix)) {
      #
      W_def <- max(8, n_preds * 0.45) + getOption("csn_ylabel_extra_width", 2.5)
      H_def <- max(6, n_path * 0.22)
      return(list(
        width  = if (is.finite(W_fix)) W_fix else W_def,
        height = if (is.finite(H_fix)) H_fix else H_def
      ))
    } else {
      return(list(
        width  = max(8, n_preds * 0.45) + getOption("csn_ylabel_extra_width", 2.5),
        height = max(6, n_path * 0.22)
      ))
    }
  } else {
    W_fix <- getOption("csn_meta_heatmap_width_in", NA_real_)
    H_fix <- getOption("csn_meta_heatmap_height_in", NA_real_)
    if (is.finite(W_fix) || is.finite(H_fix)) {
      W_def <- max(8, n_preds * 0.45) + getOption("csn_ylabel_extra_width", 2.5)
      H_def <- max(6, n_path * 0.22)
      return(list(
        width  = if (is.finite(W_fix)) W_fix else W_def,
        height = if (is.finite(H_fix)) H_fix else H_def
      ))
    } else {
      return(list(
        width  = max(8, n_preds * 0.45) + getOption("csn_ylabel_extra_width", 2.5),
        height = max(6, n_path * 0.22)
      ))
    }
  }
}


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
  L <- max(abs(df_long$NES), na.rm = TRUE)
  if (!is.finite(L) || L == 0) L <- 1
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = pathway, fill = NES)) +
    ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
    ggplot2::geom_point(
      data = df_long %>% dplyr::filter(is.finite(.data$padj), .data$padj < 0.05),
      ggplot2::aes(x = xpos, y = pathway),
      shape = 16, size = getOption("csn.dot_size_dataset", 1.6), color = "black", inherit.aes = FALSE
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
    ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = getOption("csn_ylabel_wrap", 40))) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0, size = 9, margin = ggplot2::margin(b = 8)),
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.9),
      legend.position = "right",
      plot.margin = ggplot2::margin(6, 12, 6, getOption("csn_left_margin_pts", 160)),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  n_path <- nlevels(df_long$pathway)
  size <- .choose_plot_size(n_path, length(present_preds), context = "dataset")
  list(plot = p, width = size$width, height = size$height)
}


.read_csv_safe <- function(path) {
  p <- path


  if (!file.exists(p) && .Platform$OS.type == "windows") {
    p2 <- gsub("/", "\\\\", p, fixed = TRUE)
    if (file.exists(p2)) p <- p2
  }
  if (!file.exists(p)) stop("The file does not exist:", p)

  fi <- suppressWarnings(file.info(p))
  if (isTRUE(fi$isdir)) stop("The target is not the file:", p)
  if (!is.na(fi$size) && fi$size == 0) stop("Empty file (size = 0)：", p)

  # 1)
  out <- try(suppressMessages(readr::read_csv(p, show_col_types = FALSE, progress = FALSE)), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(out)
  }

  # 2)
  out <- try(suppressMessages(readr::read_csv(
    p,
    locale = readr::locale(encoding = "UTF-8-BOM"),
    show_col_types = FALSE, progress = FALSE
  )), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(out)
  }

  # 3)
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

  # 4)
  out <- try(utils::read.csv(p, check.names = FALSE), silent = TRUE)
  if (inherits(out, "try-error")) stop(out)
  tibble::as_tibble(out)
}




.safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

# Target X-axis order
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
  if (is.na(path_col)) stop("The pathway field could not be found.", paste(path_candidates, collapse = ", "))
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
    if (!length(nes_cols)) stop("Wide table format: NES_* field not found.")
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
  if (!length(present_preds)) stop("The data does not specify any predictors.")
  df_long <- df_long %>% dplyr::filter(.data$predictor %in% present_preds)

  .coll_tok <- sub("^Summary_(.+?)_GSEA.*$", "\\1", basename(csv_file))
  if (!identical(toupper(.coll_tok), "H") && "NES_CSN_SCORE" %in% present_preds) {
    .top_n <- get0("DATASET_HEATMAP_TOP_N", ifnotfound = 25)
    .bot_n <- get0("DATASET_HEATMAP_BOTTOM_N", ifnotfound = 25)
    # optional pathway prefix filter for single-dataset heatmap
    .keep_prefix <- getOption("HEATMAP_PATHWAY_PREFIX", NULL)
    if (!is.null(.keep_prefix) && length(.keep_prefix) > 0) {
      .pat <- paste0("^(", paste(.keep_prefix, collapse = "|"), ")")
      df_long <- df_long %>% dplyr::filter(grepl(.pat, .data$pathway))
    }
    csn_tbl <- df_long %>% dplyr::filter(.data$predictor == "NES_CSN_SCORE")
    keep_up <- csn_tbl %>%
      dplyr::arrange(dplyr::desc(.data$NES)) %>%
      dplyr::slice_head(n = .top_n) %>%
      dplyr::pull(.data$pathway)
    keep_dn <- csn_tbl %>%
      dplyr::arrange(.data$NES) %>%
      dplyr::slice_head(n = .bot_n) %>%
      dplyr::pull(.data$pathway)
    keep <- unique(c(keep_up, keep_dn))
    if (length(keep)) df_long <- df_long %>% dplyr::filter(.data$pathway %in% keep)
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
  df_long <- df_long %>% dplyr::mutate(pathway = factor(.data$pathway, levels = rev(path_order)))

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

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = pathway, fill = NES)) +
    ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
    ggplot2::geom_point(
      data = df_long %>% dplyr::filter(is.finite(.data$padj), .data$padj < 0.05),
      ggplot2::aes(x = xpos, y = pathway),
      shape = 16, size = getOption("csn.dot_size_dataset", 1.6), color = "black", inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = col_low, mid = col_mid, high = col_high,
      limits = c(-L, L), oob = scales::squish, name = "NES"
    ) +
    ggplot2::scale_x_continuous(
      breaks = unname(pos_map[present_preds]),
      labels = present_preds,
      expand = ggplot2::expansion(mult = c(0.01, 0.01)),
      position = "top"
    ) +
    ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = getOption("csn_ylabel_wrap", 40))) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0, size = 9, margin = ggplot2::margin(b = 8)),
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.9),
      legend.position = "right",
      plot.margin = ggplot2::margin(6, 12, 6, getOption("csn_left_margin_pts", 160)),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  n_path <- nlevels(df_long$pathway)
  size <- .choose_plot_size(n_path, length(present_preds), context = "dataset")
  list(plot = p, width = size$width, height = size$height)
}


# Read the CSV and calculate the pathway sorting according to the rules of .make_heatmap_plot.
.compute_pathway_order_from_csv <- function(csv_file) {
  df_raw <- suppressMessages(.read_csv_safe(csv_file))
  path_candidates <- c("pathway", "Pathway", "term", "Term", "gs_name", "NAME", "set", "Set")
  path_col <- intersect(path_candidates, names(df_raw))[1]
  if (is.na(path_col)) stop("The pathway field could not be found: ", paste(path_candidates, collapse = ", "))

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
    if (!length(nes_cols)) stop("Wide table format: NES_* field not found.")
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

  ## —— Key: Use the "first 25 largest + last 25 smallest (based on CSN_SCORE)" of the ALL layer as the sole generatrix list —— ##
  .coll_tok <- sub("^Summary_(.+?)_GSEA.*$", "\\1", basename(csv_file))
  if ("NES_CSN_SCORE" %in% present_preds) {
    csn_tbl <- df_long %>%
      dplyr::filter(.data$predictor == "NES_CSN_SCORE") %>%
      dplyr::select(.data$pathway, .data$NES) %>%
      dplyr::distinct()
    .top_n <- get0("DATASET_HEATMAP_TOP_N", ifnotfound = 25L)
    .bot_n <- get0("DATASET_HEATMAP_BOTTOM_N", ifnotfound = 25L)
    # Non-H: Only keep top/bottom; H: All 50 are kept, equivalent to top/bottom = all.
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
      # H: Use all H values directly as the "ALL layer selection list"
      keep <- csn_tbl$pathway
      df_long <- df_long %>% dplyr::filter(.data$pathway %in% keep)
      csn_tbl <- csn_tbl %>% dplyr::filter(.data$pathway %in% keep)
    }
  }

  # Provide the "keep-only" order based on NES_CSN_SCORE from largest to smallest (or row average).
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


# Plot the map using the given y_order (pathway levels) (the rules are completely similar to .make_heatmap_plot).
.make_heatmap_plot_with_yorder <- function(csv_file, y_order) {
  L <- max(abs(df_long$NES), na.rm = TRUE)
  if (!is.finite(L) || L == 0) L <- 1
  col_low <- "#053061"
  col_mid <- "#FFFFFF"
  col_high <- "#67001F"

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = pathway, fill = NES)) +
    ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
    ggplot2::geom_point(
      data = df_long %>% dplyr::filter(is.finite(.data$padj), .data$padj < 0.05),
      ggplot2::aes(x = xpos, y = pathway),
      shape = 16, size = getOption("csn.dot_size_dataset", 1.6), color = "black", inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = col_low, mid = col_mid, high = col_high,
      limits = c(-L, L), oob = scales::squish, name = "NES"
    ) +
    ggplot2::scale_x_continuous(
      breaks = unname(pos_map[present_preds]),
      labels = present_preds,
      expand = ggplot2::expansion(mult = c(0.01, 0.01)),
      position = "top"
    ) +
    ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = getOption("csn_ylabel_wrap", 40))) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0, size = 9, margin = ggplot2::margin(b = 8)),
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.9),
      legend.position = "right",
      plot.margin = ggplot2::margin(6, 12, 6, getOption("csn_left_margin_pts", 160)),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  n_path <- nlevels(df_long$pathway)
  size <- .choose_plot_size(n_path, length(present_preds), context = "dataset")
  list(plot = p, width = size$width, height = size$height)
}


# Output only four formats "in-place" (for _ordered version)
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


  vidx <- which(parts_lc %in% c("raw", "batchadj"))
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
  keep <- grepl(paste0(sep_pat, "raw", sep_pat), files, ignore.case = TRUE) |
    grepl(paste0(sep_pat, "batchadj", sep_pat), files, ignore.case = TRUE)
  files <- files[keep]


  if (length(files)) {
    finfo <- suppressWarnings(file.info(files))
    files <- files[is.finite(finfo$size) & finfo$size > 0 & !finfo$isdir]
  }
  ## ====================================================

  files
}


# Main program: Run all datasets
run_gsea_heatmaps_for_all_datasets <- function(dataset_dirs_map = NULL,
                                               collect_root = "single_dataset_GSEA_heatmap",
                                               width_in = NULL, height_in = NULL, dpi = NULL) {
  old_w <- getOption("csn_ds_heatmap_width_in", NA_real_)
  old_h <- getOption("csn_ds_heatmap_height_in", NA_real_)
  old_dp <- getOption("csn_heatmap_dpi", 600)
  if (!is.null(width_in)) options(csn_ds_heatmap_width_in = width_in)
  if (!is.null(height_in)) options(csn_ds_heatmap_height_in = height_in)
  if (!is.null(dpi)) options(csn_heatmap_dpi = dpi)
  on.exit(options(
    csn_ds_heatmap_width_in = old_w,
    csn_ds_heatmap_height_in = old_h,
    csn_heatmap_dpi = old_dp
  ), add = TRUE)

  if (is.null(dataset_dirs_map)) {
    if (exists("dataset_dirs_run")) {
      dataset_dirs_map <- get("dataset_dirs_run", inherits = TRUE)
    } else if (exists("dataset_dirs")) {
      dataset_dirs_map <- get("dataset_dirs", inherits = TRUE)
    } else {
      stop("Please provide dataset_dirs_map or define dataset_dirs_run / dataset_dirs in your environment.")
    }
  }
  for (ds in names(dataset_dirs_map)) {
    ds_dir <- dataset_dirs_map[[ds]]
    scan_dir <- if (is.null(COMBO_PREFIX)) ds_dir else COMBO_PREFIX
    csvs <- .find_target_csvs_for_dataset(scan_dir)
    if (!is.null(COMBO_PREFIX)) {
      ds_name <- basename(ds_dir)
      sep_pat <- "[/\\\\]"
      pat_ds <- paste0(sep_pat, ds_name, sep_pat)
      csvs <- csvs[grepl(pat_ds, csvs, ignore.case = TRUE)]
      csvs <- csvs[grepl(paste0(sep_pat, "summary", sep_pat), tolower(csvs))]
      .grp_labels <- names(get0("genesets_by_group_ptm", ifnotfound = list()))
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
      if (!length(meta$dataset) || is.na(meta$dataset) || meta$dataset == "") meta$dataset <- ds
      if (grepl("GSEA_limma_interaction", basename(csv), fixed = TRUE) &&
        !identical(meta$stratum, "ALL")) {
        next
      }
      h <- try(.make_heatmap_plot(csv), silent = TRUE)
      if (inherits(h, "try-error")) {
        message("[heatmap] failed:", csv, " | ", as.character(h))
        next
      }
      .save_both_places(h$plot, h$width, h$height, csv, meta,
        collect_root = collect_root,
        dpi = getOption("csn_heatmap_dpi", 600)
      )
      message("[heatmap] Completed:", csv)

      if (grepl("gsea_limma_interaction", tolower(basename(csv)), fixed = TRUE) &&
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
              .save_near_only(h2$plot, h2$width, h2$height, csv, meta,
                suffix = "_ordered",
                dpi = getOption("csn_heatmap_dpi", 600)
              )
              message("[heatmap][interaction-ordered] Done:", csv)
            }
          }
        }
      }

      if (tolower(meta$stratum) %in% c("tp53_mutant", "tp53_wild_type")) {
        bn_lc <- tolower(basename(csv))
        is_cont <- grepl("gsea_limma_t_cont", bn_lc, fixed = TRUE)
        is_cont <- grepl("limma.*cont", bn_lc, fixed = FALSE)
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
                  .save_near_only(h2$plot, h2$width, h2$height, csv, meta,
                    suffix = "_ordered",
                    dpi = getOption("csn_heatmap_dpi", 600)
                  )
                  message("[heatmap][ordered] Complete (in ALL order): ", csv)
                }
              }
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

options(HEATMAP_PATHWAY_PREFIX = c("KINASE-PSP", "PATH-NP"))
## options(HEATMAP_PATHWAY_PREFIX = NULL)

options(csn.dot_size_dataset = 2.0) # Size of dots in a single dataset heatmap

## ---- Execution (Example) ----
run_gsea_heatmaps_for_all_datasets(
  dataset_dirs_map = if (exists("dataset_dirs_run")) dataset_dirs_run else NULL,
  collect_root = if (is.null(COMBO_PREFIX)) {
    "single_dataset_GSEA_heatmap"
  } else {
    file.path(COMBO_PREFIX, "single_dataset_GSEA_heatmap")
  },
  width_in = 10, height_in = 14, dpi = 600
)


## =========================================================
## [NEW] Summarize meta-FDR (Stouffer → BH) across subunits
##   - Input:  csn_gsea_pan_summary_TP53/meta_fdr/<STRATUM>/BatchAdj/<SUBUNIT>/<GROUP>/GSEA_limma_t_cont_meta_fdr.csv
##   - Output: .../summary/<STRATUM>/BatchAdj/<GROUP>/GSEA_limma_t_cont_meta_fdr/
##              ├─ Summary_<GROUP>_GSEA_limma_t_cont_meta_fdr_ALL.csv
##              ├─ Summary_<GROUP>_GSEA_limma_t_cont_meta_fdr_padjLT0.05.csv
##              └─ Summary_<GROUP>_GSEA_limma_t_cont_meta_fdr_padjLT0.25.csv
##   - Columns merged per subunit: Z_<SUBUNIT>, padj_meta_<SUBUNIT>
##   - Counting uses padj_meta thresholds; direction by sign(Z).
## =========================================================

.read_meta_fdr_table <- function(meta_root, stratum, version, subunit,
                                 group_name, stat_tag = "GSEA_limma_t_cont_meta_fdr") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first.")
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
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first.")
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install openxlsx")
  if (is.null(df) || !nrow(df)) {
    if (exists("log_msg", mode = "function")) try(log_msg("[Output omitted] {group_name} | {stat_tag} No results available"), silent = TRUE)
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

  ## ---- XLSX ----
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
    try(log_msg("[Completed Output] {group_name} | {stat_tag} -> {dirname(base)} (.csv + .xlsx)"), silent = TRUE)
  }
  invisible(NULL)
}


summarize_meta_fdr_across_subunits <- function(
  meta_root = if (is.null(COMBO_PREFIX)) "PTMsigDB/meta_fdr" else file.path(COMBO_PREFIX, "PTMsigDB", "meta_fdr"),
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  versions = c("BatchAdj"),
  genesets_by_group = genesets_by_group_ptm,
  stat_tag = "GSEA_limma_t_cont_meta_fdr"
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first.")
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


posthoc_summary_meta_fdr <- function() {
  summarize_meta_fdr_across_subunits(
    meta_root = if (is.null(COMBO_PREFIX)) "PTMsigDB/meta_fdr" else file.path(COMBO_PREFIX, "PTMsigDB", "meta_fdr"),
    strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
    versions = c("BatchAdj"),
    genesets_by_group = genesets_by_group_ptm,
    stat_tag = "GSEA_limma_t_cont_meta_fdr"
  )
  invisible(TRUE)
}


posthoc_summary_meta_fdr()




## =========================================================
## [NEW] Summarize meta-FDR for limma interaction across subunits
##   - Input:
##     csn_gsea_pan_summary_TP53/meta_fdr/<STRATUM>/BatchAdj/<SUBUNIT>/<GROUP>/GSEA_limma_interaction_meta_fdr.csv
##   - Output:
##     csn_gsea_pan_summary_TP53/meta_fdr/summary/<STRATUM>/BatchAdj/<GROUP>/GSEA_limma_interaction_meta_fdr/
##       ├─ Summary_<GROUP>_GSEA_limma_interaction_meta_fdr_ALL.csv
##       ├─ Summary_<GROUP>_GSEA_limma_interaction_meta_fdr_padjLT0.05.csv
##       ├─ Summary_<GROUP>_GSEA_limma_interaction_meta_fdr_padjLT0.25.csv
##       └─ Summary_<GROUP>_GSEA_limma_interaction_meta_fdr.xlsx
## =========================================================

posthoc_summary_meta_fdr_interaction <- function() {
  summarize_meta_fdr_across_subunits(
    meta_root = if (is.null(COMBO_PREFIX)) "PTMsigDB/meta_fdr" else file.path(COMBO_PREFIX, "PTMsigDB", "meta_fdr"),
    strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
    versions = c("BatchAdj"),
    genesets_by_group = genesets_by_group_ptm,
    stat_tag = "GSEA_limma_interaction_meta_fdr"
  )
  invisible(TRUE)
}

posthoc_summary_meta_fdr_interaction()






## =========================================================
## [NEW] Heatmaps for meta-FDR summary CSVs (Z + padj_meta)
##   - Inputs:
##     1) ALL   : Summary_H_GSEA_limma_t_cont_meta_fdr_ALL.csv
##     2) ALL   : Summary_H_GSEA_limma_interaction_meta_fdr_ALL.csv
##     3) TP53_mutant    : Summary_H_GSEA_limma_t_cont_meta_fdr_ALL.csv
##     4) TP53_wild_type : Summary_H_GSEA_limma_t_cont_meta_fdr_ALL.csv
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
    if (!is.na(fi$size) && fi$size == 0) stop("Empty file (size = 0):", p)

    # 1)
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

    # 3)
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

    # 4)
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
  if (is.na(path_col)) stop("The pathway field cannot be found at: ", csv_file)

  z_cols <- grep("^Z_", names(df), value = TRUE)
  padj_cols <- grep("^padj_meta_", names(df), value = TRUE)
  if (!length(z_cols)) stop("The Z_* column cannot be found at:", csv_file)

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

# Calculate the y-axis using meta file (Z) order: prioritize using Z_CSN_SCORE in descending order, otherwise use the row average Z.
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
  if (!length(present)) stop("The data does not specify any predictors: ", csv_file)
  df_long <- dplyr::filter(df_long, .data$predictor %in% present)

  .coll_tok <- sub("^Summary_(.+?)_GSEA.*$", "\\1", basename(csv_file))
  if (is.null(y_order)) {
    if (!identical(toupper(.coll_tok), "H") && "Z_CSN_SCORE" %in% present) {
      .top_n <- get0("PAN_HEATMAP_TOP_N", ifnotfound = 25)
      .bot_n <- get0("PAN_HEATMAP_BOTTOM_N", ifnotfound = 25)
      csn_tbl <- df_long %>% dplyr::filter(.data$predictor == "Z_CSN_SCORE")
      keep_up <- csn_tbl %>%
        dplyr::arrange(dplyr::desc(.data$Z)) %>%
        dplyr::slice_head(n = .top_n) %>%
        dplyr::pull(.data$pathway)
      keep_dn <- csn_tbl %>%
        dplyr::arrange(.data$Z) %>%
        dplyr::slice_head(n = .bot_n) %>%
        dplyr::pull(.data$pathway)
      keep <- unique(c(keep_up, keep_dn))
      if (length(keep)) df_long <- dplyr::filter(df_long, .data$pathway %in% keep)
    }
  }
  # optional pathway prefix filter for meta heatmap
  .keep_prefix <- getOption("HEATMAP_PATHWAY_PREFIX", NULL)
  if (!is.null(.keep_prefix) && length(.keep_prefix) > 0) {
    .pat <- paste0("^(", paste(.keep_prefix, collapse = "|"), ")")
    df_long <- df_long %>% dplyr::filter(grepl(.pat, .data$pathway))
  }
  if (is.null(y_order)) y_order <- .meta_compute_y_order(csv_file)
  df_long <- df_long %>% dplyr::filter(.data$pathway %in% y_order)
  y_levels <- y_order[y_order %in% df_long$pathway]
  if (!length(y_levels)) stop("y-order has no overlap with the data: ", csv_file)
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

  L <- max(abs(df_long$Z), na.rm = TRUE)
  if (!is.finite(L) || L == 0) L <- 1

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = pathway, fill = Z)) +
    ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
    ggplot2::geom_point(
      data = dplyr::filter(df_long, is.finite(.data$padj_meta), .data$padj_meta < 0.05),
      ggplot2::aes(x = xpos, y = pathway),
      shape = 16, size = getOption("csn.dot_size_meta", 1.6), color = "black", inherit.aes = FALSE
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
    ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = getOption("csn_ylabel_wrap", 40))) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0, size = 9, margin = ggplot2::margin(b = 8)),
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.9),
      legend.position = "right",
      plot.margin = ggplot2::margin(6, 12, 6, getOption("csn_left_margin_pts", 160)),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  n_path <- nlevels(df_long$pathway)
  size <- .choose_plot_size(n_path, length(present), context = "meta")
  list(plot = p, width = size$width, height = size$height)
}


.meta_parse_summary_path <- function(csv_file) {
  parts <- strsplit(normalizePath(csv_file, winslash = "/"), "/")[[1]]
  parts_lc <- tolower(parts)
  ver_idx <- which(parts_lc %in% c("raw", "batchadj"))
  version <- if (length(ver_idx)) parts[ver_idx[1]] else NA_character_
  st_idx <- which(parts_lc %in% c("all", "tp53_mutant", "tp53_wild_type"))
  stratum <- if (length(st_idx)) parts[st_idx[length(st_idx)]] else NA_character_
  list(version = version, stratum = stratum)
}


.meta_save_near <- function(p, width, height, csv_file, meta, suffix = "", dpi = 600) {
  out_dir <- dirname(csv_file)
  bn <- basename(csv_file)


  grp_tok <- sub("^Summary_([^_]+)_GSEA.*$", "\\1", bn, perl = TRUE)
  grp_tok <- .safe_fs(grp_tok)
  method_tok <- if (grepl("limma_interaction", bn, ignore.case = TRUE)) {
    "interaction"
  } else if (grepl("limma_t_cont", bn, ignore.case = TRUE)) {
    "tcont"
  } else if (grepl("limma_interaction", bn, ignore.case = TRUE)) {
    "interaction"
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


.meta_save_near_v2 <- function(p, width, height, csv_file, meta, suffix = "", dpi = 600) {
  out_dir <- dirname(csv_file)
  bn <- basename(csv_file) # e.g. Summary_C2__CP_KEGG_LEGACY_GSEA_limma_t_cont_meta_fdr_ALL.csv


  grp_tok <- sub("^Summary_(.+?)_GSEA.*$", "\\1", bn, perl = TRUE)
  grp_tok <- .safe_fs(grp_tok)

  method_tok <- if (grepl("limma_interaction", bn, ignore.case = TRUE)) {
    "interaction"
  } else if (grepl("limma_t_cont", bn, ignore.case = TRUE)) {
    "tcont"
  } else if (grepl("limma_interaction", bn, ignore.case = TRUE)) {
    "interaction"
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

# ---- v3: Absolutely short filenames & dedicated save function to avoid overwriting ----
.meta_save_near_v3_shortname <- function(p, width, height, csv_file, meta, suffix = "", dpi = 600) {
  out_dir <- dirname(csv_file)
  bn <- basename(csv_file)

  grp_tok <- sub("^Summary_(.+?)_GSEA.*$", "\\1", bn, perl = TRUE)
  grp_tok <- .safe_fs(grp_tok)
  method_tok <- if (grepl("limma_interaction", bn, ignore.case = TRUE)) {
    "interaction"
  } else if (grepl("limma_t_cont", bn, ignore.case = TRUE)) {
    "tcont"
  } else if (grepl("limma_interaction", bn, ignore.case = TRUE)) {
    "interaction"
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


.meta_find_target_csvs <- function(root = "phospho_csn_gsea_pan_summary_TP53/meta_fdr/summary") {
  if (!dir.exists(root)) {
    return(character(0))
  }
  pat_all_t <- "^Summary_.+_GSEA_limma_t_cont_meta_fdr_ALL\\.csv$"
  pat_all_i <- "^Summary_.+_GSEA_limma_interaction_meta_fdr_ALL\\.csv$"
  pat_mut_t <- "^Summary_.+_GSEA_limma_t_cont_meta_fdr_ALL\\.csv$"
  pat_wt_t <- "^Summary_.+_GSEA_limma_t_cont_meta_fdr_ALL\\.csv$"


  files <- list.files(root, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)

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

  # 剔除 size=0
  if (length(tgt)) {
    finfo <- suppressWarnings(file.info(tgt))
    tgt <- tgt[is.finite(finfo$size) & finfo$size > 0 & !finfo$isdir]
  }
  unique(tgt)
}


if (!exists(".find_target_csvs_for_dataset", mode = "function")) {
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
    keep <- grepl(paste0(sep_pat, "raw", sep_pat), files, ignore.case = TRUE) |
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


    vidx <- which(parts_lc %in% c("raw", "batchadj"))
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
run_meta_fdr_heatmaps <- function(
  root = if (is.null(COMBO_PREFIX)) {
    "PTMsigDB/meta_fdr/summary"
  } else {
    file.path(COMBO_PREFIX, "PTMsigDB", "meta_fdr", "summary")
  },
  width_in = NULL, height_in = NULL, dpi = NULL
) {
  old_w <- getOption("csn_meta_heatmap_width_in", NA_real_)
  old_h <- getOption("csn_meta_heatmap_height_in", NA_real_)
  old_dp <- getOption("csn_meta_heatmap_dpi", 600)
  if (!is.null(width_in)) options(csn_meta_heatmap_width_in = width_in)
  if (!is.null(height_in)) options(csn_meta_heatmap_height_in = height_in)
  if (!is.null(dpi)) options(csn_meta_heatmap_dpi = dpi)
  on.exit(options(
    csn_meta_heatmap_width_in = old_w,
    csn_meta_heatmap_height_in = old_h,
    csn_meta_heatmap_dpi = old_dp
  ), add = TRUE)

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
    message("[meta-heatmap] Target not found, skipped.")
    return(invisible(TRUE))
  }

  for (csv in csvs) {
    meta <- .meta_parse_summary_path(csv)
    bn_lc <- tolower(basename(csv))
    is_interaction <- grepl("gsea_limma_interaction_meta_fdr", bn_lc, fixed = TRUE)

    pal_cell <- NULL
    if (is_interaction) {
      pal0 <- .get_interaction_palette()
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
      message("[meta-heatmap] failed to build:", csv)
      next
    }

    .meta_save_near_v3_shortname(h$plot, h$width, h$height, csv, meta,
      suffix = "",
      dpi = getOption("csn_meta_heatmap_dpi", 600)
    )


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
            d50 <- .meta_read_summary_long(cont_csv) %>% dplyr::filter(.data$predictor == "Z_CSN_SCORE")
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
            .meta_save_near_v3_shortname(h2$plot, h2$width, h2$height, csv, meta,
              suffix = "_ordered",
              dpi = getOption("csn_meta_heatmap_dpi", 600)
            )
          }
        }
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
            .coll_tok2 <- sub("^Summary_(.+?)_GSEA.*$", "\\1", basename(all_csv))
            if (!identical(toupper(.coll_tok2), "H")) {
              .top_n <- get0("PAN_HEATMAP_TOP_N", ifnotfound = 25)
              .bot_n <- get0("PAN_HEATMAP_BOTTOM_N", ifnotfound = 25)
              d50 <- .meta_read_summary_long(all_csv) %>% dplyr::filter(.data$predictor == "Z_CSN_SCORE")
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
              .meta_save_near_v3_shortname(h3$plot, h3$width, h3$height, csv, meta,
                suffix = "_ordered",
                dpi = getOption("csn_meta_heatmap_dpi", 600)
              )
            }
          }
        }
      }
    }

    try(closeAllConnections(), silent = TRUE)
    invisible(gc(FALSE))
    message("[meta-heatmap] Completed:", csv)
  }
  invisible(TRUE)
}


options(HEATMAP_PATHWAY_PREFIX = c("KINASE-PSP", "PATH-NP"))
## options(HEATMAP_PATHWAY_PREFIX = NULL)

options(csn.dot_size_meta = 2.0)

## ---- Execution (Example) ----
run_meta_fdr_heatmaps(
  root = if (is.null(COMBO_PREFIX)) {
    "PTMsigDB/meta_fdr/summary"
  } else {
    file.path(COMBO_PREFIX, #"PTMsigDB", 
              "meta_fdr", "summary")
  },
  width_in = 18, height_in = 10, dpi = 600
)
