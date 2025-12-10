# Path to built-in MAGMA gene-loc files
#
# @param species One of "maize", "sorghum", "arabidopsis".
# @return Path to the built-in gene-loc file.
# @keywords internal
builtin_geneloc_path <- function(species) {
  species <- match.arg(
    tolower(species),
    c("maize", "sorghum", "arabidopsis")
  )

  fname <- switch(
    species,
    maize       = "maize.genes.loc",
    sorghum     = "sorghum.genes.loc",
    arabidopsis = "arabidopsis.genes.loc"
  )

  path <- system.file("extdata", fname, package = "MAGCAT")

  if (path == "") {
    stop("Built-in gene-loc file not found for species = '", species,
         "'. Did you put ", fname, " into inst/extdata/?",
         call. = FALSE)
  }

  path
}


#' Run MAGMA annotation starting from GWAS summary stats
#'
#' This function takes a GWAS summary statistics file, builds a SNP location
#' file (SNP, CHR, BP) internally, and runs MAGMA annotation.
#'
#' @param stats_file GWAS summary statistics file containing at least
#'   chromosome, SNP ID, and base-pair position columns.
#' @param rename_columns Named character vector mapping standard names to
#'   column names in `stats_file`, e.g.
#'   `c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald")`.
#'    NMISS = "n_miss")` or `c(..., NOBS = "Nobs")`
#'   For annotation, only `CHR`, `SNP`, and `POS` are used.
#' @param gene_loc Path to gene location file (MAGMA gene-loc format),
#'   or NULL if using `species`.
#' @param out_prefix Output prefix for MAGMA (file name, without directory).
#' @param out_dir Optional output directory. If not NULL, the final prefix
#'   will be `file.path(out_dir, out_prefix)`. The directory will be created
#'   if it does not exist.
#' @param window Optional length-2 numeric vector (kb upstream, kb downstream).
#'   For example, `window = c(25, 25)` gives "window=25,25".
#' @param species Optional; one of "maize", "sorghum", "arabidopsis"
#'   to use a built-in gene-loc file stored in the package.
#'
#' @export
magma_annotate <- function(stats_file,
                           rename_columns,
                           gene_loc   = NULL,
                           out_prefix,
                           out_dir    = NULL,
                           window     = NULL,
                           species    = NULL) {
  mp <- magma_path()

  if (!file.exists(stats_file)) {
    stop("stats_file does not exist: ", stats_file, call. = FALSE)
  }

  # ---- handle output directory + final out prefix ----
  base_prefix <- out_prefix
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    prefix_full <- file.path(out_dir, base_prefix)
  } else {
    prefix_full <- base_prefix
  }

  # ---- read stats and build SNP-loc file ----
  dat <- utils::read.table(
    stats_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # determine which columns are CHR, SNP, POS
  actual_chr <- unname(rename_columns["CHR"])
  actual_snp <- unname(rename_columns["SNP"])
  actual_pos <- unname(rename_columns["POS"])

  if (any(is.na(c(actual_chr, actual_snp, actual_pos)))) {
    stop("rename_columns must at least contain 'CHR', 'SNP', and 'POS', e.g.\n",
         "  c(CHR='chr', SNP='rs', POS='ps', PVALUE='p_wald')",
         call. = FALSE)
  }

  needed  <- c(actual_chr, actual_snp, actual_pos)
  missing <- setdiff(needed, colnames(dat))
  if (length(missing) > 0L) {
    stop("Missing required columns in stats_file for SNP-loc: ",
         paste(missing, collapse = ", "),
         call. = FALSE)
  }

  # normalise CHR: strip leading "chr" if present
  chr_vec <- as.character(dat[[actual_chr]])
  chr_vec <- sub("^chr", "", chr_vec, ignore.case = TRUE)

  snp_df <- data.frame(
    SNP = dat[[actual_snp]],
    CHR = chr_vec,
    BP  = dat[[actual_pos]],
    stringsAsFactors = FALSE
  )

  # write SNP-loc file for MAGMA
  if (!is.null(out_dir)) {
    snp_loc <- file.path(out_dir, paste0(base_prefix, ".snp.loc.txt"))
  } else {
    snp_loc <- file.path(tempdir(), paste0(base_prefix, ".snp.loc.txt"))
  }

  utils::write.table(
    snp_df,
    file      = snp_loc,
    quote     = FALSE,
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # ---- decide which gene_loc to use (species vs direct) ----
  if (!is.null(species) && !is.null(gene_loc)) {
    stop("Provide either 'gene_loc' OR 'species', not both.", call. = FALSE)
  }

  if (!is.null(species)) {
    gene_loc <- builtin_geneloc_path(species)
  }

  if (is.null(gene_loc)) {
    stop("You must provide either:\n",
         "  * gene_loc = 'path/to/file.genes.loc', OR\n",
         "  * species = 'maize' | 'sorghum' | 'arabidopsis'.",
         call. = FALSE)
  }

  # ---- build MAGMA args and run ----
  args <- c("--annotate")

  if (!is.null(window)) {
    if (length(window) != 2L) {
      stop("window must be NULL or a length-2 numeric vector, e.g. c(25, 25).")
    }
    args <- c(args, paste0("window=", window[1], ",", window[2]))
  }

  args <- c(
    args,
    "--snp-loc",  snp_loc,
    "--gene-loc", gene_loc,
    "--out",      prefix_full
  )

  system2(mp, args)
}



#' Run MAGMA gene analysis on SNP p-values
#'
#' You can either:
#' * supply a MAGMA-ready p-value file via `pval_file`, with columns
#'   specified by `snp_col` and `p_col`, OR
#' * supply a GWAS summary statistics file via `stats_file`, in which case
#'   a MAGMA p-value file is created internally using the mapping given
#'   by `rename_columns`.
#'
#' Per-SNP N handling:
#' * If `rename_columns` contains `N`, that column in `stats_file` is used
#'   as per-SNP N via `ncol=N` in MAGMA.
#' * Else if `rename_columns` contains `NOBS`, that column in `stats_file`
#'   is used as per-SNP N via `ncol=N`.
#' * Else if `rename_columns` contains `NMISS` and `n_total` is provided,
#'   per-SNP N is computed as `N = n_total - NMISS` and used via `ncol=N`.
#' * Otherwise, a constant sample size `n_total` is required and passed as
#'   `N=n_total`.
#'
#' Parallelisation by chromosome (stats_file mode only):
#' * If `chroms` is non-NULL and `n_threads > 1`, MAGMA is run separately
#'   for each chromosome in `chroms`, in parallel, using `n_threads`
#'   workers. Each run gets a suffix `_chr<k>` in `out_prefix`.
#' * This requires `stats_file` + `rename_columns['CHR']`.
#'
#' @param bfile PLINK prefix for LD reference.
#' @param gene_annot MAGMA gene annotation file (.genes.annot).
#' @param pval_file Optional; SNP p-value file for MAGMA's --pval.
#'   If NULL, `stats_file` must be provided.
#' @param stats_file Optional; GWAS summary statistics file containing
#'   at least SNP ID and p-value columns. If provided, `pval_file`
#'   is ignored and a MAGMA p-value file is generated internally.
#' @param n_total Total sample size (N=... in MAGMA) when using a
#'   constant N or computing per-SNP N from NMISS. Not required if
#'   `rename_columns` supplies an `N` or `NOBS` column.
#' @param snp_col Column name containing SNP IDs (used only when
#'   `pval_file` is given and `rename_columns` is NULL).
#' @param p_col Column name containing p-values (used only when
#'   `pval_file` is given and `rename_columns` is NULL).
#' @param rename_columns Optional named character vector mapping standard
#'   names to the column names in `stats_file`, e.g.
#'   `c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald",
#'      NMISS = "n_miss")` or `c(..., N = "N_eff")` or `c(..., NOBS = "Nobs")`.
#'   For `magma_gene()`, only `CHR`, `SNP`, `PVALUE`, `N`, `NOBS`, and `NMISS` are used.
#' @param out_prefix Base output prefix (file name without directory).
#' @param out_dir Optional output directory. If not NULL, the final prefix
#'   will be file.path(out_dir, out_prefix). The directory is created
#'   if it does not exist. For `stats_file` mode, the internal p-value
#'   file is also written into this directory.
#' @param gene_model Character vector of one or more MAGMA gene models,
#'   e.g. "snp-wise=top", "snp-wise=mean", "multi".
#' @param chroms Optional vector of chromosomes to analyse separately
#'   (e.g. 1:10). Only used in `stats_file` mode. If non-NULL and
#'   `n_threads > 1`, analysis is run in parallel, one job per chromosome.
#' @param n_threads Number of parallel workers to use when `chroms`
#'   is provided. Defaults to 1 (no parallelism). Will be clamped to
#'   `min(length(chroms), detectCores() - 1)`.
#'
#' @export
magma_gene <- function(bfile,
                       gene_annot,
                       pval_file  = NULL,
                       stats_file = NULL,
                       n_total,
                       snp_col   = "SNP",
                       p_col     = "P",
                       rename_columns = NULL,
                       out_prefix,
                       out_dir    = NULL,
                       gene_model = "snp-wise=top",
                       chroms     = NULL,
                       n_threads  = 1,
                       chr_keep   = NULL) {  # chr_keep = internal filter for per-chr mode
  mp <- magma_path()

  ## ---------- PARALLEL PER-CHR WRAPPER (stats_file mode only) ----------
  # If user supplied chroms + stats_file and wants >1 thread, run one job per chr
  if (!is.null(chroms) && length(chroms) > 0 &&
      !is.null(stats_file) && is.null(pval_file) &&
      (is.null(chr_keep) || length(chr_keep) == 0) &&
      n_threads > 1) {

    # work out how many threads to actually use
    cores <- parallel::detectCores()
    n_threads_eff <- max(1, min(length(chroms), cores - 1, n_threads))

    cl <- parallel::makeCluster(n_threads_eff)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # export this function + needed helpers to workers
    parallel::clusterExport(
      cl,
      varlist = c("magma_gene", "magma_path", "builtin_geneloc_path"),
      envir   = environment()
    )

    # run one magma_gene per chr, sequential inside worker, but workers in parallel
    res <- parallel::parLapply(
      cl,
      X = chroms,
      fun = function(chr) {
        magma_gene(
          bfile          = bfile,
          gene_annot     = gene_annot,
          pval_file      = NULL,
          stats_file     = stats_file,
          n_total        = n_total,
          snp_col        = snp_col,
          p_col          = p_col,
          rename_columns = rename_columns,
          out_prefix     = paste0(out_prefix, "_chr", chr),
          out_dir        = out_dir,
          gene_model     = gene_model,
          chroms         = NULL,   # avoid recursion
          n_threads      = 1,
          chr_keep       = chr     # internal filter
        )
      }
    )

    return(invisible(res))
  }
  ## --------------------------------------------------------------------

  ## ---- decide which p-value source to use ----
  if (!is.null(pval_file) && !is.null(stats_file)) {
    stop("Provide either 'pval_file' OR 'stats_file', not both.", call. = FALSE)
  }
  if (is.null(pval_file) && is.null(stats_file)) {
    stop("You must provide either 'pval_file' (MAGMA-ready) or 'stats_file' (GWAS summary).",
         call. = FALSE)
  }

  # handle output directory + full prefix for MAGMA outputs
  base_prefix <- out_prefix
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    prefix_full <- file.path(out_dir, base_prefix)
  } else {
    prefix_full <- base_prefix
  }

  ## ---- build pval file + decide N-argument (N= vs ncol=) ----
  if (!is.null(stats_file)) {
    if (!file.exists(stats_file)) {
      stop("stats_file does not exist: ", stats_file, call. = FALSE)
    }

    dat <- utils::read.table(
      stats_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    # Which columns are SNP, PVALUE (and optionally CHR, N, NOBS, NMISS)?
    if (!is.null(rename_columns)) {
      actual_snp <- unname(rename_columns["SNP"])
      actual_p   <- unname(rename_columns["PVALUE"])

      if (is.na(actual_snp) || is.na(actual_p)) {
        stop("rename_columns must at least contain 'SNP' and 'PVALUE', e.g.\n",
             "  c(CHR='chr', SNP='rs', POS='ps', PVALUE='p_wald')",
             call. = FALSE)
      }

      # optional CHR for per-chr filtering
      actual_chr <- if ("CHR" %in% names(rename_columns)) {
        unname(rename_columns["CHR"])
      } else {
        NA_character_
      }

      # optional N, NOBS and NMISS
      n_col_name <- if ("N" %in% names(rename_columns)) {
        unname(rename_columns["N"])
      } else {
        NA_character_
      }

      nobs_col_name <- if ("NOBS" %in% names(rename_columns)) {
        unname(rename_columns["NOBS"])
      } else {
        NA_character_
      }

      nmiss_col_raw <- if ("NMISS" %in% names(rename_columns)) {
        unname(rename_columns["NMISS"])
      } else {
        NA_character_
      }

    } else {
      actual_snp     <- snp_col
      actual_p       <- p_col
      actual_chr     <- NA_character_
      n_col_name     <- NA_character_
      nobs_col_name  <- NA_character_
      nmiss_col_raw  <- NA_character_
    }

    # basic SNP / P checks
    needed_basic <- c(actual_snp, actual_p)
    missing_basic <- setdiff(needed_basic, colnames(dat))
    if (length(missing_basic) > 0L) {
      stop("Missing required columns in stats_file: ",
           paste(missing_basic, collapse = ", "),
           call. = FALSE)
    }

    # optional per-chr filtering (internal, used by parallel wrapper)
    if (!is.null(chr_keep)) {
      if (is.na(actual_chr)) {
        stop("Per-chromosome filtering (chr_keep) requires rename_columns['CHR'] to be set.",
             call. = FALSE)
      }
      if (!actual_chr %in% colnames(dat)) {
        stop("rename_columns['CHR'] = '", actual_chr,
             "' but that column is not found in stats_file.", call. = FALSE)
      }

      chr_vec_raw  <- as.character(dat[[actual_chr]])
      chr_norm     <- sub("^chr", "", chr_vec_raw, ignore.case = TRUE)
      keep_chr_vec <- as.character(chr_keep)
      keep_idx     <- chr_norm %in% keep_chr_vec

      dat <- dat[keep_idx, , drop = FALSE]
      if (!nrow(dat)) {
        stop("No SNPs left after filtering to chromosomes: ",
             paste(keep_chr_vec, collapse = ", "),
             call. = FALSE)
      }
    }

    ## ---- Per-SNP N logic: N / NOBS / NMISS / constant ----

    # Case 1: per-SNP N column already present (N or NOBS)
    if (!is.na(n_col_name) || !is.na(nobs_col_name)) {

      # prefer explicit N over NOBS if both are present
      n_col_use <- if (!is.na(n_col_name)) n_col_name else nobs_col_name

      if (!n_col_use %in% colnames(dat)) {
        stop("Per-SNP N column '", n_col_use,
             "' (from rename_columns['N' or 'NOBS']) not found in stats_file.",
             call. = FALSE)
      }

      pval_df <- dat[, c(actual_snp, actual_p, n_col_use)]
      colnames(pval_df) <- c("SNP", "P", "N")
      N_arg   <- "ncol=N"
      use_snp <- "SNP"
      use_p   <- "P"

    } else if (!is.na(nmiss_col_raw)) {
      # Case 2: have NMISS, compute per-SNP N = n_total - NMISS
      if (missing(n_total) || is.null(n_total)) {
        stop("n_total must be provided when using NMISS to compute per-SNP N.",
             call. = FALSE)
      }
      if (!nmiss_col_raw %in% colnames(dat)) {
        stop("rename_columns['NMISS'] = '", nmiss_col_raw,
             "' but that column is not found in stats_file.", call. = FALSE)
      }

      nmiss <- dat[[nmiss_col_raw]]
      N_eff <- n_total - nmiss
      if (any(N_eff <= 0, na.rm = TRUE)) {
        stop("Computed per-SNP N (n_total - NMISS) has non-positive values.",
             call. = FALSE)
      }

      pval_df <- data.frame(
        SNP = dat[[actual_snp]],
        P   = dat[[actual_p]],
        N   = N_eff,
        stringsAsFactors = FALSE
      )
      N_arg   <- "ncol=N"
      use_snp <- "SNP"
      use_p   <- "P"

    } else {
      # Case 3: no per-SNP info, use constant n_total
      if (missing(n_total) || is.null(n_total)) {
        stop("n_total must be provided when no per-SNP N (N, NOBS or NMISS) is given.",
             call. = FALSE)
      }

      pval_df <- dat[, c(actual_snp, actual_p)]
      colnames(pval_df) <- c("SNP", "P")
      N_arg   <- paste0("N=", n_total)
      use_snp <- "SNP"
      use_p   <- "P"
    }

    # write pval file
    if (!is.null(out_dir)) {
      pval_path <- file.path(out_dir, paste0(base_prefix, ".pval.txt"))
    } else {
      pval_path <- file.path(tempdir(), paste0(base_prefix, ".pval.txt"))
    }

    utils::write.table(
      pval_df,
      file      = pval_path,
      quote     = FALSE,
      sep       = "\t",
      row.names = FALSE,
      col.names = TRUE
    )

    pval_file_for_magma <- pval_path

  } else {
    # pval_file mode: use MAGMA-ready file, constant N only (no per-CHR support)
    if (!file.exists(pval_file)) {
      stop("pval_file does not exist: ", pval_file, call. = FALSE)
    }
    if (!is.null(chr_keep)) {
      stop("chr_keep / chroms parallel mode is only supported in stats_file mode, ",
           "because pval_file may not contain chromosome information.",
           call. = FALSE)
    }
    if (missing(n_total) || is.null(n_total)) {
      stop("n_total must be provided when using pval_file mode.",
           call. = FALSE)
    }

    pval_file_for_magma <- pval_file
    N_arg   <- paste0("N=", n_total)
    use_snp <- snp_col
    use_p   <- p_col
  }
  ## -------------------------------------------

  # ensure gene_model is a vector
  gene_model <- as.character(gene_model)

  out_files <- vector("list", length(gene_model))

  for (i in seq_along(gene_model)) {
    gm <- gene_model[i]

    # nice suffix for output name (remove weird chars)
    gm_suffix <- gsub("[^A-Za-z0-9]+", "_", gm)
    out_i <- paste0(prefix_full, ".", gm_suffix)

    # MAGMA wants:
    #   --bfile <bfile> --gene-annot <annot> --pval <file> N=... or ncol=...
    #   use=SNP,P --gene-model <model> --out <prefix>
    args <- c(
      "--bfile",      bfile,
      "--gene-annot", gene_annot,
      "--pval",       pval_file_for_magma,
      N_arg,
      paste0("use=", use_snp, ",", use_p),
      "--gene-model", gm,
      "--out",        out_i
    )

    system2(mp, args)

    out_files[[i]] <- list(
      gene_model = gm,
      out_prefix = out_i,
      pval_file  = pval_file_for_magma
    )
  }

  invisible(out_files)
}
