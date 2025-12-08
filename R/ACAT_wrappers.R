## ---------- PMN pathway helpers ------------------------------------

#' Path to built-in PMN pathway file
#'
#' @param species Which built-in DB to use. For now:
#'   \itemize{
#'     \item{"maize"}   {corncyc\_pathways.20230103}
#'     \item{"sorghum"} {sorghumbicolorcyc\_pathways.20230103}
#'     \item{"arabidopsis"} {aracyc\_pathways.20230103}
#'     \item{"plant"}  {plantcyc\_pathways.20230103}
#'   }
#'
#' @keywords internal
magcat_pmn_file <- function(
  species = c("maize", "sorghum", "arabidopsis", "plant")
) {
  species <- match.arg(species)

  fname <- switch(
    species,
    maize       = "corncyc_pathways.20230103",
    sorghum     = "sorghumbicolorcyc_pathways.20230103",
    arabidopsis = "aracyc_pathways.20230103",
    plant       = "plantcyc_pathways.20230103"
  )

  path <- system.file("extdata/pathway", fname, package = "MAGCAT")
  if (path == "") {
    stop("Could not find ", fname, " in inst/extdata/pathway of MAGCAT.",
         call. = FALSE)
  }
  path
}

#' Load PMN pathways as a long data.frame
#'
#' Reads a Plant Metabolic Network (PMN) pathways file from
#' \code{inst/extdata/pathway} and returns a long table of
#' (pathway, gene) membership.
#'
#' @param species One of "maize", "sorghum", "arabidopsis", "plant".
#' @param gene_col Which column from the PMN file to use for gene IDs:
#'   "Gene-id" or "Gene-name".
#' @param drop_unknown If TRUE, drop rows with missing/unknown gene IDs.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{pathway_id}{PMN pathway ID (e.g., "PWY-7634")}
#'     \item{pathway_name}{Human-readable pathway name}
#'     \item{gene_id}{Gene identifier from \code{gene_col}}
#'   }
#' @export
#' Load PMN pathways as a long data.frame
#'
#' Reads a Plant Metabolic Network (PMN) pathways file from
#' `inst/extdata/pathway` and returns a long table of
#' (pathway, gene) membership.
#'
#' @param species One of "maize", "sorghum", "arabidopsis", "plant".
#' @param gene_col Which column from the PMN file to use for gene IDs.
#'   If NULL (default), the function will:
#'   \itemize{
#'     \item use "Gene-name" if present;
#'     \item otherwise fall back to "Gene-id".
#'   }
#' @param drop_unknown If TRUE, drop rows with missing/unknown gene IDs.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{pathway_id}{PMN pathway ID (e.g., "PWY-7634")}
#'     \item{pathway_name}{Human-readable pathway name}
#'     \item{gene_id}{Gene identifier from `gene_col`}
#'   }
#' @export
magcat_load_pathways <- function(
  species      = c("maize", "sorghum", "arabidopsis", "plant"),
  gene_col     = NULL,
  drop_unknown = TRUE
) {
  species <- match.arg(species)

  fpath <- magcat_pmn_file(species)

  x <- utils::read.delim(
    fpath,
    header           = TRUE,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  # 1) If user supplied gene_col, use that (but check it exists)
  # 2) Otherwise, prefer "Gene-name" if present, else "Gene-id"
  if (is.null(gene_col)) {
    candidates <- c("Gene-name", "Gene-id")
    gene_col <- intersect(candidates, names(x))[1]

    if (is.na(gene_col)) {
      stop(
        "Pathway file ", basename(fpath),
        " does not contain any of: ", paste(candidates, collapse = ", "),
        call. = FALSE
      )
    }
  } else {
    if (!gene_col %in% names(x)) {
      stop(
        "Requested gene_col = '", gene_col,
        "' is not a column in ", basename(fpath), ".",
        call. = FALSE
      )
    }
  }

  needed <- c("Pathway-id", "Pathway-name", gene_col)
  if (!all(needed %in% names(x))) {
    stop("Pathway file ", basename(fpath), " does not contain columns: ",
         paste(needed, collapse = ", "),
         call. = FALSE)
  }

  df <- data.frame(
    pathway_id   = x[["Pathway-id"]],
    pathway_name = x[["Pathway-name"]],
    gene_id      = x[[gene_col]],
    stringsAsFactors = FALSE
  )

  if (drop_unknown) {
    bad <- is.na(df$gene_id) |
      df$gene_id == "" |
      df$gene_id == "unknown"
    df <- df[!bad, , drop = FALSE]
  }

  df
}


## ---------- p-value cleaner -----------------------------------------

#' Clean p-values for ACAT
#'
#' @param p numeric vector of p-values (NA allowed).
#' @param min_p lower cap for very small p-values (default 1e-15).
#'
#' @return cleaned numeric vector (NAs removed, extremes capped).
#' @keywords internal
fix_p_for_acat <- function(p, min_p = 1e-15) {
  p <- p[!is.na(p)]
  d <- length(p)
  if (d == 0L) return(p)

  # lower cap
  p[p <= 0] <- min_p

  # upper cap: exact 1's (or >1) become 1 - 1/d
  p[p >= 1] <- 1 - 1 / d

  p
}

## ---------- Pathway-level ACAT wrapper ------------------------------

#' ACAT pathway test on MAGMA gene p-values
#'
#' @param gene_results data.frame with at least a gene column and a p-value column.
#'   Typically a MAGMA `.genes.out` file, where p-values are in column "P".
#' @param pathways either
#'   \itemize{
#'     \item a named list: each element is a character vector of gene IDs
#'       in that pathway, OR
#'     \item a data.frame with columns `pathway_id`, `gene_id`
#'       (and optional `pathway_name`).
#'   }
#' @param species optional; one of "maize", "sorghum", "arabidopsis", "plant".
#'   If provided, built-in PMN pathways are loaded via `magcat_load_pathways()`.
#'   You must provide either `pathways` OR `species`, but not both.
#' @param pmn_gene_col optional; passed to `magcat_load_pathways(gene_col=...)`
#'   when `species` is used. If NULL, PMN loader will prefer "Gene-name"
#'   and fall back to "Gene-id".
#' @param gene_col column name in `gene_results` containing gene IDs
#'   (default "GENE").
#' @param p_col column name in `gene_results` containing gene-level p-values
#'   (default "P").
#' @param min_p lower cap for very small p-values before ACAT (default 1e-15).
#' @param do_fix logical; if TRUE (default), clean p-values with
#'   `fix_p_for_acat()`.
#' @param B integer; number of permutations for empirical p-values.
#'   Default 0 (no permutations).
#' @param seed optional integer random seed for permutations (ignored if B == 0).
#' @param output logical; if TRUE, write a CSV of results to `out_dir`.
#' @param out_dir directory to write CSV when `output = TRUE`
#'   (default "magcat_acat").
#'
#' @return data.frame with columns:
#'   \describe{
#'     \item{pathway_id}{pathway identifier}
#'     \item{pathway_name}{pathway name (if available, otherwise same as id)}
#'     \item{n_genes}{number of genes used in ACAT for this pathway}
#'     \item{gene_names}{semicolon-separated gene IDs used in this pathway}
#'     \item{acat_p}{ACAT p-value}
#'     \item{perm_p}{permutation p-value (NA if B == 0)}
#'   }
#'   Sorted by `acat_p` in decreasing order (largest p first; NA at bottom).
#'   If `output = TRUE`, an attribute `"file"` is attached with the CSV path.
#' @export
magcat_acat_pathways <- function(gene_results,
                                 pathways     = NULL,
                                 species      = NULL,
                                 pmn_gene_col = NULL,
                                 gene_col     = "GENE",
                                 p_col        = "P",
                                 min_p        = 1e-15,
                                 do_fix       = TRUE,
                                 B            = 0L,
                                 seed         = NULL,
                                 output       = FALSE,
                                 out_dir      = "magcat_acat") {
  if (!requireNamespace("ACAT", quietly = TRUE)) {
    stop("Package 'ACAT' is required for magcat_acat_pathways(). ",
         "Please install it with install.packages('ACAT').",
         call. = FALSE)
  }

  ## -------- decide pathway source: pathways vs species ----------
  if (!is.null(pathways) && !is.null(species)) {
    stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
  }
  if (is.null(pathways) && is.null(species)) {
    stop("You must provide either:\n",
         "  * 'pathways' (list/data.frame), OR\n",
         "  * 'species' = 'maize' | 'sorghum' | 'arabidopsis' | 'plant'.",
         call. = FALSE)
  }

  # If user gave species, load PMN pathways
  if (is.null(pathways) && !is.null(species)) {
    if (is.null(pmn_gene_col)) {
      pathways <- magcat_load_pathways(species = species)
    } else {
      pathways <- magcat_load_pathways(
        species  = species,
        gene_col = pmn_gene_col
      )
    }
  }

  ## -------- standardize gene_results ----------
  if (!all(c(gene_col, p_col) %in% names(gene_results))) {
    stop("gene_results must contain columns '", gene_col,
         "' and '", p_col, "'.",
         call. = FALSE)
  }

  gr <- gene_results
  genes_all <- as.character(gr[[gene_col]])
  p_all     <- gr[[p_col]]

  # normalized (lowercase) gene IDs for matching
  genes_all_norm <- tolower(genes_all)
  gene_p_vec     <- stats::setNames(p_all, genes_all_norm)

  # map from normalized -> canonical gene ID (first occurrence)
  gene_map <- tapply(genes_all, genes_all_norm, function(x) x[1])

  ## -------- pathways -> named list + names ----------
  if (is.data.frame(pathways)) {
    if (!"pathway_id" %in% names(pathways) ||
        !"gene_id"    %in% names(pathways)) {
      stop("If 'pathways' is a data.frame it must have columns ",
           "'pathway_id' and 'gene_id'.",
           call. = FALSE)
    }
    if (!"pathway_name" %in% names(pathways)) {
      pathways$pathway_name <- pathways$pathway_id
    }

    p_list  <- split(pathways$gene_id, pathways$pathway_id)
    p_names <- tapply(
      pathways$pathway_name,
      pathways$pathway_id,
      FUN = function(x) x[1]
    )
  } else if (is.list(pathways)) {
    p_list  <- pathways
    p_names <- names(pathways)
    if (is.null(p_names)) {
      p_names <- paste0("PWY_", seq_along(pathways))
      names(p_list) <- p_names
    }
  } else {
    stop("'pathways' must be either a list or a data.frame.",
         call. = FALSE)
  }

  ## normalize gene_ids in pathways to lower-case too
  p_list <- lapply(p_list, function(g) tolower(as.character(g)))

  ## -------- prepare result container ----------
  n_pw <- length(p_list)
  res  <- data.frame(
    pathway_id   = names(p_list),
    pathway_name = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes      = NA_integer_,
    gene_names   = NA_character_,
    acat_p       = NA_real_,
    perm_p       = NA_real_,
    stringsAsFactors = FALSE
  )

  ## -------- permutations set-up ----------
  if (!is.null(seed) && B > 0L) {
    set.seed(seed)
  }

  if (B > 0L) {
    pool_p <- p_all[!is.na(p_all)]
    pool_p <- fix_p_for_acat(pool_p, min_p = min_p)
    n_pool <- length(pool_p)
    if (n_pool == 0L) {
      stop("No non-NA p-values in gene_results for permutations.",
           call. = FALSE)
    }
  }

  ## -------- loop over pathways ----------
  for (i in seq_len(n_pw)) {
    genes_i_norm <- p_list[[i]]
    p_i          <- gene_p_vec[genes_i_norm]

    # drop NAs (genes in pathway that don't appear in MAGMA output)
    keep <- !is.na(p_i)
    p_i  <- p_i[keep]

    d <- length(p_i)
    res$n_genes[i] <- d

    if (d == 0L) {
      res$acat_p[i]     <- NA_real_
      res$perm_p[i]     <- NA_real_
      res$gene_names[i] <- NA_character_
      next
    }

    # canonical gene IDs (original case) for the genes actually used
    used_norm <- names(p_i)
    canon_ids <- unname(gene_map[used_norm])
    canon_ids <- unique(canon_ids[!is.na(canon_ids)])
    res$gene_names[i] <- paste(canon_ids, collapse = ";")

    if (do_fix) {
      p_i <- fix_p_for_acat(p_i, min_p = min_p)
    }

    # ACAT p-value
    res$acat_p[i] <- ACAT::ACAT(Pvals = p_i)

    # permutation-based empirical p-value (optional)
    if (B > 0L) {
      perm_vals <- replicate(
        B,
        {
          p_perm <- sample(pool_p, size = d, replace = FALSE)
          ACAT::ACAT(Pvals = p_perm)
        }
      )

      # empirical p: how many permuted p <= observed
      res$perm_p[i] <- (1 + sum(perm_vals <= res$acat_p[i])) / (B + 1)
    }
  }

  ## -------- sort by acat_p (decreasing) ----------
  ord <- order(res$acat_p, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  ## -------- optional CSV output ----------
  if (output) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    species_tag <- if (is.null(species)) "custom" else species
    out_path <- file.path(
      out_dir,
      paste0("magcat_acat_pathways_", species_tag, ".csv")
    )
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
