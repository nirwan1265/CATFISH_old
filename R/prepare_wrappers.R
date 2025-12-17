## ============================================================
##  Shared "freeze" core used by Fisher / minP / softTFisher / Stouffer
## ============================================================

#' Prepare common pathway + gene maps for fast reuse
#'
#' This freezes the expensive, repeated setup steps you currently redo
#' inside every wrapper call: pathway list building, gene normalization,
#' named p-vector creation, canonical gene-id map, and permutation pools.
#'
#' Requires your existing:
#'   - magcat_load_pathways()
#'   - fix_p_for_acat()
#'
#' @keywords internal
magcat_prepare_core <- function(gene_results,
                                pathways     = NULL,
                                species      = NULL,
                                pmn_gene_col = NULL,
                                gene_col     = "GENE",
                                p_col        = "P",
                                weight_col   = NULL) {

  ## pathway source
  if (!is.null(pathways) && !is.null(species)) {
    stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
  }
  if (is.null(pathways) && is.null(species)) {
    stop("You must provide either 'pathways' or 'species'.", call. = FALSE)
  }
  if (is.null(pathways) && !is.null(species)) {
    if (!exists("magcat_load_pathways", mode = "function")) {
      stop("magcat_prepare_core(): missing magcat_load_pathways().", call. = FALSE)
    }
    pathways <- magcat_load_pathways(
      species  = species,
      gene_col = pmn_gene_col
    )
  }

  ## gene_results checks
  if (!all(c(gene_col, p_col) %in% names(gene_results))) {
    stop("gene_results must contain columns '", gene_col, "' and '", p_col, "'.",
         call. = FALSE)
  }

  gr         <- gene_results
  genes_all  <- as.character(gr[[gene_col]])
  p_all      <- as.numeric(gr[[p_col]])
  genes_norm <- tolower(genes_all)

  ## named vectors
  gene_p_vec <- stats::setNames(p_all, genes_norm)
  gene_map   <- tapply(genes_all, genes_norm, function(x) x[1])

  ## optional weights (for Stouffer)
  w_norm <- NULL
  if (!is.null(weight_col)) {
    if (!weight_col %in% names(gr)) {
      stop("weight_col = '", weight_col, "' not found in gene_results.", call. = FALSE)
    }
    w_all  <- as.numeric(gr[[weight_col]])
    w_norm <- stats::setNames(w_all, genes_norm)
  }

  ## pathways -> named list + names
  if (is.data.frame(pathways)) {
    if (!all(c("pathway_id", "gene_id") %in% names(pathways))) {
      stop("If 'pathways' is a data.frame it must have columns pathway_id and gene_id.",
           call. = FALSE)
    }
    if (!"pathway_name" %in% names(pathways)) {
      pathways$pathway_name <- pathways$pathway_id
    }
    p_list  <- split(pathways$gene_id, pathways$pathway_id)
    p_names <- tapply(pathways$pathway_name, pathways$pathway_id, function(x) x[1])
  } else if (is.list(pathways)) {
    p_list  <- pathways
    p_names <- names(pathways)
    if (is.null(p_names)) {
      p_names <- paste0("PWY_", seq_along(pathways))
      names(p_list) <- p_names
    }
  } else {
    stop("'pathways' must be either a list or a data.frame.", call. = FALSE)
  }

  ## normalize IDs inside pathways
  p_list <- lapply(p_list, function(g) tolower(as.character(g)))

  ## permutation pools (store both styles)
  idx_pool_na     <- which(!is.na(p_all))
  idx_pool_valid  <- which(is.finite(p_all) & p_all > 0 & p_all < 1)

  if (!length(idx_pool_na)) {
    stop("No non-NA gene-level p-values found in gene_results.", call. = FALSE)
  }
  if (!length(idx_pool_valid)) {
    ## some wrappers use valid01 pool (e.g. Stouffer); keep NA pool anyway
    idx_pool_valid <- idx_pool_na
  }

  structure(
    list(
      gene_col     = gene_col,
      p_col        = p_col,
      weight_col   = weight_col,
      genes_all    = genes_all,
      genes_norm   = genes_norm,
      p_all        = p_all,
      gene_p_vec   = gene_p_vec,
      gene_map     = gene_map,
      w_norm       = w_norm,
      p_list       = p_list,
      p_names      = p_names,
      idx_pool_na  = idx_pool_na,
      idx_pool_01  = idx_pool_valid
    ),
    class = "magcat_core_prep"
  )
}


## ============================================================
##  Fisher: prepare + run_prepared
## ============================================================

#' Prepare Fisher pathway test objects
#' @export
magcat_fisher_prepare <- function(gene_results,
                                  pathways     = NULL,
                                  species      = NULL,
                                  pmn_gene_col = NULL,
                                  gene_col     = "GENE",
                                  p_col        = "P",
                                  min_p        = 1e-15,
                                  do_fix       = TRUE) {

  core <- magcat_prepare_core(
    gene_results  = gene_results,
    pathways      = pathways,
    species       = species,
    pmn_gene_col  = pmn_gene_col,
    gene_col      = gene_col,
    p_col         = p_col,
    weight_col    = NULL
  )

  structure(
    c(core, list(min_p = min_p, do_fix = do_fix)),
    class = c("magcat_fisher_prep", class(core))
  )
}

#' Run Fisher pathway test from a prepared object
#' @export
magcat_fisher_run_prepared <- function(prep,
                                       output  = FALSE,
                                       out_dir = "magcat_fisher") {

  if (!inherits(prep, "magcat_fisher_prep")) {
    stop("magcat_fisher_run_prepared(): 'prep' must come from magcat_fisher_prepare().",
         call. = FALSE)
  }
  if (!exists("fix_p_for_acat", mode = "function")) {
    stop("magcat_fisher_run_prepared(): missing fix_p_for_acat().", call. = FALSE)
  }

  p_list  <- prep$p_list
  p_names <- prep$p_names
  gene_p  <- prep$gene_p_vec
  gmap    <- prep$gene_map

  n_pw <- length(p_list)
  res  <- data.frame(
    pathway_id   = names(p_list),
    pathway_name = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes      = NA_integer_,
    gene_names   = NA_character_,
    fisher_p     = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_pw)) {
    g_i <- p_list[[i]]
    p_i <- gene_p[g_i]
    keep <- !is.na(p_i)
    p_i  <- as.numeric(p_i[keep])
    g_use <- g_i[keep]

    d <- length(p_i)
    res$n_genes[i] <- d
    if (d == 0L) next

    canon <- unname(gmap[g_use])
    canon <- unique(canon[!is.na(canon)])
    res$gene_names[i] <- paste(canon, collapse = ";")

    if (prep$do_fix) {
      p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)
    }

    stat <- -2 * sum(log(p_i))
    res$fisher_p[i] <- stats::pchisq(stat, df = 2 * length(p_i), lower.tail = FALSE)
  }

  ord <- order(res$fisher_p, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  if (output) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(res, file.path(out_dir, "magcat_fisher_prepared.csv"), row.names = FALSE)
  }
  res
}


## ============================================================
##  minP: prepare + run_prepared (analytic only; perm handled in omni)
## ============================================================

#' Prepare minP pathway test objects
#' @export
magcat_minp_prepare <- function(gene_results,
                                pathways     = NULL,
                                species      = NULL,
                                pmn_gene_col = NULL,
                                gene_col     = "GENE",
                                p_col        = "P",
                                min_p        = 1e-15,
                                do_fix       = TRUE) {

  core <- magcat_prepare_core(
    gene_results  = gene_results,
    pathways      = pathways,
    species       = species,
    pmn_gene_col  = pmn_gene_col,
    gene_col      = gene_col,
    p_col         = p_col,
    weight_col    = NULL
  )

  structure(
    c(core, list(min_p = min_p, do_fix = do_fix)),
    class = c("magcat_minp_prep", class(core))
  )
}

#' Run minP pathway test from a prepared object (analytic)
#' @export
magcat_minp_run_prepared <- function(prep,
                                     output  = FALSE,
                                     out_dir = "magcat_minp") {

  if (!inherits(prep, "magcat_minp_prep")) {
    stop("magcat_minp_run_prepared(): 'prep' must come from magcat_minp_prepare().",
         call. = FALSE)
  }
  if (!exists("fix_p_for_acat", mode = "function")) {
    stop("magcat_minp_run_prepared(): missing fix_p_for_acat().", call. = FALSE)
  }

  .analytic_minp <- function(p_vec) {
    p_vec <- as.numeric(p_vec)
    p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec < 1]
    if (!length(p_vec)) return(NA_real_)
    if (requireNamespace("metap", quietly = TRUE)) {
      return(tryCatch(metap::minimump(p_vec)$p, error = function(e) NA_real_))
    }
    k <- length(p_vec)
    p_min <- min(p_vec)
    1 - (1 - p_min)^k
  }

  p_list  <- prep$p_list
  p_names <- prep$p_names
  gene_p  <- prep$gene_p_vec
  gmap    <- prep$gene_map

  n_pw <- length(p_list)
  res  <- data.frame(
    pathway_id      = names(p_list),
    pathway_name    = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes         = NA_integer_,
    gene_names      = NA_character_,
    minp_stat       = NA_real_,
    minp_p_analytic = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_pw)) {
    g_i <- p_list[[i]]
    p_i <- gene_p[g_i]
    keep <- !is.na(p_i)
    p_i  <- as.numeric(p_i[keep])
    g_use <- g_i[keep]

    d <- length(p_i)
    res$n_genes[i] <- d
    if (d == 0L) next

    canon <- unname(gmap[g_use])
    canon <- unique(canon[!is.na(canon)])
    res$gene_names[i] <- paste(canon, collapse = ";")

    if (prep$do_fix) {
      p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)
    }

    res$minp_stat[i]       <- min(p_i)
    res$minp_p_analytic[i] <- .analytic_minp(p_i)
  }

  ord <- order(res$minp_p_analytic, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  if (output) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(res, file.path(out_dir, "magcat_minp_prepared.csv"), row.names = FALSE)
  }
  res
}


## ============================================================
##  soft TFisher: prepare + run_prepared (analytic only; perm handled in omni)
## ============================================================

#' Prepare soft-TFisher pathway test objects
#' @export
magcat_soft_tfisher_prepare <- function(gene_results,
                                        pathways     = NULL,
                                        species      = NULL,
                                        pmn_gene_col = NULL,
                                        gene_col     = "GENE",
                                        p_col        = "P",
                                        tau1         = 0.05,
                                        min_p        = 1e-15,
                                        do_fix       = TRUE) {

  if (!requireNamespace("TFisher", quietly = TRUE)) {
    stop("magcat_soft_tfisher_prepare(): requires TFisher package.", call. = FALSE)
  }

  core <- magcat_prepare_core(
    gene_results  = gene_results,
    pathways      = pathways,
    species       = species,
    pmn_gene_col  = pmn_gene_col,
    gene_col      = gene_col,
    p_col         = p_col,
    weight_col    = NULL
  )

  structure(
    c(core, list(tau1 = tau1, min_p = min_p, do_fix = do_fix)),
    class = c("magcat_soft_tfisher_prep", class(core))
  )
}

#' Run soft-TFisher pathway test from a prepared object (analytic)
#' @export
magcat_soft_tfisher_run_prepared <- function(prep,
                                             output  = FALSE,
                                             out_dir = "magcat_tfisher_soft") {

  if (!inherits(prep, "magcat_soft_tfisher_prep")) {
    stop("magcat_soft_tfisher_run_prepared(): 'prep' must come from magcat_soft_tfisher_prepare().",
         call. = FALSE)
  }
  if (!exists("fix_p_for_acat", mode = "function")) {
    stop("magcat_soft_tfisher_run_prepared(): missing fix_p_for_acat().", call. = FALSE)
  }

  p_list  <- prep$p_list
  p_names <- prep$p_names
  gene_p  <- prep$gene_p_vec
  gmap    <- prep$gene_map

  n_pw <- length(p_list)
  res  <- data.frame(
    pathway_id         = names(p_list),
    pathway_name       = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes            = NA_integer_,
    gene_names         = NA_character_,
    tfisher_stat       = NA_real_,
    tfisher_p_analytic = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_pw)) {
    g_i <- p_list[[i]]
    p_i <- gene_p[g_i]
    keep <- !is.na(p_i)
    p_i  <- as.numeric(p_i[keep])
    g_use <- g_i[keep]

    d <- length(p_i)
    res$n_genes[i] <- d
    if (d == 0L) next

    canon <- unname(gmap[g_use])
    canon <- unique(canon[!is.na(canon)])
    res$gene_names[i] <- paste(canon, collapse = ";")

    if (prep$do_fix) {
      p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)
    }

    stat <- TFisher::stat.soft(p = p_i, tau1 = prep$tau1)
    res$tfisher_stat[i] <- stat

    Fq <- TFisher::p.soft(q = stat, n = length(p_i), tau1 = prep$tau1, M = NULL)
    res$tfisher_p_analytic[i] <- 1 - as.numeric(Fq)
  }

  ord <- order(res$tfisher_p_analytic, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  if (output) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(res, file.path(out_dir, "magcat_tfisher_soft_prepared.csv"), row.names = FALSE)
  }
  res
}


## ============================================================
##  Stouffer: prepare + run_prepared (analytic only; perm handled in omni)
## ============================================================

#' Prepare Stouffer pathway test objects
#' @export
magcat_stouffer_prepare <- function(gene_results,
                                    pathways     = NULL,
                                    species      = NULL,
                                    pmn_gene_col = NULL,
                                    gene_col     = "GENE",
                                    p_col        = "P",
                                    weight_col   = NULL,
                                    min_p        = 1e-15,
                                    do_fix       = TRUE) {

  if (!requireNamespace("metap", quietly = TRUE)) {
    stop("magcat_stouffer_prepare(): requires metap package.", call. = FALSE)
  }

  core <- magcat_prepare_core(
    gene_results  = gene_results,
    pathways      = pathways,
    species       = species,
    pmn_gene_col  = pmn_gene_col,
    gene_col      = gene_col,
    p_col         = p_col,
    weight_col    = weight_col
  )

  structure(
    c(core, list(min_p = min_p, do_fix = do_fix)),
    class = c("magcat_stouffer_prep", class(core))
  )
}

#' Run Stouffer pathway test from a prepared object (analytic)
#' @export
magcat_stouffer_run_prepared <- function(prep,
                                         output  = FALSE,
                                         out_dir = "magcat_stouffer") {

  if (!inherits(prep, "magcat_stouffer_prep")) {
    stop("magcat_stouffer_run_prepared(): 'prep' must come from magcat_stouffer_prepare().",
         call. = FALSE)
  }
  if (!exists("fix_p_for_acat", mode = "function")) {
    stop("magcat_stouffer_run_prepared(): missing fix_p_for_acat().", call. = FALSE)
  }

  p_list  <- prep$p_list
  p_names <- prep$p_names
  gene_p  <- prep$gene_p_vec
  gmap    <- prep$gene_map
  w_norm  <- prep$w_norm

  n_pw <- length(p_list)
  res  <- data.frame(
    pathway_id          = names(p_list),
    pathway_name        = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes             = NA_integer_,
    gene_names          = NA_character_,
    stouffer_z          = NA_real_,
    stouffer_p_analytic = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_pw)) {
    g_i <- p_list[[i]]
    p_i <- gene_p[g_i]
    keep <- !is.na(p_i)
    p_i  <- as.numeric(p_i[keep])
    g_use <- g_i[keep]

    d <- length(p_i)
    res$n_genes[i] <- d
    if (d < 2L) next

    canon <- unname(gmap[g_use])
    canon <- unique(canon[!is.na(canon)])
    res$gene_names[i] <- paste(canon, collapse = ";")

    if (prep$do_fix) {
      p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)
    }

    ## weights
    if (!is.null(w_norm)) {
      w_i <- as.numeric(w_norm[g_use])
      if (any(is.na(w_i) | w_i <= 0)) {
        w_pos <- w_i[!is.na(w_i) & w_i > 0]
        repl  <- if (length(w_pos)) stats::median(w_pos) else 1
        w_i[is.na(w_i) | w_i <= 0] <- repl
      }
    } else {
      w_i <- rep(1, length(p_i))
    }

    sz <- tryCatch(metap::sumz(p = p_i, weights = w_i), error = function(e) NULL)
    if (!is.null(sz)) {
      res$stouffer_z[i]          <- as.numeric(sz$z)
      res$stouffer_p_analytic[i] <- as.numeric(sz$p)
    }
  }

  ord <- order(res$stouffer_p_analytic, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  if (output) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(res, file.path(out_dir, "magcat_stouffer_prepared.csv"), row.names = FALSE)
  }
  res
}


## ============================================================
##  ACAT
## ============================================================

#' Run ACAT pathway test from a prepared object
#' @export
magcat_acat_run_prepared <- function(prep,
                                     output  = FALSE,
                                     out_dir = "magcat_acat") {

  if (is.null(prep$idx_list) || is.null(prep$p_all_u) || is.null(prep$pathway_id)) {
    stop("magcat_acat_run_prepared(): 'prep' does not look like output of magcat_acat_prepare().",
         call. = FALSE)
  }
  if (!requireNamespace("ACAT", quietly = TRUE)) {
    stop("magcat_acat_run_prepared(): requires ACAT package.", call. = FALSE)
  }

  idx_list <- prep$idx_list
  p_all_u  <- as.numeric(prep$p_all_u)

  # optional: original-case gene IDs (if you stored them)
  genes_all_u <- prep$genes_all_u
  if (is.null(genes_all_u)) {
    # fall back: use normalized ids
    genes_all_u <- prep$genes_norm_u
  }

  n_pw <- length(idx_list)

  res <- data.frame(
    pathway_id   = prep$pathway_id,
    pathway_name = prep$pathway_name,
    n_genes      = as.integer(prep$n_genes),
    gene_names   = NA_character_,
    acat_p       = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_pw)) {
    idx <- idx_list[[i]]
    d   <- length(idx)
    if (d < 1L) next

    p_i <- p_all_u[idx]
    # (p_all_u is already fixed in prepare if do_fix=TRUE, but safe anyway)
    p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)

    res$acat_p[i] <- ACAT::ACAT(Pvals = p_i)
    res$gene_names[i] <- paste(unique(genes_all_u[idx]), collapse = ";")
  }

  ord <- order(res$acat_p, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  if (output) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(res, file.path(out_dir, "magcat_acat_prepared.csv"), row.names = FALSE)
  }

  res
}


#' Prepare ACAT pathway test objects (fast reusable setup)
#' @export
magcat_acat_prepare <- function(gene_results,
                                pathways     = NULL,
                                species      = NULL,
                                pmn_gene_col = NULL,
                                gene_col     = "GENE",
                                p_col        = "P",
                                min_p        = 1e-15,
                                do_fix       = TRUE) {

  # pathway source
  if (!is.null(pathways) && !is.null(species)) {
    stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
  }
  if (is.null(pathways) && is.null(species)) {
    stop("You must provide either 'pathways' OR 'species'.", call. = FALSE)
  }
  if (is.null(pathways) && !is.null(species)) {
    if (!exists("magcat_load_pathways", mode = "function")) {
      stop("magcat_acat_prepare(): missing magcat_load_pathways().", call. = FALSE)
    }
    pathways <- magcat_load_pathways(
      species  = species,
      gene_col = pmn_gene_col
    )
  }

  # gene_results checks
  if (!all(c(gene_col, p_col) %in% names(gene_results))) {
    stop("gene_results must contain columns '", gene_col, "' and '", p_col, "'.", call. = FALSE)
  }
  if (!exists("fix_p_for_acat", mode = "function")) {
    stop("magcat_acat_prepare(): missing fix_p_for_acat().", call. = FALSE)
  }

  genes_all  <- as.character(gene_results[[gene_col]])
  p_all      <- as.numeric(gene_results[[p_col]])
  genes_norm <- tolower(genes_all)

  ok <- which(!is.na(genes_norm) & genes_norm != "" &
              !is.na(p_all) & is.finite(p_all) & p_all > 0 & p_all <= 1)

  genes_all  <- genes_all[ok]
  genes_norm <- genes_norm[ok]
  p_all      <- p_all[ok]

  # dedupe by normalized gene id, keep first (matches your earlier logic)
  keep_first <- !duplicated(genes_norm)
  genes_all_u  <- genes_all[keep_first]
  genes_norm_u <- genes_norm[keep_first]
  p_all_u      <- p_all[keep_first]

  if (do_fix) {
    p_all_u <- fix_p_for_acat(p_all_u, min_p = min_p)
  }

  # pathways -> list
  if (is.data.frame(pathways)) {
    if (!all(c("pathway_id","gene_id") %in% names(pathways))) {
      stop("If 'pathways' is a data.frame it must have columns pathway_id and gene_id.", call. = FALSE)
    }
    p_list <- split(pathways$gene_id, pathways$pathway_id)
    p_names <- if ("pathway_name" %in% names(pathways)) {
      tapply(pathways$pathway_name, pathways$pathway_id, function(x) x[1])
    } else {
      setNames(names(p_list), names(p_list))
    }
  } else if (is.list(pathways)) {
    p_list <- pathways
    if (is.null(names(p_list))) names(p_list) <- paste0("PWY_", seq_along(p_list))
    p_names <- setNames(names(p_list), names(p_list))
  } else {
    stop("'pathways' must be list or data.frame.", call. = FALSE)
  }

  # normalize pathway genes + precompute indices into genes_norm_u
  p_list_norm <- lapply(p_list, function(g) tolower(as.character(g)))
  idx_list <- lapply(p_list_norm, function(g) {
    idx <- match(g, genes_norm_u)
    idx <- idx[!is.na(idx)]
    unique(idx)
  })

  n_genes <- vapply(idx_list, length, integer(1))

  structure(
    list(
      genes_all_u  = genes_all_u,
      genes_norm_u = genes_norm_u,
      p_all_u      = p_all_u,
      idx_list     = idx_list,
      pathway_id   = names(idx_list),
      pathway_name = unname(p_names[names(idx_list)]),
      n_genes      = n_genes,
      min_p        = min_p,
      do_fix       = do_fix
    ),
    class = "magcat_acat_prep"
  )
}
