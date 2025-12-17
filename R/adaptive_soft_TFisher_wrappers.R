#' Pathway-level *adaptive* soft TFisher (stat.soft / p.soft) on MAGMA gene p-values
#'
#' For each pathway:
#'   * take gene-level p-values (e.g. MAGMA `.genes.out`),
#'   * compute soft-threshold TFisher statistics across a grid of `tau` values,
#'   * compute per-tau analytic p via TFisher::p.soft() under independence
#'     (right-tail p = 1 - F(q)),
#'   * define the adaptive TFisher score as:
#'       p_adapt = min_{tau in tau_grid} p_TFisher(tau),
#'     and record the argmin tau_hat,
#'   * optionally calibrate `p_adapt` with permutations:
#'       - perm_mode="resample_global": resample gene p-values from the global pool
#'       - perm_mode="mvn": simulate correlated null p-values using MAGMA gene correlation
#'
#' NOTE:
#'   - The analytic p-values across tau are *not* multiple-testing corrected;
#'     if you use adaptivity (min over tau), you should use `B_perm > 0` to
#'     obtain a valid calibrated `tfisher_p_perm`.
#'   - All non-NA p-values in the pathway are used (no hard truncation/filtering).
#'
#' @param gene_results data.frame with at least gene + p-value columns.
#'   Typically a MAGMA `.genes.out`-like table, with `gene_col` and `p_col`.
#' @param pathways either:
#'   * a named list: each element is a character vector of gene IDs
#'     in that pathway, OR
#'   * a data.frame with columns `pathway_id`, `gene_id`
#'     (and optional `pathway_name`).
#' @param species optional; one of "maize", "sorghum", "arabidopsis", "plant".
#'   If provided, built-in PMN pathways are loaded via `magcat_load_pathways()`.
#'   You must provide either `pathways` OR `species`, but not both.
#' @param pmn_gene_col optional; passed to `magcat_load_pathways(gene_col=...)`
#'   when `species` is used.
#' @param gene_col column name in `gene_results` containing gene IDs
#'   (default "GENE").
#' @param p_col column name in `gene_results` containing gene-level p-values
#'   (default "P").
#' @param tau_grid numeric vector of tau values to scan (default common grid).
#' @param min_p lower cap for very small p-values (default 1e-15).
#' @param do_fix logical; if TRUE, clean p-values with `fix_p_for_acat()`
#'   (to keep them strictly in (0,1)).
#' @param B_perm integer; number of permutations for empirical p
#'   (default 0L = no permutations).
#' @param perm_mode either "resample_global" or "mvn" (default "resample_global").
#' @param magma_genes_out required if perm_mode="mvn"; path to merged MAGMA *.genes.out
#'   containing at least gene + chromosome (CHR) columns for all genes.
#' @param magma_cor_file optional if perm_mode="mvn"; 3-column file gene1 gene2 r
#'   of gene-gene correlations (if NULL, your MAGMA-based builder may fall back
#'   to block/within-chr heuristics depending on your implementation).
#' @param make_PD logical; if TRUE, force the pathway correlation matrix to be PSD
#'   (default TRUE).
#' @param seed optional RNG seed for permutations (default NULL).
#' @param analytic_logical logical; if TRUE, compute per-tau analytic p via p.soft()
#'   and report the min across tau as `tfisher_p_analytic` (default TRUE).
#' @param output logical; if TRUE, write a CSV of results to `out_dir`.
#' @param out_dir directory to write CSV when `output = TRUE`
#'   (default "magcat_tfisher_soft_adaptive").
#'
#' @return data.frame with columns:
#'   * pathway_id
#'   * pathway_name
#'   * n_genes
#'   * gene_names
#'   * gene_pvals             (semicolon-separated p's used)
#'   * tau_hat                (tau achieving min analytic p)
#'   * tfisher_stat_hat       (stat.soft at tau_hat)
#'   * tfisher_p_analytic     (min over tau of right-tail analytic p; NA if analytic_logical=FALSE)
#'   * tfisher_p_perm         (permutation-calibrated p for the adaptive min-over-tau; NA if B_perm=0)
#'
#' Sorted by `tfisher_p_perm` ascending if `B_perm > 0`, otherwise by
#' `tfisher_p_analytic` ascending.
#' If `output = TRUE`, an attribute `"file"` is attached with the CSV path.
#'
#' @export
magcat_soft_tfisher_adaptive_pathways <- function(gene_results,
                                                  pathways         = NULL,
                                                  species          = NULL,
                                                  pmn_gene_col     = NULL,
                                                  gene_col         = "GENE",
                                                  p_col            = "P",
                                                  tau_grid         = c(0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.001),
                                                  min_p            = 1e-15,
                                                  do_fix           = TRUE,
                                                  B_perm           = 0L,
                                                  perm_mode        = c("resample_global", "mvn"),
                                                  magma_genes_out  = NULL,
                                                  magma_cor_file   = NULL,
                                                  make_PD          = TRUE,
                                                  seed             = NULL,
                                                  analytic_logical = TRUE,
                                                  output           = FALSE,
                                                  out_dir          = "magcat_tfisher_soft_adaptive") {

  if (!requireNamespace("TFisher", quietly = TRUE)) {
    stop("magcat_soft_tfisher_adaptive_pathways(): package 'TFisher' is required.", call. = FALSE)
  }

  perm_mode <- match.arg(perm_mode)

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

  if (is.null(pathways) && !is.null(species)) {
    if (is.null(pmn_gene_col)) {
      pathways <- magcat_load_pathways(species = species)
    } else {
      pathways <- magcat_load_pathways(species = species, gene_col = pmn_gene_col)
    }
  }

  ## -------- standardize gene_results ----------
  needed_cols <- c(gene_col, p_col)
  if (!all(needed_cols %in% names(gene_results))) {
    stop("gene_results must contain columns '", gene_col, "' and '", p_col, "'.", call. = FALSE)
  }

  gr         <- gene_results
  genes_all  <- as.character(gr[[gene_col]])
  p_all      <- as.numeric(gr[[p_col]])
  genes_norm <- tolower(genes_all)

  # map from normalized -> canonical gene ID
  gene_map <- tapply(genes_all, genes_norm, function(x) x[1])

  # named vector of p's by normalized ID
  gene_p_vec <- stats::setNames(p_all, genes_norm)

  ## -------- pathways -> named list + names ----------
  if (is.data.frame(pathways)) {
    if (!("pathway_id" %in% names(pathways)) || !("gene_id" %in% names(pathways))) {
      stop("If 'pathways' is a data.frame it must have columns 'pathway_id' and 'gene_id'.", call. = FALSE)
    }
    if (!("pathway_name" %in% names(pathways))) pathways$pathway_name <- pathways$pathway_id

    p_list  <- split(pathways$gene_id, pathways$pathway_id)
    p_names <- tapply(pathways$pathway_name, pathways$pathway_id, FUN = function(x) x[1])
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

  # normalize IDs inside pathways
  p_list <- lapply(p_list, function(g) tolower(as.character(g)))

  ## -------- helpers ----------
  .tfisher_righttail_p <- function(p_i, tau1) {
    stat <- TFisher::stat.soft(p = p_i, tau1 = tau1)
    Fq   <- TFisher::p.soft(q = stat, n = length(p_i), tau1 = tau1, M = NULL)
    list(stat = stat, p = 1 - as.numeric(Fq))
  }

  .adaptive_minp <- function(p_i, tau_grid) {
    out_p    <- rep(NA_real_, length(tau_grid))
    out_stat <- rep(NA_real_, length(tau_grid))

    for (j in seq_along(tau_grid)) {
      tmp <- .tfisher_righttail_p(p_i, tau1 = tau_grid[j])
      out_p[j]    <- tmp$p
      out_stat[j] <- tmp$stat
    }

    jhat <- which.min(out_p)
    list(
      p_min   = out_p[jhat],
      tau_hat = tau_grid[jhat],
      stat_hat= out_stat[jhat]
    )
  }

  ## -------- pool for resample_global permutations ----------
  idx_pool <- which(is.finite(p_all) & !is.na(p_all))
  if (!length(idx_pool)) stop("No valid gene-level p-values found in 'gene_results'.", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  ## -------- result container ----------
  n_pw <- length(p_list)
  res  <- data.frame(
    pathway_id         = names(p_list),
    pathway_name       = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes            = NA_integer_,
    gene_names         = NA_character_,
    gene_pvals         = NA_character_,
    tau_hat            = NA_real_,
    tfisher_stat_hat   = NA_real_,
    tfisher_p_analytic = NA_real_,
    tfisher_p_perm     = NA_real_,
    stringsAsFactors   = FALSE
  )

  ## -------- MVN prep (once) ----------
  genes_out_tab <- NULL
  cor_pairs     <- NULL
  if (B_perm > 0L && perm_mode == "mvn") {
    if (is.null(magma_genes_out) || !file.exists(magma_genes_out)) {
      stop("perm_mode='mvn' requires magma_genes_out = path to merged MAGMA *.genes.out", call. = FALSE)
    }
    genes_out_tab <- magma_read_genes_out(magma_genes_out, gene_col = gene_col, chr_col = "CHR")

    if (!is.null(magma_cor_file)) {
      if (!file.exists(magma_cor_file)) stop("magma_cor_file does not exist: ", magma_cor_file, call. = FALSE)
      cor_pairs <- magma_read_gene_cor_pairs(magma_cor_file, gene1_col = 1, gene2_col = 2, r_col = 3)
    }
  }

  ## -------- main loop ----------
  for (i in seq_len(n_pw)) {

    genes_i_norm <- p_list[[i]]
    p_i <- gene_p_vec[genes_i_norm]

    keep <- !is.na(p_i)
    p_i  <- as.numeric(p_i[keep])
    genes_used_norm <- genes_i_norm[keep]

    d <- length(p_i)
    res$n_genes[i] <- d
    if (d == 0L) next

    canon_ids <- unname(gene_map[genes_used_norm])
    canon_ids <- unique(canon_ids[!is.na(canon_ids)])
    res$gene_names[i] <- paste(canon_ids, collapse = ";")

    if (do_fix) p_i <- fix_p_for_acat(p_i, min_p = min_p)
    res$gene_pvals[i] <- paste(p_i, collapse = ";")

    ## ---- adaptive analytic (min over tau of right-tail p) ----
    if (analytic_logical) {
      ad <- .adaptive_minp(p_i, tau_grid = tau_grid)
      res$tfisher_p_analytic <- res$tfisher_p_analytic  # keep column alive
      res$tfisher_p_analytic[i] <- ad$p_min
      res$tau_hat[i]            <- ad$tau_hat
      res$tfisher_stat_hat[i]   <- ad$stat_hat
    } else {
      # still record tau_hat/stat_hat using the same machinery
      ad <- .adaptive_minp(p_i, tau_grid = tau_grid)
      res$tau_hat[i]          <- ad$tau_hat
      res$tfisher_stat_hat[i] <- ad$stat_hat
      res$tfisher_p_analytic[i] <- NA_real_
    }

    ## ---- permutations (adaptive-calibrated) ----
    if (B_perm > 0L && is.finite(res$tfisher_p_analytic[i])) {

      T_obs <- res$tfisher_p_analytic[i]  # use min-analytic-p as the adaptive score

      T_null <- rep(NA_real_, B_perm)

      if (perm_mode == "resample_global") {

        for (b in seq_len(B_perm)) {
          idx_perm <- sample(idx_pool, size = d, replace = FALSE)
          p_perm   <- as.numeric(p_all[idx_perm])
          p_perm   <- p_perm[is.finite(p_perm) & !is.na(p_perm)]
          if (!length(p_perm)) next
          if (do_fix) p_perm <- fix_p_for_acat(p_perm, min_p = min_p)

          ad_b <- .adaptive_minp(p_perm, tau_grid = tau_grid)
          T_null[b] <- ad_b$p_min
        }

        res$tfisher_p_perm[i] <- (1 + sum(T_null <= T_obs, na.rm = TRUE)) / (B_perm + 1)

      } else if (perm_mode == "mvn") {

        # build R only for this pathway’s genes
        genes_S <- genes_used_norm

        R_S <- magma_build_R_for_pathway(
          genes_S       = genes_S,
          genes_out_tab = genes_out_tab,
          cor_pairs     = cor_pairs,
          gene_col      = gene_col,
          chr_col       = "CHR",
          make_PD       = make_PD
        )
        if (is.null(R_S)) next

        sim <- magma_simulate_null_Zp(
          genes_S       = rownames(R_S),
          genes_out_tab = genes_out_tab,
          R_S           = R_S,
          B             = B_perm,
          gene_col      = gene_col,
          chr_col       = "CHR"
        )

        for (b in seq_len(B_perm)) {
          p_perm <- as.numeric(sim$P[b, ])
          p_perm <- p_perm[is.finite(p_perm) & !is.na(p_perm)]
          if (!length(p_perm)) next
          if (do_fix) p_perm <- fix_p_for_acat(p_perm, min_p = min_p)

          ad_b <- .adaptive_minp(p_perm, tau_grid = tau_grid)
          T_null[b] <- ad_b$p_min
        }

        res$tfisher_p_perm[i] <- (1 + sum(T_null <= T_obs, na.rm = TRUE)) / (B_perm + 1)
      }
    }
  }

  ## -------- sort ----------
  if (B_perm > 0L) {
    ord <- order(res$tfisher_p_perm, decreasing = FALSE, na.last = TRUE)
  } else {
    ord <- order(res$tfisher_p_analytic, decreasing = FALSE, na.last = TRUE)
  }
  res <- res[ord, , drop = FALSE]

  ## -------- optional CSV output ----------
  if (output) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    species_tag <- if (is.null(species)) "custom" else species
    out_path <- file.path(out_dir, paste0("magcat_tfisher_soft_adaptive_pathways_", species_tag, ".csv"))
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
