#' Adjust MAGMA gene p-values using gene length and NSNPS
#'
#' @param gene_results data.frame of MAGMA gene results (e.g., genes_all)
#' @param gene_lengths data.frame of gene lengths from get_gene_lengths()
#' @param gene_col column in gene_results with gene IDs (default "GENE")
#' @param nsnp_col column in gene_results with NSNPS (default "NSNPS")
#' @param p_col column in gene_results with p-values (default "P")
#' @param z_col column in gene_results with Z-stat (default "ZSTAT"); if NULL/missing, reconstruct from P (magnitude)
#' @param len_gene_col column in gene_lengths with gene IDs (default "gene_id")
#' @param len_col column in gene_lengths with length (default "length")
#' @param log1p_covars logical; log1p transform NSNPS + length (default TRUE)
#'
#' @return data.frame with columns: gene_id, z_adj, p_adj, fit (list column)
#'         and attribute "lm_fit"
#' @export
magcat_adjust_gene_p <- function(gene_results,
                                 gene_lengths,
                                 gene_col     = "GENE",
                                 nsnp_col     = "NSNPS",
                                 p_col        = "P",
                                 z_col        = "ZSTAT",
                                 len_gene_col = "gene_id",
                                 len_col      = "length",
                                 log1p_covars = TRUE) {

  if (!is.data.frame(gene_results)) stop("magcat_adjust_gene_p(): gene_results must be a data.frame.", call. = FALSE)
  if (!is.data.frame(gene_lengths)) stop("magcat_adjust_gene_p(): gene_lengths must be a data.frame.", call. = FALSE)

  for (nm in c(gene_col, nsnp_col, p_col)) {
    if (!(nm %in% names(gene_results))) {
      stop("magcat_adjust_gene_p(): missing column in gene_results: ", nm, call. = FALSE)
    }
  }
  for (nm in c(len_gene_col, len_col)) {
    if (!(nm %in% names(gene_lengths))) {
      stop("magcat_adjust_gene_p(): missing column in gene_lengths: ", nm, call. = FALSE)
    }
  }

  gr <- gene_results
  gl <- gene_lengths

  gr[[gene_col]] <- as.character(gr[[gene_col]])
  gl[[len_gene_col]] <- as.character(gl[[len_gene_col]])

  # fix "gene:Zm..." -> "Zm..."
  gl[[len_gene_col]] <- sub("^gene:", "", gl[[len_gene_col]])

  # merge length onto results
  m <- merge(
    gr,
    gl[, c(len_gene_col, len_col), drop = FALSE],
    by.x = gene_col,
    by.y = len_gene_col,
    all.x = TRUE,
    sort = FALSE
  )

  # get z_raw (prefer supplied ZSTAT; else reconstruct from P)
  have_z <- (!is.null(z_col)) && (z_col %in% names(m))
  if (have_z) {
    z_raw <- suppressWarnings(as.numeric(m[[z_col]]))
  } else {
    p <- suppressWarnings(as.numeric(m[[p_col]]))
    p[p <= 0] <- 1e-300
    p[p >= 1] <- 1 - 1e-16
    z_raw <- stats::qnorm(1 - p/2)  # magnitude only
  }

  # residualize magnitude
  y <- abs(z_raw)

  nsnps <- suppressWarnings(as.numeric(m[[nsnp_col]]))
  glen  <- suppressWarnings(as.numeric(m[[len_col]]))

  if (log1p_covars) {
    nsnps <- log1p(nsnps)
    glen  <- log1p(glen)
  }

  ok <- is.finite(y) & is.finite(nsnps) & is.finite(glen)

  if (!any(ok)) {
    stop("magcat_adjust_gene_p(): no complete cases after merging covariates.", call. = FALSE)
  }

  df <- data.frame(y = y[ok], nsnps = nsnps[ok], len = glen[ok])

  fit <- stats::lm(y ~ nsnps + len, data = df)

  z_adj <- rep(NA_real_, nrow(m))
  z_adj[ok] <- stats::residuals(fit)

  p_adj <- rep(NA_real_, nrow(m))
  p_adj[ok] <- 2 * stats::pnorm(-abs(z_adj[ok]))

  out <- data.frame(
    gene_id = m[[gene_col]],
    z_adj   = z_adj,
    p_adj   = p_adj,
    stringsAsFactors = FALSE
  )

  # ✅ save model both ways
  attr(out, "lm_fit") <- fit
  out$fit <- list(fit)

  out
}
