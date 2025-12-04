#' Run MAGMA annotation
#'
#' @param snp_loc Path to SNP location file.
#' @param gene_loc Path to gene location file.
#' @param out_prefix Output prefix for MAGMA.
#' @param window Optional length-2 numeric vector (kb upstream, kb downstream).
#'   For example, window = c(25, 25) gives "window=25,25".
#' @export
magma_annotate <- function(snp_loc,
                           gene_loc,
                           out_prefix,
                           window = NULL) {
  mp <- magma_path()

  # base args
  args <- c("--annotate")

  # add window if given
  if (!is.null(window)) {
    if (length(window) != 2L) {
      stop("window must be NULL or a length-2 numeric vector, e.g. c(25, 25).")
    }
    args <- c(args, paste0("window=", window[1], ",", window[2]))
  }

  # file + out options
  args <- c(
    args,
    paste0("--snp-loc=", shQuote(snp_loc)),
    paste0("--gene-loc=", shQuote(gene_loc)),
    paste0("--out=",     shQuote(out_prefix))
  )

  system2(mp, args)
}


#' Run MAGMA gene analysis on SNP p-values
#'
#' @param bfile PLINK prefix for LD reference.
#' @param gene_annot MAGMA gene annotation file (.genes.annot).
#' @param pval_file SNP p-value file.
#' @param n_total Total sample size (N=... in MAGMA).
#' @param snp_col Column name in pval_file containing SNP IDs.
#' @param p_col Column name in pval_file containing p-values.
#' @param out_prefix Base output prefix. A suffix per gene model will be added.
#' @param gene_model Character vector of one or more MAGMA gene models,
#'   e.g. "snp-wise=top", "snp-wise=mean", "multi".
#'
#' @export
magma_gene <- function(bfile,
                       gene_annot,
                       pval_file,
                       n_total,
                       snp_col   = "SNP",
                       p_col     = "P",
                       out_prefix,
                       gene_model = "snp-wise=top") {
  mp <- magma_path()

  # ensure it's a vector
  gene_model <- as.character(gene_model)

  out_files <- vector("list", length(gene_model))

  for (i in seq_along(gene_model)) {
    gm <- gene_model[i]

    # nice suffix for output name (remove weird chars)
    gm_suffix <- gsub("[^A-Za-z0-9]+", "_", gm)
    out_i <- paste0(out_prefix, ".", gm_suffix)

    args <- c(
      paste0("--bfile=",      shQuote(bfile)),
      paste0("--gene-annot=", shQuote(gene_annot)),
      paste0("--pval=",       shQuote(pval_file),
             " N=", n_total,
             " use=", snp_col, ",", p_col),
      paste0("--gene-model=", gm),
      paste0("--out=",        shQuote(out_i))
    )

    system2(mp, args)

    out_files[[i]] <- list(
      gene_model = gm,
      out_prefix = out_i
    )
  }

  invisible(out_files)
}
