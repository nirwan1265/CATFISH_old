#' Run MAGMA competitive gene-set (pathway) analysis
#' @export
magma_geneset_competitive <- function(gene_results_raw,
                                      set_annot,
                                      out_prefix,
                                      out_dir = NULL) {
  mp <- magma_path()

  if (!file.exists(gene_results_raw)) stop("gene_results_raw not found: ", gene_results_raw, call.=FALSE)
  if (!file.exists(set_annot))        stop("set_annot not found: ", set_annot, call.=FALSE)

  prefix_full <- if (is.null(out_dir)) out_prefix else {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    file.path(out_dir, out_prefix)
  }

  args <- c(
    "--gene-results", gene_results_raw,
    "--set-annot",    set_annot,
    "--out",          prefix_full
  )

  system2(mp, args)
  invisible(prefix_full)
}
