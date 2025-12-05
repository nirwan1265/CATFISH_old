#' Convert a GFF3 file to a MAGMA gene location file
#'
#' @param gff Path to a GFF3 file.
#' @param out Path to write the MAGMA gene location file (tab-delimited).
#' @param feature_type GFF3 feature type to use (default "gene").
#' @param id_fields Candidate attribute tags to use as gene IDs
#'   (searched in order), e.g. "gene_id", "ID", "Name".
#' @param drop_scaffolds Logical; if TRUE (default), drop seqids that look
#'   like scaffolds (seqid starting with "scaf", case-insensitive).
#'
#' @return Invisibly returns a data.frame with columns GENE, CHR, START, STOP.
#' @export
gff3_to_geneloc <- function(gff,
                            out,
                            feature_type = "gene",
                            id_fields = c("gene_id", "ID", "Name"),
                            drop_scaffolds = TRUE) {

  # read GFF3, skipping comment lines starting with '#'
  x <- utils::read.delim(
    gff,
    header           = FALSE,
    comment.char     = "#",
    stringsAsFactors = FALSE,
    sep              = "\t"
  )

  if (ncol(x) < 9L) {
    stop("GFF3 file must have at least 9 columns.")
  }

  colnames(x)[1:9] <- c(
    "seqid", "source", "type", "start",
    "end", "score", "strand", "phase", "attr"
  )

  # keep only the desired feature type(s)
  keep <- x$type %in% feature_type
  x <- x[keep, , drop = FALSE]

  if (!nrow(x)) {
    stop("No rows with type %s found in GFF3.", feature_type)
  }

  # optionally drop scaffold seqids
  if (drop_scaffolds) {
    not_scaf <- !grepl("^scaf", x$seqid, ignore.case = TRUE)
    x <- x[not_scaf, , drop = FALSE]

    if (!nrow(x)) {
      stop("All entries were removed by drop_scaffolds=TRUE; ",
           "no non-scaffold seqids left.")
    }
  }

  # helper to parse the attributes string and extract ID
  parse_id <- function(attr_string, id_fields) {
    if (is.na(attr_string) || attr_string == "") return(NA_character_)

    parts <- strsplit(attr_string, ";", fixed = TRUE)[[1]]
    kv <- strsplit(parts, "=", fixed = TRUE)

    keys <- vapply(kv, function(z) if (length(z) >= 1L) z[1] else "", character(1))
    vals <- vapply(kv, function(z) if (length(z) >= 2L) z[2] else NA_character_, character(1))

    for (tag in id_fields) {
      idx <- which(keys == tag)
      if (length(idx) > 0L && !is.na(vals[idx[1]])) {
        return(vals[idx[1]])
      }
    }

    NA_character_
  }

  gene_ids <- vapply(x$attr, parse_id, character(1), id_fields = id_fields)

  res <- data.frame(
    GENE  = gene_ids,
    CHR   = x$seqid,
    START = x$start,
    STOP  = x$end,
    stringsAsFactors = FALSE
  )

  # drop rows without a gene ID
  res <- res[!is.na(res$GENE) & res$GENE != "", , drop = FALSE]

  if (!nrow(res)) {
    stop("Could not extract any gene IDs from attributes using tags: ",
         paste(id_fields, collapse = ", "))
  }

  ## ----- NEW: order by CHR, START -----
  # if CHR looks like "chr1", "chr2", ..., sort numerically by that
  chr_pattern <- grepl("^chr[0-9]+$", res$CHR, ignore.case = TRUE)

  if (all(chr_pattern)) {
    chr_num <- as.integer(sub("^chr", "", tolower(res$CHR)))
    o <- order(chr_num, res$START)
  } else {
    # fallback: alphabetical CHR, then START
    o <- order(res$CHR, res$START)
  }

  res <- res[o, , drop = FALSE]
  ## ------------------------------------

  # write MAGMA gene-loc file
  utils::write.table(
    res,
    file      = out,
    quote     = FALSE,
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  invisible(res)
}
