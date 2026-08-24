# ══════════════════════════════════════════════════════════════
# SalmonFlow — idep_export.R
# Prepare pipeline outputs for downstream analysis in iDEP
# ══════════════════════════════════════════════════════════════
#
# iDEP (https://github.com/gexijin/idepGolem) expects a gene-level
# expression matrix with gene IDs in the first column and one column
# per sample. SalmonFlow's matrix is already that shape, but two
# details trip iDEP up:
#
#   1. Gene IDs carry Ensembl/GENCODE version suffixes (ENSG…​.16),
#      which lowers iDEP's gene ID -> Ensembl match rate.
#   2. tximport values are non-integer, while DESeq2 inside iDEP
#      expects integer counts.
#
# These helpers produce an iDEP-ready copy. The pipeline output
# (merged_lengthScaledTPM.csv) is never modified.

#' Prepare an iDEP-ready count matrix
#'
#' @param count_matrix data.frame — first column `gene_id`, then one
#'   numeric column per sample (as produced by run_tximport()).
#' @param strip_version Logical, drop Ensembl version suffixes.
#' @param round_counts Logical, round values to whole numbers.
#' @return data.frame in the same shape, with attribute `merged_ids`
#'   holding the number of rows collapsed by de-duplication, or NULL
#'   if the input is empty.
prepare_idep_counts <- function(count_matrix,
                                strip_version = TRUE,
                                round_counts  = TRUE) {

  if (is.null(count_matrix) || nrow(count_matrix) == 0) return(NULL)
  if (!"gene_id" %in% names(count_matrix)) return(NULL)

  sample_cols <- setdiff(names(count_matrix), "gene_id")
  if (length(sample_cols) == 0) return(NULL)

  ids <- as.character(count_matrix$gene_id)

  mat <- as.matrix(count_matrix[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  mat[is.na(mat)] <- 0

  merged_ids <- 0L

  if (strip_version) {
    # ENSG00000000003.16       -> ENSG00000000003
    # ENSG00000182378.14_PAR_Y -> ENSG00000182378_PAR_Y
    # PAR_Y entries are genuinely distinct loci, so the suffix is kept.
    ids <- sub("\\.[0-9]+(_PAR_Y)?$", "\\1", ids)

    # Stripping versions can collide rows. Duplicated gene IDs make
    # iDEP behave unpredictably, so collapse them by summing.
    dup <- duplicated(ids)
    if (any(dup)) {
      merged_ids <- sum(dup)
      mat <- rowsum(mat, group = ids, reorder = FALSE)
      ids <- rownames(mat)
    }
  }

  if (round_counts) {
    mat <- round(mat)
    mat[mat < 0] <- 0
  }

  out <- data.frame(gene_id = ids, mat,
                    check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- NULL

  # Preserve the original sample order
  out <- out[, c("gene_id", sample_cols), drop = FALSE]

  attr(out, "merged_ids") <- merged_ids
  out
}

#' Build an iDEP experiment design table from the sample table
#'
#' iDEP's optional design file lists one row per sample: the first
#' column holds sample names matching the expression matrix column
#' headers, and each further column is an experimental factor.
#'
#' SalmonFlow already collects a `group` per sample in the Samples
#' tab, but never writes it to disk — this is the only place it
#' becomes a durable artifact.
#'
#' @param samples data.frame with columns name, r1, r2, group
#'   (see detect_fastq_pairs() in helpers.R).
#' @return data.frame with columns Sample and Group, carrying
#'   attribute `has_groups` (FALSE when the user left the group
#'   column blank), or NULL if there are no samples.
build_idep_design <- function(samples) {

  if (is.null(samples) || nrow(samples) == 0) return(NULL)

  grp <- as.character(samples$group %||% rep("", nrow(samples)))
  grp[is.na(grp)] <- ""
  grp <- trimws(grp)

  has_groups <- any(nzchar(grp))
  grp[!nzchar(grp)] <- "Ungrouped"

  out <- data.frame(
    Sample = as.character(samples$name),
    Group  = grp,
    stringsAsFactors = FALSE
  )

  attr(out, "has_groups") <- has_groups
  out
}

#' Check whether a local iDEP instance is reachable
#'
#' Called from the SalmonFlow container, so it targets the compose
#' service name first, then the Docker Desktop host alias for users
#' who started iDEP outside compose (e.g. via run.sh).
#'
#' A negative result is advisory only — it cannot distinguish "iDEP
#' is down" from "this container has no route to it", so the UI
#' warns rather than blocking.
#'
#' @return TRUE if a TCP connection succeeded.
idep_reachable <- function(targets = list(
                             list(host = "idep",                 port = 3838),
                             list(host = "host.docker.internal", port = 3839)
                           )) {
  for (t in targets) {
    ok <- tryCatch({
      con <- suppressWarnings(
        socketConnection(host = t$host, port = t$port,
                         timeout = 1, open = "r+", blocking = TRUE)
      )
      close(con)
      TRUE
    }, error = function(e) FALSE)

    if (isTRUE(ok)) return(TRUE)
  }
  FALSE
}
