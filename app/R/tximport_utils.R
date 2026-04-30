# ══════════════════════════════════════════════════════════════
# SalmonFlow — tximport_utils.R
# Wrappers for tx2gene construction and tximport
# ══════════════════════════════════════════════════════════════

#' Build tx2gene data.frame from a GTF file
#' @param gtf_path Path to GTF annotation file
#' @param log_callback Function(msg, type) for live logging
#' @return data.frame with columns TXNAME, GENEID
build_tx2gene <- function(gtf_path, log_callback = NULL) {
  if (!is.null(log_callback)) log_callback("tximport: building tx2gene from GTF...", "info")

  suppressPackageStartupMessages({
    library(txdbmaker)
    library(GenomicFeatures)
  })

  txdb <- tryCatch(
    txdbmaker::makeTxDbFromGFF(gtf_path, format = "GTF"),
    error = function(e) {
      if (!is.null(log_callback)) log_callback(paste("tx2gene error:", e$message), "error")
      return(NULL)
    }
  )

  if (is.null(txdb)) return(NULL)

  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

  # Keep only the two required columns

  tx2gene <- tx2gene[, c("TXNAME", "GENEID")]
  tx2gene <- tx2gene[complete.cases(tx2gene), ]

  if (!is.null(log_callback)) {
    log_callback(
      paste("tx2gene: mapped", nrow(tx2gene), "transcripts to",
            length(unique(tx2gene$GENEID)), "genes ✓"),
      "success"
    )
  }

  tx2gene
}

#' Run tximport to produce merged count matrix
#' @param quant_dir Base directory containing per-sample Salmon quant folders
#' @param sample_names Character vector of sample names (subdirectory names)
#' @param tx2gene data.frame with TXNAME and GENEID columns
#' @param method countsFromAbundance method (default "lengthScaledTPM")
#' @param ignore_tx_version Logical, strip transcript version suffixes
#' @param output_dir Directory to save the CSV output
#' @param log_callback Function(msg, type) for live logging
#' @return data.frame with the merged count matrix, or NULL on error
run_tximport <- function(quant_dir, sample_names, tx2gene,
                         method = "lengthScaledTPM",
                         ignore_tx_version = TRUE,
                         output_dir = "/data/output",
                         log_callback = NULL) {

  if (!is.null(log_callback)) log_callback("tximport: importing Salmon quantifications...", "info")

  suppressPackageStartupMessages(library(tximport))

  # Build paths to quant.sf files
  quant_files <- file.path(quant_dir, sample_names, "quant.sf")
  names(quant_files) <- sample_names

  # Verify all quant.sf files exist
  missing <- !file.exists(quant_files)
  if (any(missing)) {
    msg <- paste("Missing quant.sf for:", paste(sample_names[missing], collapse = ", "))
    if (!is.null(log_callback)) log_callback(msg, "error")
    return(NULL)
  }

  # Run tximport
  txi <- tryCatch(
    tximport(
      quant_files,
      type = "salmon",
      tx2gene = tx2gene,
      countsFromAbundance = method,
      ignoreTxVersion = ignore_tx_version
    ),
    error = function(e) {
      if (!is.null(log_callback)) log_callback(paste("tximport error:", e$message), "error")
      return(NULL)
    }
  )

  if (is.null(txi)) return(NULL)

  # Build output data.frame
  count_matrix <- as.data.frame(txi$counts)
  count_matrix$gene_id <- rownames(count_matrix)
  count_matrix <- count_matrix[, c("gene_id", sample_names)]

  # Save to CSV
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  csv_path <- file.path(output_dir, "merged_lengthScaledTPM.csv")
  write.csv(count_matrix, csv_path, row.names = FALSE)

  if (!is.null(log_callback)) {
    log_callback(
      paste0("tximport: merged matrix — ", nrow(count_matrix), " genes × ",
             length(sample_names), " samples ✓"),
      "success"
    )
  }

  count_matrix
}
