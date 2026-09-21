# ══════════════════════════════════════════════════════════════
# SalmonFlow — tximport_utils.R
# Wrappers for tx2gene construction and tximport
# ══════════════════════════════════════════════════════════════

#' Cache file path for a GTF's tx2gene map.
#' Keyed on the annotation's path, size and mtime, so a replaced or edited GTF
#' misses the cache. Lives beside the references (an auxiliary file, under
#' 1 MB — it has no bearing on the sample-file storage policy).
#' The key is spelled into the filename rather than hashed: it stays readable,
#' needs no hashing dependency, and size+mtime is enough to notice a GTF that
#' was replaced or edited.
tx2gene_cache_path <- function(gtf_path, cache_dir = NULL) {
  info <- file.info(gtf_path)
  if (is.na(info$size)) stop("GTF not found: ", gtf_path)

  stem <- gsub("[^A-Za-z0-9._-]", "_", basename(gtf_path))
  key  <- paste0(stem, "_", info$size, "_",
                 format(info$mtime, "%Y%m%d%H%M%S"))
  dir  <- if (is.null(cache_dir)) dirname(gtf_path) else cache_dir
  file.path(dir, paste0(".tx2gene_", key, ".rds"))
}

#' Build tx2gene data.frame from a GTF file
#' @param gtf_path Path to GTF annotation file
#' @param log_callback Function(msg, type) for live logging
#' @param use_cache Logical, read/write the .rds cache beside the GTF
#' @return data.frame with columns TXNAME, GENEID
build_tx2gene <- function(gtf_path, log_callback = NULL, use_cache = TRUE) {

  # ── Cache lookup ────────────────────────────────────────────
  # makeTxDbFromGFF is single-threaded and unavoidably serial (measured at
  # ~40s on GENCODE vM34), sitting between Salmon finishing and MultiQC
  # starting with the machine otherwise idle. The parsed result is identical
  # for a given annotation, so it only needs deriving once.
  cache_file <- if (isTRUE(use_cache)) {
    tryCatch(tx2gene_cache_path(gtf_path), error = function(e) NULL)
  } else NULL

  if (!is.null(cache_file) && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    # A corrupt or truncated cache must never abort the run — fall through
    # to a full parse instead.
    if (is.data.frame(cached) && all(c("TXNAME", "GENEID") %in% names(cached))) {
      if (!is.null(log_callback)) {
        log_callback(paste("tx2gene: loaded from cache —", nrow(cached),
                           "transcripts (skipped GTF parse)"), "success")
      }
      return(cached)
    }
    if (!is.null(log_callback)) {
      log_callback("tx2gene: cache unreadable, rebuilding from GTF", "warn")
    }
  }

  if (!is.null(log_callback)) log_callback("tximport: building tx2gene from GTF...", "info")

  suppressPackageStartupMessages({
    library(txdbmaker)
    library(GenomicFeatures)
  })

  txdb <- tryCatch(
    txdbmaker::makeTxDbFromGFF(gtf_path, format = "auto"),
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

  # Strip version suffix from transcript IDs (e.g. ENST00000832824.1 -> ENST00000832824)
  # to match index FASTAs that lack version numbers
  tx2gene$TXNAME <- sub("\\.[0-9]+$", "", tx2gene$TXNAME)

  if (!is.null(log_callback)) {
    log_callback(
      paste("tx2gene: mapped", nrow(tx2gene), "transcripts to",
            length(unique(tx2gene$GENEID)), "genes"),
      "success"
    )
  }

  # ── Cache write ─────────────────────────────────────────────
  # Best-effort: a read-only references directory must not fail the run.
  if (!is.null(cache_file)) {
    ok <- tryCatch({ saveRDS(tx2gene, cache_file); TRUE },
                   error = function(e) FALSE)
    if (!is.null(log_callback)) {
      if (isTRUE(ok)) {
        log_callback(paste("tx2gene: cached for future runs —", basename(cache_file)), "info")
      } else {
        log_callback("tx2gene: could not write cache (references not writable?)", "warn")
      }
    }
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
      paste0("tximport: merged matrix — ", nrow(count_matrix), " genes x ",
             length(sample_names), " samples"),
      "success"
    )
  }

  count_matrix
}
