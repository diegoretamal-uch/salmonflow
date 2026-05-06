# ══════════════════════════════════════════════════════════════
# SalmonFlow — helpers.R
# Utility functions: validation, file helpers, log formatting
# ══════════════════════════════════════════════════════════════

#' Validate that FASTQ files exist and have correct extensions
#' @param paths Character vector of file paths
#' @return List with $valid (logical), $messages (character vector)
validate_fastq_files <- function(paths) {
  msgs <- character(0)
  valid <- TRUE

  for (p in paths) {
    if (!file.exists(p)) {
      msgs <- c(msgs, paste0("[ERROR] Not found: ", basename(p)))
      valid <- FALSE
    } else if (!grepl("\\.(fastq|fq)(\\.gz)?$", p, ignore.case = TRUE)) {
      msgs <- c(msgs, paste0("[WARN] Not a FASTQ file: ", basename(p)))
      valid <- FALSE
    } else {
      msgs <- c(msgs, paste0("[OK] ", basename(p)))
    }
  }

  list(valid = valid, messages = msgs)
}

#' Auto-detect paired FASTQ files from a directory
#' @param dir_path Path to directory containing FASTQs
#' @return data.frame with columns: name, r1, r2
detect_fastq_pairs <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.(fastq|fq)(\\.gz)?$",
                      full.names = TRUE, ignore.case = TRUE)

  if (length(files) == 0) return(data.frame(name = character(0),
                                             r1 = character(0),
                                             r2 = character(0),
                                             group = character(0),
                                             stringsAsFactors = FALSE))

  basenames <- basename(files)

  # Try to find R1/R2 pairs using common naming conventions
  r1_pattern <- "(_R1[_.]|_1\\.|_R1\\.|\\.R1[_.])"
  r2_pattern <- "(_R2[_.]|_2\\.|_R2\\.|\\.R2[_.])"

  r1_files <- files[grepl(r1_pattern, basenames)]
  r2_files <- files[grepl(r2_pattern, basenames)]

  if (length(r1_files) > 0 && length(r1_files) == length(r2_files)) {
    # Paired-end detected
    sample_names <- gsub(r1_pattern, "_", basename(r1_files))
    sample_names <- gsub("\\.(fastq|fq)(\\.gz)?$", "", sample_names, ignore.case = TRUE)
    sample_names <- gsub("_$", "", sample_names)

    # Sort to ensure matching
    r1_files <- sort(r1_files)
    r2_files <- sort(r2_files)

    data.frame(
      name  = sample_names,
      r1    = r1_files,
      r2    = r2_files,
      group = "",
      stringsAsFactors = FALSE
    )
  } else {
    # Single-end fallback
    sample_names <- gsub("\\.(fastq|fq)(\\.gz)?$", "", basenames, ignore.case = TRUE)
    data.frame(
      name  = sample_names,
      r1    = files,
      r2    = rep(NA_character_, length(files)),
      group = "",
      stringsAsFactors = FALSE
    )
  }
}

#' Parse Salmon meta_info.json for mapping rate
#' @param quant_dir Path to sample's Salmon quant output directory
#' @return Named list with mapping stats
parse_salmon_meta <- function(quant_dir) {
  meta_path <- file.path(quant_dir, "aux_info", "meta_info.json")
  if (!file.exists(meta_path)) {
    return(list(
      num_processed   = NA,
      num_mapped      = NA,
      mapping_rate    = NA,
      percent_mapped  = NA
    ))
  }

  meta <- jsonlite::fromJSON(meta_path)
  mapped    <- meta$num_mapped %||% 0
  processed <- meta$num_processed %||% 1
  rate      <- round(mapped / processed * 100, 2)

  list(
    num_processed  = processed,
    num_mapped     = mapped,
    mapping_rate   = mapped / processed,
    percent_mapped = rate
  )
}

#' Format a log message with timestamp
#' @param msg Message text
#' @param type One of "info", "success", "error", "warn"
#' @return HTML-formatted log line
timestamp_log <- function(msg, type = "info") {
  ts <- format(Sys.time(), "[%H:%M:%S]")
  css_class <- switch(type,
    success = "log-success",
    error   = "log-error",
    warn    = "log-warn",
    "log-info"
  )
  paste0(
    '<span class="log-timestamp">', ts, '</span> ',
    '<span class="', css_class, '">', htmltools::htmlEscape(msg), '</span>'
  )
}

#' Null-coalescing operator
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
