# ══════════════════════════════════════════════════════════════
# SalmonFlow — pipeline_functions.R
# Wrappers for external bioinformatics tools via processx
# ══════════════════════════════════════════════════════════════

#' Run FastQC on a set of FASTQ files
#' @param files Character vector of FASTQ file paths
#' @param outdir Output directory for FastQC results
#' @param threads Number of threads
#' @param log_callback Function(msg, type) for live logging
#' @return List with exit_status and outdir
run_fastqc <- function(files, outdir, threads = 4, log_callback = NULL) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(log_callback)) log_callback("FastQC: starting analysis...", "info")

  args <- c(files, "--outdir", outdir, "--threads", as.character(threads))

  result <- tryCatch(
    processx::run("fastqc", args = args, echo = FALSE,
                  stdout_line_callback = function(line) {
                    if (!is.null(log_callback)) log_callback(paste("FastQC:", line), "info")
                  },
                  stderr_line_callback = function(line) {
                    if (!is.null(log_callback)) log_callback(paste("FastQC:", line), "info")
                  }),
    error = function(e) list(status = 1, stderr = conditionMessage(e))
  )

  if (!is.null(log_callback)) {
    if (identical(result$status, 0L) || identical(result$status, 0)) {
      log_callback("FastQC: completed ✓", "success")
    } else {
      log_callback(paste("FastQC: error —", result$stderr %||% "unknown"), "error")
    }
  }

  list(exit_status = result$status %||% 1, outdir = outdir)
}

#' Run fastp on a single sample
#' @param r1 Path to R1 FASTQ
#' @param r2 Path to R2 FASTQ (NULL for SE)
#' @param out_dir Output directory for trimmed files
#' @param sample_name Sample name for output file naming
#' @param mode "PE" or "SE"
#' @param adapter_fasta Path to adapter FASTA (NULL = auto-detect for PE)
#' @param cut_front_quality --cut_front_mean_quality threshold
#' @param cut_tail_quality --cut_tail_mean_quality threshold
#' @param cut_right_quality --cut_right_mean_quality threshold (window = 4)
#' @param minlen --length_required threshold
#' @param threads Number of threads
#' @param log_callback Function(msg, type) for live logging
#' @return List with exit_status, r1_trimmed, r2_trimmed paths
run_fastp <- function(r1, r2 = NULL, out_dir, sample_name,
                      mode = "PE",
                      adapter_fasta = NULL,
                      cut_front_quality = 3, cut_tail_quality = 3,
                      cut_right_quality = 15, minlen = 36,
                      threads = 4, log_callback = NULL) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(log_callback)) log_callback(paste("fastp:", sample_name, "— starting..."), "info")

  if (mode == "PE") {
    r1_trimmed <- file.path(out_dir, paste0(sample_name, "_R1_trimmed.fastq.gz"))
    r2_trimmed <- file.path(out_dir, paste0(sample_name, "_R2_trimmed.fastq.gz"))
    args <- c("-i", r1, "-I", r2, "-o", r1_trimmed, "-O", r2_trimmed)
  } else {
    r1_trimmed <- file.path(out_dir, paste0(sample_name, "_trimmed.fastq.gz"))
    r2_trimmed <- NULL
    args <- c("-i", r1, "-o", r1_trimmed)
  }

  args <- c(args,
    "--cut_front",  "--cut_front_mean_quality",  as.character(cut_front_quality),
    "--cut_tail",   "--cut_tail_mean_quality",   as.character(cut_tail_quality),
    "--cut_right",  "--cut_right_window_size", "4",
                    "--cut_right_mean_quality",  as.character(cut_right_quality),
    "--length_required", as.character(minlen),
    "--thread", as.character(threads),
    "-j", file.path(out_dir, paste0(sample_name, "_fastp.json")),
    "-h", file.path(out_dir, paste0(sample_name, "_fastp.html"))
  )

  if (!is.null(adapter_fasta) && file.exists(adapter_fasta)) {
    args <- c(args, "--adapter_fasta", adapter_fasta)
  } else if (mode == "PE") {
    args <- c(args, "--detect_adapter_for_pe")
  }

  result <- tryCatch(
    processx::run("fastp", args = args, echo = FALSE,
                  stderr_line_callback = function(line) {
                    if (!is.null(log_callback)) log_callback(paste("fastp:", line), "info")
                  }),
    error = function(e) list(status = 1, stderr = conditionMessage(e))
  )

  if (!is.null(log_callback)) {
    if (identical(result$status, 0L) || identical(result$status, 0)) {
      log_callback(paste("fastp:", sample_name, "✓"), "success")
    } else {
      log_callback(paste("fastp:", sample_name, "error —", result$stderr %||% "unknown"), "error")
    }
  }

  list(
    exit_status = result$status %||% 1,
    r1_trimmed  = r1_trimmed,
    r2_trimmed  = r2_trimmed
  )
}

#' Build a Salmon index
#' @param fasta Path to transcriptome FASTA
#' @param outdir Output directory for the index
#' @param decoy Path to decoy file (optional)
#' @param kmer k-mer size (default 31)
#' @param threads Number of threads
#' @param log_callback Function(msg, type) for live logging
#' @return List with exit_status and index_dir
build_salmon_index <- function(fasta, outdir, decoy = NULL,
                               kmer = 31, threads = 4,
                               log_callback = NULL) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(log_callback)) log_callback("Salmon index: building...", "info")

  args <- c("index",
            "-t", fasta,
            "-i", outdir,
            "--threads", as.character(threads),
            "-k", as.character(kmer))

  if (!is.null(decoy) && file.exists(decoy)) {
    args <- c(args, "-d", decoy)
  }

  result <- tryCatch(
    processx::run("salmon", args = args, echo = FALSE,
                  stderr_line_callback = function(line) {
                    if (!is.null(log_callback)) log_callback(paste("Salmon index:", line), "info")
                  }),
    error = function(e) list(status = 1, stderr = conditionMessage(e))
  )

  if (!is.null(log_callback)) {
    if (identical(result$status, 0L) || identical(result$status, 0)) {
      log_callback("Salmon index: completed ✓", "success")
    } else {
      log_callback(paste("Salmon index: error —", result$stderr %||% "unknown"), "error")
    }
  }

  list(exit_status = result$status %||% 1, index_dir = outdir)
}

#' Run Salmon quant on a single sample
#' @param index_dir Path to Salmon index
#' @param r1 Path to R1 FASTQ (trimmed)
#' @param r2 Path to R2 FASTQ (trimmed, NULL for SE)
#' @param outdir Output directory for this sample
#' @param sample_name Sample identifier
#' @param lib_type Library type string (e.g. "A")
#' @param gc_bias Logical, enable gcBias correction
#' @param seq_bias Logical, enable seqBias correction
#' @param threads Number of threads
#' @param is_se Logical, TRUE for single-end
#' @param log_callback Function(msg, type) for live logging
#' @return List with exit_status and quant_dir
run_salmon_quant <- function(index_dir, r1, r2 = NULL, outdir, sample_name,
                             lib_type = "A", gc_bias = TRUE, seq_bias = TRUE,
                             threads = 4, is_se = FALSE,
                             log_callback = NULL) {

  sample_out <- file.path(outdir, sample_name)
  dir.create(sample_out, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(log_callback)) log_callback(paste("Salmon quant:", sample_name, "— running..."), "info")

  args <- c("quant",
            "-i", index_dir,
            "-l", lib_type,
            "-p", as.character(threads),
            "--validateMappings",
            "-o", sample_out)

  if (is_se) {
    args <- c(args, "-r", r1)
  } else {
    args <- c(args, "-1", r1, "-2", r2)
  }

  if (gc_bias)  args <- c(args, "--gcBias")
  if (seq_bias) args <- c(args, "--seqBias")

  result <- tryCatch(
    processx::run("salmon", args = args, echo = FALSE,
                  stderr_line_callback = function(line) {
                    if (!is.null(log_callback)) log_callback(paste("Salmon quant:", line), "info")
                  }),
    error = function(e) list(status = 1, stderr = conditionMessage(e))
  )

  # Parse mapping rate
  meta <- parse_salmon_meta(sample_out)

  if (!is.null(log_callback)) {
    if (identical(result$status, 0L) || identical(result$status, 0)) {
      rate_msg <- if (!is.na(meta$percent_mapped)) paste0("mapping ", meta$percent_mapped, "%") else ""
      log_callback(paste("Salmon quant:", sample_name, "✓", rate_msg), "success")
      if (!is.na(meta$percent_mapped) && meta$percent_mapped < 50) {
        log_callback(paste("⚠ Warning:", sample_name, "mapping rate below 50%!"), "warn")
      }
    } else {
      log_callback(paste("Salmon quant:", sample_name, "error —", result$stderr %||% "unknown"), "error")
    }
  }

  list(
    exit_status = result$status %||% 1,
    quant_dir   = sample_out,
    meta        = meta
  )
}

#' Run MultiQC on a directory tree
#' @param input_dir Directory to scan for tool outputs
#' @param outdir Output directory for MultiQC report
#' @param log_callback Function(msg, type) for live logging
#' @return List with exit_status and report_path
run_multiqc <- function(input_dir, outdir, log_callback = NULL) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(log_callback)) log_callback("MultiQC: generating report...", "info")

  args <- c(input_dir,
            "--outdir", outdir,
            "--filename", "multiqc_report",
            "--force")

  result <- tryCatch(
    processx::run("multiqc", args = args, echo = FALSE,
                  stderr_line_callback = function(line) {
                    if (!is.null(log_callback)) log_callback(paste("MultiQC:", line), "info")
                  }),
    error = function(e) list(status = 1, stderr = conditionMessage(e))
  )

  report_path <- file.path(outdir, "multiqc_report.html")

  if (!is.null(log_callback)) {
    if (file.exists(report_path)) {
      log_callback("MultiQC: report generated ✓", "success")
    } else {
      log_callback("MultiQC: report generation may have failed", "warn")
    }
  }

  list(exit_status = result$status %||% 1, report_path = report_path)
}
