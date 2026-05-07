# ══════════════════════════════════════════════════════════════
# SalmonFlow — pipeline_functions.R
# Wrappers for external bioinformatics tools via processx
# ══════════════════════════════════════════════════════════════

#' Detect GENCODE-format FASTA (headers use | as field separator)
detect_gencode_fasta <- function(fasta_path) {
  tryCatch({
    first_line <- if (grepl("\\.gz$", fasta_path, ignore.case = TRUE)) {
      con <- gzcon(file(fasta_path, "rb"))
      on.exit(close(con))
      readLines(con, n = 1, warn = FALSE)
    } else {
      readLines(fasta_path, n = 1, warn = FALSE)
    }
    grepl("^>.*\\|", first_line)
  }, error = function(e) FALSE)
}

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

  result <- processx::run("fastqc", args = args, echo = FALSE,
                          error_on_status = FALSE,
                          stdout_line_callback = function(line, proc) {
                            if (!is.null(log_callback)) log_callback(paste("FastQC:", line), "info")
                          },
                          stderr_line_callback = function(line, proc) {
                            if (!is.null(log_callback)) log_callback(paste("FastQC:", line), "info")
                          })

  if (!is.null(log_callback)) {
    if (identical(result$status, 0L) || identical(result$status, 0)) {
      log_callback("FastQC: completed", "success")
    } else {
      log_callback(paste("FastQC: failed with exit code", result$status), "error")
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

  result <- processx::run("fastp", args = args, echo = FALSE,
                          error_on_status = FALSE,
                          stderr_line_callback = function(line, proc) {
                            if (!is.null(log_callback)) log_callback(paste("fastp:", line), "info")
                          })

  if (!is.null(log_callback)) {
    if (identical(result$status, 0L) || identical(result$status, 0)) {
      log_callback(paste("fastp:", sample_name, "done"), "success")
    } else {
      log_callback(paste("fastp:", sample_name, "failed with exit code", result$status), "error")
    }
  }

  list(
    exit_status = result$status %||% 1,
    r1_trimmed  = r1_trimmed,
    r2_trimmed  = r2_trimmed
  )
}

#' Build a Salmon index
#' @param fasta Path to transcriptome FASTA (plain or .gz)
#' @param outdir Output directory for the index
#' @param decoy Path to GENOME FASTA for decoy-aware indexing (optional).
#'   The function auto-generates the required gentrome + decoys.txt.
#' @param kmer k-mer size (default 31)
#' @param threads Number of threads
#' @param log_callback Function(msg, type) for live logging
#' @return List with exit_status and index_dir
build_salmon_index <- function(fasta, outdir, decoy = NULL,
                               kmer = 31, threads = 4,
                               log_callback = NULL) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(log_callback)) log_callback("Salmon index: building...", "info")

  t_arg <- fasta
  d_arg <- NULL

  # ── Decoy-aware: build gentrome + decoys.txt automatically ──
  # Salmon requires -t <transcriptome+genome.fa> -d <decoy_names.txt>
  if (!is.null(decoy) && nchar(decoy) > 0 && file.exists(decoy)) {
    prep_dir      <- file.path(dirname(outdir), "salmon_index_prep")
    dir.create(prep_dir, showWarnings = FALSE, recursive = TRUE)
    decoys_file   <- file.path(prep_dir, "decoys.txt")
    gentrome_file <- file.path(prep_dir, "gentrome.fa")

    # 1. Extract sequence names from the genome FASTA header lines
    if (!is.null(log_callback)) log_callback("Salmon index: extracting decoy sequence names...", "info")
    g_gz <- grepl("\\.gz$", decoy, ignore.case = TRUE)
    names_cmd <- if (g_gz) {
      paste0("zcat ", shQuote(decoy), " | grep '^>' | cut -d' ' -f1 | sed 's/^>//'")
    } else {
      paste0("grep '^>' ", shQuote(decoy), " | cut -d' ' -f1 | sed 's/^>//'")
    }
    names_res <- processx::run("bash", args = c("-c", names_cmd), error_on_status = FALSE)
    if (!identical(names_res$status, 0L) && !identical(names_res$status, 0)) {
      if (!is.null(log_callback)) log_callback("Salmon index: failed to extract decoy names", "error")
      return(list(exit_status = 1, index_dir = outdir))
    }
    decoy_names <- Filter(nchar, strsplit(trimws(names_res$stdout), "\n")[[1]])
    writeLines(decoy_names, decoys_file)
    if (!is.null(log_callback)) log_callback(
      paste("Salmon index:", length(decoy_names), "decoy sequences identified"), "info")

    # 2. Concatenate transcriptome + genome into a single gentrome FASTA
    if (!is.null(log_callback)) log_callback(
      "Salmon index: concatenating transcriptome + genome into gentrome (may take a few minutes)...", "info")
    t_gz <- grepl("\\.gz$", fasta, ignore.case = TRUE)

    concat_cmd <- if (t_gz && g_gz) {
      gentrome_file <- paste0(gentrome_file, ".gz")
      paste0("cat ", shQuote(fasta), " ", shQuote(decoy), " > ", shQuote(gentrome_file))
    } else if (!t_gz && !g_gz) {
      paste0("cat ", shQuote(fasta), " ", shQuote(decoy), " > ", shQuote(gentrome_file))
    } else if (t_gz) {
      paste0("{ zcat ", shQuote(fasta), "; cat ",  shQuote(decoy), "; } > ", shQuote(gentrome_file))
    } else {
      paste0("{ cat ",  shQuote(fasta), "; zcat ", shQuote(decoy), "; } > ", shQuote(gentrome_file))
    }

    concat_res <- processx::run("bash", args = c("-c", concat_cmd), error_on_status = FALSE)
    if (!identical(concat_res$status, 0L) && !identical(concat_res$status, 0)) {
      if (!is.null(log_callback)) log_callback("Salmon index: failed to create gentrome", "error")
      return(list(exit_status = 1, index_dir = outdir))
    }
    if (!is.null(log_callback)) log_callback("Salmon index: gentrome ready", "success")

    t_arg <- gentrome_file
    d_arg <- decoys_file
  }

  # ── Run salmon index ──────────────────────────────────────
  args <- c("index",
            "-t", t_arg,
            "-i", outdir,
            "--threads", as.character(threads),
            "-k", as.character(kmer))

  if (!is.null(d_arg)) args <- c(args, "-d", d_arg)

  # Auto-detect GENCODE format from the original transcriptome FASTA headers.
  # GENCODE uses | as field separator; --gencode tells Salmon to use only the
  # first field (ENST ID) as the transcript name instead of the full header.
  if (detect_gencode_fasta(fasta)) {
    args <- c(args, "--gencode")
    if (!is.null(log_callback)) log_callback(
      "Salmon index: GENCODE format detected — adding --gencode flag", "info")
  }

  result <- processx::run("salmon", args = args, echo = FALSE,
                          error_on_status = FALSE,
                          stderr_line_callback = function(line, proc) {
                            if (!is.null(log_callback)) log_callback(paste("Salmon index:", line), "info")
                          })

  if (!is.null(log_callback)) {
    if (identical(result$status, 0L) || identical(result$status, 0)) {
      log_callback("Salmon index: completed", "success")
    } else {
      log_callback(paste("Salmon index: failed with exit code", result$status), "error")
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
                             validate = TRUE, bootstraps = 0,
                             min_score_frac = 0.65, discard_orphans = FALSE,
                             log_callback = NULL) {

  sample_out <- file.path(outdir, sample_name)
  dir.create(sample_out, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(log_callback)) log_callback(paste("Salmon quant:", sample_name, "— running..."), "info")

  args <- c("quant",
            "-i", index_dir,
            "-l", lib_type,
            "-p", as.character(threads),
            "--minScoreFraction", as.character(min_score_frac),
            "-o", sample_out)

  if (is_se) {
    args <- c(args, "-r", r1)
  } else {
    args <- c(args, "-1", r1, "-2", r2)
  }

  if (validate)        args <- c(args, "--validateMappings")
  if (gc_bias)         args <- c(args, "--gcBias")
  if (seq_bias)        args <- c(args, "--seqBias")
  if (bootstraps > 0)  args <- c(args, "--numBootstraps", as.character(bootstraps))
  if (discard_orphans) args <- c(args, "--discardOrphansQuasi")

  result <- processx::run("salmon", args = args, echo = FALSE,
                          error_on_status = FALSE,
                          stderr_line_callback = function(line, proc) {
                            if (!is.null(log_callback)) log_callback(paste("Salmon quant:", line), "info")
                          })

  meta <- parse_salmon_meta(sample_out)

  if (!is.null(log_callback)) {
    if (identical(result$status, 0L) || identical(result$status, 0)) {
      rate_msg <- if (!is.na(meta$percent_mapped)) paste0("mapping ", meta$percent_mapped, "%") else ""
      log_callback(paste("Salmon quant:", sample_name, "done", rate_msg), "success")
      if (!is.na(meta$percent_mapped) && meta$percent_mapped < 50) {
        log_callback(paste("Warning:", sample_name, "mapping rate below 50%!"), "warn")
      }
    } else {
      log_callback(paste("Salmon quant:", sample_name, "failed with exit code", result$status), "error")
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

  result <- processx::run("multiqc", args = args, echo = FALSE,
                          error_on_status = FALSE,
                          stderr_line_callback = function(line, proc) {
                            if (!is.null(log_callback)) log_callback(paste("MultiQC:", line), "info")
                          })

  report_path <- file.path(outdir, "multiqc_report.html")

  if (!is.null(log_callback)) {
    if (file.exists(report_path)) {
      log_callback("MultiQC: report generated", "success")
    } else {
      log_callback("MultiQC: report generation may have failed", "warn")
    }
  }

  list(exit_status = result$status %||% 1, report_path = report_path)
}
