#!/usr/bin/env Rscript
# SalmonFlow — background pipeline runner
# Args: <params_json> <log_file> <state_file>

# When Shiny auto-sources R/ files on startup, sys.nframe() > 0 — bail out silently.
# This script is only meant to be executed directly via Rscript.
if (sys.nframe() > 0) return(invisible(NULL))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: run_pipeline.R <params_json> <log_file> <state_file>")

params_file <- args[1]
log_file    <- args[2]
state_file  <- args[3]

# Locate this script's directory to source sibling files
script_args <- commandArgs(trailingOnly = FALSE)
file_arg    <- grep("^--file=", script_args, value = TRUE)
script_dir  <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}

suppressPackageStartupMessages({
  library(processx)
  library(jsonlite)
})

source(file.path(script_dir, "helpers.R"))
source(file.path(script_dir, "pipeline_functions.R"))
source(file.path(script_dir, "tximport_utils.R"))

# ── Load params ───────────────────────────────────────────────
p <- jsonlite::fromJSON(params_file, simplifyDataFrame = TRUE)

resume          <- isTRUE(p$resume)
run_start       <- Sys.time()
delete_originals <- isTRUE(p$delete_originals_after_trim)
delete_trimmed  <- isTRUE(p$delete_trimmed_after_quant)

# ── IPC helpers ───────────────────────────────────────────────
samples       <- as.data.frame(p$samples, stringsAsFactors = FALSE)
sample_names  <- samples$name
sample_status <- setNames(rep("pending", length(sample_names)), sample_names)

write_log <- function(msg, type = "info") {
  cat(paste0(type, "|", msg, "\n"), file = log_file, append = TRUE)
}

write_state <- function(step, total, running, ...) {
  state <- c(
    list(step = step, total = total, running = running,
         sample_status = as.list(sample_status)),
    list(...)
  )
  tmp <- paste0(state_file, ".tmp")
  writeLines(jsonlite::toJSON(state, auto_unbox = TRUE), tmp)
  file.rename(tmp, state_file)
}

# Delete an original FASTQ once its trimmed counterpart is confirmed valid.
# Guarded: only runs when the user opted in, the trimmed file exists and is
# non-empty, and the original is not the same path as the trimmed output.
remove_original <- function(orig_path, trimmed_path) {
  if (!delete_originals) return(invisible(NULL))
  if (is.null(orig_path) || length(orig_path) == 0 ||
      is.na(orig_path) || nchar(orig_path) == 0) return(invisible(NULL))
  if (!file.exists(orig_path)) return(invisible(NULL))

  # Trimmed output must exist and be non-empty before we touch the source
  if (is.null(trimmed_path) || length(trimmed_path) == 0 ||
      is.na(trimmed_path) || !file.exists(trimmed_path) ||
      isTRUE(file.info(trimmed_path)$size <= 0)) {
    write_log(paste("Keeping raw (invalid trimmed output):",
                    basename(orig_path)), "warn")
    return(invisible(NULL))
  }

  # Never delete the trimmed file itself
  if (identical(normalizePath(orig_path,    mustWork = FALSE),
                normalizePath(trimmed_path, mustWork = FALSE))) {
    return(invisible(NULL))
  }

  ok <- tryCatch(file.remove(orig_path), error = function(e) FALSE)
  if (isTRUE(ok)) {
    write_log(paste("Raw FASTQ deleted to free space:",
                    basename(orig_path)), "warn")
  } else {
    write_log(paste("Could not delete raw FASTQ:",
                    basename(orig_path)), "error")
  }
  invisible(NULL)
}

# Delete a trimmed FASTQ once its Salmon quant.sf is confirmed valid.
# Guarded: only runs when the user opted in and quant.sf exists and is
# non-empty. Safe to call after FastQC-post (step 2b) already consumed
# the trimmed files, since quant (step 4) always runs after it.
remove_trimmed <- function(trimmed_path, quant_sf_path) {
  if (!delete_trimmed) return(invisible(NULL))
  if (is.null(trimmed_path) || length(trimmed_path) == 0 ||
      is.na(trimmed_path) || nchar(trimmed_path) == 0) return(invisible(NULL))
  if (!file.exists(trimmed_path)) return(invisible(NULL))

  if (is.null(quant_sf_path) || length(quant_sf_path) == 0 ||
      is.na(quant_sf_path) || !file.exists(quant_sf_path) ||
      isTRUE(file.info(quant_sf_path)$size <= 0)) {
    write_log(paste("Keeping trimmed (invalid quant.sf):",
                    basename(trimmed_path)), "warn")
    return(invisible(NULL))
  }

  ok <- tryCatch(file.remove(trimmed_path), error = function(e) FALSE)
  if (isTRUE(ok)) {
    write_log(paste("Trimmed FASTQ deleted to free space:",
                    basename(trimmed_path)), "warn")
  } else {
    write_log(paste("Could not delete trimmed FASTQ:",
                    basename(trimmed_path)), "error")
  }
  invisible(NULL)
}

# ── Derive pipeline config ────────────────────────────────────
is_se       <- isTRUE(p$lib_type_se)
mode        <- if (is_se) "SE" else "PE"
threads     <- as.integer(p$salmon_threads %||% 4L)
output_dir  <- as.character(p$output_dir)
trim_dir    <- file.path(output_dir, "trimmed")
fastqc_dir  <- file.path(output_dir, "fastqc_pre")
fastqc_post_dir <- file.path(output_dir, "fastqc_post")
quant_dir   <- file.path(output_dir, "salmon_quant")
multiqc_dir     <- file.path(output_dir, "multiqc")
multiqc_pre_dir <- file.path(output_dir, "multiqc_pre")
multiqc_post_dir <- file.path(output_dir, "multiqc_post")

adapter_fasta <- {
  af <- p$adapter_fasta %||% ""
  if (nchar(af) > 0) af else NULL
}

index_dir <- if (isTRUE(p$build_new_index)) {
  file.path("/data/references", "salmon_index")
} else {
  af <- p$salmon_index_dir %||% ""
  if (nchar(af) > 0) af else stop("No salmon index directory specified")
}

n_samples <- nrow(samples)

# ── Calculate total steps ─────────────────────────────────────
total <- 1L  # Step 1: FastQC (pre-trimming)
if (isTRUE(p$build_new_index))  total <- total + 1L  # Step 2: Salmon index
total <- total + n_samples  # Step 3: Per-sample loop (trim + FastQC post + quant)
total <- total + 1L  # Step 4: tximport
if (isTRUE(p$trimming_enabled)) {
  total <- total + 2L  # Step 5: MultiQC pre and post
} else {
  total <- total + 1L  # Step 5: MultiQC
}
step <- 0L

write_log("=== SalmonFlow Pipeline started ===", "info")
write_log(paste("Samples:", n_samples, "| Mode:", mode), "info")
if (resume) write_log("RESUME mode active — steps with previous results will be skipped", "info")
write_state(step, total, TRUE)

# ── Helper: FastQC output stem ────────────────────────────────
fastqc_stem <- function(f) sub("\\.(fastq|fq)(\\.gz)?$", "", basename(f), ignore.case = TRUE)

# ── STEP 1: FastQC (pre-trimming) ─────────────────────────────
write_log("-- Step 1: FastQC (pre-trimming) --", "info")
all_fastqs <- samples$r1
if (!is_se) all_fastqs <- c(all_fastqs, samples$r2)
all_fastqs <- all_fastqs[!is.na(all_fastqs) & nchar(all_fastqs) > 0]

fastqc_pre_done <- resume && all(sapply(all_fastqs, function(f) {
  file.exists(file.path(fastqc_dir, paste0(fastqc_stem(f), "_fastqc.zip")))
}))

if (fastqc_pre_done) {
  write_log("FastQC (pre-trimming): skipped (previous results found)", "info")
} else {
  run_fastqc(all_fastqs, fastqc_dir, threads = threads, log_callback = write_log)
}
step <- step + 1L
write_state(step, total, TRUE)

# ── STEP 2: Salmon index ─────────────────────────────────────
if (isTRUE(p$build_new_index)) {
  write_log("-- Step 2: Build Salmon Index --", "info")

  index_done <- resume && file.exists(file.path(index_dir, "info.json"))
  if (index_done) {
    write_log("Salmon index: skipped (previous index found)", "info")
  } else {
    decoy_file <- if (isTRUE(p$decoy_aware)) {
      gf <- p$genome_fasta %||% ""; if (nchar(gf) > 0) gf else NULL
    } else NULL

    idx_result <- build_salmon_index(
      fasta        = as.character(p$transcriptome_fasta),
      outdir       = index_dir,
      decoy        = decoy_file,
      kmer         = as.integer(p$kmer_size %||% 31L),
      threads      = threads,
      sparse       = isTRUE(p$sparse_index),
      log_callback = write_log
    )

    if (idx_result$exit_status != 0) {
      write_log("Pipeline aborted: error building index", "error")
      write_state(step, total, FALSE)
      quit(status = 1, save = "no")
    }
    index_dir <- idx_result$index_dir
  }
  step <- step + 1L
  write_state(step, total, TRUE)
} else {
  write_log("-- Step 2: Salmon Index (using existing) --", "info")
}

# ── STEP 3: Merged Per-Sample Loop (trim + FastQC post + quant) ──
write_log("-- Step 3: Per-Sample Processing (trimming, FastQC post, Salmon quant) --", "info")
salmon_metas <- list()

for (i in seq_len(n_samples)) {
  sname <- sample_names[i]
  sample_status[sname] <- "running"
  write_state(step, total, TRUE)

  # Define output paths for this sample
  r1_trimmed <- file.path(trim_dir, paste0(sname, if (is_se) "_trimmed.fastq.gz" else "_R1_trimmed.fastq.gz"))
  r2_trimmed <- if (is_se) NULL else file.path(trim_dir, paste0(sname, "_R2_trimmed.fastq.gz"))
  quant_sf <- file.path(quant_dir, sname, "quant.sf")

  # 1. Check if we can resume/skip the entire sample (quant already done)
  if (resume && file.exists(quant_sf) && file.info(quant_sf)$size > 0) {
    sample_status[sname] <- "done"
    salmon_metas[[sname]] <- parse_salmon_meta(file.path(quant_dir, sname))
    write_log(paste("Salmon quant:", sname, "— skipped (previous quant.sf found)"), "info")

    # Even if skipped, try to remove original/trimmed files if requested
    if (isTRUE(p$trimming_enabled)) {
      remove_original(samples$r1[i], r1_trimmed)
      if (!is_se) remove_original(samples$r2[i], r2_trimmed)

      remove_trimmed(r1_trimmed, quant_sf)
      if (!is_se) remove_trimmed(r2_trimmed, quant_sf)
    }

    step <- step + 1L
    write_log(paste("  Sample completed (resume):", i, "/", n_samples), "info")
    write_state(step, total, TRUE)
    next
  }

  # 2. Trimming (fastp) if enabled
  if (isTRUE(p$trimming_enabled)) {
    fastp_skipped <- FALSE
    if (resume && file.exists(r1_trimmed) && (is_se || file.exists(r2_trimmed))) {
      write_log(paste("fastp:", sname, "— skipped (previous output found)"), "info")
      remove_original(samples$r1[i], r1_trimmed)
      if (!is_se) remove_original(samples$r2[i], r2_trimmed)
      fastp_skipped <- TRUE
    } else {
      r2_val <- if (is_se) NULL else {
        v <- samples$r2[i]; if (is.na(v) || nchar(v) == 0) NULL else v
      }

      trim_result <- run_fastp(
        r1                = samples$r1[i],
        r2                = r2_val,
        out_dir           = trim_dir,
        sample_name       = sname,
        mode              = mode,
        adapter_fasta     = adapter_fasta,
        cut_front_quality = as.numeric(p$fastp_cut_front   %||% 20),
        cut_tail_quality  = as.numeric(p$fastp_cut_tail    %||% 20),
        cut_right_quality = as.numeric(p$fastp_cut_right   %||% 20),
        window_size       = as.integer(p$fastp_window_size %||% 4L),
        minlen            = as.integer(p$fastp_minlen      %||% 36L),
        threads           = threads,
        log_callback      = write_log
      )

      if (trim_result$exit_status != 0) {
        sample_status[sname] <- "error"
        step <- step + 1L
        write_state(step, total, TRUE)
        next
      }

      # Successfully trimmed, remove originals if requested
      remove_original(samples$r1[i], r1_trimmed)
      if (!is_se) remove_original(r2_val, r2_trimmed)
    }

    # 3. Post-trim FastQC on this sample's trimmed files
    fq_trimmed <- r1_trimmed
    if (!is_se) fq_trimmed <- c(fq_trimmed, r2_trimmed)
    fq_trimmed <- fq_trimmed[!is.na(fq_trimmed) & nchar(fq_trimmed) > 0]

    fastqc_post_sample_done <- resume && all(sapply(fq_trimmed, function(f) {
      file.exists(file.path(fastqc_post_dir, paste0(fastqc_stem(f), "_fastqc.zip")))
    }))

    if (fastqc_post_sample_done) {
      write_log(paste("FastQC (post-trimming):", sname, "— skipped (previous results found)"), "info")
    } else {
      run_fastqc(fq_trimmed, fastqc_post_dir, threads = threads, log_callback = write_log)
    }
  }

  # 4. Salmon Quantification
  r1_input <- if (isTRUE(p$trimming_enabled)) r1_trimmed else samples$r1[i]
  r2_input <- if (isTRUE(p$trimming_enabled)) {
    r2_trimmed
  } else {
    if (is_se) NULL else {
      v <- samples$r2[i]; if (is.na(v) || nchar(v) == 0) NULL else v
    }
  }

  quant_result <- run_salmon_quant(
    index_dir       = index_dir,
    r1              = r1_input,
    r2              = r2_input,
    outdir          = quant_dir,
    sample_name     = sname,
    lib_type        = as.character(p$salmon_libtype  %||% "A"),
    gc_bias         = isTRUE(p$salmon_gcbias),
    seq_bias        = isTRUE(p$salmon_seqbias),
    threads         = threads,
    is_se           = is_se,
    validate        = isTRUE(p$salmon_validate),
    bootstraps      = as.integer(p$salmon_bootstraps    %||% 0L),
    min_score_frac  = as.numeric(p$salmon_min_score_frac %||% 0.65),
    discard_orphans = isTRUE(p$salmon_discard_orphans),
    log_callback    = write_log
  )

  if (quant_result$exit_status == 0) {
    sample_status[sname] <- "done"
    salmon_metas[[sname]] <- quant_result$meta

    # Optionally free disk by removing the trimmed FASTQs once quant succeeds
    if (isTRUE(p$trimming_enabled)) {
      remove_trimmed(r1_trimmed, quant_sf)
      if (!is_se) remove_trimmed(r2_trimmed, quant_sf)
    }
  } else {
    sample_status[sname] <- "error"
  }

  step <- step + 1L
  write_log(paste("  Sample completed:", i, "/", n_samples), "info")
  write_state(step, total, TRUE)
}

salmon_meta_path <- paste0(state_file, "_salmon_meta.json")
writeLines(jsonlite::toJSON(salmon_metas, auto_unbox = TRUE), salmon_meta_path)

# ── STEP 4: tximport ─────────────────────────────────────────
write_log("-- Step 4: tximport --", "info")

tx2gene <- build_tx2gene(as.character(p$gtf_path), log_callback = write_log)
if (is.null(tx2gene)) {
  write_log("Pipeline aborted: error building tx2gene", "error")
  write_state(step, total, FALSE)
  quit(status = 1, save = "no")
}

run_tximport(
  quant_dir         = quant_dir,
  sample_names      = sample_names,
  tx2gene           = tx2gene,
  method            = as.character(p$txi_method         %||% "lengthScaledTPM"),
  ignore_tx_version = isTRUE(p$txi_ignore_version),
  output_dir        = output_dir,
  log_callback      = write_log
)

step <- step + 1L
count_matrix_path <- file.path(output_dir, "merged_lengthScaledTPM.csv")
write_state(step, total, TRUE)

# ── Run summary log ──────────────────────────────────────────
# Human-readable recap of the run written to the output directory as a
# durable artifact (the live IPC log in /data/tmp is machine-formatted and
# transient). Returns the path, or "" if it could not be written.
summary_path <- file.path(output_dir, "run_summary.log")

write_run_summary <- function(multiqc = list()) {
  fmt_int <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return("N/A")
    format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
  }
  yn <- function(x) if (isTRUE(x)) "yes" else "no"

  dur_secs <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
  dur_str  <- sprintf("%dm %ds", dur_secs %/% 60, round(dur_secs %% 60))

  n_done  <- sum(sample_status == "done")
  n_error <- sum(sample_status == "error")

  L <- c(
    "=== SalmonFlow Run Summary ===",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Duration: ", dur_str),
    "",
    "-- Configuration --",
    paste0("Library type: ", if (is_se) "Single-end (SE)" else "Paired-end (PE)"),
    paste0("Trimming (fastp): ", if (isTRUE(p$trimming_enabled)) "ON" else "OFF")
  )
  if (isTRUE(p$trimming_enabled)) {
    L <- c(L, sprintf(
      "  cut_front Q: %s | cut_tail Q: %s | cut_right Q: %s | window: %s | min length: %s",
      p$fastp_cut_front %||% "20", p$fastp_cut_tail %||% "20",
      p$fastp_cut_right %||% "20", p$fastp_window_size %||% "4",
      p$fastp_minlen %||% "36"))
  }
  L <- c(L,
    paste0("Salmon library type: ", p$salmon_libtype %||% "A"),
    paste0("GC bias correction: ", yn(p$salmon_gcbias)),
    paste0("Sequence bias correction: ", yn(p$salmon_seqbias)),
    paste0("Bootstraps: ", p$salmon_bootstraps %||% 0L),
    paste0("Salmon index: ", index_dir,
           if (isTRUE(p$build_new_index)) " (newly built)" else " (existing)"),
    if (isTRUE(p$build_new_index)) paste0("Decoy-aware index: ", yn(p$decoy_aware)) else NULL,
    paste0("Threads: ", threads),
    paste0("tximport method: ", p$txi_method %||% "lengthScaledTPM",
           " (ignore tx version: ", yn(p$txi_ignore_version), ")"),
    "",
    sprintf("-- Samples (%d total) --", n_samples),
    sprintf("Completed: %d    Failed: %d", n_done, n_error)
  )

  for (s in sample_names) {
    st   <- sample_status[[s]] %||% "unknown"
    meta <- salmon_metas[[s]]
    if (identical(st, "done") && !is.null(meta)) {
      rate <- if (!is.null(meta$percent_mapped) && !is.na(meta$percent_mapped)) {
        paste0(meta$percent_mapped, "%")
      } else "N/A"
      L <- c(L, sprintf("  %-20s %-7s reads: %s   mapped: %s (%s)",
                        s, st, fmt_int(meta$num_processed),
                        fmt_int(meta$num_mapped), rate))
    } else {
      L <- c(L, sprintf("  %-20s %-7s —", s, st))
    }
  }

  L <- c(L, "", "-- Outputs --",
         paste0("Count matrix: ", count_matrix_path))
  for (nm in names(multiqc)) {
    if (nchar(multiqc[[nm]]) > 0) L <- c(L, paste0(nm, ": ", multiqc[[nm]]))
  }
  L <- c(L, paste0("Run summary:  ", summary_path))

  ok <- tryCatch({
    writeLines(L, summary_path); TRUE
  }, error = function(e) {
    write_log(paste("Could not write run summary:", conditionMessage(e)), "warn")
    FALSE
  })
  if (isTRUE(ok)) write_log(paste("Run summary written:", summary_path), "info")
  if (isTRUE(ok)) summary_path else ""
}

# ── STEP 5: MultiQC ──────────────────────────────────────────
if (isTRUE(p$trimming_enabled)) {
  write_log("-- Step 5a: MultiQC (pre-trimming) --", "info")
  run_multiqc(c(fastqc_dir, trim_dir), multiqc_pre_dir, log_callback = write_log)
  step <- step + 1L
  write_state(step, total, TRUE)

  write_log("-- Step 5b: MultiQC (post-trimming) --", "info")
  run_multiqc(c(fastqc_post_dir, quant_dir), multiqc_post_dir, log_callback = write_log)
  step <- step + 1L

  multiqc_pre_report_path  <- file.path(multiqc_pre_dir,  "multiqc_report.html")
  multiqc_post_report_path <- file.path(multiqc_post_dir, "multiqc_report.html")
  if (!file.exists(multiqc_pre_report_path))  multiqc_pre_report_path  <- ""
  if (!file.exists(multiqc_post_report_path)) multiqc_post_report_path <- ""

  summary_out <- write_run_summary(list(
    "MultiQC (pre-trimming)"  = multiqc_pre_report_path,
    "MultiQC (post-trimming)" = multiqc_post_report_path))

  write_log("=== Pipeline completed successfully ===", "success")
  write_state(step, total, FALSE,
              count_matrix_path        = count_matrix_path,
              multiqc_pre_report_path  = multiqc_pre_report_path,
              multiqc_post_report_path = multiqc_post_report_path,
              salmon_meta_path         = salmon_meta_path,
              summary_path             = summary_out)
} else {
  write_log("-- Step 5: MultiQC --", "info")
  run_multiqc(output_dir, multiqc_dir, log_callback = write_log)
  step <- step + 1L

  multiqc_report_path <- file.path(multiqc_dir, "multiqc_report.html")
  if (!file.exists(multiqc_report_path)) multiqc_report_path <- ""

  summary_out <- write_run_summary(list("MultiQC report" = multiqc_report_path))

  write_log("=== Pipeline completed successfully ===", "success")
  write_state(step, total, FALSE,
              count_matrix_path   = count_matrix_path,
              multiqc_report_path = multiqc_report_path,
              salmon_meta_path    = salmon_meta_path,
              summary_path        = summary_out)
}
