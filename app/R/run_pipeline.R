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

resume <- isTRUE(p$resume)

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
total <- 1L
if (isTRUE(p$trimming_enabled)) total <- total + n_samples + 1L + 1L  # fastp + post-FastQC + extra MultiQC
if (isTRUE(p$build_new_index))  total <- total + 1L
total <- total + n_samples + 2L  # quant per sample + tximport + MultiQC
step <- 0L

write_log("=== Pipeline SalmonFlow iniciado ===", "info")
write_log(paste("Muestras:", n_samples, "| Modo:", mode), "info")
if (resume) write_log("Modo REANUDAR activo — pasos con resultados previos seran omitidos", "info")
write_state(step, total, TRUE)

# ── Helper: FastQC output stem ────────────────────────────────
fastqc_stem <- function(f) sub("\\.(fastq|fq)(\\.gz)?$", "", basename(f), ignore.case = TRUE)

# ── STEP 1: FastQC (pre-trimming) ─────────────────────────────
write_log("-- Paso 1: FastQC (pre-trimming) --", "info")
all_fastqs <- samples$r1
if (!is_se) all_fastqs <- c(all_fastqs, samples$r2)
all_fastqs <- all_fastqs[!is.na(all_fastqs) & nchar(all_fastqs) > 0]

fastqc_pre_done <- resume && all(sapply(all_fastqs, function(f) {
  file.exists(file.path(fastqc_dir, paste0(fastqc_stem(f), "_fastqc.zip")))
}))

if (fastqc_pre_done) {
  write_log("FastQC (pre-trimming): omitido (resultados previos encontrados)", "info")
} else {
  run_fastqc(all_fastqs, fastqc_dir, threads = threads, log_callback = write_log)
}
step <- step + 1L
write_state(step, total, TRUE)

# ── STEP 2: fastp ────────────────────────────────────────────
trimmed_samples <- samples
if (isTRUE(p$trimming_enabled)) {
  write_log("-- Paso 2: fastp --", "info")

  for (i in seq_len(n_samples)) {
    sname <- sample_names[i]
    sample_status[sname] <- "running"
    write_state(step, total, TRUE)

    r1_out <- file.path(trim_dir, paste0(sname, "_R1_trimmed.fastq.gz"))
    r2_out <- if (is_se) NULL else file.path(trim_dir, paste0(sname, "_R2_trimmed.fastq.gz"))

    if (resume && file.exists(r1_out) && (is_se || file.exists(r2_out))) {
      trimmed_samples$r1[i] <- r1_out
      if (!is_se) trimmed_samples$r2[i] <- r2_out
      sample_status[sname] <- "done"
      write_log(paste("fastp:", sname, "— omitido (salida previa encontrada)"), "info")
      step <- step + 1L
      write_state(step, total, TRUE)
      next
    }

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
      cut_front_quality = as.numeric(p$fastp_cut_front %||% 3),
      cut_tail_quality  = as.numeric(p$fastp_cut_tail  %||% 3),
      cut_right_quality = as.numeric(p$fastp_cut_right %||% 15),
      minlen            = as.integer(p$fastp_minlen    %||% 36L),
      threads           = threads,
      log_callback      = write_log
    )

    if (trim_result$exit_status == 0) {
      trimmed_samples$r1[i] <- trim_result$r1_trimmed
      if (!is_se) trimmed_samples$r2[i] <- trim_result$r2_trimmed
      sample_status[sname] <- "done"
    } else {
      sample_status[sname] <- "error"
    }

    step <- step + 1L
    write_state(step, total, TRUE)
  }

  # Reset to pending for Salmon quant display
  for (s in sample_names) {
    if (sample_status[s] == "done") sample_status[s] <- "pending"
  }

  # ── STEP 2b: FastQC (post-trimming) ───────────────────────
  write_log("-- Paso 2b: FastQC (post-trimming) --", "info")
  trimmed_fastqs <- trimmed_samples$r1
  if (!is_se) trimmed_fastqs <- c(trimmed_fastqs, trimmed_samples$r2)
  trimmed_fastqs <- trimmed_fastqs[!is.na(trimmed_fastqs) & nchar(trimmed_fastqs) > 0]

  fastqc_post_done <- resume && all(sapply(trimmed_fastqs, function(f) {
    file.exists(file.path(fastqc_post_dir, paste0(fastqc_stem(f), "_fastqc.zip")))
  }))

  if (fastqc_post_done) {
    write_log("FastQC (post-trimming): omitido (resultados previos encontrados)", "info")
  } else {
    run_fastqc(trimmed_fastqs, fastqc_post_dir, threads = threads, log_callback = write_log)
  }
  step <- step + 1L
  write_state(step, total, TRUE)
} else {
  write_log("-- Paso 2: fastp (omitido) --", "info")
}

# ── STEP 3: Salmon index ─────────────────────────────────────
if (isTRUE(p$build_new_index)) {
  write_log("-- Paso 3: Construccion del indice Salmon --", "info")

  index_done <- resume && file.exists(file.path(index_dir, "info.json"))
  if (index_done) {
    write_log("Salmon index: omitido (indice previo encontrado)", "info")
    step <- step + 1L
    write_state(step, total, TRUE)
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
    step <- step + 1L
    write_state(step, total, TRUE)

    if (idx_result$exit_status != 0) {
      write_log("Pipeline abortado: error al construir el indice", "error")
      write_state(step, total, FALSE)
      quit(status = 1, save = "no")
    }
    index_dir <- idx_result$index_dir
  }
} else {
  write_log("-- Paso 3: Indice Salmon (usando existente) --", "info")
}

# ── STEP 4: Salmon quant ─────────────────────────────────────
write_log("-- Paso 4: Cuantificacion Salmon --", "info")
salmon_metas <- list()

for (i in seq_len(n_samples)) {
  sname <- sample_names[i]
  sample_status[sname] <- "running"
  write_state(step, total, TRUE)

  quant_sf <- file.path(quant_dir, sname, "quant.sf")

  if (resume && file.exists(quant_sf)) {
    sample_status[sname] <- "done"
    salmon_metas[[sname]] <- parse_salmon_meta(file.path(quant_dir, sname))
    write_log(paste("Salmon quant:", sname, "— omitido (quant.sf previo encontrado)"), "info")
    step <- step + 1L
    write_log(paste("  Salmon quant:", i, "/", n_samples), "info")
    write_state(step, total, TRUE)
    next
  }

  r2_val <- if (is_se) NULL else {
    v <- trimmed_samples$r2[i]; if (is.na(v) || nchar(v) == 0) NULL else v
  }

  quant_result <- run_salmon_quant(
    index_dir       = index_dir,
    r1              = trimmed_samples$r1[i],
    r2              = r2_val,
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
  } else {
    sample_status[sname] <- "error"
  }

  step <- step + 1L
  write_log(paste("  Salmon quant:", i, "/", n_samples), "info")
  write_state(step, total, TRUE)
}

salmon_meta_path <- paste0(state_file, "_salmon_meta.json")
writeLines(jsonlite::toJSON(salmon_metas, auto_unbox = TRUE), salmon_meta_path)

# ── STEP 5: tximport ─────────────────────────────────────────
write_log("-- Paso 5: tximport --", "info")

tx2gene <- build_tx2gene(as.character(p$gtf_path), log_callback = write_log)
if (is.null(tx2gene)) {
  write_log("Pipeline abortado: error al construir tx2gene", "error")
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

# ── STEP 6: MultiQC ──────────────────────────────────────────
if (isTRUE(p$trimming_enabled)) {
  write_log("-- Paso 6a: MultiQC (pre-trimming) --", "info")
  run_multiqc(c(fastqc_dir, trim_dir), multiqc_pre_dir, log_callback = write_log)
  step <- step + 1L
  write_state(step, total, TRUE)

  write_log("-- Paso 6b: MultiQC (post-trimming) --", "info")
  run_multiqc(c(fastqc_post_dir, quant_dir), multiqc_post_dir, log_callback = write_log)
  step <- step + 1L

  multiqc_pre_report_path  <- file.path(multiqc_pre_dir,  "multiqc_report.html")
  multiqc_post_report_path <- file.path(multiqc_post_dir, "multiqc_report.html")
  if (!file.exists(multiqc_pre_report_path))  multiqc_pre_report_path  <- ""
  if (!file.exists(multiqc_post_report_path)) multiqc_post_report_path <- ""

  write_log("=== Pipeline completado exitosamente ===", "success")
  write_state(step, total, FALSE,
              count_matrix_path        = count_matrix_path,
              multiqc_pre_report_path  = multiqc_pre_report_path,
              multiqc_post_report_path = multiqc_post_report_path,
              salmon_meta_path         = salmon_meta_path)
} else {
  write_log("-- Paso 6: MultiQC --", "info")
  run_multiqc(output_dir, multiqc_dir, log_callback = write_log)
  step <- step + 1L

  multiqc_report_path <- file.path(multiqc_dir, "multiqc_report.html")
  if (!file.exists(multiqc_report_path)) multiqc_report_path <- ""

  write_log("=== Pipeline completado exitosamente ===", "success")
  write_state(step, total, FALSE,
              count_matrix_path   = count_matrix_path,
              multiqc_report_path = multiqc_report_path,
              salmon_meta_path    = salmon_meta_path)
}
