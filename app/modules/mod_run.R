# ══════════════════════════════════════════════════════════════
# SalmonFlow — mod_run.R
# Tab 4: Pipeline execution with live logs and progress
# ══════════════════════════════════════════════════════════════

mod_run_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        box(
          title = "Resumen de configuración",
          status = "primary", solidHeader = FALSE, width = 12,
          collapsible = TRUE,
          uiOutput(ns("config_summary"))
        )
      )
    ),

    fluidRow(
      column(12,
        box(
          title = "Ejecución del Pipeline",
          status = "primary", solidHeader = FALSE, width = 12,

          fluidRow(
            column(4,
              actionButton(ns("run_btn"), "Iniciar análisis",
                           class = "btn-primary", icon = icon("play"),
                           style = "width:100%; font-size:16px; padding:14px;")
            ),
            column(4,
              actionButton(ns("cancel_btn"), "Cancelar",
                           class = "btn-danger", icon = icon("stop"),
                           style = "width:100%; font-size:16px; padding:14px;")
            ),
            column(4,
              tags$div(style = "padding-top:8px;",
                uiOutput(ns("pipeline_status_badge"))
              )
            )
          ),

          br(),

          tags$div(
            tags$strong("Progreso general:"),
            tags$div(class = "progress",
              tags$div(id = ns("progress_bar"),
                       class = "progress-bar",
                       role = "progressbar",
                       style = "width: 0%",
                       "0%")
            )
          ),

          br(),

          tags$strong("Progreso por muestra:"),
          uiOutput(ns("sample_status_list")),

          br(),

          tags$strong("Log en vivo:"),
          uiOutput(ns("live_log"))
        )
      )
    )
  )
}

mod_run_server <- function(id, shared) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Config summary ─────────────────────────────────────
    output$config_summary <- renderUI({
      samples <- shared$samples
      n_samples <- if (is.null(samples)) 0 else nrow(samples)
      lib <- if (shared$lib_type_se) "Single-end" else "Paired-end"
      trim <- if (shared$trimming_enabled) "ON" else "OFF"

      tags$div(class = "config-summary",
        tags$p(tags$strong("Muestras: "), n_samples, " | ",
               tags$strong("Tipo: "), lib, " | ",
               tags$strong("Trimming: "), trim),
        tags$p(tags$strong("Salmon library: "), shared$salmon_libtype, " | ",
               tags$strong("GC bias: "), ifelse(shared$salmon_gcbias, "ON", "OFF"), " | ",
               tags$strong("Seq bias: "), ifelse(shared$salmon_seqbias, "ON", "OFF"), " | ",
               tags$strong("Threads: "), shared$salmon_threads),
        tags$p(tags$strong("tximport: "), shared$txi_method, " | ",
               tags$strong("Ignorar versión: "), ifelse(shared$txi_ignore_version, "Sí", "No"))
      )
    })

    # ── Local pipeline state ───────────────────────────────
    rv <- reactiveValues(
      log_lines      = character(0),
      running        = FALSE,
      current_step   = 0,
      total_steps    = 6,
      sample_status  = list(),   # named list: sample_name -> "pending"|"running"|"done"|"error"
      cancel_flag    = FALSE
    )

    # Log helper
    add_log <- function(msg, type = "info") {
      rv$log_lines <- c(rv$log_lines, timestamp_log(msg, type))
    }

    # ── Live log output ────────────────────────────────────
    output$live_log <- renderUI({
      lines <- rv$log_lines
      if (length(lines) == 0) {
        tags$div(class = "log-panel",
                 tags$span(class = "log-info", "Esperando inicio del análisis..."))
      } else {
        tags$div(class = "log-panel",
                 HTML(paste(lines, collapse = "<br/>")))
      }
    })

    # ── Progress bar update ────────────────────────────────
    observe({
      pct <- if (rv$total_steps > 0) round(rv$current_step / rv$total_steps * 100) else 0
      shinyjs::runjs(sprintf(
        "document.getElementById('%s').style.width = '%d%%';
         document.getElementById('%s').textContent = '%d%%';",
        ns("progress_bar"), pct, ns("progress_bar"), pct
      ))
    })

    # ── Pipeline status badge ──────────────────────────────
    output$pipeline_status_badge <- renderUI({
      if (rv$running) {
        tags$span(class = "badge-running", "Pipeline en ejecución")
      } else if (rv$current_step >= rv$total_steps && rv$current_step > 0) {
        tags$span(class = "badge-success", "Pipeline completado")
      } else {
        tags$span(style = "color:#5a6373;", "Listo para ejecutar")
      }
    })

    # ── Per-sample status ──────────────────────────────────
    output$sample_status_list <- renderUI({
      statuses <- rv$sample_status
      if (length(statuses) == 0) return(tags$p(style="color:#5a6373;", "—"))

      tags$div(
        lapply(names(statuses), function(s) {
          st <- statuses[[s]]
          badge_class <- switch(st,
            done    = "badge-success",
            running = "badge-running",
            error   = "badge-error",
            "badge-default"
          )
          tags$span(style = "margin-right:12px;",
            tags$span(class = badge_class, s)
          )
        })
      )
    })

    # ── Cancel button ──────────────────────────────────────
    observeEvent(input$cancel_btn, {
      if (rv$running) {
        rv$cancel_flag <- TRUE
        add_log("Pipeline cancelado por el usuario", "warn")
        rv$running <- FALSE
        shared$pipeline_running <- FALSE
      }
    })

    # ── Main pipeline execution ────────────────────────────
    observeEvent(input$run_btn, {
      # Validate prerequisites
      samples <- shared$samples
      if (is.null(samples) || nrow(samples) == 0) {
        add_log("Error: no hay muestras cargadas. Ve a la pestaña Muestras.", "error")
        return()
      }

      if (shared$build_new_index && (is.null(shared$transcriptome_fasta) || length(shared$transcriptome_fasta) == 0)) {
        add_log("Error: no se ha seleccionado un transcriptoma FASTA. Ve a la pestaña Referencias.", "error")
        return()
      }

      if (is.null(shared$gtf_path) || length(shared$gtf_path) == 0) {
        add_log("Error: no se ha seleccionado un archivo GTF. Ve a la pestaña Referencias.", "error")
        return()
      }

      # Reset state
      rv$log_lines     <- character(0)
      rv$current_step  <- 0
      rv$cancel_flag   <- FALSE
      rv$running       <- TRUE
      shared$pipeline_running <- TRUE

      is_se       <- shared$lib_type_se
      mode        <- if (is_se) "SE" else "PE"
      threads     <- as.character(shared$salmon_threads)
      output_dir  <- shared$output_dir
      trim_dir    <- file.path(output_dir, "trimmed")
      fastqc_dir  <- file.path(output_dir, "fastqc_pre")
      quant_dir   <- file.path(output_dir, "salmon_quant")
      multiqc_dir <- file.path(output_dir, "multiqc")
      index_dir   <- if (shared$build_new_index) file.path("/data/references", "salmon_index") else shared$salmon_index_dir

      # Initialize sample statuses
      sample_names <- samples$name
      ss <- setNames(rep("pending", length(sample_names)), sample_names)
      rv$sample_status <- as.list(ss)

      # Calculate total steps
      n_samples <- nrow(samples)
      total <- 1  # FastQC
      if (shared$trimming_enabled) total <- total + n_samples  # fastp per sample
      if (shared$build_new_index)  total <- total + 1          # Salmon index
      total <- total + n_samples  # Salmon quant per sample
      total <- total + 1          # tximport
      total <- total + 1          # MultiQC
      rv$total_steps <- total

      add_log("=== Pipeline SalmonFlow iniciado ===", "info")
      add_log(paste("Muestras:", n_samples, "| Modo:", mode), "info")

      # ── STEP 1: FastQC ────────────────────────────────
      add_log("── Paso 1: FastQC (pre-trimming) ──", "info")
      all_fastqs <- samples$r1
      if (!is_se) all_fastqs <- c(all_fastqs, samples$r2)
      all_fastqs <- all_fastqs[!is.na(all_fastqs)]

      run_fastqc(all_fastqs, fastqc_dir,
                 threads = shared$salmon_threads,
                 log_callback = add_log)
      rv$current_step <- rv$current_step + 1
      if (rv$cancel_flag) return()

      # ── STEP 2: Trimmomatic ───────────────────────────
      trimmed_samples <- samples  # will update r1/r2 if trimming
      if (shared$trimming_enabled) {
        add_log("── Paso 2: fastp ──", "info")

        for (i in seq_len(nrow(samples))) {
          if (rv$cancel_flag) return()
          sname <- samples$name[i]
          rv$sample_status[[sname]] <- "running"

          r2_val <- if (is_se) NULL else samples$r2[i]

          trim_result <- run_fastp(
            r1                = samples$r1[i],
            r2                = r2_val,
            out_dir           = trim_dir,
            sample_name       = sname,
            mode              = mode,
            adapter_fasta     = shared$adapter_fasta,
            cut_front_quality = shared$fastp_cut_front,
            cut_tail_quality  = shared$fastp_cut_tail,
            cut_right_quality = shared$fastp_cut_right,
            minlen            = shared$fastp_minlen,
            threads           = shared$salmon_threads,
            log_callback      = add_log
          )

          if (trim_result$exit_status == 0) {
            trimmed_samples$r1[i] <- trim_result$r1_trimmed
            if (!is_se) trimmed_samples$r2[i] <- trim_result$r2_trimmed
            rv$sample_status[[sname]] <- "done"
          } else {
            rv$sample_status[[sname]] <- "error"
          }

          rv$current_step <- rv$current_step + 1
        }

        # Reset statuses for Salmon quant
        for (s in sample_names) {
          if (rv$sample_status[[s]] == "done") rv$sample_status[[s]] <- "pending"
        }
      } else {
        add_log("── Paso 2: fastp (omitido) ──", "info")
      }

      # ── STEP 3: Salmon Index ──────────────────────────
      if (shared$build_new_index) {
        add_log("── Paso 3: Construcción del índice Salmon ──", "info")
        if (rv$cancel_flag) return()

        decoy_file <- if (shared$decoy_aware) shared$genome_fasta else NULL

        idx_result <- build_salmon_index(
          fasta      = as.character(shared$transcriptome_fasta),
          outdir     = index_dir,
          decoy      = decoy_file,
          kmer       = shared$kmer_size,
          threads    = shared$salmon_threads,
          log_callback = add_log
        )
        rv$current_step <- rv$current_step + 1

        if (idx_result$exit_status != 0) {
          add_log("Pipeline abortado: error al construir el índice", "error")
          rv$running <- FALSE
          shared$pipeline_running <- FALSE
          return()
        }
        index_dir <- idx_result$index_dir
      } else {
        add_log("── Paso 3: Índice Salmon (usando existente) ──", "info")
        index_dir <- as.character(shared$salmon_index_dir)
      }

      # ── STEP 4: Salmon Quant ──────────────────────────
      add_log("── Paso 4: Cuantificación Salmon ──", "info")
      salmon_metas <- list()

      for (i in seq_len(nrow(trimmed_samples))) {
        if (rv$cancel_flag) return()
        sname <- trimmed_samples$name[i]
        rv$sample_status[[sname]] <- "running"

        r2_val <- if (is_se) NULL else trimmed_samples$r2[i]

        quant_result <- run_salmon_quant(
          index_dir      = index_dir,
          r1             = trimmed_samples$r1[i],
          r2             = r2_val,
          outdir         = quant_dir,
          sample_name    = sname,
          lib_type       = shared$salmon_libtype,
          gc_bias        = shared$salmon_gcbias,
          seq_bias       = shared$salmon_seqbias,
          threads        = shared$salmon_threads,
          is_se          = is_se,
          validate       = shared$salmon_validate,
          bootstraps     = shared$salmon_bootstraps,
          min_score_frac = shared$salmon_min_score_frac,
          discard_orphans = shared$salmon_discard_orphans,
          log_callback   = add_log
        )

        if (quant_result$exit_status == 0) {
          rv$sample_status[[sname]] <- "done"
          salmon_metas[[sname]] <- quant_result$meta
        } else {
          rv$sample_status[[sname]] <- "error"
        }

        rv$current_step <- rv$current_step + 1
        add_log(paste("  Salmon quant:", i, "/", nrow(trimmed_samples)), "info")
      }

      shared$salmon_meta <- salmon_metas

      # ── STEP 5: tximport ──────────────────────────────
      add_log("── Paso 5: tximport ──", "info")
      if (rv$cancel_flag) return()

      tx2gene <- build_tx2gene(as.character(shared$gtf_path), log_callback = add_log)
      if (is.null(tx2gene)) {
        add_log("Pipeline abortado: error al construir tx2gene", "error")
        rv$running <- FALSE
        shared$pipeline_running <- FALSE
        return()
      }

      count_matrix <- run_tximport(
        quant_dir          = quant_dir,
        sample_names       = sample_names,
        tx2gene            = tx2gene,
        method             = shared$txi_method,
        ignore_tx_version  = shared$txi_ignore_version,
        output_dir         = output_dir,
        log_callback       = add_log
      )

      rv$current_step <- rv$current_step + 1
      shared$count_matrix <- count_matrix

      # ── STEP 6: MultiQC ──────────────────────────────
      add_log("── Paso 6: MultiQC ──", "info")
      if (rv$cancel_flag) return()

      mqc <- run_multiqc(output_dir, multiqc_dir, log_callback = add_log)
      rv$current_step <- rv$current_step + 1

      if (file.exists(mqc$report_path)) {
        shared$multiqc_report <- mqc$report_path
      }

      # ── Done ──────────────────────────────────────────
      add_log("=== Pipeline completado exitosamente ===", "success")
      rv$running <- FALSE
      shared$pipeline_running <- FALSE
    })
  })
}
