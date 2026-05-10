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
            column(3,
              actionButton(ns("run_btn"), "Iniciar análisis",
                           class = "btn-primary", icon = icon("play"),
                           style = "width:100%; font-size:16px; padding:14px;")
            ),
            column(3,
              actionButton(ns("resume_btn"), "Reanudar análisis",
                           class = "btn-warning", icon = icon("rotate-right"),
                           style = "width:100%; font-size:16px; padding:14px;")
            ),
            column(3,
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
          # Static container — content is pushed via shinyjs::html() to avoid
          # Shiny re-rendering the element (which resets scrollTop every time).
          tags$div(
            id    = ns("log_content"),
            class = "log-panel",
            tags$span(class = "log-info", "Esperando inicio del análisis...")
          )
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

    # ── Reactive state ─────────────────────────────────────
    rv <- reactiveValues(
      log_pos       = 0L,
      running       = FALSE,
      proc          = NULL,
      log_file      = NULL,
      state_file    = NULL,
      current_step  = 0,
      total_steps   = 6,
      sample_status = list()
    )

    # Append a line to the static log div and scroll to bottom.
    # shinyjs::html() with add=TRUE appends without replacing the element,
    # so the browser never resets scrollTop.
    add_log <- function(msg, type = "info") {
      line <- timestamp_log(msg, type)
      shinyjs::html("log_content", html = paste0(line, "<br/>"), add = TRUE)
      shinyjs::runjs(sprintf(
        "(function(){ var el = document.getElementById('%s'); if(el) el.scrollTop = el.scrollHeight; })()",
        ns("log_content")
      ))
    }

    # ── Progress bar ───────────────────────────────────────
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
      if (length(statuses) == 0) return(tags$p(style = "color:#5a6373;", "—"))

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

    # ── Polling observer ───────────────────────────────────
    observe({
      if (!rv$running) return()
      invalidateLater(500, session)

      # Read new log lines from file and push directly into the DOM
      if (!is.null(rv$log_file) && file.exists(rv$log_file)) {
        all_lines <- readLines(rv$log_file, warn = FALSE)
        n_new <- length(all_lines) - rv$log_pos
        if (n_new > 0) {
          new_lines <- all_lines[(rv$log_pos + 1L):length(all_lines)]
          formatted <- vapply(new_lines, function(l) {
            sep <- regexpr("|", l, fixed = TRUE)
            if (sep > 0) {
              timestamp_log(substr(l, sep + 1L, nchar(l)),
                            substr(l, 1L,       sep - 1L))
            } else {
              timestamp_log(l, "info")
            }
          }, character(1L))
          rv$log_pos <- length(all_lines)

          new_html <- paste(formatted, collapse = "<br/>")
          shinyjs::html("log_content", html = paste0(new_html, "<br/>"), add = TRUE)
          shinyjs::runjs(sprintf(
            "(function(){ var el = document.getElementById('%s'); if(el) el.scrollTop = el.scrollHeight; })()",
            ns("log_content")
          ))
        }
      }

      # Read pipeline state
      if (!is.null(rv$state_file) && file.exists(rv$state_file)) {
        state <- tryCatch(
          jsonlite::fromJSON(rv$state_file),
          error = function(e) NULL
        )
        if (!is.null(state)) {
          rv$current_step <- state$step  %||% rv$current_step
          rv$total_steps  <- state$total %||% rv$total_steps
          ss <- state$sample_status
          if (!is.null(ss)) rv$sample_status <- as.list(ss)

          if (!isTRUE(state$running)) {
            rv$running              <- FALSE
            shared$pipeline_running <- FALSE

            # Load results written by the background script
            cm_path <- state$count_matrix_path %||% ""
            if (nchar(cm_path) > 0 && file.exists(cm_path)) {
              shared$count_matrix <- tryCatch(
                read.csv(cm_path, check.names = FALSE),
                error = function(e) NULL
              )
            }

            # Single unified MultiQC (no trimming)
            mq_path <- state$multiqc_report_path %||% ""
            if (nchar(mq_path) > 0 && file.exists(mq_path)) {
              shared$multiqc_report <- mq_path
            }

            # Dual MultiQC (trimming enabled)
            mq_pre_path <- state$multiqc_pre_report_path %||% ""
            if (nchar(mq_pre_path) > 0 && file.exists(mq_pre_path)) {
              shared$multiqc_pre_report <- mq_pre_path
            }
            mq_post_path <- state$multiqc_post_report_path %||% ""
            if (nchar(mq_post_path) > 0 && file.exists(mq_post_path)) {
              shared$multiqc_post_report <- mq_post_path
            }

            sm_path <- state$salmon_meta_path %||% ""
            if (nchar(sm_path) > 0 && file.exists(sm_path)) {
              shared$salmon_meta <- tryCatch(
                jsonlite::fromJSON(sm_path),
                error = function(e) NULL
              )
            }
          }
        }
      }

      # Catch unexpected process death
      if (!is.null(rv$proc) && !rv$proc$is_alive() && rv$running) {
        rv$running              <- FALSE
        shared$pipeline_running <- FALSE
        add_log("Pipeline terminado inesperadamente — revisa los logs", "warn")
      }
    })

    # ── Cancel button ──────────────────────────────────────
    observeEvent(input$cancel_btn, {
      if (rv$running) {
        if (!is.null(rv$proc) && rv$proc$is_alive()) rv$proc$kill()
        rv$running              <- FALSE
        shared$pipeline_running <- FALSE
        add_log("Pipeline cancelado por el usuario", "warn")
      }
    })

    # ── Shared launch logic ────────────────────────────────
    launch_pipeline <- function(resume = FALSE) {
      samples <- shared$samples
      if (is.null(samples) || nrow(samples) == 0) {
        add_log("Error: no hay muestras cargadas. Ve a la pestaña Muestras.", "error")
        return()
      }
      if (shared$build_new_index &&
          (is.null(shared$transcriptome_fasta) || length(shared$transcriptome_fasta) == 0)) {
        add_log("Error: no se ha seleccionado un transcriptoma FASTA. Ve a la pestaña Referencias.", "error")
        return()
      }
      if (is.null(shared$gtf_path) || length(shared$gtf_path) == 0) {
        add_log("Error: no se ha seleccionado un archivo GTF. Ve a la pestaña Referencias.", "error")
        return()
      }

      # Clear log panel, then show resume banner if applicable
      shinyjs::html("log_content", html = "", add = FALSE)
      if (resume) {
        add_log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", "warn")
        add_log("⟳ REANUDANDO — pasos con resultados previos seran omitidos automaticamente", "warn")
        add_log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", "warn")
      }

      # IPC file paths
      run_id     <- format(Sys.time(), "%Y%m%d_%H%M%S")
      tmp_dir    <- "/data/tmp"
      dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
      log_file   <- file.path(tmp_dir, paste0("sf_log_",    run_id, ".txt"))
      state_file <- file.path(tmp_dir, paste0("sf_state_",  run_id, ".json"))
      params_file <- file.path(tmp_dir, paste0("sf_params_", run_id, ".json"))

      # Serialize all params for the background script
      params <- list(
        samples               = samples,
        lib_type_se           = isTRUE(shared$lib_type_se),
        trimming_enabled      = isTRUE(shared$trimming_enabled),
        fastp_cut_front       = shared$fastp_cut_front,
        fastp_cut_tail        = shared$fastp_cut_tail,
        fastp_cut_right       = shared$fastp_cut_right,
        fastp_minlen          = shared$fastp_minlen,
        adapter_fasta         = shared$adapter_fasta %||% "",
        transcriptome_fasta   = shared$transcriptome_fasta %||% "",
        gtf_path              = shared$gtf_path %||% "",
        salmon_index_dir      = shared$salmon_index_dir %||% "",
        build_new_index       = isTRUE(shared$build_new_index),
        decoy_aware           = isTRUE(shared$decoy_aware),
        sparse_index          = isTRUE(shared$sparse_index),
        genome_fasta          = shared$genome_fasta %||% "",
        kmer_size             = shared$kmer_size,
        salmon_libtype        = shared$salmon_libtype,
        salmon_gcbias         = isTRUE(shared$salmon_gcbias),
        salmon_seqbias        = isTRUE(shared$salmon_seqbias),
        salmon_threads        = shared$salmon_threads,
        salmon_validate       = isTRUE(shared$salmon_validate),
        salmon_bootstraps     = shared$salmon_bootstraps,
        salmon_min_score_frac = shared$salmon_min_score_frac,
        salmon_discard_orphans = isTRUE(shared$salmon_discard_orphans),
        txi_method            = shared$txi_method,
        txi_ignore_version    = isTRUE(shared$txi_ignore_version),
        output_dir            = shared$output_dir,
        resume                = resume
      )
      jsonlite::write_json(params, params_file, auto_unbox = TRUE)

      # Reset state
      rv$log_pos      <- 0L
      rv$current_step <- 0
      rv$sample_status <- setNames(
        as.list(rep("pending", nrow(samples))), samples$name
      )
      rv$log_file   <- log_file
      rv$state_file <- state_file
      rv$running    <- TRUE
      shared$pipeline_running <- TRUE

      # Launch background Rscript
      script_path <- normalizePath("R/run_pipeline.R", mustWork = FALSE)
      stderr_file <- file.path(tmp_dir, paste0("sf_stderr_", run_id, ".txt"))

      rv$proc <- processx::process$new(
        "Rscript",
        args   = c(script_path, params_file, log_file, state_file),
        stdout = NULL,
        stderr = stderr_file
      )

      add_log(paste("Pipeline iniciado (PID:", rv$proc$get_pid(), ")"), "info")
    }

    # ── Run button ─────────────────────────────────────────
    observeEvent(input$run_btn, {
      launch_pipeline(resume = FALSE)
    })

    # ── Resume button ──────────────────────────────────────
    observeEvent(input$resume_btn, {
      launch_pipeline(resume = TRUE)
    })
  })
}
