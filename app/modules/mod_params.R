# ══════════════════════════════════════════════════════════════
# SalmonFlow — mod_params.R
# Tab 3: Pipeline parameter configuration
# ══════════════════════════════════════════════════════════════

mod_params_ui <- function(id) {
  ns <- NS(id)

  # Recommend threads based on the host machine: use all detected cores
  # but leave 1 free for the OS / UI. Falls back to 4 if detection fails.
  n_cores      <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
  rec_threads  <- if (is.na(n_cores) || n_cores < 1) 4L else max(1L, min(32L, n_cores - 1L))

  tagList(
    fluidRow(
      column(12,
        box(
          title = tagList(icon("info-circle"), "Quick Guide: Parameter Tuning"),
          status = "info", solidHeader = FALSE, width = 12,
          collapsible = TRUE, collapsed = TRUE,
          tags$div(
            style = "line-height: 1.6;",
            tags$p("This tab allows you to configure processing settings for each step of the pipeline:"),
            tags$ul(
              tags$li(tags$strong("fastp (Trimming):"), " Quality-based filtering and adapter removal. Use the sliders to set minimum quality thresholds for read ends (5'/3'), the sliding window average quality, and the minimum read length (recommended: 30-50bp) to discard short, non-specific reads."),
              tags$li(tags$strong("Temporary Files:"), " Policy for deleting intermediate FASTQ files to manage disk space. You can choose to delete raw FASTQs after trimming (saves ~50% peak space), or trimmed FASTQs after Salmon quantification."),
              tags$li(tags$strong("Salmon (Quantification):"), " Set strandedness (use 'A' for auto-detection) and enable GC/Sequence bias corrections to account for PCR and sequencing-specific biases."),
              tags$li(tags$strong("Salmon (Advanced):"), " Enable 'Validate mappings' for selective alignment (highly recommended). Set bootstraps to >= 100 to evaluate technical/quantification uncertainty. Keep at 0 for standard gene-level differential expression."),
              tags$li(tags$strong("tximport (Normalization):"), " Methods to summarize transcript-level abundance estimates to gene-level count matrices. 'lengthScaledTPM' is recommended for downstream DESeq2/EdgeR analyses. Enable 'Ignore version' to strip Ensembl decimal suffixes (e.g. '.1').")
            )
          )
        )
      )
    ),

    fluidRow(
      # ── fastp Parameters ───────────────────────────────────
      column(6,
        box(
          title = tagList("fastp", html_tooltip("Configure quality trimming, adapter removal, and minimum read length filters for the raw sequencing reads.")),
          status = "primary", solidHeader = FALSE, width = 12,

           checkboxInput(ns("trim_enabled"), "Enable trimming", value = TRUE),

          conditionalPanel(
            condition = paste0("input['", ns("trim_enabled"), "']"),

            sliderInput(ns("fastp_cut_front"), "Cut front quality (5' end)",
                        min = 0, max = 30, value = 20, step = 1),

            sliderInput(ns("fastp_cut_tail"), "Cut tail quality (3' end)",
                        min = 0, max = 30, value = 20, step = 1),

            sliderInput(ns("fastp_cut_right"), "Sliding window quality",
                        min = 5, max = 30, value = 20, step = 1),

            sliderInput(ns("fastp_window_size"), "Sliding window size (bases)",
                        min = 1, max = 20, value = 4, step = 1),

            sliderInput(ns("fastp_minlen"), "Minimum length (bp)",
                        min = 15, max = 100, value = 36, step = 1),

            helpText("Cut front/tail: minimum quality at read ends.",
                     "Sliding window: average quality in a sliding window of N bases.",
                     "Adapters are automatically detected in PE data.")
          )
        ),

        conditionalPanel(
          condition = paste0("input['", ns("trim_enabled"), "']"),

           box(
            title = tagList("Temporary Files", html_tooltip("Choose when to delete intermediate FASTQ files. This keeps the disk footprint minimal during pipeline runs.")),
            status = "primary", solidHeader = FALSE, width = 12,

            radioButtons(ns("temp_file_policy"),
                         "Disk space management during trimming",
                         choices = c(
                           "Keep all (raw + trimmed)" = "keep_both",
                           "Delete each raw FASTQ once its trim succeeds" = "delete_raw",
                           "Delete each trimmed FASTQ once its Salmon quant.sf succeeds" = "delete_trimmed"
                         ),
                         selected = "keep_both"),

            conditionalPanel(
              condition = paste0("input['", ns("temp_file_policy"), "'] != 'keep_both'"),
              tags$div(class = "alert alert-danger", style = "margin-bottom:10px;",
                icon("triangle-exclamation"),
                tags$strong(" Warning: "),
                "The corresponding FASTQ files will be permanently deleted ",
                "once the next step completes successfully. ",
                "This action cannot be undone. Make sure you have a backup copy ",
                "of your data if you need them later."
              )
            ),

            helpText("With trimming enabled, trimmed files can double disk space usage.",
                     "Deleting each raw FASTQ after its trim keeps the peak storage at ~1×",
                     "instead of 2×. Deleting each trimmed FASTQ after its quantification frees that space",
                     "as soon as Salmon no longer needs it. Both deletions only occur if the next file",
                     "was successfully generated (exists and is non-empty).")
          )
        )
      ),

      # ── Salmon Parameters ──────────────────────────────────
      column(6,
        box(
          title = tagList("Salmon", html_tooltip("Set the quantification parameters, thread count, and bias correction options for Salmon.")),
          status = "primary", solidHeader = FALSE, width = 12,

          selectInput(ns("salmon_libtype"), "Library type",
                      choices = c("A (auto-detect)" = "A",
                                  "ISR (Illumina PE, fr-firststrand)" = "ISR",
                                  "ISF (Illumina PE, fr-secondstrand)" = "ISF",
                                  "IU (Illumina PE, unstranded)" = "IU",
                                  "SR (SE, reverse)" = "SR",
                                  "SF (SE, forward)" = "SF",
                                  "U (SE, unstranded)" = "U"),
                      selected = "A"),

          checkboxInput(ns("salmon_gcbias"), "GC bias correction", value = TRUE),

          checkboxInput(ns("salmon_seqbias"), "Sequence bias correction", value = TRUE),

          sliderInput(ns("salmon_threads"), "Threads",
                      min = 1, max = 32, value = rec_threads, step = 1),

          helpText(
            if (is.na(n_cores)) {
              "Could not detect the number of system cores; using 4 by default."
            } else {
              sprintf(paste("System with %d cores detected.",
                            "Recommended: %d (leaves 1 free for the system)."),
                      n_cores, rec_threads)
            }
          )
        ),

        box(
          title = tagList("Salmon — Advanced", html_tooltip("Configure advanced Salmon options such as mapping validation, bootstrap sampling, and minimum score fraction.")),
          status = "primary", solidHeader = FALSE, width = 12,
          collapsible = TRUE, collapsed = TRUE,

          checkboxInput(ns("salmon_validate"),
                        "Validate mappings (--validateMappings)",
                        value = TRUE),

          sliderInput(ns("salmon_bootstraps"), "Bootstrap samples (--numBootstraps)",
                      min = 0, max = 200, value = 0, step = 10),

          numericInput(ns("salmon_min_score_frac"),
                       "Min score fraction (--minScoreFraction)",
                       value = 0.65, min = 0.0, max = 1.0, step = 0.05),

          if (FALSE) {
            checkboxInput(ns("salmon_discard_orphans"),
                          "Discard orphan reads (--discardOrphansQuasi)",
                          value = FALSE)
          },

          helpText(
            "Bootstraps: set >= 100 for uncertainty quantification (sleuth, DTU).",
            "Min score fraction: lower = more permissive mapping."
          )
        ),

        box(
          title = tagList("tximport", html_tooltip("Configure how transcript-level counts are aggregated to gene-level counts using tximport.")),
          status = "primary", solidHeader = FALSE, width = 12,

          selectInput(ns("txi_method"), "Normalization method",
                      choices = c("lengthScaledTPM (recommended)" = "lengthScaledTPM",
                                  "scaledTPM"                      = "scaledTPM",
                                  "Raw counts (no)"                = "no",
                                  "TPM (dTU only)"                 = "dtuScaledTPM"),
                      selected = "lengthScaledTPM"),

          checkboxInput(ns("txi_ignore_version"),
                        "Ignore transcript ID version (Ensembl)",
                        value = TRUE),

          helpText("If your IDs have a version suffix (e.g. ENST00000.5),",
                   "enable this option for compatibility with the GTF.")
        )
      )
    )
  )
}

mod_params_server <- function(id, shared) {
  moduleServer(id, function(input, output, session) {

    # Sync all parameters to shared reactive values
    observe({
      shared$trimming_enabled   <- input$trim_enabled
      shared$fastp_cut_front    <- input$fastp_cut_front
      shared$fastp_cut_tail     <- input$fastp_cut_tail
      shared$fastp_cut_right    <- input$fastp_cut_right
      shared$fastp_window_size  <- input$fastp_window_size
      shared$fastp_minlen       <- input$fastp_minlen

      shared$delete_originals_after_trim <- (input$temp_file_policy == "delete_raw")
      shared$delete_trimmed_after_quant  <- (input$temp_file_policy == "delete_trimmed")

      shared$salmon_libtype         <- input$salmon_libtype
      shared$salmon_gcbias          <- input$salmon_gcbias
      shared$salmon_seqbias         <- input$salmon_seqbias
      shared$salmon_threads         <- input$salmon_threads
      shared$salmon_validate        <- input$salmon_validate
      shared$salmon_bootstraps      <- input$salmon_bootstraps
      shared$salmon_min_score_frac  <- input$salmon_min_score_frac
      shared$salmon_discard_orphans <- input$salmon_discard_orphans

      shared$txi_method         <- input$txi_method
      shared$txi_ignore_version <- input$txi_ignore_version
    })
  })
}
