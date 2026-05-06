# ══════════════════════════════════════════════════════════════
# SalmonFlow — mod_params.R
# Tab 3: Pipeline parameter configuration
# ══════════════════════════════════════════════════════════════

mod_params_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      # ── fastp Parameters ───────────────────────────────────
      column(6,
        box(
          title = "fastp",
          status = "primary", solidHeader = FALSE, width = 12,

          checkboxInput(ns("trim_enabled"), "Activar trimming", value = TRUE),

          conditionalPanel(
            condition = paste0("input['", ns("trim_enabled"), "']"),

            sliderInput(ns("fastp_cut_front"), "Cut front quality (5' end)",
                        min = 0, max = 30, value = 3, step = 1),

            sliderInput(ns("fastp_cut_tail"), "Cut tail quality (3' end)",
                        min = 0, max = 30, value = 3, step = 1),

            sliderInput(ns("fastp_cut_right"), "Sliding window quality (ventana = 4)",
                        min = 5, max = 30, value = 15, step = 1),

            sliderInput(ns("fastp_minlen"), "Largo mínimo (bp)",
                        min = 15, max = 100, value = 36, step = 1),

            helpText("Cut front/tail: calidad mínima en los extremos del read.",
                     "Sliding window: calidad promedio en ventana deslizante de 4 bases.",
                     "Los adaptadores se detectan automáticamente en datos PE.")
          )
        )
      ),

      # ── Salmon Parameters ──────────────────────────────────
      column(6,
        box(
          title = "Salmon",
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
                      min = 1, max = 16, value = 4, step = 1)
        ),

        box(
          title = "Salmon — Advanced",
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

          checkboxInput(ns("salmon_discard_orphans"),
                        "Discard orphan reads (--discardOrphansQuasi)",
                        value = FALSE),

          helpText(
            "Bootstraps: set >= 100 for uncertainty quantification (sleuth, DTU).",
            "Min score fraction: lower = more permissive mapping.",
            "Discard orphans: stricter PE mode, discards reads whose mate did not map."
          )
        ),

        box(
          title = "tximport",
          status = "primary", solidHeader = FALSE, width = 12,

          selectInput(ns("txi_method"), "Método de normalización",
                      choices = c("lengthScaledTPM (recomendado)" = "lengthScaledTPM",
                                  "scaledTPM"                      = "scaledTPM",
                                  "Raw counts (no)"                = "no",
                                  "TPM (dTU only)"                 = "dtuScaledTPM"),
                      selected = "lengthScaledTPM"),

          checkboxInput(ns("txi_ignore_version"),
                        "Ignorar versión de transcript IDs (Ensembl)",
                        value = TRUE),

          helpText("Si tus IDs tienen sufijo de versión (ej. ENST00000.5),",
                   "activa esta opción para compatibilidad con el GTF.")
        )
      )
    )
  )
}

mod_params_server <- function(id, shared) {
  moduleServer(id, function(input, output, session) {

    # Sync all parameters to shared reactive values
    observe({
      shared$trimming_enabled <- input$trim_enabled
      shared$fastp_cut_front  <- input$fastp_cut_front
      shared$fastp_cut_tail   <- input$fastp_cut_tail
      shared$fastp_cut_right  <- input$fastp_cut_right
      shared$fastp_minlen     <- input$fastp_minlen

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
