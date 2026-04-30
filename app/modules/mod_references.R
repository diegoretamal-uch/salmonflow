# ══════════════════════════════════════════════════════════════
# SalmonFlow — mod_references.R
# Tab 2: Reference files — transcriptome, GTF, adapters, index
# ══════════════════════════════════════════════════════════════

mod_references_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        box(
          title = "🧬 Archivos de Referencia",
          status = "primary", solidHeader = FALSE, width = 12,

          fluidRow(
            column(6,
              tags$label("Transcriptoma FASTA"),
              shinyFilesButton(ns("fasta_file"), "📂 Seleccionar",
                               title = "Seleccionar archivo FASTA del transcriptoma",
                               multiple = FALSE),
              textOutput(ns("fasta_status"))
            ),
            column(6,
              tags$label("Anotación GTF"),
              shinyFilesButton(ns("gtf_file"), "📂 Seleccionar",
                               title = "Seleccionar archivo GTF de anotación",
                               multiple = FALSE),
              textOutput(ns("gtf_status"))
            )
          ),

          hr(),

          fluidRow(
            column(6,
              tags$label("Adaptadores FASTA (fastp, opcional)"),
              shinyFilesButton(ns("adapter_file"), "📂 Seleccionar",
                               title = "Seleccionar FASTA de adaptadores",
                               multiple = FALSE),
              textOutput(ns("adapter_status")),
              helpText("Opcional — fastp detecta adaptadores automáticamente en datos PE.")
            ),
            column(6,
              selectInput(ns("organism"), "Organismo (informativo)",
                          choices = c("Humano" = "human",
                                      "Ratón" = "mouse",
                                      "Otro" = "other"),
                          selected = "other")
            )
          )
        )
      )
    ),

    fluidRow(
      column(12,
        box(
          title = "🐟 Índice Salmon",
          status = "primary", solidHeader = FALSE, width = 12,

          radioButtons(ns("index_mode"), "Modo de índice",
                       choices = c("Construir nuevo índice" = "build",
                                   "Usar índice existente"  = "existing"),
                       selected = "build", inline = TRUE),

          conditionalPanel(
            condition = paste0("input['", ns("index_mode"), "'] == 'existing'"),
            shinyDirButton(ns("index_dir"), "📂 Seleccionar directorio del índice",
                           title = "Directorio del índice Salmon existente"),
            textOutput(ns("index_dir_status"))
          ),

          hr(),

          fluidRow(
            column(4,
              checkboxInput(ns("decoy_aware"), "Decoy-aware indexing", value = FALSE),
              helpText("Recomendado para mayor precisión — requiere genoma FASTA")
            ),
            column(4,
              conditionalPanel(
                condition = paste0("input['", ns("decoy_aware"), "']"),
                shinyFilesButton(ns("genome_file"), "📂 Genoma FASTA",
                                 title = "Seleccionar genoma FASTA para decoys",
                                 multiple = FALSE),
                textOutput(ns("genome_status"))
              )
            ),
            column(4,
              selectInput(ns("kmer_size"), "k-mer size",
                          choices = c(21, 23, 25, 27, 29, 31),
                          selected = 31)
            )
          )
        )
      )
    )
  )
}

mod_references_server <- function(id, shared, volumes) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── File choosers ──────────────────────────────────────
    shinyFileChoose(input, "fasta_file",   roots = volumes, session = session,
                    filetypes = c("fa", "fasta", "gz"))
    shinyFileChoose(input, "gtf_file",     roots = volumes, session = session,
                    filetypes = c("gtf", "gz"))
    shinyFileChoose(input, "adapter_file", roots = volumes, session = session,
                    filetypes = c("fa", "fasta"))
    shinyFileChoose(input, "genome_file",  roots = volumes, session = session,
                    filetypes = c("fa", "fasta", "gz"))
    shinyDirChoose(input, "index_dir",     roots = volumes, session = session)

    # ── Helper to parse file path ──────────────────────────
    parse_file <- function(input_val) {
      if (is.integer(input_val)) return(NULL)
      parseFilePaths(volumes, input_val)$datapath
    }

    parse_dir <- function(input_val) {
      if (is.integer(input_val)) return(NULL)
      parseDirPath(volumes, input_val)
    }

    # ── Reactives for status display ───────────────────────
    output$fasta_status <- renderText({
      p <- parse_file(input$fasta_file)
      if (is.null(p) || length(p) == 0) "No seleccionado"
      else paste("✓", basename(p))
    })

    output$gtf_status <- renderText({
      p <- parse_file(input$gtf_file)
      if (is.null(p) || length(p) == 0) "No seleccionado"
      else paste("✓", basename(p))
    })

    output$adapter_status <- renderText({
      p <- parse_file(input$adapter_file)
      if (is.null(p) || length(p) == 0) "Usando adaptadores por defecto"
      else paste("✓", basename(p))
    })

    output$genome_status <- renderText({
      p <- parse_file(input$genome_file)
      if (is.null(p) || length(p) == 0) "No seleccionado"
      else paste("✓", basename(p))
    })

    output$index_dir_status <- renderText({
      p <- parse_dir(input$index_dir)
      if (is.null(p) || length(p) == 0) "No seleccionado"
      else {
        # Check if it looks like a valid Salmon index
        has_info <- file.exists(file.path(p, "info.json"))
        if (has_info) paste("✓ Índice válido:", basename(p))
        else paste("⚠ Directorio seleccionado (no se encontró info.json):", basename(p))
      }
    })

    # ── Sync to shared reactive values ─────────────────────
    observe({
      shared$transcriptome_fasta <- parse_file(input$fasta_file)
      shared$gtf_path            <- parse_file(input$gtf_file)
      shared$organism            <- input$organism
      shared$kmer_size           <- as.integer(input$kmer_size)
      shared$build_new_index     <- (input$index_mode == "build")
      shared$decoy_aware         <- input$decoy_aware

      adapter <- parse_file(input$adapter_file)
      shared$adapter_fasta <- if (!is.null(adapter) && length(adapter) > 0) as.character(adapter) else NULL

      if (input$index_mode == "existing") {
        shared$salmon_index_dir <- parse_dir(input$index_dir)
      }

      if (input$decoy_aware) {
        shared$genome_fasta <- parse_file(input$genome_file)
      }
    })
  })
}
