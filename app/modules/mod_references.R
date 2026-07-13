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
          title = tagList(icon("info-circle"), "Quick Guide: Reference Files"),
          status = "info", solidHeader = FALSE, width = 12,
          collapsible = TRUE, collapsed = TRUE,
          tags$div(
            style = "line-height: 1.6;",
            tags$p("This tab configures the references required for alignment-free RNA-Seq quantification:"),
            tags$ul(
              tags$li(tags$strong("Transcriptome FASTA:"), " Required. FASTA containing transcript sequences (e.g. from GENCODE or Ensembl) to quantify."),
              tags$li(tags$strong("GTF Annotation:"), " Required. Maps transcript IDs to gene names, allowing tximport to aggregate transcript counts to the gene level."),
              tags$li(tags$strong("Adapter FASTA:"), " Optional. Overrides fastp's automatic adapter detection."),
              tags$li(tags$strong("Index Mode:"), " Use an existing pre-built Salmon index to save time, or build a new one from the transcriptome FASTA."),
              tags$li(tags$strong("Decoy-Aware Indexing:"), " Recommended. Includes the full genome as a 'decoy' to prevent genomic contamination from falsely mapping to transcripts. Requires a genome FASTA and >= 16 GB RAM."),
              tags$li(tags$strong("Sparse Index:"), " Reduces the index RAM footprint by 30-50% during build/run, but slightly slows down quantification.")
            )
          )
        )
      )
    ),

    fluidRow(
      column(12,
        box(
          title = tagList("Reference Files", html_tooltip("Select your transcriptome FASTA, GTF annotation, and optional adapter sequences.")),
          status = "primary", solidHeader = FALSE, width = 12,

          fluidRow(
            column(6,
              tags$label("Transcriptome FASTA"),
              shinyFilesButton(ns("fasta_file"), "Select",
                               title = "Select Transcriptome FASTA file",
                               multiple = FALSE),
              textOutput(ns("fasta_status"))
            ),
            column(6,
              tags$label("GTF Annotation"),
              shinyFilesButton(ns("gtf_file"), "Select",
                               title = "Select GTF annotation file",
                               multiple = FALSE),
              textOutput(ns("gtf_status"))
            )
          ),

          hr(),

          fluidRow(
            column(6,
              tags$label("Adapter FASTA (fastp, optional)"),
              shinyFilesButton(ns("adapter_file"), "Select",
                               title = "Select adapter FASTA",
                               multiple = FALSE),
              textOutput(ns("adapter_status")),
              helpText("Optional — fastp detects adapters automatically in PE data.")
            ),
            column(6,
              selectInput(ns("organism"), "Organism (informational)",
                          choices = c("Human" = "human",
                                      "Mouse" = "mouse",
                                      "Other" = "other"),
                          selected = "other")
            )
          )
        )
      )
    ),

    fluidRow(
      column(12,
        box(
          title = tagList("Salmon Index", html_tooltip("Configure whether to use an existing Salmon index or build a new one. Decoy-aware indexing improves quantification accuracy but requires more RAM.")),
          status = "primary", solidHeader = FALSE, width = 12,

          radioButtons(ns("index_mode"), "Index Mode",
                       choices = c("Use existing index"  = "existing",
                                   "Build new index" = "build"),
                       selected = "existing", inline = TRUE),

          conditionalPanel(
            condition = paste0("input['", ns("index_mode"), "'] == 'existing'"),
            shinyDirButton(ns("index_dir"), "Select index directory",
                           title = "Existing Salmon index directory"),
            textOutput(ns("index_dir_status"))
          ),

          hr(),

          fluidRow(
            column(4,
              checkboxInput(ns("decoy_aware"), "Decoy-aware indexing", value = FALSE),
              helpText("Higher accuracy — requires genome FASTA (.fa/.fa.gz). The pipeline generates gentrome and decoys.txt automatically.",
                       tags$br(),
                       tags$strong("System Requirements: "),
                       "including the full genome in the index significantly increases RAM",
                       " (and disk) usage during construction. For the human genome, ~16 GB of RAM or more is recommended.",
                       " If the system is low on memory, combine this option with a sparse index.")
            ),
            column(4,
              conditionalPanel(
                condition = paste0("input['", ns("decoy_aware"), "']"),
                tags$label("Genome FASTA"),
                tags$br(),
                shinyFilesButton(ns("genome_file"), "Select Genome FASTA",
                                 title = "Select genome FASTA for decoys",
                                 multiple = FALSE),
                textOutput(ns("genome_status"))
              )
            ),
            column(4,
              selectInput(ns("kmer_size"), "k-mer size",
                          choices = c(21, 23, 25, 27, 29, 31),
                          selected = 31),
              checkboxInput(ns("sparse_index"), "Sparse index (--sparse)", value = FALSE),
              helpText("Reduces RAM usage during index construction (~30-50%). Recommended if the system has low memory.",
                       tags$br(),
                       tags$strong("Note: "),
                       "uses a sparser hash table, making subsequent quantification",
                       " slightly slower.")
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
      if (is.null(p) || length(p) == 0) "Not selected"
      else paste("Selected:", basename(p))
    })

    output$gtf_status <- renderText({
      p <- parse_file(input$gtf_file)
      if (is.null(p) || length(p) == 0) "Not selected"
      else paste("Selected:", basename(p))
    })

    output$adapter_status <- renderText({
      p <- parse_file(input$adapter_file)
      if (is.null(p) || length(p) == 0) "Using default adapters"
      else paste("Selected:", basename(p))
    })

    output$genome_status <- renderText({
      p <- parse_file(input$genome_file)
      if (is.null(p) || length(p) == 0) "Not selected"
      else paste("Selected:", basename(p))
    })

    output$index_dir_status <- renderText({
      p <- parse_dir(input$index_dir)
      if (is.null(p) || length(p) == 0) "Not selected"
      else {
        # Check if it looks like a valid Salmon index
        has_info <- file.exists(file.path(p, "info.json"))
        if (has_info) paste("Valid index:", basename(p))
        else paste("Directory selected (info.json not found):", basename(p))
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
      shared$sparse_index        <- input$sparse_index

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
