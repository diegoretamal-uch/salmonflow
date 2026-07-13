# ══════════════════════════════════════════════════════════════
# SalmonFlow — mod_upload.R
# Tab 1: Sample FASTQ loading and validation
# ══════════════════════════════════════════════════════════════

mod_upload_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        box(
          title = tagList(icon("info-circle"), "Quick Guide: Loading Samples"),
          status = "info", solidHeader = FALSE, width = 12,
          collapsible = TRUE, collapsed = TRUE,
          tags$div(
            style = "line-height: 1.6;",
            tags$p("This tab allows you to configure your input RNA-Seq datasets:"),
            tags$ul(
              tags$li(tags$strong("Library Type:"), " Choose whether your samples are Paired-end (two files per sample, e.g., _R1 and _R2) or Single-end (one file)."),
              tags$li(tags$strong("Select FASTQs Folder:"), " Click the button to select the directory containing your raw sequencing data."),
              tags$li(tags$strong("Auto-detect Pairs:"), " Click this to automatically identify and pair your samples from the selected folder."),
              tags$li(tags$strong("Validate Samples:"), " Runs a check to verify that all specified files exist and are readable before launching.")
            )
          )
        )
      )
    ),

    fluidRow(
      column(12,
        box(
          title = tagList("Load Samples", html_tooltip("Select your FASTQ folder, configure the library type, and auto-detect or manually define the sample metadata table.")),
          status = "primary", solidHeader = FALSE, width = 12,

          fluidRow(
            column(4,
              radioButtons(ns("lib_type"), "Library Type",
                           choices  = c("Paired-end" = "PE", "Single-end" = "SE"),
                           selected = "PE", inline = TRUE)
            ),
            column(4,
              shinyDirButton(ns("fastq_dir"), "Select FASTQs Folder",
                             title = "Select folder with FASTQ files",
                             icon = icon("folder-open"))
            ),
            column(4,
              tags$div(style = "padding-top:25px;",
                verbatimTextOutput(ns("selected_dir_text"), placeholder = TRUE)
              )
            )
          ),

          hr(),

          DTOutput(ns("sample_table")),

          br(),

          fluidRow(
            column(3,
              actionButton(ns("add_row"), "+ Add Row",
                           class = "btn-default", icon = icon("plus"))
            ),
            column(3,
              actionButton(ns("remove_row"), "- Remove Selected",
                           class = "btn-default", icon = icon("minus"))
            ),
            column(3,
              actionButton(ns("auto_detect"), "Auto-detect Pairs",
                           class = "btn-info", icon = icon("search"))
            ),
            column(3,
              actionButton(ns("validate_btn"), "Validate Samples",
                           class = "btn-success", icon = icon("check"))
            )
          ),

          br(),
          uiOutput(ns("validation_msg"))
        )
      )
    )
  )
}

mod_upload_server <- function(id, shared, volumes) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── shinyFiles directory chooser ───────────────────────
    shinyDirChoose(input, "fastq_dir", roots = volumes, session = session)

    selected_dir <- reactive({
      if (is.integer(input$fastq_dir)) return(NULL)
      parseDirPath(volumes, input$fastq_dir)
    })

    output$selected_dir_text <- renderText({
      d <- selected_dir()
      if (is.null(d) || length(d) == 0) "No folder selected"
      else as.character(d)
    })

    # ── Local editable table ───────────────────────────────
    rv <- reactiveValues(
      table_data = data.frame(
        name  = character(0),
        r1    = character(0),
        r2    = character(0),
        group = character(0),
        stringsAsFactors = FALSE
      )
    )

    # Auto-detect pairs when directory is selected
    observeEvent(input$auto_detect, {
      d <- selected_dir()
      if (is.null(d) || length(d) == 0) return()
      detected <- detect_fastq_pairs(as.character(d))
      if (nrow(detected) > 0) {
        rv$table_data <- detected
      }
    })

    # Also auto-detect when folder changes
    observeEvent(selected_dir(), {
      d <- selected_dir()
      if (is.null(d) || length(d) == 0) return()
      detected <- detect_fastq_pairs(as.character(d))
      if (nrow(detected) > 0) {
        rv$table_data <- detected
      }
    })

    # Add row
    observeEvent(input$add_row, {
      rv$table_data <- rbind(rv$table_data,
        data.frame(name = "", r1 = "", r2 = "", group = "", stringsAsFactors = FALSE))
    })

    # Remove selected row
    observeEvent(input$remove_row, {
      sel <- input$sample_table_rows_selected
      if (!is.null(sel) && length(sel) > 0) {
        rv$table_data <- rv$table_data[-sel, , drop = FALSE]
      }
    })

    # Render editable table
    output$sample_table <- renderDT({
      df <- rv$table_data
      if (input$lib_type == "SE") df$r2 <- NULL

      datatable(df,
        editable = TRUE,
        selection = "single",
        options = list(
          pageLength = 20,
          dom = "t",
          scrollX = TRUE,
          columnDefs = list(list(className = "dt-center", targets = "_all"))
        ),
        rownames = FALSE
      )
    })

    # Handle cell edits
    observeEvent(input$sample_table_cell_edit, {
      info <- input$sample_table_cell_edit
      row <- info$row
      col <- info$col + 1  # DT is 0-indexed for columns
      val <- info$value

      col_names <- if (input$lib_type == "SE") c("name", "r1", "group") else c("name", "r1", "r2", "group")
      if (col <= length(col_names)) {
        rv$table_data[row, col_names[col]] <- val
      }
    })

    # Validate
    observeEvent(input$validate_btn, {
      df <- rv$table_data
      if (nrow(df) == 0) {
        output$validation_msg <- renderUI(
          tags$div(class = "validation-warn", "No samples loaded")
        )
        return()
      }

      # Collect all FASTQ paths
      paths <- df$r1
      if (input$lib_type == "PE") paths <- c(paths, df$r2)
      paths <- paths[!is.na(paths) & paths != ""]

      result <- validate_fastq_files(paths)

      shared$samples     <- df
      shared$lib_type_se <- (input$lib_type == "SE")

      output$validation_msg <- renderUI({
        msgs <- lapply(result$messages, function(m) {
          cls <- if (grepl("^\\[OK\\]", m)) "validation-ok"
                 else "validation-warn"
          tags$div(class = cls, m)
        })
        if (result$valid) {
          tagList(
            tags$div(class = "validation-ok",
                     paste("All", nrow(df), "samples validated successfully")),
            msgs
          )
        } else {
          tagList(
            tags$div(class = "validation-warn", "Some files have issues:"),
            msgs
          )
        }
      })
    })

    # Sync lib_type
    observe({
      shared$lib_type_se <- (input$lib_type == "SE")
    })
  })
}
