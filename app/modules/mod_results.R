# ══════════════════════════════════════════════════════════════
# SalmonFlow — mod_results.R
# Tab 5: Results — count matrix, QC, PCA, heatmap, MultiQC
# ══════════════════════════════════════════════════════════════

mod_results_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        box(
          title = "Matriz de Conteos (lengthScaledTPM)",
          status = "primary", solidHeader = FALSE, width = 12,
          collapsible = TRUE,

          DTOutput(ns("count_table")),

          br(),

          fluidRow(
            column(3,
              downloadButton(ns("download_csv"), "Descargar CSV",
                             class = "btn-success")
            ),
            column(3,
              downloadButton(ns("download_tsv"), "Descargar TSV",
                             class = "btn-success")
            )
          )
        )
      )
    ),

    fluidRow(
      column(6,
        box(
          title = "Salmon QC — Mapping Rates",
          status = "primary", solidHeader = FALSE, width = 12,
          DTOutput(ns("salmon_qc_table"))
        )
      ),
      column(6,
        box(
          title = "MultiQC Report",
          status = "primary", solidHeader = FALSE, width = 12,
          uiOutput(ns("multiqc_link")),
          helpText("El reporte MultiQC se abrirá en una nueva pestaña del navegador.")
        )
      )
    ),

    fluidRow(
      column(6,
        box(
          title = "PCA — Top Variable Genes",
          status = "primary", solidHeader = FALSE, width = 12,
          plotlyOutput(ns("pca_plot"), height = "450px")
        )
      ),
      column(6,
        box(
          title = "Heatmap — Top 50 Genes Más Variables",
          status = "primary", solidHeader = FALSE, width = 12,
          plotOutput(ns("heatmap_plot"), height = "500px")
        )
      )
    )
  )
}

mod_results_server <- function(id, shared) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Count matrix table ─────────────────────────────────
    output$count_table <- renderDT({
      cm <- shared$count_matrix
      if (is.null(cm)) return(NULL)

      datatable(cm,
        filter   = "top",
        rownames = FALSE,
        options  = list(
          pageLength = 25,
          scrollX    = TRUE,
          dom        = "lfrtip",
          order      = list()
        )
      ) %>%
        formatRound(columns = setdiff(names(cm), "gene_id"), digits = 2)
    })

    # ── Downloads ──────────────────────────────────────────
    output$download_csv <- downloadHandler(
      filename = function() {
        paste0("salmonflow_counts_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv")
      },
      content = function(file) {
        write.csv(shared$count_matrix, file, row.names = FALSE)
      }
    )

    output$download_tsv <- downloadHandler(
      filename = function() {
        paste0("salmonflow_counts_", format(Sys.time(), "%Y%m%d_%H%M"), ".tsv")
      },
      content = function(file) {
        write.table(shared$count_matrix, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )

    # ── Salmon QC summary table ────────────────────────────
    output$salmon_qc_table <- renderDT({
      meta <- shared$salmon_meta
      if (is.null(meta) || length(meta) == 0) return(NULL)

      qc_df <- data.frame(
        Sample          = names(meta),
        Reads_Processed = sapply(meta, function(m) m$num_processed %||% NA),
        Reads_Mapped    = sapply(meta, function(m) m$num_mapped %||% NA),
        Mapping_Rate    = sapply(meta, function(m) {
          if (!is.na(m$percent_mapped)) paste0(m$percent_mapped, "%") else "N/A"
        }),
        stringsAsFactors = FALSE
      )

      datatable(qc_df, rownames = FALSE,
        options = list(dom = "t", pageLength = 50, scrollX = TRUE)) %>%
        formatStyle("Mapping_Rate",
          color = styleInterval(c(50), c("#ef5350", "#66bb6a")),
          fontWeight = "bold")
    })

    # ── MultiQC link ───────────────────────────────────────
    output$multiqc_link <- renderUI({
      rpt <- shared$multiqc_report
      if (is.null(rpt) || !file.exists(rpt)) {
        tags$p(style = "color:#5a6373;", "No hay reporte MultiQC disponible aún.")
      } else {
        # Serve via Shiny addResourcePath
        addResourcePath("multiqc", dirname(rpt))
        tags$a(href = paste0("multiqc/", basename(rpt)),
               target = "_blank",
               class = "btn btn-primary",
               icon("external-link-alt"),
               " Abrir MultiQC Report")
      }
    })

    # ── PCA Plot ───────────────────────────────────────────
    output$pca_plot <- renderPlotly({
      cm <- shared$count_matrix
      samples <- shared$samples
      if (is.null(cm) || is.null(samples)) return(NULL)

      # Prepare numeric matrix (exclude gene_id)
      mat <- cm[, setdiff(names(cm), "gene_id"), drop = FALSE]
      mat <- as.matrix(mat)
      mat[is.na(mat)] <- 0

      # Log2 transform
      mat_log <- log2(mat + 1)

      # Keep top 500 most variable genes
      gene_vars <- apply(mat_log, 1, var)
      top_idx   <- head(order(gene_vars, decreasing = TRUE), 500)
      mat_top   <- mat_log[top_idx, , drop = FALSE]

      # PCA
      pca <- prcomp(t(mat_top), center = TRUE, scale. = TRUE)
      pca_df <- data.frame(
        PC1    = pca$x[, 1],
        PC2    = pca$x[, 2],
        Sample = colnames(mat_top),
        stringsAsFactors = FALSE
      )

      # Add group info
      if (!is.null(samples$group) && any(samples$group != "")) {
        grp_map <- setNames(samples$group, samples$name)
        pca_df$Group <- grp_map[pca_df$Sample]
      } else {
        pca_df$Group <- "all"
      }

      # Variance explained
      var_pct <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)

      p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, text = Sample)) +
        geom_point(size = 4, alpha = 0.85) +
        labs(
          x = paste0("PC1 (", var_pct[1], "%)"),
          y = paste0("PC2 (", var_pct[2], "%)"),
          title = "PCA — Top 500 Variable Genes"
        ) +
        theme_minimal(base_size = 13) +
        theme(
          plot.background  = element_rect(fill = "#ffffff", color = NA),
          panel.background = element_rect(fill = "#ffffff", color = NA),
          panel.grid       = element_line(color = "#e0e4ec"),
          text             = element_text(color = "#1f2430"),
          axis.text        = element_text(color = "#5a6373"),
          legend.background = element_rect(fill = "#ffffff"),
          legend.text       = element_text(color = "#1f2430")
        ) +
        scale_color_brewer(palette = "Set2")

      ggplotly(p, tooltip = c("text", "Group")) %>%
        layout(
          paper_bgcolor = "#ffffff",
          plot_bgcolor  = "#ffffff",
          font = list(color = "#1f2430")
        )
    })

    # ── Heatmap ────────────────────────────────────────────
    output$heatmap_plot <- renderPlot({
      cm <- shared$count_matrix
      samples <- shared$samples
      if (is.null(cm)) return(NULL)

      # Prepare numeric matrix
      mat <- cm[, setdiff(names(cm), "gene_id"), drop = FALSE]
      rownames(mat) <- cm$gene_id
      mat <- as.matrix(mat)
      mat[is.na(mat)] <- 0

      # Log2 transform
      mat_log <- log2(mat + 1)

      # Top 50 most variable genes
      gene_vars <- apply(mat_log, 1, var)
      top_idx   <- head(order(gene_vars, decreasing = TRUE), 50)
      mat_top   <- mat_log[top_idx, , drop = FALSE]

      # Annotation
      ann_col <- NULL
      if (!is.null(samples) && !is.null(samples$group) && any(samples$group != "")) {
        ann_col <- data.frame(
          Group = samples$group,
          row.names = samples$name
        )
      }

      pheatmap(mat_top,
        color          = colorRampPalette(c("#0277bd", "#ffffff", "#ef6c00"))(100),
        border_color   = NA,
        fontsize_row   = 7,
        fontsize_col   = 10,
        annotation_col = ann_col,
        main           = "Top 50 — Genes Más Variables",
        clustering_method = "ward.D2"
      )
    }, bg = "#ffffff")
  })
}
