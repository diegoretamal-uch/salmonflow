# ══════════════════════════════════════════════════════════════
# SalmonFlow — app.R (main entry point)
# ══════════════════════════════════════════════════════════════

library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyFiles)
library(DT)
library(ggplot2)
library(plotly)
library(dplyr)
library(readr)
library(tidyr)
library(processx)
library(future)
library(promises)
library(waiter)
library(jsonlite)
library(pheatmap)

# ── Set async plan ───────────────────────────────────────────
plan(multisession)

# ── Source modules & helpers ─────────────────────────────────
source("modules/mod_upload.R")
source("modules/mod_references.R")
source("modules/mod_params.R")
source("modules/mod_run.R")
source("modules/mod_results.R")
source("R/helpers.R")
source("R/pipeline_functions.R")
source("R/tximport_utils.R")

# ── Data directory roots (inside Docker) ─────────────────────
DATA_INPUT  <- "/data/input"
DATA_REF    <- "/data/references"
DATA_OUTPUT <- "/data/output"
DATA_TMP    <- "/data/tmp"
VOLUMES     <- c(home = DATA_INPUT, refs = DATA_REF, output = DATA_OUTPUT)

# ══════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════
ui <- dashboardPage(
  skin = "black",

  dashboardHeader(
    title = span(
      icon("fish"), "SalmonFlow",
      style = "font-weight:700; letter-spacing:1px;"
    ),
    titleWidth = 260
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      id = "main_tabs",
      menuItem("Muestras",    tabName = "tab_samples",    icon = icon("folder-open")),
      menuItem("Referencias", tabName = "tab_references", icon = icon("dna")),
      menuItem("Parámetros",  tabName = "tab_params",     icon = icon("sliders-h")),
      menuItem("Ejecutar",    tabName = "tab_run",        icon = icon("play")),
      menuItem("Resultados",  tabName = "tab_results",    icon = icon("chart-bar"))
    ),
    tags$div(
      style = "position:absolute; bottom:10px; width:100%; text-align:center;
               color:#888; font-size:11px;",
      "SalmonFlow v1.0"
    )
  ),

  dashboardBody(
    useShinyjs(),
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),

    tabItems(
      tabItem(tabName = "tab_samples",    mod_upload_ui("upload")),
      tabItem(tabName = "tab_references", mod_references_ui("refs")),
      tabItem(tabName = "tab_params",     mod_params_ui("params")),
      tabItem(tabName = "tab_run",        mod_run_ui("run")),
      tabItem(tabName = "tab_results",    mod_results_ui("results"))
    )
  )
)

# ══════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # ── Shared reactive values across modules ──────────────────
  shared <- reactiveValues(
    # Sample info
    samples      = NULL,   # data.frame: name, r1, r2, group
    lib_type_se  = FALSE,  # TRUE = single-end

    # Reference paths
    transcriptome_fasta = NULL,
    gtf_path            = NULL,
    adapter_fasta       = NULL,
    salmon_index_dir    = NULL,
    build_new_index     = TRUE,
    decoy_aware         = FALSE,
    genome_fasta        = NULL,
    kmer_size           = 31,
    organism            = "other",

    # Parameters
    trimming_enabled    = TRUE,
    fastp_cut_front     = 3,
    fastp_cut_tail      = 3,
    fastp_cut_right     = 15,
    fastp_minlen        = 36,
    salmon_libtype      = "A",
    salmon_gcbias       = TRUE,
    salmon_seqbias      = TRUE,
    salmon_threads      = 4,
    txi_method          = "lengthScaledTPM",
    txi_ignore_version  = TRUE,

    # Pipeline state
    pipeline_running = FALSE,
    pipeline_log     = character(0),
    pipeline_step    = 0,
    pipeline_total   = 6,
    pipeline_process = NULL,

    # Results
    count_matrix    = NULL,
    salmon_meta     = NULL,
    multiqc_report  = NULL,
    output_dir      = DATA_OUTPUT
  )

  # ── Call module servers ────────────────────────────────────
  mod_upload_server("upload",      shared, VOLUMES)
  mod_references_server("refs",    shared, VOLUMES)
  mod_params_server("params",      shared)
  mod_run_server("run",            shared)
  mod_results_server("results",    shared)
}

# ── Launch ───────────────────────────────────────────────────
shinyApp(ui, server)
