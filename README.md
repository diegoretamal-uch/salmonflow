# SalmonFlow

A fully local, Dockerized **R Shiny** application for bulk RNA-seq analysis.

**Pipeline:** FastQC → fastp → Salmon → tximport → MultiQC

---

## Prerequisites

- [Docker Desktop](https://docs.docker.com/get-docker/) installed and running
  - Windows: requires WSL2 backend (Docker Desktop installs this automatically)
  - Linux/Mac: Docker Engine is sufficient

---

## Quick Start

### Linux / macOS

```bash
# 1. Clone the repo
git clone https://github.com/diegoretamal-uch/salmonflow
cd salmonflow

# 2. Build the image (~10 min, one-time)
docker build -t salmonflow .

# 3. Run
./run.sh /path/to/fastqs /path/to/references /path/to/output

# 4. Open browser → http://localhost:3838
```

### Windows (PowerShell)

```powershell
# 1. Clone the repo
git clone https://github.com/diegoretamal-uch/salmonflow
cd salmonflow

# 2. Build the image (~10 min, one-time)
docker build -t salmonflow .

# 3. Run
.\run.ps1 C:\path\to\fastqs C:\path\to\references C:\path\to\output

# 4. Open browser → http://localhost:3838
```

> **Note for Windows users:** if PowerShell blocks the script, run once:
> `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`

### Live-reload (development, Linux/macOS)

```bash
docker run --rm -p 3838:3838 \
  -v "$(pwd)/app:/srv/shiny-server/salmonflow" \
  -v "$(pwd)/data/input:/data/input" \
  -v "$(pwd)/data/references:/data/references" \
  -v "$(pwd)/data/output:/data/output" \
  -v "$(pwd)/data/tmp:/data/tmp" \
  salmonflow \
  R -e "options(shiny.autoreload=TRUE); shiny::runApp('/srv/shiny-server/salmonflow', host='0.0.0.0', port=3838)"
```

---

## Volume Mounts

| Container path      | Purpose                                        |
|---------------------|------------------------------------------------|
| `/data/input`       | FASTQ files (PE or SE, `.fastq.gz`)            |
| `/data/references`  | Transcriptome FASTA, GTF, Salmon index         |
| `/data/output`      | Results: quant, counts, FastQC, MultiQC        |
| `/data/tmp`         | Intermediate / temporary files                 |

---

## Software Versions

| Tool        | Version         | Role                                      |
|-------------|-----------------|-------------------------------------------|
| Salmon      | 1.10.0          | Quasi-mapping quantification              |
| FastQC      | 0.12.1          | Pre-trimming QC                           |
| fastp       | apt (≥ 0.23)    | Adapter trimming and quality filtering    |
| MultiQC     | 1.34            | Aggregated QC report                      |
| R           | 4.4.1           | Shiny runtime                             |
| Bioconductor| 3.19            | Bioinformatics package ecosystem          |
| tximport    | Bioc 3.19       | Salmon → count matrix                     |
| txdbmaker   | Bioc 3.19       | GTF → tx2gene table                       |

---

## Salmon Parameters

### Exposed in UI (Parametros tab)

| Parameter              | Default | Description                                                   |
|------------------------|---------|---------------------------------------------------------------|
| Library type (`-l`)    | `A`     | Auto-detect strand orientation. Set manually if needed (e.g., `ISR`, `ISF`). |
| GC bias (`--gcBias`)   | ON      | Corrects for GC content bias in fragment sampling.            |
| Seq bias (`--seqBias`) | ON      | Corrects for sequence-specific bias at read starts.           |
| Threads (`-p`)         | 4       | Parallelism. Recommended: leave 20-30% of cores for the OS.  |

### Advanced (collapsible, Parametros tab)

| Parameter                        | Default | Description                                                                              |
|----------------------------------|---------|------------------------------------------------------------------------------------------|
| Validate mappings (`--validateMappings`) | ON | Re-scores and filters mappings for accuracy. Recommended for most workflows.        |
| Bootstraps (`--numBootstraps`)   | 0       | Enables bootstrap sampling for quantification uncertainty. Set ≥ 100 for sleuth/DTU.   |
| Min score fraction (`--minScoreFraction`) | 0.65 | Fraction of the optimal alignment score a mapping must achieve to be retained. Lower = more permissive. |
| Discard orphans (`--discardOrphansQuasi`) | OFF | Discards reads whose mate did not map. Stricter paired-end mode.                  |

### Hardcoded (not exposed)

| Parameter                | Value  | Reason                                                        |
|--------------------------|--------|---------------------------------------------------------------|
| `--validateMappings`     | ON     | Best practice; exposed in Advanced to allow toggling          |
| `--writeUnmappedNames`   | OFF    | Output overhead not needed for standard quantification        |
| `--numGibbsSamples`      | 0      | Alternative to bootstraps; not exposed to avoid confusion     |
| `-k` (k-mer size)        | 31     | Standard for reads ≥ 75 bp; configurable in Referencias tab   |

---

## Pipeline Steps

1. **FastQC** — per-file quality report (pre-trimming)
2. **fastp** — adapter detection (auto for PE), quality trimming
3. **Salmon index** — build or reuse an existing index
4. **Salmon quant** — quasi-mapping quantification per sample
5. **tximport** — merge per-sample quant.sf into a gene-level count matrix
6. **MultiQC** — aggregate FastQC + fastp + Salmon reports into one HTML

---

## Tabs

1. **Muestras** — Load and validate FASTQ files, auto-detect PE pairs
2. **Referencias** — Select transcriptome FASTA, GTF, adapters, Salmon index
3. **Parametros** — Configure fastp, Salmon (standard + advanced), and tximport
4. **Ejecutar** — Run pipeline with live logs and per-sample progress
5. **Resultados** — Count matrix, Salmon QC, PCA, heatmap, MultiQC report
