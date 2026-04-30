# 🐟 SalmonFlow

A fully local, Dockerized **R Shiny** application for bulk RNA-seq analysis.

Pipeline: **FastQC → fastp → Salmon → tximport → MultiQC**

## Prerequisites

- [Docker Desktop](https://docs.docker.com/desktop/) installed

## Quick Start

```bash
# 1. Clone the repo
git clone https://github.com/your-user/salmonflow
cd salmonflow

# 2. Build the image
docker build -t salmonflow .

# 3. Run with your data
./run.sh /path/to/fastqs /path/to/references /path/to/output

# 4. Open in browser
# http://localhost:3838
```

## Volume Mounts

| Container path      | Purpose                           |
|----------------------|-----------------------------------|
| `/data/input`        | Your FASTQ files                  |
| `/data/references`   | Transcriptome FASTA, GTF, indices |
| `/data/output`       | Pipeline results                  |
| `/data/tmp`          | Intermediate files                |

## Using Docker Compose

```bash
# Edit docker-compose.yml to set your local paths, then:
docker compose up --build
```

## Tabs

1. **📁 Muestras** — Load and validate FASTQ files
2. **🧬 Referencias** — Select transcriptome, GTF, adapters, Salmon index
3. **⚙️ Parámetros** — Configure Trimmomatic, Salmon, and tximport
4. **▶️ Ejecutar** — Run pipeline with live logs and progress
5. **📊 Resultados** — Interactive count matrix, PCA, heatmap, MultiQC report
