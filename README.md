# SalmonFlow

A fully local, Dockerized **R Shiny** application for bulk RNA-seq analysis.

**Pipeline:** FastQC → fastp → Salmon → tximport → MultiQC

## Tutorial

https://github.com/user-attachments/assets/a7072ae9-3d94-43aa-9d24-74ed035acb77


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

# 3. Optional, one-time: pull iDEP for downstream DE/pathway analysis
docker pull gexijin/idep:latest

# 4. Run — starts iDEP too, if its image is present
./run.sh /path/to/fastqs /path/to/references /path/to/output

# 5. Open browser → http://localhost:3838   (iDEP: http://localhost:3839)
```

### Windows (PowerShell)

```powershell
# 1. Clone the repo
git clone https://github.com/diegoretamal-uch/salmonflow
cd salmonflow

# 2. Build the image (~10 min, one-time)
docker build -t salmonflow .

# 3. Optional, one-time: pull iDEP for downstream DE/pathway analysis
docker pull gexijin/idep:latest

# 4. Run — starts iDEP too, if its image is present
.\run.ps1 C:\path\to\fastqs C:\path\to\references C:\path\to\output

# 5. Open browser → http://localhost:3838   (iDEP: http://localhost:3839)
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

## Data Folders

No manual setup required. The `data/` folder comes pre-created when you clone the repo:

```
salmonflow/
├── data/
│   ├── input/        ← drop your FASTQ files here
│   ├── references/   ← transcriptome FASTA, GTF, Salmon index
│   ├── output/       ← pipeline results written here
│   └── tmp/          ← intermediate files (auto-cleaned)
```

By default the scripts use these folders. You can also pass **any absolute path** on your machine — the run scripts create the directories automatically if they don't exist:

```bash
# Linux/macOS — use external paths
./run.sh /mnt/data/fastqs /mnt/refs /mnt/results

# Windows — use external paths
.\run.ps1 D:\data\fastqs D:\refs D:\results
```

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
5. **Resultados** — Count matrix, Salmon QC, iDEP export, MultiQC report

---

## Downstream Analysis with iDEP

SalmonFlow stops at the count matrix. For differential expression (DESeq2,
limma, edgeR), clustering and pathway enrichment, the Results tab hands off to
[**iDEP**](https://github.com/gexijin/idepGolem), which runs locally as a second
container — no internet access is needed at analysis time.

### Starting iDEP

iDEP is defined in `docker-compose.yml`, so it comes up with the rest of the stack:

```bash
docker compose up -d          # both services
docker compose up -d idep     # just iDEP
```

Then open **http://localhost:3839** (SalmonFlow stays on 3838), or use the
**Open iDEP** button in the Results tab.

`run.sh` and `run.ps1` also handle iDEP — **once its image is pulled**, they
start it automatically, put both containers on a shared network, and report
its state at startup:

```
  SalmonFlow
  FASTQs:     /home/you/salmonflow/data/input
  References: /home/you/salmonflow/data/references
  Output:     /home/you/salmonflow/data/output
  iDEP:       started on http://localhost:3839
```

Pull the image once and the quick-start flow is unchanged from before:

```bash
docker pull gexijin/idep:latest   # one time, large
./run.sh                          # starts both
```

The scripts never trigger that download themselves — if the image is absent
they print `iDEP: not installed` and start SalmonFlow alone. Set
`SALMONFLOW_NO_IDEP=1` to skip iDEP even when the image is present, or
`IDEP_PORT=…` to publish it somewhere other than 3839.

iDEP is left running when SalmonFlow exits, so you can finish a session in the
browser. Stop it with `docker rm -f idep`.

### Requirements

The iDEP image bundles its own annotation database for 220 plant and animal
genomes, so it is **large** — iDEP's own documentation gives 10 GB of storage
and 4 GB of RAM as the minimum, and the published image is considerably bigger
than SalmonFlow's. The first `docker compose up` will take a while. The compose
file sets `pull_policy: missing` so it is downloaded only once.

### Using it

The Results tab provides two purpose-built exports:

| File | Contents |
|------|----------|
| `salmonflow_idep_counts_*.csv` | The count matrix with Ensembl version suffixes removed (`ENSG00000000003.16` → `ENSG00000000003`) and values rounded to whole counts, as DESeq2 and edgeR expect. Duplicate gene IDs created by stripping versions are collapsed by summing. |
| `salmonflow_idep_design_*.csv` | One row per sample with its **group**, taken from the group column of the Samples tab. Fill that column in before running if you want it populated. |

Upload the count matrix first, then the design file, in iDEP's **Load Data**
tab. The original `merged_lengthScaledTPM.csv` is unchanged and remains
available through the existing CSV/TSV buttons.

### Licence

iDEP is distributed under **CC BY-NC 3.0** — free for non-commercial use.
Its authors ask to be contacted before local installation at private
institutions. See the [iDEP repository](https://github.com/gexijin/idepGolem)
and cite [Ge, Son & Yao (2018), *BMC Bioinformatics*](https://doi.org/10.1186/s12859-018-2486-6)
if you use it. SalmonFlow only references the published image on Docker Hub;
it does not redistribute it.
