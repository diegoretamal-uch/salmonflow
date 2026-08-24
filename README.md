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

---

## Updating an Existing Install

Already running SalmonFlow from an earlier build? **Do not clone again** —
update in place. Only two steps differ from a fresh install: you stop the old
container first, and you `git pull` instead of `git clone`. Steps 3-5 are
identical either way.

### Linux / macOS

```bash
cd salmonflow

# 1. Stop the running instance (it holds the app's port)
docker rm -f salmonflow 2>/dev/null
#    "No such container"? See "Finding an unnamed container" below.

# 2. Get the new code
git pull

# 3. Rebuild — fast, only the app layer changed
docker build -t salmonflow .

# 4. One-time, if you want iDEP
docker pull gexijin/idep:latest

# 5. Run as usual
./run.sh
```

### Windows (PowerShell)

```powershell
cd salmonflow

docker rm -f salmonflow 2>$null
git pull
docker build -t salmonflow .
docker pull gexijin/idep:latest
.\run.ps1
```

**Why the rebuild is fast.** The `Dockerfile` itself rarely changes, so Docker
reuses every cached layer and only re-copies `app/`. Seconds, not the ~10
minutes of a first build — R and Bioconductor are not reinstalled.

**Why step 1 matters.** An old container keeps the port bound, and the new one
refuses to start with `Bind for 0.0.0.0:3838 failed: port is already
allocated`.

### Finding an unnamed container

`run.sh` now starts its container as `salmonflow`, so `docker rm -f salmonflow`
is all you need. But instances started by an **older `run.sh`**, or by a bare
`docker run`, got a random name like `elastic_booth`. If step 1 said
`No such container`, list what is running:

```bash
docker ps
```

```
CONTAINER ID   IMAGE                 PORTS                        NAMES
741e50aa2d16   salmonflow            0.0.0.0:3838->3838/tcp       vigorous_newton   <- yours
3f8ec2d4c107   gexijin/idep:latest   0.0.0.0:3839->3838/tcp       idep              <- leave it
```

Find the row whose **IMAGE** is `salmonflow`, then remove it by the name in the
last column:

```bash
docker rm -f vigorous_newton      # use the name YOU see, not this one
```

Match on the **IMAGE** column, not on the ports: both containers show `3838` on
the right of their mapping, because that is the port *inside* every container.
Only the `salmonflow` one should be removed — deleting `idep` just means
pulling or restarting it again later.

> For the same reason, do **not** hunt for it with
> `docker ps --filter expose=3838` — that matches iDEP as well and would
> remove it too.

**Your data is safe.** Results live in bind-mounted host folders and no file
format changed — nothing to migrate.

**Old results won't appear automatically.** The Results tab reads the count
matrix from the *current session*, so a freshly started app shows "Run the
pipeline to enable the iDEP exports" even when `data/output/` is full. Use
**Resume Analysis** in the Run tab — it skips every completed step and
repopulates Results in a couple of minutes rather than re-running everything.

---

## Ports

SalmonFlow uses **3838**, iDEP uses **3839**. Both are host-side choices you
can change; inside the containers nothing moves.

### Checking what is already in use

```bash
# Are our two ports taken? Empty output = both free.
ss -tln 'sport = :3838 or sport = :3839'

# Who owns them? (sudo required — Docker's listeners belong to root)
sudo ss -tlnp 'sport = :3838 or sport = :3839'

# Which of them are ours?
docker ps --format 'table {{.Names}}\t{{.Image}}\t{{.Ports}}'
```

macOS has no `ss`; use `lsof -nP -iTCP:3838 -sTCP:LISTEN` instead.

```powershell
# Windows
Get-NetTCPConnection -LocalPort 3838,3839 -State Listen -ErrorAction SilentlyContinue
Get-Process -Id (Get-NetTCPConnection -LocalPort 3838).OwningProcess   # who owns it
```

Then read the result:

| What you see | Meaning | What to do |
|---|---|---|
| Nothing | Port is free | Carry on |
| Listed in **both** `ss` and `docker ps` | One of your own containers | `docker rm -f salmonflow` and reuse the port |
| In `ss` but **not** in `docker ps` | Another program owns it — RStudio Server and Shiny Server both default to 3838 | Use a different port (below) |

> Avoid `ss -tln | grep :3838`. The pattern also matches `:38380`, `:38381`
> and so on, so an unrelated service can look like a conflict. The
> `sport = :3838` filter above matches the port exactly.

### Running on different ports

```bash
SALMONFLOW_PORT=8080 IDEP_PORT=8081 ./run.sh
# → app at http://localhost:8080, iDEP at http://localhost:8081
```

```powershell
$env:SALMONFLOW_PORT=8080; $env:IDEP_PORT=8081; .\run.ps1
```

To make it permanent, put the exports in your shell profile, or use
`docker-compose.yml` and edit the `ports:` lines there.

**The iDEP button keeps working.** The Results tab builds the iDEP URL from
your browser's current hostname plus `IDEP_PORT`, so the two ports are
independent — SalmonFlow on 8080 with iDEP on 3839 is fine. What matters is
that `IDEP_PORT` matches the port iDEP is actually published on, which
`run.sh` handles by passing the same value to both.

> **Only the left number is yours to change.** In `-p 8080:3838` the left side
> is the host port and the right side is the container's, which is always
> 3838. Editing the right side breaks the app.

---

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
limma-voom, limma-trend), clustering and pathway enrichment, the Results tab
hands off to
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
| `salmonflow_idep_counts_*.csv` | The count matrix with Ensembl version suffixes removed (`ENSG00000000003.16` → `ENSG00000000003`) and values rounded to whole counts, as DESeq2 expects. Duplicate gene IDs created by stripping versions are collapsed by summing. |
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
