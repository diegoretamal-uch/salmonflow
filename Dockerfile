FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL maintainer="SalmonFlow" \
      description="Dockerized R Shiny app for bulk RNA-seq analysis"

# ── Extra system dependencies (wget/unzip for Salmon & FastQC) ───────
RUN apt-get update && apt-get install -y \
    wget curl unzip fastp \
    && rm -rf /var/lib/apt/lists/*

# ── Salmon 1.10.0 ───────────────────────────────────────────────────
RUN wget -q https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz \
    && mkdir -p /opt/salmon \
    && tar -xzf salmon-1.10.0_linux_x86_64.tar.gz -C /opt/salmon --strip-components=1 \
    && ln -s /opt/salmon/bin/salmon /usr/local/bin/salmon \
    && rm salmon-1.10.0_linux_x86_64.tar.gz

# ── FastQC 0.12.1 ───────────────────────────────────────────────────
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip \
    && unzip -q fastqc_v0.12.1.zip \
    && chmod +x FastQC/fastqc \
    && mv FastQC /opt/fastqc \
    && ln -s /opt/fastqc/fastqc /usr/local/bin/fastqc \
    && rm fastqc_v0.12.1.zip

# ── MultiQC ─────────────────────────────────────────────────────────
RUN pip3 install multiqc

# ── CRAN packages (install2.r exits non-zero on failure) ─────────────
RUN install2.r --error --skipinstalled --ncpus -1 \
    shiny shinydashboard shinyjs shinyFiles \
    DT ggplot2 plotly dplyr readr tidyr \
    processx future promises waiter \
    pheatmap jsonlite RColorBrewer \
    && rm -rf /tmp/downloaded_packages

# ── Bioconductor packages (base infra already in the image) ──────────
RUN R -e "BiocManager::install(c('tximport', 'GenomicFeatures', 'txdbmaker'), \
      ask=FALSE, update=FALSE); \
    pkgs <- c('tximport','GenomicFeatures','txdbmaker'); \
    missing <- pkgs[!pkgs %in% rownames(installed.packages())]; \
    if (length(missing)) stop(paste('Bioc install failed:', paste(missing, collapse=', ')))"

# ── Copy application ────────────────────────────────────────────────
COPY app/ /srv/shiny-server/salmonflow/

# ── Non-root runtime user ───────────────────────────────────────────
# Files written into bind-mounted volumes inherit the container process UID.
# Running as UID/GID 1000 (the typical first host user) keeps pipeline output
# owned by the host user instead of root. Override at build time to match a
# different host user:  docker build --build-arg UID=$(id -u) --build-arg GID=$(id -g)
ARG UID=1000
ARG GID=1000
RUN if ! getent group ${GID} >/dev/null; then groupadd -g ${GID} appuser; fi \
    && if ! getent passwd ${UID} >/dev/null; then \
         useradd -u ${UID} -g ${GID} -m -d /home/appuser -s /bin/bash appuser; \
       fi \
    && mkdir -p /home/appuser /data/input /data/references /data/output /data/tmp \
    && chown -R ${UID}:${GID} /home/appuser /data /srv/shiny-server/salmonflow

# R (and MultiQC) need a writable HOME for their user library / config / cache.
ENV HOME=/home/appuser
USER ${UID}:${GID}

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/salmonflow', host='0.0.0.0', port=3838)"]
