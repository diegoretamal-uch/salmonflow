FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL maintainer="SalmonFlow" \
      description="Dockerized R Shiny app for bulk RNA-seq analysis"

# ── Extra system dependencies (wget/unzip for Salmon, FastQC & fastp) ─
RUN apt-get update && apt-get install -y \
    wget curl unzip \
    && rm -rf /var/lib/apt/lists/*

# ── fastp 1.3.6 (upstream static binary, checksum-verified) ─────────
RUN wget -q -O /usr/local/bin/fastp http://opengene.org/fastp/fastp.1.3.6 \
    && echo "817df51647ecdacc1642daabe2438c744828ccda2d3634d90395ea91e3a7bc1f  /usr/local/bin/fastp" | sha256sum -c - \
    && chmod +x /usr/local/bin/fastp \
    && fastp --version 2>&1 | grep -q "1.3.6"

# ── Salmon 2.11.0 (Rust rewrite, upstream binary, checksum-verified) ─
RUN wget -q https://github.com/COMBINE-lab/salmon/releases/download/v2.11.0/salmon-cli-x86_64-unknown-linux-gnu.tar.xz \
    && echo "e15b1ef0c4bece9c2ce271e2acbafb6ac1dc22b4b038fe8380f849455b9d8527  salmon-cli-x86_64-unknown-linux-gnu.tar.xz" | sha256sum -c - \
    && mkdir -p /opt/salmon \
    && tar -xJf salmon-cli-x86_64-unknown-linux-gnu.tar.xz -C /opt/salmon --strip-components=1 \
    && chmod +x /opt/salmon/salmon \
    && ln -s /opt/salmon/salmon /usr/local/bin/salmon \
    && rm salmon-cli-x86_64-unknown-linux-gnu.tar.xz \
    && salmon --version 2>&1 | grep -q "2.11.0"

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
