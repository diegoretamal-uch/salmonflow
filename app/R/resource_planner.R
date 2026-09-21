# ══════════════════════════════════════════════════════════════
# SalmonFlow — resource_planner.R
# Detect the real CPU/RAM/disk budget (host + cgroup limits) and split it
# across tools that scale very differently, so no stage is handed a thread
# count it cannot use.
#
# Detection is deliberately runtime-only: nothing here assumes a particular
# machine. A container CPU quota, a cluster allocation and a bare-metal host
# all answer the same questions.
# ══════════════════════════════════════════════════════════════

#' Available RAM in bytes, honouring a container/cgroup memory limit.
#' Returns min(host MemAvailable, cgroup headroom). NA if nothing detectable.
detect_available_ram <- function() {
  host <- tryCatch({
    mi   <- readLines("/proc/meminfo", warn = FALSE)
    line <- grep("^MemAvailable:", mi, value = TRUE)
    if (length(line) == 0) NA_real_
    else as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", line)) * 1024
  }, error = function(e) NA_real_)

  # cgroup limit minus current usage = real headroom inside the container.
  cg <- tryCatch({
    limit <- NA_real_; used <- 0
    if (file.exists("/sys/fs/cgroup/memory.max")) {                    # cgroup v2
      lm <- readLines("/sys/fs/cgroup/memory.max", warn = FALSE)[1]
      if (!identical(lm, "max")) limit <- suppressWarnings(as.numeric(lm))
      if (file.exists("/sys/fs/cgroup/memory.current")) {
        used <- suppressWarnings(as.numeric(
          readLines("/sys/fs/cgroup/memory.current", warn = FALSE)[1]))
      }
    } else if (file.exists("/sys/fs/cgroup/memory/memory.limit_in_bytes")) {  # cgroup v1
      lm <- suppressWarnings(as.numeric(
        readLines("/sys/fs/cgroup/memory/memory.limit_in_bytes", warn = FALSE)[1]))
      # v1 reports a near-INT64_MAX sentinel when unlimited.
      if (!is.na(lm) && lm < 1e15) limit <- lm
      uf <- "/sys/fs/cgroup/memory/memory.usage_in_bytes"
      if (file.exists(uf)) {
        used <- suppressWarnings(as.numeric(readLines(uf, warn = FALSE)[1]))
      }
    }
    if (is.na(limit)) NA_real_ else max(0, limit - (if (is.na(used)) 0 else used))
  }, error = function(e) NA_real_)

  vals <- c(host, cg)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) NA_real_ else min(vals)
}

#' Usable CPU count, honouring a container/cgroup CPU quota.
#' parallel::detectCores() alone reports logical SMT threads and is blind to
#' `docker run --cpus` or a scheduler's allocation, so take the min of both.
#' Returns NA if nothing is detectable.
detect_available_cores <- function() {
  host <- tryCatch(parallel::detectCores(), error = function(e) NA_real_)

  cg <- tryCatch({
    if (file.exists("/sys/fs/cgroup/cpu.max")) {                       # cgroup v2
      parts <- strsplit(trimws(readLines("/sys/fs/cgroup/cpu.max", warn = FALSE)[1]),
                        "\\s+")[[1]]
      if (length(parts) == 2 && parts[1] != "max") {
        ceiling(as.numeric(parts[1]) / as.numeric(parts[2]))
      } else NA_real_
    } else if (file.exists("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")) {   # cgroup v1
      q  <- suppressWarnings(as.numeric(
        readLines("/sys/fs/cgroup/cpu/cpu.cfs_quota_us", warn = FALSE)[1]))
      pe <- suppressWarnings(as.numeric(
        readLines("/sys/fs/cgroup/cpu/cpu.cfs_period_us", warn = FALSE)[1]))
      if (!is.na(q) && q > 0 && !is.na(pe) && pe > 0) ceiling(q / pe) else NA_real_
    } else NA_real_
  }, error = function(e) NA_real_)

  vals <- c(host, cg)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) NA_integer_ else as.integer(max(1, floor(min(vals))))
}

#' Free bytes on the filesystem holding `path`.
#' Walks up to the nearest existing ancestor, since the output directory may
#' not have been created yet when a projection is requested.
detect_free_disk <- function(path) {
  tryCatch({
    p <- path
    while (nchar(p) > 1 && !dir.exists(p)) p <- dirname(p)
    out <- suppressWarnings(
      system2("df", args = c("-Pk", shQuote(p)), stdout = TRUE, stderr = FALSE))
    if (length(out) < 2) return(NA_real_)
    fields <- strsplit(trimws(out[length(out)]), "\\s+")[[1]]
    if (length(fields) < 4) return(NA_real_)
    as.numeric(fields[4]) * 1024
  }, error = function(e) NA_real_)
}

# ── Per-tool thread budget ────────────────────────────────────
# One integer cannot serve four tools that scale differently:
#   FastQC  — exactly one core per file; --threads IS the file count, so any
#             value above the number of files is wasted. It also sets the JVM
#             heap (memory x threads), so an inflated value inflates the heap.
#   fastp   — measured plateau at 4-8 worker threads; hard-capped at 16.
#   salmon  — the only tool that genuinely uses the full budget.

FASTQC_MB_PER_THREAD <- 512L   # the fastqc launcher's own default
FASTQC_MB_MINIMUM    <- 250L   # per FastQC's docs, the memory one file needs
FASTP_MAX_USEFUL     <- 8L     # measured: flat past 8, slightly worse at 16

#' Threads for a FastQC invocation: never more than the files handed to it.
fastqc_threads <- function(n_files, budget) {
  max(1L, min(as.integer(n_files), as.integer(budget)))
}

#' Pick a thread count AND heap size for FastQC that fit inside a RAM budget.
#'
#' The launcher computes -Xmx(memory x threads)m, so the two are coupled. When
#' memory is tight the lever is the THREAD COUNT, not the heap per thread:
#' starving the heap below what one file needs makes FastQC die with an
#' OutOfMemoryError rather than merely run lean. So drop threads until
#' FASTQC_MB_MINIMUM each fits, and only then trim the per-thread figure.
#'
#' @param fraction Share of available RAM FastQC may claim.
#' @return list(threads, memory_mb); memory_mb is NA to mean "launcher default".
fastqc_plan <- function(n_files, budget, ram_bytes = NA_real_, fraction = 0.25) {
  threads <- fastqc_threads(n_files, budget)

  if (is.na(ram_bytes) || ram_bytes <= 0) {
    return(list(threads = threads, memory_mb = NA_integer_))
  }

  budget_mb <- (ram_bytes * fraction) / 1024^2

  # Step 1: cut threads until each can hold the documented minimum.
  max_by_ram <- max(1L, as.integer(floor(budget_mb / FASTQC_MB_MINIMUM)))
  threads    <- max(1L, min(threads, max_by_ram))

  # Step 2: only cap the heap if the default would overshoot the budget.
  per <- floor(budget_mb / threads)
  memory_mb <- if (per >= FASTQC_MB_PER_THREAD) {
    NA_integer_                                   # default is already fine
  } else {
    # The launcher rejects anything outside 100-10000 MB.
    as.integer(max(FASTQC_MB_MINIMUM, min(10000L, per)))
  }

  list(threads = threads, memory_mb = memory_mb)
}

#' Threads for fastp: capped where the measured scaling curve goes flat.
fastp_threads <- function(budget) {
  max(1L, min(as.integer(budget), FASTP_MAX_USEFUL))
}

#' Threads for salmon quant, leaving room for a concurrently running FastQC.
salmon_quant_threads <- function(budget, reserved = 0L) {
  max(1L, as.integer(budget) - as.integer(reserved))
}

# ── Launch-time storage projection ────────────────────────────
# The storage policy already exists (mod_params.R); what was missing is that a
# user could not see what their choice costs until the run was under way.
#
# R = total raw bytes, N = samples, t = per-sample trimmed bytes, W = samples
# in flight. Trimming is roughly size-neutral (fastp recompresses at -z 4, so
# output is typically a little LARGER than input), hence TRIMMED_RATIO.
TRIMMED_RATIO <- 1.09

#' Total bytes of the raw FASTQs for a sample table (columns r1, r2).
raw_bytes_of <- function(samples) {
  if (is.null(samples) || nrow(samples) == 0) return(0)
  paths <- c(samples$r1, samples$r2)
  paths <- paths[!is.na(paths) & nchar(paths) > 0]
  if (length(paths) == 0) return(0)
  sz <- file.info(paths)$size
  sum(sz[!is.na(sz)])
}

#' Project peak SAMPLE-FILE bytes for a run.
#'
#' Only sample FASTQs are counted — indexes, QC reports and caches are
#' auxiliary and deliberately excluded.
#'
#' @param policy "keep_both", "delete_raw" or "delete_trimmed"
#' @param concurrency Samples in flight. Only "delete_trimmed" is sensitive to
#'   it: under the other two, trimmed files either accumulate regardless or are
#'   funded by freed raws, so overlapping samples costs no extra disk.
#' @return list(raw, peak, limiting) in bytes.
project_peak_storage <- function(samples, policy = "keep_both",
                                 trimming_enabled = TRUE, concurrency = 1L) {
  R <- raw_bytes_of(samples)
  N <- if (is.null(samples)) 0 else nrow(samples)

  if (!isTRUE(trimming_enabled) || N == 0) {
    return(list(raw = R, peak = R, limiting = "no trimmed files"))
  }

  t_total <- R * TRIMMED_RATIO
  t_per   <- t_total / N
  W       <- max(1L, min(as.integer(concurrency), N))

  switch(policy,
    "delete_raw" = list(
      raw = R,
      # Each raw is freed as its trimmed appears, so the total drifts from R
      # up to t_total; the W in-flight sets are already counted in that drift.
      peak = max(R, t_total) + (W - 1) * t_per,
      limiting = "trimmed files accumulate; raws freed as they go"),
    "delete_trimmed" = list(
      raw = R,
      # The only policy where trimmed files are transient, so the number live
      # at once IS the concurrency.
      peak = R + W * t_per,
      limiting = if (W > 1) sprintf("raws kept + %d trimmed set(s) live", W)
                 else "raws kept + 1 trimmed set live"),
    list(  # keep_both (default)
      raw = R,
      peak = R + t_total,
      limiting = "raws and all trimmed files kept")
  )
}

#' Format a byte count for the UI.
fmt_bytes <- function(b) {
  if (is.null(b) || length(b) == 0 || is.na(b)) return("unknown")
  u <- c("B", "KB", "MB", "GB", "TB"); i <- 1
  while (b >= 1024 && i < length(u)) { b <- b / 1024; i <- i + 1 }
  sprintf(if (b >= 100 || i <= 2) "%.0f %s" else "%.1f %s", b, u[i])
}
