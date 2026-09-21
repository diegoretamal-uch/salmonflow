# ══════════════════════════════════════════════════════════════
# SalmonFlow — cpu_monitor.R
# Record how much CPU (and RAM) each pipeline stage actually used, so the
# thread budget handed out by resource_planner.R can be checked against
# reality instead of assumed.
#
# Why an EXTERNAL sampler: processx::run() blocks the runner's single R
# thread for the whole duration of a tool, so no R-level timer can sample
# while salmon or fastp is working. A tiny bash loop started once per run
# writes a CSV of raw monotonic counters; each stage is then bracketed by
# timestamps and its numbers derived from the rows inside that window.
#
# What the numbers mean: the sampler reads the CGROUP counters when they
# exist, so the figure is the container's own usage. Where no cgroup is
# visible it falls back to /proc/stat, which is HOST-WIDE — those lines are
# labelled as such, because anything else on the machine is counted too.
# Windows can also overlap (FastQC runs alongside salmon quant by design);
# overlapping stages are marked and their shared CPU is counted in both.
# ══════════════════════════════════════════════════════════════

# 2s costs nothing (two small file reads and a printf) and keeps stages of a
# few seconds measurable; a 4-hour run still writes well under 1 MB of CSV.
CPU_SAMPLE_INTERVAL <- 2L    # seconds between samples
CPU_MIN_BLOCK_SECS  <- 5L    # shorter stages cannot be measured meaningfully

.cpu_state <- new.env(parent = emptyenv())
.cpu_state$proc   <- NULL
.cpu_state$csv    <- NULL
.cpu_state$open   <- list()   # currently open blocks, for overlap detection
.cpu_state$seq    <- 0L
.cpu_state$blocks <- list()   # finished records, for the run summary
.cpu_state$interval <- CPU_SAMPLE_INTERVAL
.cpu_state$scoped   <- NULL   # memoised "counters are ours alone" answer

# Marker for a stage that shared its window with another. Deliberately ASCII:
# R's sprintf pads %-32s by BYTES, so a multi-byte glyph here shifts the
# summary table's columns out of line.
CPU_OVERLAP_MARK <- "*"

#' The sampler loop, run by bash.
#'
#' Emits cumulative counters (never rates) and lets R do the deltas, which
#' keeps the shell trivial. It appends with >> on every tick rather than
#' holding stdout open, because a block-buffered stdout would withhold rows
#' for minutes and make short windows unmeasurable.
cpu_sampler_script <- function() {
  r"---(
set -u
out=$1
interval=$2
# The runner may be SIGKILLed (the Cancel button), leaving no chance to reap
# this loop from R. Watch the parent and exit with it rather than appending
# to the CSV forever.
ppid=$PPID
clk=$(getconf CLK_TCK 2>/dev/null || echo 100)
src=proc
v1dir=""
if [ -r /sys/fs/cgroup/cpu.stat ] && grep -q usage_usec /sys/fs/cgroup/cpu.stat 2>/dev/null; then
  src=cgroup2
else
  for d in /sys/fs/cgroup/cpuacct /sys/fs/cgroup/cpu,cpuacct; do
    if [ -r "$d/cpuacct.usage" ]; then src=cgroup1; v1dir=$d; break; fi
  done
fi
v1mem=/sys/fs/cgroup/memory/memory.usage_in_bytes

# Used bytes from /proc/meminfo. The fallback for both cgroup branches too:
# a root cgroup exposes no memory.current, and a silent 0 would be reported
# as a peak RAM of nothing.
readmem_proc() {
  mt=0
  ma=0
  while read -r k v _u; do
    case $k in
      MemTotal:) mt=$v ;;
      MemAvailable:) ma=$v; break ;;
    esac
  done < /proc/meminfo
  mem=$(( (mt - ma) * 1024 ))
}

while :; do
  kill -0 "$ppid" 2>/dev/null || exit 0
  cpu=0
  mem=0
  case $src in
    cgroup2)
      while read -r k v; do
        if [ "$k" = usage_usec ]; then cpu=$v; break; fi
      done < /sys/fs/cgroup/cpu.stat
      if [ -r /sys/fs/cgroup/memory.current ]; then
        read -r mem < /sys/fs/cgroup/memory.current
      else
        readmem_proc
      fi
      ;;
    cgroup1)
      read -r ns < "$v1dir/cpuacct.usage"
      cpu=$(( ns / 1000 ))
      if [ -r $v1mem ]; then read -r mem < $v1mem; else readmem_proc; fi
      ;;
    *)
      read -r _lbl u n s _i _w q sq st _rest < /proc/stat
      cpu=$(( (u + n + s + q + sq + st) * 1000000 / clk ))
      readmem_proc
      ;;
  esac
  [ -n "$mem" ] || mem=0
  printf '%s,%s,%s,%s\n' "${EPOCHSECONDS:-$(date +%s)}" "$cpu" "$mem" "$src" >> "$out"
  sleep "$interval"
done
)---"
}

#' Start the background sampler. Never fails the run: on any error the
#' monitor simply stays off and no CPU lines are produced. The loop also
#' watches its parent, so a killed runner takes it down with it.
#' @return TRUE if sampling started.
start_cpu_sampler <- function(csv_path, interval = CPU_SAMPLE_INTERVAL,
                              log_callback = NULL) {
  ok <- tryCatch({
    dir.create(dirname(csv_path), showWarnings = FALSE, recursive = TRUE)
    writeLines("epoch,cpu_usec,mem_bytes,source", csv_path)
    proc <- processx::process$new(
      "bash",
      args    = c("-c", cpu_sampler_script(), "sf_cpu_sampler",
                  csv_path, as.character(as.integer(interval))),
      stdout  = NULL, stderr = NULL, cleanup = TRUE)
    .cpu_state$proc     <- proc
    .cpu_state$csv      <- csv_path
    .cpu_state$interval <- as.integer(interval)
    .cpu_state$blocks   <- list()
    .cpu_state$open     <- list()
    TRUE
  }, error = function(e) {
    if (!is.null(log_callback)) {
      log_callback(paste("CPU monitor: could not start —", conditionMessage(e)), "warn")
    }
    FALSE
  })

  if (isTRUE(ok) && !is.null(log_callback)) {
    log_callback(sprintf("CPU monitor: sampling every %ds → %s",
                         as.integer(interval), csv_path), "info")
  }
  isTRUE(ok)
}

#' Stop the sampler. Safe to call more than once, and when it never started.
stop_cpu_sampler <- function(log_callback = NULL) {
  if (is.null(.cpu_state$proc)) return(invisible(NULL))
  tryCatch({
    if (.cpu_state$proc$is_alive()) .cpu_state$proc$kill()
  }, error = function(e) NULL)
  .cpu_state$proc <- NULL
  invisible(NULL)
}

#' Do the cgroup counters describe THIS container alone?
#' A bare host's root cgroup answers cpu.stat just like /proc/stat does, so a
#' cgroup reading is only our own when we are in a container, or in a cgroup
#' that carries a cpu/memory limit (a systemd slice, a scheduler allocation).
#' Memoised: the answer cannot change during a run.
cpu_scoped <- function() {
  if (!is.null(.cpu_state$scoped)) return(.cpu_state$scoped)
  res <- tryCatch({
    # In a container the counters are ours whether or not a limit is set.
    lim <- file.exists("/.dockerenv") || file.exists("/run/.containerenv")
    if (!lim && file.exists("/sys/fs/cgroup/cpu.max")) {               # cgroup v2
      q <- strsplit(trimws(readLines("/sys/fs/cgroup/cpu.max", warn = FALSE)[1]),
                    "\\s+")[[1]][1]
      if (!identical(q, "max")) lim <- TRUE
    }
    if (!lim && file.exists("/sys/fs/cgroup/memory.max")) {
      if (!identical(trimws(readLines("/sys/fs/cgroup/memory.max", warn = FALSE)[1]),
                     "max")) lim <- TRUE
    }
    if (!lim && file.exists("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")) {  # cgroup v1
      q <- suppressWarnings(as.numeric(
        readLines("/sys/fs/cgroup/cpu/cpu.cfs_quota_us", warn = FALSE)[1]))
      if (!is.na(q) && q > 0) lim <- TRUE
    }
    if (!lim && file.exists("/sys/fs/cgroup/memory/memory.limit_in_bytes")) {
      lm <- suppressWarnings(as.numeric(
        readLines("/sys/fs/cgroup/memory/memory.limit_in_bytes", warn = FALSE)[1]))
      if (!is.na(lm) && lm < 1e15) lim <- TRUE   # v1 sentinel for "unlimited"
    }
    lim
  }, error = function(e) FALSE)
  .cpu_state$scoped <- res
  res
}

#' TRUE when the figures cover more than this container's own work.
cpu_host_wide <- function(src) identical(src, "proc") || !cpu_scoped()

#' Read the timeline, dropping any row the sampler was still writing.
cpu_read_samples <- function() {
  csv <- .cpu_state$csv
  if (is.null(csv) || !file.exists(csv)) return(NULL)
  lines <- tryCatch(readLines(csv, warn = FALSE), error = function(e) character(0))
  lines <- lines[grepl("^[0-9]+,[0-9]+,[0-9]+,[a-z0-9]+$", lines)]
  if (length(lines) < 2) return(NULL)
  parts <- strsplit(lines, ",", fixed = TRUE)
  data.frame(
    epoch  = as.numeric(vapply(parts, `[`, "", 1L)),
    cpu    = as.numeric(vapply(parts, `[`, "", 2L)),
    mem    = as.numeric(vapply(parts, `[`, "", 3L)),
    source = vapply(parts, `[`, "", 4L),
    stringsAsFactors = FALSE
  )
}

#' Derive usage over [t0, t1] from the counter deltas inside that window.
#' @return list(mean_cores, peak_cores, peak_mem, n, source), or NULL when
#'   fewer than two samples fall inside the window.
cpu_window_stats <- function(t0, t1) {
  s <- cpu_read_samples()
  if (is.null(s)) return(NULL)
  s <- s[s$epoch >= as.numeric(t0) & s$epoch <= as.numeric(t1), , drop = FALSE]
  if (nrow(s) < 2) return(NULL)

  span <- s$epoch[nrow(s)] - s$epoch[1]
  if (span <= 0) return(NULL)

  d_cpu <- diff(s$cpu)
  d_t   <- diff(s$epoch)
  ok    <- d_t > 0 & d_cpu >= 0          # guards a counter reset
  inst  <- if (any(ok)) (d_cpu[ok] / 1e6) / d_t[ok] else numeric(0)

  list(
    mean_cores = (s$cpu[nrow(s)] - s$cpu[1]) / 1e6 / span,
    peak_cores = if (length(inst) > 0) max(inst) else NA_real_,
    peak_mem   = suppressWarnings(max(s$mem, na.rm = TRUE)),
    n          = nrow(s),
    source     = s$source[1]
  )
}

#' Open a measurement window.
#' @param allocated Threads this stage was given, for the efficiency figure.
#' @return A token for cpu_block_end(), or NULL when monitoring is off.
cpu_block_start <- function(label, allocated = NA_integer_) {
  if (is.null(.cpu_state$csv)) return(NULL)
  .cpu_state$seq <- .cpu_state$seq + 1L
  id <- as.character(.cpu_state$seq)

  env <- new.env(parent = emptyenv())
  env$concurrent <- FALSE
  .cpu_state$open[[id]] <- env

  # Any stage open at the same time as another shares its window; mark them
  # all so no line is read as belonging to one tool alone.
  if (length(.cpu_state$open) > 1L) {
    for (o in .cpu_state$open) o$concurrent <- TRUE
  }

  list(id = id, label = label, allocated = allocated, t0 = Sys.time(), env = env)
}

#' Close a window, record it for the summary and log one line.
#' Safe with NULL (monitoring off, or a stage that was skipped on resume).
cpu_block_end <- function(token, log_callback = NULL) {
  if (is.null(token)) return(invisible(NULL))
  t1 <- Sys.time()
  .cpu_state$open[[token$id]] <- NULL

  elapsed <- as.numeric(difftime(t1, token$t0, units = "secs"))
  stats   <- if (elapsed >= CPU_MIN_BLOCK_SECS) cpu_window_stats(token$t0, t1) else NULL

  rec <- list(label      = token$label,
              allocated  = token$allocated,
              elapsed    = elapsed,
              concurrent = isTRUE(token$env$concurrent),
              stats      = stats)
  .cpu_state$blocks[[length(.cpu_state$blocks) + 1L]] <- rec

  if (!is.null(stats) && !is.null(log_callback)) log_callback(cpu_block_line(rec), "info")
  invisible(stats)
}

#' Bracket an expression with a measurement window, returning its value.
with_cpu_block <- function(label, allocated = NA_integer_, log_callback = NULL, expr) {
  token <- cpu_block_start(label, allocated)
  on.exit(cpu_block_end(token, log_callback), add = TRUE)
  expr
}

cpu_fmt_dur <- function(secs) {
  secs <- as.numeric(secs)
  if (!is.finite(secs)) return("n/a")
  if (secs < 60) return(sprintf("%ds", round(secs)))
  sprintf("%dm %ds", floor(secs / 60), round(secs %% 60))
}

cpu_use_pct <- function(mean_cores, allocated) {
  if (is.null(allocated) || length(allocated) == 0 ||
      is.na(allocated) || allocated <= 0) return(NA_integer_)
  as.integer(round(mean_cores / allocated * 100))
}

#' One live-log line for a finished stage.
cpu_block_line <- function(rec) {
  s   <- rec$stats
  pct <- cpu_use_pct(s$mean_cores, rec$allocated)
  paste0(
    "CPU  ", if (isTRUE(rec$concurrent)) paste0(CPU_OVERLAP_MARK, " ") else "",
    rec$label, ": ",
    sprintf("%.1f cores avg", s$mean_cores),
    if (!is.na(s$peak_cores)) sprintf(", %.1f peak", s$peak_cores) else "",
    if (!is.na(pct)) sprintf(", of %d allocated (%d%%)",
                             as.integer(rec$allocated), pct) else "",
    ", peak RAM ", fmt_bytes(s$peak_mem),
    ", ", cpu_fmt_dur(rec$elapsed),
    if (cpu_host_wide(s$source)) " (host-wide)" else "")
}

#' The "-- Resource usage --" section of run_summary.log.
#' Returns character(0) when nothing was measured.
cpu_summary_lines <- function(run_start, run_end = Sys.time()) {
  if (is.null(.cpu_state$csv)) return(character(0))
  if (length(.cpu_state$blocks) == 0L) return(character(0))

  overall <- cpu_window_stats(run_start, run_end)
  src     <- if (!is.null(overall)) overall$source else "proc"
  host_wide <- cpu_host_wide(src)
  src_lbl   <- paste0(
    switch(src, cgroup2 = "cgroup v2 counters",
                cgroup1 = "cgroup v1 counters",
                "/proc/stat counters"),
    if (host_wide) " — HOST-WIDE" else " — container-scoped")

  L <- c("", "-- Resource usage --",
         paste0("Source: ", src_lbl,
                sprintf(" | sampled every %ds",
                        as.integer(.cpu_state$interval %||% CPU_SAMPLE_INTERVAL))))

  if (!is.null(overall)) {
    L <- c(L, sprintf("Run-wide: %.1f cores avg, %.1f peak | peak RAM %s",
                      overall$mean_cores,
                      if (is.na(overall$peak_cores)) 0 else overall$peak_cores,
                      fmt_bytes(overall$peak_mem)))
  }

  L <- c(L, "",
         sprintf("  %-32s %10s %10s %7s %7s %6s",
                 "Stage", "Duration", "Avg cores", "Peak", "Alloc", "Use"))

  any_concurrent <- FALSE
  for (rec in .cpu_state$blocks) {
    if (isTRUE(rec$concurrent)) any_concurrent <- TRUE
    lbl <- paste0(if (isTRUE(rec$concurrent)) paste0(CPU_OVERLAP_MARK, " ") else "",
                  rec$label)
    s   <- rec$stats
    if (is.null(s)) {
      L <- c(L, sprintf("  %-32s %10s %10s %7s %7s %6s",
                        lbl, cpu_fmt_dur(rec$elapsed), "n/a", "-", "-", "-"))
    } else {
      pct <- cpu_use_pct(s$mean_cores, rec$allocated)
      L <- c(L, sprintf("  %-32s %10s %10.1f %7s %7s %6s",
                        lbl, cpu_fmt_dur(rec$elapsed), s$mean_cores,
                        if (is.na(s$peak_cores)) "-" else sprintf("%.1f", s$peak_cores),
                        if (is.na(rec$allocated)) "-" else as.character(as.integer(rec$allocated)),
                        if (is.na(pct)) "-" else paste0(pct, "%")))
    }
  }

  if (any_concurrent) {
    L <- c(L, "",
           paste0("  ", CPU_OVERLAP_MARK,
                  " ran at the same time as another stage (FastQC overlaps"),
           "    quant by design); the CPU they shared is counted in both rows.")
  }
  if (host_wide) {
    L <- c(L, "  No cgroup limit was visible, so these figures cover the whole",
              "  machine, not SalmonFlow alone.")
  }
  L
}
