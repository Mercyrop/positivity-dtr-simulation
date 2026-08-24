#!/usr/bin/env Rscript

## Run the primary simulation grid.
##
## Usage:
##   Rscript scripts/02_run_primary.R [n_sim] [workers] [output_file]
##
## Manuscript settings:
##   n_sim   = 1000
##   n       = 500, 1500, 5000
##   scenario= benchmark, mild, moderate, severe
##   G-computation bootstrap CIs are run separately by 03_run_gcomp_bootstrap.R.

args <- commandArgs(trailingOnly = TRUE)
n_sim <- if (length(args) >= 1L) as.integer(args[1]) else 1000L
workers <- if (length(args) >= 2L) as.integer(args[2]) else 1L

locate_repo_root <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)

  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1L])
    script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
    return(normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
  }

  source_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(source_path) && length(source_path) && !is.na(source_path) && nzchar(source_path)) {
    script_dir <- dirname(normalizePath(source_path, mustWork = TRUE))
    return(normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
  }

  wd <- normalizePath(getwd(), mustWork = TRUE)
  if (dir.exists(file.path(wd, "R")) && dir.exists(file.path(wd, "scripts"))) {
    return(wd)
  }
  if (basename(wd) == "scripts" && dir.exists(file.path(dirname(wd), "R"))) {
    return(normalizePath(dirname(wd), mustWork = TRUE))
  }

  stop(
    "Could not locate the repository root. Run this script from the repository root ",
    "or execute it with Rscript."
  )
}

repo_root <- locate_repo_root()

output_file <- if (length(args) >= 3L) {
  args[3]
} else {
  file.path(repo_root, "results", "primary_simulation.rds")
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_root, "R", "01_data_generating_mechanism.R"))
source(file.path(repo_root, "R", "02_estimators.R"))

n_values <- c(500L, 1500L, 5000L)
scenarios <- c("benchmark", "mild", "moderate", "severe")

grid <- expand.grid(
  n = n_values,
  scenario = scenarios,
  sim = seq_len(n_sim),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid <- grid[order(grid$n, match(grid$scenario, scenarios), grid$sim), ]
grid$seed <- as.integer(mapply(
  make_seed_v7,
  n = grid$n,
  scenario = grid$scenario,
  sim = grid$sim
))

run_one <- function(i) {
  row <- grid[i, ]
  est <- tryCatch(
    one_sim_all_estimators_v7(
      n = row$n,
      seed = row$seed,
      scenario = row$scenario,
      gcomp_B = 0L
    ),
    error = function(e) {
      message(sprintf(
        "Replicate failed: n=%d scenario=%s sim=%d: %s",
        row$n, row$scenario, row$sim, conditionMessage(e)
      ))
      NULL
    }
  )

  if (is.null(est)) {
    return(data.frame(
      n = row$n, scenario = row$scenario, sim = row$sim, seed = row$seed,
      failed = TRUE, error_message = "replicate-level failure",
      stringsAsFactors = FALSE
    ))
  }

  cbind(
    data.frame(
      n = row$n, scenario = row$scenario, sim = row$sim, seed = row$seed,
      failed = FALSE, error_message = NA_character_,
      stringsAsFactors = FALSE
    ),
    as.data.frame(as.list(est), check.names = FALSE)
  )
}

cat(sprintf("Primary grid: %d rows (%d simulations per cell).\n", nrow(grid), n_sim))
cat(sprintf("Workers: %d\n", workers))

if (workers > 1L) {
  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Install packages 'future' and 'future.apply' to use workers > 1.")
  }
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  out <- future.apply::future_lapply(
    seq_len(nrow(grid)), run_one, future.seed = TRUE
  )
} else {
  out <- lapply(seq_len(nrow(grid)), run_one)
}

res <- dplyr::bind_rows(out)
saveRDS(res, output_file)
write.csv(res, sub("\\.rds$", ".csv", output_file), row.names = FALSE)

cat(sprintf("Saved %d rows to %s\n", nrow(res), output_file))
cat(sprintf("Replicate-level failures: %d\n", sum(res$failed %in% TRUE, na.rm = TRUE)))
