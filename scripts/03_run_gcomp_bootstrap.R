#!/usr/bin/env Rscript

## G-computation bootstrap for confidence intervals.
##
## Usage:
##   Rscript scripts/03_run_gcomp_bootstrap.R [n|all] [scenario|all] [n_sim] [B] [workers] [output_file]
##
## Manuscript settings: n_sim=1000 and B=200.
## Running one n x scenario cell at a time is recommended on an HPC system.

args <- commandArgs(trailingOnly = TRUE)
n_arg <- if (length(args) >= 1L) args[1] else "all"
scenario_arg <- if (length(args) >= 2L) args[2] else "all"
n_sim <- if (length(args) >= 3L) as.integer(args[3]) else 1000L
B <- if (length(args) >= 4L) as.integer(args[4]) else 200L
workers <- if (length(args) >= 5L) as.integer(args[5]) else 1L

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

source(file.path(repo_root, "R", "01_data_generating_mechanism.R"))
source(file.path(repo_root, "R", "02_estimators.R"))

all_n <- c(500L, 1500L, 5000L)
all_scenarios <- c("benchmark", "mild", "moderate", "severe")

n_values <- if (identical(tolower(n_arg), "all")) all_n else as.integer(n_arg)
scenarios <- if (identical(tolower(scenario_arg), "all")) all_scenarios else scenario_arg

if (any(!n_values %in% all_n)) stop("n must be one of 500, 1500, 5000, or 'all'.")
if (any(!scenarios %in% all_scenarios)) stop("Unknown scenario.")

suffix <- paste0(
  if (length(n_values) == 1L) paste0("n", n_values) else "alln",
  "_",
  if (length(scenarios) == 1L) scenarios else "allscenarios"
)

output_file <- if (length(args) >= 6L) {
  args[6]
} else {
  file.path(repo_root, "results", paste0("gcomp_bootstrap_", suffix, ".rds"))
}
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

grid <- expand.grid(
  n = n_values,
  scenario = scenarios,
  sim = seq_len(n_sim),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid <- grid[order(grid$n, match(grid$scenario, all_scenarios), grid$sim), ]
grid$seed <- as.integer(mapply(
  make_seed_v7,
  n = grid$n,
  scenario = grid$scenario,
  sim = grid$sim
))

prefix_result <- function(x, prefix) {
  names(x) <- paste0(prefix, "_", names(x))
  x
}

run_one <- function(i) {
  row <- grid[i, ]
  tryCatch({
    dat <- simulate_hiv_violation_v6(
      n = row$n,
      seed = row$seed,
      scenario = row$scenario
    )
    wide <- prepare_wide_estimation_v7(dat)

    par <- bootstrap_gcomp_oracle_v7(
      wide,
      use_sl = FALSE,
      B = B,
      seed = row$seed,
      n_quad = DEFAULT_N_QUAD
    )
    sl <- bootstrap_gcomp_oracle_v7(
      wide,
      use_sl = TRUE,
      SL.library = stable_sl_library_v7(),
      B = B,
      seed = row$seed,
      n_quad = DEFAULT_N_QUAD
    )

    cbind(
      data.frame(
        n = row$n, scenario = row$scenario, sim = row$sim, seed = row$seed,
        B = B, failed = FALSE, error_message = NA_character_,
        stringsAsFactors = FALSE
      ),
      as.data.frame(as.list(c(
        prefix_result(par, "GCOMP_PAR"),
        prefix_result(sl, "GCOMP_SL")
      )), check.names = FALSE)
    )
  }, error = function(e) {
    data.frame(
      n = row$n, scenario = row$scenario, sim = row$sim, seed = row$seed,
      B = B, failed = TRUE, error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

cat(sprintf("G-computation bootstrap: %d outer replicates; B=%d; workers=%d\n",
            nrow(grid), B, workers))

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
