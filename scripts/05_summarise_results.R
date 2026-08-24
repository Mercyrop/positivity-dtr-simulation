#!/usr/bin/env Rscript

## Combine simulation outputs and calculate ADEMP performance measures.
##
## Run this after the primary grid. If G-computation bootstrap files are present,
## their SE/CI/coverage results replace the point-only G-computation rows from
## the primary grid. If sensitivity files are present, a separate sensitivity
## summary is also written.

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
results_dir <- file.path(repo_root, "results")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_root, "R", "01_data_generating_mechanism.R"))
source(file.path(repo_root, "R", "02_estimators.R"))

read_many <- function(pattern) {
  files <- list.files(results_dir, pattern = pattern, full.names = TRUE)
  files <- files[!grepl("_test\\.rds$", basename(files))]
  if (!length(files)) return(NULL)
  dplyr::bind_rows(lapply(files, readRDS))
}

sanitize_numeric <- function(x, protect = c("n", "sim", "seed", "B")) {
  for (nm in setdiff(names(x), protect)) {
    if (is.numeric(x[[nm]])) x[[nm]][!is.finite(x[[nm]])] <- NA_real_
  }
  x
}

check_unique <- function(x, keys, label) {
  dup <- x |>
    dplyr::count(dplyr::across(dplyr::all_of(keys)), name = "copies") |>
    dplyr::filter(copies > 1L)
  if (nrow(dup)) {
    stop(label, " contains duplicated keys. Resolve duplicates before summarising.")
  }
}

primary <- read_many("^primary_simulation.*\\.rds$")
if (is.null(primary) || !nrow(primary)) {
  stop("No results/primary_simulation*.rds files were found.")
}
primary <- sanitize_numeric(primary)
check_unique(primary, c("n", "scenario", "sim"), "Primary results")

cat(sprintf("Primary results: %d rows.\n", nrow(primary)))

main_summary <- summarise_results_grid_v7c(
  primary,
  true_rr = TRUE_RR,
  gcomp_ci_expected = FALSE
)

gcomp <- read_many("^gcomp_bootstrap_.*\\.rds$")

if (!is.null(gcomp) && nrow(gcomp)) {
  gcomp <- sanitize_numeric(gcomp)
  check_unique(gcomp, c("n", "scenario", "sim"), "G-computation bootstrap results")

  ## Confirm the separately bootstrapped point estimates reproduce the same
  ## point estimates as the primary run for every matched outer replicate.
  check <- primary |>
    dplyr::select(
      n, scenario, sim,
      main_par = GCOMP_PAR_est,
      main_sl = GCOMP_SL_est
    ) |>
    dplyr::inner_join(
      gcomp |>
        dplyr::select(
          n, scenario, sim,
          boot_par = GCOMP_PAR_est,
          boot_sl = GCOMP_SL_est
        ),
      by = c("n", "scenario", "sim")
    )

  max_diff <- function(a, b) {
    ok <- is.finite(a) & is.finite(b)
    if (!any(ok)) return(NA_real_)
    max(abs(a[ok] - b[ok]))
  }

  d_par <- max_diff(check$main_par, check$boot_par)
  d_sl <- max_diff(check$main_sl, check$boot_sl)
  cat(sprintf("GCOMP point-estimate max |difference|: PAR=%g, SL=%g\n", d_par, d_sl))

  if (is.finite(d_par) && d_par > 1e-10) stop("Parametric GCOMP point estimates do not match.")
  if (is.finite(d_sl) && d_sl > 1e-8) stop("SL GCOMP point estimates do not match.")

  gcomp_summary <- summarise_results_grid_v7c(
    gcomp,
    true_rr = TRUE_RR,
    gcomp_ci_expected = TRUE
  ) |>
    dplyr::filter(grepl("^GCOMP_", as.character(estimator)))

  final_summary <- dplyr::bind_rows(
    main_summary |>
      dplyr::filter(!grepl("^GCOMP_", as.character(estimator))),
    gcomp_summary
  ) |>
    dplyr::mutate(
      scenario = factor(scenario, levels = c("benchmark", "mild", "moderate", "severe"))
    ) |>
    dplyr::arrange(n, scenario, estimator)
} else {
  cat("No G-computation bootstrap files found; GCOMP interval metrics remain NA.\n")
  final_summary <- main_summary
}

saveRDS(final_summary, file.path(results_dir, "final_summary.rds"))
write.csv(final_summary, file.path(results_dir, "final_summary.csv"), row.names = FALSE)

## Sensitivity analyses -------------------------------------------------------
sensitivity <- read_many("^sensitivity_.*\\.rds$")
if (!is.null(sensitivity) && nrow(sensitivity)) {
  sensitivity <- sanitize_numeric(sensitivity)
  check_unique(sensitivity, c("n", "scenario", "sim", "setting"), "Sensitivity results")

  settings <- c("Correct", "Pmis", "Qmis", "Both")
  summaries <- lapply(settings, function(st) {
    z <- sensitivity[sensitivity$setting == st, , drop = FALSE]
    if (!nrow(z)) return(NULL)
    summarise_results_grid_v7c(z, true_rr = TRUE_RR, gcomp_ci_expected = FALSE) |>
      dplyr::mutate(setting = st)
  })
  sensitivity_summary <- dplyr::bind_rows(summaries)

  keep_pair <- with(sensitivity_summary,
    (grepl("^IPTW_", as.character(estimator)) & setting %in% c("Correct", "Pmis")) |
    (grepl("^AIPTW_", as.character(estimator)) & setting %in% settings) |
    (grepl("^TMLE_", as.character(estimator)) & setting %in% settings) |
    (grepl("^GCOMP_", as.character(estimator)) & setting %in% c("Correct", "Qmis"))
  )
  sensitivity_summary <- sensitivity_summary[keep_pair, , drop = FALSE]

  saveRDS(sensitivity_summary, file.path(results_dir, "sensitivity_summary.rds"))
  write.csv(
    sensitivity_summary,
    file.path(results_dir, "sensitivity_summary.csv"),
    row.names = FALSE
  )
  cat(sprintf("Sensitivity summary: %d rows.\n", nrow(sensitivity_summary)))
}

cat(sprintf("Final summary written to %s\n", file.path(results_dir, "final_summary.csv")))
