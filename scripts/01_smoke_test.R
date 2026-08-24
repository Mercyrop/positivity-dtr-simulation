#!/usr/bin/env Rscript

## Fast end-to-end check of the public release.

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

stopifnot(abs(TRUE_RR - 1.469604947) < 1e-12)
stopifnot(all(vapply(
  c("simulate_hiv_violation_v6", "make_seed_v7", "one_sim_all_estimators_v7"),
  exists, logical(1), mode = "function"
)))

n <- 500L
scenario <- "benchmark"
sim <- 1L
seed <- make_seed_v7(n = n, scenario = scenario, sim = sim)

cat("Running one smoke-test replicate...\n")
cat(sprintf("n=%d, scenario=%s, sim=%d, seed=%d\n", n, scenario, sim, seed))

ans <- one_sim_all_estimators_v7(
  n = n,
  seed = seed,
  scenario = scenario,
  gcomp_B = 0L
)

primary <- c("IPTW_PAR_est", "AIPTW_PAR_est", "TMLE_PAR_est", "GCOMP_PAR_est")
print(ans[primary])

if (!all(is.finite(ans[primary]))) {
  stop("Smoke test failed: at least one primary parametric estimate is non-finite.")
}

cat("\nSmoke test passed.\n")
cat("Locked true RR:", TRUE_RR, "\n")
cat("R version:", R.version.string, "\n")
