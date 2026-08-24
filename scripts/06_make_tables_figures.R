#!/usr/bin/env Rscript

## Recreate the core manuscript tables and figures from the public simulation
## outputs. CSV tables are intentionally plain and journal-agnostic.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

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
figures_dir <- file.path(repo_root, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(results_dir, "final_summary.csv")
if (!file.exists(summary_path)) stop("Run scripts/05_summarise_results.R first.")
sumdat <- read.csv(summary_path, stringsAsFactors = FALSE)

primary_files <- list.files(
  results_dir,
  pattern = "^primary_simulation.*\\.rds$",
  full.names = TRUE
)
primary_files <- primary_files[!grepl("_test\\.rds$", basename(primary_files))]
if (!length(primary_files)) stop("Primary simulation RDS files are required for positivity diagnostics.")
primary <- bind_rows(lapply(primary_files, readRDS))

scenario_levels <- c("benchmark", "mild", "moderate", "severe")
scenario_labels <- c(
  benchmark = "Benchmark",
  mild = "Mild",
  moderate = "Moderate",
  severe = "Severe"
)

## Table 1: regime-specific positivity diagnostics ---------------------------
table1 <- primary |>
  group_by(n, scenario) |>
  summarise(
    dtr1_ess = median(IPTW_PAR_ess_d1, na.rm = TRUE),
    dtr1_ess_pct = 100 * dtr1_ess / first(n),
    dtr1_w99 = median(IPTW_PAR_w99_d1, na.rm = TRUE),
    dtr2_ess = median(IPTW_PAR_ess_d2, na.rm = TRUE),
    dtr2_ess_pct = 100 * dtr2_ess / first(n),
    dtr2_w99 = median(IPTW_PAR_w99_d2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    Scenario = unname(scenario_labels[as.character(scenario)])
  ) |>
  arrange(n, scenario) |>
  select(
    n, Scenario,
    DTR1_ESS_percent = dtr1_ess_pct,
    DTR1_ESS_n = dtr1_ess,
    DTR1_weight_p99 = dtr1_w99,
    DTR2_ESS_percent = dtr2_ess_pct,
    DTR2_ESS_n = dtr2_ess,
    DTR2_weight_p99 = dtr2_w99
  )
write.csv(table1, file.path(results_dir, "Table1_positivity_diagnostics.csv"), row.names = FALSE)

## Primary estimator rows -----------------------------------------------------
primary_estimators <- c(
  "IPTW_PAR", "IPTW_SL",
  "AIPTW_PAR", "AIPTW_SL",
  "TMLE_PAR", "TMLE_SL",
  "GCOMP_PAR", "GCOMP_SL"
)

main <- sumdat |>
  filter(estimator %in% primary_estimators) |>
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    estimator_family = case_when(
      grepl("^IPTW_", estimator) ~ "IPTW-MSM",
      grepl("^AIPTW_", estimator) ~ "AIPTW",
      grepl("^TMLE_", estimator) ~ "TMLE",
      grepl("^GCOMP_", estimator) ~ "G-computation",
      TRUE ~ estimator
    ),
    nuisance = ifelse(grepl("_SL$", estimator), "Super Learner", "Parametric")
  )

cols <- c(
  "n", "estimator_family", "nuisance", "mean_rr", "variance",
  "rel_bias_pct", "empirical_sd", "mse", "se_to_empirical_sd",
  "mean_ci_width", "coverage_pct", "n_valid", "n_ci_valid"
)

write.csv(
  main |> filter(scenario == "benchmark") |> select(all_of(cols)),
  file.path(results_dir, "Table2_benchmark_performance.csv"),
  row.names = FALSE
)
write.csv(
  main |> filter(scenario == "severe") |> select(all_of(cols)),
  file.path(results_dir, "Table3_severe_performance.csv"),
  row.names = FALSE
)

## Sensitivity table, if available -------------------------------------------
sens_path <- file.path(results_dir, "sensitivity_summary.csv")
if (file.exists(sens_path)) {
  sens <- read.csv(sens_path, stringsAsFactors = FALSE)
  write.csv(
    sens,
    file.path(results_dir, "TableS_sensitivity_full.csv"),
    row.names = FALSE
  )
}

## Figure data ---------------------------------------------------------------
plotdat <- main |>
  mutate(
    scenario_label = factor(
      unname(scenario_labels[as.character(scenario)]),
      levels = unname(scenario_labels[scenario_levels])
    ),
    n_label = factor(
      paste0("n = ", format(n, big.mark = ",")),
      levels = c("n = 500", "n = 1,500", "n = 5,000")
    )
  )

base_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

p1 <- ggplot(
  plotdat,
  aes(
    x = scenario_label, y = rel_bias_pct,
    group = interaction(estimator_family, nuisance),
    colour = estimator_family, linetype = nuisance, shape = nuisance
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~n_label, nrow = 1) +
  labs(
    x = "Positivity scenario", y = "Relative bias (%)",
    colour = "Estimator", linetype = "Nuisance model", shape = "Nuisance model"
  ) +
  base_theme

p2 <- ggplot(
  plotdat,
  aes(
    x = scenario_label, y = mse,
    group = interaction(estimator_family, nuisance),
    colour = estimator_family, linetype = nuisance, shape = nuisance
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~n_label, nrow = 1) +
  labs(
    x = "Positivity scenario", y = "Mean squared error (log scale)",
    colour = "Estimator", linetype = "Nuisance model", shape = "Nuisance model"
  ) +
  base_theme

p3 <- ggplot(
  plotdat |> filter(is.finite(coverage_pct)),
  aes(
    x = scenario_label, y = coverage_pct,
    group = interaction(estimator_family, nuisance),
    colour = estimator_family, linetype = nuisance, shape = nuisance
  )
) +
  geom_hline(yintercept = 95, linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~n_label, nrow = 1) +
  labs(
    x = "Positivity scenario", y = "95% CI coverage (%)",
    colour = "Estimator", linetype = "Nuisance model", shape = "Nuisance model"
  ) +
  base_theme

save_plot <- function(plot, stem) {
  ggsave(file.path(figures_dir, paste0(stem, ".png")), plot,
         width = 11, height = 4.8, dpi = 600)
  ggsave(file.path(figures_dir, paste0(stem, ".pdf")), plot,
         width = 11, height = 4.8)
}

save_plot(p1, "Figure1_relative_bias")
save_plot(p2, "Figure2_mse")
save_plot(p3, "Figure3_coverage")

cat("Created Tables 1-3 and Figures 1-3.\n")
