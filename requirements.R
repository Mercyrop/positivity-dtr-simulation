#!/usr/bin/env Rscript

required_packages <- c(
  "dplyr",
  "tidyr",
  "SuperLearner",
  "survey",
  "gam",
  "future",
  "future.apply",
  "ggplot2"
)

installed <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)

cat("Required R packages:\n")
print(data.frame(package = required_packages, installed = installed), row.names = FALSE)

if (!all(installed)) {
  cat("\nMissing packages:\n  ", paste(required_packages[!installed], collapse = ", "), "\n", sep = "")
  cat("Install them using your institution's preferred R package workflow.\n")
} else {
  cat("\nAll required packages are available.\n")
}

cat("\nSession information:\n")
print(sessionInfo())
