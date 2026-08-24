# Results directory

Simulation outputs are generated here by the scripts in `../scripts/`.

Large replicate-level files are intentionally excluded from Git by default. The main generated summaries are:

- `final_summary.csv` — primary estimator performance across sample sizes and positivity scenarios;
- `sensitivity_summary.csv` — nuisance-model and IPTW weight-modification sensitivity results;
- `Table1_positivity_diagnostics.csv`;
- `Table2_benchmark_performance.csv`; and
- `Table3_severe_performance.csv`.

The manuscript production analysis used 1000 Monte Carlo replicates per primary cell and 200 bootstrap resamples per outer replicate for g-computation confidence intervals.
