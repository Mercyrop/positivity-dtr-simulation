.PHONY: requirements smoke quick-primary summarise figures

requirements:
	Rscript requirements.R

smoke:
	Rscript scripts/01_smoke_test.R

quick-primary:
	Rscript scripts/02_run_primary.R 5 1 results/primary_simulation_test.rds

summarise:
	Rscript scripts/05_summarise_results.R

figures:
	Rscript scripts/06_make_tables_figures.R
