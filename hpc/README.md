# Optional HPC execution

The scripts in `../scripts/` are scheduler-independent. For large production runs, the recommended approach is to submit one sample-size/scenario cell per job rather than running every cell in a single process.

The example PBS files in this directory assume an environment-modules system and R 4.5.1. Edit queue/account/resource directives for the target cluster before use.

A typical workflow is:

```bash
qsub -v N_VAL=500,SCENARIO=benchmark,N_SIM=1000,NCPUS=24 hpc/run_primary_cell.pbs.example
qsub -v N_VAL=500,SCENARIO=benchmark,N_SIM=1000,B_VAL=200,NCPUS=24 hpc/run_gcomp_cell.pbs.example
```

These examples write cell-specific RDS files into `results/`. Once all cells are present, run:

```bash
Rscript scripts/05_summarise_results.R
Rscript scripts/06_make_tables_figures.R
```

The original study used additional chunking to accommodate wall-time limits. Chunk boundaries do not affect the simulated data because seeds are deterministic functions of `(n, scenario, sim)`.
