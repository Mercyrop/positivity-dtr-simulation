# Simulation code for practical positivity violations in longitudinal dynamic treatment regimes

This repository contains the reproducible simulation code accompanying a methodological study of causal estimators for longitudinal dynamic treatment regimes (DTRs) under progressively worsening **practical positivity violations** in an HIV treatment-switching setting.

The public release is organised so that the scientific workflow can be run without the original project-specific HPC directory structure. The core estimator implementation is retained from the frozen analysis code; lightweight runner scripts provide a simpler interface for reproduction and manuscript review.

## Study design

The simulation follows individuals from baseline through treatment-decision visits at months 3, 6 and 9, with viral suppression assessed at month 12.

- **DTR1 (early switch):** switch when viral load is at least 50 copies/mL and remain switched thereafter.
- **DTR2 (delayed switch):** switch when viral load exceeds 1000 copies/mL and remain switched thereafter.
- **Outcome:** viral load below 50 copies/mL at month 12.
- **Primary estimand:** marginal risk ratio, `E[Y^d1] / E[Y^d2]`.
- **Locked simulation truth:** `RR = 1.469604947`, with `E[Y^d1] = 0.5716074` and `E[Y^d2] = 0.3889531`.
- **Sample sizes:** 500, 1500 and 5000.
- **Positivity scenarios:** benchmark, mild, moderate and severe.
- **Monte Carlo replicates:** 1000 per sample-size/scenario cell in the manuscript analysis.

Positivity deterioration is generated through a confirmation barrier in the observational treatment process. As the scenario becomes more severe, switching after an unconfirmed detectable-but-not-failing viral-load result becomes increasingly uncommon, while switching after virologic failure remains supported. The biological data-generating process and the two intervention definitions are unchanged across scenarios.

## Estimators

The primary comparison includes:

- inverse-probability-of-treatment weighted marginal structural models (IPTW-MSM);
- augmented IPTW (AIPTW);
- targeted maximum likelihood estimation (TMLE); and
- parametric g-computation.

Each primary estimator is implemented using parametric nuisance models and, where applicable, Super Learner. The primary treatment probabilities are bounded at 0.01 and 0.99; deterministic post-switch treatment probabilities remain exactly one. No additional percentile truncation of final weights is used in the primary IPTW analysis.

Sensitivity analyses evaluate 1st/99th percentile weight truncation, stabilized weights, stabilized plus truncated weights, treatment-model misspecification, outcome-model misspecification, and joint nuisance-model misspecification.

## Repository structure

```text
positivity_simulation_code/
├── README.md
├── CODE_VERSION.md
├── requirements.R
├── CITATION.cff.template
├── .gitignore
├── R/
│   ├── 01_data_generating_mechanism.R
│   ├── 02_estimators.R
│   └── 03_misspecification.R
├── scripts/
│   ├── 01_smoke_test.R
│   ├── 02_run_primary.R
│   ├── 03_run_gcomp_bootstrap.R
│   ├── 04_run_sensitivity.R
│   ├── 05_summarise_results.R
│   └── 06_make_tables_figures.R
├── hpc/
│   ├── README.md
│   ├── run_primary_cell.pbs.example
│   └── run_gcomp_cell.pbs.example
├── results/
└── figures/
```

## Software requirements

The manuscript analyses were run in R 4.5.1. Required packages are listed in `requirements.R`:

- `dplyr`
- `tidyr`
- `SuperLearner`
- `survey`
- `gam`
- `future`
- `future.apply`
- `ggplot2`

Run the following to check your local installation:

```bash
Rscript requirements.R
```

The script reports missing packages and prints `sessionInfo()`; it does not install packages automatically.

## Quick reproducibility check

From the repository root:

```bash
Rscript scripts/01_smoke_test.R
```

In an interactive R/RStudio session opened at the repository root, the equivalent command is:

```r
source("scripts/01_smoke_test.R")
```

This runs a single `n = 500` benchmark replicate and checks that the four primary parametric estimators return finite risk-ratio estimates. The public release has also passed this smoke test under R 4.4.1 on Windows 11; the manuscript production environment used R 4.5.1.

## Reproducing the manuscript simulations

### 1. Primary simulation

The defaults reproduce the manuscript grid of 1000 simulations per cell:

```bash
Rscript scripts/02_run_primary.R 1000 1
```

The second argument is the number of local workers. For example, to use 12 workers:

```bash
Rscript scripts/02_run_primary.R 1000 12
```

Outputs are written to `results/primary_simulation.rds` and `results/primary_simulation.csv`.

For a quick test before the full run:

```bash
Rscript scripts/02_run_primary.R 5 1 results/primary_simulation_test.rds
```

### 2. G-computation bootstrap confidence intervals

G-computation confidence intervals were obtained separately using **200 nonparametric bootstrap samples per Monte Carlo replicate**. Because this is computationally expensive, running one sample-size/scenario cell at a time is recommended.

For example:

```bash
Rscript scripts/03_run_gcomp_bootstrap.R 500 benchmark 1000 200 12
Rscript scripts/03_run_gcomp_bootstrap.R 500 mild      1000 200 12
Rscript scripts/03_run_gcomp_bootstrap.R 500 moderate  1000 200 12
Rscript scripts/03_run_gcomp_bootstrap.R 500 severe    1000 200 12
```

Repeat for `n = 1500` and `n = 5000`. The summarisation script automatically combines files named `gcomp_bootstrap_*.rds`.

A small plumbing test can be run with:

```bash
Rscript scripts/03_run_gcomp_bootstrap.R 500 benchmark 2 5 1 results/gcomp_bootstrap_test.rds
```

### 3. Sensitivity analyses

The manuscript sensitivity analyses use `n = 500` and `n = 5000` with 1000 simulations per cell:

```bash
Rscript scripts/04_run_sensitivity.R all all 1000 12
```

Individual cells can also be run separately, for example:

```bash
Rscript scripts/04_run_sensitivity.R 500 severe 1000 12
```

### 4. Summarise estimator performance

After the primary grid and, where required, the G-computation bootstrap have completed:

```bash
Rscript scripts/05_summarise_results.R
```

The script checks that G-computation point estimates from the bootstrap run agree with the corresponding primary-run point estimates before combining the results. It then writes `results/final_summary.csv` and `results/final_summary.rds`.

If sensitivity files are present, `results/sensitivity_summary.csv` is also created.

### 5. Recreate core tables and figures

```bash
Rscript scripts/06_make_tables_figures.R
```

This creates plain CSV versions of the core tables in `results/` and figures in both PNG and PDF format in `figures/`.

## Parallel and HPC execution

The public runner scripts are scheduler-independent. They can be run serially, with local `future` workers, or wrapped in a cluster scheduler. Example PBS wrappers are provided in `hpc/`. The manuscript production run used a chunked HPC workflow, but those project-specific orchestration details are not required to understand or reproduce the statistical analysis.

## Reproducibility and random-number generation

Monte Carlo seeds are deterministic functions of sample size, positivity scenario and replicate number through `make_seed_v7()`. Therefore, a given `(n, scenario, sim)` combination is reproducible independently of submission order or chunk boundaries. The primary grid and the separate G-computation bootstrap use the same outer-replicate seed, which permits direct point-estimate consistency checks before bootstrap interval results are attached.

## Code provenance

The manuscript results were finalised under the internal analysis label `v7_frozen_2026-08-04_v7`, with a locked true marginal risk ratio of `1.469604947`. The archived estimator snapshot had MD5 checksum `30476166f4b698161bc2304cef98e2c0`. A byte-identical copy is retained under `archive/estimators_rr_v7_FROZEN.R`; `R/02_estimators.R` contains the same executable statements with development-oriented comments cleaned for public release.

`R/01_data_generating_mechanism.R` is a public-facing extraction of the active v6 data generator and positivity-diagnostic functions. Development-only and deprecated sections were omitted for clarity, so its checksum is intentionally different from the internal project archive. See `CODE_VERSION.md` for details.

## Results files

Large replicate-level simulation outputs are not tracked by default. The `.gitignore` excludes generated `.rds`, `.csv`, `.png` and `.pdf` files under `results/` and `figures/`, while preserving the directory structure. For manuscript submission, summary tables may be committed separately if required by the journal.

## Citation and licence

A `CITATION.cff.template` is included so the manuscript authors can add the final article title, author list, DOI and repository URL once these are fixed. No software licence has been imposed in this release because the appropriate licence should be selected by the authors/institution before making the repository public.
