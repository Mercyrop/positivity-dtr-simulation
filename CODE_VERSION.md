# Code version and provenance

## Manuscript analysis snapshot

The final manuscript tables and results archive identify the analysis as:

- **analysis label:** `v7_frozen_2026-08-04_v7`
- **true E[Y^d1]:** `0.5716074`
- **true E[Y^d2]:** `0.3889531`
- **true RR:** `1.469604947`
- **archived DGP MD5 recorded with the final results:** `0e1081a5ec0d9a2c5f8edd1bc33fb64b`
- **archived estimator MD5:** `30476166f4b698161bc2304cef98e2c0`

## Public release

This repository separates the scientific core from the original HPC orchestration.

A byte-identical copy of the frozen estimator snapshot is retained as `archive/estimators_rr_v7_FROZEN.R`, with MD5 `30476166f4b698161bc2304cef98e2c0`. The public-facing `R/02_estimators.R` preserves the same executable statements but cleans development-oriented comments, so its checksum differs.

`R/01_data_generating_mechanism.R` is a cleaned public extraction of the active `simulate_hiv_violation_v6()` generator, the switch-barrier diagnostics, and the locked truth constants. Historical/deprecated development functions and comments that are not required to reproduce the analysis were removed. Consequently, the public DGP file is not expected to have the same byte-level checksum as the internal archived DGP. The calibrated scenario parameters and active simulation equations were preserved from the frozen v6 code lineage.

The public runner scripts are newly written orchestration wrappers. They do not redefine the DGP or estimator mathematics; their purpose is to provide a transparent, scheduler-independent path from simulation to summary tables and figures.

## Reproducibility principle

For scientific reproduction, the relevant invariants are:

1. the locked estimand and truth values above;
2. the four final switch-barrier scenarios;
3. the deterministic `make_seed_v7()` mapping;
4. the final estimator functions in `R/02_estimators.R`;
5. the nuisance-model specifications and treatment-probability bounds;
6. 1000 Monte Carlo replicates per primary cell; and
7. 200 bootstrap samples per outer replicate for g-computation interval estimation.

## Public-release runtime validation

The released code passed the end-to-end smoke test under R 4.4.1 on Windows 11 after package-availability checks succeeded. The manuscript production analyses used R 4.5.1. Runner scripts include repository-root detection for both `Rscript` execution and interactive `source()` use.
