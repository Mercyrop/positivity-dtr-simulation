# Reference outputs from the final manuscript analysis

These values provide convenient checks that a full reproduction is on the expected scale. Monte Carlo summaries will not match every displayed decimal unless the same full 1000-replicate analysis is reproduced.

## Locked causal truth

- `E[Y^d1] = 0.5716074`
- `E[Y^d2] = 0.3889531`
- `RR = 1.469604947`

## Positivity diagnostics

For `n = 5000`, the final analysis reported the following median regime-specific diagnostics:

| Scenario | DTR1 relative ESS | DTR1 99th-percentile weight | DTR2 relative ESS | DTR2 99th-percentile weight |
|---|---:|---:|---:|---:|
| Benchmark | 45.2% | 3.6 | 36.8% | 7.8 |
| Mild | 32.9% | 5.5 | 44.4% | 6.2 |
| Moderate | 24.5% | 7.8 | 53.7% | 4.4 |
| Severe | 8.4% | 23.6 | 68.9% | 2.8 |

The intended pattern is asymmetric: support for the early-switch regime deteriorates as the confirmation barrier strengthens, whereas support for the delayed-switch regime improves.

## Benchmark estimator check at n = 500

Selected parametric results from the final 1000-replicate analysis were approximately:

| Estimator | Mean RR | Relative bias | Empirical SD | MSE | Coverage |
|---|---:|---:|---:|---:|---:|
| IPTW-MSM | 1.4719 | 0.16% | 0.0939 | 0.00881 | 99.5% |
| AIPTW | 1.4760 | 0.43% | 0.0778 | 0.00609 | 94.7% |
| TMLE | approximately 1.476 | approximately 0.43% | approximately 0.0776 | 0.00606 | 94.6% |
| G-computation | approximately 1.471 | 0.09% | 0.0521 | 0.00272 | 94.4% |

These are reference checks, not hard-coded analysis inputs.

## Public-release smoke-test validation

The public repository smoke test was also run successfully under **R 4.4.1 on Windows 11**. For the deterministic test cell
`n = 500`, `scenario = "benchmark"`, `sim = 1`, `seed = 101010001`, the primary parametric estimates were:

| Estimator | RR estimate |
|---|---:|
| IPTW-MSM | 1.409674 |
| AIPTW | 1.358813 |
| TMLE | 1.360698 |
| Parametric g-computation | 1.368523 |

These are single-replicate values, not performance summaries. Their purpose is to provide a practical reference that the code executes end-to-end and produces finite estimates. Small platform- or package-version numerical differences are not, by themselves, evidence of a reproduction failure.
