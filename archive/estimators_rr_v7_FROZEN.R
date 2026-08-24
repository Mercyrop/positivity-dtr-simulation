## ============================================================================
## ESTIMATORS v7e -- v7 PRIMARY + SOFTENED Pmis + IPTW RAW/TRUNC01_99/
##                    STAB/STAB_TRUNC01_99, modified-Poisson stabilized MSM
##
## Pmis is the INTERMEDIATE ("softened") specification: det_only + fail +
## persist + age_c + male (drops deep persistence and the conf/unconf split,
## keeps ordinary persist). This is NOT the earlier "strong" Pmis (which
## additionally dropped persist and produced 84.5%/57.9% coverage even at
## benchmark) -- that version is retired. Confirm this header matches the
## actual ps_xvars_v7() code below before trusting it; version-label/code
## mismatches have caused real confusion on HPC before.
##
## Qmis (transition_x_v7) is the ACCEPTED moderated specification (drops
## only A:slope_B) -- matches the validated ~-4.8% to -7.5% GCOMP bias
## pattern. Do not change this while updating IPTW-related code.
##
## Primary estimand:
##   RR = E[Y^d1] / E[Y^d2]
##
## DTR1: switch when VL >= 50 copies/mL, then remain switched.
## DTR2: switch when VL > 1000 copies/mL, then remain switched.
## Outcome: VL < 50 copies/mL at month 12.
##
## DESIGN PRINCIPLES
##   1. The treatment model is correctly specified for the v6 switch-barrier DGP.
##   2. The outcome model is the correctly specified Gaussian viral-load
##      transition model used by the DGP. Its coefficients and residual SD are
##      re-estimated in every simulated sample.
##   3. G-computation, AIPTW and TMLE use the same fitted outcome-transition
##      engine and the same dynamic-regime definitions.
##   4. IPTW, AIPTW and TMLE use the same estimated treatment probabilities.
##   5. Stochastic treatment probabilities are bounded at 0.01 and 0.99.
##      Deterministic post-switch probabilities remain exactly 1.
##   6. No additional percentile truncation of the final weights is applied in
##      the primary analysis. The 0.01/0.99 bound is applied to node-specific
##      stochastic treatment probabilities.
##
## REQUIRED DGP
##   Source simulate_hiv_all_v6.R before using this file. The active generator
##   must be simulate_hiv_violation_v6().
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(SuperLearner)
  library(survey)
})

## *** DO NOT TRUST THESE UNTIL compute_true_values_v7() HAS BEEN RUN ***
## These stored constants are internally inconsistent: 0.5717/0.3888 = 1.4704,
## not 1.4678. The discrepancy (0.18%) is too large to be rounding noise --
## most likely these were computed as a mean of five seed-specific RR ratios
## rather than the ratio of pooled marginal risks, which are not the same
## quantity in general. Every bias/coverage number computed against TRUE_RR
## anywhere in this project inherits this inconsistency until it's fixed.
## Run compute_true_values_v7() below and replace these three lines with its
## full-precision output before trusting any bias interpretation.
TRUE_RR    <- 1.469604947
TRUE_EY_D1 <- 0.5716074
TRUE_EY_D2 <- 0.3889531

## ----------------------------------------------------------------------------
## Recomputes the true risks and RR as the ratio of pooled marginal risks
## (not a mean of per-seed ratios), using common random numbers across the
## two forced regimes so the Monte Carlo SE of log(RR) can be computed
## correctly via the paired covariance, rather than assuming independence.
## ----------------------------------------------------------------------------
compute_true_values_v7 <- function(n_chunk = 200000L, n_chunks = 100L, seed_start = 900000L) {
  total_n <- 0L
  sum_y1  <- 0
  sum_y2  <- 0
  sum_y12 <- 0

  for (j in seq_len(n_chunks)) {
    current_seed <- seed_start + j - 1L

    ## Same seed for both regimes -- common random numbers, so the paired
    ## covariance below is a real, estimable quantity, not assumed away.
    d1 <- simulate_hiv_violation_v6(n = n_chunk, seed = current_seed, scenario = "benchmark", regime = "d1")
    d2 <- simulate_hiv_violation_v6(n = n_chunk, seed = current_seed, scenario = "benchmark", regime = "d2")

    y1 <- as.integer(d1$viral_load[d1$month == 12] < 50)
    y2 <- as.integer(d2$viral_load[d2$month == 12] < 50)
    stopifnot(length(y1) == length(y2))

    total_n <- total_n + length(y1)
    sum_y1  <- sum_y1 + sum(y1)
    sum_y2  <- sum_y2 + sum(y2)
    sum_y12 <- sum_y12 + sum(y1 * y2)
  }

  EY_d1 <- sum_y1 / total_n
  EY_d2 <- sum_y2 / total_n
  rr    <- EY_d1 / EY_d2
  rd    <- EY_d1 - EY_d2

  cov_y1_y2 <- sum_y12 / total_n - EY_d1 * EY_d2
  var_log_rr <-
    EY_d1 * (1 - EY_d1) / (total_n * EY_d1^2) +
    EY_d2 * (1 - EY_d2) / (total_n * EY_d2^2) -
    2 * cov_y1_y2 / (total_n * EY_d1 * EY_d2)

  se_log_rr <- sqrt(max(var_log_rr, 0))
  se_rr     <- rr * se_log_rr

  data.frame(
    n_per_regime = total_n, EY_d1 = EY_d1, EY_d2 = EY_d2, RD = rd, RR = rr,
    MCSE_RR = se_rr,
    RR_LCL_MC = exp(log(rr) - 1.96 * se_log_rr),
    RR_UCL_MC = exp(log(rr) + 1.96 * se_log_rr)
  )
}

## ----------------------------------------------------------------------------
## Deterministic seed as a pure function of (n, scenario, sim) -- NOT of job
## submission order or chunk boundaries. Results are reproducible regardless
## of how jobs are split, resubmitted, or reordered, and this lets the main
## grid and the GCOMP bootstrap use the SAME underlying simulated data for a
## given outer replicate, rather than two independently-seeded sets.
## ----------------------------------------------------------------------------
SCENARIO_LEVELS_V7 <- c("benchmark", "mild", "moderate", "severe")

make_seed_v7 <- function(n, scenario, sim) {
  n_id <- match(n, c(500, 1500, 5000))
  scenario_id <- match(scenario, SCENARIO_LEVELS_V7)
  if (is.na(n_id) || is.na(scenario_id)) stop("make_seed_v7: unrecognized n or scenario")
  as.integer(100000000L + n_id * 1000000L + scenario_id * 10000L + sim)
}

DEFAULT_G_BOUNDS <- c(0.01, 0.99)
DEFAULT_N_QUAD   <- 15L
PROB_EPS         <- 1e-8

## ============================================================================
## 1. SMALL NUMERICAL HELPERS
## ============================================================================

clamp_prob <- function(x, eps = PROB_EPS) {
  pmin(pmax(x, eps), 1 - eps)
}

expit <- function(x) {
  1 / (1 + exp(-x))
}

logit <- function(p) {
  qlogis(clamp_prob(p))
}

apply_logit_shift <- function(p, epsilon = 0) {
  expit(logit(p) + epsilon)
}

safe_ratio_result <- function() {
  c(
    est = NA_real_, se = NA_real_, lcl = NA_real_, ucl = NA_real_,
    EY_d1 = NA_real_, EY_d2 = NA_real_,
    ess_d1 = NA_real_, ess_d2 = NA_real_,
    w99_d1 = NA_real_, w99_d2 = NA_real_,
    se_log_rr = NA_real_
  )
}

## Standard-normal Gauss-Hermite quadrature using the Golub-Welsch algorithm.
## For Z ~ N(0,1), E[f(Z)] is approximated by sum_j weight_j f(node_j).
normal_quadrature <- function(n_quad = DEFAULT_N_QUAD) {
  n_quad <- as.integer(n_quad)
  if (n_quad < 3L) stop("n_quad must be at least 3")

  J <- matrix(0, nrow = n_quad, ncol = n_quad)
  off <- sqrt(seq_len(n_quad - 1L))
  J[cbind(seq_len(n_quad - 1L), 2:n_quad)] <- off
  J[cbind(2:n_quad, seq_len(n_quad - 1L))] <- off

  ee <- eigen(J, symmetric = TRUE)
  ord <- order(ee$values)

  list(
    nodes = ee$values[ord],
    weights = ee$vectors[1L, ord]^2
  )
}

## ============================================================================
## 2. PREPARE ONE-ROW-PER-SUBJECT DATA
## ============================================================================

prepare_wide_estimation_v7 <- function(
    dat_obs,
    dtr1_threshold = 50,
    dtr2_threshold = 1000
) {
  dat_obs %>%
    filter(month %in% c(3, 6, 9, 12)) %>%
    mutate(kk = match(month, c(3, 6, 9, 12))) %>%
    select(
      id, age, female, male, age_c, traj_type,
      slope_A, slope_B, baseline_log10_vl,
      kk, A_next, det, fail, det_only,
      persist_det_only, deep_persist_det_only,
      confirmed_rebound, unconfirmed_det_only,
      viral_load
    ) %>%
    pivot_wider(
      id_cols = c(
        id, age, female, male, age_c, traj_type,
        slope_A, slope_B, baseline_log10_vl
      ),
      names_from = kk,
      values_from = c(
        A_next, det, fail, det_only,
        persist_det_only, deep_persist_det_only,
        confirmed_rebound, unconfirmed_det_only,
        viral_load
      ),
      names_sep = "_k"
    ) %>%
    rename(
      A1 = A_next_k1,
      A2 = A_next_k2,
      A3 = A_next_k3,

      det1 = det_k1,
      det2 = det_k2,
      det3 = det_k3,

      fail1 = fail_k1,
      fail2 = fail_k2,
      fail3 = fail_k3,

      det_only1 = det_only_k1,
      det_only2 = det_only_k2,
      det_only3 = det_only_k3,

      persist1 = persist_det_only_k1,
      persist2 = persist_det_only_k2,
      persist3 = persist_det_only_k3,

      deep1 = deep_persist_det_only_k1,
      deep2 = deep_persist_det_only_k2,
      deep3 = deep_persist_det_only_k3,

      conf1 = confirmed_rebound_k1,
      conf2 = confirmed_rebound_k2,
      conf3 = confirmed_rebound_k3,

      unconf1 = unconfirmed_det_only_k1,
      unconf2 = unconfirmed_det_only_k2,
      unconf3 = unconfirmed_det_only_k3,

      vl1 = viral_load_k1,
      vl2 = viral_load_k2,
      vl3 = viral_load_k3,
      vl4 = viral_load_k4
    ) %>%
    mutate(
      analysis_id = row_number(),
      Y4 = as.integer(vl4 < dtr1_threshold),

      Astar1_d1 = as.integer(vl1 >= dtr1_threshold),
      Astar2_d1 = pmax(Astar1_d1, as.integer(vl2 >= dtr1_threshold)),
      Astar3_d1 = pmax(Astar2_d1, as.integer(vl3 >= dtr1_threshold)),

      Astar1_d2 = as.integer(vl1 > dtr2_threshold),
      Astar2_d2 = pmax(Astar1_d2, as.integer(vl2 > dtr2_threshold)),
      Astar3_d2 = pmax(Astar2_d2, as.integer(vl3 > dtr2_threshold)),

      logvl1 = log10(pmax(vl1, 1)),
      logvl2 = log10(pmax(vl2, 1)),
      logvl3 = log10(pmax(vl3, 1)),
      logvl4 = log10(pmax(vl4, 1))
    )
}

## Backward-compatible alias for existing runners.
prepare_wide_estimation_v5 <- prepare_wide_estimation_v7

## ============================================================================
## 3. CORRECTLY SPECIFIED TREATMENT MODEL
## ============================================================================

## det_only is intentionally omitted because, by construction,
## det_only = confirmed_rebound + unconfirmed_det_only.
## Using conf and unconf separately gives an exact nonredundant reparameterisation
## of the v6 treatment law.
ps_xvars_v7 <- function(df, k, misspecified = FALSE) {
  if (!misspecified) {
    data.frame(
      fail = df[[paste0("fail", k)]],
      conf = df[[paste0("conf", k)]],
      unconf = df[[paste0("unconf", k)]],
      persist = df[[paste0("persist", k)]],
      deep = df[[paste0("deep", k)]],
      age_c = df$age_c,
      male = df$male
    )
  } else {
    ## Softened, clinically interpretable propensity-model misspecification.
    ## The analyst observes the current detectable/failure categories and
    ## ordinary one-visit persistence, but does not distinguish confirmed
    ## from unconfirmed rebound and does not retain deeper persistence history.
    ## This is deliberately intermediate between the original weak Pmis and
    ## the v7b strong Pmis. The DGP itself remains unchanged.
    data.frame(
      det_only = df[[paste0("det_only", k)]],
      fail = df[[paste0("fail", k)]],
      persist = df[[paste0("persist", k)]],
      age_c = df$age_c,
      male = df$male
    )
  }
}

fit_binary_nuisance_v7 <- function(
    Y,
    X,
    use_sl = FALSE,
    SL.library = c("SL.glm", "SL.gam")
) {
  X <- as.data.frame(X)
  keep <- vapply(
    X,
    function(x) length(unique(x[is.finite(x)])) > 1L,
    logical(1)
  )

  if (!any(keep)) {
    X <- data.frame(intercept_proxy = rep(0, length(Y)))
    keep_names <- "intercept_proxy"
  } else {
    X <- X[, keep, drop = FALSE]
    keep_names <- names(X)
  }

  if (!use_sl) {
    fit <- glm(Y ~ ., data = cbind(Y = Y, X), family = binomial())
  } else {
    fit <- SuperLearner(
      Y = Y,
      X = X,
      family = binomial(),
      SL.library = SL.library
    )
  }

  list(
    fit = fit,
    use_sl = use_sl,
    keep_names = keep_names
  )
}

predict_binary_nuisance_v7 <- function(object, Xnew) {
  Xnew <- as.data.frame(Xnew)

  if (identical(object$keep_names, "intercept_proxy")) {
    Xnew <- data.frame(intercept_proxy = rep(0, nrow(Xnew)))
  } else {
    Xnew <- Xnew[, object$keep_names, drop = FALSE]
  }

  if (!object$use_sl) {
    as.numeric(predict(object$fit, newdata = Xnew, type = "response"))
  } else {
    as.numeric(
      predict(object$fit, newdata = Xnew, onlySL = TRUE)$pred[, 1]
    )
  }
}

fit_ps_v7 <- function(
    wide_df,
    use_sl = FALSE,
    misspecified = FALSE,
    SL.library = c("SL.glm", "SL.gam"),
    g_bounds = DEFAULT_G_BOUNDS,
    seed = 12345L
) {
  if (length(g_bounds) != 2L || g_bounds[1] <= 0 || g_bounds[2] >= 1 ||
      g_bounds[1] >= g_bounds[2]) {
    stop("g_bounds must satisfy 0 < lower < upper < 1")
  }

  fit_one <- function(df, Ycol, k) {
    ## Without this, SL's internal CV-fold assignment depends on whatever the
    ## ambient RNG state happens to be, which is not reproducible and --
    ## specifically for the misspecification comparison -- means Correct and
    ## Misspecified propensity fits could end up using DIFFERENT random fold
    ## splits, adding noise to a comparison that should differ only because
    ## of model specification, not because of unrelated CV-split randomness.
    if (use_sl) set.seed(as.integer(seed + 100L * k))
    fit_binary_nuisance_v7(
      Y = df[[Ycol]],
      X = ps_xvars_v7(df, k, misspecified = misspecified),
      use_sl = use_sl,
      SL.library = SL.library
    )
  }

  bound_stochastic <- function(x) {
    pmin(pmax(x, g_bounds[1]), g_bounds[2])
  }

  f1 <- fit_one(wide_df, "A1", 1L)

  idx2 <- wide_df$A1 == 0
  f2 <- fit_one(wide_df[idx2, , drop = FALSE], "A2", 2L)

  idx3 <- wide_df$A1 == 0 & wide_df$A2 == 0
  f3 <- fit_one(wide_df[idx3, , drop = FALSE], "A3", 3L)

  g1 <- bound_stochastic(
    predict_binary_nuisance_v7(f1, ps_xvars_v7(wide_df, 1L, misspecified = misspecified))
  )

  ## Once switched, treatment is deterministic and its probability must remain 1.
  g2 <- rep(1, nrow(wide_df))
  g2[idx2] <- bound_stochastic(
    predict_binary_nuisance_v7(
      f2,
      ps_xvars_v7(wide_df[idx2, , drop = FALSE], 2L, misspecified = misspecified)
    )
  )

  g3 <- rep(1, nrow(wide_df))
  g3[idx3] <- bound_stochastic(
    predict_binary_nuisance_v7(
      f3,
      ps_xvars_v7(wide_df[idx3, , drop = FALSE], 3L, misspecified = misspecified)
    )
  )

  stopifnot(all(g2[wide_df$A1 == 1] == 1))
  stopifnot(all(g3[wide_df$A1 == 1 | wide_df$A2 == 1] == 1))

  list(g1 = g1, g2 = g2, g3 = g3, fits = list(f1, f2, f3))
}

## Backward-compatible alias.
fit_ps <- fit_ps_v7

regime_probs_v7 <- function(wide_df, g1, g2, g3, regime) {
  regime <- match.arg(regime, c("d1", "d2"))

  if (regime == "d1") {
    A1s <- wide_df$Astar1_d1
    A2s <- wide_df$Astar2_d1
    A3s <- wide_df$Astar3_d1
  } else {
    A1s <- wide_df$Astar1_d2
    A2s <- wide_df$Astar2_d2
    A3s <- wide_df$Astar3_d2
  }

  p1 <- ifelse(A1s == 1, g1, 1 - g1)

  ## Later treatment probabilities contribute only while the observed subject
  ## remains consistent with the regime. After an earlier deviation the later
  ## factors are set to 1, because those decisions are outside that subject's
  ## regime-following path. This also prevents undefined 0/0 terms when an
  ## observed early switch makes later treatment deterministic but the regime
  ## had prescribed no switch at the earlier visit.
  c1_obs <- as.integer(wide_df$A1 == A1s)
  p2_raw <- ifelse(A2s == 1, g2, 1 - g2)
  p2 <- ifelse(c1_obs == 1, p2_raw, 1)

  c2_obs <- c1_obs * as.integer(wide_df$A2 == A2s)
  p3_raw <- ifelse(A3s == 1, g3, 1 - g3)
  p3 <- ifelse(c2_obs == 1, p3_raw, 1)

  list(
    A1s = A1s,
    A2s = A2s,
    A3s = A3s,
    p1 = p1,
    p2 = p2,
    p3 = p3,
    gbar1 = p1,
    gbar2 = p1 * p2,
    gbar3 = p1 * p2 * p3
  )
}

regime_probs <- regime_probs_v7

## ============================================================================
## 4. CORRECT ORACLE VIRAL-LOAD TRANSITION MODEL
## ============================================================================

build_transition_long_v7 <- function(wide_df) {
  bind_rows(
    wide_df %>%
      transmute(
        analysis_id,
        interval = 1L,
        logvl_current = logvl1,
        A = A1,
        slope_A,
        slope_B,
        logvl_next = logvl2
      ),
    wide_df %>%
      transmute(
        analysis_id,
        interval = 2L,
        logvl_current = logvl2,
        A = A2,
        slope_A,
        slope_B,
        logvl_next = logvl3
      ),
    wide_df %>%
      transmute(
        analysis_id,
        interval = 3L,
        logvl_current = logvl3,
        A = A3,
        slope_A,
        slope_B,
        logvl_next = logvl4
      )
  ) %>%
    mutate(
      A_slope_A = A * slope_A,
      A_slope_B = A * slope_B
    )
}

transition_x_v7 <- function(logvl_current, A, slope_A, slope_B, misspecified = FALSE) {
  if (!misspecified) {
    data.frame(
      logvl_current = logvl_current,
      slope_A = slope_A,
      A_slope_A = A * slope_A,
      A_slope_B = A * slope_B
    )
  } else {
    ## Moderated misspecification: correctly represents the first-line
    ## trajectory (retains slope_A and its treatment interaction) but
    ## replaces heterogeneous second-line response with one common additive
    ## switch effect -- omits ONLY the A:slope_B interaction, not both
    ## interactions. This is a more realistic degree of analyst error than
    ## dropping both interactions (which produced ~28-45% GCOMP bias, an
    ## extreme rather than moderate misspecification).
    data.frame(
      logvl_current = logvl_current,
      slope_A = slope_A,
      A_slope_A = A * slope_A,
      A = A
    )
  }
}

make_subject_level_valid_rows_v7 <- function(ids, V = 5L, seed = 12345L) {
  V <- max(2L, min(as.integer(V), length(unique(ids))))
  set.seed(seed)

  uid <- unique(ids)
  fold_id <- sample(rep(seq_len(V), length.out = length(uid)))
  fold_map <- setNames(fold_id, uid)
  row_fold <- unname(fold_map[as.character(ids)])

  lapply(seq_len(V), function(v) which(row_fold == v))
}

fit_oracle_transition_v7 <- function(
    wide_df,
    use_sl = FALSE,
    misspecified = FALSE,
    SL.library = c("SL.glm", "SL.gam"),
    V = 5L,
    seed = 12345L
) {
  dlong <- build_transition_long_v7(wide_df)
  X <- transition_x_v7(
    dlong$logvl_current,
    dlong$A,
    dlong$slope_A,
    dlong$slope_B,
    misspecified = misspecified
  )
  Y <- dlong$logvl_next

  if (!use_sl) {
    fit <- lm(Y ~ ., data = cbind(Y = Y, X))
    sigma_hat <- summary(fit)$sigma

    return(list(
      fit = fit,
      use_sl = FALSE,
      misspecified = misspecified,
      sigma = max(as.numeric(sigma_hat), 1e-6),
      library = "lm"
    ))
  }

  validRows <- make_subject_level_valid_rows_v7(
    ids = dlong$analysis_id,
    V = V,
    seed = seed
  )

  fit <- SuperLearner(
    Y = Y,
    X = X,
    family = gaussian(),
    SL.library = SL.library,
    cvControl = list(V = length(validRows), validRows = validRows)
  )

  ## Estimate the Gaussian innovation SD from cross-validated ensemble
  ## predictions when available, avoiding optimistic in-sample residuals.
  pred_cv <- NULL
  if (!is.null(fit$Z) && !is.null(fit$coef) &&
      nrow(fit$Z) == length(Y) && ncol(fit$Z) == length(fit$coef)) {
    pred_cv <- as.numeric(fit$Z %*% fit$coef)
  }

  if (is.null(pred_cv) || any(!is.finite(pred_cv))) {
    pred_cv <- as.numeric(
      predict(fit, newdata = X, onlySL = TRUE)$pred[, 1]
    )
  }

  sigma_hat <- sqrt(mean((Y - pred_cv)^2, na.rm = TRUE))

  list(
    fit = fit,
    use_sl = TRUE,
    misspecified = misspecified,
    sigma = max(as.numeric(sigma_hat), 1e-6),
    library = SL.library
  )
}

predict_transition_mu_v7 <- function(
    object,
    logvl_current,
    A,
    slope_A,
    slope_B,
    chunk_size = 100000L
) {
  n <- length(logvl_current)
  if (!all(c(length(A), length(slope_A), length(slope_B)) == n)) {
    stop("Transition prediction inputs must have equal lengths")
  }

  out <- numeric(n)
  starts <- seq.int(1L, n, by = as.integer(chunk_size))

  for (s in starts) {
    e <- min(n, s + as.integer(chunk_size) - 1L)
    idx <- s:e

    Xnew <- transition_x_v7(
      logvl_current[idx],
      A[idx],
      slope_A[idx],
      slope_B[idx],
      misspecified = isTRUE(object$misspecified)
    )

    if (!object$use_sl) {
      out[idx] <- as.numeric(
        predict(object$fit, newdata = Xnew)
      )
    } else {
      out[idx] <- as.numeric(
        predict(object$fit, newdata = Xnew, onlySL = TRUE)$pred[, 1]
      )
    }
  }

  out
}

## ============================================================================
## 5. ORACLE DYNAMIC-REGIME Q ENGINE
## ============================================================================

regime_action_v7 <- function(
    logvl,
    prior_switched,
    regime,
    dtr1_threshold = 50,
    dtr2_threshold = 1000
) {
  regime <- match.arg(regime, c("d1", "d2"))

  new_switch <- if (regime == "d1") {
    logvl >= log10(dtr1_threshold)
  } else {
    logvl > log10(dtr2_threshold)
  }

  as.integer(prior_switched == 1 | new_switch)
}

oracle_q3_state_v7 <- function(
    logvl3,
    prior_switched,
    slope_A,
    slope_B,
    transition_fit,
    regime,
    epsilon3 = 0,
    outcome_threshold = 50
) {
  A3d <- regime_action_v7(
    logvl = logvl3,
    prior_switched = prior_switched,
    regime = regime
  )

  mu4 <- predict_transition_mu_v7(
    transition_fit,
    logvl_current = logvl3,
    A = A3d,
    slope_A = slope_A,
    slope_B = slope_B
  )

  q3 <- pnorm(
    (log10(outcome_threshold) - mu4) / transition_fit$sigma
  )

  apply_logit_shift(q3, epsilon3)
}

oracle_q2_state_v7 <- function(
    logvl2,
    prior_switched,
    slope_A,
    slope_B,
    transition_fit,
    regime,
    quadrature,
    epsilon3 = 0,
    epsilon2 = 0,
    chunk_size = 500L
) {
  n <- length(logvl2)
  out <- numeric(n)
  K <- length(quadrature$nodes)

  starts <- seq.int(1L, n, by = as.integer(chunk_size))

  for (s in starts) {
    e <- min(n, s + as.integer(chunk_size) - 1L)
    idx <- s:e
    m <- length(idx)

    A2d <- regime_action_v7(
      logvl = logvl2[idx],
      prior_switched = prior_switched[idx],
      regime = regime
    )

    mu3 <- predict_transition_mu_v7(
      transition_fit,
      logvl_current = logvl2[idx],
      A = A2d,
      slope_A = slope_A[idx],
      slope_B = slope_B[idx]
    )

    logvl3_grid <-
      rep(mu3, each = K) +
      transition_fit$sigma * rep(quadrature$nodes, times = m)

    q3_grid <- oracle_q3_state_v7(
      logvl3 = logvl3_grid,
      prior_switched = rep(A2d, each = K),
      slope_A = rep(slope_A[idx], each = K),
      slope_B = rep(slope_B[idx], each = K),
      transition_fit = transition_fit,
      regime = regime,
      epsilon3 = epsilon3
    )

    q3_matrix <- matrix(q3_grid, nrow = m, ncol = K, byrow = TRUE)
    q2 <- as.numeric(q3_matrix %*% quadrature$weights)

    out[idx] <- apply_logit_shift(q2, epsilon2)
  }

  out
}

oracle_q1_state_v7 <- function(
    logvl1,
    slope_A,
    slope_B,
    transition_fit,
    regime,
    quadrature,
    epsilon3 = 0,
    epsilon2 = 0,
    epsilon1 = 0,
    chunk_size = 150L
) {
  n <- length(logvl1)
  out <- numeric(n)
  K <- length(quadrature$nodes)

  starts <- seq.int(1L, n, by = as.integer(chunk_size))

  for (s in starts) {
    e <- min(n, s + as.integer(chunk_size) - 1L)
    idx <- s:e
    m <- length(idx)

    A1d <- regime_action_v7(
      logvl = logvl1[idx],
      prior_switched = rep(0L, m),
      regime = regime
    )

    mu2 <- predict_transition_mu_v7(
      transition_fit,
      logvl_current = logvl1[idx],
      A = A1d,
      slope_A = slope_A[idx],
      slope_B = slope_B[idx]
    )

    logvl2_grid <-
      rep(mu2, each = K) +
      transition_fit$sigma * rep(quadrature$nodes, times = m)

    q2_grid <- oracle_q2_state_v7(
      logvl2 = logvl2_grid,
      prior_switched = rep(A1d, each = K),
      slope_A = rep(slope_A[idx], each = K),
      slope_B = rep(slope_B[idx], each = K),
      transition_fit = transition_fit,
      regime = regime,
      quadrature = quadrature,
      epsilon3 = epsilon3,
      epsilon2 = epsilon2,
      chunk_size = max(100L, chunk_size * K)
    )

    q2_matrix <- matrix(q2_grid, nrow = m, ncol = K, byrow = TRUE)
    q1 <- as.numeric(q2_matrix %*% quadrature$weights)

    out[idx] <- apply_logit_shift(q1, epsilon1)
  }

  out
}

oracle_q_bundle_v7 <- function(
    wide_df,
    transition_fit,
    regime,
    quadrature = normal_quadrature(DEFAULT_N_QUAD)
) {
  regime <- match.arg(regime, c("d1", "d2"))

  A1d <- regime_action_v7(
    wide_df$logvl1,
    prior_switched = rep(0L, nrow(wide_df)),
    regime = regime
  )

  A2d <- regime_action_v7(
    wide_df$logvl2,
    prior_switched = A1d,
    regime = regime
  )

  q3 <- oracle_q3_state_v7(
    logvl3 = wide_df$logvl3,
    prior_switched = A2d,
    slope_A = wide_df$slope_A,
    slope_B = wide_df$slope_B,
    transition_fit = transition_fit,
    regime = regime
  )

  q2 <- oracle_q2_state_v7(
    logvl2 = wide_df$logvl2,
    prior_switched = A1d,
    slope_A = wide_df$slope_A,
    slope_B = wide_df$slope_B,
    transition_fit = transition_fit,
    regime = regime,
    quadrature = quadrature
  )

  q1 <- oracle_q1_state_v7(
    logvl1 = wide_df$logvl1,
    slope_A = wide_df$slope_A,
    slope_B = wide_df$slope_B,
    transition_fit = transition_fit,
    regime = regime,
    quadrature = quadrature
  )

  list(q1 = q1, q2 = q2, q3 = q3)
}

## ============================================================================
## 6. IPTW-MSM WITH NODE-SPECIFIC 0.01/0.99 PROBABILITY BOUNDS
## ============================================================================

iptw_msm_v7 <- function(wide_df, g1, g2, g3) {
  build_long <- function(regime) {
    rp <- regime_probs_v7(wide_df, g1, g2, g3, regime)

    I <-
      as.integer(wide_df$A1 == rp$A1s) *
      as.integer(wide_df$A2 == rp$A2s) *
      as.integer(wide_df$A3 == rp$A3s)

    denom <- rp$gbar3
    w <- I / denom

    data.frame(
      analysis_id = wide_df$analysis_id,
      d1 = as.integer(regime == "d1"),
      Y = wide_df$Y4,
      weight = w
    )
  }

  msm <- bind_rows(build_long("d1"), build_long("d2"))

  w1 <- msm$weight[msm$d1 == 1]
  w2 <- msm$weight[msm$d1 == 0]

  ess <- function(w) {
    if (!any(w > 0)) return(NA_real_)
    sum(w)^2 / sum(w^2)
  }

  w99 <- function(w) {
    if (!any(w > 0)) return(NA_real_)
    as.numeric(quantile(w[w > 0], 0.99, na.rm = TRUE))
  }

  extra1 <- extra_weight_diagnostics_v7(w1)
  extra2 <- extra_weight_diagnostics_v7(w2)

  ## Zero-weight rows contribute nothing to estimation or variance -- removing
  ## them prevents thousands of harmless "zero weight" dispersion warnings
  ## from flooding HPC logs across a 12,000-replicate run. ESS/w99 above are
  ## still computed from the full (unfiltered) msm, since those diagnostics
  ## need to see the zero-weight rows to be meaningful.
  msm_fit <- msm[is.finite(msm$weight) & msm$weight > 0, , drop = FALSE]

  des <- survey::svydesign(
    ids = ~analysis_id,
    weights = ~weight,
    data = msm_fit
  )

  fit <- tryCatch(
    survey::svyglm(
      Y ~ d1,
      design = des,
      family = quasibinomial(link = "log")
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    ans <- safe_ratio_result()
    ans["ess_d1"] <- ess(w1)
    ans["ess_d2"] <- ess(w2)
    ans["w99_d1"] <- w99(w1)
    ans["w99_d2"] <- w99(w2)
    ans["mean_positive_weight_d1"] <- extra1["mean_positive"]; ans["mean_positive_weight_d2"] <- extra2["mean_positive"]
    ans["sd_positive_weight_d1"] <- extra1["sd_positive"]; ans["sd_positive_weight_d2"] <- extra2["sd_positive"]
    ans["cv_d1"] <- extra1["cv_positive"]; ans["cv_d2"] <- extra2["cv_positive"]
    ans["max_weight_d1"] <- extra1["max_positive"]; ans["max_weight_d2"] <- extra2["max_positive"]
    ans["pct_positive_d1"] <- extra1["pct_positive"]; ans["pct_positive_d2"] <- extra2["pct_positive"]
    return(ans)
  }

  co <- coef(summary(fit))
  log_rr <- unname(co["d1", "Estimate"])
  se_log <- unname(co["d1", "Std. Error"])
  rr <- exp(log_rr)
  EY_d2 <- exp(unname(co["(Intercept)", "Estimate"]))
  EY_d1 <- EY_d2 * rr

  c(
    est = rr,
    se = rr * se_log,
    lcl = exp(log_rr - 1.96 * se_log),
    ucl = exp(log_rr + 1.96 * se_log),
    EY_d1 = EY_d1,
    EY_d2 = EY_d2,
    ess_d1 = ess(w1),
    ess_d2 = ess(w2),
    w99_d1 = w99(w1),
    w99_d2 = w99(w2),
    mean_positive_weight_d1 = unname(extra1["mean_positive"]), mean_positive_weight_d2 = unname(extra2["mean_positive"]),
    sd_positive_weight_d1 = unname(extra1["sd_positive"]), sd_positive_weight_d2 = unname(extra2["sd_positive"]),
    cv_d1 = unname(extra1["cv_positive"]), cv_d2 = unname(extra2["cv_positive"]),
    max_weight_d1 = unname(extra1["max_positive"]), max_weight_d2 = unname(extra2["max_positive"]),
    pct_positive_d1 = unname(extra1["pct_positive"]), pct_positive_d2 = unname(extra2["pct_positive"]),
    se_log_rr = se_log
  )
}


## ----------------------------------------------------------------------------
## IPTW sensitivity estimator with two-sided final-weight truncation
##
## Positive inverse-probability weights are winsorized separately within each
## regime at the 1st and 99th percentiles. Zero weights remain zero because
## they identify individuals who did not follow that regime. The original
## iptw_msm_v7() above remains the raw, untruncated primary IPTW estimator.
## ----------------------------------------------------------------------------

.iptw_positive_ess_v7b <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L || sum(w^2) <= 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

.iptw_positive_quantile_v7b <- function(w, p) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L) return(NA_real_)
  as.numeric(stats::quantile(w, probs = p, na.rm = TRUE,
                             names = FALSE, type = 7))
}

safe_iptw_trunc01_99_result_v7b <- function() {
  c(
    est = NA_real_, se = NA_real_, lcl = NA_real_, ucl = NA_real_,
    EY_d1 = NA_real_, EY_d2 = NA_real_,

    ess_d1 = NA_real_, w01_d1 = NA_real_, w99_d1 = NA_real_,
    raw_ess_d1 = NA_real_, raw_w01_d1 = NA_real_, raw_w99_d1 = NA_real_,
    lower_cutoff_d1 = NA_real_, upper_cutoff_d1 = NA_real_,
    pct_lower_truncated_d1 = NA_real_, pct_upper_truncated_d1 = NA_real_,
    pct_any_truncated_d1 = NA_real_,

    ess_d2 = NA_real_, w01_d2 = NA_real_, w99_d2 = NA_real_,
    raw_ess_d2 = NA_real_, raw_w01_d2 = NA_real_, raw_w99_d2 = NA_real_,
    lower_cutoff_d2 = NA_real_, upper_cutoff_d2 = NA_real_,
    pct_lower_truncated_d2 = NA_real_, pct_upper_truncated_d2 = NA_real_,
    pct_any_truncated_d2 = NA_real_,

    se_log_rr = NA_real_
  )
}

iptw_msm_trunc01_99_v7 <- function(
    wide_df,
    g1,
    g2,
    g3,
    lower_quantile = 0.01,
    upper_quantile = 0.99
) {
  if (length(lower_quantile) != 1L || length(upper_quantile) != 1L ||
      !is.finite(lower_quantile) || !is.finite(upper_quantile) ||
      lower_quantile < 0 || lower_quantile >= 0.5 ||
      upper_quantile <= 0.5 || upper_quantile > 1 ||
      lower_quantile >= upper_quantile) {
    stop("Require 0 <= lower_quantile < 0.5 < upper_quantile <= 1.")
  }

  build_long <- function(regime) {
    rp <- regime_probs_v7(wide_df, g1, g2, g3, regime)

    I <-
      as.integer(wide_df$A1 == rp$A1s) *
      as.integer(wide_df$A2 == rp$A2s) *
      as.integer(wide_df$A3 == rp$A3s)

    denom <- rp$gbar3
    raw_weight <- ifelse(
      is.finite(denom) & denom > 0,
      I / denom,
      NA_real_
    )

    data.frame(
      analysis_id = wide_df$analysis_id,
      d1 = as.integer(regime == "d1"),
      Y = wide_df$Y4,
      raw_weight = raw_weight,
      weight = raw_weight
    )
  }

  msm <- bind_rows(build_long("d1"), build_long("d2"))
  diagnostics <- list()

  for (regime_value in c(1L, 0L)) {
    suffix <- if (regime_value == 1L) "d1" else "d2"
    idx <- which(
      msm$d1 == regime_value &
        is.finite(msm$raw_weight) & msm$raw_weight > 0
    )

    lower_cutoff <- NA_real_
    upper_cutoff <- NA_real_
    pct_lower <- NA_real_
    pct_upper <- NA_real_
    pct_any <- NA_real_

    if (length(idx) > 0L) {
      raw_positive <- msm$raw_weight[idx]
      lower_cutoff <- .iptw_positive_quantile_v7b(
        raw_positive, lower_quantile
      )
      upper_cutoff <- .iptw_positive_quantile_v7b(
        raw_positive, upper_quantile
      )

      below <- raw_positive < lower_cutoff
      above <- raw_positive > upper_cutoff

      msm$weight[idx] <- pmin(
        pmax(raw_positive, lower_cutoff),
        upper_cutoff
      )

      pct_lower <- 100 * mean(below)
      pct_upper <- 100 * mean(above)
      pct_any <- 100 * mean(below | above)
    }

    used_w <- msm$weight[msm$d1 == regime_value]
    raw_w <- msm$raw_weight[msm$d1 == regime_value]

    diagnostics[[paste0("ess_", suffix)]] <-
      .iptw_positive_ess_v7b(used_w)
    diagnostics[[paste0("w01_", suffix)]] <-
      .iptw_positive_quantile_v7b(used_w, 0.01)
    diagnostics[[paste0("w99_", suffix)]] <-
      .iptw_positive_quantile_v7b(used_w, 0.99)

    diagnostics[[paste0("raw_ess_", suffix)]] <-
      .iptw_positive_ess_v7b(raw_w)
    diagnostics[[paste0("raw_w01_", suffix)]] <-
      .iptw_positive_quantile_v7b(raw_w, 0.01)
    diagnostics[[paste0("raw_w99_", suffix)]] <-
      .iptw_positive_quantile_v7b(raw_w, 0.99)

    extra_used <- extra_weight_diagnostics_v7(used_w)
    diagnostics[[paste0("mean_positive_weight_", suffix)]] <- unname(extra_used["mean_positive"])
    diagnostics[[paste0("sd_positive_weight_", suffix)]] <- unname(extra_used["sd_positive"])
    diagnostics[[paste0("cv_", suffix)]] <- unname(extra_used["cv_positive"])
    diagnostics[[paste0("max_weight_", suffix)]] <- unname(extra_used["max_positive"])
    diagnostics[[paste0("pct_positive_", suffix)]] <- unname(extra_used["pct_positive"])

    ## Pre-truncation versions -- previously only ESS/w01/w99 were saved for
    ## the raw (pre-truncation) weights; mean/SD/CV/max were missing.
    extra_raw <- extra_weight_diagnostics_v7(raw_w)
    diagnostics[[paste0("pretrunc_mean_positive_weight_", suffix)]] <- unname(extra_raw["mean_positive"])
    diagnostics[[paste0("pretrunc_sd_positive_weight_", suffix)]] <- unname(extra_raw["sd_positive"])
    diagnostics[[paste0("pretrunc_cv_", suffix)]] <- unname(extra_raw["cv_positive"])
    diagnostics[[paste0("pretrunc_max_weight_", suffix)]] <- unname(extra_raw["max_positive"])
    ## pct_positive not repeated pre/post -- truncation (winsorizing) doesn't
    ## change which weights are structurally zero, so this value is identical
    ## before and after; storing it once (post-truncation, above) is sufficient.

    diagnostics[[paste0("lower_cutoff_", suffix)]] <- lower_cutoff
    diagnostics[[paste0("upper_cutoff_", suffix)]] <- upper_cutoff
    diagnostics[[paste0("pct_lower_truncated_", suffix)]] <- pct_lower
    diagnostics[[paste0("pct_upper_truncated_", suffix)]] <- pct_upper
    diagnostics[[paste0("pct_any_truncated_", suffix)]] <- pct_any
  }

  diagnostic_values <- unlist(diagnostics, use.names = TRUE)
  msm_fit <- msm[
    is.finite(msm$weight) & msm$weight > 0 &
      is.finite(msm$Y) & is.finite(msm$d1),
    , drop = FALSE
  ]

  if (nrow(msm_fit) == 0L || length(unique(msm_fit$d1)) < 2L) {
    ans <- safe_iptw_trunc01_99_result_v7b()
    ans[names(diagnostic_values)] <- diagnostic_values
    return(ans)
  }

  des <- survey::svydesign(
    ids = ~analysis_id,
    weights = ~weight,
    data = msm_fit
  )

  fit <- tryCatch(
    survey::svyglm(
      Y ~ d1,
      design = des,
      family = quasibinomial(link = "log")
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    ans <- safe_iptw_trunc01_99_result_v7b()
    ans[names(diagnostic_values)] <- diagnostic_values
    return(ans)
  }

  co <- coef(summary(fit))
  if (!all(c("(Intercept)", "d1") %in% rownames(co))) {
    ans <- safe_iptw_trunc01_99_result_v7b()
    ans[names(diagnostic_values)] <- diagnostic_values
    return(ans)
  }

  log_rr <- unname(co["d1", "Estimate"])
  se_log <- unname(co["d1", "Std. Error"])
  log_mu_d2 <- unname(co["(Intercept)", "Estimate"])

  if (!all(is.finite(c(log_rr, se_log, log_mu_d2)))) {
    ans <- safe_iptw_trunc01_99_result_v7b()
    ans[names(diagnostic_values)] <- diagnostic_values
    return(ans)
  }

  rr <- exp(log_rr)
  EY_d2 <- exp(log_mu_d2)
  EY_d1 <- EY_d2 * rr

  c(
    est = rr,
    se = rr * se_log,
    lcl = exp(log_rr - 1.96 * se_log),
    ucl = exp(log_rr + 1.96 * se_log),
    EY_d1 = EY_d1,
    EY_d2 = EY_d2,
    diagnostic_values,
    se_log_rr = se_log
  )
}

## Backward-compatible name. The old stabilization/99th-weight-quantile
## arguments are deliberately unsupported in the v7 primary estimator.
iptw_msm <- function(
    wide_df,
    g1,
    g2,
    g3,
    stabilize = FALSE,
    trunc_q = NULL
) {
  if (isTRUE(stabilize) || !is.null(trunc_q)) {
    warning(
      "v7 primary IPTW uses node-specific g bounds [0.01,0.99] and no ",
      "additional weight-percentile truncation; stabilize/trunc_q ignored."
    )
  }
  iptw_msm_v7(wide_df, g1, g2, g3)
}

## ============================================================================
## 7. ORACLE G-COMPUTATION
## ============================================================================

gcomp_one_regime_v7 <- function(
    wide_df,
    transition_fit,
    regime,
    quadrature
) {
  q1 <- oracle_q1_state_v7(
    logvl1 = wide_df$logvl1,
    slope_A = wide_df$slope_A,
    slope_B = wide_df$slope_B,
    transition_fit = transition_fit,
    regime = regime,
    quadrature = quadrature
  )

  mean(q1)
}

gcomp_rr_v7 <- function(
    wide_df,
    transition_fit,
    quadrature = normal_quadrature(DEFAULT_N_QUAD)
) {
  EY_d1 <- gcomp_one_regime_v7(
    wide_df, transition_fit, "d1", quadrature
  )
  EY_d2 <- gcomp_one_regime_v7(
    wide_df, transition_fit, "d2", quadrature
  )

  c(
    est = EY_d1 / EY_d2,
    se = NA_real_,
    lcl = NA_real_,
    ucl = NA_real_,
    EY_d1 = EY_d1,
    EY_d2 = EY_d2,
    ess_d1 = NA_real_,
    ess_d2 = NA_real_,
    w99_d1 = NA_real_,
    w99_d2 = NA_real_,
    se_log_rr = NA_real_
  )
}

bootstrap_gcomp_oracle_v7 <- function(
    wide_df,
    use_sl = FALSE,
    misspecified = FALSE,
    SL.library = c("SL.glm", "SL.gam"),
    B = 0L,
    seed = 1L,
    n_quad = DEFAULT_N_QUAD,
    min_valid_fraction = 0.80
) {
  quadrature <- normal_quadrature(n_quad)

  fit0 <- fit_oracle_transition_v7(
    wide_df,
    use_sl = use_sl,
    misspecified = misspecified,
    SL.library = SL.library,
    seed = seed
  )
  point <- gcomp_rr_v7(wide_df, fit0, quadrature)

  if (B <= 0L || !is.finite(point["est"])) {
    return(point)
  }

  set.seed(seed)
  n <- nrow(wide_df)

  ## Generate ALL resamples before any SL fit changes the RNG state -- avoids
  ## the resampling draws and each replicate's internal SL-CV-fold seed from
  ## being entangled in one shared, hard-to-reason-about RNG sequence.
  boot_indices <- replicate(B, sample.int(n, n, replace = TRUE), simplify = FALSE)
  boot_log_rr <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    wb <- wide_df[boot_indices[[b]], , drop = FALSE]
    ## Deliberately keep the ORIGINAL analysis_id, not a fresh sequential one.
    ## Repeated copies of the same original subject (inevitable under
    ## resampling with replacement) must stay in the same subject-level SL
    ## validation fold -- otherwise a duplicate row can appear in both the
    ## training and validation folds simultaneously, since it's literally
    ## identical data, producing overly optimistic cross-validated performance.

    est_b <- tryCatch({
      fit_b <- fit_oracle_transition_v7(
        wb,
        use_sl = use_sl,
        misspecified = misspecified,
        SL.library = SL.library,
        seed = seed + 100000L + b  ## offset well clear of the resampling seed range
      )
      unname(gcomp_rr_v7(wb, fit_b, quadrature)["est"])
    }, error = function(e) NA_real_)

    if (is.finite(est_b) && est_b > 0) {
      boot_log_rr[b] <- log(est_b)
    }
  }

  valid <- is.finite(boot_log_rr)
  n_valid <- sum(valid)
  required <- max(2L, ceiling(min_valid_fraction * B))

  if (n_valid < required) {
    return(c(point, n_boot_valid = n_valid))
  }

  se_log <- sd(boot_log_rr[valid])
  rr <- unname(point["est"])

  point["se"] <- rr * se_log
  point["lcl"] <- exp(log(rr) - 1.96 * se_log)
  point["ucl"] <- exp(log(rr) + 1.96 * se_log)
  point["se_log_rr"] <- se_log

  c(point, n_boot_valid = n_valid)
}

## ============================================================================
## 8. ORACLE LONGITUDINAL AIPTW
## ============================================================================

aiptw_one_regime_v7 <- function(
    wide_df,
    g1,
    g2,
    g3,
    transition_fit,
    regime,
    quadrature
) {
  rp <- regime_probs_v7(wide_df, g1, g2, g3, regime)

  c1 <- as.integer(wide_df$A1 == rp$A1s)
  c2 <- c1 * as.integer(wide_df$A2 == rp$A2s)
  c3 <- c2 * as.integer(wide_df$A3 == rp$A3s)

  q <- oracle_q_bundle_v7(
    wide_df,
    transition_fit,
    regime,
    quadrature
  )

  phi <-
    q$q1 +
    (c1 / rp$gbar1) * (q$q2 - q$q1) +
    (c2 / rp$gbar2) * (q$q3 - q$q2) +
    (c3 / rp$gbar3) * (wide_df$Y4 - q$q3)

  psi <- mean(phi)

  list(
    EY_hat = psi,
    IC = phi - psi,
    q = q
  )
}

aiptw_rr_v7 <- function(
    wide_df,
    g1,
    g2,
    g3,
    transition_fit,
    quadrature = normal_quadrature(DEFAULT_N_QUAD)
) {
  r1 <- aiptw_one_regime_v7(
    wide_df, g1, g2, g3, transition_fit, "d1", quadrature
  )
  r2 <- aiptw_one_regime_v7(
    wide_df, g1, g2, g3, transition_fit, "d2", quadrature
  )

  EY_d1 <- r1$EY_hat
  EY_d2 <- r2$EY_hat
  rr <- EY_d1 / EY_d2

  IC_rr <-
    r1$IC / EY_d2 -
    EY_d1 * r2$IC / EY_d2^2

  se_rr <- sqrt(var(IC_rr) / nrow(wide_df))
  se_log <- se_rr / rr

  c(
    est = rr,
    se = se_rr,
    lcl = exp(log(rr) - 1.96 * se_log),
    ucl = exp(log(rr) + 1.96 * se_log),
    EY_d1 = EY_d1,
    EY_d2 = EY_d2,
    ess_d1 = NA_real_,
    ess_d2 = NA_real_,
    w99_d1 = NA_real_,
    w99_d2 = NA_real_,
    se_log_rr = se_log
  )
}

## ============================================================================
## 9. ORACLE LONGITUDINAL TMLE
## ============================================================================

## Solves the TMLE fluctuation's score equation directly via uniroot() rather
## than fitting a glm() with an offset. Two real advantages over the glm
## approach: (1) uniroot's bisection search doesn't share glm/IRLS's
## convergence failure modes under pathological data (e.g. a badly wrong
## initial Q from a misspecified transition model combined with extreme
## weights under severe positivity) -- this was very likely the actual cause
## of the billion-scale TMLE estimates seen under Qmis/Both; (2) on failure
## this returns NA, which the caller must treat as an invalid replicate --
## NOT silently fall back to epsilon=0, which would quietly produce an
## untargeted, biased "TMLE" estimate with no indication anything went wrong.
fit_fluctuation_intercept_v7 <- function(y, qinit, weights, search_limit = 40) {
  use <-
    is.finite(y) &
    is.finite(qinit) &
    is.finite(weights) &
    weights > 0

  if (!any(use)) return(NA_real_)

  y_use <- y[use]
  q_use <- clamp_prob(qinit[use])

  w <- weights[use]
  w <- w / mean(w)

  eta0 <- qlogis(q_use)

  ## The TMLE fluctuation score equation: sum of weighted residuals under the
  ## one-parameter logistic submodel Q(epsilon) = expit(eta0 + epsilon).
  score <- function(epsilon) {
    sum(w * (y_use - plogis(eta0 + epsilon)))
  }

  lower_score <- score(-search_limit)
  upper_score <- score(search_limit)

  if (!is.finite(lower_score) || !is.finite(upper_score) || lower_score * upper_score > 0) {
    return(NA_real_)
  }

  uniroot(score, interval = c(-search_limit, search_limit), tol = 1e-10)$root
}

tmle_one_regime_v7 <- function(
    wide_df,
    g1,
    g2,
    g3,
    transition_fit,
    regime,
    quadrature
) {
  rp <- regime_probs_v7(wide_df, g1, g2, g3, regime)

  c1 <- as.integer(wide_df$A1 == rp$A1s)
  c2 <- c1 * as.integer(wide_df$A2 == rp$A2s)
  c3 <- c2 * as.integer(wide_df$A3 == rp$A3s)

  A1d <- regime_action_v7(
    wide_df$logvl1,
    prior_switched = rep(0L, nrow(wide_df)),
    regime = regime
  )
  A2d <- regime_action_v7(
    wide_df$logvl2,
    prior_switched = A1d,
    regime = regime
  )

  ## Stage 3
  q3_init <- oracle_q3_state_v7(
    logvl3 = wide_df$logvl3,
    prior_switched = A2d,
    slope_A = wide_df$slope_A,
    slope_B = wide_df$slope_B,
    transition_fit = transition_fit,
    regime = regime,
    epsilon3 = 0
  )

  eps3 <- fit_fluctuation_intercept_v7(
    y = wide_df$Y4,
    qinit = q3_init,
    weights = c3 / rp$gbar3
  )
  if (!is.finite(eps3)) stop("TMLE stage-3 fluctuation failed")
  q3_star <- apply_logit_shift(q3_init, eps3)

  ## Stage 2
  q2_init <- oracle_q2_state_v7(
    logvl2 = wide_df$logvl2,
    prior_switched = A1d,
    slope_A = wide_df$slope_A,
    slope_B = wide_df$slope_B,
    transition_fit = transition_fit,
    regime = regime,
    quadrature = quadrature,
    epsilon3 = eps3,
    epsilon2 = 0
  )

  eps2 <- fit_fluctuation_intercept_v7(
    y = q3_star,
    qinit = q2_init,
    weights = c2 / rp$gbar2
  )
  if (!is.finite(eps2)) stop("TMLE stage-2 fluctuation failed")
  q2_star <- apply_logit_shift(q2_init, eps2)

  ## Stage 1
  q1_init <- oracle_q1_state_v7(
    logvl1 = wide_df$logvl1,
    slope_A = wide_df$slope_A,
    slope_B = wide_df$slope_B,
    transition_fit = transition_fit,
    regime = regime,
    quadrature = quadrature,
    epsilon3 = eps3,
    epsilon2 = eps2,
    epsilon1 = 0
  )

  eps1 <- fit_fluctuation_intercept_v7(
    y = q2_star,
    qinit = q1_init,
    weights = c1 / rp$gbar1
  )
  if (!is.finite(eps1)) stop("TMLE stage-1 fluctuation failed")
  q1_star <- apply_logit_shift(q1_init, eps1)

  psi <- mean(q1_star)

  IC <-
    (c3 / rp$gbar3) * (wide_df$Y4 - q3_star) +
    (c2 / rp$gbar2) * (q3_star - q2_star) +
    (c1 / rp$gbar1) * (q2_star - q1_star) +
    q1_star - psi

  list(
    EY_hat = psi,
    IC = IC,
    eps = c(eps1 = eps1, eps2 = eps2, eps3 = eps3)
  )
}

tmle_rr_v7 <- function(
    wide_df,
    g1,
    g2,
    g3,
    transition_fit,
    quadrature = normal_quadrature(DEFAULT_N_QUAD)
) {
  r1 <- tmle_one_regime_v7(
    wide_df, g1, g2, g3, transition_fit, "d1", quadrature
  )
  r2 <- tmle_one_regime_v7(
    wide_df, g1, g2, g3, transition_fit, "d2", quadrature
  )

  EY_d1 <- r1$EY_hat
  EY_d2 <- r2$EY_hat
  rr <- EY_d1 / EY_d2

  IC_rr <-
    r1$IC / EY_d2 -
    EY_d1 * r2$IC / EY_d2^2

  se_rr <- sqrt(var(IC_rr) / nrow(wide_df))
  se_log <- se_rr / rr

  c(
    est = rr,
    se = se_rr,
    lcl = exp(log(rr) - 1.96 * se_log),
    ucl = exp(log(rr) + 1.96 * se_log),
    EY_d1 = EY_d1,
    EY_d2 = EY_d2,
    ess_d1 = NA_real_,
    ess_d2 = NA_real_,
    w99_d1 = NA_real_,
    w99_d2 = NA_real_,
    se_log_rr = se_log
  )
}

## ============================================================================
## 10. SUPER LEARNER LIBRARIES
## ============================================================================

stable_sl_library_v7 <- function() {
  ## The exact GLM is included. GAM supplies a flexible alternative.
  c("SL.glm", "SL.gam")
}

stable_sl_library_v6 <- stable_sl_library_v7

## ============================================================================
## 11. ONE SIMULATION REPLICATE
## ============================================================================

## ============================================================================
## v7c ADD-ON: STABILIZED IPTW FOR DYNAMIC TREATMENT REGIMES
## (merged into v7e -- previously a separate sourced-after file)
##
## Adds two estimators:
##   IPTW_STAB               : stabilized dynamic-regime IPTW
##   IPTW_STAB_TRUNC01_99    : the same stabilized weights, winsorized at the
##                             1st and 99th percentiles within each regime
##
## Stabilization is implemented through the artificial-censoring/cloning view
## of dynamic regimes. The denominator is the existing probability of following
## the regime given full observed history. The numerator is the probability of
## remaining adherent given only regime and baseline covariates. Because those
## baseline variables enter the numerator, they are also included in the
## weighted outcome model, and the regime-specific risks are standardized back
## to the empirical baseline population.
## ============================================================================

safe_iptw_stab_result_v7c <- function(truncated = FALSE) {
  out <- c(
    est = NA_real_, se = NA_real_, lcl = NA_real_, ucl = NA_real_,
    EY_d1 = NA_real_, EY_d2 = NA_real_,
    max_predicted_risk_d1 = NA_real_, max_predicted_risk_d2 = NA_real_,
    pct_predicted_risk_gt1_d1 = NA_real_, pct_predicted_risk_gt1_d2 = NA_real_,
    ess_d1 = NA_real_, ess_d2 = NA_real_,
    w99_d1 = NA_real_, w99_d2 = NA_real_,
    se_log_rr = NA_real_,
    w01_d1 = NA_real_, w01_d2 = NA_real_,
    mean_positive_weight_d1 = NA_real_,
    mean_positive_weight_d2 = NA_real_
  )

  if (isTRUE(truncated)) {
    out <- c(
      out,
      pretrunc_ess_d1 = NA_real_, pretrunc_ess_d2 = NA_real_,
      pretrunc_w01_d1 = NA_real_, pretrunc_w01_d2 = NA_real_,
      pretrunc_w99_d1 = NA_real_, pretrunc_w99_d2 = NA_real_,
      lower_cutoff_d1 = NA_real_, lower_cutoff_d2 = NA_real_,
      upper_cutoff_d1 = NA_real_, upper_cutoff_d2 = NA_real_,
      pct_lower_truncated_d1 = NA_real_, pct_lower_truncated_d2 = NA_real_,
      pct_upper_truncated_d1 = NA_real_, pct_upper_truncated_d2 = NA_real_,
      pct_any_truncated_d1 = NA_real_, pct_any_truncated_d2 = NA_real_
    )
  }

  out
}

weight_ess_v7c <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

weight_quantile_v7c <- function(w, prob) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L) return(NA_real_)
  as.numeric(stats::quantile(w, probs = prob, na.rm = TRUE, names = FALSE))
}

weight_mean_positive_v7c <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L) return(NA_real_)
  mean(w)
}

## Consolidated weight diagnostics -- mean, SD, CV, max, and % positive all
## computed together from one pass, replacing the earlier split where mean
## came from a separate helper (weight_mean_positive_v7c) and CV/max/
## pct_positive came from here. This is a targeted, explicitly-requested
## consolidation of these specific fields, not the broader helper-function
## refactor deferred until after the current HPC run.
extra_weight_diagnostics_v7 <- function(w_full) {
  finite <- w_full[is.finite(w_full)]
  pos <- finite[finite > 0]

  mean_pos <- if (length(pos) > 0L) mean(pos) else NA_real_
  sd_pos <- if (length(pos) > 1L) stats::sd(pos) else NA_real_

  c(
    mean_positive = mean_pos,
    sd_positive = sd_pos,
    cv_positive = if (is.finite(mean_pos) && mean_pos > 0 && is.finite(sd_pos)) {
      sd_pos / mean_pos
    } else {
      NA_real_
    },
    max_positive = if (length(pos) > 0L) max(pos) else NA_real_,
    pct_positive = if (length(finite) > 0L) 100 * mean(finite > 0) else NA_real_
  )
}

## Build two clones per subject, one for each dynamic regime, and identify
## whether each clone remains adherent at each decision time.
build_regime_clones_v7c <- function(wide_df) {
  build_one <- function(regime) {
    if (regime == "d1") {
      a1s <- wide_df$Astar1_d1
      a2s <- wide_df$Astar2_d1
      a3s <- wide_df$Astar3_d1
    } else {
      a1s <- wide_df$Astar1_d2
      a2s <- wide_df$Astar2_d2
      a3s <- wide_df$Astar3_d2
    }

    adhere1 <- as.integer(wide_df$A1 == a1s)
    adhere2_now <- as.integer(wide_df$A2 == a2s)
    adhere3_now <- as.integer(wide_df$A3 == a3s)

    c1 <- adhere1
    c2 <- c1 * adhere2_now
    c3 <- c2 * adhere3_now

    data.frame(
      row_id = seq_len(nrow(wide_df)),
      analysis_id = wide_df$analysis_id,
      regime = regime,
      d1 = as.integer(regime == "d1"),
      Y = wide_df$Y4,
      age_c = wide_df$age_c,
      male = wide_df$male,
      baseline_log10_vl = wide_df$baseline_log10_vl,
      A1s = a1s,
      A2s = a2s,
      A3s = a3s,
      adhere1 = adhere1,
      adhere2_now = adhere2_now,
      adhere3_now = adhere3_now,
      c1 = c1,
      c2 = c2,
      c3 = c3,
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(build_one("d1"), build_one("d2"))
}

## Fit the numerator of the stabilized artificial-censoring weights.
fit_stabilization_numerator_v7c <- function(
    wide_df,
    ## Purely numerical safeguard (avoid literal 0/1 in floating point), NOT
    ## a substantive positivity bound -- the denominator already carries the
    ## deliberate 0.01/0.99 positivity-management bound (DEFAULT_G_BOUNDS).
    ## Applying that same substantive bound to the numerator would add a
    ## second, redundant positivity intervention specifically on the
    ## stabilization component, which is supposed to only reduce variance,
    ## not manage positivity itself.
    bounds = c(1e-6, 1 - 1e-6)
) {
  if (length(bounds) != 2L || bounds[1] <= 0 || bounds[2] >= 1 ||
      bounds[1] >= bounds[2]) {
    stop("bounds must satisfy 0 < lower < upper < 1")
  }

  clones <- build_regime_clones_v7c(wide_df)

  numerator_x <- function(df) {
    data.frame(
      d1 = df$d1,
      age_c = df$age_c,
      male = df$male,
      baseline_log10_vl = df$baseline_log10_vl
    )
  }

  bound_prob <- function(p) pmin(pmax(p, bounds[1]), bounds[2])

  fit1 <- fit_binary_nuisance_v7(
    Y = clones$adhere1,
    X = numerator_x(clones),
    use_sl = FALSE
  )
  n1 <- bound_prob(
    predict_binary_nuisance_v7(fit1, numerator_x(clones))
  )

  ## Risk sets condition ONLY on remaining uncensored through the preceding
  ## visit (c1==1L / c2==1L), NOT on A1s/A2s (whether the regime has already
  ## prescribed switching). Cain et al.'s CCW numerator must depend only on
  ## regime and baseline covariates -- A1s/A2s are functions of the clone's
  ## own evolving VL trajectory (time-varying, regime-evolving information),
  ## so restricting the numerator's ESTIMATION SAMPLE on them (even without
  ## using them as an explicit regression covariate) violates that
  ## specification. Already-switched, still-adherent clones are now included
  ## in the same fit -- since treatment is absorbing, adhere2_now/adhere3_now
  ## should be deterministically 1 for them, so this adds well-behaved,
  ## easily-fit observations rather than corrupting the model.
  risk2 <- clones$c1 == 1L
  n2 <- rep(1, nrow(clones))
  if (any(risk2)) {
    fit2 <- fit_binary_nuisance_v7(
      Y = clones$adhere2_now[risk2],
      X = numerator_x(clones[risk2, , drop = FALSE]),
      use_sl = FALSE
    )
    n2[risk2] <- bound_prob(
      predict_binary_nuisance_v7(
        fit2,
        numerator_x(clones[risk2, , drop = FALSE])
      )
    )
  } else {
    fit2 <- NULL
  }

  risk3 <- clones$c2 == 1L
  n3 <- rep(1, nrow(clones))
  if (any(risk3)) {
    fit3 <- fit_binary_nuisance_v7(
      Y = clones$adhere3_now[risk3],
      X = numerator_x(clones[risk3, , drop = FALSE]),
      use_sl = FALSE
    )
    n3[risk3] <- bound_prob(
      predict_binary_nuisance_v7(
        fit3,
        numerator_x(clones[risk3, , drop = FALSE])
      )
    )
  } else {
    fit3 <- NULL
  }

  ## Validation: absorbing treatment means already-switched adherent clones
  ## must show adhere2_now/adhere3_now == 1 deterministically. If this fails,
  ## the DGP's absorbing-treatment assumption isn't holding as expected
  ## somewhere upstream, and that's a more serious problem than this function.
  stopifnot(all(
    clones$adhere2_now[clones$c1 == 1L & clones$A1s == 1L] == 1L
  ))
  stopifnot(all(
    clones$adhere3_now[clones$c2 == 1L & clones$A2s == 1L] == 1L
  ))

  clones$num1 <- n1
  clones$num2 <- n2
  clones$num3 <- n3
  clones$numbar3 <- n1 * n2 * n3

  list(
    clones = clones,
    fits = list(fit1 = fit1, fit2 = fit2, fit3 = fit3),
    bounds = bounds
  )
}

fit_standardized_weighted_msm_v7c <- function(msm, wide_df) {
  msm_fit <- msm[
    is.finite(msm$weight) & msm$weight > 0 & is.finite(msm$Y),
    ,
    drop = FALSE
  ]

  if (nrow(msm_fit) == 0L || length(unique(msm_fit$d1)) < 2L) {
    return(safe_iptw_stab_result_v7c(FALSE))
  }

  des <- survey::svydesign(
    ids = ~analysis_id,
    weights = ~weight,
    data = msm_fit
  )

  ## Modified Poisson (quasipoisson, log link) rather than log-binomial:
  ## log-binomial's constrained parameter space (linear predictor must imply
  ## probabilities in [0,1]) becomes numerically fragile once continuous
  ## covariates and estimated weights are both in play -- modified Poisson
  ## (Zou 2004) is the standard epidemiological fix for exactly this, and
  ## svyglm's design-based variance is already robust/sandwich-type, so it
  ## correctly handles the "wrong" Poisson variance function underneath.
  ##
  ## Regime-by-baseline interactions: since the numerator already makes
  ## weights depend on age/male/baseline_log10_vl, allowing the regime
  ## effect to vary by those same covariates avoids assuming a constant risk
  ## ratio across baseline strata that the weighting itself doesn't assume.
  fit_formula <- Y ~ d1 * (age_c + male + baseline_log10_vl)
  fit <- tryCatch(
    survey::svyglm(
      formula = fit_formula,
      design = des,
      family = quasipoisson(link = "log")
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(safe_iptw_stab_result_v7c(FALSE))

  beta <- stats::coef(fit)
  vc <- stats::vcov(fit)
  if (any(!is.finite(beta)) || any(!is.finite(vc))) {
    return(safe_iptw_stab_result_v7c(FALSE))
  }

  rhs_formula <- ~ d1 * (age_c + male + baseline_log10_vl)
  new1 <- data.frame(
    d1 = 1,
    age_c = wide_df$age_c,
    male = wide_df$male,
    baseline_log10_vl = wide_df$baseline_log10_vl
  )
  new0 <- new1
  new0$d1 <- 0

  X1 <- stats::model.matrix(rhs_formula, data = new1)
  X0 <- stats::model.matrix(rhs_formula, data = new0)
  X1 <- X1[, names(beta), drop = FALSE]
  X0 <- X0[, names(beta), drop = FALSE]

  mu1_i <- exp(drop(X1 %*% beta))
  mu0_i <- exp(drop(X0 %*% beta))

  ## Diagnostic: modified Poisson doesn't structurally bound predictions to
  ## [0,1] the way log-binomial does. Unlikely to bite here given regime
  ## risks are ~0.39-0.57, but worth tracking regime-specifically rather
  ## than only a single combined maximum.
  max_predicted_risk_d1 <- suppressWarnings(max(mu1_i, na.rm = TRUE))
  max_predicted_risk_d2 <- suppressWarnings(max(mu0_i, na.rm = TRUE))
  pct_predicted_risk_gt1_d1 <- 100 * mean(mu1_i > 1, na.rm = TRUE)
  pct_predicted_risk_gt1_d2 <- 100 * mean(mu0_i > 1, na.rm = TRUE)

  EY_d1 <- mean(mu1_i)
  EY_d2 <- mean(mu0_i)
  if (!is.finite(EY_d1) || !is.finite(EY_d2) || EY_d1 <= 0 || EY_d2 <= 0) {
    return(safe_iptw_stab_result_v7c(FALSE))
  }

  rr <- EY_d1 / EY_d2
  grad1 <- colMeans(X1 * mu1_i)
  grad0 <- colMeans(X0 * mu0_i)
  grad_log_rr <- grad1 / EY_d1 - grad0 / EY_d2

  var_log_rr <- as.numeric(
    t(grad_log_rr) %*% vc %*% grad_log_rr
  )
  if (!is.finite(var_log_rr)) {
    return(safe_iptw_stab_result_v7c(FALSE))
  }
  se_log <- sqrt(max(var_log_rr, 0))

  c(
    est = rr,
    se = rr * se_log,
    lcl = exp(log(rr) - 1.96 * se_log),
    ucl = exp(log(rr) + 1.96 * se_log),
    EY_d1 = EY_d1,
    EY_d2 = EY_d2,
    max_predicted_risk_d1 = max_predicted_risk_d1,
    max_predicted_risk_d2 = max_predicted_risk_d2,
    pct_predicted_risk_gt1_d1 = pct_predicted_risk_gt1_d1,
    pct_predicted_risk_gt1_d2 = pct_predicted_risk_gt1_d2,
    se_log_rr = se_log
  )
}

build_stabilized_msm_v7c <- function(
    wide_df,
    g1,
    g2,
    g3,
    numerator_fit
) {
  clones <- numerator_fit$clones

  build_one <- function(regime) {
    rp <- regime_probs_v7(wide_df, g1, g2, g3, regime)
    clone_r <- clones[clones$regime == regime, , drop = FALSE]
    clone_r <- clone_r[order(clone_r$row_id), , drop = FALSE]

    I <- clone_r$c3
    denom <- rp$gbar3
    weight <- I * clone_r$numbar3 / denom

    data.frame(
      analysis_id = clone_r$analysis_id,
      d1 = clone_r$d1,
      Y = clone_r$Y,
      age_c = clone_r$age_c,
      male = clone_r$male,
      baseline_log10_vl = clone_r$baseline_log10_vl,
      regime = clone_r$regime,
      weight = weight,
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(build_one("d1"), build_one("d2"))
}

iptw_msm_stabilized_v7c <- function(
    wide_df,
    g1,
    g2,
    g3,
    numerator_fit = fit_stabilization_numerator_v7c(wide_df)
) {
  msm <- build_stabilized_msm_v7c(
    wide_df, g1, g2, g3, numerator_fit
  )

  fit_result <- fit_standardized_weighted_msm_v7c(msm, wide_df)
  if (!is.finite(fit_result["est"])) {
    ans <- safe_iptw_stab_result_v7c(FALSE)
  } else {
    ans <- safe_iptw_stab_result_v7c(FALSE)
    ans[names(fit_result)] <- fit_result
  }

  w1 <- msm$weight[msm$d1 == 1]
  w2 <- msm$weight[msm$d1 == 0]
  ans["ess_d1"] <- weight_ess_v7c(w1)
  ans["ess_d2"] <- weight_ess_v7c(w2)
  ans["w01_d1"] <- weight_quantile_v7c(w1, 0.01)
  ans["w01_d2"] <- weight_quantile_v7c(w2, 0.01)
  ans["w99_d1"] <- weight_quantile_v7c(w1, 0.99)
  ans["w99_d2"] <- weight_quantile_v7c(w2, 0.99)
  extra1 <- extra_weight_diagnostics_v7(w1)
  extra2 <- extra_weight_diagnostics_v7(w2)
  ans["mean_positive_weight_d1"] <- extra1["mean_positive"]; ans["mean_positive_weight_d2"] <- extra2["mean_positive"]
  ans["sd_positive_weight_d1"] <- extra1["sd_positive"]; ans["sd_positive_weight_d2"] <- extra2["sd_positive"]
  ans["cv_d1"] <- extra1["cv_positive"]; ans["cv_d2"] <- extra2["cv_positive"]
  ans["max_weight_d1"] <- extra1["max_positive"]; ans["max_weight_d2"] <- extra2["max_positive"]
  ans["pct_positive_d1"] <- extra1["pct_positive"]; ans["pct_positive_d2"] <- extra2["pct_positive"]
  ans
}

truncate_positive_weights_v7c <- function(w, lower = 0.01, upper = 0.99) {
  out <- w
  pos <- which(is.finite(w) & w > 0)
  if (length(pos) == 0L) {
    return(list(
      weights = out,
      lower_cutoff = NA_real_,
      upper_cutoff = NA_real_,
      pct_lower = NA_real_,
      pct_upper = NA_real_,
      pct_any = NA_real_
    ))
  }

  lo <- as.numeric(stats::quantile(w[pos], lower, na.rm = TRUE, names = FALSE))
  hi <- as.numeric(stats::quantile(w[pos], upper, na.rm = TRUE, names = FALSE))
  lower_hit <- w[pos] < lo
  upper_hit <- w[pos] > hi
  out[pos] <- pmin(pmax(w[pos], lo), hi)

  list(
    weights = out,
    lower_cutoff = lo,
    upper_cutoff = hi,
    pct_lower = 100 * mean(lower_hit),
    pct_upper = 100 * mean(upper_hit),
    pct_any = 100 * mean(lower_hit | upper_hit)
  )
}

iptw_msm_stabilized_trunc_v7c <- function(
    wide_df,
    g1,
    g2,
    g3,
    numerator_fit = fit_stabilization_numerator_v7c(wide_df),
    lower_quantile = 0.01,
    upper_quantile = 0.99
) {
  if (lower_quantile < 0 || upper_quantile > 1 ||
      lower_quantile >= upper_quantile) {
    stop("Require 0 <= lower_quantile < upper_quantile <= 1")
  }

  msm <- build_stabilized_msm_v7c(
    wide_df, g1, g2, g3, numerator_fit
  )

  raw_w1 <- msm$weight[msm$d1 == 1]
  raw_w2 <- msm$weight[msm$d1 == 0]

  t1 <- truncate_positive_weights_v7c(
    raw_w1, lower_quantile, upper_quantile
  )
  t2 <- truncate_positive_weights_v7c(
    raw_w2, lower_quantile, upper_quantile
  )

  msm$weight[msm$d1 == 1] <- t1$weights
  msm$weight[msm$d1 == 0] <- t2$weights

  fit_result <- fit_standardized_weighted_msm_v7c(msm, wide_df)
  ans <- safe_iptw_stab_result_v7c(TRUE)
  if (is.finite(fit_result["est"])) {
    ans[names(fit_result)] <- fit_result
  }

  w1 <- msm$weight[msm$d1 == 1]
  w2 <- msm$weight[msm$d1 == 0]

  ans["ess_d1"] <- weight_ess_v7c(w1)
  ans["ess_d2"] <- weight_ess_v7c(w2)
  ans["w01_d1"] <- weight_quantile_v7c(w1, 0.01)
  ans["w01_d2"] <- weight_quantile_v7c(w2, 0.01)
  ans["w99_d1"] <- weight_quantile_v7c(w1, 0.99)
  ans["w99_d2"] <- weight_quantile_v7c(w2, 0.99)
  extra1 <- extra_weight_diagnostics_v7(w1)
  extra2 <- extra_weight_diagnostics_v7(w2)
  ans["mean_positive_weight_d1"] <- extra1["mean_positive"]; ans["mean_positive_weight_d2"] <- extra2["mean_positive"]
  ans["sd_positive_weight_d1"] <- extra1["sd_positive"]; ans["sd_positive_weight_d2"] <- extra2["sd_positive"]
  ans["cv_d1"] <- extra1["cv_positive"]; ans["cv_d2"] <- extra2["cv_positive"]
  ans["max_weight_d1"] <- extra1["max_positive"]; ans["max_weight_d2"] <- extra2["max_positive"]
  ans["pct_positive_d1"] <- extra1["pct_positive"]; ans["pct_positive_d2"] <- extra2["pct_positive"]

  ans["pretrunc_ess_d1"] <- weight_ess_v7c(raw_w1)
  ans["pretrunc_ess_d2"] <- weight_ess_v7c(raw_w2)
  ans["pretrunc_w01_d1"] <- weight_quantile_v7c(raw_w1, 0.01)
  ans["pretrunc_w01_d2"] <- weight_quantile_v7c(raw_w2, 0.01)
  ans["pretrunc_w99_d1"] <- weight_quantile_v7c(raw_w1, 0.99)
  ans["pretrunc_w99_d2"] <- weight_quantile_v7c(raw_w2, 0.99)
  extra_raw1 <- extra_weight_diagnostics_v7(raw_w1)
  extra_raw2 <- extra_weight_diagnostics_v7(raw_w2)
  ans["pretrunc_mean_positive_weight_d1"] <- extra_raw1["mean_positive"]; ans["pretrunc_mean_positive_weight_d2"] <- extra_raw2["mean_positive"]
  ans["pretrunc_sd_positive_weight_d1"] <- extra_raw1["sd_positive"]; ans["pretrunc_sd_positive_weight_d2"] <- extra_raw2["sd_positive"]
  ans["pretrunc_cv_d1"] <- extra_raw1["cv_positive"]; ans["pretrunc_cv_d2"] <- extra_raw2["cv_positive"]
  ans["pretrunc_max_weight_d1"] <- extra_raw1["max_positive"]; ans["pretrunc_max_weight_d2"] <- extra_raw2["max_positive"]
  ## pct_positive not repeated pre/post -- winsorizing doesn't change which
  ## weights are structurally zero, already saved once above.

  ans["lower_cutoff_d1"] <- t1$lower_cutoff
  ans["lower_cutoff_d2"] <- t2$lower_cutoff
  ans["upper_cutoff_d1"] <- t1$upper_cutoff
  ans["upper_cutoff_d2"] <- t2$upper_cutoff
  ans["pct_lower_truncated_d1"] <- t1$pct_lower
  ans["pct_lower_truncated_d2"] <- t2$pct_lower
  ans["pct_upper_truncated_d1"] <- t1$pct_upper
  ans["pct_upper_truncated_d2"] <- t2$pct_upper
  ans["pct_any_truncated_d1"] <- t1$pct_any
  ans["pct_any_truncated_d2"] <- t2$pct_any
  ans
}


one_sim_all_estimators_v7 <- function(
    n,
    seed,
    scenario,
    dtr1_threshold = 50,
    dtr2_threshold = 1000,
    g_bounds = DEFAULT_G_BOUNDS,
    n_quad = DEFAULT_N_QUAD,
    sl_lib = stable_sl_library_v7(),
    gcomp_B = 0L
) {
  if (!exists("simulate_hiv_violation_v6", mode = "function")) {
    stop("Source simulate_hiv_all_v6.R before running the estimators.")
  }

  dat <- simulate_hiv_violation_v6(
    n = n,
    seed = seed,
    scenario = scenario,
    dtr1_threshold = dtr1_threshold,
    dtr2_threshold = dtr2_threshold
  )

  wide <- prepare_wide_estimation_v7(
    dat,
    dtr1_threshold = dtr1_threshold,
    dtr2_threshold = dtr2_threshold
  )

  quadrature <- normal_quadrature(n_quad)

  ps_par <- fit_ps_v7(
    wide,
    use_sl = FALSE,
    g_bounds = g_bounds,
    seed = seed + 1000L
  )

  ps_sl <- suppressWarnings(
    fit_ps_v7(
      wide,
      use_sl = TRUE,
      SL.library = sl_lib,
      g_bounds = g_bounds,
      seed = seed + 1000L
    )
  )

  q_par <- fit_oracle_transition_v7(
    wide,
    use_sl = FALSE,
    seed = seed
  )

  q_sl <- suppressWarnings(
    fit_oracle_transition_v7(
      wide,
      use_sl = TRUE,
      SL.library = sl_lib,
      seed = seed
    )
  )

  ## Computed once per PS branch and reused by both iptw_msm_stabilized_v7c()
  ## and iptw_msm_stabilized_trunc_v7c() below -- avoids refitting the same
  ## three numerator models twice per branch.
  numerator_fit_par <- tryCatch(
    fit_stabilization_numerator_v7c(wide),
    error = function(e) NULL
  )
  numerator_fit_sl <- numerator_fit_par  ## numerator models are parametric by
  ## design (fit_binary_nuisance_v7(..., use_sl=FALSE) is hardcoded inside
  ## fit_stabilization_numerator_v7c) -- identical for both PS branches, so
  ## reused rather than recomputed.

  safe <- function(nm, expr) {
    out <- tryCatch(
      expr,
      error = function(e) {
        message(sprintf("[safe] %s failed: %s", nm, conditionMessage(e)))
        safe_ratio_result()
      }
    )
    names(out) <- paste0(nm, "_", names(out))
    out
  }

  c(
    safe(
      "IPTW_PAR",
      iptw_msm_v7(wide, ps_par$g1, ps_par$g2, ps_par$g3)
    ),
    safe(
      "IPTW_SL",
      iptw_msm_v7(wide, ps_sl$g1, ps_sl$g2, ps_sl$g3)
    ),
    safe(
      "IPTW_TRUNC01_99_PAR",
      iptw_msm_trunc01_99_v7(wide, ps_par$g1, ps_par$g2, ps_par$g3)
    ),
    safe(
      "IPTW_TRUNC01_99_SL",
      iptw_msm_trunc01_99_v7(wide, ps_sl$g1, ps_sl$g2, ps_sl$g3)
    ),
    safe(
      "IPTW_STAB_PAR",
      iptw_msm_stabilized_v7c(wide, ps_par$g1, ps_par$g2, ps_par$g3,
                                numerator_fit = numerator_fit_par)
    ),
    safe(
      "IPTW_STAB_SL",
      iptw_msm_stabilized_v7c(wide, ps_sl$g1, ps_sl$g2, ps_sl$g3,
                                numerator_fit = numerator_fit_sl)
    ),
    safe(
      "IPTW_STAB_TRUNC01_99_PAR",
      iptw_msm_stabilized_trunc_v7c(wide, ps_par$g1, ps_par$g2, ps_par$g3,
                                      numerator_fit = numerator_fit_par)
    ),
    safe(
      "IPTW_STAB_TRUNC01_99_SL",
      iptw_msm_stabilized_trunc_v7c(wide, ps_sl$g1, ps_sl$g2, ps_sl$g3,
                                      numerator_fit = numerator_fit_sl)
    ),
    safe(
      "AIPTW_PAR",
      aiptw_rr_v7(
        wide, ps_par$g1, ps_par$g2, ps_par$g3,
        q_par, quadrature
      )
    ),
    safe(
      "AIPTW_SL",
      aiptw_rr_v7(
        wide, ps_sl$g1, ps_sl$g2, ps_sl$g3,
        q_sl, quadrature
      )
    ),
    safe(
      "TMLE_PAR",
      tmle_rr_v7(
        wide, ps_par$g1, ps_par$g2, ps_par$g3,
        q_par, quadrature
      )
    ),
    safe(
      "TMLE_SL",
      tmle_rr_v7(
        wide, ps_sl$g1, ps_sl$g2, ps_sl$g3,
        q_sl, quadrature
      )
    ),
    safe(
      "GCOMP_PAR",
      if (gcomp_B > 0L) {
        bootstrap_gcomp_oracle_v7(
          wide,
          use_sl = FALSE,
          B = gcomp_B,
          seed = seed,
          n_quad = n_quad
        )
      } else {
        gcomp_rr_v7(wide, q_par, quadrature)
      }
    ),
    safe(
      "GCOMP_SL",
      if (gcomp_B > 0L) {
        bootstrap_gcomp_oracle_v7(
          wide,
          use_sl = TRUE,
          SL.library = sl_lib,
          B = gcomp_B,
          seed = seed,
          n_quad = n_quad
        )
      } else {
        gcomp_rr_v7(wide, q_sl, quadrature)
      }
    )
  )
}

## Backward-compatible alias used by older HPC runner scripts.
one_sim_all_estimators_v5 <- one_sim_all_estimators_v7

## ============================================================================
## 12. GRID RUNNER AND SUMMARY
## ============================================================================

run_simulation_grid_v7 <- function(
    n_values = c(500, 1500, 5000),
    scenarios = c("benchmark", "mild", "moderate", "severe"),
    n_sim = 50L,
    ...
) {
  grid <- expand.grid(
    n = n_values,
    scenario = scenarios,
    sim = seq_len(n_sim),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    arrange(n, scenario, sim) %>%
    mutate(seed = as.integer(mapply(make_seed_v7, n = n, scenario = scenario, sim = sim)))

  out <- lapply(seq_len(nrow(grid)), function(i) {
    row <- grid[i, ]
    est <- tryCatch(
      one_sim_all_estimators_v7(
        n = row$n,
        seed = row$seed,
        scenario = row$scenario,
        ...
      ),
      error = function(e) {
        message(sprintf("Simulation n=%d scenario=%s sim=%d (seed=%d) failed entirely: %s",
                         row$n, row$scenario, row$sim, row$seed, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(est)) return(cbind(row, failed = TRUE))
    cbind(row, failed = FALSE, as.data.frame(as.list(est), check.names = FALSE))
  })

  bind_rows(out)
}

run_simulation_grid_v5 <- run_simulation_grid_v7

## Original summarise_results_grid_v7() (8-estimator list) is superseded
## by summarise_results_grid_v7c() (12-estimator list, includes the
## truncated/stabilized IPTW variants added above) -- aliased under the
## ORIGINAL name so existing HPC scripts (combine_hpc_results_v7.R etc.)
## that call summarise_results_grid_v7() by name get the complete,
## up-to-date estimator list with no changes needed on their end.
## Summary used after the v7c add-on is merged with the existing main-grid data.
summarise_results_grid_v7c <- function(
    sim_res,
    true_rr = TRUE_RR,
    gcomp_ci_expected = FALSE
) {
  estimator_names <- c(
    "IPTW_PAR", "IPTW_SL",
    "IPTW_TRUNC01_99_PAR", "IPTW_TRUNC01_99_SL",
    "IPTW_STAB_PAR", "IPTW_STAB_SL",
    "IPTW_STAB_TRUNC01_99_PAR", "IPTW_STAB_TRUNC01_99_SL",
    "AIPTW_PAR", "AIPTW_SL",
    "TMLE_PAR", "TMLE_SL",
    "GCOMP_PAR", "GCOMP_SL"
  )

  dplyr::bind_rows(lapply(estimator_names, function(est) {
    ec <- paste0(est, "_est")
    sc <- paste0(est, "_se")
    lc <- paste0(est, "_lcl")
    uc <- paste0(est, "_ucl")
    e1 <- paste0(est, "_EY_d1")
    e2 <- paste0(est, "_EY_d2")
    es1 <- paste0(est, "_ess_d1")
    es2 <- paste0(est, "_ess_d2")
    w1 <- paste0(est, "_w99_d1")
    w2 <- paste0(est, "_w99_d2")

    if (!ec %in% names(sim_res)) return(NULL)

    ci_expected <- !grepl("^GCOMP_", est) || isTRUE(gcomp_ci_expected)

    sim_res |>
      dplyr::group_by(n, scenario) |>
      dplyr::summarise(
        estimator = est,
        mean_rr = mean(.data[[ec]], na.rm = TRUE),
        bias = mean(.data[[ec]], na.rm = TRUE) - true_rr,
        rel_bias_pct = 100 *
          (mean(.data[[ec]], na.rm = TRUE) - true_rr) / true_rr,
        empirical_sd = stats::sd(.data[[ec]], na.rm = TRUE),
        variance = stats::var(.data[[ec]], na.rm = TRUE),
        mse = mean((.data[[ec]] - true_rr)^2, na.rm = TRUE),
        mean_se = if (sc %in% names(sim_res)) {
          mean(.data[[sc]], na.rm = TRUE)
        } else {
          NA_real_
        },
        se_to_empirical_sd = if (sc %in% names(sim_res)) {
          mean(.data[[sc]], na.rm = TRUE) /
            stats::sd(.data[[ec]], na.rm = TRUE)
        } else {
          NA_real_
        },
        mean_ci_width = if (all(c(lc, uc) %in% names(sim_res))) {
          mean(.data[[uc]] - .data[[lc]], na.rm = TRUE)
        } else {
          NA_real_
        },
        n_valid = sum(is.finite(.data[[ec]])),
        n_ci_valid = if (!ci_expected) {
          NA_integer_
        } else if (all(c(lc, uc) %in% names(sim_res))) {
          sum(is.finite(.data[[lc]]) & is.finite(.data[[uc]]))
        } else {
          0L
        },
        coverage_pct = if (!ci_expected) {
          NA_real_
        } else if (all(c(lc, uc) %in% names(sim_res))) {
          ok <- is.finite(.data[[lc]]) & is.finite(.data[[uc]])
          if (any(ok)) {
            mean(
              .data[[lc]][ok] <= true_rr &
                .data[[uc]][ok] >= true_rr
            ) * 100
          } else {
            NA_real_
          }
        } else {
          NA_real_
        },
        ci_failure_pct = if (!ci_expected) {
          NA_real_
        } else if (all(c(lc, uc) %in% names(sim_res))) {
          mean(!(
            is.finite(.data[[lc]]) & is.finite(.data[[uc]])
          )) * 100
        } else {
          100
        },
        mean_EY_d1 = if (e1 %in% names(sim_res)) {
          mean(.data[[e1]], na.rm = TRUE)
        } else {
          NA_real_
        },
        mean_EY_d2 = if (e2 %in% names(sim_res)) {
          mean(.data[[e2]], na.rm = TRUE)
        } else {
          NA_real_
        },
        ess_d1 = if (es1 %in% names(sim_res)) {
          stats::median(.data[[es1]], na.rm = TRUE)
        } else {
          NA_real_
        },
        ess_d2 = if (es2 %in% names(sim_res)) {
          stats::median(.data[[es2]], na.rm = TRUE)
        } else {
          NA_real_
        },
        w99_d1 = if (w1 %in% names(sim_res)) {
          stats::median(.data[[w1]], na.rm = TRUE)
        } else {
          NA_real_
        },
        w99_d2 = if (w2 %in% names(sim_res)) {
          stats::median(.data[[w2]], na.rm = TRUE)
        } else {
          NA_real_
        },
        .groups = "drop"
      )
  })) |>
    dplyr::mutate(
      scenario = factor(
        scenario,
        levels = c("benchmark", "mild", "moderate", "severe")
      ),
      estimator = factor(estimator, levels = estimator_names)
    ) |>
    dplyr::arrange(n, scenario, estimator)
}

summarise_results_grid_v7 <- summarise_results_grid_v7c
summarise_results_grid_v5 <- summarise_results_grid_v7

## ============================================================================
## 13. REQUIRED VALIDATION CHECKS BEFORE THE FULL HPC RUN
## ============================================================================

validate_transition_model_v7 <- function(
    n = 50000,
    seed = 999,
    scenario = "benchmark"
) {
  dat <- simulate_hiv_violation_v6(
    n = n,
    seed = seed,
    scenario = scenario
  )
  wide <- prepare_wide_estimation_v7(dat)
  fit <- fit_oracle_transition_v7(wide, use_sl = FALSE)

  list(
    coefficients = coef(fit$fit),
    expected = c(
      `(Intercept)` = 0,
      logvl_current = 1,
      slope_A = 3,
      A_slope_A = -3,
      A_slope_B = 3
    ),
    sigma_hat = fit$sigma,
    sigma_expected = 0.08
  )
}

validate_g_determinism_v7 <- function(
    n = 5000,
    seed = 999,
    scenario = "benchmark"
) {
  dat <- simulate_hiv_violation_v6(
    n = n,
    seed = seed,
    scenario = scenario
  )
  wide <- prepare_wide_estimation_v7(dat)
  ps <- fit_ps_v7(wide, use_sl = FALSE)

  c(
    all_g2_post_A1_equal_1 = all(ps$g2[wide$A1 == 1] == 1),
    all_g3_post_switch_equal_1 = all(
      ps$g3[wide$A1 == 1 | wide$A2 == 1] == 1
    ),
    min_stochastic_g = min(
      c(
        ps$g1,
        ps$g2[wide$A1 == 0],
        ps$g3[wide$A1 == 0 & wide$A2 == 0]
      )
    ),
    max_stochastic_g = max(
      c(
        ps$g1,
        ps$g2[wide$A1 == 0],
        ps$g3[wide$A1 == 0 & wide$A2 == 0]
      )
    )
  )
}

check_quadrature_stability_v7 <- function(
    n = 5000,
    seed = 999,
    scenario = "benchmark",
    n_quad_small = 15L,
    n_quad_large = 25L
) {
  dat <- simulate_hiv_violation_v6(
    n = n,
    seed = seed,
    scenario = scenario
  )
  wide <- prepare_wide_estimation_v7(dat)
  fit <- fit_oracle_transition_v7(wide, use_sl = FALSE)

  small <- gcomp_rr_v7(
    wide,
    fit,
    normal_quadrature(n_quad_small)
  )
  large <- gcomp_rr_v7(
    wide,
    fit,
    normal_quadrature(n_quad_large)
  )

  c(
    rr_small = unname(small["est"]),
    rr_large = unname(large["est"]),
    absolute_difference = abs(
      unname(small["est"] - large["est"])
    )
  )
}

check_primary_estimators_v7 <- function(
    n = 5000,
    n_sim = 20L,
    scenario = "benchmark",
    n_quad = DEFAULT_N_QUAD
) {
  res <- lapply(seq_len(n_sim), function(j) {
    one <- one_sim_all_estimators_v7(
      n = n,
      seed = make_seed_v7(n = n, scenario = scenario, sim = j),
      scenario = scenario,
      n_quad = n_quad,
      gcomp_B = 0L
    )

    data.frame(
      sim = j,
      IPTW_PAR = unname(one["IPTW_PAR_est"]),
      AIPTW_PAR = unname(one["AIPTW_PAR_est"]),
      TMLE_PAR = unname(one["TMLE_PAR_est"]),
      GCOMP_PAR = unname(one["GCOMP_PAR_est"])
    )
  }) %>% bind_rows()

  summary <- bind_rows(lapply(names(res)[-1], function(nm) {
    x <- res[[nm]]
    data.frame(
      estimator = nm,
      n_valid = sum(is.finite(x)),
      mean_rr = mean(x, na.rm = TRUE),
      bias = mean(x, na.rm = TRUE) - TRUE_RR,
      rel_bias_pct = 100 * (mean(x, na.rm = TRUE) - TRUE_RR) / TRUE_RR,
      sd_rr = sd(x, na.rm = TRUE),
      mcse_mean = sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
    )
  }))

  list(results = res, summary = summary)
}

## ============================================================================
## END OF FILE
## ============================================================================
