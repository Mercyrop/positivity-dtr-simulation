## ============================================================================
## 01_data_generating_mechanism.R
##
## Public release of the active data-generating mechanism for the longitudinal
## HIV dynamic-treatment-regime positivity simulation study.
##
## DTR1: switch when viral load >= 50 copies/mL, then remain switched.
## DTR2: switch when viral load > 1000 copies/mL, then remain switched.
## Outcome: viral suppression (viral load < 50 copies/mL) at month 12.
##
## Positivity deterioration is induced only through the observational switching
## process. The viral-load biology and the intervention definitions are held
## fixed across benchmark, mild, moderate, and severe scenarios.
##
## Locked Monte Carlo truth used in the manuscript:
##   E[Y^d1] = 0.5716074
##   E[Y^d2] = 0.3889531
##   RR      = 1.469604947
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

simulate_hiv_violation_v6 <- function(

    ## Sample size and random seed
    n    = 1000,
    seed = 123,

    ## Intervention (same as benchmark)
    regime = NULL,

    ## Positivity scenario
    scenario = c("benchmark", "mild", "moderate", "severe"),

    ## DTR thresholds (copies/mL) -- fixed across scenarios
    dtr1_threshold = 50,
    dtr2_threshold = 1000,

    ## Biology parameters -- identical to benchmark, do not change
    p_stable     = 0.68,
    p_slow_fail  = 0.12,
    p_rapid_fail = 0.20,
    baseline_log10_mean = log10(30),
    baseline_log10_sd   = 0.08,
    slope_A_stable_mean = 0.010,  slope_A_slow_mean = 0.090,  slope_A_rapid_mean = 0.320,
    slope_A_stable_sd   = 0.020,  slope_A_slow_sd   = 0.040,  slope_A_rapid_sd   = 0.060,
    age_coef  = 0.006,
    male_coef = 0.012,
    slope_B_stable_mean = 0.000,  slope_B_slow_mean = -0.050, slope_B_rapid_mean = -0.020,
    slope_B_stable_sd   = 0.020,  slope_B_slow_sd   =  0.040, slope_B_rapid_sd   =  0.080,
    slope_B_rho    = 0.25,
    sigma_interval = 0.080

) {

  scenario <- match.arg(scenario)
  set.seed(seed)
  sigmoid <- function(x) 1 / (1 + exp(-x))

  ## ---- Scenario-specific switch-barrier parameters ----
  ##
  ## The biology is unchanged. Only the observed treatment mechanism changes.
  ##
  ## confirm_prob:
  ##   Probability that confirmation evidence is available among patients
  ##   with current detectable-but-not-failing VL and prior detectable VL.
  ##
  ## no_confirm_penalty:
  ##   Additional penalty against switching when VL is detectable-but-not-failing
  ##   but rebound is not confirmed.
  ##
  ## The failure coefficient is kept high and stable so that support for the
  ## delayed-switch regime is preserved in the primary design.
  ##
  ## Final calibrated switch-barrier settings. Biological parameters and the
  ## coefficients in the switching model are held fixed; only confirmation
  ## availability and the penalty for switching without confirmation vary.
  ## This creates a monotone loss of support for the early-switch regime while
  ## preserving the failure-triggered delayed-switch pathway.

  pars <- switch(scenario,

    benchmark = list(
      intercept = -2.20,
      det_only  =  2.40,
      persist   =  0.10,
      deep      =  0.20,
      confirm   =  0.30,
      fail      =  4.00,
      confirm_prob       = 0.95,
      no_confirm_penalty = 0.20
    ),

    mild = list(
      intercept = -2.20,
      det_only  =  2.40,
      persist   =  0.10,
      deep      =  0.20,
      confirm   =  0.30,
      fail      =  4.00,
      confirm_prob       = 0.85,
      no_confirm_penalty = 1.00
    ),

    moderate = list(
      intercept = -2.20,
      det_only  =  2.40,
      persist   =  0.10,
      deep      =  0.20,
      confirm   =  0.30,
      fail      =  4.00,
      confirm_prob       = 0.55,
      no_confirm_penalty = 1.90
    ),

    severe = list(
      intercept = -2.20,
      det_only  =  2.40,
      persist   =  0.10,
      deep      =  0.20,
      confirm   =  0.30,
      fail      =  4.00,
      confirm_prob       = 0.25,
      no_confirm_penalty = 3.20
    )
  )

  ## ---- Biology (identical to benchmark) ----
  visits <- c(0, 3, 6, 9, 12)
  K      <- length(visits) - 1
  id     <- seq_len(n)
  age    <- pmax(round(rnorm(n, 38, 12)), 15)
  female <- rbinom(n, 1, 0.65)
  male   <- 1 - female
  age_c  <- (age - 38) / 12

  traj_type <- factor(
    sample(c("stable", "slow_fail", "rapid_fail"), n, replace = TRUE,
           prob = c(p_stable, p_slow_fail, p_rapid_fail)),
    levels = c("stable", "slow_fail", "rapid_fail"))

  sA_mean <- ifelse(traj_type == "stable",    slope_A_stable_mean,
             ifelse(traj_type == "slow_fail", slope_A_slow_mean,
                                              slope_A_rapid_mean))
  sA_sd   <- ifelse(traj_type == "stable",    slope_A_stable_sd,
             ifelse(traj_type == "slow_fail", slope_A_slow_sd,
                                              slope_A_rapid_sd))
  slope_A <- rnorm(n, sA_mean + age_coef * age_c + male_coef * male, sA_sd)

  sB_mean <- ifelse(traj_type == "stable",    slope_B_stable_mean,
             ifelse(traj_type == "slow_fail", slope_B_slow_mean,
                                              slope_B_rapid_mean))
  sB_sd   <- ifelse(traj_type == "stable",    slope_B_stable_sd,
             ifelse(traj_type == "slow_fail", slope_B_slow_sd,
                                              slope_B_rapid_sd))
  sAc     <- ave(slope_A, traj_type, FUN = function(x) x - mean(x))
  slope_B <- sB_mean + slope_B_rho * sAc + rnorm(n, 0, sB_sd)

  baseline_log10_vl <- pmin(
    rnorm(n, baseline_log10_mean, baseline_log10_sd),
    log10(dtr1_threshold - 1))

  lv_prev  <- 10^baseline_log10_vl
  log_prev <- baseline_log10_vl
  A_prev   <- rep(0L, n)
  out      <- vector("list", K + 1)

  out[[1]] <- data.frame(
    id     = id,  k = 0,  month = 0,
    age    = age,  female = female,  male = male,  age_c = age_c,
    traj_type = traj_type,
    slope_A = slope_A,  slope_B = slope_B,  baseline_log10_vl = baseline_log10_vl,
    A_prev = 0L,  A_next = 0L,
    log10_vl            = baseline_log10_vl,
    viral_load          = lv_prev,
    viral_load_reported = pmax(round(lv_prev), dtr1_threshold - 1),
    det = 0L,  fail = 0L,  det_only = 0L,
    lag_det = 0L,  lag_fail = 0L,  lag_det_only = 0L,
    persist_det_only      = 0L,
    deep_persist_det_only = 0L,
    confirmation_available = 0L,
    confirmed_rebound      = 0L,
    unconfirmed_det_only   = 0L,
    pA = NA_real_)

  for (k in 1:K) {

    month_k   <- visits[k + 1]
    eff_slope <- ifelse(A_prev == 1L, slope_B, slope_A)
    log10_vl  <- log_prev + 3 * eff_slope + rnorm(n, 0, sigma_interval)
    latent_vl <- 10^log10_vl

    det      <- as.integer(latent_vl >= dtr1_threshold)
    fail     <- as.integer(latent_vl >  dtr2_threshold)
    det_only <- as.integer(det == 1L & fail == 0L)

    lag_det      <- as.integer(lv_prev >= dtr1_threshold)
    lag_fail     <- as.integer(lv_prev >  dtr2_threshold)
    lag_det_only <- as.integer(lag_det == 1L & lag_fail == 0L)

    persist_det_only <- as.integer(det_only == 1L & lag_det_only == 1L)

    deep_persist_det_only <- if (k == 1) {
      rep(0L, n)
    } else {
      as.integer(det_only == 1L & lag_det_only == 1L & out[[k - 1]]$det_only == 1L)
    }

    ## ---- Switch-barrier variables ----
    ##
    ## Confirmation is relevant only for current detectable-but-not-failing VL.
    ## It requires prior detectable VL history plus availability of confirmation.
    ## We initialize here so that forced regimes and month 12 remain well-defined.

    confirmation_available <- rep(0L, n)
    confirmed_rebound      <- rep(0L, n)
    unconfirmed_det_only   <- rep(0L, n)

    if (k <= 3) {

      if (is.null(regime)) {

        ## OBSERVATIONAL: scenario-specific switch-barrier law

        ## Confirmation eligibility:
        ## current VL is detectable-but-not-failing AND there was prior detectable VL.
        eligible_confirmation <- as.integer(det_only == 1L & lag_det == 1L)

        ## Confirmation availability becomes less common as scenarios worsen.
        confirm_draw <- rbinom(n, 1, pars$confirm_prob)

        confirmation_available <- as.integer(
          eligible_confirmation == 1L & confirm_draw == 1L
        )

        ## Confirmed rebound is an observed care-process/history variable.
        confirmed_rebound <- as.integer(
          det_only == 1L & confirmation_available == 1L
        )

        ## The switch barrier applies only to detectable-but-not-failing histories
        ## without confirmation. It does NOT penalize failure histories.
        unconfirmed_det_only <- as.integer(
          det_only == 1L & confirmed_rebound == 0L
        )

        logitA <- pars$intercept +
          pars$det_only * det_only +
          pars$fail     * fail +
          pars$persist  * persist_det_only +
          pars$deep     * deep_persist_det_only +
          pars$confirm  * confirmed_rebound -
          pars$no_confirm_penalty * unconfirmed_det_only +
          0.15 * age_c +
          0.10 * male

        pA     <- sigmoid(logitA)
        A_next <- rbinom(n, 1, pA)
        A_next[A_prev == 1L] <- 1L

      } else {

        ## FORCED REGIME: switch according to DTR definition
        pA     <- NA_real_
        target <- if (regime == "d1") det else fail
        A_next <- target
        A_next[A_prev == 1L] <- 1L

      }

    } else {
      ## Month 12: no switching decision, just record outcome
      pA     <- NA_real_
      A_next <- A_prev
    }

    out[[k + 1]] <- data.frame(
      id     = id,  k = k,  month = month_k,
      age    = age,  female = female,  male = male,  age_c = age_c,
      traj_type = traj_type,
      slope_A = slope_A,  slope_B = slope_B,  baseline_log10_vl = baseline_log10_vl,
      A_prev = A_prev,  A_next = A_next,
      log10_vl            = log10_vl,
      viral_load          = latent_vl,
      viral_load_reported = pmax(round(latent_vl), dtr1_threshold - 1),
      det = det,  fail = fail,  det_only = det_only,
      lag_det = lag_det,  lag_fail = lag_fail,  lag_det_only = lag_det_only,
      persist_det_only      = persist_det_only,
      deep_persist_det_only = deep_persist_det_only,
      confirmation_available = confirmation_available,
      confirmed_rebound      = confirmed_rebound,
      unconfirmed_det_only   = unconfirmed_det_only,
      pA = pA)

    A_prev   <- A_next
    log_prev <- log10_vl
    lv_prev  <- latent_vl
  }

  bind_rows(out) %>% arrange(id, k)
}


## ============================================================================
## 2c.  FIRST CALIBRATION CHECK FOR v6
##
##      Band-level median switching probability, broken out by whether
##      det_only rebound was confirmed. Run this before handing v6 to any
##      estimator. What you want to see:
##        - unconfirmed_det_only median pA decreases progressively
##          (benchmark > mild > moderate > severe)
##        - confirmed_det_only stays meaningfully higher than unconfirmed
##          within each scenario
##        - failure stays high and stable across all scenarios
##
## ============================================================================

check_switch_barrier_v6 <- function(
    n         = 20000,
    seed      = 999,
    scenarios = c("benchmark", "mild", "moderate", "severe")) {

  for (sc in scenarios) {
    dat <- simulate_hiv_violation_v6(n = n, seed = seed, scenario = sc)

    dec <- dat |>
      dplyr::filter(month %in% c(3, 6, 9), A_prev == 0)

    cat("\n\nScenario:", sc, "\n")

    print(
      dec |>
        dplyr::mutate(
          band = dplyr::case_when(
            fail == 1 ~ "failure",
            det_only == 1 & confirmed_rebound == 1 ~ "confirmed_det_only",
            det_only == 1 & unconfirmed_det_only == 1 ~ "unconfirmed_det_only",
            TRUE ~ "suppressed"
          )
        ) |>
        dplyr::group_by(band) |>
        dplyr::summarise(
          n = dplyr::n(),
          median_pA = round(median(pA, na.rm = TRUE), 3),
          pct_pA_lt_0_10 = round(100 * mean(pA < 0.10, na.rm = TRUE), 1),
          .groups = "drop"
        )
    )
  }

  invisible(NULL)
}


## ============================================================================
## 2d.  FULL DIAGNOSTIC SUITE FOR v6 (band + path + final level)
##
##      Adapted from verify_violation_scenarios_v5() in simulate_hiv_all_v5.R,
##      fixing two issues found when that function is pointed at v6 data:
##
##      1. HARDCODED DGP CALL. The v5 version calls simulate_hiv_violation_v5()
##         by name. This version takes dgp_fun (default simulate_hiv_violation_v6)
##         so it works for either DGP, or any future variant with the same
##         calling convention.
##
##      2. ABSORBING-TREATMENT PATH-PROBABILITY BUG. The v5 version computes
##         p2 <- ifelse(A2s==1, pA_k2, 1-pA_k2) unconditionally -- it never
##         checks whether the patient had ALREADY switched at a prior visit.
##         pA_k2/pA_k3 are computed for every patient at every visit in the
##         DGP (before the absorbing A_next[A_prev==1]<-1 override is
##         applied), so for an already-switched patient pA_k2 is a leftover
##         counterfactual value for a decision that was never actually made.
##         Once switched, the next visit's action is deterministic and must
##         contribute probability 1 to the path, not pA_k2. Failing to gate
##         on the observed prior action (A_next_k{k-1}) distorts the path
##         probability disproportionately for early switchers -- i.e.
##         disproportionately for DTR1 -- and biases ESS/w99 for that regime.
## ============================================================================

verify_violation_scenarios_v6 <- function(
    n         = 20000,
    seed      = 999,
    scenarios = c("benchmark", "mild", "moderate", "severe"),
    dgp_fun   = simulate_hiv_violation_v6,
    ...) {

  get_regime_probs <- function(dat, regime) {

    wide <- dat %>%
      filter(month %in% c(3, 6, 9)) %>%
      mutate(kk = match(month, c(3, 6, 9))) %>%
      select(id, kk, A_next, det, fail, pA, unconfirmed_det_only) %>%
      pivot_wider(id_cols = id, names_from = kk,
                  values_from = c(A_next, det, fail, pA, unconfirmed_det_only),
                  names_sep = "_k")

    if (regime == "d1") {
      A1s <- as.integer(wide$det_k1 == 1)
      A2s <- pmax(A1s, as.integer(wide$det_k2 == 1))
      A3s <- pmax(A2s, as.integer(wide$det_k3 == 1))
    } else {
      A1s <- as.integer(wide$fail_k1 == 1)
      A2s <- pmax(A1s, as.integer(wide$fail_k2 == 1))
      A3s <- pmax(A2s, as.integer(wide$fail_k3 == 1))
    }

    ## Per-visit probability of the prescribed action.
    ## FIX: once a patient has already switched at a prior visit
    ## (A_next_k{k-1}==1), the current visit's action is deterministic
    ## (absorbing) and contributes probability 1 to the path -- not
    ## wide$pA_k{k}, which was computed for a decision never actually made.
    p1 <- ifelse(A1s == 1, wide$pA_k1, 1 - wide$pA_k1)

    entering2 <- wide$A_next_k1
    p2 <- ifelse(entering2 == 1, 1,
                 ifelse(A2s == 1, wide$pA_k2, 1 - wide$pA_k2))

    entering3 <- wide$A_next_k2
    p3 <- ifelse(entering3 == 1, 1,
                 ifelse(A3s == 1, wide$pA_k3, 1 - wide$pA_k3))

    f1   <- as.integer(wide$A_next_k1 == A1s)
    f2   <- as.integer(wide$A_next_k2 == A2s)
    f3   <- as.integer(wide$A_next_k3 == A3s)
    fall <- f1 * f2 * f3

    pp <- p1 * p2 * p3
    w  <- fall / pp
    wp <- w[fall == 1 & is.finite(w) & w > 0]

    ## Subgroup diagnostic: patients who were EVER in the unconfirmed
    ## det-only state at any of the three decision visits, regardless of
    ## which regime is being evaluated. This is the subgroup the switch-
    ## barrier mechanism actually targets -- most patients never become
    ## det_only at all and are structurally unaffected by the barrier, which
    ## dilutes whole-population summaries (see calibration note above: the
    ## whole-population p1p2p3 %<0.10 undershot hand-derived targets across
    ## every recalibration round even as the targeted mechanism clearly
    ## strengthened -- this subgroup split is why).
    ever_unconfirmed <- as.integer(
      wide$unconfirmed_det_only_k1 == 1 |
      wide$unconfirmed_det_only_k2 == 1 |
      wide$unconfirmed_det_only_k3 == 1
    )

    stats_for <- function(idx) {
      wp_i <- w[idx & fall == 1 & is.finite(w) & w > 0]
      if (length(wp_i) == 0) {
        return(c(ess = NA_real_, w99 = NA_real_, support = NA_real_))
      }
      c(ess     = (sum(wp_i)^2) / sum(wp_i^2),
        w99     = as.numeric(quantile(wp_i, 0.99, na.rm = TRUE)),
        support = mean(fall[idx], na.rm = TRUE))
    }
    whole <- stats_for(rep(TRUE, length(fall)))
    sub   <- stats_for(ever_unconfirmed == 1)

    list(
      p1      = p1,
      p12     = p1 * p2,
      p123    = pp,
      followed= fall,
      ess     = whole["ess"],  w99 = whole["w99"],  support = whole["support"],
      ever_unconfirmed        = ever_unconfirmed,
      n_ever_unconfirmed      = sum(ever_unconfirmed == 1),
      p123_sub_median         = round(median(pp[ever_unconfirmed == 1], na.rm = TRUE), 3),
      p123_sub_pct_lt_010     = round(100 * mean(pp[ever_unconfirmed == 1] < 0.10, na.rm = TRUE), 1),
      ess_sub = sub["ess"],  w99_sub = sub["w99"],  support_sub = sub["support"])
  }

  sp <- function(x) c(
    median     = round(median(x, na.rm = TRUE), 3),
    p05        = round(as.numeric(quantile(x, 0.05, na.rm = TRUE)), 3),
    p01        = round(as.numeric(quantile(x, 0.01, na.rm = TRUE)), 3),
    pct_lt_010 = round(100 * mean(x < 0.10, na.rm = TRUE), 1),
    pct_lt_005 = round(100 * mean(x < 0.05, na.rm = TRUE), 1))

  out <- lapply(scenarios, function(sc) {

    dat <- dgp_fun(n = n, seed = seed, scenario = sc, ...)
    dec <- dat %>% filter(month %in% c(3, 6, 9), A_prev == 0)

    ## v6 bands: confirmation status is the primary split, since that's
    ## what the switch-barrier mechanism actually acts on (persist/deep
    ## remain available as covariates but are no longer the band definition).
    band <- dplyr::case_when(
      dec$fail == 1 ~ "failure",
      dec$det_only == 1 & dec$confirmed_rebound   == 1 ~ "confirmed_det_only",
      dec$det_only == 1 & dec$unconfirmed_det_only == 1 ~ "unconfirmed_det_only",
      TRUE ~ "suppressed")

    med_by_band <- function(b) {
      x <- dec$pA[band == b]
      if (length(x) == 0) return(c(n = 0L, median = NA_real_, pct_lt_010 = NA_real_))
      c(n = length(x),
        median = round(median(x, na.rm = TRUE), 3),
        pct_lt_010 = round(100 * mean(x < 0.10, na.rm = TRUE), 1))
    }
    b_supp   <- med_by_band("suppressed")
    b_unconf <- med_by_band("unconfirmed_det_only")
    b_conf   <- med_by_band("confirmed_det_only")
    b_fail   <- med_by_band("failure")

    d1 <- get_regime_probs(dat, "d1")
    d2 <- get_regime_probs(dat, "d2")

    list(
      visit_level = data.frame(
        scenario               = sc,
        n_suppressed           = b_supp["n"],   median_pA_suppressed   = b_supp["median"],
        n_unconfirmed_det_only = b_unconf["n"], median_pA_unconfirmed  = b_unconf["median"],
        pct_lt10_unconfirmed   = b_unconf["pct_lt_010"],
        n_confirmed_det_only   = b_conf["n"],   median_pA_confirmed    = b_conf["median"],
        n_failure              = b_fail["n"],   median_pA_failure      = b_fail["median"]),

      path_level = bind_rows(
        data.frame(scenario = sc, regime = "d1", stage = "p1",     t(sp(d1$p1))),
        data.frame(scenario = sc, regime = "d1", stage = "p1p2",   t(sp(d1$p12))),
        data.frame(scenario = sc, regime = "d1", stage = "p1p2p3", t(sp(d1$p123))),
        data.frame(scenario = sc, regime = "d2", stage = "p1",     t(sp(d2$p1))),
        data.frame(scenario = sc, regime = "d2", stage = "p1p2",   t(sp(d2$p12))),
        data.frame(scenario = sc, regime = "d2", stage = "p1p2p3", t(sp(d2$p123)))),

      final_level = data.frame(
        scenario   = sc,
        support_d1 = round(d1$support, 3), ess_d1 = round(d1$ess, 1), w99_d1 = round(d1$w99, 1),
        support_d2 = round(d2$support, 3), ess_d2 = round(d2$ess, 1), w99_d2 = round(d2$w99, 1)),

      ## Subgroup diagnostic: DTR1/DTR2 path-probability behavior restricted
      ## to patients who ever entered the unconfirmed det-only state at any
      ## decision visit -- the subgroup the switch-barrier mechanism targets,
      ## as distinct from the whole population (most of whom never become
      ## det_only and are structurally unaffected). Compare *_whole against
      ## *_sub columns: the subgroup numbers show the targeted violation far
      ## more clearly than the whole-population final_level table above.
      subgroup_level = data.frame(
        scenario                = sc,
        n_ever_unconfirmed      = d1$n_ever_unconfirmed,
        pct_ever_unconfirmed    = round(100 * d1$n_ever_unconfirmed / nrow(dec %>% dplyr::distinct(id)), 1),
        d1_p123_median_whole    = round(median(d1$p123, na.rm = TRUE), 3),
        d1_p123_median_sub      = d1$p123_sub_median,
        d1_pct_lt10_whole       = round(100 * mean(d1$p123 < 0.10, na.rm = TRUE), 1),
        d1_pct_lt10_sub         = d1$p123_sub_pct_lt_010,
        d1_ess_whole            = round(d1$ess, 1),
        d1_ess_sub              = round(d1$ess_sub, 1),
        d1_w99_sub              = round(d1$w99_sub, 1),
        d2_p123_median_sub      = d2$p123_sub_median,
        d2_pct_lt10_sub         = d2$p123_sub_pct_lt_010,
        d2_ess_sub              = round(d2$ess_sub, 1)))
  })

  list(
    visit_level    = bind_rows(lapply(out, `[[`, "visit_level")),
    path_level     = bind_rows(lapply(out, `[[`, "path_level")),
    final_level    = bind_rows(lapply(out, `[[`, "final_level")),
    subgroup_level = bind_rows(lapply(out, `[[`, "subgroup_level")))
}


## ============================================================================
## 3.  TRUE RR UNDER FORCED REGIMES
##
##     Computed by bypassing the switching law entirely.
##     All patients are assigned to DTR1 (or DTR2) deterministically.
##     This is identical across all scenarios because the biology is fixed.
##     TRUE_RR should be stable regardless of the violation scenario used.
## ============================================================================

TRUE_RR    <- 1.469604947
TRUE_EY_D1 <- 0.5716074
TRUE_EY_D2 <- 0.3889531

