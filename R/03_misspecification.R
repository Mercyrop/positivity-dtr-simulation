## ============================================================================
## 03_misspecification.R
##
## Prespecified nuisance-model misspecification analysis.
##
## Settings:
##   Correct : correctly specified treatment and outcome models
##   Pmis    : misspecified treatment model only
##   Qmis    : misspecified outcome-transition model only
##   Both    : both nuisance models misspecified
##
## Pmis retains det_only, fail, persist, age_c, and male, while omitting
## deep persistence and the confirmed/unconfirmed split.
## Qmis omits heterogeneous second-line response through A:slope_B.
##
## Source 01_data_generating_mechanism.R and 02_estimators.R first.
## ============================================================================

options(cli.dynamic = FALSE, cli.progress_show_after = Inf)

MISSPEC_SETTINGS_V7E <- data.frame(
  setting = c("Correct", "Pmis", "Qmis", "Both"),
  g_mis   = c(FALSE, TRUE,  FALSE, TRUE),
  q_mis   = c(FALSE, FALSE, TRUE,  TRUE),
  stringsAsFactors = FALSE
)

required_misspec_functions_v7e <- c(
  "simulate_hiv_violation_v6",
  "prepare_wide_estimation_v7",
  "make_seed_v7",
  "normal_quadrature",
  "fit_ps_v7",
  "fit_oracle_transition_v7",
  "fit_stabilization_numerator_v7c",
  "iptw_msm_v7",
  "iptw_msm_trunc01_99_v7",
  "iptw_msm_stabilized_v7c",
  "iptw_msm_stabilized_trunc_v7c",
  "aiptw_rr_v7",
  "tmle_rr_v7",
  "gcomp_rr_v7",
  "safe_ratio_result",
  "safe_iptw_trunc01_99_result_v7b",
  "safe_iptw_stab_result_v7c",
  "ps_xvars_v7",
  "transition_x_v7"
)

validate_misspec_code_v7e <- function() {
  missing_functions <- required_misspec_functions_v7e[
    !vapply(required_misspec_functions_v7e, exists, logical(1), mode = "function")
  ]
  if (length(missing_functions) > 0L) {
    stop(
      "Required v7e function(s) are missing:\n  ",
      paste(missing_functions, collapse = "\n  "),
      "\nSource the supplied DGP and estimator files first."
    )
  }

  if (!exists("TRUE_RR", inherits = TRUE) ||
      abs(get("TRUE_RR", inherits = TRUE) - 1.469604947) > 1e-10) {
    stop("The estimator file does not contain the final TRUE_RR = 1.469604947.")
  }

  dummy <- data.frame(
    fail1 = 0, conf1 = 0, unconf1 = 0, persist1 = 0, deep1 = 0,
    det_only1 = 0, age_c = 0, male = 0
  )

  correct_g_names <- names(ps_xvars_v7(dummy, 1L, misspecified = FALSE))
  pmis_g_names <- names(ps_xvars_v7(dummy, 1L, misspecified = TRUE))

  expected_correct_g <- c("fail", "conf", "unconf", "persist", "deep", "age_c", "male")
  expected_pmis_g <- c("det_only", "fail", "persist", "age_c", "male")

  if (!identical(correct_g_names, expected_correct_g)) {
    stop(
      "Correct treatment-model specification has drifted.\nExpected: ",
      paste(expected_correct_g, collapse = ", "),
      "\nFound: ", paste(correct_g_names, collapse = ", ")
    )
  }
  if (!identical(pmis_g_names, expected_pmis_g)) {
    stop(
      "Pmis specification is not the accepted softened model.\nExpected: ",
      paste(expected_pmis_g, collapse = ", "),
      "\nFound: ", paste(pmis_g_names, collapse = ", ")
    )
  }

  correct_q_names <- names(
    transition_x_v7(0, 0, 0, 0, misspecified = FALSE)
  )
  qmis_q_names <- names(
    transition_x_v7(0, 0, 0, 0, misspecified = TRUE)
  )

  expected_correct_q <- c("logvl_current", "slope_A", "A_slope_A", "A_slope_B")
  expected_qmis_q <- c("logvl_current", "slope_A", "A_slope_A", "A")

  if (!identical(correct_q_names, expected_correct_q)) {
    stop(
      "Correct Q-model specification has drifted.\nExpected: ",
      paste(expected_correct_q, collapse = ", "),
      "\nFound: ", paste(correct_q_names, collapse = ", ")
    )
  }
  if (!identical(qmis_q_names, expected_qmis_q)) {
    stop(
      "Qmis is not the accepted moderated specification.\nExpected: ",
      paste(expected_qmis_q, collapse = ", "),
      "\nFound: ", paste(qmis_q_names, collapse = ", ")
    )
  }

  invisible(TRUE)
}

safe_named_v7e <- function(prefix, expr, fallback) {
  out <- tryCatch(
    expr,
    error = function(e) {
      message("[", prefix, "] ", conditionMessage(e))
      fallback
    }
  )
  names(out) <- paste0(prefix, "_", names(out))
  out
}

compute_iptw_branch_v7e <- function(wide, ps, numerator_fit) {
  c(
    safe_named_v7e(
      "IPTW_PAR",
      iptw_msm_v7(wide, ps$g1, ps$g2, ps$g3),
      safe_ratio_result()
    ),
    safe_named_v7e(
      "IPTW_TRUNC01_99_PAR",
      iptw_msm_trunc01_99_v7(wide, ps$g1, ps$g2, ps$g3),
      safe_iptw_trunc01_99_result_v7b()
    ),
    safe_named_v7e(
      "IPTW_STAB_PAR",
      iptw_msm_stabilized_v7c(
        wide, ps$g1, ps$g2, ps$g3,
        numerator_fit = numerator_fit
      ),
      safe_iptw_stab_result_v7c(FALSE)
    ),
    safe_named_v7e(
      "IPTW_STAB_TRUNC01_99_PAR",
      iptw_msm_stabilized_trunc_v7c(
        wide, ps$g1, ps$g2, ps$g3,
        numerator_fit = numerator_fit
      ),
      safe_iptw_stab_result_v7c(TRUE)
    )
  )
}

compute_dr_branch_v7e <- function(wide, ps, qfit, quadrature) {
  c(
    safe_named_v7e(
      "AIPTW_PAR",
      aiptw_rr_v7(
        wide, ps$g1, ps$g2, ps$g3, qfit, quadrature
      ),
      safe_ratio_result()
    ),
    safe_named_v7e(
      "TMLE_PAR",
      tmle_rr_v7(
        wide, ps$g1, ps$g2, ps$g3, qfit, quadrature
      ),
      safe_ratio_result()
    )
  )
}

compute_gcomp_branch_v7e <- function(wide, qfit, quadrature) {
  safe_named_v7e(
    "GCOMP_PAR",
    gcomp_rr_v7(wide, qfit, quadrature),
    safe_ratio_result()
  )
}

one_sim_misspec_v7e <- function(
    n,
    scenario,
    sim,
    n_quad = 15L,
    dtr1_threshold = 50,
    dtr2_threshold = 1000,
    g_bounds = c(0.01, 0.99)
) {
  validate_misspec_code_v7e()

  n <- as.integer(n)
  sim <- as.integer(sim)
  n_quad <- as.integer(n_quad)
  seed_i <- make_seed_v7(n = n, scenario = scenario, sim = sim)

  dat <- simulate_hiv_violation_v6(
    n = n,
    seed = seed_i,
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

  ## Fit each distinct nuisance model once and reuse it across settings.
  ps_correct <- fit_ps_v7(
    wide,
    use_sl = FALSE,
    misspecified = FALSE,
    g_bounds = g_bounds,
    seed = seed_i + 1000L
  )
  ps_pmis <- fit_ps_v7(
    wide,
    use_sl = FALSE,
    misspecified = TRUE,
    g_bounds = g_bounds,
    seed = seed_i + 1000L
  )

  q_correct <- fit_oracle_transition_v7(
    wide,
    use_sl = FALSE,
    misspecified = FALSE,
    seed = seed_i
  )
  q_qmis <- fit_oracle_transition_v7(
    wide,
    use_sl = FALSE,
    misspecified = TRUE,
    seed = seed_i
  )

  numerator_fit <- tryCatch(
    fit_stabilization_numerator_v7c(wide),
    error = function(e) {
      message("[stabilization numerator] ", conditionMessage(e))
      NULL
    }
  )

  iptw_branches <- list(
    Correct = compute_iptw_branch_v7e(wide, ps_correct, numerator_fit),
    Pmis    = compute_iptw_branch_v7e(wide, ps_pmis, numerator_fit)
  )

  q_branches <- list(
    Correct = q_correct,
    Qmis    = q_qmis
  )

  gcomp_branches <- list(
    Correct = compute_gcomp_branch_v7e(wide, q_correct, quadrature),
    Qmis    = compute_gcomp_branch_v7e(wide, q_qmis, quadrature)
  )

  setting_rows <- lapply(seq_len(nrow(MISSPEC_SETTINGS_V7E)), function(i) {
    spec <- MISSPEC_SETTINGS_V7E[i, ]
    setting_i <- spec$setting
    g_key <- if (isTRUE(spec$g_mis)) "Pmis" else "Correct"
    q_key <- if (isTRUE(spec$q_mis)) "Qmis" else "Correct"
    ps_i <- if (g_key == "Pmis") ps_pmis else ps_correct
    q_i <- q_branches[[q_key]]

    est <- c(
      iptw_branches[[g_key]],
      compute_dr_branch_v7e(wide, ps_i, q_i, quadrature),
      gcomp_branches[[q_key]]
    )

    est_columns <- grepl("_est$", names(est))
    estimator_failure_count <- sum(
      !is.finite(as.numeric(est[est_columns]))
    )

    cbind(
      data.frame(
        code_release = "github_public_v1",
        n = n,
        scenario = scenario,
        sim = sim,
        seed = seed_i,
        setting = setting_i,
        g_mis = isTRUE(spec$g_mis),
        q_mis = isTRUE(spec$q_mis),
        failed = FALSE,
        estimator_failure_count = estimator_failure_count,
        error_message = NA_character_,
        n_quad = n_quad,
        stringsAsFactors = FALSE
      ),
      as.data.frame(as.list(est), check.names = FALSE)
    )
  })

  dplyr::bind_rows(setting_rows)
}

failed_misspec_rows_v7e <- function(n, scenario, sim, error_message, n_quad = NA_integer_) {
  seed_i <- tryCatch(
    make_seed_v7(n = n, scenario = scenario, sim = sim),
    error = function(e) NA_integer_
  )

  dplyr::bind_rows(lapply(seq_len(nrow(MISSPEC_SETTINGS_V7E)), function(i) {
    spec <- MISSPEC_SETTINGS_V7E[i, ]
    data.frame(
      code_release = "github_public_v1",
      n = as.integer(n),
      scenario = as.character(scenario),
      sim = as.integer(sim),
      seed = as.integer(seed_i),
      setting = spec$setting,
      g_mis = isTRUE(spec$g_mis),
      q_mis = isTRUE(spec$q_mis),
      failed = TRUE,
      estimator_failure_count = NA_integer_,
      error_message = as.character(error_message),
      n_quad = as.integer(n_quad),
      stringsAsFactors = FALSE
    )
  }))
}
