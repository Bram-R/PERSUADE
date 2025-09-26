library(testthat)

# Tests for extended PERSUADE functionality:
#
# These tests include:
# 1) Integration-style checks of helper functions (f_hazard, f_cum_hazard, f_tp)
#    to ensure consistent shapes, valid ranges, and expected column structures.
# 2) Full-run evaluations of f_PERSUADE on larger survival datasets
#    (e.g., survival::lung), verifying parametric model predictions,
#    group-specific outputs, and numerical plausibility of survival/hazard values.
# 3) Validation of Excel-export tables (f_surv_model_excel) for correct
#    presence of model metadata such as Distnames, Parnames, and covariates.
# 4) Checks of truncation logic in f_surv_model_pred_tp_gr, ensuring values
#    beyond thresholds are correctly blanked/NA.
# 5) Tests of input transformations and misc attributes when advanced options
#    (spline_mod, cure_mod, cure_link) are enabled, confirming consistent
#    propagation through input and misc slots.
#
# These are more integration-like and will be skipped when dependencies are missing.

skip_on_cran()

# These tests are more integration-like and will be skipped when dependencies are missing.

test_that("helper functions f_hazard, f_cum_hazard, f_tp produce consistent shapes", {
  skip_if_not_installed("muhaz")
  skip_if_not_installed("sft")

  years_small <- c(1, 2, 3, 4, 5, 6)
  status_small <- c(1, 1, 0, 1, 0, 1)
  group_A <- factor(rep("A", length(years_small)))

  # f_hazard output basic checks
  haz_out <- f_hazard(years_small, status_small, group_A, ngroups = 1)
  expect_type(haz_out, "list")
  expect_true(all(c("hazards", "names", "max") %in% names(haz_out)))
  expect_true(is.data.frame(haz_out$max))
  expect_true(is.numeric(haz_out$max$time) || is.numeric(haz_out$max$t))

  # f_cum_hazard returns a data.frame with expected columns and non-negative H/var
  time_pred <- seq(0, 6, by = 1)
  cumh <- f_cum_hazard(years_small, status_small, group_A, ngroups = 1, time_pred = time_pred, time_unit = 1)
  expect_s3_class(cumh, "data.frame")
  expect_true(all(c("group", "time", "H", "var") %in% names(cumh)))
  expect_true(all(cumh$H >= 0))
  expect_true(all(cumh$var >= 0))

  # f_tp returns smoothed probabilities with bounds [0,1] (except possible leading NA)
  tp_out <- f_tp(ngroups = 1, cum_haz = cumh, time_unit = 1)
  expect_type(tp_out, "list")
  expect_true("gr_1" %in% names(tp_out))
  expect_true(is.data.frame(tp_out$gr_1))
  if (nrow(tp_out$gr_1) > 1) {
    vals <- unlist(tp_out$gr_1[-1, intersect(c("smooth", "smooth_upper", "smooth_lower"), colnames(tp_out$gr_1))], use.names = FALSE)
    vals <- vals[!is.na(vals)]
    if (length(vals)) expect_true(all(vals >= 0 & vals <= 1))
  }
})

test_that("f_PERSUADE full run: parametric model predictions are sensible", {
  # require packages for parametric fits and predictions
  required_pkgs <- c("flexsurv", "muhaz", "sft", "rms", "data.table")
  miss <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) skip(paste("Missing packages:", paste(miss, collapse = ", ")))

  yrs <- survival::lung$time
  sts <- survival::lung$status
  grp <- factor(survival::lung$sex) # two groups

  out <- f_PERSUADE(
    name = "exhaustive_check",
    years = yrs,
    status = sts,
    group = grp,
    strata = FALSE,
    spline_mod = FALSE,   # parametric models only to keep runtime reasonable
    cure_mod = FALSE,
    time_unit = 365.25 / 12,
    time_horizon = 2000,
    time_pred_surv_table = seq(0, 2000, by = 365.25)
  )

  expect_s3_class(out, "PERSUADE")

  # 1) Survival estimates between 0 and 1 for parametric models (if present)
  if (!is.null(out$surv_pred$model) && length(out$surv_pred$model) > 0) {
    for (mod in names(out$surv_pred$model)) {
      model_entry <- out$surv_pred$model[[mod]]
      if (!is.null(model_entry$surv)) {
        pred_df <- as.data.frame(model_entry$surv)
        if (ncol(pred_df) > 1) {
          surv_vals <- as.numeric(unlist(pred_df[, -1, drop = FALSE]))
          surv_vals <- surv_vals[!is.na(surv_vals)]
          if (length(surv_vals)) expect_true(all(surv_vals >= 0 & surv_vals <= 1))
        }
      }
      if (!is.null(model_entry$hazard)) {
        haz_df <- as.data.frame(model_entry$hazard)
        if (ncol(haz_df) > 1) {
          haz_vals <- as.numeric(unlist(haz_df[, -1, drop = FALSE]))
          haz_vals <- haz_vals[!is.na(haz_vals)]
          if (length(haz_vals)) expect_true(all(haz_vals >= 0))
        }
      }
    }
  }

  # 2) Grouped predictions should exist for each group
  expect_true(is.list(out$surv_pred$gr))
  expect_true(length(out$surv_pred$gr) == out$misc$ngroups)
  for (gname in names(out$surv_pred$gr)) {
    tbl <- out$surv_pred$gr[[gname]]
    expect_true(is.data.frame(tbl) || is.matrix(tbl))
    expect_true("time" %in% colnames(tbl) || "Time" %in% colnames(tbl))
  }
})

test_that("f_surv_model_excel returns a table with Distnames and Parnames rows when models run", {
  required_pkgs <- c("flexsurv", "muhaz", "sft", "rms")
  miss <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) skip(paste("Missing packages:", paste(miss, collapse = ", ")))

  yrs <- survival::lung$time
  sts <- survival::lung$status
  grp <- factor(survival::lung$sex)

  out <- f_PERSUADE(
    name = "excel_check",
    years = yrs,
    status = sts,
    group = grp,
    strata = FALSE,
    spline_mod = FALSE,
    cure_mod = FALSE,
    time_unit = 365.25 / 12,
    time_horizon = 2000,
    time_pred_surv_table = seq(0, 2000, by = 365.25)
  )

  excel_df <- out$surv_model_excel
  expect_s3_class(excel_df, "data.frame")
  rn <- rownames(excel_df)
  expect_true(any(c("Distnames", "Parnames") %in% rn))
  expect_true(any(grepl("^Cov_", rn)))
})

test_that("f_surv_model_pred_tp_gr truncation blanks values after exceedance threshold", {
  # This test does not require heavy modelling packages; it verifies truncation logic.
  group_tbl <- data.frame(
    time = 0:10,
    VerySharpModel_1 = c(1.0, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.0001, 0.0001, 0.0001),
    VerySharpModel_2 = c(1.0, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0, 0.0001)
  )

  surv_model_pred_gr <- list(gr_1 = group_tbl)
  res <- f_surv_model_pred_tp_gr(ngroups = 1, time_pred = 0:3, time_unit = 1, surv_model_pred_gr = surv_model_pred_gr, cols_tp = 3)

  expect_true(is.list(res))
  expect_true("gr_1" %in% names(res))
  df_res <- res$gr_1
  probcols <- setdiff(colnames(df_res), "Time")
  found_na <- FALSE
  for (pc in probcols) {
    if (any(is.na(df_res[[pc]]))) found_na <- TRUE
  }
  expect_true(found_na)
})

test_that("input transformations and misc values are consistent when spline and cure options are used", {
  required_pkgs <- c("muhaz", "sft", "rms", "flexsurv")
  miss <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) skip(paste("Missing packages:", paste(miss, collapse = ", ")))

  # Minimal working dataset
  set.seed(123)
  n_per_group <- 50   # 50 observations per group
  years <- rep(seq(1, 30, length.out = n_per_group), 2)   # Continuous years between 1 and 30
  status <- sample(0:1, 2 * n_per_group, replace = TRUE)   # Random events, balanced
  group <- factor(rep(c("A", "B"), each = n_per_group))   # Two groups

  out <- f_PERSUADE(
    name = "input_transform_check",
    years = years,
    status = status,
    group = group,
    strata = FALSE,
    spline_mod = TRUE,
    cure_mod = TRUE,
    cure_link = "logistic",
    time_unit = 5,
    time_horizon = 30,
    time_pred_surv_table = seq(0, 30, by = 10)
  )

  # time_pred should be numeric and cover horizon
  expect_true(is.numeric(out$input$time_pred))
  expect_true(max(out$input$time_pred) >= out$input$time_horizon)

  # time_pred_surv_table in input should be numeric and sorted
  expect_true(is.numeric(out$input$time_pred_surv_table))
  expect_equal(sort(out$input$time_pred_surv_table), out$input$time_pred_surv_table)

  # cure_link should be preserved when cure_mod = TRUE
  expect_true(!is.null(out$input$cure_link))
  expect_true(out$input$cure_link == "logistic")

  # cols_tp should reflect additions for spline and cure (numeric)
  expect_true(is.numeric(out$misc$cols_tp))
  expect_true(out$misc$cols_tp >= 8)
})
