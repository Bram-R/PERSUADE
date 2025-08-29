library(testthat)

# Tests for plotting and output helpers defined in R/PERSUADE_output_functions.R:
#
# These tests include:
# 1) Verification of f_summary on simple numeric data.frames, checking
#    descriptive statistics and output structure.
# 2) Error handling checks for f_generate_report when the required
#    Rmd template is missing.
# 3) Lightweight plotting tests (e.g., f_plot_km_survival_base,
#    f_plot_log_cumhaz) using minimal synthetic PERSUADE-like objects.
# 4) Conditional tests of survminer-based KM plotting and Schoenfeld residual
#    plotting with fitted Cox models, guarded by skip_on_cran /
#    skip_if_not_installed.
# 5) Smoothed hazard plotting with minimal hazard data structures.
# 6) Extensive overlay and extrapolation plot functions, exercised with
#    synthetic PERSUADE objects containing parametric, spline, and cure models,
#    ensuring all plotting helpers run without error.
# 7) Integration-style checks of plot(PERSUADE, type=...) on a full
#    f_PERSUADE run (flexsurv::bc), confirming end-to-end compatibility.
#
# Heavy/integration tests are guarded with skip_on_cran()/skip_if_not_installed().

test_that("f_summary computes descriptive stats for numeric data.frame", {
  df <- mtcars[, c("mpg", "hp", "wt")]
  res <- f_summary(df)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("Mean", "Std.Dev", "Min", "Q1", "Median", "Q3", "Max", "IQR") %in% colnames(res)))
  # Rownames correspond to variables
  expect_true(all(rownames(res) == colnames(df)))
  # Values numeric
  expect_true(all(vapply(res, is.numeric, logical(1))))
})

test_that("f_generate_report errors when Rmd template is missing", {
  # Minimal PERSUADE object
  PERS <- list(name = "tmp_report")
  class(PERS) <- "PERSUADE"

  # The function checks system.file(...) == "" and stops with a specific message
  expect_error(f_generate_report(PERS), "The PERSUADE Rmd template was not found in the package")
})

test_that("f_plot_km_survival_base and f_plot_log_cumhaz run with a simple survfit object", {
  skip_on_cran()
  skip_if_not_installed("survival")

  # Create a tiny dataset and survfit for single-group plotting
  df <- data.frame(
    time = c(1,2,3,4,5,6),
    status = c(1,1,0,1,0,1),
    group = factor(c("A","A","A","A","A","A"))
  )
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = df)

  PERS <- list(
    input = list(years = df$time, status = df$status, group = df$group),
    misc = list(ngroups = 1, group_names = "A"),
    surv_obs = list(
      km = km,
      km_names = rep(1, length(km$time))
    )
  )
  class(PERS) <- "PERSUADE"

  # base KM plot should run without error (draws to current device)
  expect_silent(f_plot_km_survival_base(PERS))

  # log-cumulative-hazard diagnostic expects km with surv and time
  expect_silent(f_plot_log_cumhaz(PERS))
})

test_that("f_plot_km_survival (survminer) runs when survminer is available", {
  skip_on_cran()
  skip_if_not_installed("survminer")
  skip_if_not_installed("survival")
  # build simple two-group dataset
  df <- data.frame(
    years = c(1,2,3,4,5,6,7,8),
    status = c(1,1,0,1,0,1,0,1),
    group = factor(c("A","A","A","A","B","B","B","B"))
  )
  # build PERSUADE-like object
  PERS <- list(
    input = list(years = df$years, status = df$status, group = df$group),
    misc = list(ngroups = 2, group_names = c("A","B"), form = survival::Surv(years, status) ~ group)
  )
  class(PERS) <- "PERSUADE"

  expect_silent(f_plot_km_survival(PERS))
})

test_that("f_plot_schoenfeld_residuals works with a fitted cox model (2 groups)", {
  skip_on_cran()
  skip_if_not_installed("survival")

  df <- data.frame(
    time = c(1,2,3,4,5,6,7,8),
    status = c(1,1,0,1,0,1,0,1),
    group = factor(c("A","A","A","B","B","B","B","A"))
  )

  # Fit a simple Cox model with group as covariate
  cox_fit <- survival::coxph(survival::Surv(time, status) ~ group, data = df)

  # Build a PERSUADE-like object
  km <- survival::survfit(survival::Surv(time, status) ~ group, data = df)
  PERS <- list(
    surv_obs = list(cox_reg = cox_fit, km = km, km_names = rep(1, length(km$time))),
    misc = list(ngroups = 2, group_names = c("A", "B"))
  )
  class(PERS) <- "PERSUADE"

  # Should run (produces base plots)
  expect_silent(f_plot_schoenfeld_residuals(PERS))
})

test_that("f_plot_smoothed_hazard runs with synthetic hazard data", {
  skip_on_cran()

  # Build minimal smoothed hazard structure for one group
  haz_df <- data.frame(est.grid = seq(0, 10, 1), haz.est = seq(0.1, 1.1, length.out = 11))
  surv_obs <- list(haz = list(max = data.frame(time = 10, smooth = max(haz_df$haz.est)),
                              hazards = list(smooth_gr1 = haz_df),
                              names = 1))
  PERS <- list(surv_obs = surv_obs, misc = list(ngroups = 1, group_names = "G1"))
  class(PERS) <- "PERSUADE"

  expect_silent(f_plot_smoothed_hazard(PERS))
})

test_that("model-overlay plot functions run with minimal model prediction objects", {
  skip_on_cran()
  skip_if_not_installed("data.table")

  # Construct minimal PERSUADE object with required structures for overlay plotting
  time_pred <- seq(0, 10, by = 1)
  # helper to create model matrices (time + single group column)
  make_model_matrices <- function() {
    surv_mat <- cbind(time_pred, matrix(seq(1, by = -0.05, length.out = length(time_pred)), ncol = 1))
    hazard_mat <- cbind(time_pred, matrix(seq(0.2, by = 0.02, length.out = length(time_pred)), ncol = 1))
    list(surv = as.data.frame(surv_mat), hazard = as.data.frame(hazard_mat))
  }

  # Create param/spline/cure model name lists to match internal mapping expectations
  param_names  <- c("expo","weib","gom","lnorm","llog","gam","ggam")
  spline_names <- c("spline_hazard_1","spline_odds_1")
  cure_names   <- c("weib_mixture_1","lnorm_nonmix_1")

  surv_pred <- list(model = list())
  for (nm in param_names) surv_pred$model[[nm]] <- make_model_matrices()
  surv_pred$model$spline <- list()
  for (nm in spline_names) surv_pred$model$spline[[nm]] <- make_model_matrices()
  surv_pred$model$cure <- list()
  for (nm in cure_names) surv_pred$model$cure[[nm]] <- make_model_matrices()

  surv_model <- list(
    param_models = setNames(vector("list", length(param_names)), param_names),
    spline_models = setNames(vector("list", length(spline_names)), spline_names),
    cure_models = setNames(vector("list", length(cure_names)), cure_names)
  )
  # Provide knots for spline models (used when drawing vertical lines)
  for (nm in spline_names) surv_model$spline_models[[nm]] <- list(knots = c(log(1), log(3), log(10)))

  # smoothed hazard to serve as base
  haz_df <- data.frame(est.grid = time_pred, haz.est = seq(0.1, 1.1, length.out = length(time_pred)))
  surv_obs <- list(
    haz = list(max = data.frame(time = max(time_pred), smooth = max(haz_df$haz.est)),
               hazards = list(smooth_gr1 = haz_df, smooth_gr2 = haz_df),
               names = c(1,2)),
    km = structure(list(time = time_pred, surv = seq(1, 0.5, length.out = length(time_pred))),
                   class = "survfit"),
    km_names = rep(1, length(time_pred)),
    tp = list(gr_1 = data.frame(time = seq(0,5), smooth = runif(6), smooth_upper = runif(6, .6, 1), smooth_lower = runif(6, 0, .4)),
              gr_2 = data.frame(time = seq(0,5), smooth = runif(6), smooth_upper = runif(6, .6, 1), smooth_lower = runif(6, 0, .4))),
    tp_gr = list()
  )

  # Build tp_gr as data.table with many columns to satisfy column-indexed selection
  dt_cols <- data.table::data.table(Time = seq(0,5))
  for (j in 2:23) dt_cols[[paste0("V", j)]] <- runif(6)
  surv_pred$tp_gr <- list(gr_1 = data.table::copy(dt_cols), gr_2 = data.table::copy(dt_cols))

  PERS <- list(
    input = list(time_pred = time_pred, spline_mod = TRUE, cure_mod = TRUE, time_horizon = 10),
    misc = list(ngroups = 2, group_names = c("G1","G2"), lbls = paste0("M", seq_along(param_names)),
                lbls_spline = paste0("S", seq_along(spline_names)), lbls_cure = paste0("C", seq_along(cure_names))),
    surv_obs = surv_obs,
    surv_pred = surv_pred,
    surv_model = surv_model
  )
  class(PERS) <- "PERSUADE"

  # Call plotting functions to ensure they run without error.
  expect_silent(f_plot_hazard_with_models(PERS))
  expect_silent(f_plot_param_surv_model(PERS, model_index = 1))
  expect_silent(f_plot_spline_surv_model(PERS, model_index = 1))
  expect_silent(f_plot_cure_surv_model(PERS, model_index = 1))

  expect_silent(f_plot_diag_param_surv_model(PERS, model_index = 1))
  expect_silent(f_plot_diag_spline_surv_model(PERS, model_index = 1))
  expect_silent(f_plot_diag_cure_surv_model(PERS, model_index = 1))

  expect_silent(f_plot_tp_param_surv_model(PERS, model_index = 1))
  expect_silent(f_plot_tp_spline_surv_model(PERS, model_index = 1))
  expect_silent(f_plot_tp_cure_surv_model(PERS, model_index = 1))

  expect_silent(f_plot_param_surv_extrap(PERS))
  expect_silent(f_plot_spline_surv_extrap(PERS))
  expect_silent(f_plot_cure_surv_extrap(PERS))
  expect_silent(f_plot_tp_param_surv_extrap(PERS))
  expect_silent(f_plot_tp_spline_surv_extrap(PERS))
  expect_silent(f_plot_tp_cure_surv_extrap(PERS))
  expect_silent(f_plot_hazard_parametric_extrap(PERS))
  expect_silent(f_plot_hazard_spline_extrap(PERS))
  expect_silent(f_plot_hazard_cure_extrap(PERS))
})

test_that("plot functions run without error (integration with f_PERSUADE and flexsurv::bc)", {
  skip_on_cran()
  skip_if_not_installed("flexsurv")

  data <- flexsurv::bc
  PERSUADE <- f_PERSUADE(
    name = "test_plots",
    years = data$recyrs,
    status = data$censrec,
    group = data$group,
    strata = TRUE,
    spline_mod = FALSE,
    cure_mod = FALSE,
    time_unit = 1/12,
    time_horizon = 5,
    time_pred_surv_table = c(0, 1, 2, 3, 4, 5)
  )

  expect_silent(plot(PERSUADE, type = "km"))
  expect_silent(plot(PERSUADE, type = "ph"))
  expect_silent(plot(PERSUADE, type = "hr"))
})
