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
  skip_if_not_installed("kableExtra")
  skip_if_not_installed("knitr")

  # Minimal PERSUADE object
  PERS <- list(name = "tmp_report")
  class(PERS) <- "PERSUADE"

  # The function checks system.file(...) == "" and stops with a specific message
  expect_error(f_generate_report(PERS, template_path = ""), "The PERSUADE Rmd template was not found in the package")
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

  expect_error(suppressWarnings(f_plot_km_survival(PERS)), NA) # so the test robust and doesn’t fail on cosmetic ggplot2 warnings
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

test_that("model-overlay plot functions run with minimal model prediction objects (fixed)", {
  skip_on_cran()
  skip_if_not_installed("survival")

  # Use a small positive start to avoid log(0) in diagnostic transforms
  time_pred <- seq(0.01, 10, by = 1)      # first element > 0
  pred_rows <- length(time_pred) - 1
  ngroups <- 2L

  eps <- 1e-6   # clip survival probabilities into (eps, 1-eps)

  # helper to create model matrices (time + one column per group)
  make_model_matrices <- function(ng, times) {
    surv_cols <- sapply(seq_len(ng), function(k)
      seq(from = 0.95 - (k - 1) * 0.05, by = -0.05, length.out = length(times))
    )
    hazard_cols <- sapply(seq_len(ng), function(k)
      seq(from = 0.2 + (k - 1) * 0.02, by = 0.02, length.out = length(times))
    )
    surv_mat <- cbind(time = times, as.data.frame(surv_cols))
    hazard_mat <- cbind(time = times, as.data.frame(hazard_cols))
    colnames(surv_mat) <- c("time", paste0("G", seq_len(ng)))
    colnames(hazard_mat) <- c("time", paste0("G", seq_len(ng)))
    # clip survival probabilities into (eps, 1 - eps)
    surv_mat[, -1] <- pmin(pmax(as.matrix(surv_mat[, -1]), eps), 1 - eps)
    list(surv = as.data.frame(surv_mat), hazard = as.data.frame(hazard_mat))
  }

  # Model names expected by code
  param_names  <- c("expo","weib","gom","lnorm","llog","gam","ggam")
  # use package naming convention for spline models
  spline_names <- c("spl_hazard_1","spl_odds_1")
  cure_names   <- c("weib_mixture_1","lnorm_nonmix_1")

  surv_pred <- list(model = list())
  for (nm in param_names) surv_pred$model[[nm]] <- make_model_matrices(ngroups, time_pred)
  surv_pred$model$spline <- list()
  for (nm in spline_names) surv_pred$model$spline[[nm]] <- make_model_matrices(ngroups, time_pred)
  surv_pred$model$cure <- list()
  for (nm in cure_names) surv_pred$model$cure[[nm]] <- make_model_matrices(ngroups, time_pred)

  surv_model <- list(
    param_models = stats::setNames(vector("list", length(param_names)), param_names),
    spline_models = stats::setNames(vector("list", length(spline_names)), spline_names),
    cure_models = stats::setNames(vector("list", length(cure_names)), cure_names)
  )
  for (nm in spline_names) surv_model$spline_models[[nm]] <- list(knots = c(log(1), log(3), log(10)))

  # smoothed hazard base data (one series per group expected by plotting code)
  haz_df <- data.frame(est.grid = time_pred, haz.est = seq(0.1, 1.1, length.out = length(time_pred)))

  # --- Create a real survfit object but ensure no survival = 1 or 0 that breaks diagnostics ---
  group_levels <- paste0("G", seq_len(ngroups))
  df_fake <- data.frame(
    time = rep(time_pred, times = ngroups),
    status = rep(1L, length(time_pred) * ngroups),
    group = factor(rep(group_levels, each = length(time_pred)), levels = group_levels)
  )

  km_fit <- survival::survfit(survival::Surv(time, status) ~ group, data = df_fake,
                              conf.type = "log", conf.int = 0.95)
  # ensure maxtime present and > 0
  km_fit$maxtime <- max(time_pred)
  # Clip km survival values into (eps, 1 - eps) to avoid log(-log(1)) and similar issues
  if (!is.null(km_fit$surv)) {
    km_fit$surv <- pmin(pmax(km_fit$surv, eps), 1 - eps)
  }
  if (!is.null(km_fit$lower)) km_fit$lower <- pmin(pmax(km_fit$lower, eps), 1 - eps)
  if (!is.null(km_fit$upper)) km_fit$upper <- pmin(pmax(km_fit$upper, eps), 1 - eps)

  # km_names mapping used by some helpers (strata lengths are safe because we used df_fake)
  strata_lengths <- if (!is.null(attr(km_fit, "strata"))) as.integer(attr(km_fit, "strata")) else rep(length(time_pred), ngroups)
  km_names <- rep(seq_len(length(strata_lengths)), times = strata_lengths)

  surv_obs <- list(
    haz = list(
      max = data.frame(time = max(time_pred), smooth = max(haz_df$haz.est)),
      hazards = list(smooth_gr1 = haz_df, smooth_gr2 = haz_df),
      names = seq_len(ngroups)
    ),
    km = km_fit,
    km_names = km_names,
    tp = list(
      gr_1 = data.frame(time = time_pred[-1], smooth = runif(pred_rows, 0.01, 0.9),
                        smooth_upper = runif(pred_rows, .6, 1 - eps), smooth_lower = runif(pred_rows, eps, .4)),
      gr_2 = data.frame(time = time_pred[-1], smooth = runif(pred_rows, 0.01, 0.9),
                        smooth_upper = runif(pred_rows, .6, 1 - eps), smooth_lower = runif(pred_rows, eps, .4)),
      max = 1
    ),
    tp_gr = list()
  )

  # Build tp_gr as data.frames with many columns to satisfy column-indexed selection.
  n_cols_needed <- 23
  tp_df <- as.data.frame(matrix(runif(pred_rows * (n_cols_needed - 1), min = 0.01, max = 0.99),
                                nrow = pred_rows, ncol = n_cols_needed - 1))
  colnames(tp_df) <- paste0("V", 2:n_cols_needed)
  tp_df <- cbind(Time = time_pred[-1], tp_df)
  surv_pred$tp_gr <- list(gr_1 = tp_df, gr_2 = tp_df)

  PERS <- list(
    input = list(time_pred = time_pred, time_unit = 1, spline_mod = TRUE, cure_mod = TRUE, time_horizon = max(time_pred)),
    misc = list(ngroups = ngroups, group_names = group_levels,
                lbls = paste0("M", seq_along(param_names)),
                lbls_spline = paste0("S", seq_along(spline_names)),
                lbls_cure = paste0("C", seq_along(cure_names)),
                cols_tp = n_cols_needed),
    surv_obs = surv_obs,
    surv_pred = surv_pred,
    surv_model = surv_model
  )
  class(PERS) <- "PERSUADE"

  # Additional defensive clipping: ensure all spline predicted survival columns are in (eps, 1 - eps)
  for (nm in names(PERS$surv_pred$model$spline)) {
    mat <- PERS$surv_pred$model$spline[[nm]]$surv
    if (is.matrix(mat) || is.data.frame(mat)) {
      mat[, -1] <- apply(as.matrix(mat[, -1, drop = FALSE]), 2, function(col) pmin(pmax(col, eps), 1 - eps))
      PERS$surv_pred$model$spline[[nm]]$surv <- as.data.frame(mat)
    }
  }

  # Ensure diagnostic functions won't see time == 0 anywhere important
  expect_true(all(PERS$input$time_pred > 0))
  expect_true(all(sapply(PERS$surv_pred$model$spline, function(x) all(x$surv[,1] > 0))))

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
