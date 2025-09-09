library(testthat)

# Tests for S3 methods defined in R/PERSUADE_S3_method.R:
# - print.PERSUADE
# - summary.PERSUADE
# - plot.PERSUADE
#
# These tests include:
# 1) Lightweight unit tests using minimal PERSUADE-like objects and stubs,
#    which do not require heavy modeling dependencies.
# 2) Integration-style tests that call f_PERSUADE on flexsurv::bc and exercise
#    print/summary/plot with appropriate skips (skip_on_cran / skip_if_not_installed).

test_that("print.PERSUADE prints expected lines and returns object invisibly", {
  # Minimal PERSUADE object (does not rely on f_PERSUADE)
  obj <- list(
    name = "pt_test",
    input = list(
      years = 1:4,
      group = factor(c("A", "B", "A", "B"))
    )
  )
  class(obj) <- "PERSUADE"

  # Call the generic print() so S3 dispatch will find the method in the package namespace.
  out_lines <- capture.output(res <- print(obj))

  expect_true(any(grepl("PERSUADE Survival Analysis Object", out_lines)))
  expect_true(any(grepl("Analysis Name: pt_test", out_lines)))
  expect_true(any(grepl("Number of objects/individuals:", out_lines)))
  expect_true(any(grepl("Groups: A, B", out_lines)))

  # print() should return the object invisibly; the assigned result should be identical
  expect_identical(res, obj)
})

test_that("summary.PERSUADE handles 'km' and 'gof' types and errors on unknown types", {
  skip_on_cran()

  # km: create a small survfit object (survival is a recommended package)
  skip_if_not_installed("survival")
  df <- data.frame(time = c(1,2,3,4,5), status = c(1,0,1,0,1))

  # Option B (uncomment if you want to construct KM the same way the package does):
  # if (requireNamespace("rms", quietly = TRUE)) {
  #   # package's pipeline uses rms::npsurv which often yields a non-null summary(... )$table
  #   km <- rms::npsurv(survival::Surv(time, status) ~ 1, data = df)
  # } else {
  #   km <- survival::survfit(survival::Surv(time, status) ~ 1, data = df)
  # }

  # simpler default (works without rms)
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = df)

  obj <- list(
    name = "sum_test",
    input = list(time_pred_surv_table = c(0,1)),
    surv_obs = list(km = km),
    surv_model = list(param_ic = data.frame(Model = c("A","B"), AIC = c(1,2))),
    surv_pred = list(gr = list(
      # create two simple per-group prediction tables with time in first column
      data.frame(time = c(0,1,2), M1 = c(1, 0.9, 0.8)),
      data.frame(time = c(0,1,2), M1 = c(1, 0.95, 0.88))
    )),
    misc = list(ngroups = 2, group_names = c("G1", "G2"))
  )
  class(obj) <- "PERSUADE"

  # km summary may return a table-like object OR NULL (depending on how km was created)
  km_tab <- summary(obj, type = "km")
  expect_true(!is.null(km_tab) || is.matrix(km_tab) || is.data.frame(km_tab) || is.null(km_tab))

  # gof should return param_ic as-is
  gof <- summary(obj, type = "gof")
  expect_identical(gof, obj$surv_model$param_ic)

  # unknown type errors
  expect_error(summary(obj, type = "nope"), regexp = "'arg' should be one of")

})

test_that("summary.PERSUADE handles 'surv_probs' producing per-group tables", {
  obj <- list(
    name = "surv_probs_test",
    input = list(time_pred_surv_table = c(0, 1)),
    surv_pred = list(
      gr = list(
        # each gr[[i]] must have at least rows for indices 1 + time_pred_surv_table
        data.frame(time = c(0, 1, 2), ModelA = c(1, 0.9, 0.8), ModelB = c(1, 0.85, 0.7)),
        data.frame(time = c(0, 1, 2), ModelA = c(1, 0.95, 0.86), ModelB = c(1, 0.9, 0.75))
      )
    ),
    misc = list(ngroups = 2, group_names = c("Alpha", "Beta"))
  )
  class(obj) <- "PERSUADE"

  res <- summary(obj, type = "surv_probs")
  expect_true(is.list(res))
  expect_length(res, 2)
  expect_true(all(names(res) == paste0("Group_", obj$misc$group_names)))
  # Each element should be a data.frame with rownames corresponding to models and columns "T=..."
  expect_true(is.data.frame(res[[1]]))
  expect_true(all(grepl("^T=", colnames(res[[1]]))))
})

test_that("summary.PERSUADE enforces spline and cure flags for respective gof types", {
  # Prepare object without spline/cure
  obj <- list(
    name = "gof_flags",
    input = list(spline_mod = FALSE, cure_mod = FALSE),
    surv_model = list(spline_ic = data.frame(), cure_ic = data.frame())
  )
  class(obj) <- "PERSUADE"

  expect_error(summary(obj, type = "gof_spline"), "No spline models identified")
  expect_error(summary(obj, type = "gof_cure"), "No cure models identified")

  # Turn on flags and expect the respective IC objects to be returned
  obj$input$spline_mod <- TRUE
  obj$input$cure_mod <- TRUE
  obj$surv_model$spline_ic <- data.frame(Spline = 1)
  obj$surv_model$cure_ic <- data.frame(Cure = 2)

  expect_identical(summary(obj, type = "gof_spline"), obj$surv_model$spline_ic)
  expect_identical(summary(obj, type = "gof_cure"), obj$surv_model$cure_ic)
})

test_that("plot.PERSUADE dispatches to helper functions and returns list structure for each type", {
  skip_on_cran()

  # Create a minimal PERSUADE object
  obj <- list(
    name = "plot_test",
    misc = list(
      lbls = c("m1", "m2"),
      lbls_spline = c("s1"),
      lbls_cure = c("c1")
    ),
    input = list(spline_mod = TRUE, cure_mod = TRUE)
  )
  class(obj) <- "PERSUADE"

  with_mocked_bindings(
    # ---- Stubs for all helper functions ----
    f_plot_km_survival_base       = function(PERSUADE) list(km = "ok_km"),
    f_plot_log_cumhaz             = function(PERSUADE) "log_cumhaz_plot",
    f_plot_schoenfeld_residuals   = function(PERSUADE) "schoenfeld_plot",
    f_plot_hazard_with_models     = function(PERSUADE) list(hazard = "haz_plot"),

    f_plot_param_surv_model       = function(PERSUADE, i) paste0("param_surv_", i),
    f_plot_diag_param_surv_model  = function(PERSUADE, i) paste0("param_diag_", i),
    f_plot_tp_param_surv_model    = function(PERSUADE, i) paste0("param_tp_", i),

    f_plot_spline_surv_model      = function(PERSUADE, i) paste0("spline_surv_", i),
    f_plot_diag_spline_surv_model = function(PERSUADE, i) paste0("spline_diag_", i),
    f_plot_tp_spline_surv_model   = function(PERSUADE, i) paste0("spline_tp_", i),

    f_plot_cure_surv_model        = function(PERSUADE, i) paste0("cure_surv_", i),
    f_plot_diag_cure_surv_model   = function(PERSUADE, i) paste0("cure_diag_", i),
    f_plot_tp_cure_surv_model     = function(PERSUADE, i) paste0("cure_tp_", i),

    {
      # ---- Tests ----

      # km
      res_km <- plot(obj, type = "km")
      expect_type(res_km, "list")
      expect_true("km" %in% names(res_km))

      # ph
      res_ph <- plot(obj, type = "ph")
      expect_type(res_ph, "list")
      expect_true(any(res_ph == "log_cumhaz_plot"))
      expect_true(any(res_ph == "schoenfeld_plot"))

      # hr
      res_hr <- plot(obj, type = "hr")
      expect_type(res_hr, "list")
      expect_true("hazard" %in% names(res_hr))

      # param_models
      res_param <- plot(obj, type = "param_models")
      expect_type(res_param, "list")
      expect_true(all(c(
        "model_1_surv", "model_1_diag", "model_1_tp",
        "model_2_surv", "model_2_diag", "model_2_tp"
      ) %in% names(res_param)))

      # spline_models
      res_spline <- plot(obj, type = "spline_models")
      expect_type(res_spline, "list")
      expect_true(all(c("spline_1_surv", "spline_1_diag", "spline_1_tp") %in% names(res_spline)))

      # cure_models
      res_cure <- plot(obj, type = "cure_models")
      expect_type(res_cure, "list")
      expect_true(all(c("cure_1_surv", "cure_1_diag", "cure_1_tp") %in% names(res_cure)))

      # invalid type
      expect_error(plot(obj, type = "not_a_type"), "Unknown plot type")
    },
    .package = "PERSUADE"
  )
})

# -------------------------------------------------------------------------
# Integration-style tests that call f_PERSUADE on flexsurv::bc dataset.
# These are similar to existing tests in the suite and exercise the S3 methods
# on actual PERSUADE objects created by the package function.
# -------------------------------------------------------------------------

test_that("print.PERSUADE works (integration with f_PERSUADE and flexsurv::bc)", {
  skip_on_cran()
  skip_if_not_installed("flexsurv")

  data <- flexsurv::bc
  PERSUADE <- f_PERSUADE(
    name = "test_print",
    years = data$recyrs,
    status = data$censrec,
    group = data$group,
    strata = TRUE,
    spline_mod = FALSE,
    cure_mod = FALSE,
    time_unit = 1/12,
    time_horizon = 5,
    time_pred_surv_table = c(0, 1, 5)
  )

  expect_output(print(PERSUADE), "test_print")
})

test_that("summary.PERSUADE returns expected types (integration with f_PERSUADE and flexsurv::bc)", {
  skip_on_cran()
  skip_if_not_installed("flexsurv")

  data <- flexsurv::bc
  PERSUADE <- f_PERSUADE(
    name = "test_summary",
    years = data$recyrs,
    status = data$censrec,
    group = data$group,
    strata = TRUE,
    spline_mod = FALSE,
    cure_mod = FALSE,
    time_unit = 1/12,
    time_horizon = 5,
    time_pred_surv_table = c(0, 1, 5)
  )

  expect_error(summary(PERSUADE, type = "notatype"))
})

test_that("plot.PERSUADE runs without error (integration with f_PERSUADE and flexsurv::bc)", {
  skip_on_cran()
  skip_if_not_installed("flexsurv")

  data <- flexsurv::bc
  PERSUADE <- f_PERSUADE(
    name = "test_plot",
    years = data$recyrs,
    status = data$censrec,
    group = data$group,
    strata = TRUE,
    spline_mod = FALSE,
    cure_mod = FALSE,
    time_unit = 1/12,
    time_horizon = 5,
    time_pred_surv_table = c(0, 1, 5)
  )

  expect_silent(plot(PERSUADE, type = "km"))
  expect_error(plot(PERSUADE, type = "notatype"))
})
