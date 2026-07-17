library(testthat)

# Regression test for the fix where a single spline model variant failing
# to fit (e.g. flexsurv::flexsurvspline(..., k = 1, scale = "normal")
# throwing "non-finite finite-difference value") used to abort the entire
# f_surv_model() / f_PERSUADE() call. It should now be skipped (stored as
# NULL, with NA AIC/BIC) while the other 8 spline variants (and all
# standard models) still fit and are returned normally.

test_that("f_surv_model() tolerates a single spline variant failing without aborting", {
  skip_if_not_installed("flexsurv")
  skip_on_cran()

  years <- survival::lung$time
  status <- survival::lung$status
  group <- factor(survival::lung$sex)
  form <- stats::as.formula(survival::Surv(years, status) ~ group)

  real_flexsurvspline <- flexsurv::flexsurvspline
  local_mocked_bindings(
    flexsurvspline = function(formula, k, scale, ...) {
      if (identical(k, 1) && identical(scale, "normal")) {
        stop("simulated fitting failure")
      }
      real_flexsurvspline(formula, k = k, scale = scale, ...)
    },
    .package = "flexsurv"
  )

  result <- suppressWarnings(f_surv_model(
    years = years, status = status, group = group, strata = FALSE,
    ngroups = nlevels(group), form = form, spline_mod = TRUE, cure_mod = FALSE,
    cure_link = "logistic", group_names = levels(group)
  ))

  other_models <- result$spline_models[setdiff(names(result$spline_models), "spl_normal_1")]

  expect_null(result$spline_models[["spl_normal_1"]])
  expect_true(all(!vapply(other_models, is.null, logical(1))))

  expect_true(is.na(result$spline_ic$AIC[result$spline_ic$Model == "Normal 1 knots"]))
  expect_true(all(!is.na(result$spline_ic$AIC[result$spline_ic$Model != "Normal 1 knots"])))
})

test_that("f_surv_model() warns (not errors) when a spline variant fails", {
  skip_if_not_installed("flexsurv")
  skip_on_cran()

  years <- survival::lung$time
  status <- survival::lung$status
  group <- factor(survival::lung$sex)
  form <- stats::as.formula(survival::Surv(years, status) ~ group)

  local_mocked_bindings(
    flexsurvspline = function(...) stop("simulated fitting failure"),
    .package = "flexsurv"
  )

  expect_warning(
    f_surv_model(
      years = years, status = status, group = group, strata = FALSE,
      ngroups = nlevels(group), form = form, spline_mod = TRUE, cure_mod = FALSE,
      cure_link = "logistic", group_names = levels(group)
    ),
    "failed to fit and was skipped"
  )
})
