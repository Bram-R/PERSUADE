library(testthat)

# Tests for the bundled example dataset `persuade_trial_data`, used as the
# default dataset in the PERSUADE dashboard (f_run_dashboard()) and in
# examples/vignettes for users who don't have their own survival data yet.

test_that("persuade_trial_data is bundled and well-formed", {
  expect_true(data("persuade_trial_data", package = "PERSUADE", envir = environment())[1] == "persuade_trial_data")

  expect_s3_class(persuade_trial_data, "data.frame")
  expect_named(persuade_trial_data, c("years", "status", "group"))

  expect_true(is.numeric(persuade_trial_data$years))
  expect_true(all(persuade_trial_data$years > 0))

  expect_true(all(persuade_trial_data$status %in% c(0, 1)))

  expect_s3_class(persuade_trial_data$group, "factor")
  expect_true(nlevels(persuade_trial_data$group) >= 2 && nlevels(persuade_trial_data$group) <= 3)

  expect_false(anyNA(persuade_trial_data))
  expect_gte(nrow(persuade_trial_data), 100)

  # Each group should contain both events and censored observations, so the
  # dataset is actually useful for fitting parametric survival models.
  by_group <- split(persuade_trial_data$status, persuade_trial_data$group)
  for (g in by_group) {
    expect_true(any(g == 1))
    expect_true(any(g == 0))
  }
})

test_that("persuade_trial_data runs through the full f_PERSUADE workflow", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("flexsurvcure")
  skip_if_not_installed("muhaz")
  skip_if_not_installed("sft")
  skip_on_cran()

  data("persuade_trial_data", package = "PERSUADE", envir = environment())

  result <- f_PERSUADE(
    name = "dashboard_example",
    years = persuade_trial_data$years,
    status = persuade_trial_data$status,
    group = persuade_trial_data$group,
    strata = TRUE,
    spline_mod = FALSE,
    cure_mod = FALSE,
    time_unit = 1 / 12,
    time_horizon = 10,
    time_pred_surv_table = c(0, 1, 5, 10)
  )

  expect_s3_class(result, "PERSUADE")
})
