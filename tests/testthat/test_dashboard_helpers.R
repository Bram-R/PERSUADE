library(testthat)

# Tests for PERSUADE functions defined in R/PERSUADE_dashboard.R:
#
# These are the pure-logic helpers behind the Shiny dashboard
# (f_run_dashboard()): validating/cleaning user-supplied data, building a
# sensible prediction time grid, listing the built-in example datasets, and
# picking the best-fitting model from a goodness-of-fit table. They are kept
# free of flexsurv/flexsurvcure/muhaz/sft so they can be unit tested without
# actually fitting survival models.

## f_dashboard_validate_data ----

test_that("f_dashboard_validate_data cleans and summarizes valid input", {
  df <- data.frame(
    t = c(1, 2, 3, 4, NA, 6),
    ev = c(1, 0, 1, 1, 1, 0),
    arm = c("A", "A", "B", "B", "A", "B")
  )

  out <- f_dashboard_validate_data(df, years_col = "t", status_col = "ev", group_col = "arm")

  expect_type(out, "list")
  expect_s3_class(out$data, "data.frame")
  expect_equal(nrow(out$data), 5) # the NA row is dropped
  expect_true(all(c("years", "status", "group") %in% names(out$data)))
  expect_s3_class(out$data$group, "factor")
  expect_equal(out$n, 5)
  expect_equal(out$n_events, 3)
  expect_equal(out$n_groups, 2)
  expect_setequal(out$group_names, c("A", "B"))
})

test_that("f_dashboard_validate_data errors when a required column is missing", {
  df <- data.frame(t = 1:5, ev = c(1, 0, 1, 0, 1))
  expect_error(
    f_dashboard_validate_data(df, years_col = "t", status_col = "ev", group_col = "missing_col"),
    "group"
  )
})

test_that("f_dashboard_validate_data errors on non-positive years", {
  df <- data.frame(t = c(1, 0, -2), ev = c(1, 0, 1), arm = c("A", "A", "B"))
  expect_error(
    f_dashboard_validate_data(df, years_col = "t", status_col = "ev", group_col = "arm"),
    "positive"
  )
})

test_that("f_dashboard_validate_data errors on status values outside 0/1", {
  df <- data.frame(t = c(1, 2, 3), ev = c(1, 2, 0), arm = c("A", "A", "B"))
  expect_error(
    f_dashboard_validate_data(df, years_col = "t", status_col = "ev", group_col = "arm"),
    "status"
  )
})

test_that("f_dashboard_validate_data errors when more than 3 groups are present", {
  df <- data.frame(
    t = 1:8,
    ev = rep(c(1, 0), 4),
    arm = c("A", "A", "B", "B", "C", "C", "D", "D")
  )
  expect_error(
    f_dashboard_validate_data(df, years_col = "t", status_col = "ev", group_col = "arm"),
    "1 to 3 groups"
  )
})

test_that("f_dashboard_validate_data drops unused group levels", {
  df <- data.frame(
    t = 1:6,
    ev = c(1, 0, 1, 0, 1, 0),
    arm = factor(c("A", "A", "B", "B", "B", "B"), levels = c("A", "B", "C"))
  )
  out <- f_dashboard_validate_data(df, years_col = "t", status_col = "ev", group_col = "arm")
  expect_equal(out$n_groups, 2)
  expect_setequal(levels(out$data$group), c("A", "B"))
})

## f_dashboard_time_grid ----

test_that("f_dashboard_time_grid returns an increasing grid spanning [0, time_horizon]", {
  grid <- f_dashboard_time_grid(40)
  expect_true(is.numeric(grid))
  expect_equal(grid[1], 0)
  expect_equal(grid[length(grid)], 40)
  expect_true(all(diff(grid) > 0))
})

test_that("f_dashboard_time_grid respects a small time horizon", {
  grid <- f_dashboard_time_grid(2, n_points = 4)
  expect_equal(grid[1], 0)
  expect_equal(grid[length(grid)], 2)
  expect_true(all(grid >= 0 & grid <= 2))
})

test_that("f_dashboard_time_grid rejects invalid input", {
  expect_error(f_dashboard_time_grid(-1))
  expect_error(f_dashboard_time_grid(0))
  expect_error(f_dashboard_time_grid("a"))
})

## f_dashboard_example_datasets ----

test_that("f_dashboard_example_datasets returns a well-formed registry", {
  reg <- f_dashboard_example_datasets()
  expect_type(reg, "list")
  expect_true(length(reg) >= 1)

  for (nm in names(reg)) {
    entry <- reg[[nm]]
    expect_true(all(c("label", "years_col", "status_col", "group_col", "loader") %in% names(entry)))
    expect_true(is.function(entry$loader))
  }

  expect_true("persuade_trial_data" %in% names(reg))
  bundled <- reg[["persuade_trial_data"]]$loader()
  expect_s3_class(bundled, "data.frame")
  expect_true(reg[["persuade_trial_data"]]$years_col %in% names(bundled))
})

## f_dashboard_best_model ----

test_that("f_dashboard_best_model picks the row with the lowest AIC", {
  gof <- data.frame(Model = c("Exponential", "Weibull", "Gompertz"), AIC = c(120, 95, 130), BIC = c(125, 100, 135))
  best <- f_dashboard_best_model(gof)
  expect_equal(best$Model, "Weibull")
  expect_equal(best$AIC, 95)
})

test_that("f_dashboard_best_model errors on a table without an AIC column", {
  expect_error(f_dashboard_best_model(data.frame(Model = "A")), "AIC")
})

## f_dashboard_run_persuade ----

test_that("f_dashboard_run_persuade returns ok = TRUE and the model on success", {
  local_mocked_bindings(f_PERSUADE = function(...) structure(list(name = "stub"), class = "PERSUADE"))
  out <- f_dashboard_run_persuade(name = "x", years = 1:3, status = c(1, 0, 1), group = factor(c("A", "A", "B")))
  expect_true(out$ok)
  expect_s3_class(out$value, "PERSUADE")
  expect_null(out$error)
})

test_that("f_dashboard_run_persuade captures errors instead of throwing", {
  local_mocked_bindings(f_PERSUADE = function(...) stop("simulated model-fitting failure"))
  out <- f_dashboard_run_persuade(name = "x", years = 1:3, status = c(1, 0, 1), group = factor(c("A", "A", "B")))
  expect_false(out$ok)
  expect_null(out$value)
  expect_match(out$error, "simulated model-fitting failure")
})
