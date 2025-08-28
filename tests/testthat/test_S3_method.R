test_that("print.PERSUADE works", {
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

test_that("summary.PERSUADE returns expected types", {
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

test_that("plot.PERSUADE runs without error", {
  skip_on_cran()
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
