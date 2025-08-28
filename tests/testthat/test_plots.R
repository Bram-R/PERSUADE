test_that("plot functions run without error", {
  skip_on_cran()

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
