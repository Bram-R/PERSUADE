test_that("f_PERSUADE runs and returns expected structure", {
  skip_on_cran()  # avoid long runs on CRAN

  data <- flexsurv::bc

  PERSUADE <- f_PERSUADE(
    name = "test_run",
    years = data$recyrs,
    status = data$censrec,
    group = data$group,
    strata = TRUE,
    spline_mod = FALSE,
    cure_mod = FALSE,
    time_unit = 1/12,
    time_horizon = 20,
    time_pred_surv_table = c(0, 1, 2, 3, 4, 5, 10, 15, 20)
  )

  # Object type and structure
  expect_s3_class(PERSUADE, "PERSUADE")
  expect_type(PERSUADE, "list")

  # Contains expected elements
  expect_true("surv_obs" %in% names(PERSUADE))
  expect_true("surv_model" %in% names(PERSUADE))
  expect_true("misc" %in% names(PERSUADE))

  # Check misc info
  expect_equal(PERSUADE$name, "test_run")
  expect_true(is.numeric(PERSUADE$misc$ngroups))
})
