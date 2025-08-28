test_that("f_PERSUADE runs and returns expected structure", {
  skip_on_cran()

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
    time_horizon = 10,
    time_pred_surv_table = c(0, 1, 5, 10)
  )

  # Object type and structure
  expect_s3_class(PERSUADE, "PERSUADE")
  expect_type(PERSUADE, "list")

  # Contains expected elements
  expect_true("surv_obs" %in% names(PERSUADE))
  expect_true("surv_model" %in% names(PERSUADE))
  expect_true("misc" %in% names(PERSUADE))

  # Metadata checks
  expect_equal(PERSUADE$name, "test_run")
  expect_true(is.numeric(PERSUADE$misc$ngroups))
})

test_that("f_PERSUADE errors with wrong input types", {
  data <- flexsurv::bc

  expect_error(f_PERSUADE("bad", "not numeric", "wrong", "oops"))
})

