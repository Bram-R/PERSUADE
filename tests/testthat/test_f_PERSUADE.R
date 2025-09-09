library(testthat)

# Tests for PERSUADE functions defined in R/PERSUADE_function.R:
#
# These tests include:
# 1) Core functionality checks of f_PERSUADE to ensure it runs successfully
#    on flexsurv::bc and returns a valid PERSUADE S3 object with the expected
#    structure and elements.
# 2) Input validation tests that confirm informative errors are raised when
#    arguments of the wrong type, class, or value are supplied.
# 3) Lightweight runs on small synthetic datasets to verify output consistency
#    without invoking heavy modeling options (e.g., spline_mod, cure_mod).
# 4) Metadata and content checks for top-level list elements, input slots,
#    survival predictions, and group-related attributes.
# 5) Flexibility checks on optional fields (e.g., cure_link) to ensure correct
#    defaults are applied when certain modeling features are disabled.

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

test_that("f_PERSUADE input validation errors for clearly bad inputs", {
  # Use generic expect_error to avoid depending on exact error messages
  years_chr <- as.character(1:10)
  status_num <- as.numeric(0:9 %% 2)
  group_char <- as.character(rep("A", 10))
  group_fac <- factor(rep("A", 10))

  expect_error(
    f_PERSUADE(
      name = "bad_years",
      years = years_chr,
      status = status_num,
      group = group_fac,
      strata = FALSE,
      spline_mod = FALSE,
      cure_mod = FALSE,
      time_unit = 1,
      time_horizon = 10,
      time_pred_surv_table = seq(0, 10, 1)
    )
  )

  expect_error(
    f_PERSUADE(
      name = "bad_group",
      years = as.numeric(1:10),
      status = status_num,
      group = group_char, # not a factor
      strata = FALSE,
      spline_mod = FALSE,
      cure_mod = FALSE,
      time_unit = 1,
      time_horizon = 10,
      time_pred_surv_table = seq(0, 10, 1)
    )
  )

  expect_error(
    f_PERSUADE(
      name = "bad_time_unit",
      years = as.numeric(1:10),
      status = status_num,
      group = group_fac,
      strata = FALSE,
      spline_mod = FALSE,
      cure_mod = FALSE,
      time_unit = 0, # invalid
      time_horizon = 10,
      time_pred_surv_table = seq(0, 10, 1)
    )
  )

  expect_error(
    f_PERSUADE(
      name = "bad_time_horizon",
      years = as.numeric(1:10),
      status = status_num,
      group = group_fac,
      strata = FALSE,
      spline_mod = FALSE,
      cure_mod = FALSE,
      time_unit = 1,
      time_horizon = -10, # invalid
      time_pred_surv_table = seq(0, 10, 1)
    )
  )
})

test_that("f_PERSUADE returns expected structure for single-group minimal run", {
  skip_on_cran()

  # Use a small synthetic dataset; avoid heavy modelling options
  years <- c(5, 10, 15, 20, 25, 30, 35, 40)
  status <- c(1, 1, 0, 1, 0, 1, 0, 1)
  group <- factor(rep("A", length(years)))

  out <- f_PERSUADE(
    name = "single_group_test",
    years = years,
    status = status,
    group = group,
    strata = FALSE,
    spline_mod = FALSE,
    cure_mod = FALSE,
    time_unit = 5,
    time_horizon = 40,
    time_pred_surv_table = seq(0, 40, by = 10)
  )

  expect_s3_class(out, "PERSUADE")

  # top-level names
  expect_true(all(c("name", "input", "surv_obs", "surv_model", "surv_pred", "surv_model_excel", "misc") %in% names(out)))

  # input slot basic checks
  expect_true(is.list(out$input))
  expect_equal(out$name, "single_group_test")
  expect_true(is.numeric(out$input$time_pred))
  expect_true(min(out$input$time_pred) == 0)
  expect_true(max(out$input$time_pred) >= out$input$time_horizon)

  # time_pred_surv_table should be numeric and sorted
  expect_true(is.numeric(out$input$time_pred_surv_table))
  expect_equal(sort(out$input$time_pred_surv_table), out$input$time_pred_surv_table)

  # cure_link should be NA when cure_mod = FALSE (or missing) - accept NA or NULL behavior flexibly
  if (!is.null(out$input$cure_link)) {
    expect_true(is.na(out$input$cure_link) || identical(out$input$cure_link, ""))
  }

  # misc contains group_names and ngroups
  expect_true(is.list(out$misc))
  expect_equal(out$misc$ngroups, 1)
  expect_equal(out$misc$group_names, levels(group))

  # cols_tp is numeric
  expect_true(is.numeric(out$misc$cols_tp))
})
