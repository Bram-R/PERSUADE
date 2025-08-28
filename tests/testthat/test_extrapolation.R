test_that("Extrapolation outputs look sensible", {
  skip_on_cran()
  data <- flexsurv::bc
  PERSUADE <- f_PERSUADE(
    name = "test_extrap",
    years = data$recyrs,
    status = data$censrec,
    group = data$group,
    strata = TRUE,
    spline_mod = FALSE,
    cure_mod = FALSE,
    time_unit = 1/12,
    time_horizon = 20,
    time_pred_surv_table = c(0, 1, 5, 10, 20)
  )

  # Survival at time 0 should be 1
  surv_tab <- summary(PERSUADE, type = "surv_probs")

  # For each group, check all rows in column "T= 0" are equal to 1
  for (g in surv_tab) {
    expect_true(all(abs(g[["T= 0"]] - 1) < 1e-8))
  }
})
