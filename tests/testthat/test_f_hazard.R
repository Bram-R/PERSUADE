test_that("f_hazard works with 2 groups", {
  data <- flexsurv::bc
  haz <- f_hazard(
    years = data$recyrs,
    status = data$censrec,
    group = factor(data$group),
    ngroups = 2
  )

  expect_type(haz, "list")
  expect_true("hazards" %in% names(haz))
  expect_true(all(haz$hazards[[1]]$haz.est >= 0))
})

test_that("f_hazard errors with wrong inputs", {
  data <- flexsurv::bc

  # Non-factor group
  expect_error(f_hazard(data$recyrs, data$censrec, data$group, 2))

  # Too many groups
  expect_error(f_hazard(data$recyrs, data$censrec, factor(data$group), 5))
})
