test_that("Package loads correctly", {
  expect_true("PERSUADE" %in% (.packages(all.available = TRUE)))
})
