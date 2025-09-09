library(testthat)

# Basic package load test:
#
# These tests include:
# 1) Verification that the PERSUADE package is correctly installed and
#    discoverable via .packages(all.available = TRUE).
# 2) Provides an early-fail check before running any functional or
#    integration-style tests.

test_that("Package loads correctly", {
  expect_true("PERSUADE" %in% (.packages(all.available = TRUE)))
})
