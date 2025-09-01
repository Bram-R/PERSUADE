library(testthat)

# Major Update Regression Test
#
# Purpose:
#   This is a long-running, exhaustive integration/regression test suite designed
#   to validate that the PERSUADE package continues to function correctly after
#   major updates. It systematically exercises the main modelling function
#   `f_PERSUADE` across a wide grid of realistic combinations:
#
#     • Grouping: 1-group, 2-groups, 3-groups
#     • Strata: enabled/disabled
#     • Spline-based model: enabled/disabled
#     • Cure model: enabled/disabled, with multiple link functions
#
#   The test ensures:
#     • The package loads correctly.
#     • `f_PERSUADE` runs silently without error across all combinations.
#     • Returned objects have the expected S3 class and internal structure.
#     • Model predictions (survival probabilities) are always within [0, 1].
#     • Reporting functionality integrates cleanly:
#         – First via `f_generate_report()` if the package template is present.
#         – Otherwise via a fallback minimal RMarkdown render to HTML.
#
# Scope:
#   This test is *not* intended for CRAN or normal pull request CI runs because
#   it is computationally expensive and may require LaTeX or rmarkdown setup.
#   It should be executed locally or in nightly/QA regression jobs after major
#   refactors, dependency upgrades, or changes to reporting code.
# ------------------------------------------------------------------------------

skip_on_cran()
skip_if_not_installed("flexsurv")
skip_if_not_installed("rmarkdown")

# --- Exhaustive grid test -----------------------------------------------------
test_that("f_PERSUADE works across 1-3 groups, strata/spline/cure options and integrates with reporting", {

  bc_sample <- flexsurv::bc
  bc_sample$group2 <- factor(ifelse(bc_sample$group == "Good", "Good", "NotGood"))   # Create a 2-group variable

  group_opts <- list(
    "1group" = NULL,
    "2groups" = "group2",   # new 2-level grouping
    "3groups" = "group"     # original 3-level grouping
  )

  strata_opts <- c(FALSE, TRUE)
  spline_opts <- c(FALSE, TRUE)
  cure_opts <- c(FALSE, TRUE)
  cure_links <- c("logistic", "identity", "loglog")

  all_fits <- list()
  runs <- 0L

  for (gname in names(group_opts)) {
    for (strata in strata_opts) {
      for (spline in spline_opts) {
        for (cure in cure_opts) {
          links <- if (cure) cure_links else NA
          for (link in links) {
            runs <- runs + 1L
            msg <- paste("combo", runs, gname,
                         "strata=", strata,
                         "spline=", spline,
                         "cure=", cure,
                         if (!is.na(link)) paste("cure_link=", link) else "")

            yrs <- bc_sample$recyrs
            sts <- bc_sample$censrec
            grp <- if (is.null(group_opts[[gname]])) factor(rep("A", length(yrs)))
            else bc_sample[[group_opts[[gname]]]]
            if (!is.factor(grp)) grp <- factor(grp)

            fit <- NULL
            # Suppress warnings here, but still fail on actual errors
            expect_error(fit <- suppressWarnings(f_PERSUADE(
              name = paste0("int_check_", runs),
              years = yrs,
              status = sts,
              group = grp,
              strata = strata,
              spline_mod = spline,
              cure_mod = cure,
              cure_link = if (cure) link else NA,
              time_unit = 1,
              time_horizon = 100,
              time_pred_surv_table = seq(0, 100, by = 10)
            )), NA)

            expect_s3_class(fit, "PERSUADE")
            expect_true(is.list(fit$input), info = msg)
            expect_true(is.list(fit$misc), info = msg)

            # Check survival predictions are probabilities
            if (!is.null(fit$surv_pred) && length(fit$surv_pred) > 0) {
              model_list <- fit$surv_pred$model
              if (is.list(model_list) && length(model_list) > 0) {
                for (m in names(model_list)) {
                  mobj <- model_list[[m]]
                  if (!is.null(mobj$surv)) {
                    vals <- as.numeric(unlist(mobj$surv[, -1, drop = FALSE]))
                    vals <- vals[is.finite(vals)]
                    expect_true(all(vals >= 0 & vals <= 1),
                                info = paste(msg, "model", m))
                  }
                }
              }
            }

            all_fits[[msg]] <- fit
          }
        }
      }
    }
  }

  # --- Report rendering for all fits ------------------------------------------
  tmpdir <- tempfile("persuade_report_")
  dir.create(tmpdir, recursive = TRUE)

  for (i in seq_along(all_fits)) {
    fit <- all_fits[[i]]
    html_out <- file.path(tmpdir, paste0("integration_report_", i, ".html"))

    expect_silent({
      out_path <- suppressMessages(suppressWarnings(
        f_generate_report(fit)
      ))
      if (is.character(out_path) && nzchar(out_path)) {
        expect_true(file.exists(out_path))
        expect_gt(file.info(out_path)$size, 1000)  # check file is not empty
      }
    })
  }



  last_fit <- tail(all_fits, 1)[[1]]

  f_generate_report(last_fit)

  tmpdir <- tempfile("persuade_report_")
  dir.create(tmpdir, recursive = TRUE)
  html_out <- file.path(tmpdir, "integration_report.html")

  # 1) Try package Rmd template via f_generate_report if available
  pkg_rmd <- system.file("rmd", "PERSUADE_output.Rmd", package = "PERSUADE")
  if (nzchar(pkg_rmd)) {
    expect_silent({
      out_path <- suppressMessages(suppressWarnings(
        f_generate_report(last_fit)
      ))
      if (is.character(out_path) && nzchar(out_path)) {
        expect_true(file.exists(out_path))
      }
    })
  }

  # 2) Fallback minimal Rmd render
  rmd_file <- file.path(tmpdir, "minimal_report.Rmd")
  rmd_contents <- c(
    "---",
    "title: \"Minimal PERSUADE Integration Report\"",
    "output: html_document",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "```",
    "",
    "## PERSUADE object summary",
    "```{r}",
    "print(class(PERSUADE))",
    "print(names(PERSUADE)[1:10])",
    "if (!is.null(PERSUADE$input)) print(str(PERSUADE$input, max.level = 1))",
    "```"
  )
  writeLines(rmd_contents, con = rmd_file)

  expect_silent({
    rmarkdown::render(
      input = rmd_file,
      output_file = html_out,
      envir = list2env(list(PERSUADE = last_fit), parent = globalenv()),
      quiet = TRUE
    )
  })

  last_fit <- tail(all_fits, 1)[[1]]

  f_generate_report(last_fit)

  expect_true(file.exists(html_out))
  expect_gt(file.info(html_out)$size, 1000)
})
