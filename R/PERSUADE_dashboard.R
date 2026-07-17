#### PERSUADE interactive dashboard ----
#### Helper functions backing the Shiny dashboard launched by f_run_dashboard().
#### Kept free of flexsurv/flexsurvcure/muhaz/sft so this file's logic can be
#### unit tested independently of the (heavier) model-fitting stack.

#' Validate and Clean Data for the PERSUADE Dashboard
#'
#' Checks that the requested columns exist in `data`, coerces them to the
#' types expected by [f_PERSUADE()], drops rows with missing values in any
#' of the three columns, and enforces the 1-3 group limit required by
#' [f_hazard()] / [f_cum_hazard()].
#'
#' @param data A data frame containing time-to-event data.
#' @param years_col Character. Name of the time-to-event column in `data`.
#' @param status_col Character. Name of the event/censoring indicator column
#'   in `data` (values must be codeable as `0`/`1`).
#' @param group_col Character. Name of the grouping column in `data`.
#'
#' @return A list with elements:
#'   - `data`: cleaned data frame with columns `years`, `status`, `group`.
#'   - `n`: number of retained rows.
#'   - `n_events`: number of events (`status == 1`).
#'   - `n_groups`: number of groups.
#'   - `group_names`: character vector of group levels.
#'
#' @examples
#' df <- data.frame(t = c(1, 2, 3, 4), ev = c(1, 0, 1, 0), arm = c("A", "A", "B", "B"))
#' f_dashboard_validate_data(df, "t", "ev", "arm")
#'
#' @export
f_dashboard_validate_data <- function(data, years_col, status_col, group_col) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.")

  roles <- c(years = years_col, status = status_col, group = group_col)
  missing_roles <- roles[!roles %in% names(data)]
  if (length(missing_roles) > 0) {
    stop(
      "Column(s) not found in data: ",
      paste(sprintf("%s (%s)", missing_roles, names(missing_roles)), collapse = ", ")
    )
  }

  cleaned <- data.frame(
    years = suppressWarnings(as.numeric(data[[years_col]])),
    status = suppressWarnings(as.numeric(data[[status_col]])),
    group = as.character(data[[group_col]]),
    stringsAsFactors = FALSE
  )
  cleaned <- cleaned[stats::complete.cases(cleaned), , drop = FALSE]

  if (nrow(cleaned) == 0) stop("No complete rows remain after removing missing values.")

  if (!all(cleaned$years > 0)) {
    stop("`", years_col, "` must contain only positive time values.")
  }

  if (!all(cleaned$status %in% c(0, 1))) {
    stop("`", status_col, "` (status) must contain only 0 (censored) or 1 (event) values.")
  }

  cleaned$group <- droplevels(as.factor(cleaned$group))
  n_groups <- nlevels(cleaned$group)

  if (n_groups < 1 || n_groups > 3) {
    stop("PERSUADE supports 1 to 3 groups; found ", n_groups, " group(s) in `", group_col, "`.")
  }

  rownames(cleaned) <- NULL

  list(
    data = cleaned,
    n = nrow(cleaned),
    n_events = sum(cleaned$status == 1),
    n_groups = n_groups,
    group_names = levels(cleaned$group)
  )
}

#' Build a Prediction Time Grid for the Dashboard
#'
#' Produces a sensible, increasing set of time points spanning `[0,
#' time_horizon]` for use as `time_pred_surv_table` in [f_PERSUADE()],
#' using [pretty()] to pick round numbers.
#'
#' @param time_horizon Numeric. Positive maximum follow-up/prediction time.
#' @param n_points Integer. Approximate number of interior grid points
#'   (default: 8).
#'
#' @return A numeric vector starting at `0` and ending at `time_horizon`.
#'
#' @examples
#' f_dashboard_time_grid(40)
#'
#' @export
f_dashboard_time_grid <- function(time_horizon, n_points = 8) {
  if (!is.numeric(time_horizon) || length(time_horizon) != 1 || time_horizon <= 0) {
    stop("`time_horizon` must be a single positive numeric value.")
  }
  if (!is.numeric(n_points) || n_points < 2) {
    stop("`n_points` must be numeric and at least 2.")
  }

  pts <- pretty(c(0, time_horizon), n = n_points)
  pts <- pts[pts >= 0 & pts <= time_horizon]
  sort(unique(c(0, pts, time_horizon)))
}

#' Registry of Example Datasets for the PERSUADE Dashboard
#'
#' Lists the datasets offered in the dashboard's "Data source" selector: the
#' bundled [persuade_trial_data] (always available) and, when installed,
#' `flexsurv::bc`. Each entry's `loader` is a zero-argument function so
#' datasets are only materialized when actually selected.
#'
#' @return A named list of entries, each with elements `label`,
#'   `years_col`, `status_col`, `group_col`, and `loader`.
#'
#' @examples
#' reg <- f_dashboard_example_datasets()
#' reg$persuade_trial_data$loader()
#'
#' @export
f_dashboard_example_datasets <- function() {
  registry <- list(
    persuade_trial_data = list(
      label = "Simulated 3-arm trial (bundled with PERSUADE)",
      years_col = "years",
      status_col = "status",
      group_col = "group",
      loader = function() {
        env <- new.env()
        utils::data("persuade_trial_data", package = "PERSUADE", envir = env)
        env[["persuade_trial_data"]]
      }
    )
  )

  if (requireNamespace("flexsurv", quietly = TRUE)) {
    registry$flexsurv_bc <- list(
      label = "Breast cancer recurrence (flexsurv::bc)",
      years_col = "recyrs",
      status_col = "censrec",
      group_col = "group",
      loader = function() {
        env <- new.env()
        utils::data("bc", package = "flexsurv", envir = env)
        env[["bc"]]
      }
    )
  }

  registry
}

#' Identify the Best-Fitting Model from a Goodness-of-Fit Table
#'
#' Returns the row of `gof` with the lowest AIC, as produced by
#' `summary.PERSUADE(object, type = "gof" / "gof_spline" / "gof_cure")`.
#'
#' @param gof A data frame containing at least `Model` and `AIC` columns.
#'
#' @return A one-row data frame: the best-fitting model.
#'
#' @examples
#' gof <- data.frame(Model = c("Exponential", "Weibull"), AIC = c(120, 95))
#' f_dashboard_best_model(gof)
#'
#' @export
f_dashboard_best_model <- function(gof) {
  if (!is.data.frame(gof) || !"AIC" %in% names(gof)) {
    stop("`gof` must be a data frame containing an `AIC` column.")
  }
  gof[which.min(gof$AIC), , drop = FALSE]
}

#' Run f_PERSUADE() Without Throwing (Dashboard-Safe Wrapper)
#'
#' Calls [f_PERSUADE()] and captures any error as a structured result instead
#' of letting it propagate, so the Shiny dashboard can surface a friendly
#' message rather than crashing.
#'
#' @param ... Arguments passed on to [f_PERSUADE()].
#'
#' @return A list with elements `ok` (logical), `value` (the `PERSUADE`
#'   object, or `NULL` on failure), and `error` (character message, or
#'   `NULL` on success).
#'
#' @examples
#' \donttest{
#' f_dashboard_run_persuade(
#'   name = "Example",
#'   years = survival::lung$time,
#'   status = survival::lung$status,
#'   group = factor(survival::lung$sex),
#'   strata = FALSE,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE,
#'   time_unit = 365.25 / 12,
#'   time_horizon = 2000,
#'   time_pred_surv_table = seq(0, 2000, 365.25)
#' )
#' }
#'
#' @export
f_dashboard_run_persuade <- function(...) {
  tryCatch(
    list(ok = TRUE, value = f_PERSUADE(...), error = NULL),
    error = function(e) list(ok = FALSE, value = NULL, error = conditionMessage(e))
  )
}

#' Launch the PERSUADE Shiny Dashboard
#'
#' Starts an interactive Shiny dashboard for running [f_PERSUADE()] and
#' exploring its outputs (Kaplan-Meier, PH diagnostics, hazards, parametric/
#' spline/cure model fits, GOF tables, exports, and PDF reports) without
#' writing any R code. Requires the Suggested packages `shiny`,
#' `shinydashboard`, `DT`, `bslib`, and `shinycssloaders`.
#'
#' @param ... Additional arguments passed on to [shiny::runApp()] (e.g.
#'   `launch.browser`, `port`).
#'
#' @return Does not return; runs the Shiny application (see
#'   [shiny::runApp()]).
#'
#' @examples
#' \dontrun{
#' f_run_dashboard()
#' }
#'
#' @export
f_run_dashboard <- function(...) {
  required_pkgs <- c("shiny", "shinydashboard", "DT", "bslib", "shinycssloaders")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are required for f_run_dashboard() but are not installed: ",
      paste(missing_pkgs, collapse = ", "),
      ".\nPlease install them with:\n  install.packages(c(\"",
      paste(missing_pkgs, collapse = "\", \""), "\"))"
    )
  }

  app_dir <- system.file("shiny-dashboard", package = "PERSUADE")
  if (app_dir == "" || !file.exists(file.path(app_dir, "app.R"))) {
    stop("The PERSUADE dashboard app was not found in the installed package.")
  }

  shiny::runApp(app_dir, ...)
}
