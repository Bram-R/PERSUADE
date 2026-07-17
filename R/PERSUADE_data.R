#### Bundled example data ----
#' Simulated Three-Arm Trial Survival Data
#'
#' A simulated time-to-event dataset (e.g. representing overall survival in a
#' 3-arm oncology trial) bundled with PERSUADE so the package - and the
#' interactive dashboard launched via [f_run_dashboard()] - can be explored
#' without requiring the user's own data. Event times are Weibull-distributed
#' and differ by arm; censoring combines administrative censoring at 5 years
#' of follow-up with random dropout, so every arm contains both events and
#' censored observations.
#'
#' @format A data frame with 450 rows and 3 columns:
#' \describe{
#'   \item{years}{Numeric. Time to event or censoring, in years.}
#'   \item{status}{Numeric. Event indicator (`1` = event, `0` = censored).}
#'   \item{group}{Factor with 3 levels, `"Control"`, `"Low dose"`,
#'     `"High dose"`, identifying the trial arm.}
#' }
#'
#' @source Simulated data; see `data-raw/create_persuade_trial_data.R` for
#'   the generating code.
#'
#' @examples
#' data(persuade_trial_data)
#' head(persuade_trial_data)
#'
#' \donttest{
#' f_PERSUADE(
#'   name = "dashboard_example",
#'   years = persuade_trial_data$years,
#'   status = persuade_trial_data$status,
#'   group = persuade_trial_data$group,
#'   strata = TRUE,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE,
#'   time_unit = 1 / 12,
#'   time_horizon = 10,
#'   time_pred_surv_table = c(0, 1, 5, 10)
#' )
#' }
"persuade_trial_data"
