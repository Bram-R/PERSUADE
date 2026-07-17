#### Generate the bundled `persuade_trial_data` example dataset ----
#### Used as the default dataset in the PERSUADE dashboard (f_run_dashboard())
#### and for examples that do not require the `flexsurv` package.
####
#### Simulates a 3-arm time-to-event trial (e.g. oncology overall survival)
#### as a mixture-cure process: each arm has a genuine long-term-survivor
#### ("cured") fraction that never experiences the event, plus a susceptible
#### fraction with Weibull-distributed event times. Staggered enrollment
#### (accrual) plus a fixed calendar-time data cutoff gives administrative
#### censoring that varies continuously per subject, on top of random
#### dropout censoring.
####
#### Two numerically-motivated design choices, both empirically verified
#### against real fits (not just reasoned about) - see
#### tests/testthat/test_persuade_trial_data.R and MAP.md:
####
#### 1. A real cured fraction (i.e. a Kaplan-Meier plateau above 0) is
####    included: fitting flexsurvcure's cure models to data with no actual
####    plateau is a known source of unstable/boundary optimization.
#### 2. Administrative censoring is staggered (varies per subject via random
####    enrollment time) rather than a single fixed cutoff shared by every
####    subject. An earlier version used one fixed cutoff, which ties a
####    large mass of subjects to an identical censoring time; that
####    destabilized one specific spline model (flexsurvspline(..., k = 1,
####    scale = "normal")) with "non-finite finite-difference value [2]" -
####    reproduced and confirmed fixed by iterating with a live tester who
####    had the full flexsurv/flexsurvcure/muhaz/sft stack installed (this
####    sandbox does not).

set.seed(20240607)

n_per_group <- 200
groups <- c("Control", "Low dose", "High dose")

# Weibull shape/scale (in years, for the *susceptible* fraction) and cure
# fraction per arm: increasing survival benefit and increasing cure rate
weibull_par <- list(
  "Control"   = list(shape = 1.20, scale = 2.2),
  "Low dose"  = list(shape = 1.15, scale = 3.0),
  "High dose" = list(shape = 1.10, scale = 4.0)
)
cure_fraction <- c("Control" = 0.10, "Low dose" = 0.20, "High dose" = 0.32)

accrual_period <- 2   # years over which subjects are enrolled (staggered entry)
study_duration <- 5   # calendar years from first enrollment to data cutoff
dropout_rate <- 0.05  # random dropout hazard (per year)

sim_group <- function(group_name, n) {
  par <- weibull_par[[group_name]]
  cured <- stats::rbinom(n, size = 1, prob = cure_fraction[[group_name]]) == 1

  # Cured subjects never experience the event within (or beyond) the study
  event_time <- ifelse(cured, Inf, stats::rweibull(n, shape = par$shape, scale = par$scale))

  # Staggered enrollment: each subject's administrative censoring time (time
  # from their own enrollment to the shared calendar cutoff) varies
  # continuously instead of every subject sharing one identical cutoff.
  enroll_time <- stats::runif(n, min = 0, max = accrual_period)
  admin_censor_time <- study_duration - enroll_time

  dropout_time <- stats::rexp(n, rate = dropout_rate)
  censor_time <- pmin(admin_censor_time, dropout_time)

  years <- pmin(event_time, censor_time)
  status <- as.numeric(event_time <= censor_time)

  data.frame(
    years = round(years, 3),
    status = status,
    group = group_name,
    stringsAsFactors = FALSE
  )
}

persuade_trial_data <- do.call(
  rbind,
  lapply(groups, sim_group, n = n_per_group)
)
persuade_trial_data$group <- factor(persuade_trial_data$group, levels = groups)
rownames(persuade_trial_data) <- NULL

# Guard: fitting several PERSUADE models needs a positive minimum time
stopifnot(all(persuade_trial_data$years > 0))

save(persuade_trial_data, file = file.path("data", "persuade_trial_data.rda"), compress = "xz")
