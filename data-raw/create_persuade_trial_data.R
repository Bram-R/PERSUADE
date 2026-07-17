#### Generate the bundled `persuade_trial_data` example dataset ----
#### Used as the default dataset in the PERSUADE dashboard (f_run_dashboard())
#### and for examples that do not require the `flexsurv` package.
####
#### Simulates a 3-arm time-to-event trial (e.g. oncology overall survival)
#### with Weibull-distributed event times that differ by arm, administrative
#### censoring at the end of follow-up, and random dropout censoring -
#### producing data with both events and censoring in every group, suitable
#### for fitting the standard, spline, and cure models supported by
#### f_PERSUADE().

set.seed(20240607)

n_per_group <- 150
groups <- c("Control", "Low dose", "High dose")

# Weibull shape/scale (in years) per arm: increasing survival benefit
weibull_par <- list(
  "Control"   = list(shape = 1.20, scale = 2.2),
  "Low dose"  = list(shape = 1.15, scale = 3.0),
  "High dose" = list(shape = 1.10, scale = 4.0)
)

admin_censor_time <- 5   # study end (years)
dropout_rate <- 0.05     # random dropout hazard (per year)

sim_group <- function(group_name, n) {
  par <- weibull_par[[group_name]]
  event_time <- stats::rweibull(n, shape = par$shape, scale = par$scale)
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
