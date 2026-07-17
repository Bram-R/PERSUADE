#### Generate the bundled `persuade_trial_data` example dataset ----
#### Used as the default dataset in the PERSUADE dashboard (f_run_dashboard())
#### and for examples that do not require the `flexsurv` package.
####
#### Simulates a 3-arm time-to-event trial (e.g. oncology overall survival)
#### as a mixture-cure process: each arm has a genuine long-term-survivor
#### ("cured") fraction that never experiences the event, plus a susceptible
#### fraction with Weibull-distributed event times. Administrative censoring
#### at the end of follow-up and random dropout censoring are layered on top.
####
#### A real cured fraction (i.e. a Kaplan-Meier plateau above 0) is
#### deliberately included: fitting flexsurvcure's cure models to data with
#### no actual plateau is a known source of unstable/boundary optimization
#### (surfacing as e.g. "non-finite finite-difference value" errors from the
#### numerical Hessian), and a demo dataset for a package built around cure
#### models should showcase them working, not fail on them.

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

admin_censor_time <- 5   # study end (years)
dropout_rate <- 0.05     # random dropout hazard (per year)

sim_group <- function(group_name, n) {
  par <- weibull_par[[group_name]]
  cured <- stats::rbinom(n, size = 1, prob = cure_fraction[[group_name]]) == 1

  # Cured subjects never experience the event within (or beyond) the study
  event_time <- ifelse(cured, Inf, stats::rweibull(n, shape = par$shape, scale = par$scale))

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
