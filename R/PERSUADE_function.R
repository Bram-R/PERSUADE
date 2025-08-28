#### Main PERSUADE function ----
#' Main PERSUADE Function
#'
#' Executes the PERSUADE workflow for parametric survival analysis, including Kaplan–Meier,
#' parametric, spline, and cure models. Produces outputs for visualization, prediction,
#' and Excel export.
#'
#' @param name Character. Name identifier for the analysis (default: "no_name").
#' @param years Numeric vector of time-to-event data.
#' @param status Numeric vector indicating event occurrence (1 = event, 0 = censoring).
#' @param group Factor indicating group membership.
#' @param strata Logical. Whether to stratify models by group.
#' @param spline_mod Logical. Whether spline models should be fitted.
#' @param cure_mod Logical. Whether cure models should be fitted.
#' @param cure_link Character string specifying the link function for cure models
#'   ("logistic", "loglog", "identity", "probit"; default = "logistic").
#' @param time_unit Numeric. The unit of time for annualization.
#' @param time_horizon Numeric. The maximum prediction time horizon.
#' @param time_pred_surv_table Numeric vector of time points for survival table predictions.
#'
#' @return A list of class `"PERSUADE"` containing:
#'   - `input`: Input arguments used in the analysis.
#'   - `surv_obs`: Observed survival results (Kaplan–Meier, hazards, Cox model).
#'   - `surv_model`: Fitted parametric/spline/cure models.
#'   - `surv_pred`: Model predictions.
#'   - `surv_model_excel`: Excel-ready parameter table.
#'   - `misc`: Auxiliary results (labels, number of groups, etc.).
#'
#' @details The workflow proceeds in three main stages:
#'   1. Observed data (Kaplan–Meier, hazards, Cox regression).
#'   2. Parametric, spline, and cure model fitting.
#'   3. Prediction and export of results.
#'
#' @seealso [f_hazard()], [f_cum_hazard()], [f_surv_model()]
#'
#' @examples
#' years <- survival::lung$time
#' status <-  survival::lung$status
#' group <- factor(survival::lung$sex)
#' f_PERSUADE(
#'   name = "Example",
#'   years = years,
#'   status = status,
#'   group = group,
#'   strata = FALSE,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE,
#'   time_unit = 365.25/12,
#'   time_horizon = 5000,
#'   time_pred_surv_table = seq(0, 5000, 100)
#' )
#' @export
# Validate inputs
f_PERSUADE <- function(name = "no_name", years, status, group,
                       strata = FALSE, spline_mod = FALSE, cure_mod = FALSE,
                       cure_link = "logistic", time_unit, time_horizon,
                       time_pred_surv_table) {
  if (!is.numeric(years) || !is.numeric(status)) stop("`years` and `status` must be numeric vectors.")
  if (!is.factor(group)) stop("`group` must be a factor.")
  if (!is.numeric(time_unit) || time_unit <= 0) stop("`time_unit` must be a positive numeric value.")
  if (!is.numeric(time_horizon) || time_horizon <= 0) stop("`time_horizon` must be a positive numeric value.")

  # strata <- TRUE # for validation purposes
  # spline_mod <- TRUE # for validation purposes
  # cure_mod <- TRUE # for validation purposes
  # cure_link = "logistic" # for validation purposes

  # Inputs
  years <- as.numeric(years)
  status <- as.numeric(status)
  group <- droplevels(as.factor(group))
  if (!cure_mod) cure_link <- NA

  # Determine groups and time variables
  ngroups <- nlevels(group)
  group_names <- levels(group)
  if (ngroups == 1) strata <- FALSE
  time_horizon <- max(time_pred_surv_table, time_horizon)
  time_pred_surv_table <- time_pred_surv_table / time_unit
  time_pred <- seq(0, time_horizon, by = time_unit)

  # Define survival formula
  form <- if (ngroups == 1) as.formula(survival::Surv(years, status) ~ 1) else as.formula(survival::Surv(years, status) ~ group)

  # Step 1: Observed data
  km <- rms::npsurv(form)
  km_names <- if (ngroups == 1) {
    rep(1, length(km$time))
  } else {
    unlist(lapply(seq_len(ngroups), function(i) rep(i, km$strata[i])))
  }

  haz <- f_hazard(years, status, group, ngroups)
  cum_haz <- f_cum_hazard(years, status, group, ngroups, time_pred, time_unit)
  tp <- f_tp(ngroups, cum_haz, time_unit)
  cox_reg <- survival::coxph(form)

  # Step 2: Parametric models
  surv_model <- f_surv_model(
    years, status, group, strata, ngroups, form, spline_mod, cure_mod, cure_link, group_names
  )

  surv_model_pred <- f_surv_model_pred(
    ngroups, time_pred, surv_model, spline_mod, cure_mod, group_names
  )

  surv_model_pred_gr <- f_surv_model_pred_gr(
    ngroups, surv_model, surv_model_pred, spline_mod, cure_mod
  )

  cols_tp <- 8 + if (spline_mod) 9 else 0
  cols_tp <- cols_tp + if (cure_mod) 6 else 0

  surv_model_pred_tp_gr <- f_surv_model_pred_tp_gr(
    ngroups, time_pred, time_unit, surv_model_pred_gr, cols_tp
  )

  surv_model_excel <- f_surv_model_excel(
    ngroups, strata, surv_model, spline_mod, cure_mod
  )

  # Output
  output <- list(
    name = name,
    input = list(
      years = years,
      status = status,
      group = group,
      strata = strata,
      spline_mod = spline_mod,
      cure_mod = cure_mod,
      cure_link = cure_link,
      time_unit = time_unit,
      time_horizon = time_horizon,
      time_pred_surv_table = time_pred_surv_table,
      time_pred = time_pred
    ),
    surv_obs = list(
      km = km,
      km_names = km_names,
      cum_haz = cum_haz,
      haz = haz,
      tp = tp,
      cox_reg = cox_reg
    ),
    surv_model = surv_model,
    surv_pred = list(
      model = surv_model_pred,
      gr = surv_model_pred_gr,
      tp_gr = surv_model_pred_tp_gr
    ),
    surv_model_excel = surv_model_excel,
    misc = c(
      list(
        form = form,
        group_names = group_names,
        ngroups = ngroups,
        lbls = surv_model$param_ic$Model
      ),
      if (spline_mod) list(lbls_spline = surv_model$spline_ic$Model),
      if (cure_mod) list(lbls_cure = surv_model$cure_ic$Model),
      list(cols_tp = cols_tp)
    )
  )
  class(output) <- "PERSUADE" # for implementing S3 class
  return(output)
}

#### Functions used in the main PERSUADE function ----
#' Calculate Smoothed Hazard Estimates
#'
#' Computes smoothed hazard estimates for up to three groups using
#' the \pkg{muhaz} package.
#'
#' @inheritParams f_PERSUADE
#' @param ngroups Integer. Number of groups (1–3).
#'
#' @return A list with elements:
#'   - `hazards`: List of hazard objects (one per group).
#'   - `names`: Vector of group identifiers for hazard values.
#'   - `max`: Data frame with maximum time and hazard values.
#'
#' @examples
#' years <- survival::lung$time
#' status <-  survival::lung$status
#' group <- factor(survival::lung$sex)
#' f_hazard(
#'   years = years,
#'   status = status,
#'   group = group,
#'   ngroups = nlevels(group)
#' )
#'
#' @export
f_hazard <- function(years, status, group, ngroups) {
  if (!is.numeric(years) || !is.numeric(status)) {
    stop("`years` and `status` must be numeric vectors.")
  }
  if (!is.factor(group)) {
    stop("`group` must be a factor.")
  }
  if (!is.numeric(ngroups) || ngroups < 1 || ngroups > 3) {
    stop("`ngroups` must be an integer between 1 and 3.")
  }

  # Helper function to compute smoothed hazard for a group
  compute_hazard <- function(gr_idx) {
    if (gr_idx <= ngroups) {
      return(muhaz::muhaz(years, status, group == levels(group)[gr_idx], max.time = max(years)))
    }
    NULL
  }

  # Compute hazards for each group
  hazards <- lapply(1:3, compute_hazard)
  names(hazards) <- paste0("smooth_gr", 1:3)

  # Construct group names and maximum values
  haz_names <- unlist(
    lapply(seq_len(ngroups), function(i) rep(i, length(hazards[[i]]$est.grid))),
    use.names = FALSE
  )
  haz_max <- data.frame(
    time = ceiling(do.call(max, lapply(hazards, function(h) h$est.grid))),
    smooth = ceiling(do.call(max, lapply(hazards, function(h) h$haz.est)))
  )

  list(hazards = hazards, names = haz_names, max = haz_max)
}

#' Calculate Cumulative Hazard Estimates
#'
#' Computes cumulative hazard estimates for up to three groups along with
#' variance and confidence intervals, using the \pkg{estimateNAH} package.
#'
#' @inheritParams f_PERSUADE
#' @param ngroups Integer. Number of groups (1–3).
#' @param time_pred Numeric vector of prediction times.
#' @param time_unit Numeric. Time unit length for scaling.
#'
#' @return A data frame with columns:
#'   - `group`: Group identifier.
#'   - `time`: Prediction times.
#'   - `H`: Cumulative hazard values.
#'   - `var`: Variance estimates.
#'   - `H_upper`, `H_lower`: 95% confidence interval bounds.
#'   - `H_delta`, `H_upper_delta`, `H_lower_delta`: Differences between time steps.
#'
#' @examples
#' years <- survival::lung$time
#' status <-  survival::lung$status
#' group <- factor(survival::lung$sex)
#' f_cum_hazard(
#'   years = years,
#'   status = status,
#'   group = group,
#'   ngroups = nlevels(group),
#'   time_pred = seq(0, 5000, 100),
#'   time_unit = 30
#' )
#'
#' @export
f_cum_hazard <- function(years, status, group, ngroups, time_pred, time_unit) {
  if (!is.numeric(years) || !is.numeric(status)) {
    stop("`years` and `status` must be numeric vectors.")
  }
  if (!is.factor(group)) {
    stop("`group` must be a factor.")
  }
  if (!is.numeric(ngroups) || ngroups < 1 || ngroups > 3) {
    stop("`ngroups` must be an integer between 1 and 3.")
  }
  if (!is.numeric(time_pred) || !is.numeric(time_unit)) {
    stop("`time_pred` and `time_unit` must be a numeric vectors")
  }

  # Helper function to calculate hazard for a specific group
  calculate_group_hazard <- function(gr_idx) {
    valid_times <- time_pred[time_pred <= max(years[group == levels(group)[gr_idx]])]
    estimate <- sft::estimateNAH(years[group == levels(group)[gr_idx]],
                                 status[group == levels(group)[gr_idx]])
    data.frame(
      time = valid_times,
      H = estimate$H(valid_times),
      var = estimate$Var(valid_times)
    )
  }

  # Generate cumulative hazard data for each group
  group_data <- lapply(1:ngroups, calculate_group_hazard)
  group_labels <- unlist(lapply(seq_len(ngroups), function(i) rep(i, nrow(group_data[[i]]))),
                         use.names = FALSE)
  cumulative_data <- do.call(rbind, group_data)
  cumulative_data$group <- group_labels

  # Calculate confidence intervals and deltas
  cumulative_data$H_upper <- pmax(0, cumulative_data$H + sqrt(cumulative_data$var) * qnorm(0.975))
  cumulative_data$H_lower <- pmax(0, cumulative_data$H - sqrt(cumulative_data$var) * qnorm(0.975))

  calculate_deltas <- function(column, group_col) {
    unlist(lapply(split(column, group_col), function(x) c(NA, diff(x))), use.names = FALSE)
  }

  cumulative_data$H_delta <- -calculate_deltas(cumulative_data$H, cumulative_data$group)
  cumulative_data$H_upper_delta <- -calculate_deltas(cumulative_data$H_upper, cumulative_data$group)
  cumulative_data$H_lower_delta <- -calculate_deltas(cumulative_data$H_lower, cumulative_data$group)

  return(cumulative_data)
}

#' Calculate Transition Probabilities
#'
#' Derives annualized transition probabilities (and confidence bounds)
#' from cumulative hazard estimates, smoothed with LOESS.
#'
#' @param ngroups Integer. Number of groups (1–3).
#' @param cum_haz Data frame from [f_cum_hazard()] with columns
#'   `group`, `time`, `H_delta`, `H_upper_delta`, `H_lower_delta`.
#' @param time_unit Numeric. Time unit for annualization.
#'
#' @return A list with:
#'   - `gr_1`, `gr_2`, `gr_3`: Data frames of smoothed probabilities per group.
#'   - `max`: Maximum upper bound across all groups.
#'
#' @examples
#' years <- survival::lung$time
#' status <-  survival::lung$status
#' group <- factor(survival::lung$sex)
#' cum_haz <- f_cum_hazard(
#'   years = years,
#'   status = status,
#'   group = group,
#'   ngroups = nlevels(group),
#'   time_pred = seq(0, 5000, 100),
#'   time_unit = 30
#' )
#' f_tp(ngroups = nlevels(group), cum_haz = cum_haz, time_unit = 30)
#'
#' @export
f_tp <- function(ngroups, cum_haz, time_unit) {
  if (!is.numeric(ngroups) || ngroups < 1 || ngroups > 3) {
    stop("`ngroups` must be numeric between 1 and 3.")
  }
  if (!all(c("group", "time", "H_delta", "H_upper_delta", "H_lower_delta") %in% names(cum_haz))) {
    stop("`cum_haz` must contain columns: group, time, H_delta, H_upper_delta, H_lower_delta.")
  }
  if (!is.numeric(time_unit) || time_unit <= 0) {
    stop("`time_unit` must be a positive numeric value.")
  }

  # Helper function to calculate transition probabilities for a group
  calculate_tp_group <- function(group_id) {
    group_data <- cum_haz[cum_haz$group == group_id, ]
    tp <- 1 - exp(group_data$H_delta / time_unit)
    tp_upper <- 1 - exp(group_data$H_upper_delta / time_unit)
    tp_lower <- 1 - exp(group_data$H_lower_delta / time_unit)

    # Smooth probabilities using LOESS
    smooth_tp <- loess(tp ~ group_data$time)$fitted
    smooth_tp_upper <- loess(tp_upper ~ group_data$time)$fitted
    smooth_tp_lower <- loess(tp_lower ~ group_data$time)$fitted

    data.frame(
      time = group_data$time,
      smooth = c(NA, pmax(0, pmin(1, smooth_tp))),
      smooth_lower = c(NA, pmax(0, pmin(1, smooth_tp_lower))),
      smooth_upper = c(NA, pmax(0, pmin(1, smooth_tp_upper)))
    )
  }

  # Calculate transition probabilities for each group
  results <- list()
  for (i in seq_len(ngroups)) {
    results[[paste0("gr_", i)]] <- calculate_tp_group(i)
  }

  # Find the maximum smoothed upper bound across all groups
  max_upper <- max(unlist(lapply(results, function(gr) gr$smooth_upper)), na.rm = TRUE)
  results$max <- max_upper

  return(results)
}

#' Fit Parametric Survival Models
#'
#' Fits standard parametric models, spline models, and cure models
#' using the \pkg{flexsurv} package.
#'
#' @inheritParams f_PERSUADE
#' @param form A survival model formula (e.g., `Surv(years, status) ~ group`).
#' @param group_names Character vector of group labels (for cure fractions).
#'
#' @return A list containing:
#'   - `param_models`, `param_ic`: Parametric models and information criteria.
#'   - `spline_models`, `spline_ic`: Spline models and IC (if fitted).
#'   - `cure_models`, `cure_ic`: Cure models and IC (if fitted).
#'
#' @details
#' Models fitted include Exponential, Weibull, Gompertz, Log-normal,
#' Log-logistic, Gamma, Generalised Gamma. Optional spline models
#' (1–3 knots, scales: hazard, odds, normal) and cure models
#' (Weibull, Log-normal, Log-logistic with logistic/probit/etc. link).
#'
#' @examples
#' years <- survival::lung$time
#' status <-  survival::lung$status
#' group <- factor(survival::lung$sex)
#' form <- as.formula(survival::Surv(years, status) ~ group)
#' f_surv_model(
#'   years = years,
#'   status = status,
#'   group = group,
#'   strata = FALSE,
#'   ngroups = nlevels(group),
#'   form = form,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE,
#'   cure_link = "logistic",
#'   group_names = levels(group)
#' )
#'
#' @export
f_surv_model <- function(years, status, group, strata, ngroups, form, spline_mod, cure_mod, cure_link, group_names) {
  if (!is.numeric(years) || !is.numeric(status)) {
    stop("`years` and `status` must be numeric vectors.")
  }
  if (!is.factor(group)) {
    stop("`group` must be a factor.")
  }
  if (!is.numeric(ngroups) || ngroups < 1 || ngroups != as.integer(ngroups)) {
    stop("`ngroups` must be a positive integer.")
  }

  # Helper function to fit parametric models
  fit_parametric_model <- function(dist, anc = NULL) {
    flexsurv::flexsurvreg(form, dist = dist, anc = anc)
  }

  # Fit standard parametric models
  param_models <- list(
    expo = fit_parametric_model("exp"),
    weib = if (strata) fit_parametric_model("weibull", anc = list(shape = ~group)) else fit_parametric_model("weibull"),
    gom = if (strata) fit_parametric_model("gompertz", anc = list(shape = ~group)) else fit_parametric_model("gompertz"),
    lnorm = if (strata) fit_parametric_model("lnorm", anc = list(sdlog = ~group)) else fit_parametric_model("lnorm"),
    llog = if (strata) fit_parametric_model("llogis", anc = list(shape = ~group)) else fit_parametric_model("llogis"),
    gam = if (strata) fit_parametric_model("gamma", anc = list(shape = ~group)) else fit_parametric_model("gamma"),
    ggam = if (strata) fit_parametric_model("gengamma", anc = list(sigma = ~group, Q = ~group)) else fit_parametric_model("gengamma")
  )

  # Compile information criteria for parametric models
  param_labels <- c("Exponential", "Weibull", "Gompertz", "Log-normal", "Log-logistic", "Gamma", "Generalised Gamma")
  param_ic <- data.frame(
    Model = param_labels,
    AIC = sapply(param_models, function(model) model$AIC),
    BIC = sapply(param_models, function(model) BIC(model))
  )

  # Fit spline models if requested
  spline_models <- list()
  spline_ic <- NULL
  if (spline_mod) {
    for (k in 1:3) {
      anc_list <- if (strata) {
        if (k >= 3) {
          list(gamma1 = ~group, gamma2 = ~group, gamma3 = ~group, gamma4 = ~group)
        } else if (k >= 2) {
          list(gamma1 = ~group, gamma2 = ~group, gamma3 = ~group)
        } else {
          list(gamma1 = ~group, gamma2 = ~group)
        }
      } else {
        NULL
      }
      for (scale in c("hazard", "odds", "normal")) {
        spline_models[[paste0("spl_", scale, "_", k)]] <- flexsurv::flexsurvspline(
          form, k = k, scale = scale, anc = anc_list
        )
      }
    }
    model_order <- with(expand.grid(k = 1:3, scale = c("hazard", "odds", "normal")),
                        paste0("spl_", scale, "_", k))
    spline_models <- spline_models[model_order]

    spline_labels <- expand.grid(k = 1:3, scale = c("Hazard", "Odds", "Normal"))
    spline_labels <- apply(spline_labels, 1, function(x) paste(x[2], x[1], "knots"))
    spline_ic <- data.frame(
      Model = spline_labels,
      AIC = sapply(spline_models, function(model) model$AIC),
      BIC = sapply(spline_models, function(model) BIC(model))
    )
  }

  # Fit cure models if requested
  cure_models <- list()
  cure_ic <- NULL
  if (cure_mod) {
    cure_dists <- c("weibull", "lnorm", "llogis")
    for (dist in cure_dists) {
      for (mixture in c(TRUE, FALSE)) {
        cure_models[[paste0("cure_", dist, "_", ifelse(mixture, "mix", "nmix"))]] <- lapply(seq_len(ngroups), function(i) {
          flexsurvcure::flexsurvcure(
            survival::Surv(years, status) ~ 1,
            data = data.frame(years, status)[group == levels(group)[i], ],
            link = cure_link,
            dist = dist,
            mixture = mixture
          )
        })
      }
    }
    cure_labels <- expand.grid(type = c("Mixture", "Non-mixture"), dist = c("Weibull", "Log-normal", "Log-logistic"))
    cure_labels <- apply(cure_labels, 1, paste, collapse = " cure ")
    cure_ic <- data.frame(
      Model = cure_labels,
      AIC = sapply(cure_models, function(models) {
        sum(sapply(models, function(model) model$AIC))
      }),
      CureFraction = sapply(cure_models, function(models) {
        cf_strings <- sapply(models, function(model) {
          est <- round(model$res[1, 1] * 100, 0)
          lci <- round(model$res[1, 2] * 100, 0)
          uci <- round(model$res[1, 3] * 100, 0)
          paste0(est, "% (", lci, "%–", uci, "%)")
        })
        paste(cf_strings, collapse = ", ")  # join if multiple per group
      })
    )
  }

  # Combine and return results
  results <- list(
    param_models = param_models,
    param_ic = param_ic
  )
  if (spline_mod) results$spline_models <- spline_models
  if (spline_mod) results$spline_ic <- spline_ic
  if (cure_mod) results$cure_models <- cure_models
  if (cure_mod) results$cure_ic <- cure_ic
  return(results)
}

#' Predict from Survival Models
#'
#' Generates predicted survival and hazard values from fitted parametric,
#' spline, and cure models.
#'
#' @param ngroups Integer. Number of groups.
#' @param time_pred Numeric vector of prediction times.
#' @param surv_model List of fitted survival models from [f_surv_model()].
#' @param spline_mod Logical. Whether spline models were fitted.
#' @param cure_mod Logical. Whether cure models were fitted.
#' @param group_names Character vector of group labels.
#'
#' @return A list of predictions containing:
#'   - `param_models`: Survival & hazard predictions for standard models.
#'   - `spline`: Predictions for spline models (if fitted).
#'   - `cure`: Predictions for cure models (if fitted).
#'
#' @examples
#' years <- survival::lung$time
#' status <-  survival::lung$status
#' group <- factor(survival::lung$sex)
#' form <- as.formula(survival::Surv(years, status) ~ group)
#' surv_model <- f_surv_model(
#'   years = years,
#'   status = status,
#'   group = group,
#'   strata = FALSE,
#'   ngroups = nlevels(group),
#'   form = form,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE,
#'   cure_link = "logistic",
#'   group_names = levels(group)
#' )
#' f_surv_model_pred(
#'   ngroups = nlevels(group),
#'   time_pred = seq(0, 5000, 100),
#'   surv_model = surv_model,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE,
#'   group_names = levels(group)
#' )
#' @export
f_surv_model_pred <- function(ngroups, time_pred, surv_model, spline_mod, cure_mod, group_names) {
  if (!is.numeric(time_pred)) stop("`time_pred` must be a numeric vector.")
  if (!is.list(surv_model)) stop("`surv_model` must be a list of fitted models.")
  if (!is.character(group_names)) stop("`group_names` must be a character vector.")

  # Helper to get predictions
  extract_predictions <- function(model, time, ngroups, type = "survival") {
    est <- summary(model, t = time, type = type)
    if (ngroups == 1) names(est) <- paste("group=", group_names, sep = "")
    cbind(time, sapply(seq_len(ngroups), function(x) est[[paste("group=", group_names[x], sep = "")]]$est))
    # out <- data.frame(time = est[[1]]$time, sapply(est, function(x) x$est))
    # colnames(out) <- c("time", group_names)
  }

  # Parametric models predictions
  param_models <- c("expo", "weib", "gom", "lnorm", "llog", "gam", "ggam")
  predictions <- lapply(param_models, function(model) {
    list(
      surv = extract_predictions(surv_model$param_models[[model]], time_pred, ngroups),
      hazard = extract_predictions(surv_model$param_models[[model]], time_pred, ngroups, type = "hazard")
    )
  })

  names(predictions) <- param_models

  # Spline models
  if (spline_mod) {
    spline_models <- c("spl_hazard_1", "spl_hazard_2", "spl_hazard_3",
                       "spl_odds_1", "spl_odds_2", "spl_odds_3",
                       "spl_normal_1", "spl_normal_2", "spl_normal_3")
    predictions$spline <- lapply(spline_models, function(model) {
      list(
        surv = extract_predictions(surv_model$spline_models[[model]], time_pred, ngroups),
        hazard = extract_predictions(surv_model$spline_models[[model]], time_pred, ngroups, type = "hazard")
      )
    })
    names(predictions$spline) <- spline_models
  }

  # Cure models
  if (cure_mod) {
    cure_models <- c("cure_weibull_mix", "cure_weibull_nmix", "cure_lnorm_mix",
                     "cure_lnorm_nmix", "cure_llogis_mix", "cure_llogis_nmix")
    predictions$cure <- lapply(cure_models, function(model) {
      list(
        surv = cbind(time_pred, sapply(seq_len(ngroups), function(x) {
          data.frame(summary(surv_model$cure_models[[model]][[x]], t = time_pred))[, 2]
        })),
        hazard = cbind(time_pred, sapply(seq_len(ngroups), function(x) {
          data.frame(summary(surv_model$cure_models[[model]][[x]], t = time_pred, type = "hazard"))[, 2]
        }))
      )
    })
    names(predictions$cure) <- cure_models
  }

  return(predictions)
}

#' Group Predictions by Survival Model
#'
#' Consolidates predictions from [f_surv_model_pred()] into
#' group-specific data frames.
#'
#' @param ngroups Integer. Number of groups.
#' @param surv_model List of survival models from [f_surv_model()].
#' @param surv_model_pred List of predictions from [f_surv_model_pred()].
#' @param spline_mod Logical. Whether spline models were fitted.
#' @param cure_mod Logical. Whether cure models were fitted.
#'
#' @return A list of length `ngroups`, each a data frame with columns:
#'   - `time`
#'   - survival predictions for all models (parametric, spline, cure).
#'
#' @examples
#' years <- survival::lung$time
#' status <-  survival::lung$status
#' group <- factor(survival::lung$sex)
#' form <- as.formula(survival::Surv(years, status) ~ group)
#' surv_model <- f_surv_model(
#'   years = years,
#'   status = status,
#'   group = group,
#'   strata = FALSE,
#'   ngroups = nlevels(group),
#'   form = form,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE,
#'   cure_link = "logistic",
#'   group_names = levels(group)
#' )
#' surv_model_pred <- f_surv_model_pred(
#'   ngroups = nlevels(group),
#'   time_pred = seq(0, 5000, 100),
#'   surv_model = surv_model,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE,
#'   group_names = levels(group)
#' )
#' f_surv_model_pred_gr(
#'   ngroups = nlevels(group),
#'   surv_model = surv_model,
#'   surv_model_pred = surv_model_pred,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE
#' )
#'
#' @export
f_surv_model_pred_gr <- function(ngroups, surv_model, surv_model_pred, spline_mod, cure_mod) {
  if (!is.list(surv_model_pred)) stop("`surv_model_pred` must be a list.")

  # Helper to extract group data
  extract_group_data <- function(group_idx) {
    cbind(
      surv_model_pred[[1]][[1]][, 1],
      sapply(1:7, function(x) surv_model_pred[[x]][[1]][, group_idx + 1]),
      if (spline_mod) sapply(1:9, function(x) surv_model_pred$spline[[x]][[1]][, group_idx + 1]),
      if (cure_mod) sapply(1:6, function(x) surv_model_pred$cure[[x]][[1]][, group_idx + 1])
    )
  }

  groups <- lapply(seq_len(ngroups), extract_group_data)
  names(groups) <- paste0("gr_", seq_len(ngroups))

  # Add column labels
  lbls_all <- c("time", surv_model$param_ic$Model,
                if (spline_mod) surv_model$spline_ic$Model,
                if (cure_mod) surv_model$cure_ic$Model)
  lapply(groups, function(gr) {
    colnames(gr) <- lbls_all
    gr
  })
}

#' Compute Transition Probabilities for Survival Model Predictions
#'
#' @param ngroups Integer, number of groups.
#' @param time_pred Numeric vector of prediction times (currently unused).
#' @param time_unit Numeric, time unit for transition probability calculation.
#' @param surv_model_pred_gr List of group predictions; each element coercible to data.table.
#'   Each group's table should have a time column in column 1 and survival-related columns from 2:cols_tp.
#' @param cols_tp Integer, index of the last survival-related column (i.e., use columns 2:cols_tp).
#'
#' @return Named list of data.frames with transition probabilities (truncated after threshold).
#'
#' @examples
#' years <- survival::lung$time
#' status <-  survival::lung$status
#' group <- factor(survival::lung$sex)
#' form <- as.formula(survival::Surv(years, status) ~ group)
#' surv_model <- f_surv_model(
#'   years = years,
#'   status = status,
#'   group = group,
#'   strata = FALSE,
#'   ngroups = nlevels(group),
#'   form = form,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE,
#'   cure_link = "logistic",
#'   group_names = levels(group)
#' )
#' surv_model_pred <- f_surv_model_pred(
#'   ngroups = nlevels(group),
#'   time_pred = seq(0, 5000, 100),
#'   surv_model = surv_model,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE,
#'   group_names = levels(group)
#' )
#' f_surv_model_pred_gr(
#'   ngroups = nlevels(group),
#'   surv_model = surv_model,
#'   surv_model_pred = surv_model_pred,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE
#' )
#' f_surv_model_pred_tp_gr(
#'   ngroups = nlevels(group),
#'   time_pred = seq(0, 5000, 100),
#'   time_unit = 365.25/12,
#'   surv_model_pred_gr = surv_model_pred_gr,
#'   cols_tp = 8
#')
#'
#' @export
f_surv_model_pred_tp_gr <- function(ngroups, time_pred, time_unit, surv_model_pred_gr, cols_tp) {
  if (!is.numeric(time_pred)) stop("`time_pred` must be numeric.")
  if (!is.list(surv_model_pred_gr)) stop("`surv_model_pred_gr` must be a list.")
  if (!is.numeric(cols_tp) || length(cols_tp) != 1L) stop("`cols_tp` must be a single integer index.")

  # --- Helper: compute transition probabilities for a single group ---
  compute_tp <- function(group) {
    dt <- data.table::as.data.table(group)

    # columns 2:cols_tp are the survival/probability-like columns
    col_idx <- 2:cols_tp
    surv_cols <- dt[, col_idx, with = FALSE]      # keep as data.table for shift()
    shifted   <- data.table::shift(surv_cols, 1L, type = "lag")

    # elementwise formula
    tp <- 1 - (1 - (shifted - surv_cols) / shifted)^(1 / time_unit)

    # drop first row (where shift produced NA)
    tp_df <- as.data.frame(tp[-1, , drop = FALSE])
    Time  <- dt[-1, 1][[1]]

    # make into a proper data.frame
    res <- data.frame(Time = Time, tp_df, check.names = FALSE, row.names = NULL)
    return(res)
  }

  # --- Helper: truncate after first exceedance of threshold in each column ---
  truncate_after_threshold <- function(groups, threshold = 0.95,
                                       include_exceed = FALSE, time_col = 1L) {
    lapply(groups, function(g) {
      dt <- data.table::copy(data.table::as.data.table(g))
      cols <- setdiff(seq_len(ncol(dt)), time_col)  # all non-time columns

      # find first index > threshold per column
      first_idx <- lapply(dt[, cols, with = FALSE], function(col) which(col > threshold)[1])
      first_idx <- as.integer(first_idx)

      # blank later values with NA_real_
      for (k in seq_along(cols)) {
        idx <- first_idx[k]
        if (!is.na(idx)) {
          start <- idx + if (include_exceed) 0L else 1L
          if (start <= nrow(dt)) {
            data.table::set(dt, i = start:nrow(dt), j = cols[k], value = NA_real_)
          }
        }
      }
      return(dt)
    })
  }

  # compute per-group
  groups <- lapply(seq_len(ngroups), function(i) {
    compute_tp(surv_model_pred_gr[[paste0("gr_", i)]])
  })

  # truncate spikes (e.g., Gompertz)
  groups_trunc <- truncate_after_threshold(groups, threshold = 0.95, include_exceed = FALSE)

  names(groups_trunc) <- paste0("gr_", seq_len(ngroups))
  return(groups_trunc)
}

#' Prepare Excel-Ready Survival Model Output
#'
#' Formats model parameters (including spline knots) into a
#' table suitable for export to Excel.
#'
#' @param ngroups Integer. Number of groups.
#' @param strata Logical. Whether stratified models were used.
#' @param surv_model List of fitted models from [f_surv_model()].
#' @param spline_mod Logical. Whether spline models were included.
#' @param cure_mod Logical. Whether cure models were included.
#'
#' @return A transposed data frame containing:
#'   - Distribution names
#'   - Parameter names
#'   - Estimates, SE, CI
#'   - Knot values (if splines fitted)
#'   - Covariance matrix
#'
#' @examples
#' years <- survival::lung$time
#' status <-  survival::lung$status
#' group <- factor(survival::lung$sex)
#' form <- as.formula(survival::Surv(years, status) ~ group)
#' surv_model <- f_surv_model(
#'   years = years,
#'   status = status,
#'   group = group,
#'   strata = FALSE,
#'   ngroups = nlevels(group),
#'   form = form,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE,
#'   cure_link = "logistic",
#'   group_names = levels(group)
#' )
#' f_surv_model_excel(
#'   ngroups = nlevels(group),
#'   strata = FALSE,
#'   surv_model = surv_model,
#'   spline_mod = FALSE,
#'   cure_mod = FALSE
#' )
#'
#' @export
# Helper function: Pad covariance matrices for consistent dimensions
f_surv_model_excel <- function(ngroups, strata, surv_model, spline_mod, cure_mod) {
  pad_matrix <- function(matrix, nrows, ncols) {
    matrix <- as.matrix(matrix)
    if (ncol(matrix) < ncols) {
      matrix <- cbind(matrix, matrix(0, nrow = nrow(matrix), ncol = ncols - ncol(matrix)))
    }
    matrix
  }

  # Collect all model names and categorize them
  param_models <- c("expo", "weib", "gom", "lnorm", "llog", "gam", "ggam")
  spline_models <- c(
    "spl_hazard_1", "spl_odds_1", "spl_normal_1",
    "spl_hazard_2", "spl_odds_2", "spl_normal_2",
    "spl_hazard_3", "spl_odds_3", "spl_normal_3"
  )

  # Combine models based on whether spline_mod is TRUE
  names_models_combined <- c(param_models, if (spline_mod) spline_models)
  surv_models_combined <- c(surv_model$param_models, if (spline_mod) surv_model$spline_models)

  # Prepare data: Distributions and parameters
  distnames <- unlist(lapply(names_models_combined, function(model) {
    if (!is.null(surv_models_combined[[model]])) {
      rep(model, nrow(surv_models_combined[[model]]$res.t))
    } else {
      NULL
    }
  }))

  parnames <- unlist(lapply(names_models_combined, function(model) {
    if (!is.null(surv_models_combined[[model]])) {
      rownames(surv_models_combined[[model]]$res.t)
    } else {
      NULL
    }
  }))

  # Collect model results
  res <- do.call(rbind, lapply(names_models_combined, function(model) {
    if (!is.null(surv_models_combined[[model]])) surv_models_combined[[model]]$res.t else NULL
  }))

  # Determine maximum covariance matrix size
  max_cols <- max(sapply(names_models_combined, function(model) {
    if (!is.null(surv_models_combined[[model]])) ncol(surv_models_combined[[model]]$cov) else 0
  }))

  # Collect and pad covariance matrices
  cov <- do.call(rbind, lapply(names_models_combined, function(model) {
    if (!is.null(surv_models_combined[[model]])) {
      pad_matrix(surv_models_combined[[model]]$cov, nrow(surv_models_combined[[model]]$cov), max_cols)
    } else {
      NULL
    }
  }))

  # Initialize output with blank "Knots" column
  empty <- rep("", nrow(res))
  knots <- empty

  # Add knots for spline models
  if (spline_mod) {
    for (spline_model in spline_models) {
      if (!is.null(surv_model$spline_models[[spline_model]]$knots)) {
        spline_knots <- surv_model$spline_models[[spline_model]]$knots
        for (i in seq_along(spline_knots)) {
          row_idx <- which(
            distnames == spline_model & parnames == paste0("gamma", i - 1)
          )
          if (length(row_idx) > 0) {
            knots[row_idx] <- spline_knots[i]
          }
        }
      }
    }
  }

  # Assemble the final dataframe
  survmod <- cbind(
    Distnames = distnames,
    Parnames = parnames,
    res,
    empty = empty,
    Knots = knots,
    empty = empty,
    Cov_Matrix = cov
  )

  # Rename columns for clarity
  colnames(survmod) <- c(
    "Distnames", "Parnames", colnames(res), "", "Knots", "",
    paste0("Cov_", seq_len(ncol(cov)))
  )

  # Transpose for Excel-friendly format
  return(as.data.frame(t(survmod), stringsAsFactors = FALSE))
}
