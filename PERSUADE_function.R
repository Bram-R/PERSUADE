#### Main PERSUADE function ----
f_PERSUADE <- function(name = "no_name", years, status, group, 
                       strata = FALSE, spline_mod = FALSE, cure_mod = FALSE, 
                       cure_link = "logistic", time_unit, time_horizon, 
                       time_pred_surv_table) {
  #' Main PERSUADE Function
  #'
  #' Executes the PERSUADE workflow for parametric survival analysis, producing outputs for visualization,
  #' prediction, and Excel export.
  #'
  #' @param name Character, name identifier for the analysis (default: "no_name").
  #' @param years Numeric vector of time-to-event data.
  #' @param status Numeric vector indicating event occurrence (1 for event, 0 for censoring).
  #' @param group Factor indicating group membership.
  #' @param strata Logical, whether to stratify models by group.
  #' @param spline_mod Logical, whether spline models are included.
  #' @param cure_mod Logical, whether cure models are included.
  #' @param cure_link Character string specifying the link function for cure models (default: "logistic").  Link options: "logistic", "loglog", "identity", "probit"
  #' @param time_unit Numeric, the unit of time for annualization.
  #' @param time_horizon Numeric, the time horizon for predictions.
  #' @param time_pred_surv_table Numeric vector, time points for survival table predictions.
  #'
  #' @return A list of results including survival observations, model predictions, and Excel-ready data.
  #'
  #' @examples
  #' \dontrun{
  # f_PERSUADE(
  #   name = "Example",
  #   years = bc$recyrs,
  #   status = bc$censrec,
  #   group = factor(bc$group),
  #   strata = FALSE,
  #   spline_mod = TRUE,
  #   cure_mod = FALSE,
  #   cure_link = "logistic",
  #   time_unit = 1,
  #   time_horizon = 10,
  #   time_pred_surv_table = seq(1, 10, 1)
  # )
  #' }
  #' @export
  # Validate inputs
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
  form <- if (ngroups == 1) as.formula(Surv(years, status) ~ 1) else as.formula(Surv(years, status) ~ group)
  
  # Step 1: Observed data
  km <- npsurv(form)
  km_names <- if (ngroups == 1) {
    rep(1, length(km$time))
  } else {
    unlist(lapply(seq_len(ngroups), function(i) rep(i, km$strata[i])))
  }
  
  haz <- f_hazard(years, status, group, ngroups)
  cum_haz <- f_cum_hazard(years, status, group, ngroups, time_pred, time_unit)
  tp <- f_tp(ngroups, cum_haz, time_unit)
  cox_reg <- coxph(form)
  
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

f_hazard <- function(years, status, group, ngroups) {
  #' Calculate Smoothed Hazard Estimates for Groups
  #'
  #' This function computes smoothed hazard estimates for up to three groups 
  #' using the `muhaz` package and returns a consolidated list of results.
  #'
  #' @param years Numeric vector of times to event or censoring.
  #' @param status Numeric vector indicating event occurrence (1 for event, 0 for censoring).
  #' @param group Factor indicating the group each observation belongs to.
  #' @param ngroups Integer specifying the number of groups to analyze (up to 3).
  #'
  #' @return A list containing smoothed hazard estimates for each group, 
  #'         a vector of group identifiers, and the maximum observed hazard values.
  #'         The list includes:
  #'         - `smooth_grX`: Smoothed hazard results for each group.
  #'         - `names`: Vector indicating group membership of hazard values.
  #'         - `max`: Data frame with the maximum time and hazard values.
  #'
  #' @importFrom muhaz muhaz
  #' @examples
  #' \dontrun{
  #' f_hazard(years = bc$recyrs, status = bc$censrec,
  #'          group = factor(bc$group), ngroups = nlevels(bc$group))
  #' }
  #' @export
  
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

f_cum_hazard <- function(years, status, group, ngroups, time_pred, time_unit) {
  #' Calculate Cumulative Hazard Estimates for Groups
  #'
  #' This function computes cumulative hazard estimates for up to three groups
  #' along with variance and confidence intervals using the `estimateNAH` package.
  #'
  #' @param years Numeric vector of times to event or censoring.
  #' @param status Numeric vector indicating event occurrence (1 for event, 0 for censoring).
  #' @param group Factor indicating the group each observation belongs to.
  #' @param ngroups Integer specifying the number of groups to analyze (up to 3).
  #' @param time_pred Numeric vector of prediction times for cumulative hazard calculation.
  #' @param time_unit Numeric vector specifying the time unit length.
  #'
  #' @return A data frame containing:
  #'         - `group`: Group identifier for each observation.
  #'         - `time`: Prediction times.
  #'         - `H`: Cumulative hazard values.
  #'         - `var`: Variance of cumulative hazard.
  #'         - `H_upper`: Upper bound of the confidence interval.
  #'         - `H_lower`: Lower bound of the confidence interval.
  #'         - `H_delta`: Differences in cumulative hazard values (lagged).
  #'         - `H_upper_delta`: Differences in upper confidence bound (lagged).
  #'         - `H_lower_delta`: Differences in lower confidence bound (lagged).
  #'
  #' @importFrom estimateNAH estimateNAH
  #' @importFrom stats qnorm
  #' @examples
  #' \dontrun{
  #' f_cum_hazard(years = bc$recyrs, status = bc$censrec,
  #'              group =  factor(bc$group), ngroups = nlevels(bc$group),
  #'              time_pred = seq(0, 5, by = 1/12), time_unit = 1/12)
  #' }
  #' @export
  
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
    estimate <- estimateNAH(years[group == levels(group)[gr_idx]], 
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

f_tp <- function(ngroups, cum_haz, time_unit) {
  #' Calculate Annual Transition Probabilities from Observed Data
  #'
  #' This function calculates annual transition probabilities and their confidence intervals 
  #' for up to three groups based on cumulative hazard data. The probabilities are smoothed 
  #' using LOESS regression for each group.
  #'
  #' @param ngroups Integer specifying the number of groups to analyze (up to 3).
  #' @param cum_haz Data frame containing cumulative hazard data. Must include columns:
  #'        `group`, `time`, `H_delta`, `H_upper_delta`, `H_lower_delta`.
  #' @param time_unit Numeric specifying the time unit for annualization (e.g., 1 for years).
  #'
  #' @return A list containing:
  #'         - `gr_1`, `gr_2`, `gr_3`: Data frames for each group with columns:
  #'           `time`, `smooth` (smoothed transition probability),
  #'           `smooth_lower` (lower confidence bound),
  #'           `smooth_upper` (upper confidence bound).
  #'         - `max`: Maximum smoothed upper bound across all groups.
  #'
  #' @importFrom stats loess
  #' @examples
  #' \dontrun{
  #' f_tp(ngroups = nlevels(bc$group), cum_haz = f_cum_hazard() output, time_unit = 1/12)
  #' }
  #' @export
  
  if (!is.integer(ngroups) || ngroups < 1 || ngroups > 3) {
    stop("`ngroups` must be an integer between 1 and 3.")
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

f_surv_model <- function(years, status, group, strata, ngroups, form, spline_mod, cure_mod, cure_link, group_names) {
  #' Fit Parametric Survival Models
  #'
  #' This function fits a variety of parametric survival models, including spline and cure models,
  #' and computes AIC and BIC for model comparison.
  #'
  #' @param years Numeric vector of survival times.
  #' @param status Numeric vector indicating event occurrence (1 for event, 0 for censoring).
  #' @param group Factor indicating the group each observation belongs to.
  #' @param strata Logical indicating whether to use stratified models.
  #' @param ngroups Integer specifying the number of groups.
  #' @param form Formula specifying the survival model.
  #' @param spline_mod Logical, whether to fit spline models.
  #' @param cure_mod Logical, whether to fit cure models.
  #' @param cure_link Character string specifying the link function for cure models.
  #' @param group_names Character vector of group names for cure fraction output.
  #'
  #' @return A list containing fitted models, information criteria (AIC/BIC), and cure fractions.
  #'
  #' @importFrom flexsurv flexsurvreg flexsurvspline flexsurvcure
  #' @importFrom stats BIC
  #' @examples
  #' \dontrun{
  #' f_surv_model(years = bc$recyrs, status = bc$censrec,
  #'              group = factor(bc$group), strata = FALSE,
  #'              ngroups = nlevels(bc$group), form = as.formula(Surv(bc$recyrs, bc$censrec) ~ factor(bc$group)),
  #'              spline_mod = TRUE, cure_mod = FALSE,
  #'              cure_link = "logistic", group_names = levels(bc$group))
  #' }
  #' @export
  
  if (!is.numeric(years) || !is.numeric(status)) {
    stop("`years` and `status` must be numeric vectors.")
  }
  if (!is.factor(group)) {
    stop("`group` must be a factor.")
  }
  if (!is.integer(ngroups) || ngroups < 1) {
    stop("`ngroups` must be a positive integer.")
  }
  
  # Helper function to fit parametric models
  fit_parametric_model <- function(dist, anc = NULL) {
    flexsurvreg(form, dist = dist, anc = anc)
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
        spline_models[[paste0("spl_", scale, "_", k)]] <- flexsurvspline(
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
          flexsurvcure(
            Surv(years, status) ~ 1,
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

f_surv_model_pred <- function(ngroups, time_pred, surv_model, spline_mod, cure_mod, group_names) {
  #' Predict Values from Parametric Survival Models (Including Cure Models)
  #'
  #' Generates predictions (survival and hazard) for parametric survival models,
  #' spline models, and cure models.
  #'
  #' @param ngroups Integer specifying the number of groups.
  #' @param time_pred Numeric vector of prediction times.
  #' @param surv_model List of fitted survival models.
  #' @param spline_mod Logical, whether spline models were fitted.
  #' @param cure_mod Logical, whether cure models were fitted.
  #' @param group_names Character vector of group names.
  #'
  #' @return A list of predicted survival probabilities and hazards for all models.
  #'
  #' @importFrom stats summary
  #' @examples
  #' \dontrun{
  #' f_surv_model_pred(
  #'   ngroups = nlevels(bc$group),
  #'   time_pred = seq(0, 5, 1/12),
  #'   surv_model = f_surv_model() output,
  #'   spline_mod = TRUE,
  #'   cure_mod = FALSE,
  #'   group_names = levels(bc$group)
  #' )
  #' }
  #' @export
  
  if (!is.numeric(time_pred)) stop("`time_pred` must be a numeric vector.")
  if (!is.list(surv_model)) stop("`surv_model` must be a list of fitted models.")
  if (!is.character(group_names)) stop("`group_names` must be a character vector.")
  
  # Helper to get predictions
  extract_predictions <- function(model, time, ngroups, type = "survival") {
    est <- summary(model, t = time, type = type)
    if (ngroups == 1) names(est) <- paste("group=", group_names, sep = "")
    cbind(time, sapply(seq_len(ngroups), function(x) est[[paste("group=", group_names[x], sep = "")]]$est))
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

f_surv_model_pred_gr <- function(ngroups, surv_model, surv_model_pred, spline_mod, cure_mod) {
  #' Grouped Predictions from Parametric Survival Models
  #'
  #' Organizes predictions for each group into a consolidated data frame.
  #'
  #' @param ngroups Integer specifying the number of groups.
  #' @param surv_model List of survival models.
  #' @param surv_model_pred List of predictions from `f_surv_model_pred`.
  #' @param spline_mod Logical, whether spline models were fitted.
  #' @param cure_mod Logical, whether cure models were fitted.
  #'
  #' @return A list of grouped predictions.
  #'
  #' @examples
  #' \dontrun{
  #' f_surv_model_pred_gr(ngroups = nlevels(bc$group), surv_model = f_surv_model() output,
  #'                      surv_model_pred = f_surv_model_pred() output, spline_mod = TRUE, cure_mod = FALSE)
  #' }
  #' @export
  
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

f_surv_model_pred_tp_gr <- function(ngroups, time_pred, time_unit, surv_model_pred_gr, cols_tp) {
  #' Compute Annual Transition Probabilities from Survival Model Predictions
  #'
  #' Computes annual transition probabilities for each group based on survival model predictions.
  #'
  #' @param ngroups Integer specifying the number of groups.
  #' @param time_pred Numeric vector of prediction times.
  #' @param time_unit Numeric, time unit for annualization.
  #' @param surv_model_pred_gr List of grouped survival model predictions.
  #' @param cols_tp Integer, number of transition probability columns.
  #'
  #' @return A list of grouped transition probabilities.
  #'
  #' @importFrom data.table shift
  #' @examples
  #' \dontrun{
  #' f_surv_model_pred_tp_gr(
  #'   ngroups = nlevels(bc$group),
  #'   time_pred = seq(0, 5, 1/12),
  #'   time_unit = 1/12,
  #'   surv_model_pred_gr = surv_model_pred_gr() output,
  #'   cols_tp = 10
  #' )
  #' }
  #' @export
  
  if (!is.numeric(time_pred)) stop("`time_pred` must be numeric.")
  if (!is.list(surv_model_pred_gr)) stop("`surv_model_pred_gr` must be a list.")
  
  # Helper function to compute transition probabilities for a single group
  compute_tp <- function(group) {
    # Convert to data.table for proper handling with shift
    dt <- as.data.table(group)
    tp <- 1 - (1 - (shift(dt[, 2:cols_tp], 1L, type = "lag") - dt[, 2:cols_tp]) /
                 shift(dt[, 2:cols_tp], 1L, type = "lag"))^(1 / time_unit)
    tp <- cbind(Time = dt[-1, 1], tp[-1])  # Remove the first row and attach the time column
    
    return(tp)
  }
  
  # Helper function to truncate each group's TP table after first exceedance of 0.95 (to prevent spike artefacts particularly for the Gompertz)
  truncate_after_threshold <- function(groups, threshold = 0.95,
                                       include_exceed = FALSE,   # TRUE to also blank the exceedance point
                                       time_col = 1L) {          # index of time column
    lapply(groups, function(dt) {
      dt <- copy(as.data.table(dt))          # don't modify input by reference
      cols <- setdiff(seq_len(ncol(dt)), time_col)
      
      # find first index > threshold for each column (returns integer or NA)
      first_idx <- dt[, lapply(.SD, function(col) which(col > threshold)[1]), .SDcols = cols]
      first_idx <- as.integer(first_idx)
      
      # blank later values with NA (vectorised over columns, fast with data.table::set)
      for (k in seq_along(cols)) {
        idx <- first_idx[k]
        if (!is.na(idx)) {
          start <- idx + if (include_exceed) 0L else 1L  # "after" the exceedance by default
          if (start <= nrow(dt)) set(dt, i = start:nrow(dt), j = cols[k], value = NA_real_)
        }
      }
      dt
    })
  }
  
  # Compute transition probabilities for each group
  groups <- lapply(seq_len(ngroups), function(i) {
    compute_tp(surv_model_pred_gr[[paste0("gr_", i)]])
  })
  
  # Set all values after first exceedance of 0.95 to NA for each column (to prevent spike artefacts particularly for the Gompertz)
  groups_trunc <- truncate_after_threshold(groups, threshold = 0.95, include_exceed = FALSE)
  
  # Assign names to groups
  names(groups_trunc) <- paste0("gr_", seq_len(ngroups))
  
  return(groups_trunc)
}

f_surv_model_excel <- function(ngroups, strata, surv_model, spline_mod, cure_mod) {
  #' Generate Output Dataframe for MS Excel 
  #'
  #' Creates a structured dataframe of survival model parameters, including knot values
  #' for spline models, to be exported to Excel.
  #'
  #' @param ngroups Integer specifying the number of groups.
  #' @param strata Logical, whether stratified models are used.
  #' @param surv_model List of fitted survival models.
  #' @param spline_mod Logical, whether spline models are included.
  #' @param cure_mod Logical, whether cure models are included.
  #'
  #' @return A dataframe formatted for MS Excel export, including spline knot values.
  #'
  #' @export
  
  # Helper function: Pad covariance matrices for consistent dimensions
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

