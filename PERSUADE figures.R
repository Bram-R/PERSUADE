#### Functions used to create Figures in the R Markdown file ----
f_plot_km_survival <- function(PERSUADE) {
  #' Plot Kaplan-Meier Survival Curves from PERSUADE Object
  #'
  #' Generates Kaplan-Meier survival plots using ggsurvplot, adapting automatically 
  #' to the number of groups in the PERSUADE object.
  #'
  #' @param PERSUADE A PERSUADE object created by `f_PERSUADE()`.
  #'
  #' @return A `ggsurvplot` object showing KM curves with risk table and optional CI/rug.
  #' @export
  #'
  #' @examples
  #' f_plot_km_survival(PERSUADE)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  ngroups <- misc$ngroups
  form <- misc$form
  group_names <- misc$group_names
  
  # Build the data frame used for survival fit
  df_km <- if (ngroups > 1) {
    data.frame(
      years = input$years,
      status = input$status,
      group = input$group
    )
  } else {
    data.frame(
      years = input$years,
      status = input$status
    )
  }
  
  # Fit KM model
  surv_object <- surv_fit(formula = form, data = df_km)
  
  # Plot
  plot_obj <- ggsurvplot(
    fit = surv_object,
    data = df_km,
    legend.labs = group_names,
    risk.table = TRUE,
    conf.int = TRUE,
    surv.median.line = "hv",
    ncensor.plot = ngroups > 1,
    censor = ngroups == 1,
    ggtheme = theme_light(),
    color = if (ngroups > 1) "strata" else "black",
    linetype = if (ngroups > 1) as.integer(c(1, "2222", "5212")),
    size = 0.5,
    legend = "top"
  )
  
  return(plot_obj)
}

f_plot_log_cumhaz <- function(PERSUADE) {
  #' Plot Log(-log(Survival)) vs Log(Time)
  #'
  #' Creates a base R plot for the log cumulative hazard (ln(-ln(S(t)))) against ln(time),
  #' commonly used for checking proportional hazards.
  #'
  #' @param PERSUADE A PERSUADE object created by `f_PERSUADE()`.
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_log_cumhaz(PERSUADE)
  surv_obs <- PERSUADE$surv_obs
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  
  km <- surv_obs$km
  km_names <- surv_obs$km_names
  
  log_time <- log(km$time)
  log_hazard <- log(-log(km$surv))
  
  valid_idx <- is.finite(log_time) & is.finite(log_hazard) & km$surv != 1
  
  x_vals <- log_time[valid_idx]
  y_vals <- log_hazard[valid_idx]
  group_vals <- km_names[valid_idx]
  
  xlim_range <- range(x_vals)
  
  pch_type <- c(1, 8, 9)[group_vals]
  line_color <- c("black", "lightgrey", "darkgrey")[group_vals]
  
  plot(
    x = x_vals, y = y_vals,
    main = "A: LN(cumulative hazard)",
    xlab = "LN(time)",
    ylab = "LN(cumulative hazard)",
    xlim = xlim_range,
    cex = 0.6,
    pch = pch_type,
    col = line_color
  )
  
  if (misc$ngroups > 1) {
    legend(
      "topleft",
      legend = misc$group_names,
      col = c("black", "lightgrey", "darkgrey")[unique(km_names)],
      lty = 1,
      bty = "n"
    )
  }
}

f_plot_schoenfeld_residuals <- function(PERSUADE) {
  #' Plot Schoenfeld Residuals for Proportional Hazards Diagnostics
  #'
  #' Plots scaled Schoenfeld residuals and regression lines from Cox PH model to assess PH assumptions.
  #'
  #' @param PERSUADE A PERSUADE object created by `f_PERSUADE()`.
  #'
  #' @return One or more base R plots.
  #' @export
  #'
  #' @examples
  #' f_plot_schoenfeld_residuals(PERSUADE)
  cox_fit <- PERSUADE$surv_obs$cox_reg
  misc <- PERSUADE$misc
  
  zph <- cox.zph(cox_fit, terms = FALSE)
  zph_data <- zph$y
  zph_time <- zph$x
  ngroups <- misc$ngroups
  
  for (i in 1:(ngroups - 1)) {
    plot(
      zph[i], col = "5", lwd = 2,
      main = paste0(ifelse(ngroups > 2, LETTERS[i + 1], ""), ": Scaled Schoenfeld residuals"),
      xlab = "time",
      ylab = paste("Beta(t) for", misc$group_names[i + 1], "vs", misc$group_names[1])
    )
    
    abline(h = 0, lty = 3)
    
    lm_fit <- lm(zph_data[, i] ~ zph_time)
    abline(lm_fit$coefficients, col = "7", lty = 1, lwd = 2)
    
    legend("bottomleft", legend = c(
      "smoothed line (natural spline with df = 4)",
      "regression line"
    ),
    col = c("5", "7"), lty = 1, cex = 0.8, bty = "n"
    )
  }
}

f_plot_smoothed_hazard <- function(PERSUADE) {
  #' Plot Smoothed Hazard Functions by Group
  #'
  #' Creates a base R line plot showing smoothed hazard functions for all groups in the PERSUADE object.
  #'
  #' @param PERSUADE A PERSUADE object created by `f_PERSUADE()`.
  #'
  #' @return A base R plot.
  #' @export
  #'
  #' @examples
  #' f_plot_smoothed_hazard(PERSUADE)
  surv_obs <- PERSUADE$surv_obs
  misc <- PERSUADE$misc
  ngroups <- misc$ngroups
  
  # Define plotting colors and line types (extendable, will recycle if needed)
  haz_line_color <- c("black", "lightgrey", "darkgrey")
  line_type <- as.integer(c(1, "2222", "5212", "3313", "1144", "42"))
  
  # Extract plot limits
  xlim_range <- c(0, surv_obs$haz$max$time)
  ylim_range <- c(0, surv_obs$haz$max$smooth)
  
  # Loop through groups to extract smoothed hazard data
  for (i in seq_len(ngroups)) {
    group_name <- paste0("smooth_gr", i)
    haz_data <- surv_obs$haz$hazards[[group_name]]
    
    if (i == 1) {
      # Start plot with group 1
      with(haz_data, {
        plot(
          x = est.grid,
          y = haz.est,
          type = "l",
          col = haz_line_color[(i - 1) %% length(haz_line_color) + 1],
          lty = line_type[(i - 1) %% length(line_type) + 1],
          lwd = 2,
          xlab = "time",
          ylab = "smoothed hazard rate",
          xlim = xlim_range,
          ylim = ylim_range
        )
      })
    } else {
      # Add lines for remaining groups
      with(haz_data, {
        lines(
          x = est.grid,
          y = haz.est,
          col = haz_line_color[(i - 1) %% length(haz_line_color) + 1],
          lty = line_type[(i - 1) %% length(line_type) + 1],
          lwd = 2
        )
      })
    }
  }
  
  # Add legend
  legend(
    "topleft",
    legend = misc$group_names,
    col = haz_line_color[seq_len(ngroups)],
    lty = line_type[seq_len(ngroups)],
    cex = 0.8,
    bty = "n"
  )
}

f_plot_hazard_with_models <- function(PERSUADE) {
  #' Plot Observed and Parametric Hazard Estimates by Group
  #'
  #' Plots smoothed hazard estimates from observed data with optional overlays from
  #' parametric, spline-based, and cure models for each group.
  #'
  #' @param PERSUADE A PERSUADE object created by `f_PERSUADE()`.
  #'
  #' @return Base R plots by group with multiple hazard overlays.
  #' @export
  #'
  #' @examples
  #' f_plot_hazard_with_models(PERSUADE)
  input <- PERSUADE$input
  surv_pred <- PERSUADE$surv_pred
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  ngroups <- misc$ngroups
  
  # Define base plot styles
  haz_line_color <- c("black", "lightgrey", "darkgrey")
  line_color <- 1:9
  line_type <- c(1, 2222, 5212)
  
  # Define model groups to iterate over
  model_types <- list(
    base = names(PERSUADE$surv_model$param_models),
    spline = names(PERSUADE$surv_model$spline_models),
    cure = names(PERSUADE$surv_model$cure_models)
  )
  
  # Determine cure column per group
  cure_col_index <- c(2, 3, 4)  # assumes group 1: col 2, group 2: col 3, etc.
  
  for (i in seq_len(ngroups)) {
    group_label <- paste("Group", misc$group_names[i])
    smooth_name <- paste0("smooth_gr", i)
    smooth_data <- surv_obs$haz$hazards[[smooth_name]]
    
    plot(
      smooth_data$est.grid, smooth_data$haz.est,
      main = group_label,
      type = "l",
      col = haz_line_color[i],
      lty = line_type[i],
      lwd = 2,
      xlab = "time",
      ylab = "smoothed hazard rate",
      xlim = c(0, surv_obs$haz$max$time),
      ylim = c(0, surv_obs$haz$max$smooth)
    )
    
    # ── Parametric Model Overlays ─────────────────────────────────────────────
    for (j in seq_along(model_types$base)) {
      lines(
        surv_pred$model[[model_types$base[j]]]$hazard[,1],
        surv_pred$model[[model_types$base[j]]]$hazard[,i + 1],
        col = line_color[j],
        lty = line_type[i],
        lwd = 1
      )
    }
    legend("topleft", legend = misc$lbls, col = line_color[1:7], lty = line_type[i], cex = 0.8, bty = "n")
    
    # ── Spline Models ─────────────────────────────────────────────────────────
    if (isTRUE(input$spline_mod)) {
      plot(
        smooth_data$est.grid, smooth_data$haz.est,
        main = group_label,
        type = "l",
        col = haz_line_color[i],
        lty = line_type[i],
        lwd = 2,
        xlab = "time",
        ylab = "smoothed hazard rate",
        xlim = c(0, surv_obs$haz$max$time),
        ylim = c(0, surv_obs$haz$max$smooth)
      )
      for (j in seq_along(model_types$spline)) {
        lines(
          surv_pred$model$spline[[model_types$spline[j]]]$hazard[,1],
          surv_pred$model$spline[[model_types$spline[j]]]$hazard[,i + 1],
          col = line_color[j],
          lty = line_type[i],
          lwd = 1
        )
      }
      legend("topleft", legend = misc$lbls_spline, col = line_color[1:9], lty = line_type[i], cex = 0.8, bty = "n")
    }
    
    # ── Cure Models ───────────────────────────────────────────────────────────
    if (isTRUE(input$cure_mod)) {
      plot(
        smooth_data$est.grid, smooth_data$haz.est,
        main = group_label,
        type = "l",
        col = haz_line_color[i],
        lty = line_type[i],
        lwd = 2,
        xlab = "time",
        ylab = "smoothed hazard rate",
        xlim = c(0, surv_obs$haz$max$time),
        ylim = c(0, surv_obs$haz$max$smooth)
      )
      for (j in seq_along(model_types$cure)) {
        col_idx <- if (length(cure_col_index) >= i) cure_col_index[i] else 2
        lines(
          surv_pred$model$cure[[model_types$cure[j]]]$hazard[,1],
          surv_pred$model$cure[[model_types$cure[j]]]$hazard[,i + 1],
          col = line_color[j],
          lty = line_type[i],
          lwd = 1
        )
      }
      legend("topleft", legend = misc$lbls_cure, col = line_color[1:6], lty = line_type[i], cex = 0.8, bty = "n")
    }
  }
}

f_plot_param_surv_model <- function(PERSUADE, model_key = "expo", model_label = "Exponential", col_index = 1) {
  #' Plot Kaplan-Meier and Fitted Parametric Survival Model
  #'
  #' Dynamically overlays a parametric survival model on top of KM curves
  #' with shaded confidence bands per group.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_key Character. Name of the model as stored in `surv_pred$model`.
  #' @param model_label Character. Main title of the plot.
  #' @param col_index Integer. Base color index for fitted curves (1–9).
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_param_surv_model(PERSUADE, "weib", "Weibull", 2)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "2222", "5212"))
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Base KM plot
  plot(
    surv_obs$km, lwd = 2, col = km_line_color[1:ngroups],
    main = paste("A: Kaplan-Meier (", model_label, ")", sep = ""),
    xlab = "time", ylab = "survival"
  )
  
  # Add shaded CI per group
  for (i in seq_len(ngroups)) {
    idx <- surv_obs$km_names == i
    polygon(
      na.omit(data.frame(
        x = c(surv_obs$km$time[idx], rev(surv_obs$km$time[idx])),
        y = c(surv_obs$km$lower[idx], rev(surv_obs$km$upper[idx]))
      )),
      col = adjustcolor(km_line_color[i], alpha.f = 0.3),
      border = NA
    )
  }
  
  # Overlay model predictions
  for (i in seq_len(ngroups)) {
    lines(
      x = input$time_pred,
      y = surv_pred$model[[model_key]]$surv[, i + 1],
      col = as.character(col_index),
      lty = line_type[i],
      lwd = 2
    )
  }
  
  # Add legend with model and KM lines
  legend(
    "bottomleft",
    legend = c(group_names, rep("", ngroups)),
    col = c(rep(as.character(col_index), ngroups), km_line_color[1:ngroups]),
    lty = c(line_type[1:ngroups], rep(1, ngroups)),
    cex = 0.8, ncol = 2, bty = "n"
  )
}

f_plot_spline_surv_model <- function(PERSUADE, model_key, model_label, col_index = 1) {
  #' Plot Kaplan-Meier and Fitted Spline-Based Survival Model
  #'
  #' Dynamically overlays a spline-based survival model on top of KM curves
  #' with shaded confidence bands per group, and vertical lines at knot locations.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_key Character. Name of the model as stored in `surv_pred$model` (e.g., "spl_hazard1").
  #' @param model_label Character. Main title of the plot.
  #' @param col_index Integer. Base color index for fitted curves (1–9).
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_spline_surv_model(PERSUADE, "spl_hazard_1", "Spline, 1 knot, hazard scale", 1)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  surv_model <- PERSUADE$surv_model
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "2222", "5212"))
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Base KM plot
  plot(
    surv_obs$km, lwd = 2, col = km_line_color[1:ngroups],
    main = paste("A: Kaplan-Meier (", model_label, ")", sep = ""),
    xlab = "time", ylab = "survival"
  )
  
  # Add shaded CI per group
  for (i in seq_len(ngroups)) {
    idx <- surv_obs$km_names == i
    polygon(
      na.omit(data.frame(
        x = c(surv_obs$km$time[idx], rev(surv_obs$km$time[idx])),
        y = c(surv_obs$km$lower[idx], rev(surv_obs$km$upper[idx]))
      )),
      col = adjustcolor(km_line_color[i], alpha.f = 0.3),
      border = NA
    )
  }
  
  # Overlay model predictions for each group
  for (i in seq_len(ngroups)) {
    lines(
      x = input$time_pred,
      y = surv_pred$model$spline[[model_key]]$surv[, i + 1],
      col = as.character(col_index),
      lty = line_type[i],
      lwd = 2
    )
  }
  
  # Add vertical lines for knots
  knots <- exp(surv_model$spline_models[[model_key]]$knots)
  if (length(knots) >= 2) {
    # Middle knots (internal): solid dashed line
    internal_knots <- knots[2:(length(knots) - 1)]
    abline(v = internal_knots, lty = 2, lwd = 1.5)
    # Boundary knots (first & last): lighter dashed line
    abline(v = c(knots[1], knots[length(knots)]), lty = 3, lwd = 1)
  }
  
  # Add legend
  legend(
    "bottomleft",
    legend = c(group_names, rep("", ngroups)),
    col = c(rep(as.character(col_index), ngroups), km_line_color[1:ngroups]),
    lty = c(line_type[1:ngroups], rep(1, ngroups)),
    cex = 0.8, ncol = 2, bty = "n"
  )
}

f_plot_cure_surv_model <- function(PERSUADE, model_key, model_label, col_index = 1) {
  #' Plot Kaplan-Meier and Fitted Cure-Based Survival Model
  #'
  #' Dynamically overlays a cure-based survival model on top of KM curves
  #' with shaded confidence bands per group.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_key Character. Name of the model as stored in `surv_pred$model$cure`.
  #' @param model_label Character. Main title of the plot.
  #' @param col_index Integer. Base color index for fitted curves (1–9).
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_cure_surv_model(PERSUADE, "cure_weibull_mix", "Weibull mixture cure", 1)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "2222", "5212"))
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Base KM plot
  plot(
    surv_obs$km, lwd = 2, col = km_line_color[1:ngroups],
    main = paste("A: Kaplan-Meier (", model_label, ")", sep = ""),
    xlab = "time", ylab = "survival"
  )
  
  # Add shaded CI per group
  for (i in seq_len(ngroups)) {
    idx <- surv_obs$km_names == i
    polygon(
      na.omit(data.frame(
        x = c(surv_obs$km$time[idx], rev(surv_obs$km$time[idx])),
        y = c(surv_obs$km$lower[idx], rev(surv_obs$km$upper[idx]))
      )),
      col = adjustcolor(km_line_color[i], alpha.f = 0.3),
      border = NA
    )
  }
  
  # Overlay cure model predictions
  for (i in seq_len(ngroups)) {
    lines(
      x = input$time_pred,
      y = surv_pred$model$cure[[model_key]]$surv[, i + 1],
      col = as.character(col_index),
      lty = line_type[i],
      lwd = 2
    )
  }
  
  # Add legend
  legend(
    "bottomleft",
    legend = c(group_names, rep("", ngroups)),
    col = c(rep(as.character(col_index), ngroups), km_line_color[1:ngroups]),
    lty = c(line_type[1:ngroups], rep(1, ngroups)),
    cex = 0.8, ncol = 2, bty = "n"
  )
}

f_plot_diag_param_surv_model <- function(PERSUADE, model_key, model_label, col_index = 1) {
  #' Plot Diagnostic Figures for Standard Parametric Survival Models
  #'
  #' Produces diagnostic plots for various parametric survival models,
  #' using transformations specific to each model family.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_key Character. Name of the model as stored in `surv_pred$model`.
  #' @param model_label Character. Main title of the plot.
  #' @param col_index Integer. Base color index for fitted curves (1–9).
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_diag_param_surv_model(PERSUADE, "weib", "Weibull", 2)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "2222", "5212"))
  point_shape <- c(1, 8, 9)
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Define transformation logic for each model
  transform_funs <- list(
    expo  = list(
      obs  = function(time, surv) cbind(time, -log(surv)),
      pred = function(time, surv) cbind(time, -log(surv))
    ),
    weib  = list(
      obs  = function(time, surv) cbind(log(time), log(-log(surv))),
      pred = function(time, surv) cbind(log(time), log(-log(surv)))
    ),
    gom   = list(
      obs  = function(time, haz) cbind(time, log(haz)),
      pred = function(time, haz) cbind(time, log(haz))
    ),
    lnorm = list(
      obs  = function(time, surv) cbind(log(time), -qnorm(surv)),
      pred = function(time, surv) cbind(log(time), -qnorm(surv))
    ),
    llog  = list(
      obs  = function(time, surv) cbind(log(time), -log(surv / (1 - surv))),
      pred = function(time, surv) cbind(log(time), -log(surv / (1 - surv)))
    ),
    gam   = list(
      obs  = function(time, surv) cbind(log(time), log(-log(surv))),
      pred = function(time, surv) cbind(log(time), log(-log(surv)))
    ),
    ggam  = list(
      obs  = function(time, surv) cbind(log(time), log(-log(surv))),
      pred = function(time, surv) cbind(log(time), log(-log(surv)))
    )
  )
  
  if (!model_key %in% names(transform_funs)) {
    stop("Model key not recognized in transform_funs")
  }
  
  tf <- transform_funs[[model_key]]
  
  # Prepare observed data
  if (model_key == "gom") {
    obs_data <- tf$obs(
      time = c(surv_obs$haz$hazards$smooth_gr1$est.grid,
               if (ngroups > 1) surv_obs$haz$hazards$smooth_gr2$est.grid else NULL,
               if (ngroups > 2) surv_obs$haz$hazards$smooth_gr3$est.grid else NULL),
      haz = c(surv_obs$haz$hazards$smooth_gr1$haz.est,
              if (ngroups > 1) surv_obs$haz$hazards$smooth_gr2$haz.est else NULL,
              if (ngroups > 2) surv_obs$haz$hazards$smooth_gr3$haz.est else NULL)
    )
    obs_col <- km_line_color[surv_obs$haz$names]
    obs_pch <- point_shape[surv_obs$haz$names]
  } else {
    obs_data <- tf$obs(
      time = surv_obs$km$time,
      surv = surv_obs$km$surv
    )
    obs_col <- km_line_color[surv_obs$km_names]
    obs_pch <- point_shape[surv_obs$km_names]
  }
  
  # Base plot
  plot(
    obs_data, cex = 0.6, pch = obs_pch, col = obs_col,
    main = paste("B: Diagnostic plot (", model_label, ")", sep = ""),
    xlab = switch(model_key,
                  expo  = "time",
                  weib  = "LN(time)",
                  gom   = "time",
                  lnorm = "LN(time)",
                  llog  = "LN(time)",
                  gam   = "LN(time)",
                  ggam  = "LN(time)"),
    ylab = switch(model_key,
                  expo  = "cumulative hazard",
                  weib  = "LN(cumulative hazard)",
                  gom   = "LN(smoothed hazard)",
                  lnorm = "- standard normal quantiles",
                  llog  = "-LN(survival odds)",
                  gam   = "LN(cumulative hazard)",
                  ggam  = "LN(cumulative hazard)")
  )
  
  # Add predicted lines
  for (i in seq_len(ngroups)) {
    if (model_key == "gom") {
      pred_data <- tf$pred(
        time = surv_pred$model$gom$hazard[, 1],
        haz = surv_pred$model$gom$hazard[, i + 1]
      )
    } else {
      pred_data <- tf$pred(
        time = input$time_pred,
        surv = surv_pred$model[[model_key]]$surv[, i + 1]
      )
    }
    lines(pred_data, col = as.character(col_index), lty = line_type[i], lwd = 2)
  }
  
  # Legend
  legend(
    "topleft",
    legend = c(group_names, rep("", ngroups)),
    col = c(rep(as.character(col_index), ngroups), km_line_color[1:ngroups]),
    lty = c(line_type[1:ngroups], rep(NA, ngroups)),
    pch = c(rep(NA, ngroups), point_shape[1:ngroups]),
    cex = 0.8, ncol = 2, bty = "n"
  )
}

f_plot_diag_spline_surv_model <- function(PERSUADE, model_key, model_label, col_index = 1) {
  #' Plot Diagnostic Figures for Spline-based Survival Models
  #'
  #' Produces diagnostic plots for spline-based survival models,
  #' with transformations depending on the spline scale (hazard, odds, normal).
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_key Character. Name of the spline model in `surv_pred$model`.
  #' @param model_label Character. Main title of the plot.
  #' @param col_index Integer. Base color index for fitted curves (1–9).
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_diag_spline_surv_model(PERSUADE, "spl_hazard_1", "Spline, 1 knot, hazard scale", 1)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  surv_model <- PERSUADE$surv_model
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "2222", "5212"))
  point_shape <- c(1, 8, 9)
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Determine spline type from model_key
  spline_type <- if (grepl("hazard", model_key)) {
    "hazard"
  } else if (grepl("odds", model_key)) {
    "odds"
  } else if (grepl("normal", model_key)) {
    "normal"
  } else {
    stop("Model key not recognized as spline hazard/odds/normal")
  }
  
  # Define transformations
  transform_funs <- list(
    hazard = list(
      obs  = function(time, surv) cbind(log(time), log(-log(surv))),
      pred = function(time, surv) cbind(log(time), log(-log(surv))),
      xlab = "LN(time)",
      ylab = "LN(cumulative hazard)"
    ),
    odds = list(
      obs  = function(time, surv) cbind(log(time), -log(surv / (1 - surv))),
      pred = function(time, surv) cbind(log(time), -log(surv / (1 - surv))),
      xlab = "LN(time)",
      ylab = "-LN(survival odds)"
    ),
    normal = list(
      obs  = function(time, surv) cbind(log(time), -qnorm(surv)),
      pred = function(time, surv) cbind(log(time), -qnorm(surv)),
      xlab = "LN(time)",
      ylab = "- standard normal quartiles"
    )
  )
  
  tf <- transform_funs[[spline_type]]
  
  # Observed data
  obs_data <- tf$obs(
    time = surv_obs$km$time,
    surv = surv_obs$km$surv
  )
  
  obs_col <- km_line_color[surv_obs$km_names]
  obs_pch <- point_shape[surv_obs$km_names]
  
  # Base plot
  plot(
    obs_data, cex = 0.6, pch = obs_pch, col = obs_col,
    main = paste("B: Diagnostic plot (", model_label, ")", sep = ""),
    xlab = tf$xlab,
    ylab = tf$ylab,
    xlim = c(
      floor(min(log(surv_obs$km$time)[surv_obs$km$surv != 1])),
      ceiling(max(log(surv_obs$km$time)))
    )
  )
  
  # Predicted lines for all groups
  for (i in seq_len(ngroups)) {
    pred_data <- tf$pred(
      time = input$time_pred,
      surv = surv_pred$model$spline[[model_key]]$surv[, i + 1]
    )
    lines(pred_data, col = as.character(col_index), lty = line_type[i], lwd = 2)
  }
  
  # Draw knots
  knots <- surv_model$spline_models[[model_key]]$knots
  if (!is.null(knots)) {
    if (length(knots) >= 2) {
      # boundary knots
      abline(v = knots[c(1, length(knots))], lty = 3, lwd = 1)
      # internal knots
      if (length(knots) > 2) {
        abline(v = knots[2:(length(knots) - 1)], lty = 2, lwd = 1.5)
      }
    }
  }
  
  # Legend
  legend(
    "topleft",
    legend = c(group_names, rep("", ngroups)),
    col = c(rep(as.character(col_index), ngroups), km_line_color[1:ngroups]),
    lty = c(line_type[1:ngroups], rep(NA, ngroups)),
    pch = c(rep(NA, ngroups), point_shape[1:ngroups]),
    cex = 0.8, ncol = 2, bty = "n"
  )
}

f_plot_diag_cure_surv_model <- function(PERSUADE, model_key, model_label, col_index = 1) {
  #' Plot Diagnostic Figures for Cure Survival Models
  #'
  #' Produces diagnostic plots for mixture and non-mixture cure models,
  #' with transformations depending on the distribution (Weibull, log-normal, log-logistic).
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_key Character. Name of the cure model in `surv_pred$model`.
  #' @param model_label Character. Main title of the plot.
  #' @param col_index Integer. Base color index for fitted curves (1–9).
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_diag_cure_surv_model(PERSUADE, "cure_weibull_mix", "Weibull mixture cure", 1)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "2222", "5212"))
  point_shape <- c(1, 8, 9)
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Detect cure model distribution
  dist_type <- if (grepl("weib", model_key, ignore.case = TRUE)) {
    "weibull"
  } else if (grepl("lnorm", model_key, ignore.case = TRUE)) {
    "lognormal"
  } else if (grepl("llog", model_key, ignore.case = TRUE)) {
    "loglogistic"
  } else {
    stop("Model key not recognized as Weibull, log-normal, or log-logistic cure model")
  }
  
  # Define transformations
  transform_funs <- list(
    weibull = list(
      obs  = function(time, surv) cbind(log(time), log(-log(surv))),
      pred = function(time, surv) cbind(log(time), log(-log(surv))),
      xlab = "LN(time)",
      ylab = "LN(cumulative hazard)"
    ),
    lognormal = list(
      obs  = function(time, surv) cbind(log(time), -qnorm(surv)),
      pred = function(time, surv) cbind(log(time), -qnorm(surv)),
      xlab = "LN(time)",
      ylab = "- standard normal quartiles"
    ),
    loglogistic = list(
      obs  = function(time, surv) cbind(log(time), -log(surv / (1 - surv))),
      pred = function(time, surv) cbind(log(time), -log(surv / (1 - surv))),
      xlab = "LN(time)",
      ylab = "-LN(survival odds)"
    )
  )
  
  tf <- transform_funs[[dist_type]]
  
  # Observed data
  obs_data <- tf$obs(
    time = surv_obs$km$time,
    surv = surv_obs$km$surv
  )
  
  obs_col <- km_line_color[surv_obs$km_names]
  obs_pch <- point_shape[surv_obs$km_names]
  
  # Base plot
  plot(
    obs_data, cex = 0.6, pch = obs_pch, col = obs_col,
    main = paste("B: Diagnostic plot (", model_label, ")", sep = ""),
    xlab = tf$xlab,
    ylab = tf$ylab,
    xlim = c(
      floor(min(log(surv_obs$km$time)[surv_obs$km$surv != 1])),
      ceiling(max(log(surv_obs$km$time)))
    )
  )
  
  # Predicted lines for all groups
  for (i in seq_len(ngroups)) {
    pred_data <- tf$pred(
      time = input$time_pred,
      surv = surv_pred$model$cure[[model_key]]$surv[, i + 1]
    )
    lines(pred_data, col = as.character(col_index), lty = line_type[i], lwd = 2)
  }
  
  # Legend
  legend(
    "topleft",
    legend = c(group_names, rep("", ngroups)),
    col = c(rep(as.character(col_index), ngroups), km_line_color[1:ngroups]),
    lty = c(line_type[1:ngroups], rep(NA, ngroups)),
    pch = c(rep(NA, ngroups), point_shape[1:ngroups]),
    cex = 0.8, ncol = 2, bty = "n"
  )
}

f_plot_tp_param_surv_model <- function(PERSUADE, model_key, model_label, col_index = 1) {
  #' Plot Annual Transition Probabilities for Parametric Survival Models
  #'
  #' Plots smoothed observed annual transition probabilities alongside
  #' model-predicted probabilities for each group, with shaded confidence intervals.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_key Character. Name of the model as stored in `surv_pred$tp_gr`.
  #' @param model_label Character. Main title of the plot.
  #' @param col_index Integer. Base color index for fitted curves (1–9).
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_tp_param_surv_model(PERSUADE, "weib", "Weibull", 2)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "2222", "5212"))
  point_shape <- c(1, 8, 9)
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Map model_key to column index in tp_gr predictions
  model_map <- list(
    expo  = 2,
    weib  = 3,
    gom   = 4,
    lnorm = 5,
    llog  = 6,
    gam   = 7,
    ggam  = 8
  )
  if (!model_key %in% names(model_map)) {
    stop("Model key not recognized in model_map")
  }
  pred_col_index <- model_map[[model_key]]
  
  # Base plot — group 1
  plot(
    x = surv_obs$tp$gr_1$time, y = surv_obs$tp$gr_1$smooth,
    type = "b", lwd = 1, lty = 1, cex = 0.6, pch = point_shape[1],
    col = km_line_color[1],
    main = paste("C: Annual transition probability (", model_label, ")", sep = ""),
    xlab = "time", ylab = "annual transition probability (smoothed)",
    ylim = c(0, ceiling(surv_obs$tp$max * 10) / 10),
    xlim = c(0, ceiling(surv_obs$km$maxtime * 10) / 10)
  )
  lines(
    x = input$time_pred[-1],
    y = unlist(surv_pred$tp_gr$gr_1[, ..pred_col_index]),
    col = as.character(col_index), lty = line_type[1], lwd = 2
  )
  
  # Additional groups
  for (i in 2:ngroups) {
    gr_name <- paste0("gr_", i)
    lines(
      x = surv_obs$tp[[gr_name]]$time,
      y = surv_obs$tp[[gr_name]]$smooth,
      col = km_line_color[i], type = "b", lwd = 1,
      cex = 0.6, pch = point_shape[i]
    )
    lines(
      x = input$time_pred[-1],
      y = unlist(surv_pred$tp_gr[[gr_name]][, ..pred_col_index]),
      col = as.character(col_index), lty = line_type[i], lwd = 2
    )
  }
  
  # Shaded CIs
  for (i in seq_len(ngroups)) {
    gr_name <- paste0("gr_", i)
    polygon(
      x = c(surv_obs$tp[[gr_name]]$time, rev(surv_obs$tp[[gr_name]]$time)),
      y = c(surv_obs$tp[[gr_name]]$smooth_lower,
            rev(surv_obs$tp[[gr_name]]$smooth_upper)),
      col = adjustcolor(km_line_color[i], alpha.f = 0.3),
      border = NA
    )
  }
  
  # Legend
  legend(
    "topleft",
    legend = c(group_names, rep("", ngroups)),
    col = c(rep(as.character(col_index), ngroups), km_line_color[1:ngroups]),
    lty = c(line_type[1:ngroups], rep(NA, ngroups)),
    pch = c(rep(NA, ngroups), point_shape[1:ngroups]),
    cex = 0.8, ncol = 2, bty = "n"
  )
}

f_plot_tp_spline_surv_model <- function(PERSUADE, model_key, model_label, col_index = 1, pred_col_index) {
  #' Plot Annual Transition Probabilities for Spline Survival Models
  #'
  #' Plots smoothed observed annual transition probabilities alongside
  #' spline model predictions, with shaded confidence intervals and knot lines.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_key Character. Name of the spline model in `surv_model` (e.g., "spl_hazard1").
  #' @param model_label Character. Main title of the plot.
  #' @param col_index Integer. Base colour index for fitted curves (1–9).
  #' @param pred_col_index Integer. Column index in `surv_pred$tp_gr` with model predictions.
  #'
  #' @return A base R plot.
  #' @export
  #'
  #' @examples
  #' f_plot_tp_spline_surv_model(PERSUADE, "spl_hazard_1", "Spline, 1 knot, hazard scale", 1)
  input      <- PERSUADE$input
  misc       <- PERSUADE$misc
  surv_obs   <- PERSUADE$surv_obs
  surv_pred  <- PERSUADE$surv_pred
  surv_model <- PERSUADE$surv_model
  
  ngroups       <- misc$ngroups
  group_names   <- misc$group_names
  line_type     <- as.integer(c(1, "2222", "5212"))
  point_shape   <- c(1, 8, 9)
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Map model_key to column index in tp_gr predictions
  model_map <- list(
    spl_hazard_1  = 9,
    spl_hazard_2  = 10,
    spl_hazard_3   = 11,
    spl_odds_1 = 12,
    spl_odds_2  = 13,
    spl_odds_3   = 14,
    spl_normal_1  = 15,
    spl_normal_2  = 16,
    spl_normal_3  = 17
  )
  if (!model_key %in% names(model_map)) {
    stop("Model key not recognized in model_map")
  }
  pred_col_index <- model_map[[model_key]]
  
  # Base plot: group 1 observed
  plot(
    x = surv_obs$tp$gr_1$time,
    y = surv_obs$tp$gr_1$smooth,
    type = "b", lwd = 1, lty = 1, cex = 0.6, pch = point_shape[1],
    col = km_line_color[1],
    main = paste("C: Annual transition probability (", model_label, ")", sep = ""),
    xlab = "time",
    ylab = "annual transition probability (smoothed)",
    ylim = c(0, ceiling(surv_obs$tp$max * 10) / 10),
    xlim = c(0, ceiling(surv_obs$km$maxtime * 10) / 10)
  )
  lines(
    x = input$time_pred[-1],
    y = unlist(surv_pred$tp_gr$gr_1[, ..pred_col_index]),
    col = as.character(col_index),
    lty = line_type[1], lwd = 2
  )
  
  # Additional groups
  for (i in 2:ngroups) {
    gr_name <- paste0("gr_", i)
    lines(
      x = surv_obs$tp[[gr_name]]$time,
      y = surv_obs$tp[[gr_name]]$smooth,
      col = km_line_color[i], type = "b", lwd = 1,
      cex = 0.6, pch = point_shape[i]
    )
    lines(
      x = input$time_pred[-1],
      y = unlist(surv_pred$tp_gr[[gr_name]][, ..pred_col_index]),
      col = as.character(col_index),
      lty = line_type[i], lwd = 2
    )
  }
  
  # Shaded CIs
  for (i in seq_len(ngroups)) {
    gr_name <- paste0("gr_", i)
    polygon(
      x = c(surv_obs$tp[[gr_name]]$time, rev(surv_obs$tp[[gr_name]]$time)),
      y = c(surv_obs$tp[[gr_name]]$smooth_lower,
            rev(surv_obs$tp[[gr_name]]$smooth_upper)),
      col = adjustcolor(km_line_color[i], alpha.f = 0.3),
      border = NA
    )
  }
  
  # Vertical lines for knots
  knots <- exp(surv_model$spline_models[[model_key]]$knots)
  if (length(knots) >= 2) {
    # Internal knots (middle ones) get lty = 2, thicker
    internal <- knots[2:(length(knots) - 1)]
    boundary <- knots[c(1, length(knots))]
    abline(v = internal, lty = 2, lwd = 1.5)
    abline(v = boundary, lty = 3, lwd = 1)
  }
  
  # Legend
  legend(
    "topleft",
    legend = c(group_names, rep("", ngroups)),
    col = c(rep(as.character(col_index), ngroups), km_line_color[1:ngroups]),
    lty = c(line_type[1:ngroups], rep(NA, ngroups)),
    pch = c(rep(NA, ngroups), point_shape[1:ngroups]),
    cex = 0.8, ncol = 2, bty = "n"
  )
}

f_plot_tp_cure_surv_model <- function(PERSUADE, model_key, model_label, col_index = 1) {
  #' Plot Annual Transition Probabilities for Cure Survival Models
  #'
  #' Produces smoothed observed annual transition probabilities with shaded CIs,
  #' overlaid with predictions from a cure survival model (mixture or non-mixture).
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_key Character. Name of the cure model in `surv_model$cure`.
  #' @param model_label Character. Main title of the plot.
  #' @param col_index Integer. Base colour index for fitted curves (1–9).
  #'
  #' @return A base R plot.
  #' @export
  #'
  #' @examples
  #' f_plot_tp_cure_surv_model(PERSUADE, "cure_weibull_mix", "Weibull mixture cure", 1)
  input     <- PERSUADE$input
  misc      <- PERSUADE$misc
  surv_obs  <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  ngroups       <- misc$ngroups
  group_names   <- misc$group_names
  line_type     <- as.integer(c(1, "2222", "5212"))
  point_shape   <- c(1, 8, 9)
  km_line_color <- c("black", "lightgrey", "darkgrey")

  # Map model_key to column index in tp_gr predictions
  model_map <- list(
    cure_weibull_mix  = 18,
    cure_weibull_nmix  = 19,
    cure_lnorm_mix   = 20,
    cure_lnorm_nmix = 21,
    cure_llogis_mix  = 22,
    cure_llogis_nmix   = 23
  )
  if (!model_key %in% names(model_map)) {
    stop("Model key not recognized in model_map")
  }
  pred_col_index <- model_map[[model_key]]
  pred_col_index <- pred_col_index - if (input$spline_mod == FALSE) { 9 } else { 0 }
  
  # Base plot: group 1 observed
  plot(
    x = surv_obs$tp$gr_1$time,
    y = surv_obs$tp$gr_1$smooth,
    type = "b", lwd = 1, lty = 1, cex = 0.6, pch = point_shape[1],
    col = km_line_color[1],
    main = paste("C: Annual transition probability (", model_label, ")", sep = ""),
    xlab = "time",
    ylab = "annual transition probability (smoothed)",
    ylim = c(0, ceiling(surv_obs$tp$max * 10) / 10),
    xlim = c(0, ceiling(surv_obs$km$maxtime * 10) / 10)
  )
  lines(
    x = input$time_pred[-1],
    y = unlist(surv_pred$tp_gr$gr_1[, ..pred_col_index]),
    col = as.character(col_index),
    lty = line_type[1], lwd = 2
  )
  
  # Remaining groups
  for (i in 2:ngroups) {
    gr_name <- paste0("gr_", i)
    lines(
      x = surv_obs$tp[[gr_name]]$time,
      y = surv_obs$tp[[gr_name]]$smooth,
      col = km_line_color[i], type = "b", lwd = 1,
      cex = 0.6, pch = point_shape[i]
    )
    lines(
      x = input$time_pred[-1],
      y = unlist(surv_pred$tp_gr[[gr_name]][, ..pred_col_index]),
      col = as.character(col_index),
      lty = line_type[i], lwd = 2
    )
  }
  
  # Shaded CIs
  for (i in seq_len(ngroups)) {
    gr_name <- paste0("gr_", i)
    polygon(
      x = c(surv_obs$tp[[gr_name]]$time, rev(surv_obs$tp[[gr_name]]$time)),
      y = c(surv_obs$tp[[gr_name]]$smooth_lower,
            rev(surv_obs$tp[[gr_name]]$smooth_upper)),
      col = adjustcolor(km_line_color[i], alpha.f = 0.3),
      border = NA
    )
  }
  
  # Legend
  legend(
    "topleft",
    legend = c(group_names, rep("", ngroups)),
    col = c(rep(as.character(col_index), ngroups), km_line_color[1:ngroups]),
    lty = c(line_type[1:ngroups], rep(NA, ngroups)),
    pch = c(rep(NA, ngroups), point_shape[1:ngroups]),
    cex = 0.8, ncol = 2, bty = "n"
  )
}
