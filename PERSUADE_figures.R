#### Functions used to create Figures and other output (e.g. for use in an R Markdown file) ----
f_plot_km_survival <- function(PERSUADE) {
  #' Plot Kaplan-Meier Survival Curves (ggsurvplot)
  #'
  #' Generates Kaplan-Meier survival plots from a PERSUADE object using
  #' [survminer::ggsurvplot()], automatically adapting to the number of groups.
  #'
  #' @param PERSUADE A PERSUADE object created by [f_PERSUADE()].
  #'
  #' @return A `ggsurvplot` object with KM curves, risk table, CI bands, and optional censor marks.
  #'
  #' @examples
  #' \dontrun{
  #' f_plot_km_survival(PERSUADE)
  #' }
  #' @export
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
    linetype = if (ngroups > 1) as.integer(c(1, "3333", "5212")),
    size = 0.5,
    legend = "top"
  )
  
  return(plot_obj)
}

f_plot_km_survival_base <- function(PERSUADE) {
  #' Plot Kaplan-Meier Survival Curves (Base R)
  #'
  #' Generates Kaplan-Meier survival plots from a PERSUADE object using base R graphics,
  #' with shaded confidence intervals and group-specific legends.
  #'
  #' @param PERSUADE A PERSUADE object created by [f_PERSUADE()].
  #'
  #' @return A base R plot.
  #'
  #' @examples
  #' \dontrun{
  #' f_plot_km_survival_base(PERSUADE)
  #' }
  #' @export
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "3333", "5212"))
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Base KM plot
  plot(
    surv_obs$km, lwd = 2, col = km_line_color[1:ngroups],
    main = "Kaplan-Meier", lty = line_type,
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
  
  # Add legend with model and KM lines
  legend(
    "bottomleft",
    legend = group_names,
    col = km_line_color[1:ngroups],
    lty = c(line_type[1:ngroups], rep(1, ngroups)),
    cex = 0.8, ncol = 2, bty = "n"
  )
}

f_plot_log_cumhaz <- function(PERSUADE) {
  #' Log-Log Survival Diagnostic Plot
  #'
  #' Creates a log(-log(S(t))) vs log(time) plot to visually assess proportional hazards.
  #'
  #' @param PERSUADE A PERSUADE object created by [f_PERSUADE()].
  #'
  #' @return A base R plot showing ln(-ln(S(t))) against ln(time).
  #'
  #' @examples
  #' \dontrun{
  #' f_plot_log_cumhaz(PERSUADE)
  #' }
  #' @export
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
  #' Schoenfeld Residuals Plot
  #'
  #' Produces scaled Schoenfeld residual plots with fitted regression lines
  #' to evaluate Cox proportional hazards assumptions.
  #'
  #' @param PERSUADE A PERSUADE object created by [f_PERSUADE()].
  #'
  #' @return One or more base R plots, one per group comparison.
  #'
  #' @examples
  #' \dontrun{
  #' f_plot_schoenfeld_residuals(PERSUADE)
  #' }
  #' @export
  if (PERSUADE$misc$ngroups < 2) {
    stop("Scaled Schoenfeld residuals not available for a single group")
  }
  
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
  #' Smoothed Hazard Function Plot
  #'
  #' Plots smoothed hazard estimates for each group in the PERSUADE object.
  #'
  #' @param PERSUADE A PERSUADE object created by [f_PERSUADE()].
  #'
  #' @return A base R plot of smoothed hazards by group.
  #'
  #' @examples
  #' \dontrun{
  #' f_plot_smoothed_hazard(PERSUADE)
  #' }
  #' @export
  surv_obs <- PERSUADE$surv_obs
  misc <- PERSUADE$misc
  ngroups <- misc$ngroups
  
  # Define plotting colors and line types (extendable, will recycle if needed)
  haz_line_color <- c("black", "lightgrey", "darkgrey")
  line_type <- as.integer(c(1, "3333", "5212", "3313", "1144", "42"))
  
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
  #' Hazard Plot with Model Overlays
  #'
  #' Plots observed smoothed hazard estimates together with hazard predictions
  #' from parametric, spline, and cure survival models (if fitted).
  #'
  #' @param PERSUADE A PERSUADE object created by [f_PERSUADE()].
  #'
  #' @return A series of base R plots, one per group, with hazard overlays by model family.
  #'
  #' @examples
  #' \dontrun{
  #' f_plot_hazard_with_models(PERSUADE)
  #' }
  #' @export
  input <- PERSUADE$input
  surv_pred <- PERSUADE$surv_pred
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  ngroups <- misc$ngroups
  
  # Define base plot styles
  haz_line_color <- c("black", "lightgrey", "darkgrey")
  line_color <- 1:9
  line_type <- c(1, 3333, 5212)
  
  # Define model groups to iterate over
  model_types <- list(
    base = names(PERSUADE$surv_model$param_models),
    spline = names(PERSUADE$surv_model$spline_models),
    cure = names(PERSUADE$surv_model$cure_models)
  )
  
  # Determine cure column per group
  cure_col_index <- c(2, 3, 4)  # assumes group 1: col 2, group 2: col 3, etc.
  
  for (i in seq_len(ngroups)) {
    group_label <- paste("Group:", misc$group_names[i])
    smooth_name <- paste0("smooth_gr", i)
    smooth_data <- surv_obs$haz$hazards[[smooth_name]]
    
    plot(
      smooth_data$est.grid, smooth_data$haz.est,
      main = paste0("Standard parametric models, ", group_label),
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
        main = paste0("Spline models, ", group_label),
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
        main = paste0("Cure models, ", group_label),
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

f_plot_param_surv_model <- function(PERSUADE, model_index = 1) {
  #' Parametric Survival Model Overlay
  #'
  #' Overlays a fitted parametric survival model on top of KM curves, including
  #' shaded KM confidence bands per group.
  #'
  #' @param PERSUADE A PERSUADE object created by [f_PERSUADE()].
  #' @param model_index Integer. Index of the parametric model in `PERSUADE$surv_pred$model`.
  #'
  #' @return A base R plot of KM curves with parametric model overlays.
  #'
  #' @examples
  #' \dontrun{
  #' f_plot_param_surv_model(PERSUADE, model_index = 1)
  #' }
  #' @export
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  model_key = names(surv_pred$model)[model_index]
  model_label = misc$lbls[model_index]
  col_index = model_index
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "3333", "5212"))
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

f_plot_spline_surv_model <- function(PERSUADE, model_index = 1) {
  #' Spline Survival Model Overlay
  #'
  #' Overlays a spline-based survival model on KM curves, including shaded KM
  #' confidence bands and vertical lines for knot positions.
  #'
  #' @param PERSUADE A PERSUADE object created by [f_PERSUADE()].
  #' @param model_index Integer. Index of the spline model in `PERSUADE$surv_pred$model$spline`.
  #'
  #' @return A base R plot of KM curves with spline model overlays and knots.
  #'
  #' @examples
  #' \dontrun{
  #' f_plot_spline_surv_model(PERSUADE, model_index = 1)
  #' }
  #' @export
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  surv_model <- PERSUADE$surv_model
  
  model_key = names(surv_pred$model$spline)[model_index]
  model_label = misc$lbls_spline[model_index]
  col_index = model_index
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "3333", "5212"))
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

f_plot_cure_surv_model <- function(PERSUADE, model_index = 1) {
  #' Cure Survival Model Overlay
  #'
  #' Overlays a fitted cure survival model on KM curves, including shaded KM
  #' confidence bands per group.
  #'
  #' @param PERSUADE A PERSUADE object created by [f_PERSUADE()].
  #' @param model_index Integer. Index of the cure model in `PERSUADE$surv_pred$model$cure`.
  #'
  #' @return A base R plot of KM curves with cure model overlays.
  #'
  #' @examples
  #' \dontrun{
  #' f_plot_cure_surv_model(PERSUADE, model_index = 1)
  #' }
  #' @export
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  model_key = names(surv_pred$model$cure)[model_index]
  model_label = misc$lbls_cure[model_index]
  col_index = model_index
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "3333", "5212"))
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

f_plot_diag_param_surv_model <- function(PERSUADE, model_index = 1) {
  #' Plot Diagnostic Figures for Standard Parametric Survival Models
  #'
  #' Produces diagnostic plots for various parametric survival models,
  #' using transformations specific to each model family.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_index Integer. Index for model selection.
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_diag_param_surv_model(PERSUADE, 1)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  model_key = names(surv_pred$model)[model_index]
  model_label = misc$lbls[model_index]
  col_index = model_index
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "3333", "5212"))
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

f_plot_diag_spline_surv_model <- function(PERSUADE, model_index = 1) {
  #' Plot Diagnostic Figures for Spline-based Survival Models
  #'
  #' Produces diagnostic plots for spline-based survival models,
  #' with transformations depending on the spline scale (hazard, odds, normal).
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_index Integer. Index for model selection.
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_diag_spline_surv_model(PERSUADE, 1)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  surv_model <- PERSUADE$surv_model
  
  model_key = names(surv_pred$model$spline)[model_index]
  model_label = misc$lbls_spline[model_index]
  col_index = model_index
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "3333", "5212"))
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

f_plot_diag_cure_surv_model <- function(PERSUADE, model_index = 1) {
  #' Plot Diagnostic Figures for Cure Survival Models
  #'
  #' Produces diagnostic plots for mixture and non-mixture cure models,
  #' with transformations depending on the distribution (Weibull, log-normal, log-logistic).
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_index Integer. Index for model selection.
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_diag_cure_surv_model(PERSUADE, 1)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  model_key = names(surv_pred$model$cure)[model_index]
  model_label = misc$lbls_cure[model_index]
  col_index = model_index
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "3333", "5212"))
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

f_plot_tp_param_surv_model <- function(PERSUADE, model_index = 1) {
  #' Plot Annual Transition Probabilities for Parametric Survival Models
  #'
  #' Plots smoothed observed annual transition probabilities alongside
  #' model-predicted probabilities for each group, with shaded confidence intervals.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_index Integer. Index for model selection.
  #'
  #' @return A base R plot.
  #' @export
  #' 
  #' @examples
  #' f_plot_tp_param_surv_model(PERSUADE, 1)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  model_key = names(surv_pred$model)[model_index]
  model_label = misc$lbls[model_index]
  col_index = model_index
  
  ngroups <- misc$ngroups
  group_names <- misc$group_names
  line_type <- as.integer(c(1, "3333", "5212"))
  point_shape <- c(1, 8, 9)
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Map model_key to column index in tp_gr predictions
  model_map <- setNames(as.list(2:8), names(PERSUADE$surv_model$param_models))
  
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
  if (ngroups > 1) {
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

f_plot_tp_spline_surv_model <- function(PERSUADE, model_index = 1) {
  #' Plot Annual Transition Probabilities for Spline Survival Models
  #'
  #' Plots smoothed observed annual transition probabilities alongside
  #' spline model predictions, with shaded confidence intervals and knot lines.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_index Integer. Index for model selection.
  #'
  #' @return A base R plot.
  #' @export
  #'
  #' @examples
  #' f_plot_tp_spline_surv_model(PERSUADE, 1)
  input      <- PERSUADE$input
  misc       <- PERSUADE$misc
  surv_obs   <- PERSUADE$surv_obs
  surv_pred  <- PERSUADE$surv_pred
  surv_model <- PERSUADE$surv_model
  
  model_key = names(surv_pred$model$spline)[model_index]
  model_label = misc$lbls_spline[model_index]
  col_index = model_index
  
  ngroups       <- misc$ngroups
  group_names   <- misc$group_names
  line_type     <- as.integer(c(1, "3333", "5212"))
  point_shape   <- c(1, 8, 9)
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Map model_key to column index in tp_gr predictions
  model_map <- setNames(as.list(9:17), names(PERSUADE$surv_model$spline_models))
  
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
  if (ngroups > 1) {
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

f_plot_tp_cure_surv_model <- function(PERSUADE, model_index = 1) {
  #' Plot Annual Transition Probabilities for Cure Survival Models
  #'
  #' Produces smoothed observed annual transition probabilities with shaded CIs,
  #' overlaid with predictions from a cure survival model (mixture or non-mixture).
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param model_index Integer. Index for model selection.
  #'
  #' @return A base R plot.
  #' @export
  #'
  #' @examples
  #' f_plot_tp_cure_surv_model(PERSUADE, 1)
  input     <- PERSUADE$input
  misc      <- PERSUADE$misc
  surv_obs  <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  model_key = names(surv_pred$model$cure)[model_index]
  model_label = misc$lbls_cure[model_index]
  col_index = model_index
  
  ngroups       <- misc$ngroups
  group_names   <- misc$group_names
  line_type     <- as.integer(c(1, "3333", "5212"))
  point_shape   <- c(1, 8, 9)
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  # Map model_key to column index in tp_gr predictions
  model_map <- setNames(as.list(18:23), names(PERSUADE$surv_model$cure_models))
  
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
  
  # Additional groups
  if (ngroups > 1) {
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

f_plot_param_surv_extrap <- function(PERSUADE) {
  #' Plot Extrapolated Parametric Survival Models per Group
  #'
  #' Produces Kaplan-Meier curves per group overlaid with fitted parametric
  #' survival models, extrapolated up to the defined time horizon.
  #' The models include: Exponential, Weibull, Gompertz, Log-normal,
  #' Log-logistic, Gamma, and Generalized Gamma.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #'
  #' @return A base R plot per group.
  #' @export
  #'
  #' @examples
  #' f_plot_param_surv_extrap(PERSUADE)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  model_names <- names(PERSUADE$surv_model$param_models)
  
  line_type <- as.integer(c(1, "3333", "5212"))
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  for (i in seq_len(misc$ngroups)) {
    # Base KM plot
    plot(
      surv_obs$km[i], lwd = 2, col = km_line_color[i], conf.int = FALSE, lty = line_type[i],
      main = paste0("A: Kaplan-Meier (parametric curves), Group: ", misc$group_names[i]),
      xlab = "time", ylab = "survival",
      xlim = c(0, input$time_horizon)
    )
    
    # Add shaded CI per group
    idx <- surv_obs$km_names == i
    polygon(
      na.omit(data.frame(
        x = c(surv_obs$km$time[idx], rev(surv_obs$km$time[idx])),
        y = c(surv_obs$km$lower[idx], rev(surv_obs$km$upper[idx]))
      )),
      col = adjustcolor(km_line_color[i], alpha.f = 0.3),
      border = NA
    )
    
    # Add parametric curves
    for (j in seq_along(model_names)) {
      lines(
        x = input$time_pred,
        y = surv_pred$model[[model_names[j]]]$surv[, i + 1],
        col = j,
        lty = line_type[i],
        lwd = 1
      )
    }
    
    # Add legend  
    legend("topright", legend = misc$lbls, col = 1:length(model_names),
           lty = line_type[i], cex = 0.8)
  }
}

f_plot_spline_surv_extrap <- function(PERSUADE) {
  #' Plot Extrapolated Spline Survival Models per Group
  #'
  #' Produces Kaplan-Meier curves per group overlaid with fitted spline
  #' survival models, extrapolated up to the defined time horizon.
  #' Models include: spline hazard, odds, and normal forms.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #'
  #' @return A base R plot per group.
  #' @export
  #'
  #' @examples
  #' f_plot_spline_surv_extrap(PERSUADE)
  if (!isTRUE(PERSUADE$input$spline_mod)) return(invisible())
  
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  model_names <- names(PERSUADE$surv_model$spline_models)
  
  line_type <- as.integer(c(1, "3333", "5212"))
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  for (i in seq_len(misc$ngroups)) {
    # Base KM plot
    plot(
      surv_obs$km[i], lwd = 2, col = km_line_color[i], conf.int = FALSE, lty = line_type[i],
      main = paste0("A: Kaplan-Meier (spline curves), Group: ", misc$group_names[i]),
      xlab = "time", ylab = "survival",
      xlim = c(0, input$time_horizon)
    )
    
    # Add shaded CI per group
    idx <- surv_obs$km_names == i
    polygon(
      na.omit(data.frame(
        x = c(surv_obs$km$time[idx], rev(surv_obs$km$time[idx])),
        y = c(surv_obs$km$lower[idx], rev(surv_obs$km$upper[idx]))
      )),
      col = adjustcolor(km_line_color[i], alpha.f = 0.3),
      border = NA
    )
    
    # Add parametric curves
    for (j in seq_along(model_names)) {
      lines(
        x = input$time_pred,
        y = surv_pred$model$spline[[model_names[j]]]$surv[, i + 1],
        col = j,
        lty = line_type[i],
        lwd = 1
      )
    }
    
    # Add legend 
    legend("topright", legend = misc$lbls_spline, col = 1:length(model_names),
           lty = line_type[i], cex = 0.8)
  }
}

f_plot_cure_surv_extrap <- function(PERSUADE) {
  #' Plot Extrapolated Cure Survival Models per Group
  #'
  #' Produces Kaplan-Meier curves per group overlaid with fitted cure
  #' survival models, extrapolated up to the defined time horizon.
  #' Models include: Weibull, Log-normal, and Log-logistic (mixture and non-mixture).
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #'
  #' @return A base R plot per group.
  #' @export
  #'
  #' @examples
  #' f_plot_cure_surv_extrap(PERSUADE)
  if (!isTRUE(PERSUADE$input$cure_mod)) return(invisible())
  
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  model_names <- names(PERSUADE$surv_model$cure_models)
  
  line_type <- as.integer(c(1, "3333", "5212"))
  km_line_color <- c("black", "lightgrey", "darkgrey")
  
  for (i in seq_len(misc$ngroups)) {
    # Base KM plot
    plot(
      surv_obs$km[i], lwd = 2, col = km_line_color[i], conf.int = FALSE, lty = line_type[i],
      main = paste0("A: Kaplan-Meier (cure curves), Group: ", misc$group_names[i]),
      xlab = "time", ylab = "survival",
      xlim = c(0, input$time_horizon)
    )
    
    # Add shaded CI per group
    idx <- surv_obs$km_names == i
    polygon(
      na.omit(data.frame(
        x = c(surv_obs$km$time[idx], rev(surv_obs$km$time[idx])),
        y = c(surv_obs$km$lower[idx], rev(surv_obs$km$upper[idx]))
      )),
      col = adjustcolor(km_line_color[i], alpha.f = 0.3),
      border = NA
    )
    
    # Add parametric curves
    for (j in seq_along(model_names)) {
      lines(
        x = input$time_pred,
        y = surv_pred$model$cure[[model_names[j]]]$surv[, i + 1],
        col = j,
        lty = line_type[i],
        lwd = 1
      )
    }
    
    # Add legend  
    legend("topright", legend = misc$lbls_cure, col = 1:length(model_names),
           lty = line_type[i], cex = 0.8)
  }
}

f_plot_tp_param_surv_extrap <- function(PERSUADE) {
  #' Plot Extrapolated Annual Transition Probabilities (Parametric Models)
  #'
  #' Produces smoothed observed annual transition probabilities with shaded CIs,
  #' overlaid with predictions from parametric survival models.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #'
  #' @return A base R plot per group.
  #' @export
  #'
  #' @examples
  #' f_plot_tp_param_surv_extrap(PERSUADE)
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  line_color <- c("black", "lightgrey", "darkgrey")
  line_type <- as.integer(c(1, "3333", "5212"))
  
  for (i in seq_len(misc$ngroups)) {
    tp_obs <- surv_obs$tp[[paste0("gr_", i)]]
    tp_pred <- surv_pred$tp_gr[[paste0("gr_", i)]]
    
    plot(tp_obs$time, tp_obs$smooth,
         main =  paste0("B: Annual transition probability (parametric curves), Group: ", misc$group_names[i]),
         xlab = "time", ylab = "annual transition probability",
         ylim = c(0, 1), xlim = c(0, input$time_horizon),
         lty = line_type[i], col = line_color[i], type = "l", lwd = 2)
    
    # Add shaded CI per group
    polygon(
      na.omit(data.frame(
        x = c(tp_obs$time, rev(tp_obs$time)),
        y = c(tp_obs$smooth_upper, rev(tp_obs$smooth_lower))
      )),
      col = adjustcolor(line_color[i], alpha.f = 0.3),
      border = NA
    )
    
    for (j in 2:8) {
      lines(input$time_pred[-1], unlist(tp_pred[, ..j]), col = j - 1, lty = line_type[i], lwd = 1)
    }
    
    legend("topleft", legend = misc$lbls, col = 1:7, lty = line_type[i], cex = 0.8)
  }
}

f_plot_tp_spline_surv_extrap <- function(PERSUADE) {
  #' Plot Extrapolated Annual Transition Probabilities (Spline Models)
  #'
  #' Produces smoothed observed annual transition probabilities with shaded CIs,
  #' overlaid with predictions from spline survival models.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #'
  #' @return A base R plot per group.
  #' @export
  #'
  #' @examples
  #' f_plot_tp_spline_surv_extrap(PERSUADE)
  if (!isTRUE(PERSUADE$input$spline_mod)) return(invisible())
  
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  line_color <- c("black", "lightgrey", "darkgrey")
  line_type <- as.integer(c(1, "3333", "5212"))
  
  for (i in seq_len(misc$ngroups)) {
    tp_obs <- surv_obs$tp[[paste0("gr_", i)]]
    tp_pred <- surv_pred$tp_gr[[paste0("gr_", i)]]
    
    plot(tp_obs$time, tp_obs$smooth,
         main =  paste0("B: Annual transition probability (spline curves), Group: ", misc$group_names[i]),
         xlab = "time", ylab = "annual transition probability",
         ylim = c(0, 1), xlim = c(0, input$time_horizon),
         lty = line_type[i], col = line_color[i], type = "l", lwd = 2)
    
    # Add shaded CI per group
    polygon(
      na.omit(data.frame(
        x = c(tp_obs$time, rev(tp_obs$time)),
        y = c(tp_obs$smooth_upper, rev(tp_obs$smooth_lower))
      )),
      col = adjustcolor(line_color[i], alpha.f = 0.3),
      border = NA
    )
    
    for (j in 9:17) {
      lines(input$time_pred[-1], unlist(tp_pred[, ..j]), col = j - 8, lty = line_type[i], lwd = 1)
    }
    
    legend("topleft", legend = misc$lbls_spline, col = 1:9, lty = line_type[i], cex = 0.8)
  }
}

f_plot_tp_cure_surv_extrap <- function(PERSUADE) {
  #' Plot Extrapolated Annual Transition Probabilities (Cure Models)
  #'
  #' Produces smoothed observed annual transition probabilities with shaded CIs,
  #' overlaid with predictions from cure survival models.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #'
  #' @return A base R plot per group.
  #' @export
  #'
  #' @examples
  #' f_plot_tp_cure_surv_extrap(PERSUADE)
  if (!isTRUE(PERSUADE$input$cure_mod)) return(invisible())
  
  input <- PERSUADE$input
  misc <- PERSUADE$misc
  surv_obs <- PERSUADE$surv_obs
  surv_pred <- PERSUADE$surv_pred
  
  line_color <- c("black", "lightgrey", "darkgrey")
  line_type <- as.integer(c(1, "3333", "5212"))
  offset <- if (isTRUE(input$spline_mod)) 0 else -9
  
  for (i in seq_len(misc$ngroups)) {
    tp_obs <- surv_obs$tp[[paste0("gr_", i)]]
    tp_pred <- surv_pred$tp_gr[[paste0("gr_", i)]]
    
    plot(tp_obs$time, tp_obs$smooth,
         main =  paste0("B: Annual transition probability (cure curves), Group: ", misc$group_names[i]),
         xlab = "time", ylab = "annual transition probability",
         ylim = c(0, 1), xlim = c(0, input$time_horizon),
         lty = line_type[i], col = line_color[i], type = "l", lwd = 2)
    
    # Add shaded CI per group
    polygon(
      na.omit(data.frame(
        x = c(tp_obs$time, rev(tp_obs$time)),
        y = c(tp_obs$smooth_upper, rev(tp_obs$smooth_lower))
      )),
      col = adjustcolor(line_color[i], alpha.f = 0.3),
      border = NA
    )
    
    for (j in 18:23) {
      j_offset <- j + offset
      lines(input$time_pred[-1], unlist(tp_pred[, ..j_offset]), col = j - 17, lty = line_type[i], lwd = 1)
    }
    
    legend("topleft", legend = misc$lbls_cure, col = 1:6, lty = line_type[i], cex = 0.8)
  }
}

f_plot_hazard_parametric_extrap <- function(PERSUADE) {
  #' Plot Extrapolated Hazard Functions (Parametric Models)
  #'
  #' Produces observed smoothed hazard rates overlaid with extrapolated hazard functions
  #' from parametric survival models.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #'
  #' @return Base R plots, one per group.
  #' @export
  #'
  #' @examples
  #' f_plot_hazard_parametric_extrap(PERSUADE)
  models <- names(PERSUADE$surv_model$param_models)
  misc <- PERSUADE$misc
  
  cols <- seq_along(models)
  line_color <- c("black", "lightgrey", "darkgrey")
  line_type <- as.integer(c(1, "3333", "5212"))
  
  for (i in seq_len(misc$ngroups)) {
    obs_data <- PERSUADE$surv_obs$haz$hazards[[paste0("smooth_gr", i)]]
    plot(cbind(obs_data$est.grid, obs_data$haz.est),
         main = paste0("C: Hazard function (parametric curves), Group: ", misc$group_names[i]),
         type = "l", col = line_color[i], lty = line_type[i], lwd = 2,
         xlab = "time", ylab = "smoothed hazard rate",
         xlim = c(0, PERSUADE$input$time_horizon),
         ylim = c(0, PERSUADE$surv_obs$haz$max$smooth))
    for (j in seq_along(models)) {
      lines(cbind(PERSUADE$surv_pred$model[[models[j]]]$hazard[, 1],
                  PERSUADE$surv_pred$model[[models[j]]]$hazard[, i + 1]),
            col = cols[j], lty = line_type[i], lwd = 1)
    }
    legend("topleft", legend = PERSUADE$misc$lbls,
           col = cols, lty = line_type[i], cex = 0.8)
  }
}

f_plot_hazard_spline_extrap <- function(PERSUADE) {
  #' Plot Extrapolated Hazard Functions (Spline Models)
  #'
  #' Produces observed smoothed hazard rates overlaid with extrapolated hazard functions
  #' from spline-based survival models.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #'
  #' @return Base R plots, one per group (only if `input$spline_mod` is TRUE).
  #' @export
  #'
  #' @examples
  #' f_plot_hazard_spline_extrap(PERSUADE)
  if (!isTRUE(PERSUADE$input$spline_mod)) return(invisible(NULL))
  
  models <- names(PERSUADE$surv_model$spline_models)
  misc <- PERSUADE$misc
  
  cols <- seq_along(models)
  line_color <- c("black", "lightgrey", "darkgrey")
  line_type <- as.integer(c(1, "3333", "5212"))
  
  for (i in seq_len(misc$ngroups)) {
    obs_data <- PERSUADE$surv_obs$haz$hazards[[paste0("smooth_gr", i)]]
    plot(cbind(obs_data$est.grid, obs_data$haz.est),
         main =  paste0("C: Hazard function (spline curves), Group: ", misc$group_names[i]),
         type = "l", col = line_color[i], lty = line_type[i], lwd = 2,
         xlab = "time", ylab = "smoothed hazard rate",
         xlim = c(0, PERSUADE$input$time_horizon),
         ylim = c(0, PERSUADE$surv_obs$haz$max$smooth))
    for (j in seq_along(models)) {
      lines(cbind(PERSUADE$surv_pred$model$spline[[models[j]]]$hazard[, 1],
                  PERSUADE$surv_pred$model$spline[[models[j]]]$hazard[, i + 1]),
            col = cols[j], lty = line_type[i], lwd = 1)
    }
    legend("topleft", legend = PERSUADE$misc$lbls_spline,
           col = cols, lty = line_type[i], cex = 0.8)
  }
}

f_plot_hazard_cure_extrap <- function(PERSUADE) {
  #' Plot Extrapolated Hazard Functions (Cure Models)
  #'
  #' Produces observed smoothed hazard rates overlaid with extrapolated hazard functions
  #' from cure survival models (mixture and non-mixture).
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #'
  #' @return Base R plots, one per group (only if `input$cure_mod` is TRUE).
  #' @export
  #'
  #' @examples
  #' f_plot_hazard_cure_extrap(PERSUADE)
  if (!isTRUE(PERSUADE$input$cure_mod)) return(invisible(NULL))
  
  models <- names(PERSUADE$surv_model$cure_models)
  misc <- PERSUADE$misc
  
  cols <- seq_along(models)
  line_color <- c("black", "lightgrey", "darkgrey")
  line_type <- as.integer(c(1, "3333", "5212"))
  
  for (i in seq_len(misc$ngroups)) {
    obs_data <- PERSUADE$surv_obs$haz$hazards[[paste0("smooth_gr", i)]]
    
    plot(cbind(obs_data$est.grid, obs_data$haz.est),
         main =  paste0("C: Hazard function (cure curves), Group: ", misc$group_names[i]),
         type = "l", col = line_color[i], lty = line_type[i], lwd = 2,
         xlab = "time", ylab = "smoothed hazard rate",
         xlim = c(0, PERSUADE$input$time_horizon),
         ylim = c(0, PERSUADE$surv_obs$haz$max$smooth))
    
    for (j in seq_along(models)) {
      lines(cbind(PERSUADE$surv_pred$model$cure[[models[j]]]$hazard[, 1],
                  PERSUADE$surv_pred$model$cure[[models[j]]]$hazard[, i + 1]),
            col = cols[j], lty = line_type[i], lwd = 1)
    }
    
    legend("topleft", legend = PERSUADE$misc$lbls_cure,
           col = cols, lty = line_type[i], cex = 0.8)
  }
}

f_summary <- function(df) {
  #' Compute Summary Statistics for Numeric Variables
  #'
  #' Produces descriptive statistics for each numeric variable in a data frame.
  #' The statistics include mean, standard deviation, minimum, first quartile (Q1),
  #' median, third quartile (Q3), maximum, and interquartile range (IQR),
  #' rounded to three decimal places.
  #'
  #' @param df A data frame containing numeric variables to summarise.
  #'
  #' @return A data frame with one row per variable and columns for each
  #'   summary statistic.
  #' @export
  #'
  #' @examples
  #' f_summary(mtcars)
  res <- t(sapply(df, function(x) {
    q <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE)
    round(c(
      Mean   = mean(x, na.rm = TRUE),
      Std.Dev= sd(x, na.rm = TRUE),
      Min    = q[1],
      Q1     = q[2],
      Median = q[3],
      Q3     = q[4],
      Max    = q[5],
      IQR    = q[4] - q[2]
    ), 3)
  }))
  as.data.frame(res, check.names = FALSE)
}

f_generate_report <- function(PERSUADE) {
  #' Generate PDF report for PERSUADE analysis
  #'
  #' Saves the PERSUADE object and generates a PDF report using `PERSUADE_output.Rmd`.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #'
  #' @return Path to generated PDF (invisible).
  #' @export
  #' 
  #' @examples
  #' f_generate_report(PERSUADE)
  name <- PERSUADE$name
  
  # Create output directories
  output_dir <- file.path(getwd(), paste0(name, "_output"))
  fig_dir    <- file.path(output_dir, "Images")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE)
  
  # Save PERSUADE object inside the output directory
  save(PERSUADE, file = file.path(output_dir, "PERSUADE.RData"))
  
  # Render R Markdown into PDF
  rmd_file <- "PERSUADE_output.Rmd"
  if (!file.exists(rmd_file)) stop("The file '", rmd_file, "' was not found.")
  
  # Ensure knit options inside Rmd can use this figure directory
  rmarkdown::render(
    input = rmd_file,
    output_file = paste0(name, ".pdf"),
    output_dir = output_dir,
    intermediates_dir = output_dir,
    knit_root_dir = output_dir,
    params = list(fig_dir = fig_dir), # pass fig_dir into the Rmd
    envir = list2env(list(PERSUADE = PERSUADE), parent = globalenv()),
    clean = TRUE
  )
  
  invisible(file.path(output_dir, paste0(name, ".pdf")))
}

