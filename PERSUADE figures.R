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
  #' @return A base R plot of smoothed hazard estimates by group.
  #' @export
  #'
  #' @examples
  #' f_plot_smoothed_hazard(PERSUADE)
  surv_obs <- PERSUADE$surv_obs
  misc <- PERSUADE$misc
  ngroups <- misc$ngroups
  
  # Define plotting colors and line types (extendable, will recycle if needed)
  line_color <- c("red", "green", "blue", "purple", "orange", "brown")
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
          col = line_color[(i - 1) %% length(line_color) + 1],
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
          col = line_color[(i - 1) %% length(line_color) + 1],
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
    col = line_color[seq_len(ngroups)],
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
      col = "black",
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
        col = "black",
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
        col = "black",
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

