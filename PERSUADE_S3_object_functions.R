#### S3 object functions for PERSUADE ----

print.PERSUADE <- function(PERSUADE, ...) {
  #' Print Method for PERSUADE Objects
  #'
  #' Displays a brief summary of the PERSUADE object in the console.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param ... Additional arguments (currently unused).
  #'
  #' @return Invisibly returns the PERSUADE object.
  #' 
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' print(PERSUADE)
  #' }
  cat("PERSUADE Survival Analysis Object\n")
  cat("Analysis Name:", PERSUADE$name, "\n")
  cat("Number of Observations:", length(PERSUADE$input$years), "\n")
  cat("Groups:", paste(levels(PERSUADE$input$group), collapse = ", "), "\n")
  invisible(PERSUADE)
}

summary.PERSUADE <- function(object, type = "km", ...) {
  #' Summary Method for PERSUADE Objects
  #'
  #' The `type` argument controls which summary is produced:
  #'   - `"km"`: Kaplan–Meier estimates (default).
  #'   - `"surv_probs"`: Survival probabilities at specified prediction times for each group.
  #'   - `"gof"`: Goodness-of-fit statistics for standard parametric models.
  #'   - `"gof_spline"`: Goodness-of-fit statistics for spline models.
  #'   - `"gof_cure"`: Goodness-of-fit statistics for cure models (including cure fraction).
  #'
  #' @param object A PERSUADE object from `f_PERSUADE()`.
  #' @param type Character. Type of summary to return.
  #' @param ... Additional arguments passed to methods (currently unused).
  #'
  #' @return A data frame or list of data frames depending on `type`.
  #'
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' summary(PERSUADE, type = "surv_probs")
  #' }
  if (type == "km") {
    return(summary(object$surv_obs$km)$table)
    
  } else if  (type == "surv_probs") {
    n_groups <- object$misc$ngroups
    surv_tables <- vector("list", n_groups)
    names(surv_tables) <- paste0("Group_", object$misc$group_names)
    
    for (i in seq_len(n_groups)) {
      surv_mat <- object$surv_pred$gr[[i]][1 + object$input$time_pred_surv_table, ]
      surv_mat <- t(round(surv_mat, 3))[-1, ]
      colnames(surv_mat) <- paste("T=", object$surv_pred$gr[[i]][1 + object$input$time_pred_surv_table, ][, 1])
      surv_tables[[i]] <- as.data.frame(surv_mat)
    }
    return(surv_tables)
    
  } else if  (type == "gof") {
    return(object$surv_model$param_ic)
    
  } else if  (type == "gof_spline") {
    if (!isTRUE(object$input$spline_mod)) stop("No spline models identified")
    return(object$surv_model$spline_ic)
    
  } else if  (type == "gof_cure") {
    if (!isTRUE(object$input$cure_mod)) stop("No cure models identified")
    return(object$surv_model$cure_ic) 
    
  } else {
    stop("Unknown summary type: ", type)
  }
}

plot.PERSUADE <- function(PERSUADE, type = "km", ...) {
  #' Plot Method for PERSUADE Objects
  #'
  #' Generates diagnostic and model fit plots for PERSUADE survival analysis objects.
  #' The `type` argument controls which plot(s) are produced:
  #'   - `"km"`: Kaplan–Meier survival curves.
  #'   - `"ph"`: Proportional hazards diagnostics.
  #'   - `"hr"`: Hazard function with fitted models.
  #'   - `"param_models"`: Fitted parametric survival models with diagnostics and transition probability plots.
  #'   - `"spline_models"`: Fitted spline-based survival models with diagnostics and transition probability plots.
  #'   - `"cure_models"`: Fitted cure survival models with diagnostics and transition probability plots.
  #'
  #' @param PERSUADE A PERSUADE object from `f_PERSUADE()`.
  #' @param type Character. The type of plot to produce.
  #' @param ... Additional arguments (currently unused).
  #'
  #' @return Invisibly returns a list of results from the plotting functions.
  #'   Also produces base R plots as side effects.
  #' 
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' plot(PERSUADE, "km")
  #' }
  plots <- list()
  
  if (type == "km") {
    plots <- f_plot_km_survival_base(PERSUADE)
    
  } else if (type == "ph") {
    plots <- list(
      f_plot_log_cumhaz(PERSUADE),
      f_plot_schoenfeld_residuals(PERSUADE)
    )
    
  } else if (type == "hr") {
    plots <- f_plot_hazard_with_models(PERSUADE)
    
  } else if (type == "param_models") {
    n_models <- length(PERSUADE$misc$lbls)
    for (i in seq_len(n_models)) {
      plots[[paste0("model_", i, "_surv")]] <- f_plot_param_surv_model(PERSUADE, i)
      plots[[paste0("model_", i, "_diag")]] <- f_plot_diag_param_surv_model(PERSUADE, i)
      plots[[paste0("model_", i, "_tp")]]   <- f_plot_tp_param_surv_model(PERSUADE, i)
    }
    
  } else if (type == "spline_models") {
    n_models <- length(PERSUADE$misc$lbls_spline)
    for (i in seq_len(n_models)) {
      plots[[paste0("spline_", i, "_surv")]] <- f_plot_spline_surv_model(PERSUADE, i)
      plots[[paste0("spline_", i, "_diag")]] <- f_plot_diag_spline_surv_model(PERSUADE, i)
      plots[[paste0("spline_", i, "_tp")]]   <- f_plot_tp_spline_surv_model(PERSUADE, i)
    }
    
  } else if (type == "cure_models") {
    n_models <- length(PERSUADE$misc$lbls_cure)
    for (i in seq_len(n_models)) {
      plots[[paste0("cure_", i, "_surv")]] <- f_plot_cure_surv_model(PERSUADE, i)
      plots[[paste0("cure_", i, "_diag")]] <- f_plot_diag_cure_surv_model(PERSUADE, i)
      plots[[paste0("cure_", i, "_tp")]]   <- f_plot_tp_cure_surv_model(PERSUADE, i)
    }
    
  } else {
    stop("Unknown plot type: ", type)
  }
  
  invisible(plots)
}
