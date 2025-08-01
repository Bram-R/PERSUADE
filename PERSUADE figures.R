f_plot_km_survival <- function(PERSUADE, size = 0.5, legend_position = "top") {
  #' Plot Kaplan-Meier Survival Curves from PERSUADE Object
  #'
  #' Generates Kaplan-Meier survival plots using ggsurvplot, adapting automatically 
  #' to the number of groups in the PERSUADE object.
  #'
  #' @param PERSUADE A list object returned from the `f_PERSUADE()` function.
  #' @param size Numeric. Line size for survival curves.
  #' @param legend_position Character. Position of the legend ("right", "bottom", etc.).
  #'
  #' @return A `ggsurvplot` object showing KM curves with risk table and optional CI/rug.
  #' @export
  #'
  #' @examples
  #' plot_km_survival(PERSUADE)
  if (!is.list(PERSUADE)) {
    stop("`PERSUADE` must be a list created by f_PERSUADE().")
  }
  
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
    size = size,
    legend = legend_position
  )
  
  return(plot_obj)
}
