## main PERSUADE function
f_PERSUADE <- function(name = "no_name", years, status, group, 
                       strata = FALSE, spline_mod = FALSE, cure_mod = FALSE, cure_link = "logistic",
                       time_unit, time_horizon, time_pred_surv_table) {
  
  # strata <- TRUE # for validation purposes
  # spline_mod <- TRUE # for validation purposes
  # cure_mod <- TRUE # for validation purposes
  # cure_link = "logistic" # for validation purposes
  
  #input
  years <- as.numeric(years)  # time variable should be numeric
  status <- as.numeric(status)  # status / event variable should be numeric
  group <- as.factor(group)  # group variable should be a factor 
  if (cure_mod == FALSE) {cure_link <- NA}
  
  # number of groups
  group <- droplevels(group)  # drop unused levels
  ngroups <- length(levels(group))  # number of groups
  group_names <- levels(group)  # name groups
  if (ngroups == 1) {strata <- FALSE}
  
  # time variables
  time_horizon <- max(time_pred_surv_table, time_horizon)
  time_pred_surv_table <- time_pred_surv_table/time_unit
  time_pred <- seq(from = 0, to = time_horizon, by = time_unit)
  
  # form
  if (ngroups == 1) {form <- as.formula(Surv(years, status) ~ 1)} else {
    form <- as.formula(Surv(years, status) ~ group)}
  
  # km
  km <- npsurv(form)
  km_names <- if (ngroups == 1) {rep(1, length(km$time))} else {
    if (ngroups > 1) {c(rep(1, km$strata[1]), rep(2, km$strata[2]), if (ngroups == 3) {rep(3, km$strata[3])})}
  }
  
  # hazard rate
  haz <- f_hazard(years = years, status = status, group = group, ngroups = ngroups)
  
  # cumulative hazard
  cum_haz <- f_cum_hazard(years = years, status = status, group = group, ngroups = ngroups, time_pred = time_pred, time_unit = time_unit)
  
  # annual transition probability based on observed data (km) 
  tp <- f_tp(ngroups = ngroups, cum_haz = cum_haz)
  
  # cox (for Scaled Schoenfeld residuals)
  cox_reg <- coxph(form)
  
  # fit parameteric models
  surv_model <- f_surv_model(years = years, status = status, group = group, ngroups = ngroups, form = form, strata = strata, spline_mod = spline_mod, cure_mod = cure_mod, cure_link = cure_link, group_names = group_names)
  
  # predicted values parametric models
  surv_model_pred <- f_surv_model_pred(ngroups = ngroups, time_pred = time_pred, surv_model = surv_model, spline_mod = spline_mod, cure_mod = cure_mod, group_names = group_names) 
  
  # predicted values parametric models per group
  surv_model_pred_gr <- f_surv_model_pred_gr(ngroups = ngroups, surv_model = surv_model, surv_model_pred = surv_model_pred, spline_mod = spline_mod, cure_mod = cure_mod) 
  
  # predicted annual transition probability based on parametric survival models
  cols_tp <- 8 + (if (spline_mod == TRUE) {9} else {0}) + (if (cure_mod == TRUE) {6} else {0})
  surv_model_pred_tp_gr <- f_surv_model_pred_tp_gr(ngroups = ngroups, time_pred = time_pred, time_unit = time_unit, surv_model_pred_gr = surv_model_pred_gr, cols_tp = cols_tp)
  
  # create output dataframe for use in MS Excel 
  surv_model_excel <- f_surv_model_excel(ngroups = ngroups, strata = strata, surv_model = surv_model, spline_mod = spline_mod, cure_mod = cure_mod) 
  
  # output
  input <- list(years = years, status = status, group = group, strata = strata, spline_mod = spline_mod, 
                cure_mod = cure_mod, cure_link = cure_link, time_unit = time_unit, time_horizon = time_horizon, 
                time_pred_surv_table = time_pred_surv_table, time_pred = time_pred)
  
  surv_obs <- list(km = km, km_names = km_names, cum_haz = cum_haz, haz = haz, tp = tp, cox_reg = cox_reg)
  
  surv_pred <- list(model = surv_model_pred, gr = surv_model_pred_gr, tp_gr = surv_model_pred_tp_gr)
  
  misc <- c(list(form = form, group_names = group_names, ngroups = ngroups, lbls = surv_model$IC[,1]), 
            if (spline_mod == TRUE) {list(lbls_spline = surv_model$IC_spl[,1])}, 
            if (cure_mod == TRUE) {list(lbls_cure = surv_model$IC_cure[,1])}, 
            list(cols_tp = cols_tp))
  
  output <- list(name = name, input = input, surv_obs = surv_obs, surv_model = surv_model, surv_pred = surv_pred, surv_model_excel = surv_model_excel, misc = misc)
  
  return(output) 
}

## functions used in the main PERSUADE function

# function to obtain hazard
f_hazard <- function(years, status, group, ngroups) {
  haz <- c(list(smooth_gr1 = muhaz(years, status, group == levels(group)[1])),
           if (ngroups > 1) {list(smooth_gr2 = muhaz(years, status, group == levels(group)[2]))},
           if (ngroups > 2) {list(smooth_gr3 = muhaz(years, status, group == levels(group)[3]))})
  haz_names <- c(rep(1, length(haz$smooth_gr1$est.grid)), 
                 if (ngroups > 1) {rep(2, length(haz$smooth_gr2$est.grid))} else {NA}, 
                 if (ngroups > 2) {rep(3, length(haz$smooth_gr3$est.grid))})
  haz_max <- data.frame(time = ceiling(max(haz$smooth_gr1$est.grid, 
                                           if (ngroups > 1) {haz$smooth_gr2$est.grid} else {NA}, 
                                           if (ngroups > 2) {haz$smooth_gr3$est.grid} else {NA}, na.rm = TRUE)), 
                        smooth = ceiling(max(haz$smooth_gr1$haz.est, 
                                             if (ngroups > 1) {haz$smooth_gr2$haz.est} else {NA}, 
                                             if (ngroups > 2) {haz$smooth_gr3$haz.est} else {NA}, na.rm = TRUE)))
  
  return(c(haz, list(names = haz_names, max = haz_max))) 
}

# function to obtain cumulative hazard
f_cum_hazard <- function(years, status, group, ngroups, time_pred, time_unit) {
  cum_haz <- data.frame(group = c(rep(1, length(time_pred[time_pred <= max(years[group == levels(group)[1]])])), 
                                  if (ngroups > 1) {rep(2, length(time_pred[time_pred <= max(years[group == levels(group)[2]])]))}, 
                                  if (ngroups > 2) {rep(3, length(time_pred[time_pred <= max(years[group == levels(group)[3]])]))}))
  cum_haz$time <- c(time_pred[0:length(cum_haz$group[cum_haz$group == 1])], 
                    if (ngroups > 1) {time_pred[0:length(cum_haz$group[cum_haz$group == 2])]}, 
                    if (ngroups > 2) {time_pred[0:length(cum_haz$group[cum_haz$group == 3])]})
  cum_haz$H <- c(estimateNAH(years[group == levels(group)[1]], status[group == levels(group)[1]])$H(cum_haz$time[cum_haz$group == 1]), 
                 if (ngroups > 1) {estimateNAH(years[group == levels(group)[2]], status[group == levels(group)[2]])$H(cum_haz$time[cum_haz$group == 2])}, 
                 if (ngroups > 2) {estimateNAH(years[group == levels(group)[3]], status[group == levels(group)[3]])$H(cum_haz$time[cum_haz$group == 3])})
  cum_haz$var <- c(estimateNAH(years[group == levels(group)[1]], status[group == levels(group)[1]])$Var(cum_haz$time[cum_haz$group == 1]), 
                   if (ngroups > 1) {estimateNAH(years[group == levels(group)[2]], status[group == levels(group)[2]])$Var(cum_haz$time[cum_haz$group == 2])}, 
                   if (ngroups > 2) {estimateNAH(years[group == levels(group)[3]], status[group == levels(group)[3]])$Var(cum_haz$time[cum_haz$group == 3])})
  cum_haz$H_upper <- pmax(0, cum_haz$H + sqrt(cum_haz$var) * qnorm(p = 0.975))
  cum_haz$H_lower <- pmax(0, cum_haz$H - sqrt(cum_haz$var) * qnorm(p = 0.975))
  cum_haz$H_delta <- c(shift(cum_haz$H[cum_haz$group == 1], 1L, type = "lag") - cum_haz$H[cum_haz$group == 1], 
                       if (ngroups > 1) {shift(cum_haz$H[cum_haz$group == 2], 1L, type = "lag") - cum_haz$H[cum_haz$group == 2]}, 
                       if (ngroups > 2) {shift(cum_haz$H[cum_haz$group == 3], 1L, type = "lag") - cum_haz$H[cum_haz$group == 3]})
  cum_haz$H_upper_delta <- c(shift(cum_haz$H_upper[cum_haz$group == 1], 1L, type = "lag") - cum_haz$H_upper[cum_haz$group == 1], 
                             if (ngroups > 1) {shift(cum_haz$H_upper[cum_haz$group == 2], 1L, type = "lag") - cum_haz$H_upper[cum_haz$group == 2]}, 
                             if (ngroups > 2) {shift(cum_haz$H_upper[cum_haz$group == 3], 1L, type = "lag") - cum_haz$H_upper[cum_haz$group == 3]})
  cum_haz$H_lower_delta <- c(shift(cum_haz$H_lower[cum_haz$group == 1], 1L, type = "lag") - cum_haz$H_lower[cum_haz$group == 1], 
                             if (ngroups > 1) {shift(cum_haz$H_lower[cum_haz$group == 2], 1L, type = "lag") - cum_haz$H_lower[cum_haz$group == 2]},
                             if (ngroups > 2) {shift(cum_haz$H_lower[cum_haz$group == 3], 1L, type = "lag") - cum_haz$H_lower[cum_haz$group == 3]})
  
  return(cum_haz) 
}

# function to obtain annual transition probability based on observed data (km)
f_tp <- function(ngroups, cum_haz) { 
  
  tp <- 1 - exp(cum_haz$H_delta)^(1/time_unit)
  tp_upper <- 1 - exp(cum_haz$H_upper_delta)^(1/time_unit)
  tp_lower <- 1 - exp(cum_haz$H_lower_delta)^(1/time_unit)
  tp_smooth <- pmax(0, pmin(1, c(NA, loess(tp[cum_haz$group == 1] ~ cum_haz$time[cum_haz$group ==  1])$fitted, 
                                 if (ngroups > 1) {c(NA, loess(tp[cum_haz$group == 2] ~ cum_haz$time[cum_haz$group == 2])$fitted)},
                                 if (ngroups > 2) {c(NA, loess(tp[cum_haz$group == 3] ~ cum_haz$time[cum_haz$group == 3])$fitted)})))
  tp_upper_smooth <- pmax(0, pmin(1, c(NA, loess(tp_upper[cum_haz$group == 1] ~ cum_haz$time[cum_haz$group == 1])$fitted, 
                                       if (ngroups > 1) {c(NA, loess(tp_upper[cum_haz$group == 2] ~ cum_haz$time[cum_haz$group == 2])$fitted)},
                                       if (ngroups > 2) {c(NA, loess(tp_upper[cum_haz$group == 3] ~ cum_haz$time[cum_haz$group == 3])$fitted)})))
  tp_lower_smooth <- pmax(0, pmin(1, c(NA, loess(tp_lower[cum_haz$group == 1] ~ cum_haz$time[cum_haz$group == 1])$fitted, 
                                       if (ngroups > 1) {c(NA, loess(tp_lower[cum_haz$group == 2] ~ cum_haz$time[cum_haz$group == 2])$fitted)}, 
                                       if (ngroups > 2) {c(NA, loess(tp_lower[cum_haz$group == 3] ~ cum_haz$time[cum_haz$group == 3])$fitted)})))
  
  km_tp_gr_1 <- data.frame(
    time = cum_haz$time[cum_haz$group==1],
    smooth = tp_smooth[cum_haz$group==1],
    smooth_lower = tp_lower_smooth[cum_haz$group==1],
    smooth_upper = tp_upper_smooth[cum_haz$group==1]
  )
  
  if (ngroups > 1) {
    km_tp_gr_2 <- data.frame(
      time = cum_haz$time[cum_haz$group==2],
      smooth = tp_smooth[cum_haz$group==2],
      smooth_lower = tp_lower_smooth[cum_haz$group==2],
      smooth_upper = tp_upper_smooth[cum_haz$group==2]
    )
  }
  if (ngroups > 2) {
    km_tp_gr_3 <- data.frame(
      time = cum_haz$time[cum_haz$group==3],
      smooth = tp_smooth[cum_haz$group==3],
      smooth_lower = tp_lower_smooth[cum_haz$group==3],
      smooth_upper = tp_upper_smooth[cum_haz$group==3]
    )
  }
  
  km_tp_max <- max(c(km_tp_gr_1$smooth_upper, 
                     if (ngroups > 1) {km_tp_gr_2$smooth_upper}, 
                     if (ngroups > 2) {km_tp_gr_3$smooth_upper}), 
                   na.rm = TRUE)
  
  return(c(list(gr_1 = km_tp_gr_1), 
           if (ngroups > 1) {list(gr_2 = km_tp_gr_2)}, 
           if (ngroups > 2) {list(gr_3 = km_tp_gr_3)}, 
           list(max = km_tp_max)))
}


# function to obtain parametric survival models
f_surv_model <- function(years, status, group, strata, ngroups, form, spline_mod, cure_mod, cure_link, group_names) {
  
  # fit parameteric models
  expo <- flexsurvreg(form, dist = "exp")
  weib <- if (strata == FALSE) 
    flexsurvreg(form, dist = "weibull") else {
      flexsurvreg(form, anc = list(shape = ~group), dist = "weibull")
    }
  gom <- if (strata == FALSE) 
    flexsurvreg(form, dist = "gompertz") else {
      flexsurvreg(form, anc = list(shape = ~group), dist = "gompertz")
    }
  lnorm <- if (strata == FALSE) 
    flexsurvreg(form, dist = "lnorm") else {
      flexsurvreg(form, anc = list(sdlog = ~group), dist = "lnorm")
    }
  llog <- if (strata == FALSE) 
    flexsurvreg(form, dist = "llogis") else {
      flexsurvreg(form, anc = list(shape = ~group), dist = "llogis")
    }
  gam <- if (strata == FALSE) 
    flexsurvreg(form, dist = "gamma") else {
      flexsurvreg(form, anc = list(shape = ~group), dist = "gamma")
    }
  ggam <- if (strata == FALSE) 
    flexsurvreg(form, dist = "gengamma") else {
      flexsurvreg(form, anc = list(sigma = ~group, Q = ~group), dist = "gengamma")
    }
  
  lbls <- c(" 1. Exponential", " 2. Weibull", " 3. Gompertz", " 4. Log-normal", " 5. Log-logistic", " 6. Gamma", 
            " 7. Generalised Gamma")
  
  # calculate AIC and BIC
  AIC <- c(expo$AIC, weib$AIC, gom$AIC, lnorm$AIC, llog$AIC, gam$AIC, ggam$AIC)
  BIC <- BIC(expo, weib, gom, lnorm, llog, gam, ggam)[2]
  IC <- data.frame(lbls, AIC, BIC, row.names = NULL)
  colnames(IC) <- c("Model", "AIC", "BIC")
  
  # fit spline models
  if (spline_mod == TRUE) {
    k <- 1  #number of knots
    spl_hazard1 <- if (strata == FALSE) {flexsurvspline(form, k = k, scale = "hazard")} else {
      flexsurvspline(form, k = k, scale = "hazard", anc = list(gamma1 = ~group, gamma2 = ~group))
    }
    spl_odds1 <- if (strata == FALSE) {flexsurvspline(form, k = k, scale = "odds")} else {
      flexsurvspline(form, k = k, scale = "odds", anc = list(gamma1 = ~group, gamma2 = ~group))
    }
    spl_normal1 <- if (strata == FALSE) {flexsurvspline(form, k = k, scale = "normal")} else {
      flexsurvspline(form, k = k, scale = "normal", anc = list(gamma1 = ~group, gamma2 = ~group))
    }
    
    k <- 2  #number of knots
    spl_hazard2 <- if (strata == FALSE) {flexsurvspline(form, k = k, scale = "hazard")} else {
      flexsurvspline(form, k = k, scale = "hazard", anc = list(gamma1 = ~group, gamma2 = ~group, gamma3 = ~group))
    }
    spl_odds2 <- if (strata == FALSE) {flexsurvspline(form, k = k, scale = "odds")} else {
      flexsurvspline(form, k = k, scale = "odds", anc = list(gamma1 = ~group, gamma2 = ~group, gamma3 = ~group))
    }
    spl_normal2 <- if (strata == FALSE) {flexsurvspline(form, k = k, scale = "normal")} else {
      flexsurvspline(form, k = k, scale = "normal", anc = list(gamma1 = ~group, gamma2 = ~group, gamma3 = ~group))
    }
    
    k <- 3  #number of knots
    spl_hazard3 <- if (strata == FALSE) {flexsurvspline(form, k = k, scale = "hazard")} else {
      flexsurvspline(form, k = k, scale = "hazard", anc = list(gamma1 = ~group, gamma2 = ~group, gamma3 = ~group, gamma4 = ~group))
    }
    spl_odds3 <- if (strata == FALSE) {flexsurvspline(form, k = k, scale = "odds")} else {
      flexsurvspline(form, k = k, scale = "odds", anc = list(gamma1 = ~group, gamma2 = ~group, gamma3 = ~group, gamma4 = ~group))
    }
    spl_normal3 <- if (strata == FALSE) {flexsurvspline(form, k = k, scale = "normal")} else {
      flexsurvspline(form, k = k, scale = "normal", anc = list(gamma1 = ~group, gamma2 = ~group, gamma3 = ~group, gamma4 = ~group))
    }
    
    lbls_spline <- c(" 8. Spline 1 knot hazard", " 9. Spline 2 knots hazard", "10. Spline 3 knots hazard", "11. Spline 1 knot odds", 
                     "12. Spline 2 knots odds", "13. Spline 3 knots odds", "14. Spline 1 knot normal", "15. Spline 2 knots normal", 
                     "16. Spline 3 knots normal")
    
    # calculate AIC and BIC
    AIC_spl <- c(spl_hazard1$AIC, spl_hazard2$AIC, spl_hazard3$AIC, spl_odds1$AIC, spl_odds2$AIC, spl_odds3$AIC, spl_normal1$AIC, 
                 spl_normal2$AIC, spl_normal3$AIC)
    BIC_spl <- BIC(spl_hazard1, spl_hazard2, spl_hazard3, spl_odds1, spl_odds2, spl_odds3, spl_normal1, spl_normal2, spl_normal3)[2]
    IC_spl <- data.frame(lbls_spline, AIC_spl, BIC_spl)
    colnames(IC_spl) <- c("Model", "AIC", "BIC")
  }
  
  # fit cure models
  if (cure_mod == TRUE) {
    cure_weib_mix <- lapply(c(1:ngroups), function(x)
      flexsurvcure(Surv(years, status) ~ 1, data = data.frame(years, status)[group == levels(group)[x],], 
                   link = cure_link, dist = "weibull", mixture = TRUE)
    )
    cure_weib_nmix <- lapply(c(1:ngroups), function(x)
      flexsurvcure(Surv(years, status) ~ 1, data = data.frame(years, status)[group == levels(group)[x],], 
                   link = cure_link, dist = "weibull", mixture = FALSE)
    )
    cure_lnorm_mix <- lapply(c(1:ngroups), function(x)
      flexsurvcure(Surv(years, status) ~ 1, data = data.frame(years, status)[group == levels(group)[x],], 
                   link = cure_link, dist = "lnorm", mixture = TRUE)
    )
    cure_lnorm_nmix <- lapply(c(1:ngroups), function(x)
      flexsurvcure(Surv(years, status) ~ 1, data = data.frame(years, status)[group == levels(group)[x],], 
                   link = cure_link, dist = "lnorm", mixture = FALSE)
    )
    cure_llog_mix <- lapply(c(1:ngroups), function(x)
      flexsurvcure(Surv(years, status) ~ 1, data = data.frame(years, status)[group == levels(group)[x],], 
                   link = cure_link, dist = "llogis", mixture = TRUE)
    )
    cure_llog_nmix <- lapply(c(1:ngroups), function(x)
      flexsurvcure(Surv(years, status) ~ 1, data = data.frame(years, status)[group == levels(group)[x],], 
                   link = cure_link, dist = "llogis", mixture = FALSE)
    )
    
    lbls_cure <- c("17. Mixture cure Weibull", "18. Non-mixture cure Weibull", 
                   "19. Mixture cure Log-normal", "20. Non-mixture cure Log-normal", 
                   "21. Mixture cure Log-logistic", "22. Non-mixture cure Log-logistic") 
    
    # calculate AIC 
    AIC_cure <- c(sum(sapply(c(1:ngroups), function(x) cure_weib_mix[[x]]$AIC)),
                  sum(sapply(c(1:ngroups), function(x) cure_weib_nmix[[x]]$AIC)),
                  sum(sapply(c(1:ngroups), function(x) cure_lnorm_mix[[x]]$AIC)),
                  sum(sapply(c(1:ngroups), function(x) cure_lnorm_nmix[[x]]$AIC)),
                  sum(sapply(c(1:ngroups), function(x) cure_llog_mix[[x]]$AIC)),
                  sum(sapply(c(1:ngroups), function(x) cure_llog_nmix[[x]]$AIC)))
    cure_fraction <- rbind(sapply(c(1:ngroups), function(x) cure_weib_mix[[x]]$res[1,1:3]),
                           sapply(c(1:ngroups), function(x) cure_weib_nmix[[x]]$res[1,1:3]),
                           sapply(c(1:ngroups), function(x) cure_lnorm_mix[[x]]$res[1,1:3]),
                           sapply(c(1:ngroups), function(x) cure_lnorm_nmix[[x]]$res[1,1:3]),
                           sapply(c(1:ngroups), function(x) cure_llog_mix[[x]]$res[1,1:3]),
                           sapply(c(1:ngroups), function(x) cure_llog_nmix[[x]]$res[1,1:3])) 
    
    cure_fraction <- matrix(ncol = ngroups, data = 
                              sapply(c(1:ngroups), function(x) 
                                sapply(seq(1, 18, by = 3), function(y) 
                                  paste0(round(cure_fraction[y, x],3)*100, "% (",round(cure_fraction[y + 1, x],3)*100,"% - ",round(cure_fraction[y + 2, x],3)*100,"%)")
                                )
                              )
    )
    IC_cure <- data.frame(lbls_cure, AIC_cure, cure_fraction)
    colnames(IC_cure) <- c("Model", "AIC", paste0("Cure fraction ", group_names))
  }
  
  return(c(list(expo = expo, weib = weib, gom = gom, lnorm = lnorm, llog = llog, gam = gam, ggam = ggam, IC = IC), 
           if (spline_mod == TRUE) {list(spl_hazard1 = spl_hazard1, spl_hazard2 = spl_hazard2, spl_hazard3 = spl_hazard3, 
                                         spl_odds1 = spl_odds1, spl_odds2 = spl_odds2,  spl_odds3 = spl_odds3, 
                                         spl_normal1 = spl_normal1, spl_normal2 = spl_normal2, spl_normal3 = spl_normal3, 
                                         IC_spl = IC_spl)}, 
           if (cure_mod == TRUE) {list(cure_weib_mix = cure_weib_mix, cure_weib_nmix = cure_weib_nmix, 
                                       cure_lnorm_mix = cure_lnorm_mix, cure_lnorm_nmix = cure_lnorm_nmix, 
                                       cure_llog_mix = cure_llog_mix,  cure_llog_nmix = cure_llog_nmix, 
                                       IC_cure = IC_cure)})) 
}

# function to obtain predicted values from parametric survival models
f_surv_model_pred <- function(ngroups, time_pred, surv_model, spline_mod, cure_mod, group_names) { 
  
  column_names <- c("Time", group_names)
  
  expo_est <- summary(surv_model$expo, t = time_pred)
  weib_est <- summary(surv_model$weib, t = time_pred)
  gom_est <- summary(surv_model$gom, t = time_pred)
  lnorm_est <- summary(surv_model$lnorm, t = time_pred)
  llog_est <- summary(surv_model$llog, t = time_pred)
  gam_est <- summary(surv_model$gam, t = time_pred)
  ggam_est <- summary(surv_model$ggam, t = time_pred)
  
  expo_pred_h <- summary(surv_model$expo, t = time_pred, type = "hazard")
  weib_pred_h <- summary(surv_model$weib, t = time_pred, type = "hazard")
  gom_pred_h <- summary(surv_model$gom, t = time_pred, type = "hazard")
  lnorm_pred_h <- summary(surv_model$lnorm, t = time_pred, type = "hazard")
  llog_pred_h <- summary(surv_model$llog, t = time_pred, type = "hazard")
  gam_pred_h <- summary(surv_model$gam, t = time_pred, type = "hazard")
  ggam_pred_h <- summary(surv_model$ggam, t = time_pred, type = "hazard")
  
  if (ngroups == 1) {
    names(expo_est) <- names(weib_est) <- names(gom_est) <- names(lnorm_est) <- names(llog_est) <- 
      names(gam_est) <- names(ggam_est) <- paste("group=", group_names, sep = "")
  } # name has to be changed if there is only one group
  
  # extract prediction for each group
  expo_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
    expo_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est)) 
  weib_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
    weib_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
  gom_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
    gom_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
  lnorm_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
    lnorm_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
  llog_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
    llog_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
  gam_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
    gam_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
  ggam_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
    ggam_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
  
  gom_pred[, -1][gom_pred[, -1] < 1e-15] <- 0  # prevent rounding errors for predicted transition probabilities
  
  colnames(expo_pred) <- colnames(weib_pred) <- colnames(gom_pred) <- colnames(lnorm_pred) <- colnames(llog_pred) <- colnames(gam_pred) <- 
    colnames(ggam_pred) <- column_names
  
  if (spline_mod == TRUE) {
    spl_hazard1_est <- summary(surv_model$spl_hazard1, t = time_pred)
    spl_hazard2_est <- summary(surv_model$spl_hazard2, t = time_pred)
    spl_hazard3_est <- summary(surv_model$spl_hazard3, t = time_pred)
    spl_odds1_est <- summary(surv_model$spl_odds1, t = time_pred)
    spl_odds2_est <- summary(surv_model$spl_odds2, t = time_pred)
    spl_odds3_est <- summary(surv_model$spl_odds3, t = time_pred)
    spl_normal1_est <- summary(surv_model$spl_normal1, t = time_pred)
    spl_normal2_est <- summary(surv_model$spl_normal2, t = time_pred)
    spl_normal3_est <- summary(surv_model$spl_normal3, t = time_pred)
    
    spl_hazard1_pred_h <- summary(surv_model$spl_hazard1, t = time_pred, type = "hazard")
    spl_hazard2_pred_h <- summary(surv_model$spl_hazard2, t = time_pred, type = "hazard")
    spl_hazard3_pred_h <- summary(surv_model$spl_hazard3, t = time_pred, type = "hazard")
    spl_odds1_pred_h <- summary(surv_model$spl_odds1, t = time_pred, type = "hazard")
    spl_odds2_pred_h <- summary(surv_model$spl_odds2, t = time_pred, type = "hazard")
    spl_odds3_pred_h <- summary(surv_model$spl_odds3, t = time_pred, type = "hazard")
    spl_normal1_pred_h <- summary(surv_model$spl_normal1, t = time_pred, type = "hazard")
    spl_normal2_pred_h <- summary(surv_model$spl_normal2, t = time_pred, type = "hazard")
    spl_normal3_pred_h <- summary(surv_model$spl_normal3, t = time_pred, type = "hazard")
    
    if (ngroups == 1) {
      names(spl_hazard1_est) <- names(spl_hazard2_est) <- names(spl_hazard3_est) <- 
        names(spl_odds1_est) <- names(spl_odds2_est) <- names(spl_odds3_est) <- 
        names(spl_normal1_est) <- names(spl_normal2_est) <- names(spl_normal3_est) <- paste("group=", group_names, sep = "")
    } # name has to be changed if there is only one group
    
    # extract prediction for each group
    spl_hazard1_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      spl_hazard1_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    spl_hazard2_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      spl_hazard2_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    spl_hazard3_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      spl_hazard3_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    spl_odds1_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      spl_odds1_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    spl_odds2_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      spl_odds2_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    spl_odds3_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      spl_odds3_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    spl_normal1_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      spl_normal1_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    spl_normal2_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      spl_normal2_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    spl_normal3_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      spl_normal3_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    
    colnames(spl_hazard1_pred) <- colnames(spl_hazard2_pred) <- colnames(spl_hazard3_pred)<- colnames(spl_odds1_pred) <- 
      colnames(spl_odds2_pred) <- colnames(spl_odds3_pred) <- colnames(spl_normal1_pred) <- colnames(spl_normal2_pred) <- 
      colnames(spl_normal3_pred) <- column_names
  }
  
  if (cure_mod == TRUE) {## BR do this similarly for the other models (1 instead of 2 steps)
    cure_weib_mix_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_weib_mix[[x]], t = time_pred))[,2]))
    cure_weib_nmix_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_weib_nmix[[x]], t = time_pred))[,2]))
    cure_lnorm_mix_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_lnorm_mix[[x]], t = time_pred))[,2]))
    cure_lnorm_nmix_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_lnorm_nmix[[x]], t = time_pred))[,2]))
    cure_llog_mix_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_llog_mix[[x]], t = time_pred))[,2]))
    cure_llog_nmix_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_llog_nmix[[x]], t = time_pred))[,2]))
    
    cure_weib_mix_pred_h <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_weib_mix[[x]], t = time_pred, type = "hazard"))[,2]))
    cure_weib_nmix_pred_h <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_weib_nmix[[x]], t = time_pred, type = "hazard"))[,2]))
    cure_lnorm_mix_pred_h <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_lnorm_mix[[x]], t = time_pred, type = "hazard"))[,2]))
    cure_lnorm_nmix_pred_h <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_lnorm_nmix[[x]], t = time_pred, type = "hazard"))[,2]))
    cure_llog_mix_pred_h <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_llog_mix[[x]], t = time_pred, type = "hazard"))[,2]))
    cure_llog_nmix_pred_h <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
      data.frame(summary(surv_model$cure_llog_nmix[[x]], t = time_pred, type = "hazard"))[,2]))
    
    if (ngroups == 1) {
      names(cure_weib_mix_pred) <- names(cure_weib_nmix_pred) <- names(cure_lnorm_mix_pred) <- 
        names(cure_lnorm_nmix_pred) <- names(cure_llog_mix_pred) <- names(cure_llog_nmix_pred) <- 
        names(cure_weib_mix_pred_h) <- names(cure_weib_nmix_pred_h) <- names(cure_lnorm_mix_pred_h) <- 
        names(cure_lnorm_nmix_pred_h) <- names(cure_llog_mix_pred_h) <- names(cure_llog_nmix_pred_h) <- 
        paste("group=", group_names, sep = "")
    } # name has to be changed if there is only one group
    
    colnames(cure_weib_mix_pred) <- colnames(cure_weib_nmix_pred) <- colnames(cure_lnorm_mix_pred) <- 
      colnames(cure_lnorm_nmix_pred) <- colnames(cure_llog_mix_pred) <- colnames(cure_llog_nmix_pred) <- 
      colnames(cure_weib_mix_pred_h) <- colnames(cure_weib_nmix_pred_h) <- colnames(cure_lnorm_mix_pred_h) <- 
      colnames(cure_lnorm_nmix_pred_h) <- colnames(cure_llog_mix_pred_h) <- colnames(cure_llog_nmix_pred_h) <- 
      column_names
  }
  
  return(c(list(expo = expo_pred, weib = weib_pred, gom = gom_pred, lnorm = lnorm_pred, llog = llog_pred, 
                gam = gam_pred, ggam = ggam_pred, expo_h = expo_pred_h, weib_h = weib_pred_h, 
                gom_h = gom_pred_h, lnorm_h = lnorm_pred_h, llog_h = llog_pred_h, gam_h = gam_pred_h, 
                ggam_h = ggam_pred_h), 
           if (spline_mod == TRUE) {list(spl_hazard1 = spl_hazard1_pred, spl_hazard2 = spl_hazard2_pred, spl_hazard3 = spl_hazard3_pred, 
                                         spl_odds1 = spl_odds1_pred, spl_odds2 = spl_odds2_pred, spl_odds3 = spl_odds3_pred, 
                                         spl_normal1 = spl_normal1_pred, spl_normal2 = spl_normal2_pred, spl_normal3 = spl_normal3_pred, 
                                         spl_hazard1_h = spl_hazard1_pred_h, spl_hazard2_h = spl_hazard2_pred_h, spl_hazard3_h = spl_hazard3_pred_h, 
                                         spl_odds1_h = spl_odds1_pred_h, spl_odds2_h = spl_odds2_pred_h, spl_odds3_h = spl_odds3_pred_h,
                                         spl_normal1_h = spl_normal1_pred_h, spl_normal2_h = spl_normal2_pred_h, spl_normal3_h = spl_normal3_pred_h)},
           if (cure_mod == TRUE) {list(cure_weib_mix = cure_weib_mix_pred, cure_weib_nmix = cure_weib_nmix_pred, 
                                       cure_lnorm_mix = cure_lnorm_mix_pred, cure_lnorm_nmix = cure_lnorm_nmix_pred, 
                                       cure_llog_mix = cure_llog_mix_pred,  cure_llog_nmix = cure_llog_nmix_pred,
                                       cure_weib_mix_h = cure_weib_mix_pred_h, cure_weib_nmix_h = cure_weib_nmix_pred_h, 
                                       cure_lnorm_mix_h = cure_lnorm_mix_pred_h, cure_lnorm_nmix_h = cure_lnorm_nmix_pred_h, 
                                       cure_llog_mix_h = cure_llog_mix_pred_h,  cure_llog_nmix_h = cure_llog_nmix_pred_h)}))
}

# function to obtain predicted values from parametric survival models
f_surv_model_pred_gr <- function(ngroups, surv_model, surv_model_pred, spline_mod, cure_mod) { 
  
  surv_gr_1 <- cbind(surv_model_pred[[1]][,1],
                     sapply(c(1:7), function(x) surv_model_pred[[x]][,2]),
                     if (spline_mod == TRUE) {sapply(c(15:23), function(x) surv_model_pred[[x]][,2])},
                     if (cure_mod == TRUE) {sapply(c(33:38), function(x) surv_model_pred[[x]][,2])}
  )
  
  if (ngroups > 1) {
    surv_gr_2 <- cbind(surv_model_pred[[1]][,1],
                       sapply(c(1:7), function(x) surv_model_pred[[x]][,3]),
                       if (spline_mod == TRUE) {sapply(c(15:23), function(x) surv_model_pred[[x]][,3])},
                       if (cure_mod == TRUE) {sapply(c(33:38), function(x) surv_model_pred[[x]][,3])}
    )
  }
  
  if (ngroups > 2) {
    surv_gr_3 <- cbind(surv_model_pred[[1]][,1],
                       sapply(c(1:7), function(x) surv_model_pred[[x]][,4]),
                       if (spline_mod == TRUE) {sapply(c(15:23), function(x) surv_model_pred[[x]][,4])},
                       if (cure_mod == TRUE) {sapply(c(33:38), function(x) surv_model_pred[[x]][,4])}
    )
  }
  
  lbls_all <- c("time", surv_model$IC$Model, if (spline_mod == TRUE) {surv_model$IC_spl$Model}, if (cure_mod == TRUE) {surv_model$IC_cure$Model})
  
  colnames(surv_gr_1) <- lbls_all
  if (ngroups > 1) {colnames(surv_gr_2) <- lbls_all}
  if (ngroups > 2) {colnames(surv_gr_3) <- lbls_all}
  
  return(c(list(gr_1 = surv_gr_1), 
           if (ngroups > 1) {list(gr_2 = surv_gr_2)}, 
           if (ngroups > 2) {list(gr_3 = surv_gr_3)}))
}

# function to obtain predicted annual transition probability based on parametric survival models
f_surv_model_pred_tp_gr <- function(ngroups, time_pred, time_unit, surv_model_pred_gr, cols_tp) { 
  
  tp_gr_1 <- cbind(time_pred, 1 - (1 - (shift(surv_model_pred_gr$gr_1[, 2:cols_tp], 1L, type = "lag") - surv_model_pred_gr$gr_1[, 2:cols_tp])/
                                     shift(surv_model_pred_gr$gr_1[, 2:cols_tp], 1L, type = "lag"))^(1/time_unit))[-1, ]
  
  if (ngroups > 1) {
    tp_gr_2 <- cbind(time_pred, 1 - (1 - (shift(surv_model_pred_gr$gr_2[, 2:cols_tp], 1L, type = "lag") - surv_model_pred_gr$gr_2[, 2:cols_tp])/
                                       shift(surv_model_pred_gr$gr_2[, 2:cols_tp], 1L, type = "lag"))^(1/time_unit))[-1, ]
  }
  
  if (ngroups > 2) {
    tp_gr_3 <- cbind(time_pred, 1 - (1 - (shift(surv_model_pred_gr$gr_3[, 2:cols_tp], 1L, type = "lag") - surv_model_pred_gr$gr_3[, 2:cols_tp])/
                                       shift(surv_model_pred_gr$gr_3[, 2:cols_tp], 1L, type = "lag"))^(1/time_unit))[-1, ]
  }
  
  return(c(list(gr_1 = tp_gr_1), 
           if (ngroups > 1) {list(gr_2 = tp_gr_2)}, 
           if (ngroups > 2) {list(gr_3 = tp_gr_3)}))
}

# function tocreate output dataframe for use in MS Excel 
f_surv_model_excel <- function(ngroups, strata, surv_model, spline_mod, cure_mod) { # BR wellicht bestaande labels gebruiken? (surv_model$IC[,1] surv_model$IC_spl[,1] surv_model$IC_cure[,1])
  
  # distributions names
  distnames <- if (spline_mod == TRUE) {
    c(rep("1. Exponential", nrow(surv_model$expo$res.t)), rep("2. Weibull", nrow(surv_model$weib$res.t)), rep("3. Gompertz", nrow(surv_model$gom$res.t)), 
      rep("4. Log-normal", nrow(surv_model$lnorm$res.t)), rep("5. Log-logistic", nrow(surv_model$llog$res.t)), rep("6. Gamma", nrow(surv_model$gam$res.t)), 
      rep("7. Generalisedgamma", nrow(surv_model$ggam$res.t)), rep("8. 1-knot spline hazard", nrow(surv_model$spl_hazard1$res.t)), 
      rep("9. 1-knot spline odds", nrow(surv_model$spl_odds1$res.t)), rep("10. 1-knot spline normal", nrow(surv_model$spl_normal1$res.t)), 
      rep("11. 2-knot spline hazard", nrow(surv_model$spl_hazard2$res.t)), rep("12. 2-knot spline odds", nrow(surv_model$spl_odds2$res.t)), 
      rep("13. 2-knot spline normal", nrow(surv_model$spl_normal2$res.t)))} else { ##XP: overwegen 'else' op nieuwe lijn te zetten zodat het zichtbaarder is, ik moest even zoeken
        c(rep("1. Exponential", nrow(surv_model$expo$res.t)), rep("2. Weibull", nrow(surv_model$weib$res.t)), rep("3. Gompertz", nrow(surv_model$gom$res.t)), 
          rep("4. Log-normal", nrow(surv_model$lnorm$res.t)), rep("5. Log-logistic", nrow(surv_model$llog$res.t)), rep("6. Gamma", nrow(surv_model$gam$res.t)), 
          rep("7. Generalisedgamma", nrow(surv_model$ggam$res.t)))
      }
  
  # parameters' names
  parnames <- if (spline_mod == TRUE) {
    c(rownames(surv_model$expo$res.t), rownames(surv_model$weib$res.t), rownames(surv_model$gom$res.t), rownames(surv_model$lnorm$res.t), rownames(surv_model$llog$res.t), 
      rownames(surv_model$gam$res.t), rownames(surv_model$ggam$res.t), rownames(surv_model$spl_hazard1$res.t), rownames(surv_model$spl_odds1$res.t), 
      rownames(surv_model$spl_normal1$res.t), rownames(surv_model$spl_hazard2$res.t), rownames(surv_model$spl_odds2$res.t), rownames(surv_model$spl_normal2$res.t))} else { ##XP: idem
        c(rownames(surv_model$expo$res.t), rownames(surv_model$weib$res.t), rownames(surv_model$gom$res.t), rownames(surv_model$lnorm$res.t), rownames(surv_model$llog$res.t), 
          rownames(surv_model$gam$res.t), rownames(surv_model$ggam$res.t))
      }
  
  # extract parameters for each distribution
  res <- if (spline_mod == TRUE) {
    rbind(surv_model$expo$res.t, surv_model$weib$res.t, surv_model$gom$res.t, surv_model$lnorm$res.t, surv_model$llog$res.t, surv_model$gam$res.t, surv_model$ggam$res.t, surv_model$spl_hazard1$res.t, 
          surv_model$spl_odds1$res.t, surv_model$spl_normal1$res.t, surv_model$spl_hazard2$res.t, surv_model$spl_odds2$res.t, surv_model$spl_normal2$res.t)} else { ##XP: idem
            rbind(surv_model$expo$res.t, surv_model$weib$res.t, surv_model$gom$res.t, surv_model$lnorm$res.t, surv_model$llog$res.t, surv_model$gam$res.t, surv_model$ggam$res.t)
          }
  
  empty <- rep("", nrow(res))  # create vector of length of the parameters in order to separate the parameters from other outputs
  
  # compute number of additional rows
  addrows <- if (strata == FALSE) {1} else {ifelse(ngroups == 2, 2, 3)} ##XP: kan dit geautomatiseerd worden als we nog meer modellen willen toevoegen?
  addcols <- if (strata == FALSE) {1} else {ifelse(ngroups == 2, 2, 3)}
  
  # create covariance matrices of equal lengths
  cov <- if (spline_mod == TRUE) { ##XP: hier nog splines met 3 knots toevoegen?
    rbind(cbind(surv_model$expo$cov, 
                matrix(0, nrow = ngroups, ncol = if (strata == FALSE) {3} else {ifelse(ngroups == 2, 6, 9)})), 
          cbind(surv_model$weib$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
          cbind(surv_model$gom$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
          cbind(surv_model$lnorm$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
          cbind(surv_model$llog$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
          cbind(surv_model$gam$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
          cbind(surv_model$ggam$cov, matrix(0, nrow = ngroups + 2 * addrows, ncol = 1 * addcols)), 
          cbind(surv_model$spl_hazard1$cov, matrix(0, nrow = ngroups + 2 * addrows, ncol = 1 * addcols)), 
          cbind(surv_model$spl_odds1$cov, matrix(0, nrow = ngroups + 2 * addrows, ncol = 1 * addcols)), 
          cbind(surv_model$spl_normal1$cov, matrix(0, nrow = ngroups + 2 * addrows, ncol = 1 * addcols)), 
          surv_model$spl_hazard2$cov, surv_model$spl_odds2$cov, surv_model$spl_normal2$cov)} else {
            rbind(cbind(surv_model$expo$cov, matrix(0, nrow = ngroups, ncol = if (strata == FALSE) {3} else {
              ifelse(ngroups == 2, 6, 9)})), 
              cbind(surv_model$weib$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
              cbind(surv_model$gom$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
              cbind(surv_model$lnorm$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
              cbind(surv_model$llog$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
              cbind(surv_model$gam$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
              cbind(surv_model$ggam$cov, matrix(0, nrow = ngroups + 2 * addrows, ncol = 1 * addcols)))
          }
  
  survmod <- cbind(distnames, parnames, res, empty, empty, empty, cov)
  
  # rename the columns and transpose it
  colnames(survmod) <- c("Distnames", "Parnames", colnames(surv_model$expo$res.t), "", "Knots", "Cov_matrix", 
                         c(1:if (spline_mod == TRUE) {ncol(surv_model$spl_hazard2$cov)} else {
                           abs(ncol(surv_model$ggam$cov) + 
                                 abs(length(colnames(survmod)) - 
                                       length(c("Distnames", "Parnames", colnames(surv_model$expo$res.t), "", "Knots", "Cov_matrix", c(1:ncol(surv_model$ggam$cov))))))
                         }))
  
  survmod <- t(survmod)
  
  # add the knots
  if (spline_mod == TRUE) {##XP: hier nog splines met 3 knots toevoegen?
    # 1-knot splines
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "8. 1-knot spline hazard" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "9. 1-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "10. 1-knot spline normal") & 
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma0")] <- c(surv_model$spl_hazard1$knots[1], surv_model$spl_odds1$knots[1], 
                                                                                         surv_model$spl_normal1$knots[1])
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "8. 1-knot spline hazard" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "9. 1-knot spline odds" |
                     survmod[which(rownames(survmod) == "Distnames"), ] == "10. 1-knot spline normal") & 
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma1")] <- c(surv_model$spl_hazard1$knots[2], surv_model$spl_odds1$knots[2], 
                                                                                         surv_model$spl_normal1$knots[2])
    
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "8. 1-knot spline hazard" |
                     survmod[which(rownames(survmod) == "Distnames"), ] == "9. 1-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "10. 1-knot spline normal") &
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma2")] <- c(surv_model$spl_hazard1$knots[3], surv_model$spl_odds1$knots[3], 
                                                                                         surv_model$spl_normal1$knots[3])
    # 2-knot splines
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "11. 2-knot spline hazard" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "12. 2-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "13. 2-knot spline normal") & 
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma0")] <- c(surv_model$spl_hazard2$knots[1], surv_model$spl_odds2$knots[1], 
                                                                                         surv_model$spl_normal2$knots[1])
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "11. 2-knot spline hazard" |
                     survmod[which(rownames(survmod) == "Distnames"), ] == "12. 2-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "13. 2-knot spline normal") & 
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma1")] <- c(surv_model$spl_hazard2$knots[2], surv_model$spl_odds2$knots[2], 
                                                                                         surv_model$spl_normal2$knots[2])
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "11. 2-knot spline hazard" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "12. 2-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "13. 2-knot spline normal") & 
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma2")] <- c(surv_model$spl_hazard2$knots[3], surv_model$spl_odds2$knots[3], 
                                                                                         surv_model$spl_normal2$knots[3])
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "11. 2-knot spline hazard" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "12. 2-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "13. 2-knot spline normal") &
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma3")] <- c(surv_model$spl_hazard2$knots[4], surv_model$spl_odds2$knots[4], 
                                                                                         surv_model$spl_normal2$knots[4])
  }
  
  # remove column names
  colnames(survmod) <- c("Time-to-event models parameters", rep("", abs(ncol(survmod)) - 1))
  
  return(survmod)
}
