# code is formatted using formatR::tidy_source(width.cutoff = 100)
f_PERSUADE <- function(years, status, group, strata = FALSE, spline_mod = FALSE, 
                       time_unit, time_horizon, time_pred_surv_table) {
  
  # number of groups
  group <- droplevels(group)  #drop unused levels
  ngroups <- length(levels(group))  #number of groups
  group_names <- levels(group)  #name groups
  
  if (ngroups == 1) {
    strata <- FALSE
  }
  if (spline_mod == TRUE) {
    show_spline <- TRUE
  } else {
    show_spline <- FALSE
  }
  
  # time variables
  time_horizon <- max(time_pred_surv_table, time_horizon)
  time_pred_surv_table <- time_pred_surv_table/time_unit
  time_pred <- seq(from = 0, to = time_horizon, by = time_unit)
  
  # form
  if (ngroups == 1) {
    form <- as.formula(Surv(years, status) ~ 1)
  } else {
    form <- as.formula(Surv(years, status) ~ group)
  }
  
  # km
  km <- npsurv(form)
  km_names <- if (ngroups == 1) {
    rep(1, length(km$time))
  } else {
    if (ngroups >= 2) {
      c(rep(1, km$strata[1]), rep(2, km$strata[2]), if (ngroups == 3) {
        rep(3, km$strata[3])
      })
    }
  }
  
  # hazard rate
  hr_smooth1 <- muhaz(years, status, group == levels(group)[1])
  if (ngroups > 1) {
    hr_smooth2 <- muhaz(years, status, group == levels(group)[2])
  }
  if (ngroups > 2) {
    hr_smooth3 <- muhaz(years, status, group == levels(group)[3])
  }
  hr_names <- c(rep(1, length(hr_smooth1$est.grid)), if (ngroups > 1) {
    rep(2, length(hr_smooth2$est.grid))
  } else {
    NA
  }, if (ngroups > 2) {
    rep(3, length(hr_smooth3$est.grid))
  })
  hr_max <- data.frame(time = ceiling(max(hr_smooth1$est.grid, if (ngroups > 1) {
    hr_smooth2$est.grid
  } else {
    NA
  }, if (ngroups > 2) {
    hr_smooth3$est.grid
  } else {
    NA
  }, na.rm = TRUE)), hr_smooth = ceiling(max(hr_smooth1$haz.est, if (ngroups > 1) {
    hr_smooth2$haz.est
  } else {
    NA
  }, if (ngroups > 2) {
    hr_smooth3$haz.est
  } else {
    NA
  }, na.rm = TRUE)))
  
  # cumulative hazard
  cum_haz <- data.frame(group = c(rep(1, length(time_pred[time_pred <= max(years[group == levels(group)[1]])])), 
                                  if (ngroups > 1) {
                                    rep(2, length(time_pred[time_pred <= max(years[group == levels(group)[2]])]))
                                  }, if (ngroups > 2) {
                                    rep(3, length(time_pred[time_pred <= max(years[group == levels(group)[3]])]))
                                  }))
  cum_haz$time <- c(time_pred[0:length(cum_haz$group[cum_haz$group == 1])], if (ngroups > 1) {
    time_pred[0:length(cum_haz$group[cum_haz$group == 2])]
  }, if (ngroups > 2) {
    time_pred[0:length(cum_haz$group[cum_haz$group == 3])]
  })
  cum_haz$H <- c(estimateNAH(years[group == levels(group)[1]], status[group == levels(group)[1]])$H(cum_haz$time[cum_haz$group == 1]), 
                 if (ngroups > 1) {
                   estimateNAH(years[group == levels(group)[2]], status[group == levels(group)[2]])$H(cum_haz$time[cum_haz$group == 2])
                 }, if (ngroups > 2) {
                   estimateNAH(years[group == levels(group)[3]], status[group == levels(group)[3]])$H(cum_haz$time[cum_haz$group == 3])
                 })
  cum_haz$var <- c(estimateNAH(years[group == levels(group)[1]], status[group == levels(group)[1]])$Var(cum_haz$time[cum_haz$group == 1]), 
                   if (ngroups > 1) {
                     estimateNAH(years[group == levels(group)[2]], status[group == levels(group)[2]])$Var(cum_haz$time[cum_haz$group == 2])
                   }, if (ngroups > 2) {
                     estimateNAH(years[group == levels(group)[3]], status[group == levels(group)[3]])$Var(cum_haz$time[cum_haz$group == 3])
                   })
  cum_haz$H_upper <- pmax(0, cum_haz$H + sqrt(cum_haz$var) * qnorm(p = 0.975))
  cum_haz$H_lower <- pmax(0, cum_haz$H - sqrt(cum_haz$var) * qnorm(p = 0.975))
  cum_haz$H_delta <- c(shift(cum_haz$H[cum_haz$group == 1], 1L, type = "lag") - cum_haz$H[cum_haz$group == 1], 
                       if (ngroups > 1) {
                         shift(cum_haz$H[cum_haz$group == 2], 1L, type = "lag") - cum_haz$H[cum_haz$group == 2]
                       }, if (ngroups > 2) {
                         shift(cum_haz$H[cum_haz$group == 3], 1L, type = "lag") - cum_haz$H[cum_haz$group == 3]
                       })
  cum_haz$H_upper_delta <- c(shift(cum_haz$H_upper[cum_haz$group == 1], 1L, type = "lag") - cum_haz$H_upper[cum_haz$group == 1], 
                             if (ngroups > 1) {
                               shift(cum_haz$H_upper[cum_haz$group == 2], 1L, type = "lag") - cum_haz$H_upper[cum_haz$group == 2]
                             }, if (ngroups > 2) {
                               shift(cum_haz$H_upper[cum_haz$group == 3], 1L, type = "lag") - cum_haz$H_upper[cum_haz$group == 3]
                             })
  cum_haz$H_lower_delta <- c(shift(cum_haz$H_lower[cum_haz$group == 1], 1L, type = "lag") - cum_haz$H_lower[cum_haz$group == 1], 
                             if (ngroups > 1) {
                               shift(cum_haz$H_lower[cum_haz$group == 2], 1L, type = "lag") - cum_haz$H_lower[cum_haz$group == 2]
                             }, if (ngroups > 2) {
                               shift(cum_haz$H_lower[cum_haz$group == 3], 1L, type = "lag") - cum_haz$H_lower[cum_haz$group == 3]
                             })
  cum_haz$tp <- 1 - exp(cum_haz$H_delta)^(1/time_unit)
  cum_haz$tp_lower <- 1 - exp(cum_haz$H_lower_delta)^(1/time_unit)
  cum_haz$tp_upper <- 1 - exp(cum_haz$H_upper_delta)^(1/time_unit)
  cum_haz$tp_smooth <- pmax(0, pmin(1, c(NA, loess(cum_haz$tp[cum_haz$group == 1] ~ cum_haz$time[cum_haz$group ==  1])$fitted, 
                                         if (ngroups > 1) {
                                           c(NA, loess(cum_haz$tp[cum_haz$group == 2] ~ cum_haz$time[cum_haz$group == 2])$fitted)
                                         }, if (ngroups > 2) {
                                           c(NA, loess(cum_haz$tp[cum_haz$group == 3] ~ cum_haz$time[cum_haz$group == 3])$fitted)
                                         })))
  cum_haz$tp_upper_smooth <- pmax(0, pmin(1, c(NA, loess(cum_haz$tp_upper[cum_haz$group == 1] ~ cum_haz$time[cum_haz$group == 1])$fitted, 
                                               if (ngroups > 1) {
                                                 c(NA, loess(cum_haz$tp_upper[cum_haz$group == 2] ~ cum_haz$time[cum_haz$group == 2])$fitted)
                                               }, if (ngroups > 2) {
                                                 c(NA, loess(cum_haz$tp_upper[cum_haz$group == 3] ~ cum_haz$time[cum_haz$group == 3])$fitted)
                                               })))
  cum_haz$tp_lower_smooth <- pmax(0, pmin(1, c(NA, loess(cum_haz$tp_lower[cum_haz$group == 1] ~ cum_haz$time[cum_haz$group == 1])$fitted, 
                                               if (ngroups > 1) {
                                                 c(NA, loess(cum_haz$tp_lower[cum_haz$group == 2] ~ cum_haz$time[cum_haz$group == 2])$fitted)
                                               }, if (ngroups > 2) {
                                                 c(NA, loess(cum_haz$tp_lower[cum_haz$group == 3] ~ cum_haz$time[cum_haz$group == 3])$fitted)
                                               })))
  
  # cox (for Scaled Schoenfeld residuals)
  cox_reg <- coxph(form)
  
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
  IC <- IC
  
  # fit spline models
  if (spline_mod == TRUE) {
    k <- 1  #number of splines
    spl_hazard1 <- if (strata == FALSE) {
      flexsurvspline(form, k = k, scale = "hazard")
    } else {
      flexsurvspline(form, k = k, scale = "hazard", anc = list(gamma1 = ~group, gamma2 = ~group))
    }
    spl_odds1 <- if (strata == FALSE) {
      flexsurvspline(form, k = k, scale = "odds")
    } else {
      flexsurvspline(form, k = k, scale = "odds", anc = list(gamma1 = ~group, gamma2 = ~group))
    }
    spl_normal1 <- if (strata == FALSE) {
      flexsurvspline(form, k = k, scale = "normal")
    } else {
      flexsurvspline(form, k = k, scale = "normal", anc = list(gamma1 = ~group, gamma2 = ~group))
    }
    
    k <- 2  #number of splines
    spl_hazard2 <- if (strata == FALSE) {
      flexsurvspline(form, k = k, scale = "hazard")
    } else {
      flexsurvspline(form, k = k, scale = "hazard", anc = list(gamma1 = ~group, gamma2 = ~group, 
                                                               gamma3 = ~group))
    }
    spl_odds2 <- if (strata == FALSE) {
      flexsurvspline(form, k = k, scale = "odds")
    } else {
      flexsurvspline(form, k = k, scale = "odds", anc = list(gamma1 = ~group, gamma2 = ~group, 
                                                             gamma3 = ~group))
    }
    spl_normal2 <- if (strata == FALSE) {
      flexsurvspline(form, k = k, scale = "normal")
    } else {
      flexsurvspline(form, k = k, scale = "normal", anc = list(gamma1 = ~group, gamma2 = ~group, 
                                                               gamma3 = ~group))
    }
    
    lbls_spline <- c(" 8. Spline 1 knot hazard", " 9. Spline 2 knots hazard", "10. Spline 1 knot odds", 
                     "11. Spline 2 knots odds", "12. Spline 1 knot normal", "13. Spline 2 knots normal")
    
    # calculate AIC and BIC
    AIC_spl <- c(spl_hazard1$AIC, spl_hazard2$AIC, spl_odds1$AIC, spl_odds2$AIC, spl_normal1$AIC, 
                 spl_normal2$AIC)
    BIC_spl <- BIC(spl_hazard1, spl_hazard2, spl_odds1, spl_odds2, spl_normal1, spl_normal2)[2]
    IC_spl <- data.frame(lbls_spline, AIC_spl, BIC_spl)
    colnames(IC_spl) <- c("Model", "AIC", "BIC")
    IC_spl <- IC_spl
  }
  # predicted values parametric models
  column_names <- c("Time", group_names)
  
  expo_est <- summary(expo, t = time_pred)
  weib_est <- summary(weib, t = time_pred)
  gom_est <- summary(gom, t = time_pred)
  lnorm_est <- summary(lnorm, t = time_pred)
  llog_est <- summary(llog, t = time_pred)
  gam_est <- summary(gam, t = time_pred)
  ggam_est <- summary(ggam, t = time_pred)
  
  gom_pred_h <- summary(gom, t = time_pred, type = "hazard")  #used for diagnostics for Gompertz distribution
  
  if (ngroups == 1) {
    expo_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) expo_est[[1]]$est))
    weib_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) weib_est[[1]]$est))
    gom_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) gom_est[[1]]$est))
    lnorm_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) lnorm_est[[1]]$est))
    llog_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) llog_est[[1]]$est))
    gam_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) gam_est[[1]]$est))
    ggam_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) ggam_est[[1]]$est))
  }
  
  if (ngroups > 1) {
    expo_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) expo_est[[which(names(expo_est) == 
                                                                                     paste("group=", group_names[x], sep = ""))]]$est))
    weib_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) weib_est[[which(names(expo_est) == 
                                                                                     paste("group=", group_names[x], sep = ""))]]$est))
    gom_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) gom_est[[which(names(expo_est) == 
                                                                                   paste("group=", group_names[x], sep = ""))]]$est))
    lnorm_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) lnorm_est[[which(names(expo_est) == 
                                                                                       paste("group=", group_names[x], sep = ""))]]$est))
    llog_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) llog_est[[which(names(expo_est) == 
                                                                                     paste("group=", group_names[x], sep = ""))]]$est))
    gam_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) gam_est[[which(names(expo_est) == 
                                                                                   paste("group=", group_names[x], sep = ""))]]$est))
    ggam_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) ggam_est[[which(names(expo_est) == 
                                                                                     paste("group=", group_names[x], sep = ""))]]$est))
  }
  
  gom_pred[, -1][gom_pred[, -1] < 1e-15] <- 0  #prevent rounding errors for predicted transition probabilities
  
  colnames(expo_pred) <- colnames(weib_pred) <- colnames(gom_pred) <- colnames(lnorm_pred) <- colnames(llog_pred) <- colnames(gam_pred) <- 
    colnames(ggam_pred) <- column_names
  
  if (spline_mod == TRUE) {
    spl_hazard1_est <- summary(spl_hazard1, t = time_pred)
    spl_hazard2_est <- summary(spl_hazard2, t = time_pred)
    spl_odds1_est <- summary(spl_odds1, t = time_pred)
    spl_odds2_est <- summary(spl_odds2, t = time_pred)
    spl_normal1_est <- summary(spl_normal1, t = time_pred)
    spl_normal2_est <- summary(spl_normal2, t = time_pred)
    
    if (ngroups == 1) {
      spl_hazard1_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) spl_hazard1_est[[1]]$est))
      spl_hazard2_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) spl_hazard2_est[[1]]$est))
      spl_odds1_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) spl_odds1_est[[1]]$est))
      spl_odds2_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) spl_odds2_est[[1]]$est))
      spl_normal1_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) spl_normal1_est[[1]]$est))
      spl_normal2_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) spl_normal2_est[[1]]$est))
    }
    
    if (ngroups > 1) {
      spl_hazard1_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
        spl_hazard1_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
      spl_hazard2_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
        spl_hazard2_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
      spl_odds1_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
        spl_odds1_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
      spl_odds2_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
        spl_odds2_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
      spl_normal1_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
        spl_normal1_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
      spl_normal2_pred <- cbind(time_pred, sapply(c(1:ngroups), function(x) 
        spl_normal2_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    }
    
    colnames(spl_hazard1_pred) <- colnames(spl_hazard2_pred) <- colnames(spl_odds1_pred) <- colnames(spl_odds2_pred) <- 
      colnames(spl_normal1_pred) <- colnames(spl_normal2_pred) <- column_names
  }
  
  # predicted survival
  surv_gr_1 <- if (spline_mod == TRUE) {
    cbind(time_pred, expo_pred[, 2], weib_pred[, 2], gom_pred[, 2], lnorm_pred[, 2], llog_pred[, 2], 
          gam_pred[, 2], ggam_pred[, 2], spl_hazard1_pred[, 2], spl_hazard2_pred[, 2], spl_odds1_pred[, 2], 
          spl_odds2_pred[, 2], spl_normal1_pred[, 2], spl_normal2_pred[, 2])
  } else {
    cbind(time_pred, expo_pred[, 2], weib_pred[, 2], gom_pred[, 2], lnorm_pred[, 2], llog_pred[, 2],
          gam_pred[, 2], ggam_pred[, 2])
  }
  if (ngroups > 1) {
    surv_gr_2 <- if (spline_mod == TRUE) {
      cbind(time_pred, expo_pred[, 3], weib_pred[, 3], gom_pred[, 3], lnorm_pred[, 3], llog_pred[, 3],
            gam_pred[, 3], ggam_pred[, 3], spl_hazard1_pred[, 3], spl_hazard2_pred[, 3], spl_odds1_pred[, 3], 
            spl_odds2_pred[, 3], spl_normal1_pred[, 3], spl_normal2_pred[, 3])
    } else {
      cbind(time_pred, expo_pred[, 3], weib_pred[, 3], gom_pred[, 3], lnorm_pred[, 3], llog_pred[, 3],
            gam_pred[, 3], ggam_pred[, 3])
    }
  }
  
  if (ngroups > 2) {
    surv_gr_3 <- if (spline_mod == TRUE) {
      cbind(time_pred, expo_pred[, 4], weib_pred[, 4], gom_pred[, 4], lnorm_pred[, 4], llog_pred[, 4],
            gam_pred[, 4], ggam_pred[, 4], spl_hazard1_pred[, 4], spl_hazard2_pred[, 4], spl_odds1_pred[, 4],
            spl_odds2_pred[, 4], spl_normal1_pred[, 4], spl_normal2_pred[, 4])
    } else {
      cbind(time_pred, expo_pred[, 4], weib_pred[, 4], gom_pred[, 4], lnorm_pred[, 4], llog_pred[,4],
            gam_pred[, 4], ggam_pred[, 4])
    }
  }
  
  lbls_all <- if (spline_mod == TRUE) {
    c("Time", " 1. Exponential", " 2. Weibull", " 3. Gompertz", " 4. Log-normal", " 5. Log-logistic", 
      " 6. Gamma", " 7. Generalised Gamma", " 8. Spline 1 knot hazard", " 9. Spline 2 knots hazard", 
      "10. Spline 1 knot odds", "11. Spline 2 knots odds", "12. Spline 1 knot normal", "13. Spline 2 knots normal")
  } else {
    c("Time", "1. Exponential", "2. Weibull", "3. Gompertz", "4. Log-normal", "5. Log-logistic", "6. Gamma", 
      "7. Generalised Gamma")
  }
  
  colnames(surv_gr_1) <- lbls_all
  if (ngroups > 1) {
    colnames(surv_gr_2) <- lbls_all
  }
  if (ngroups > 2) {
    colnames(surv_gr_3) <- lbls_all
  }
  
  # calculate annual transition probability based on observed data (km)
  km_tp_gr_1 <- data.frame(
    time = cum_haz$time[cum_haz$group==1],
    tp_smooth = cum_haz$tp_smooth[cum_haz$group==1],
    tp_smooth_lower = cum_haz$tp_lower_smooth[cum_haz$group==1],
    tp_smooth_upper = cum_haz$tp_upper_smooth[cum_haz$group==1]
  )
  
  if (ngroups > 1) {
    km_tp_gr_2 <- data.frame(
      time = cum_haz$time[cum_haz$group==2],
      tp_smooth = cum_haz$tp_smooth[cum_haz$group==2],
      tp_smooth_lower = cum_haz$tp_lower_smooth[cum_haz$group==2],
      tp_smooth_upper = cum_haz$tp_upper_smooth[cum_haz$group==2]
    )
  }
  if (ngroups > 2) {
    km_tp_gr_3 <- data.frame(
      time = cum_haz$time[cum_haz$group==3],
      tp_smooth = cum_haz$tp_smooth[cum_haz$group==3],
      tp_smooth_lower = cum_haz$tp_lower_smooth[cum_haz$group==3],
      tp_smooth_upper = cum_haz$tp_upper_smooth[cum_haz$group==3]
    )
  }
  
  km_tp_max <- max(c(km_tp_gr_1$tp_smooth_upper, if (ngroups > 1) {
    km_tp_gr_2$tp_smooth_upper
  }, if (ngroups > 2) {
    km_tp_gr_3$tp_smooth_upper
  }), na.rm = TRUE)
  
  hr_names <- c(rep(1, length(hr_smooth1$est.grid)), if (ngroups > 1) {
    rep(2, length(hr_smooth2$est.grid))
  } else {
    NA
  }, if (ngroups > 2) {
    rep(3, length(hr_smooth3$est.grid))
  })
  
  # parametric survival models
  cols_extr <- ifelse(spline_mod == TRUE, 14, 8)  # define data frame width for the annual TP calculations
  
  tp_gr_1 <- cbind(time_pred, 1 - (1 - (shift(surv_gr_1[, 2:cols_extr], 1L, type = "lag") - surv_gr_1[, 2:cols_extr])/
                                     shift(surv_gr_1[, 2:cols_extr], 1L, type = "lag"))^(1/time_unit))[-1, ]
  
  if (ngroups > 1) {
    tp_gr_2 <- cbind(time_pred, 1 - (1 - (shift(surv_gr_2[, 2:cols_extr], 1L, type = "lag") - surv_gr_2[, 2:cols_extr])/
                                       shift(surv_gr_2[, 2:cols_extr], 1L, type = "lag"))^(1/time_unit))[-1, ]
  }
  
  if (ngroups > 2) {
    tp_gr_3 <- cbind(time_pred, 1 - (1 - (shift(surv_gr_3[, 2:cols_extr], 1L, type = "lag") - surv_gr_3[, 2:cols_extr])/
                                       shift(surv_gr_3[, 2:cols_extr], 1L, type = "lag"))^(1/time_unit))[-1, ]
  }
  
  # create output dataframe containing each distributions' name, the parameters' name, the parameters,
  # the covariance matrix and the knots (if applicable)
  
  # distributions names
  distnames <- if (spline_mod == TRUE) {
    c(rep("1. Exponential", nrow(expo$res.t)), rep("2. Weibull", nrow(weib$res.t)), rep("3. Gompertz", nrow(gom$res.t)), 
      rep("4. Log-normal", nrow(lnorm$res.t)), rep("5. Log-logistic", nrow(llog$res.t)), rep("6. Gamma", nrow(gam$res.t)), 
      rep("7. Generalisedgamma", nrow(ggam$res.t)), rep("8. 1-knot spline hazard", nrow(spl_hazard1$res.t)), 
      rep("9. 1-knot spline odds", nrow(spl_odds1$res.t)), rep("10. 1-knot spline normal", nrow(spl_normal1$res.t)), 
      rep("11. 2-knot spline hazard", nrow(spl_hazard2$res.t)), rep("12. 2-knot spline odds", nrow(spl_odds2$res.t)), 
      rep("13. 2-knot spline normal", nrow(spl_normal2$res.t)))
  } else {
    c(rep("1. Exponential", nrow(expo$res.t)), rep("2. Weibull", nrow(weib$res.t)), rep("3. Gompertz", nrow(gom$res.t)), 
      rep("4. Log-normal", nrow(lnorm$res.t)), rep("5. Log-logistic", nrow(llog$res.t)), rep("6. Gamma", nrow(gam$res.t)), 
      rep("7. Generalisedgamma", nrow(ggam$res.t)))
  }
  
  # parameters' names
  parnames <- if (spline_mod == TRUE) {
    c(rownames(expo$res.t), rownames(weib$res.t), rownames(gom$res.t), rownames(lnorm$res.t), rownames(llog$res.t), 
      rownames(gam$res.t), rownames(ggam$res.t), rownames(spl_hazard1$res.t), rownames(spl_odds1$res.t), 
      rownames(spl_normal1$res.t), rownames(spl_hazard2$res.t), rownames(spl_odds2$res.t), rownames(spl_normal2$res.t))
  } else {
    c(rownames(expo$res.t), rownames(weib$res.t), rownames(gom$res.t), rownames(lnorm$res.t), rownames(llog$res.t), 
      rownames(gam$res.t), rownames(ggam$res.t))
  }
  
  # extract parameters of each distribution
  res <- if (spline_mod == TRUE) {
    rbind(expo$res.t, weib$res.t, gom$res.t, lnorm$res.t, llog$res.t, gam$res.t, ggam$res.t, spl_hazard1$res.t, 
          spl_odds1$res.t, spl_normal1$res.t, spl_hazard2$res.t, spl_odds2$res.t, spl_normal2$res.t)
  } else {
    rbind(expo$res.t, weib$res.t, gom$res.t, lnorm$res.t, llog$res.t, gam$res.t, ggam$res.t)
  }
  
  empty <- rep("", nrow(res))  # create vector of length of the parameters in order to separate the parameters from other outputs
  
  # compute number of additional rows
  addrows <- if (strata == FALSE) {
    1
  } else {
    ifelse(ngroups == 2, 2, 3)
  }
  addcols <- if (strata == FALSE) {
    1
  } else {
    ifelse(ngroups == 2, 2, 3)
  }
  
  # create covariance matrices of equal lengths
  cov <- if (spline_mod == TRUE) {
    rbind(cbind(expo$cov, matrix(0, nrow = ngroups, ncol = if (strata == FALSE) {
      3
    } else {
      ifelse(ngroups == 2, 6, 9)
    })), cbind(weib$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
    cbind(gom$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
    cbind(lnorm$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
    cbind(llog$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
    cbind(gam$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
    cbind(ggam$cov, matrix(0, nrow = ngroups + 2 * addrows, ncol = 1 * addcols)), 
    cbind(spl_hazard1$cov, matrix(0, nrow = ngroups + 2 * addrows, ncol = 1 * addcols)), 
    cbind(spl_odds1$cov, matrix(0, nrow = ngroups + 2 * addrows, ncol = 1 * addcols)), 
    cbind(spl_normal1$cov, matrix(0, nrow = ngroups + 2 * addrows, ncol = 1 * addcols)), 
    spl_hazard2$cov, spl_odds2$cov, spl_normal2$cov)
  } else {
    rbind(cbind(expo$cov, matrix(0, nrow = ngroups, ncol = if (strata == FALSE) {
      3
    } else {
      ifelse(ngroups == 2, 6, 9)
    })), cbind(weib$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
    cbind(gom$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
    cbind(lnorm$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
    cbind(llog$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
    cbind(gam$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
    cbind(ggam$cov, matrix(0, nrow = ngroups + 2 * addrows, ncol = 1 * addcols)))
  }
  
  survmod <- cbind(distnames, parnames, res, empty, empty, empty, cov)
  
  # rename the columns and transpose it
  colnames(survmod) <- c("Distnames", "Parnames", colnames(expo$res.t), "", "Knots", "Cov_matrix", 
                         c(1:if (spline_mod == TRUE) {
                           ncol(spl_hazard2$cov)
                         } else {
                           abs(ncol(ggam$cov) + abs(length(colnames(survmod)) - length(c("Distnames", "Parnames", colnames(expo$res.t), 
                                                                                         "", "Knots", "Cov_matrix", c(1:ncol(ggam$cov))))))
                         }))
  
  survmod <- t(survmod)
  
  # add the knots
  if (spline_mod == TRUE) {
    # 1-knot splines
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "8. 1-knot spline hazard" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "9. 1-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "10. 1-knot spline normal") & 
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma0")] <- c(spl_hazard1$knots[1], spl_odds1$knots[1], 
                                                                                         spl_normal1$knots[1])
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "8. 1-knot spline hazard" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "9. 1-knot spline odds" |
                     survmod[which(rownames(survmod) == "Distnames"), ] == "10. 1-knot spline normal") & 
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma1")] <- c(spl_hazard1$knots[2], spl_odds1$knots[2], 
                                                                                         spl_normal1$knots[2])
    
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "8. 1-knot spline hazard" |
                     survmod[which(rownames(survmod) == "Distnames"), ] == "9. 1-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "10. 1-knot spline normal") &
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma2")] <- c(spl_hazard1$knots[3], spl_odds1$knots[3], 
                                                                                         spl_normal1$knots[3])
    # 2-knot splines
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "11. 2-knot spline hazard" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "12. 2-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "13. 2-knot spline normal") & 
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma0")] <- c(spl_hazard2$knots[1], spl_odds2$knots[1], 
                                                                                         spl_normal2$knots[1])
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "11. 2-knot spline hazard" |
                     survmod[which(rownames(survmod) == "Distnames"), ] == "12. 2-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "13. 2-knot spline normal") & 
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma1")] <- c(spl_hazard2$knots[2], spl_odds2$knots[2], 
                                                                                         spl_normal2$knots[2])
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "11. 2-knot spline hazard" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "12. 2-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "13. 2-knot spline normal") & 
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma2")] <- c(spl_hazard2$knots[3], spl_odds2$knots[3], 
                                                                                         spl_normal2$knots[3])
    survmod[which(rownames(survmod) == "Knots"), 
            which((survmod[which(rownames(survmod) == "Distnames"), ] == "11. 2-knot spline hazard" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "12. 2-knot spline odds" | 
                     survmod[which(rownames(survmod) == "Distnames"), ] == "13. 2-knot spline normal") &
                    survmod[which(rownames(survmod) == "Parnames"), ] == "gamma3")] <- c(spl_hazard2$knots[4], spl_odds2$knots[4], 
                                                                                         spl_normal2$knots[4])
  }
  
  # remove column names
  colnames(survmod) <- c("Time-to-event models parameters", rep("", abs(ncol(survmod)) - 1))
  
  # export to global environment **** NEEdS ADJUSTING **** cluster (e.g. input, km, hr, cum_h, models, model_pred, misc)
  l1 <- list(years = years, status = status, group = group, ngroups = ngroups, group_names = group_names, show_spline = show_spline, 
             time_horizon = time_horizon, time_pred_surv_table = time_pred_surv_table, time_pred = time_pred, form = form, 
             km = km, km_names = km_names, hr_smooth1 = hr_smooth1, hr_names = hr_names, hr_max = hr_max, cum_haz = cum_haz,
             cox_reg = cox_reg, expo = expo, weib = weib, gom = gom, lnorm = lnorm, llog = llog, gam = gam, ggam = ggam, 
             expo_pred = expo_pred, weib_pred = weib_pred, gom_pred = gom_pred, gom_pred_h = gom_pred_h, lnorm_pred = lnorm_pred, 
             llog_pred = llog_pred, gam_pred = gam_pred, ggam_pred = ggam_pred, lbls = lbls, IC = IC, survmod = survmod,
             surv_gr_1 = surv_gr_1, km_tp_max = km_tp_max, km_tp_gr_1 = km_tp_gr_1, tp_gr_1 = tp_gr_1, cols_extr = cols_extr)
  if (ngroups > 1) {
    l2 <- list(hr_smooth2 = hr_smooth2, surv_gr_2 = surv_gr_2, km_tp_gr_2 = km_tp_gr_2, 
               tp_gr_2 = tp_gr_2)
  } else {
    l2 <- NA
  }
  if (ngroups > 2) {
    l3 <- list(hr_smooth3 = hr_smooth3, surv_gr_3 = surv_gr_3, km_tp_gr_3 = km_tp_gr_3, 
               tp_gr_3 = tp_gr_3)
  } else {
    l3 <- NA
  }
  if (spline_mod == TRUE) {
    l4 <- list(spl_hazard1 = spl_hazard1, spl_hazard2 = spl_hazard2, spl_odds1 = spl_odds1, spl_odds2 = spl_odds2, 
               spl_normal1 = spl_normal1, spl_normal2 = spl_normal2, spl_hazard1_pred = spl_hazard1_pred, 
               spl_hazard2_pred = spl_hazard2_pred, spl_odds1_pred = spl_odds1_pred, spl_odds2_pred = spl_odds2_pred, 
               spl_normal1_pred = spl_normal1_pred, spl_normal2_pred = spl_normal2_pred, IC_spl = IC_spl, 
               lbls_spline = lbls_spline)
  } else {
    l4 <- NA
  }
  
  output <- c(l1, l2, l3, l4)
  output <- output[names(output) != ""]
  output <- output[order(names(output))]
  return(output)
}