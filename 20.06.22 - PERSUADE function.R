PERSUADE <- function(years, status, group, strata = FALSE, time_unit, 
                                time_horizon, time_pred_surv_table, spline_mod = FALSE, csv_semicolon=FALSE,
                                csv_comma=FALSE, clipboard=FALSE) {
  
  #number of groups
  group <<- droplevels(group)                #drop unused levels
  group <- droplevels(group)                 #drop unused levels
  ngroups <<- length(levels(group))          #number of groups
  group_names <<- levels(group)              #name groups
  
  if (ngroups==1) {strata <- FALSE}
  if (spline_mod == TRUE) {show_spline <<- TRUE} else {show_spline <<- FALSE}
  
  #time variables
  time_horizon <<- max(time_pred_surv_table, time_horizon)
  time_pred_surv_table <<- time_pred_surv_table / time_unit
  time_pred <<- seq(from = 0, to = time_horizon, by = time_unit)
  
  #form
  if (ngroups==1) {form <<- as.formula(Surv(years, status) ~ 1)} else {
    form <<- as.formula(Surv(years, status) ~ group)}
  
  #km
  km  <<- npsurv(form)
  
  #hr
  hr_smooth1  <<- muhaz(years, status, group==levels(group)[1])
  if (ngroups>1) {hr_smooth2  <<- muhaz(years, status, group==levels(group)[2])} 
  if (ngroups>2) {hr_smooth3  <<- muhaz(years, status, group==levels(group)[3])} 
  
  #cox (for Scaled Schoenfeld residuals)
  cox_reg  <<- coxph(form)
  
  #fit parameteric models
  expo  <<- flexsurvreg(form, dist="exp")
  
  weib  <<- if (strata == FALSE) flexsurvreg(form, dist="weibull") else {
    flexsurvreg(form, anc= list(shape= ~ group), dist="weibull")}
  
  gom  <<- if (strata == FALSE) flexsurvreg(form, dist="gompertz") else {
    flexsurvreg(form, anc= list(shape= ~ group), dist="gompertz")}
  
  lnorm  <<- if (strata == FALSE) flexsurvreg(form, dist="lnorm") else {
    flexsurvreg(form, anc= list(sdlog= ~ group),dist="lnorm")}
  
  llog  <<- if (strata == FALSE) flexsurvreg(form, dist="llogis") else {
    flexsurvreg(form, anc= list(shape= ~ group), dist="llogis")}
  
  gam  <<- if (strata == FALSE) flexsurvreg(form, dist="gamma") else {
    flexsurvreg(form, anc= list(shape= ~ group), dist="gamma")}
  
  ggam  <<- if (strata == FALSE) flexsurvreg(form, dist="gengamma") else {
    flexsurvreg(form, anc= list(sigma= ~ group, Q= ~ group ), dist="gengamma")}
  
  lbls  <<- c("Exponential", "Weibull", "Gompertz", "Lognormal", "Loglogistic", "Gamma", 
              "Generalised Gamma")
  
  #calculate AIC and BIC
  AIC <- c(expo$AIC, weib$AIC, gom$AIC, lnorm$AIC, llog$AIC, gam$AIC, ggam$AIC)
  BIC <- BIC(expo, weib, gom, lnorm, llog, gam, ggam)[2]
  IC <- data.frame(c("exp", "weib", "gom", "lnorm", "llog", "gam", "ggam"), AIC, BIC,
                   row.names=NULL)
  colnames(IC) <- c("Distribution", "AIC", "BIC")
  IC <<- IC
  
  #fit spline models
  if (spline_mod == TRUE){
    k <- 1 #number of splines
    spl_hazard1  <<- if (strata == FALSE) {flexsurvspline(form, k=k, scale="hazard")} else { 
      flexsurvspline(form, k=k, scale="hazard", anc= list(gamma1 = ~ group, gamma2 = ~ group))}
    
    spl_odds1  <<- if (strata == FALSE) {flexsurvspline(form, k=k, scale="odds")} else { 
      flexsurvspline(form, k=k, scale="odds", anc= list(gamma1 = ~ group, gamma2 = ~ group))}
    
    spl_normal1  <<- if (strata == FALSE) {flexsurvspline(form, k=k, scale="normal")} else { 
      flexsurvspline(form, k=k, scale="normal", anc= list(gamma1 = ~ group, gamma2 = ~ group))}     
    
    k <- 2 #number of splines
    spl_hazard2  <<- if (strata == FALSE) {flexsurvspline(form, k=k, scale="hazard")} else { 
      flexsurvspline(form, k=k, scale="hazard", anc= list(gamma1 = ~ group, gamma2 = ~ group, gamma3 = ~ group))}
    
    spl_odds2  <<- if (strata == FALSE) {flexsurvspline(form, k=k, scale="odds")} else { 
      flexsurvspline(form, k=k, scale="odds", anc= list(gamma1 = ~ group, gamma2 = ~ group, gamma3 = ~ group))}
    
    spl_normal2  <<- if (strata == FALSE) {flexsurvspline(form, k=k, scale="normal")} else { 
      flexsurvspline(form, k=k, scale="normal", anc= list(gamma1 = ~ group, gamma2 = ~ group, gamma3 = ~ group))} 
    
    lbls_spline  <<- c("Spline 1 knot hazard", "Spline 2 knots hazard", "Spline 1 knot odds", 
                       "Spline 2 knots odds", "Spline 1 knot normal", "Spline 2 knots normal")
    
    #calculate AIC and BIC
    AIC_spl <- c(spl_hazard1$AIC, spl_odds1$AIC, spl_normal1$AIC, spl_hazard2$AIC, spl_odds2$AIC, spl_normal2$AIC)
    BIC_spl <- BIC(spl_hazard1, spl_odds1, spl_normal1, spl_hazard2, spl_odds2, spl_normal2)[2]
    IC_spl <- data.frame(rep(c("hazard", "odds", "normal"),2), c(rep(1,3),rep(2,3)), AIC_spl, BIC_spl,
                         row.names=NULL)
    colnames(IC_spl) <- c("Scale", "knots", "AIC", "BIC")
    IC_spl <<- IC_spl
  }
  #predicted values parametric models
  column_names <- c("Time", group_names)

  expo_est  <<- summary(expo, t=time_pred)
  weib_est  <<- summary(weib, t=time_pred)
  gom_est  <<- summary(gom, t=time_pred)
  gom_est_h  <<- summary(gom, t=time_pred, type="hazard") #used for diagnostics for Gompertz distribution
  lnorm_est  <<- summary(lnorm, t=time_pred)
  llog_est  <<- summary(llog, t=time_pred)
  gam_est  <<- summary(gam, t=time_pred)
  ggam_est  <<- summary(ggam, t=time_pred)

  if (ngroups==1) {
    expo_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) expo_est[[1]]$est))
    weib_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) weib_est[[1]]$est))
    gom_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) gom_est[[1]]$est))
    lnorm_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) lnorm_est[[1]]$est))
    llog_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) llog_est[[1]]$est))
    gam_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) gam_est[[1]]$est))
    ggam_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) ggam_est[[1]]$est))
  }
  
  if (ngroups>1) {
    expo_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) expo_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    weib_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) weib_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    gom_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) gom_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    lnorm_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) lnorm_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    llog_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) llog_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    gam_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) gam_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    ggam_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) ggam_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
  }

  colnames(expo_pred) <- colnames(weib_pred) <- colnames(gom_pred) <- colnames(lnorm_pred) <- colnames(llog_pred) <-
    colnames(gam_pred) <- colnames(ggam_pred) <- column_names

  if (spline_mod == TRUE) {
    spl_hazard1_est <<- summary(spl_hazard1, t=time_pred)
    spl_hazard2_est <<- summary(spl_hazard2, t=time_pred)
    spl_odds1_est <<- summary(spl_odds1, t=time_pred)
    spl_odds2_est <<- summary(spl_odds2, t=time_pred)
    spl_normal1_est <<- summary(spl_normal1, t=time_pred)
    spl_normal2_est <<- summary(spl_normal2, t=time_pred)
    
    if (ngroups==1) {
      spl_hazard1_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_hazard1_est[[1]]$est))
      spl_hazard2_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_hazard2_est[[1]]$est))
      spl_odds1_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_odds1_est[[1]]$est))
      spl_odds2_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_odds2_est[[1]]$est))
      spl_normal1_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_normal1_est[[1]]$est))
      spl_normal2_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_normal2_est[[1]]$est))
    }
    
    if (ngroups>1) {
      spl_hazard1_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_hazard1_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
      spl_hazard2_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_hazard2_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
      spl_odds1_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_odds1_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
      spl_odds2_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_odds2_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
      spl_normal1_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_normal1_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
      spl_normal2_pred <<- cbind(time_pred, sapply(c(1:ngroups), function (x) spl_normal2_est[[which(names(expo_est) == paste("group=", group_names[x], sep = ""))]]$est))
    }
    
    colnames(spl_hazard1_pred) <- colnames(spl_hazard2_pred) <- colnames(spl_odds1_pred) <- colnames(spl_odds2_pred) <-
      colnames(spl_normal1_pred) <- colnames(spl_normal2_pred) <- column_names
    
  }
  #predicted survival
  extrapolation_gr_1 <- if (spline_mod == TRUE){cbind(time_pred, expo_pred[,2], weib_pred[,2], gom_pred[,2], lnorm_pred[,2], llog_pred[,2],
                                                      gam_pred[,2], ggam_pred[,2], spl_hazard1_pred[,2], spl_hazard2_pred[,2], 
                                                      spl_odds1_pred[,2], spl_odds2_pred[,2], spl_normal1_pred[,2], 
                                                      spl_normal2_pred[,2])} else {
                                                        cbind(time_pred, expo_pred[,2], weib_pred[,2], gom_pred[,2], lnorm_pred[,2], llog_pred[,2],
                                                              gam_pred[,2], ggam_pred[,2])}
  if (ngroups>1) {
    extrapolation_gr_2 <- if (spline_mod == TRUE){cbind(time_pred, expo_pred[,3], weib_pred[,3], gom_pred[,3], lnorm_pred[,3], llog_pred[,3], 
                                                        gam_pred[,3], ggam_pred[,3], spl_hazard1_pred[,3], spl_hazard2_pred[,3],
                                                        spl_odds1_pred[,3], spl_odds2_pred[,3], spl_normal1_pred[,3], spl_normal2_pred[,3])} else{
                                                          cbind(time_pred, expo_pred[,3], weib_pred[,3], gom_pred[,3], lnorm_pred[,3], llog_pred[,3], 
                                                                gam_pred[,3], ggam_pred[,3])}
  }
  
  if (ngroups>2) {
    extrapolation_gr_3 <- if (spline_mod == TRUE){cbind(time_pred, expo_pred[,4], weib_pred[,4], gom_pred[,4], lnorm_pred[,4], llog_pred[,4], 
                                                        gam_pred[,4], ggam_pred[,4], spl_hazard1_pred[,4], spl_hazard2_pred[,4], 
                                                        spl_odds1_pred[,4], spl_odds2_pred[,4], spl_normal1_pred[,4],
                                                        spl_normal2_pred[,4])} else {
                                                          cbind(time_pred, expo_pred[,4], weib_pred[,4], gom_pred[,4], lnorm_pred[,4], llog_pred[,4], 
                                                                gam_pred[,4], ggam_pred[,4])}
  }
  
  lbls_all <- if (spline_mod == TRUE){c("Time", "Exponential", "Weibull", "Gompertz", "Lognormal", "Loglogistic", 
                                        "Gamma", "Generalised Gamma", "Spline 1 knot hazard", "Spline 2 knots hazard",
                                        "Spline 1 knot odds", "Spline 2 knots odds", "Spline 1 knot normal",
                                        "Spline 2 knots normal")} else {
                                          c("Time", "Exponential", "Weibull", "Gompertz", "Lognormal", "Loglogistic", 
                                            "Gamma", "Generalised Gamma")
                                        }
  
  colnames(extrapolation_gr_1) <- lbls_all
  if (ngroups>1) {colnames(extrapolation_gr_2) <- lbls_all}
  if (ngroups>2) {colnames(extrapolation_gr_3) <- lbls_all}
  
  extrapolation_gr_1  <<- extrapolation_gr_1
  if (ngroups>1) {extrapolation_gr_2  <<- extrapolation_gr_2}
  if (ngroups>2) {extrapolation_gr_3  <<- extrapolation_gr_3}
  
  #calculate annual transition probability
  cols_extr <<- ifelse(spline_mod == TRUE, 14, 8) # define data frame width for the annual TP calculations, based on whether spline models are asked
  
  tp_gr_1  <<- cbind(time_pred, 1 - (1 - 
                                       (shift(extrapolation_gr_1[, 2:cols_extr], 1L, type="lag") - extrapolation_gr_1[, 2:cols_extr]) 
                                     / shift(extrapolation_gr_1[, 2:cols_extr], 1L, type="lag")) ^ (1 / time_unit))[-1,]
  
  if (ngroups>1) {
    tp_gr_2  <<- cbind(time_pred, 1 - (1 - 
                                         (shift(extrapolation_gr_2[, 2:cols_extr], 1L, type="lag") - extrapolation_gr_2[, 2:cols_extr]) 
                                       / shift(extrapolation_gr_2[, 2:cols_extr], 1L, type="lag")) ^ (1 / time_unit))[-1,]}
  
  if (ngroups>2) {
    tp_gr_3  <<- cbind(time_pred, 1 - (1 - 
                                         (shift(extrapolation_gr_3[, 2:cols_extr], 1L, type="lag") - extrapolation_gr_3[, 2:cols_extr]) 
                                       / shift(extrapolation_gr_3[, 2:cols_extr], 1L, type="lag")) ^ (1 / time_unit))[-1,]}
  
  #create output dataframe containing each distributions' name, the parameters' name, the parameters, the covariance matrix and the knots (if applicable)
  
  #distributions names
  distnames <- if (spline_mod == TRUE) {c(rep("Exponential", nrow(expo$res.t)), rep("Weibull",nrow(weib$res.t)), 
                                          rep("Gompertz",nrow(gom$res.t)), rep("Lognormal",nrow(lnorm$res.t)),
                                          rep("Loglogistic",nrow(llog$res.t)), rep("Gamma",nrow(gam$res.t)),
                                          rep("Generalisedgamma",nrow(ggam$res.t)), rep("1-knot spline hazard", nrow(spl_hazard1$res.t)),
                                          rep("1-knot spline odds", nrow(spl_odds1$res.t)), 
                                          rep("1-knot spline normal", nrow(spl_normal1$res.t)),
                                          rep("2-knot spline hazard", nrow(spl_hazard2$res.t)), 
                                          rep("2-knot spline odds", nrow(spl_odds2$res.t)),
                                          rep("2-knot spline normal", nrow(spl_normal2$res.t)))} else {
                                            c(rep("Exponential", nrow(expo$res.t)), rep("Weibull",nrow(weib$res.t)), 
                                              rep("Gompertz",nrow(gom$res.t)), rep("Lognormal",nrow(lnorm$res.t)),
                                              rep("Loglogistic",nrow(llog$res.t)), rep("Gamma",nrow(gam$res.t)),
                                              rep("Generalisedgamma",nrow(ggam$res.t)))
                                          }
  
  #parameters' names
  parnames <- if (spline_mod == TRUE) {c(rownames(expo$res.t), rownames(weib$res.t), rownames(gom$res.t), rownames(lnorm$res.t), 
                                         rownames(llog$res.t), rownames(gam$res.t), rownames(ggam$res.t), rownames(spl_hazard1$res.t),
                                         rownames(spl_odds1$res.t), rownames(spl_normal1$res.t), rownames(spl_hazard2$res.t), 
                                         rownames(spl_odds2$res.t), rownames(spl_normal2$res.t))} else {
                                           c(rownames(expo$res.t), rownames(weib$res.t), rownames(gom$res.t), rownames(lnorm$res.t), 
                                             rownames(llog$res.t), rownames(gam$res.t), rownames(ggam$res.t))
                                         }
  
  #extract parameters of each distribution
  res <- if (spline_mod == TRUE) {rbind(expo$res.t, weib$res.t, gom$res.t, lnorm$res.t, llog$res.t, gam$res.t, ggam$res.t, spl_hazard1$res.t,
                                        spl_odds1$res.t, spl_normal1$res.t, spl_hazard2$res.t, spl_odds2$res.t, spl_normal2$res.t)} else {
                                          rbind(expo$res.t, weib$res.t, gom$res.t, lnorm$res.t, llog$res.t, gam$res.t, ggam$res.t)
                                        }
  
  empty <- rep("", nrow(res)) # create vector of length of the parameters in order to separate the parameters from other outputs
  
  #compute number of additional rows
  addrows <- if (strata == FALSE) {1} else {ifelse(ngroups == 2, 2, 3)}
  addcols <- if (strata == FALSE) {1} else {ifelse(ngroups == 2, 2, 3)}
  
  #create covariance matrices of equal lengths
  cov <- if (spline_mod == TRUE) {rbind(cbind(expo$cov, matrix(0, nrow = ngroups, ncol = if (strata == FALSE) {3} else {ifelse(ngroups == 2, 6, 9)})), 
                                        cbind(weib$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
                                        cbind(gom$cov,  matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)),  
                                        cbind(lnorm$cov,  matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)),
                                        cbind(llog$cov,  matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
                                        cbind(gam$cov,  matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
                                        cbind(ggam$cov,  matrix(0, nrow = ngroups + 2 * addrows, ncol = 1* addcols)), 
                                        cbind(spl_hazard1$cov,  matrix(0, nrow = ngroups + 2 * addrows, ncol = 1* addcols)), 
                                        cbind(spl_odds1$cov,  matrix(0, nrow = ngroups + 2 * addrows, ncol = 1* addcols)), 
                                        cbind(spl_normal1$cov,  matrix(0, nrow = ngroups + 2 * addrows, ncol = 1* addcols)), 
                                        spl_hazard2$cov, spl_odds2$cov, spl_normal2$cov)} else { # THIS ONE NEEDS TO BE CHECKED BECAUSE I THINK THIS IS NOT GOING TO GO WELL
                                          rbind(cbind(expo$cov, matrix(0, nrow = ngroups, ncol = if (strata == FALSE) {3} else {ifelse(ngroups == 2, 6, 9)})), 
                                                cbind(weib$cov, matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
                                                cbind(gom$cov,  matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)),  
                                                cbind(lnorm$cov,  matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)),
                                                cbind(llog$cov,  matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
                                                cbind(gam$cov,  matrix(0, nrow = ngroups + 1 * addrows, ncol = 2 * addcols)), 
                                                cbind(ggam$cov,  matrix(0, nrow = ngroups + 2 * addrows, ncol = 1* addcols)))
                                          
                                        }
  
  Survmod <- cbind(distnames, parnames, res, empty, empty, empty, cov)
  
  #rename the columns and transpose it
  colnames(Survmod) <- c("Distnames","Parnames",colnames(expo$res.t),"","Knots", 
                         "Cov_matrix",c(1:if (spline_mod == TRUE) {ncol(spl_hazard2$cov)} else {
                           abs(ncol(ggam$cov) + abs(length(colnames(Survmod)) - length(c("Distnames","Parnames",colnames(expo$res.t),"","Knots", 
                                                                                         "Cov_matrix",c(1:ncol(ggam$cov))))))
                         })) 
  
  Survmod <- t(Survmod)
  
  #add the knots
  if (spline_mod == TRUE) {
    #1-knot splines
    Survmod[which(rownames(Survmod)=="Knots"), 
            which((Survmod[which(rownames(Survmod)=="Distnames"),]=="1-knot spline hazard"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="1-knot spline odds"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="1-knot spline normal")&
                    Survmod[which(rownames(Survmod)=="Parnames"),]=="gamma0")] <- c(spl_hazard1$knots[1], 
                                                                                    spl_odds1$knots[1],spl_normal1$knots[1])
    Survmod[which(rownames(Survmod)=="Knots"), 
            which((Survmod[which(rownames(Survmod)=="Distnames"),]=="1-knot spline hazard"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="1-knot spline odds"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="1-knot spline normal")&
                    Survmod[which(rownames(Survmod)=="Parnames"),]=="gamma1")] <- c(spl_hazard1$knots[2], 
                                                                                    spl_odds1$knots[2],spl_normal1$knots[2])
    
    Survmod[which(rownames(Survmod)=="Knots"), 
            which((Survmod[which(rownames(Survmod)=="Distnames"),]=="1-knot spline hazard"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="1-knot spline odds"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="1-knot spline normal")&
                    Survmod[which(rownames(Survmod)=="Parnames"),]=="gamma2")] <- c(spl_hazard1$knots[3], 
                                                                                    spl_odds1$knots[3],spl_normal1$knots[3])
    #2-knot splines
    Survmod[which(rownames(Survmod)=="Knots"), 
            which((Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline hazard"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline odds"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline normal")&
                    Survmod[which(rownames(Survmod)=="Parnames"),]=="gamma0")] <- c(spl_hazard2$knots[1], 
                                                                                    spl_odds2$knots[1],spl_normal2$knots[1])
    Survmod[which(rownames(Survmod)=="Knots"), 
            which((Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline hazard"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline odds"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline normal")&
                    Survmod[which(rownames(Survmod)=="Parnames"),]=="gamma1")] <- c(spl_hazard2$knots[2], 
                                                                                    spl_odds2$knots[2],spl_normal2$knots[2])
    Survmod[which(rownames(Survmod)=="Knots"), 
            which((Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline hazard"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline odds"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline normal")&
                    Survmod[which(rownames(Survmod)=="Parnames"),]=="gamma2")] <- c(spl_hazard2$knots[3], 
                                                                                    spl_odds2$knots[3],spl_normal2$knots[3])
    Survmod[which(rownames(Survmod)=="Knots"), 
            which((Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline hazard"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline odds"|
                     Survmod[which(rownames(Survmod)=="Distnames"),]=="2-knot spline normal")&
                    Survmod[which(rownames(Survmod)=="Parnames"),]=="gamma3")] <- c(spl_hazard2$knots[4], 
                                                                                    spl_odds2$knots[4],spl_normal2$knots[4])
  }
  #remove column names
  colnames(Survmod) <- c("Time-to-event models parameters",rep("", abs(ncol(Survmod))-1))
  
  #Export to clipboard and .csv
  if (csv_semicolon==TRUE) {write.csv2(Survmod,"PERSUADE_Time-to-event_models_parameters_semicolon.csv")}
  if (csv_comma==TRUE) {write.csv(Survmod,"PERSUADE_Time-to-event_models_parameters_comma.csv")}
  if (clipboard==TRUE) {write.table(Survmod, "clipboard-128", sep="\t")}
}
