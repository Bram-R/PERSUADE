rm(list = setdiff(ls(), lsf.str()))  #remove all objects except functions

####### LOAD AND INSTALL PACKAGES AND FUNCTIONS #######
packages <- c("rms", "survival", "flexsurv", "muhaz", "survminer", "ggplot2", "data.table", "summarytools", 
              "knitr", "kableExtra", "sft", "flexsurvcure")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]  #check for new packages
if (length(new.packages)) install.packages(new.packages)  #install new packages (if needed)
suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE))  #load packages

source("PERSUADE function.R")

####### INPUT DATA ####### (**please adjust to add your input data**)
name <- "BC_OS"  #name (will be printed on PDF as well as used to name output directory and PDF file)
#bc <- bc[bc$group==levels(bc$group)[1],] # used for validation purposes
#bc <- bc[bc$group==levels(bc$group)[-1],] # used for validation purposes
years <- bc$recyrs  #time in years
status <- bc$censrec  #status / event variable
group <- bc$group  #grouping variable (with a maximum of 3 levels), please amend the names of the groups here if you wish

# INPUT TIME POINTS FOR PREDICTED SURVIVAL TABLE
time_pred_surv_table <- c(0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35)  #time points (in years) for predicted survival Table
time_unit <- 1/12  #time unit (in years) for predicted survival Figures
time_horizon <- 40  #time horizon (in years) for predicted survival Figures

####### PERSUADE ####### 
# PERSUADE FUNCTION (**please adjust TRUE/ FALSE if necessary**)
PERSUADE <- f_PERSUADE(name = name, years = years, status = status, group = group, 
                       strata = TRUE, spline_mod = FALSE, cure_mod = TRUE, cure_link = "logistic", # link can either be "logistic", "loglog", "identity" or "probit"
                       time_unit = time_unit, time_horizon = time_horizon,
                       time_pred_surv_table = time_pred_surv_table)

# PERSUADE RMARKDOWN (create PDF file and images) 
dir.create(paste0(name, "_output"), showWarnings = FALSE) # create output directory
save(PERSUADE, file = paste0(name, "_output/PERSUADE.RData"))  # save PERSUADE in output directory (so it can be loaded in the RMD script)

xfun::Rscript_call( #Rscript_call renders the Rmd document in a new R session (similar to clicking the Knit button in RStudio)
  rmarkdown::render,
  list(input = "PERSUADE output.Rmd", output_file = paste0(name, ".pdf"), output_dir = paste0(name, "_output"), 
       intermediates_dir = paste0(name, "_output"))
)

# export parametric survival models to clipboard and .csv (e.g. for use in accompanying Excel file)
write.table(PERSUADE$surv_model$survmod, "clipboard-128", sep = "\t")
write.csv(PERSUADE$surv_model$survmod, paste0(name, "_output/PERSUADE_Time-to-event_models_parameters_comma.csv"))
write.csv2(PERSUADE$surv_model$survmod, paste0(name, "_output/PERSUADE_Time-to-event_models_parameters_semicolon.csv"))



