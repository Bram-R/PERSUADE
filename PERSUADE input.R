rm(list = setdiff(ls(), lsf.str()))  #remove all objects except functions

####### LOAD AND INSTALL PACKAGES AND FUNCTIONS #######
packages <- c("rms", "survival", "flexsurv", "muhaz", "survminer", "ggplot2", "data.table", "summarytools", 
              "knitr", "kableExtra")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]  #check for new packages
if (length(new.packages)) install.packages(new.packages)  #install new packages (if needed)
suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE))  #load packages

source("PERSUADE function.R")

####### INPUT DATA #######
years <- bc$recyrs  #time in years
status <- bc$censrec  #status / event variable
group <- bc$group  #grouping variable, please amend the names of the groups here if you wish

years <- as.numeric(years)  #time variable should be numeric
status <- as.numeric(status)  #status / event variable should be numeric
group <- as.factor(group)  #group variable should be a factor

# INPUT TIME POINTS FOR PREDICTED SURVIVAL TABLE
time_pred_surv_table <- c(0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35)  #time points (in years) for predicted survival Table
time_unit <- 6/12  #time unit (in years) for predicted survival Figures
time_horizon <- 40  #time horizon (in years) for predicted survival Figures

####### PERSUADE #######
# PERSUADE FUNCTION (please adjust TRUE/ FALSE if necessary)
PERSUADE <- PERSUADE(years = years, status = status, group = group, strata = TRUE, time_unit = time_unit, 
                     time_horizon = time_horizon, time_pred_surv_table = time_pred_surv_table, spline_mod = TRUE, csv_semicolon = FALSE, 
                     csv_comma = TRUE, clipboard = FALSE)

# PERSUADE RMARKDOWN (PDF file)
evalq({
  rmarkdown::render(input = "PERSUADE main - bc case.Rmd", envir = new.env())  # please adjust name of RMD file if necessary
}, envir = PERSUADE)