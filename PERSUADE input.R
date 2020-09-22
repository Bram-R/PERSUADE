# code is formatted using formatR::tidy_source(width.cutoff = 100)
rm(list = setdiff(ls(), lsf.str()))  #remove all objects except functions

####### LOAD AND INSTALL PACKAGES AND FUNCTIONS #######
packages <- c("rms", "survival", "flexsurv", "muhaz", "survminer", "ggplot2", "data.table", "summarytools", 
              "knitr", "kableExtra", "sft")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]  #check for new packages
if (length(new.packages)) install.packages(new.packages)  #install new packages (if needed)
suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE))  #load packages

source("PERSUADE function.R")

####### INPUT DATA ####### (**please adjust to add your input data**)
name <- "BC_OS"  #name (will be printed on PDF as well as used to name output directory and PDF file)
# bc <- bc[bc$group==levels(bc$group)[3],] #used for validation purposes
# bc <- bc[bc$group==levels(bc$group)[-1],] #used for validation purposes
years <- bc$recyrs  #time in years
status <- bc$censrec  #status / event variable
group <- bc$group  #grouping variable, please amend the names of the groups here if you wish

years <- as.numeric(years)  #time variable should be numeric
status <- as.numeric(status)  #status / event variable should be numeric
group <- as.factor(group)  #group variable should be a factor (with a maximum of 3 levels)

# INPUT TIME POINTS FOR PREDICTED SURVIVAL TABLE
time_pred_surv_table <- c(0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35)  #time points (in years) for predicted survival Table
time_unit <- 1/12  #time unit (in years) for predicted survival Figures
time_horizon <- 40  #time horizon (in years) for predicted survival Figures

####### PERSUADE ####### 
# PERSUADE FUNCTION (**please adjust TRUE/ FALSE if necessary**)
PERSUADE <- f_PERSUADE(years = years, status = status, group = group, strata = TRUE, time_unit = time_unit, 
                           time_horizon = time_horizon, time_pred_surv_table = time_pred_surv_table, 
                           spline_mod = TRUE, csv_semicolon = FALSE, csv_comma = TRUE, clipboard = FALSE)

PERSUADE$name <- name
save(PERSUADE, file = "PERSUADE.RData")  # save PERSUADE (so it can be loaded in the RMARKDOWN script)

# PERSUADE RMARKDOWN (create PDF file and images) 
xfun::Rscript_call( #Rscript_call renders the Rmd document in a new R session (similar to clicking the Knit button in RStudio)
  rmarkdown::render,
  list(input = 'PERSUADE output.Rmd', output_file = paste0(name, '.pdf'), output_dir = paste0(name, '_output'), 
       intermediates_dir = paste0(name, '_output'))
)

# file management (move PERSUADE object to output directory)
file.copy(from = "PERSUADE.RData", to = paste0(name, "_output"))  #copy PERSUADE object to output directory
file.remove("PERSUADE.RData") #remove PERSUADE object from working directory



