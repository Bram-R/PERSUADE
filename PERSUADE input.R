#### Description: The ParamEtRic SUrvivAl moDel sElection (PERSUADE) template ####
#### This is a standardised survival analysis template to support the selection of parametric survival models and their implementation in decision analytic models. ####

# General settings
options(scipen = 999, max.print = 10000, digits = 4)

# Load and install necessary libraries
required_packages <- c(
  "rms", "survival", "flexsurv", "muhaz", "survminer", "ggplot2", 
  "data.table", "summarytools", "knitr", "kableExtra", "sft", 
  "flexsurvcure", "docstring"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages) 
suppressPackageStartupMessages(lapply(required_packages, require, character.only = TRUE))

# Clear workspace
rm(list = ls())

#### LOAD PERSUADE FUNCTIONS ----
if (!file.exists("PERSUADE function.R")) stop("The file 'PERSUADE function.R' was not found.")
source("PERSUADE function.R")
source("PERSUADE figures.R")

# check PERSUADE functions using docstring()
# docstring(f_PERSUADE)
# docstring(f_cum_hazard)
# docstring(f_tp)
# docstring(f_surv_model)
# docstring(f_surv_model_pred)
# docstring(f_surv_model_pred_gr)
# docstring(f_surv_model_pred_tp_gr)
# docstring(f_surv_model_excel)

# check Figure functions using docstring()
# docstring(f_plot_km_survival)

#### INPUT DATA ----
name <- "BC_OS" # Analysis name

# Input variables
years <- bc$recyrs  # Time to event
status <- bc$censrec  # Event status
group <- bc$group  # Grouping variable (max 3 levels)

# Predicted survival table times (in years)
time_pred_surv_table <- c(0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35)
time_unit <- 1 / 12  # Time unit in years (monthly)
time_horizon <- 40  # Time horizon in years

#### RUN PERSUADE ----
PERSUADE <- f_PERSUADE(
  name = name, 
  years = years, 
  status = status, 
  group = group, 
  strata = TRUE, 
  spline_mod = TRUE, 
  cure_mod = TRUE, 
  cure_link = "logistic",  # Link options: "logistic", "loglog", "identity", "probit"
  time_unit = time_unit, 
  time_horizon = time_horizon, 
  time_pred_surv_table = time_pred_surv_table
)

#### EXPORT RESULTS ----
# Create output directory
output_dir <- paste0(name, "_output")
dir.create(output_dir, showWarnings = FALSE)

# Save PERSUADE object
save(PERSUADE, file = file.path(getwd(), "PERSUADE.RData"))

# Generate PDF report using R Markdown
if (!file.exists("PERSUADE output.Rmd")) stop("The file 'PERSUADE output.Rmd' was not found.")
xfun::Rscript_call(
  rmarkdown::render,
  list(
    input = "PERSUADE output.Rmd",
    output_file = paste0(name, ".pdf"),
    output_dir = output_dir,
    intermediates_dir = output_dir
  )
)

# Export parametric survival models to clipboard and CSV files
write.table(PERSUADE$surv_model_excel, "clipboard-128", sep = "\t", col.names = FALSE)
write.csv(PERSUADE$surv_model_excel, file.path(output_dir, "PERSUADE_Time-to-event_models_parameters_comma.csv"))
write.csv2(PERSUADE$surv_model_excel, file.path(output_dir, "PERSUADE_Time-to-event_models_parameters_semicolon.csv"))
