#### Description: Example workflow for the ParamEtRic SUrvivAl moDel sElection (PERSUADE) communication tool ####
#### This is a standardised survival analysis tool to support the selection and communication of parametric survival models and their implementation in decision analytic models. ####

# General settings
options(scipen = 999, max.print = 10000, digits = 4)

# Clear workspace
rm(list = ls())

#### LOAD PERSUADE ----
devtools::install_github("Bram-R/PERSUADE", quiet = TRUE, upgrade = "never") # To install the development version of PERSUADE
library("PERSUADE")

# Colour palette for Figures
n <- 9  #number of different colors (to be used for palette)
palette(rainbow(n = n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1))

#### INPUT DATA ----
name <- "BC_OS" # Analysis name
name <- file.path("inst", "example", name) # Analysis name including path

# bc1 <- bc[bc$group=="Medium",] # 1 group data set (for testing purposes)
# bc2 <- bc[bc$group!="Medium",] # 2 group data set (for testing purposes)

# Input variables
years <- flexsurv::bc$recyrs  # Time to event (alternative example: factor(survival::lung$time))
status <- flexsurv::bc$censrec  # Event status (alternative example: factor(survival::lung$status))
group <- flexsurv::bc$group  # Grouping variable (alternative example: factor(survival::lung$sex))

# Predicted survival table times (in years)
time_pred_surv_table <- c(0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35)
time_unit <- 1 / 12  # Time unit in years (monthly), e.g. use 365.25/12 when time unit is days
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

#### RESULTS ----
# S3 functionality
print(PERSUADE)

summary(PERSUADE,type = "km")
summary(PERSUADE,type = "surv_probs")
summary(PERSUADE,type = "gof")
summary(PERSUADE,type = "gof_spline")
summary(PERSUADE,type = "gof_cure")

plot(PERSUADE, type = "km")
plot(PERSUADE, type = "ph")
plot(PERSUADE, type = "hr")
plot(PERSUADE, type = "param_models")
plot(PERSUADE, type = "spline_models")
plot(PERSUADE, type = "cure_models")

# Create report
f_generate_report(PERSUADE) # check RMD file: system.file("rmd", "PERSUADE_output.Rmd", package = "PERSUADE")

# Export parametric survival models to clipboard and CSV files
write.table(PERSUADE$surv_model_excel, "clipboard-128", sep = "\t", col.names = FALSE)
write.csv(PERSUADE$surv_model_excel, file.path(file.path(getwd(), paste0(name, "_output")), "PERSUADE_Time-to-event_models_parameters_comma.csv"))
write.csv2(PERSUADE$surv_model_excel, file.path(file.path(getwd(), paste0(name, "_output")), "PERSUADE_Time-to-event_models_parameters_semicolon.csv"))
