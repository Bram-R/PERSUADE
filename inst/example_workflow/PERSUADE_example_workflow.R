#### Description: Example workflow for the ParamEtRic SUrvivAl moDel sElection (PERSUADE) communication tool ####
#### This is a standardised survival analysis tool to support the selection and communication of parametric survival models and their implementation in decision analytic models. ####

#### LOAD PERSUADE ---- (development version can be installed by: devtools::install_github("Bram-R/PERSUADE", quiet = TRUE, upgrade = "never")
library("PERSUADE")

#### INPUT DATA ----
# Example dataset
dat <- subset(PERSUADE::nsclc_ipd, dataset_id == "Immunotherapy_OAK_PDL1_unselected")

years  <- dat$time_months / 12  # time to event
status <- dat$event             # event indicator
group  <- dat$arm               # grouping variable

name <- "NSCLC_OS"        # project name
time_pred_surv_table <- c(0, 1, 2, 3, 4, 5, 10, 15, 20)  # survival table times
time_unit    <- 1/12   # monthly data in years
time_horizon <- 20     # analysis horizon (years)

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
print(PERSUADE)

summary(PERSUADE,type = "km")
summary(PERSUADE,type = "surv_probs")
summary(PERSUADE,type = "gof")
summary(PERSUADE,type = "gof_spline")
summary(PERSUADE,type = "gof_cure")

palette(rainbow(n = 9, s = 1, v = 1, start = 0, end = max(1, 9 - 1)/9, alpha = 1)) # Set colour palette for Figures
plot(PERSUADE, type = "km")
plot(PERSUADE, type = "ph")
plot(PERSUADE, type = "hr")
plot(PERSUADE, type = "param_models")
plot(PERSUADE, type = "spline_models")
plot(PERSUADE, type = "cure_models")
palette("default") # Set colour palette to default

# Create report
f_generate_report(PERSUADE, open = TRUE) # check RMD file: system.file("rmd", "PERSUADE_output.Rmd", package = "PERSUADE")

# Export parametric survival models to clipboard and CSV files (change `tempdir()` into `getwd()` for copying to working directory)
write.table(PERSUADE$surv_model_excel, "clipboard-128", sep = "\t", col.names = FALSE)
write.csv(PERSUADE$surv_model_excel, file.path(tempdir(), paste0(name, "_output"), "PERSUADE_Time-to-event_models_parameters_comma.csv"))
write.csv2(PERSUADE$surv_model_excel, file.path(tempdir(), paste0(name, "_output"), "PERSUADE_Time-to-event_models_parameters_semicolon.csv"))

# Obtain Excel template to incorporate model parameters in decision-analytic model
f_get_excel_template()
