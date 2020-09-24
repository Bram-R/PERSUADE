# PERSUADE

The ParamEtRic SUrvivAl moDel sElection (PERSUADE) template is a standardised survival analysis template to support the selection of parametric survival models and their implementation in decision analytic models. 

Specifically, The intention of the template was to provide information and support the decisions concerning parametric survival model selection:
1)	Visualisation of time-to-event data
2)	Examine the proportional hazard assumption
3)	Examine the hazard over time
4)	Visually examine the goodness of fit by a. visualisation of the parametric survival models (estimated survival function as well as transition probabilities) and b. examine diagnostics of the parametric survival models
5)	Examine statistical goodness of fit
6)	Examine extrapolation

How to use PERSUADE
1)  Make sure the the files "PERSUADE input.R", "PERSUADE function.R" and "PERSUADE output.Rmd" are saved in the R working directory
2)  Open "PERSUADE input.R" 
3)  Specify the input data and change the outcome/ project name using the "name" parameter
4)  Adjust TRUE/ FALSE in the PERSUADE function "f_PERSUADE(...)" if necessary
5)  Run all code in "PERSUADE input.R", the results will be stored in a subdirectory using the outcome/ project name 

