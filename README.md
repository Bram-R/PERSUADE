# PERSUADE  
**ParamEtRic SUrvivAl moDel sElection**

A standardized R template designed to support the selection of parametric survival models and their implementation in decision-analytic models.

---

## Overview

PERSUADE is crafted to guide analysts through a robust, reproducible workflow for selecting the most appropriate **parametric survival model** for time-to-event data. It offers:

- Visual tools for exploring time-to-event data  
- Diagnostics for proportional hazards and hazard patterns over time  
- Goodness-of-fit evaluation via visual and statistical methods  
- Model extrapolation insights; essential for long-term projections  

This template streamlines complex survival analysis within decision modeling frameworks.

---

## Features & Workflow

The framework supports the analyst in:

1. **Visualizing** raw time-to-event data  
2. **Assessing** the proportional hazards assumption  
3. **Examining** hazard rates over time  
4. **Comparing** parametric survival model fits:
   - Visualization of survival curves and transition probabilities  
   - Model diagnostics  
5. **Quantifying** fit with statistical metrics  
6. **Exploring** model extrapolation beyond observed data  

---

## Repository Structure
├── PERSUADE_input.R # Main script to configure inputs and initiate the analysis
├── PERSUADE_function.R # Core function implementation: f_PERSUADE(...)
├── PERSUADE_figures.R # Visualization utilities
├── PERSUADE_S3_object_functions.R # S3 methods and object manipulation tools
├── PERSUADE_output.Rmd # Markdown template for generating reports
├── PERSUADE_figure_process.png/.vsdx # Visual diagram explaining the workflow
├── PERSUADE Excel template 2020_05_15.xltx # Optional spreadsheet template
├── README.md # This readme file
├── .gitignore
├── PERSUADE.Rproj # Optional RStudio project file
├── PERSUADE.enl # Auxiliary file
└── BC_OS_output # Output artifacts (generated results)

---

## Getting Started

### **Prerequisites**
- R (version compatible with tidyverse, survival, and other dependencies)
- Recommended: RStudio for better environment handling

### **Setup Instructions**

1. **Place files in your working directory:**  
   Ensure the following are present:
   - `PERSUADE_input.R`
   - `PERSUADE_function.R`  
   - `PERSUADE_figures.R`  
   - `PERSUADE_S3_object_functions.R`  
   - `PERSUADE_output.Rmd`  

2. **Open `PERSUADE_input.R`:**  
   - Specify your input dataset (e.g., CSV or R object)  
   - Set the project name and outcome identifier using the `name` parameter  

3. **Customize Analysis Options:**  
   Within the call to `f_PERSUADE(...)`, toggle options (`TRUE`/`FALSE`) for the analyses you want to run  

4. **Run the Workflow:**  
   Execute all code in `PERSUADE_input.R`. Outputs—including figures, results, and R Markdown reports—will be saved into a subdirectory named after the `name` parameter.

---

## Outputs & Interpretation

Upon execution, expect:

- A structured output directory tied to your project name  
- Visualizations: Kaplan–Meier curves, hazard plots, model-fit overlays  
- Diagnostic outputs: diagnostic plots, goodness-of-fit comparisons  
- An R Markdown report ready for rendering (PDF, HTML, etc.) detailing your analysis  

These outputs help critically evaluate model fit, especially for extrapolating survival beyond the observed timeline.

---

## Workflow Diagram

Visualize the process via the included diagram:  

![PERSUADE Workflow Overview](PERSUADE_figure_process.png)  

*(Alternatively, open the `.vsdx` version if preferred.)*

---

## Extending or Customizing

- **Add new model types:** Incorporate additional parametric distributions or non-parametric methods by extending `f_PERSUADE()`.  
- **Refine plots:** Customize visualizations in `PERSUADE_figures.R` or adopt themes for publication-quality graphics.  
- **Scale to new data formats:** Build data preprocessing modules for commonly used formats (e.g., survival objects, breakpoints, long-format datasets).  
- **Automate Reporting:** Enhance the R Markdown template for production-ready reports with headers, metadata, or version control.  

---

## Acknowledgments & History

- **Origins:** This repository represents an ongoing effort to standardize survival model workflows in decision analysis.  
- **Contributors:** Developed by **Bram-R** and **Xa4P**; contributions are welcome via pull requests or issues.  

---

## License & Contact

*(If applicable, include:)*  
Licensed under **[LICENSE NAME]**.  
For questions, enhancements, or collaboration, reach out via GitHub Discussions or Issues.  

---


