## Check results

0 errors | 0 warnings | 1 note (unable to verify current time)

- This is a new release.
- This NOTE is caused by running on Windows where R CMD check sometimes cannot verify file timestamps. 
  It is unrelated to package functionality and can be safely ignored.

## Main changes compared with PERSUADE v0.1.1

- Added f_get_excel_template() to copy the Excel template `PERSUADE_Excel_template.xltx` to a user specified directory (default = `tempdir()`), this template provides a convenient structure for transferring survival model outputs into health economic models.
- The example workflow `PERSUADE_example_workflow.R` can be inspected by `file.edit(system.file("example_workflow", "PERSUADE_example_workflow.R", package = "PERSUADE"))`
- Minor updates to vignette and README files
- Added single quotes for software names and API (e.g. `rmarkdown`)
- Added references to DESCRIPTION file
- Replaced \dontrun with \donttest (and tested by `devtools::check(args = c("--run-donttest"))`), `except for f_generate_report()` as it requires LaTeX to be installed
- Set default output directories in `f_generate_report()` and `f_get_excel_template()` (and examples and unit tests) to `tempdir()` to avoid writing to the user’s working directory.
- Removed `options(scipen = 999, max.print = 10000, digits = 4)` from "R/examples/PERSUADE_example_workflow.R"
- Removed `rm(list = ls())` from "R/examples/PERSUADE_example_workflow.R"
- Removed `devtools::install_github("Bram-R/PERSUADE", quiet = TRUE, upgrade = "never")` from `R/examples/PERSUADE_example_workflow.R`
- Removed `options(warn = -1)` from `PERSUADE_output.Rmd`



