## Submission comments

- This is the submission of PERSUADE 0.1.2. The package provides a standardized framework to support parametric survival model selection for decision-analytic models.
- All comments from the previous CRAN review have been addressed.
- Compared with version 0.1.1, this release introduces a new helper function (`f_get_excel_template()`), refines examples and vignettes to avoid writing to the user’s working directory, and includes minor documentation and workflow updates (details below).
- All local checks pass with 0 errors, 0 warnings, 0 notes (both `devtools::check()` and `devtools::check(args = "--as-cran")` have been run).
- The only NOTE observed on `devtools::check_win_devel()` (possibly misspelled words) is due to cited author names: “Grambsch”, “Therneau”, “Ishak”, “Royston”, “Parmar”.

## Check results using `devtools::check(args = "--as-cran")`

- 0 errors | 0 warnings | 0 notes

## Main changes compared with PERSUADE v0.1.1

- Added `f_get_excel_template()` to copy the Excel template `PERSUADE_Excel_template.xltx` to a user specified directory (default = `tempdir()`), this template provides a convenient structure for transferring survival model outputs into health economic models.
- The example workflow `PERSUADE_example_workflow.R` can be inspected by `file.edit(system.file("example_workflow", "PERSUADE_example_workflow.R", package = "PERSUADE"))`
- Minor updates to vignette and README files
- Added single quotes for software names and API (e.g. `rmarkdown`)
- Added references to DESCRIPTION file
- Replaced \dontrun with \donttest (and tested by `devtools::check(args = c("--run-donttest"))`), except for `f_generate_report()` as it requires LaTeX to be installed
- Set default output directories in `f_generate_report()` and `f_get_excel_template()` (and examples and unit tests) to `tempdir()` to avoid writing to the user’s working directory.
- Cleaned example workflow (`R/examples/PERSUADE_example_workflow.R`) by: 
   -- removing `options(scipen = 999, max.print = 10000, digits = 4)`
   -- removing `rm(list = ls())` 
   -- commenting `devtools::install_github("Bram-R/PERSUADE", quiet = TRUE, upgrade = "never")` 
- Removed `options(warn = -1)` from `PERSUADE_output.Rmd`

