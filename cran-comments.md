## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Changes compared with PERSUADE v0.1.1

- Added f_get_excel_template() to copy the Excel template `PERSUADE_Excel_template.xltx` to the current working directory, this template provides a convenient structure for transferring survival model outputs into health economic models.
- Added example workflow `PERSUADE_example_workflow.R` that can be inspected by `file.edit(system.file("example_workflow", "PERSUADE_example_workflow.R", package = "PERSUADE"))`
- Minor updates to vignette and README file
- Added single quotes for software names and API (e.g. 'rmarkdown')
- Added references to DESCRIPTION file
- Replaced \dontrun with \donttest



