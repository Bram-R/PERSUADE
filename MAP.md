# MAP — Development Log: PERSUADE Dashboard

This file documents the steps taken to add an interactive dashboard to the
PERSUADE R package, in the order they happened, including the constraints
of the development sandbox and how they shaped the approach.

## 1. Understand the existing package

PERSUADE is an R package (not a web app) that provides a standardized
workflow for parametric survival model selection:

- `f_PERSUADE()` — the main entry point. Fits Kaplan-Meier, Cox PH,
  7 standard parametric models, 9 spline models (`flexsurv`/`flexsurvcure`),
  and 12 cure models, returning an S3 object of class `"PERSUADE"`.
- `print.PERSUADE()`, `summary.PERSUADE()`, `plot.PERSUADE()` — S3 methods
  for inspecting results (`R/PERSUADE_S3_method.R`).
- ~25 `f_plot_*()` functions producing base-R diagnostic/extrapolation
  plots (`R/PERSUADE_output_functions.R`).
- `f_generate_report()` — renders an RMarkdown → PDF report.
- `f_get_excel_template()` — copies a bundled Excel template for exporting
  model parameters into decision-analytic models.

No data was bundled with the package (no `data/` directory); the README,
vignette, and example workflow script all use `flexsurv::bc` (breast
cancer recurrence data) as example data. No Shiny/web code existed
anywhere in the repository.

## 2. Clarify ambiguous instructions

Two items in the task were ambiguous and worth resolving before writing
code rather than guessing:

1. "Create an AGENTS.md file with this content" — no content was actually
   attached to the message.
2. The dashboard's implementation technology was unspecified.

Asked the user directly (`AskUserQuestion`); both times the recommended
option was chosen:

- AGENTS.md: generate the standard `agents.md`-style template, tailored to
  this repository.
- Dashboard: an R Shiny app living inside the package (`inst/shiny-dashboard/`),
  launched via a new exported `f_run_dashboard()` function — the idiomatic
  approach for an R package, keeping everything in one language/repo.

## 3. Sandbox constraints discovered (important — read this before running tests)

This development session runs in a sandboxed environment with **no network
route to CRAN** (`cran.r-project.org`, `cloud.r-project.org`, and
`packagemanager.posit.co` all reject the CONNECT tunnel). Ubuntu's `apt`
mirrors *are* reachable, so `r-base` plus everything packaged as
`r-cran-*` in the Ubuntu repositories was installed successfully:
`survival`, `rms`, `survminer`, `data.table`, `ggplot2`, `rmarkdown`,
`testthat`, `shiny`, `shinydashboard`, `DT`, `bslib`, `shinycssloaders`,
`kableExtra`, `knitr`, `roxygen2`.

However, **`flexsurv`, `flexsurvcure`, `muhaz`, and `sft`** — all hard
dependencies in `DESCRIPTION`'s `Imports:` — are not packaged for Ubuntu
and could not be installed. Consequences, discovered empirically:

- `R CMD INSTALL .` fails outright ("dependencies ... are not available").
- `pkgload::load_all(".")` fails the same way.
- `roxygen2::roxygenise()` fails the same way (it calls `pkgload::load_all()`
  internally).
- The existing test harness (`tests/testthat.R`) does `library(PERSUADE)`,
  which therefore cannot run in this sandbox either — this is a pre-existing
  property of the repository, not something introduced here.

Practical effect: **nothing that calls `f_PERSUADE()`, `f_hazard()`,
`f_cum_hazard()`, or `f_surv_model()` (all of which need `flexsurv`/
`muhaz`/`sft`) could be executed or verified locally.** The project's own
CI (`.github/workflows/rhub.yaml`, via r-hub) has full CRAN access and is
the real verification path for those code paths — this was already true
of the existing test suite before this change.

## 4. Red/green TDD, adapted to the constraint

Work was scoped so the genuinely testable parts (new code, no `flexsurv`/
`muhaz`/`sft` dependency) got real red→green cycles locally, and the parts
that need the full model-fitting stack were written carefully, statically
checked, and left for CI/the user's machine to execute.

**a. Bundled example dataset (`persuade_trial_data`)**

- Wrote `tests/testthat/test_persuade_trial_data.R` first (structural
  checks always runnable; an `f_PERSUADE()` round-trip test guarded by
  `skip_if_not_installed("flexsurv")` etc. for CI).
- Confirmed red via a scratch harness that loads the `.rda` directly
  (`data("...", package = "PERSUADE")` itself needs the package installed,
  which isn't possible here — see §3 — so the harness substitutes a direct
  `load()` for local verification only; the committed test is the one that
  actually runs in CI).
- Implemented `data-raw/create_persuade_trial_data.R` (simulates a 3-arm
  trial with `stats::rweibull`/`stats::rexp` — no `flexsurv` needed to
  *simulate* data, only to *fit* models to it) and `R/PERSUADE_data.R`.
- Reran → green.

**b. Dashboard helper functions (`R/PERSUADE_dashboard.R`)**

Extracted the dashboard's business logic into small, pure(ish) functions
specifically so they could be unit tested without the heavy model-fitting
stack: `f_dashboard_validate_data()`, `f_dashboard_time_grid()`,
`f_dashboard_example_datasets()`, `f_dashboard_best_model()`, and
`f_dashboard_run_persuade()` (a `tryCatch` wrapper around `f_PERSUADE()`
so the dashboard can show a friendly error instead of crashing).

- Wrote `tests/testthat/test_dashboard_helpers.R` first (17 `test_that()`
  blocks). The committed version uses `testthat::local_mocked_bindings()`
  to stub `f_PERSUADE()` for the wrapper's success/failure paths — the
  idiomatic approach, but it requires the package namespace to be loaded,
  which (§3) isn't possible here.
- For local verification only, mirrored the same assertions in a scratch
  harness that sources the file directly and stubs `f_PERSUADE` by
  reassigning it in the sourced function's environment instead of via
  `local_mocked_bindings()`.
- Confirmed red (file didn't exist) → implemented → green (all 41
  assertions passing), including a real bug caught along the way: the
  first draft's error message for a missing column didn't name which
  logical role (years/status/group) was missing; the test asserted on
  that, caught it, and the implementation was fixed to report roles.

**c. The Shiny app itself (`inst/shiny-dashboard/app.R`)**

Not unit-testable in the traditional sense (it's UI wiring), so verified
by actually running it:

1. Parsed for syntax errors, then loaded with `source()` (not `Rscript
   app.R` directly, which would auto-print and launch the server) to
   confirm the `ui`/`server` objects construct without error.
2. Launched it for real (`shiny::runApp(..., port = 8791)`) and drove it
   with a headless Chromium via Playwright (pre-installed in this
   environment): loaded the page, clicked "Load / validate data", clicked
   "Run PERSUADE analysis", navigated every tab, and captured screenshots.
3. This caught a real bug: the "Data source" `radioButtons()` used the
   value `"example"`, but the dataset registry (`f_dashboard_example_datasets()`)
   keys its bundled entry as `"persuade_trial_data"`. The mismatch made
   `f_dashboard_validate_data()` receive `NULL` column names, producing a
   cryptic `attempt to select less than one element in get1index` error —
   only visible by actually clicking the button in a browser, not from
   reading the code. Fixed by aligning the radio value with the registry
   key.
4. Because `f_PERSUADE()` itself needs `flexsurv` (§3), the one thing that
   *cannot* be verified end-to-end here is an actual model fit completing
   successfully inside the running app. What *was* verified: data loads,
   validates, and previews correctly; stat tiles and the auto-suggested
   time horizon/grid populate correctly; clicking "Run" calls
   `f_dashboard_run_persuade()`, which — even when the underlying call
   fails (here, because `f_PERSUADE` isn't loadable at all in this
   sandbox) — is caught and surfaced as a non-blocking notification
   instead of crashing the app, which is exactly the intended behavior.

## 5. Docs and package metadata

`roxygen2::roxygenise()` couldn't run (§3), so `NAMESPACE` and the new
`man/*.Rd` files were written by hand, mirroring the exact format/style of
the existing (roxygen2-generated) files, and validated with
`tools::parse_Rd()`. `DESCRIPTION` gained `shiny`, `shinydashboard`, `DT`,
`bslib`, `shinycssloaders` under `Suggests` (not `Imports`, so the base
package stays lightweight and only needs these when someone actually runs
the dashboard — `f_run_dashboard()` checks for them via `requireNamespace()`
and gives an actionable error if missing, following the same pattern
already used by `f_generate_report()` for its own Suggested packages).
`README.Rmd`/`README.md` and `NEWS.md` were updated to document
`f_run_dashboard()`.

## 6. Post-review fix: `library(PERSUADE)` missing from `app.R`

After the initial push, the user ran `f_run_dashboard()` for real (on a
machine with the full dependency stack) and hit:

```
Warning: Error in f_dashboard_example_datasets: could not find function "f_dashboard_example_datasets"
```

Root cause: `shiny::runApp()` sources `app.R` into a **new environment
whose parent is `globalenv()`**, not the PERSUADE package namespace. `app.R`
attached `shiny`/`shinydashboard`/`DT`/`shinycssloaders` via `library()` but
never attached PERSUADE itself, so every unqualified call to a PERSUADE
function (`f_PERSUADE()`, all `f_dashboard_*()`, all `f_plot_*()`,
`f_generate_report()`, `f_get_excel_template()`) was unresolvable.

This slipped through local verification (§4c) because the sandbox
verification harness called `source("R/PERSUADE_dashboard.R")` *before*
`shiny::runApp(...)` in the same R session — which defines those functions
directly in `.GlobalEnv`, and since `app.R`'s sourcing environment's parent
*is* `.GlobalEnv`, they resolved by accident. That path only exists in the
sandbox harness; a real user calling `f_run_dashboard()` never sources
`PERSUADE_dashboard.R` directly, so the functions were never actually
reachable. Lesson: a "does it run" check that pre-populates the global
environment with the code under test isn't a substitute for reproducing the
real call path.

Fix: added `library(PERSUADE)` to the top of `inst/shiny-dashboard/app.R`,
alongside the UI-library `library()` calls. Verified the mechanism (not
just patched blindly) with a minimal reproduction: sourcing a stand-in
"app.R" that calls an unqualified function defined only in a detached
environment reproduces the exact `could not find function` error; attaching
that environment to the search path first (what `library()` does) resolves
it. See `tests/` — no new automated test was added for this specific
regression because it is purely about *environment attachment*, which
requires the real installed package (unavailable in this sandbox, see §3)
to exercise meaningfully; a `shiny::runApp()`-based integration check
against an actually-installed PERSUADE is the right place to catch this
class of bug going forward, ideally in CI once dependencies are available.

## 7. What to verify next (on a machine with CRAN access)

- `R CMD check` end-to-end (impossible here per §3).
- Run the full `testthat` suite via `library(PERSUADE); testthat::test_check("PERSUADE")`
  and confirm the `skip_if_not_installed()`-guarded tests actually execute
  (they were only exercised via the scratch-harness substitutes here).
- Click through the dashboard with `spline_mod`/`cure_mod` enabled against
  real data and confirm a fit completes (the 12 cure models in particular
  are the slowest part of `f_PERSUADE()` — worth timing on a realistic
  dataset size to judge whether the dashboard needs an async/progress
  affordance beyond the current `withProgress()` spinner).
- Try `f_generate_report()`'s PDF download from the dashboard on a machine
  with pandoc + a LaTeX distribution installed (untestable here).
