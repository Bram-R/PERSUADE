
\name{nsclc_trials}
\alias{nsclc_trials}
\docType{data}

\title{
Trial Characteristics and References for Reconstructed NSCLC Survival Data
}

\description{
Trial-level metadata accompanying the reconstructed
pseudo-individual patient data (pseudo-IPD) in
\code{nsclc_ipd}. The dataset provides information
on the original clinical trials, patient populations,
treatment comparisons, publication references,
and data reconstruction provenance.
}

\usage{nsclc_trials}

\format{
A data frame containing 33 trial-population datasets
with the following 13 variables:

\describe{
  \item{\code{dataset_id}}{
    Character. Unique identifier linking each
    trial-population dataset to the corresponding
    reconstructed observations in \code{nsclc_ipd}.
  }
  \item{\code{trial_id}}{
    Character. Name or identifier of the original
    clinical trial.
  }
  \item{\code{population}}{
    Character. Patient population or subgroup
    represented by the reconstructed survival data.
    Different populations from the same clinical
    trial may have separate dataset identifiers.
  }
  \item{\code{source_group}}{
    Character. Broad clinical or treatment category
    used to organize the source data, including
    EGFR-targeted therapies, ALK-targeted therapies,
    and immunotherapies.
  }
  \item{\code{control_treatment}}{
    Character. Treatment or treatment combination
    assigned to the control arm.
  }
  \item{\code{intervention_treatment}}{
    Character. Treatment or treatment combination
    assigned to the intervention arm.
  }
  \item{\code{included_in_paper_analysis}}{
    Logical. Indicates whether the dataset was included
    in the quantitative survival analysis reported by
    Simons et al. (2020). Additional datasets available
    in the source archive are also represented.
  }
  \item{\code{trial_publication}}{
    Character. Short reference identifying the original
    clinical trial publication or publications.
  }
  \item{\code{trial_doi}}{
    Character. Digital Object Identifier (DOI) of the
    original clinical trial publication, where available.
    Missing or empty entries do not necessarily indicate
    that a publication has no DOI.
  }
  \item{\code{trial_bibliography}}{
    Character. Bibliographic information for the
    original clinical trial publication or publications
    associated with the reconstructed survival data.
  }
  \item{\code{source_notes}}{
    Character. Additional information concerning
    the source data, population definitions, treatment
    assignments, or data reconstruction.
  }
  \item{\code{swap_source_arms}}{
    Numeric. Indicator specifying whether the original
    source-file arm labels were reversed when assembling
    the reconstructed dataset:
    0 = original labels retained;
    1 = source labels reversed.
  }
  \item{\code{trial_url}}{
  Character. Direct URL to the original clinical trial
  publication, constructed from \code{trial_doi}.
  The URL resolves through https://doi.org/.
}
}
}

\details{
The metadata accompany reconstructed overall survival
data from clinical trials investigating targeted
therapies and immunotherapies in patients with
locally advanced or metastatic non-small cell
lung cancer (NSCLC).

The original source data were collected and digitized
as part of the systematic review by Simons et al.
(2020). Published Kaplan-Meier curves and numbers
at risk were used to reconstruct pseudo-individual
patient time-to-event data using the methodology
described by Guyot et al. (2012).

Each row represents a unique trial-population dataset
containing a control and intervention arm. Some
clinical trials contribute multiple datasets
representing distinct or potentially overlapping
patient subgroups.

The variable \code{dataset_id} provides the link
between this metadata and the reconstructed
observations in \code{nsclc_ipd}.

The treatment variables describe the control and
intervention assignments in the assembled dataset.
The variable \code{swap_source_arms} documents
whether these assignments required reversing the
original source-file labels. The corresponding
\code{arm} and \code{source_arm} variables are
available in \code{nsclc_ipd}.

The variable \code{included_in_paper_analysis}
distinguishes datasets included in the quantitative
analysis by Simons et al. (2020) from additional
datasets available in the digitized source archive.
This indicator should not be interpreted as an
assessment of data quality or reconstruction accuracy.

Bibliographic information is provided to enable
users to identify and cite the original clinical
trial publications alongside the Simons dataset
and reconstruction methodology.

The reconstructed observations are not original
patient-level trial data. They are intended primarily
for educational and methodological research,
including survival analysis, model selection,
and survival extrapolation.
}

\source{
The reconstructed pseudo-individual patient data were
generated from digitized Kaplan-Meier survival curves
and numbers-at-risk tables made publicly available by
Simons et al. (2020):

Simons, M., et al. (2020).
Digitized patient level time to event data of overall survival.
Zenodo.
\href{https://doi.org/10.5281/zenodo.4288605}{
doi:10.5281/zenodo.4288605}

The original digitized data are distributed under the
Creative Commons Attribution 4.0 International
(CC BY 4.0) license.
}

\references{
Simons, M., Ramaekers, B., Peeters, A., Mankor, J.,
Paats, M., Aerts, J., van Harten, W., Retel, V.,
and Joore, M. (2020).
Observed versus modelled lifetime overall survival
of targeted therapies and immunotherapies for
advanced non-small cell lung cancer patients -
A systematic review.
Critical Reviews in Oncology/Hematology, 153, 103035.
\href{https://doi.org/10.1016/j.critrevonc.2020.103035}{doi:10.1016/j.critrevonc.2020.103035}

Guyot, P., Ades, A. E., Ouwens, M. J. N. M.,
and Welton, N. J. (2012).
Enhanced secondary analysis of survival data:
reconstructing the data from published
Kaplan-Meier survival curves.
BMC Medical Research Methodology, 12, 9.
\href{https://doi.org/10.1186/1471-2288-12-9}{doi:10.1186/1471-2288-12-9}

References to the original clinical trial
publications are included in the dataset.
}

\seealso{
\code{\link{nsclc_ipd}}
}

\examples{
data("nsclc_trials")

# Inspect trial metadata
head(nsclc_trials)

# List available clinical trials
unique(nsclc_trials$trial_id)

# Select LUX-Lung 3 metadata
luxlung3_info <- subset(
  nsclc_trials,
  dataset_id == "EGFR_Lux-Lung3"
)

# View treatment comparison and references
luxlung3_info[, c(
  "trial_id",
  "control_treatment",
  "intervention_treatment",
  "trial_bibliography"
)]

# Identify datasets included in the published analysis
published_trials <- subset(
  nsclc_trials,
  included_in_paper_analysis
)

nrow(published_trials)
}

\keyword{datasets}
