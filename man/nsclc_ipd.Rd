\name{nsclc_ipd}
\alias{nsclc_ipd}
\docType{data}

\title{
Reconstructed Individual Patient Data for Overall Survival in NSCLC Trials
}

\description{
Reconstructed pseudo-individual patient data (pseudo-IPD) for overall
survival (OS) from published Kaplan-Meier curves in trials of targeted
therapies and immunotherapies for advanced non-small cell lung cancer
(NSCLC). Reconstruction used digitized survival coordinates and, where
available, published numbers at risk and total event counts, following
the approach described by Guyot et al. (2012).

A review of all 33 trial-population datasets found that the reconstructed
curves, survival medians, and treatment-effect estimates generally agree
closely with their published sources. Known discrepancies and limitations
are documented below; the records are not the original patient data.
}

\usage{nsclc_ipd}

\format{
A data frame with 11,572 reconstructed patient records from 33
trial-population datasets and 66 treatment arms, with 8 variables:
\describe{
  \item{\code{dataset_id}}{
    Character. Identifier for a trial-population comparison. One trial can
    contribute more than one population or subgroup.
  }
  \item{\code{trial_id}}{
    Character. Name or identifier of the source clinical trial.
  }
  \item{\code{population}}{
    Character. Population or subgroup represented by the published OS curves.
  }
  \item{\code{arm}}{
    Factor with levels \code{control} and \code{intervention}. Harmonized
    treatment assignment used in this data frame. These labels identify
    the comparison, not necessarily the reference category of a published HR.
  }
  \item{\code{source_arm}}{
    Factor with levels \code{control} and \code{intervention}. Arm label
    used in the digitized source archive. It can differ from \code{arm}
    where source-arm labels have been reassigned.
  }
  \item{\code{patient_id}}{
    Numeric. Artificial identifier within each dataset and arm. These
    identifiers do not identify real patients, including across overlapping
    subgroups reconstructed from the same trial.
  }
  \item{\code{time_months}}{
    Numeric. Reconstructed time from randomization to death or right-censoring,
    expressed in months.
  }
  \item{\code{event}}{
    Numeric. Event indicator: \code{1} for death and \code{0} for
    right-censoring.
  }
}
}

\details{
\strong{Origin and purpose.}

These data build on the published survival curves and digitized source data
from Simons et al. (2020). Kaplan-Meier coordinates were converted into
pseudo-IPD using an inversion/reconstruction procedure based on Guyot et al.
(2012). Where reported, numbers at risk and total deaths informed the
reconstruction. The resulting observations are estimates consistent with
published aggregate survival information; they are \emph{not} the original
individual patient records.

The collection contains comparisons involving EGFR-targeted treatment,
ALK-targeted treatment, and immunotherapy. Each \code{dataset_id} identifies
one two-arm survival comparison. Population eligibility, line of treatment,
treatment sequence, analysis cutoff, and survival maturity differ across
trials. Some datasets are biomarker-defined subgroups rather than the whole
randomized trial. In some trials, different PD-L1 or mutation-defined
datasets overlap. Such datasets must not be analysed or pooled as independent
patient cohorts.

\strong{Overall validation conclusion.}

All 33 datasets were reviewed against the corresponding published OS
Kaplan-Meier figures and reported survival outcomes. Across the collection,
reconstructed curve shapes, median OS point estimates, and Cox hazard ratios
generally correspond very closely to the published results and to the original
recreated estimates in Simons et al. (2020). Several comparisons reproduce
published arm-specific death totals exactly. For other comparisons,
reconstructed death totals differ by a few events, and some original
publications do not report arm-specific event totals. The most pronounced
event-count discrepancy concerns CTONG-0901.

The validation is a comparison with \emph{published aggregate information},
not patient-level validation against unavailable original IPD. Visual
agreement, agreement in the HR, or a matching median does not establish
correct individual event or censoring histories. The trial-specific notes
below disclose important differences without implying that they must be
resolved before the datasets can be used for methodological research.

\strong{PICO and reconstructed event proportions.}

The trial-specific section reports the population, intervention, comparator,
and outcome (PICO) for each comparison, together with the reconstructed
sample sizes and event proportions for both arms. Because OS is the common
outcome for all datasets, it is consistently defined as death from any
cause. Event proportions describe the reconstructed observations rather than
validated published event totals; they should not be treated as direct
measures of follow-up maturity or survival-tail precision. They can be
calculated directly from \code{nsclc_ipd}, as shown in the examples.

\strong{Interpreting reported quantities.}

Unless stated otherwise, reconstructed HRs compare \code{intervention}
with \code{control} in an unstratified Cox model, and reconstructed median
OS values come from the reconstructed Kaplan-Meier estimator. Published
HRs can be stratified, adjusted, or expressed in the opposite direction;
published p-values can additionally reflect prespecified interim analyses
and multiplicity adjustment. Therefore, differences in HRs, confidence
intervals, and p-values need not be reconstruction errors. A confidence
interval for median OS may be unreached even when its point estimate is
available; estimates near the 50 percent KM threshold are especially
sensitive to small differences in event and censoring allocation.

\strong{Scope of the data.}

The pseudo-IPD contain survival times, censoring indicators, trial-population
identifiers, and treatment assignments, but not individual mutation status,
PD-L1 expression, other prognostic factors, subsequent therapy, crossover,
or treatment-switching history. Where crossover or subsequent therapies
occurred, the published OS curves generally describe survival under the
randomized treatment strategies, not a counterfactual estimate adjusted
for switching. Reconstructed subgroup records cannot be matched to the
same individuals in a separately reconstructed whole-trial dataset.

\strong{Finding original KM curves.}

The companion data documentation \code{\link{nsclc_trials}} (see
\file{nsclc_trials.Rd}) describes the trial-level metadata, including
\code{dataset_id}, treatment labels, publication details, DOI, and the
\code{trial_url} column. A convenient way to open the original publication
containing a survival curve is:
\preformatted{browseURL(nsclc_trials$trial_url[28])}
With the supplied ordering, row 28 is the PD-L1 >=50 percent KEYNOTE-042
comparison. The URL resolves to a publication page; the KM figure may need
to be opened within the article or its PDF, and publisher access conditions
may apply. For reproducible code, selecting a row by \code{dataset_id}
is safer than relying on its numerical row position.
}

\section{Trial-specific validation findings}{
The descriptions below summarize a manual comparison of the reconstructed
pseudo-IPD with the relevant published OS figure and statistics. Values
headed as published or reported refer to the particular source analysis,
not necessarily a later follow-up publication. Death-count discrepancies
are stated where there is an identifiable published or independently
sourced reference; unknown or unspecified counts are not treated as zero.
Each entry provides a compact PICO description (population, intervention,
comparator, outcome), the reconstructed sample size, number of deaths, and
event proportion for each harmonized treatment arm, followed by the original
validation findings. Arm-level figures are calculated from the current
\code{nsclc_ipd} records; they are not necessarily the deaths reported in
the corresponding publication. \emph{Intervention} and \emph{control} refer
to the harmonized \code{arm} labels in this package, not necessarily the
reference contrast in the original study. Treatment and population details
refer to the particular reconstructed comparison, which may be a subgroup
or two selected arms from a larger randomized trial.

The event proportion is the reconstructed number of deaths divided by the
reconstructed number of patients in that arm (multiplied by 100). It is a
simple descriptor of observed event information, sometimes informally used
as a proxy for survival-data maturity, but \emph{not} a direct measure of
follow-up completeness, precision of tail estimates, or the percentage of
patients who will eventually experience the event. The event proportion
may differ between arms because of a treatment effect, different censoring
patterns, or reconstruction discrepancies. In particular, the high
reconstructed event proportions in CTONG-0901 should not be interpreted as
verified published OS maturity. Comparisons with very small sample sizes,
such as TITAN, remain highly uncertain irrespective of the percentage.

Abbreviations: HR, hazard ratio; NR, not reached; ITT, intention to treat;
TPS, tumour proportion score; WT, wild type.

\describe{
  \item{\code{ALK_ALEX}}{
    \strong{Population (P):} Previously untreated advanced ALK-positive NSCLC
    in the ALEX trial.

    \strong{Intervention (I):} Alectinib.

    \strong{Comparator (C):} Crizotinib.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=152,
    deaths=36 (23.7 percent); control: N=151, deaths=40 (26.5 percent).

    \strong{Validation and limitations:}
    Alectinib versus crizotinib in ALK-positive NSCLC. Both OS medians are NR.
    The reconstructed HR is approximately 0.81 versus a published HR of
    approximately 0.76. The reconstruction has 36 alectinib deaths versus
    35 specified in the source reconstruction information, while the 40
    crizotinib deaths agree. OS was relatively immature.

  }
  \item{\code{ALK_ALTA-L1}}{
    \strong{Population (P):} Previously untreated advanced ALK-positive NSCLC
    in ALTA-1L.

    \strong{Intervention (I):} Brigatinib.

    \strong{Comparator (C):} Crizotinib.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=137,
    deaths=17 (12.4 percent); control: N=138, deaths=17 (12.3 percent).

    \strong{Validation and limitations:}
    Brigatinib versus crizotinib, ALK-positive population. Both arms have
    17 reconstructed deaths, matching reported counts, and neither median
    is reached. Reconstructed HR 1.01 is close to the published 0.98;
    only 34 deaths in 275 patients inform this early OS comparison.

  }
  \item{\code{ALK_ASCEND-4}}{
    \strong{Population (P):} Previously untreated advanced ALK-positive NSCLC
    in ASCEND-4.

    \strong{Intervention (I):} Ceritinib.

    \strong{Comparator (C):} Platinum-based chemotherapy plus pemetrexed.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=189,
    deaths=47 (24.9 percent); control: N=187, deaths=59 (31.6 percent).

    \strong{Validation and limitations:}
    Ceritinib versus platinum/pemetrexed in ALK-positive NSCLC. The
    reconstructed chemotherapy median is 26.33 months versus 26.2 published;
    ceritinib median is NR. The reconstructed HR is 0.71 versus 0.73
    published. There is one fewer reconstructed ceritinib death (47 versus
    48); the chemotherapy event count (59) matches.

  }
  \item{\code{ALK_PROFILE1007}}{
    \strong{Population (P):} Previously treated advanced ALK-positive NSCLC in
    PROFILE 1007.

    \strong{Intervention (I):} Crizotinib.

    \strong{Comparator (C):} Pemetrexed or docetaxel chemotherapy.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=173,
    deaths=48 (27.7 percent); control: N=174, deaths=47 (27.0 percent).

    \strong{Validation and limitations:}
    Crizotinib versus pemetrexed or docetaxel in previously treated
    ALK-positive NSCLC. The reconstructed HR is approximately 1.02,
    matching the published 1.02; medians are approximately 20.5 and 23.0
    months versus 20.3 and 22.8 published. The reconstruction has 48
    crizotinib deaths versus 49 specified, with 47 chemotherapy deaths
    reproduced. Subsequent crizotinib treatment in the chemotherapy arm
    complicates interpretation of the ITT OS comparison.

  }
  \item{\code{ALK_PROFILE1014}}{
    \strong{Population (P):} Previously untreated advanced ALK-positive NSCLC
    in PROFILE 1014; mature OS analysis.

    \strong{Intervention (I):} Crizotinib.

    \strong{Comparator (C):} Platinum-based chemotherapy plus pemetrexed.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=172,
    deaths=72 (41.9 percent); control: N=171, deaths=81 (47.4 percent).

    \strong{Validation and limitations:}
    Crizotinib versus platinum/pemetrexed in ALK-positive NSCLC. The
    digitized curves correspond to the more mature 2018 OS report, not
    the earlier 2014 initial report. Reconstructed chemotherapy median
    43.91 months differs from 47.5 months published, largely in a KM
    region near the 50 percent threshold. Crizotinib median is NR;
    reconstructed HR 0.82 agrees closely with the Simons-recreated HR.
    There is one additional reconstructed crizotinib death (72 versus 71).

  }
  \item{\code{EGFR_CTONG-0901}}{
    \strong{Population (P):} Advanced EGFR-mutated NSCLC in the CTONG-0901
    comparison.

    \strong{Intervention (I):} Erlotinib.

    \strong{Comparator (C):} Gefitinib.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=128,
    deaths=117 (91.4 percent); control: N=128, deaths=121 (94.5 percent).

    \strong{Validation and limitations:}
    Erlotinib versus gefitinib in EGFR-mutated NSCLC. The original digitized
    source-arm labels were reversed relative to these treatment names;
    \code{arm} gives the harmonized assignment. Median OS and the
    reconstructed HR (approximately 0.86) match the publication closely.
    \emph{However}, the reconstruction contains 121 gefitinib and 117
    erlotinib deaths versus 92 deaths reported for each arm. This major
    discrepancy indicates that the event/censoring histories, despite their
    similar KM shapes and HR, should not be treated as faithfully recovered.

  }
  \item{\code{EGFR_DELTA}}{
    \strong{Population (P):} Previously treated advanced NSCLC with an
    activating EGFR mutation; small DELTA subgroup.

    \strong{Intervention (I):} Erlotinib.

    \strong{Comparator (C):} Docetaxel.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=21, deaths=4
    (19.0 percent); control: N=30, deaths=12 (40.0 percent).

    \strong{Validation and limitations:}
    EGFR-mutated subgroup of DELTA: erlotinib (21 patients) versus
    docetaxel (30). The reconstructed HR 0.45 resembles the published 0.43;
    the docetaxel median is approximately 27.8 months and the erlotinib
    median is NR. Only four and 12 deaths, respectively, are reconstructed.
    Original arm-specific total-event constraints were recorded as NA;
    the small sample and sparse tail limit precision.

  }
  \item{\code{EGFR_ENSURE}}{
    \strong{Population (P):} Previously untreated advanced EGFR-mutated NSCLC
    in ENSURE.

    \strong{Intervention (I):} Erlotinib.

    \strong{Comparator (C):} Cisplatin plus gemcitabine.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=110,
    deaths=57 (51.8 percent); control: N=107, deaths=56 (52.3 percent).

    \strong{Validation and limitations:}
    Erlotinib versus cisplatin/gemcitabine. Reconstructed medians are
    26.42 and 25.66 months, closely matching published medians 26.3 and
    25.5 months; reconstructed HR 0.94 versus 0.91 published. The
    chemotherapy reconstruction contains 56 rather than 54 specified
    deaths; the 57 erlotinib deaths agree.

  }
  \item{\code{EGFR_EURTAC}}{
    \strong{Population (P):} Previously untreated advanced EGFR
    mutation-positive NSCLC in EURTAC.

    \strong{Intervention (I):} Erlotinib.

    \strong{Comparator (C):} Platinum-based chemotherapy.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=86, deaths=38
    (44.2 percent); control: N=87, deaths=31 (35.6 percent).

    \strong{Validation and limitations:}
    Erlotinib versus platinum-based chemotherapy in EGFR-mutated NSCLC.
    Reconstructed median OS (19.37 versus 19.56 months) and HR (1.05)
    are close to the published estimates (19.3 versus 19.5 months and
    HR 1.04). The specified death totals, 38 and 31, are reproduced;
    confidence limits for median survival are not identical.

  }
  \item{\code{EGFR_FLAURA}}{
    \strong{Population (P):} Previously untreated advanced EGFR-mutated NSCLC
    (exon 19 deletion or L858R); initial FLAURA OS analysis.

    \strong{Intervention (I):} Osimertinib.

    \strong{Comparator (C):} Gefitinib or erlotinib.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=279,
    deaths=58 (20.8 percent); control: N=277, deaths=83 (30.0 percent).

    \strong{Validation and limitations:}
    First-line osimertinib versus gefitinib or erlotinib. Both median
    OS estimates are NR at the early analysis cutoff. Published and
    reconstructed event counts match (58 osimertinib, 83 control), and
    reconstructed HR 0.62 closely matches the published 0.63. These
    pseudo-IPD do not represent the later mature OS analysis.

  }
  \item{\code{EGFR_IMPRESS}}{
    \strong{Population (P):} Advanced EGFR mutation-positive NSCLC following
    progression on first-line gefitinib.

    \strong{Intervention (I):} Continued gefitinib plus cisplatin and
    pemetrexed.

    \strong{Comparator (C):} Placebo plus cisplatin and pemetrexed.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=133,
    deaths=94 (70.7 percent); control: N=132, deaths=80 (60.6 percent).

    \strong{Validation and limitations:}
    Gefitinib plus cisplatin/pemetrexed versus placebo plus the same
    chemotherapy following progression on gefitinib. Reconstructed medians
    are 13.40 and 19.85 months versus published 13.4 and 19.5; reconstructed
    HR 1.46 agrees with 1.44 published. The chemotherapy-only arm has
    80 versus 82 specified deaths; 94 deaths in the gefitinib-containing
    arm are reproduced.

  }
  \item{\code{EGFR_IPASS}}{
    \strong{Population (P):} EGFR mutation-positive subgroup of the first-line
    IPASS trial, not its full population.

    \strong{Intervention (I):} Gefitinib.

    \strong{Comparator (C):} Carboplatin plus paclitaxel.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=132,
    deaths=102 (77.3 percent); control: N=129, deaths=95 (73.6 percent).

    \strong{Validation and limitations:}
    EGFR mutation-positive subgroup, not the full IPASS randomized
    population. Gefitinib is compared with carboplatin/paclitaxel.
    Reconstructed medians are 21.02 and 22.19 months, close to published
    21.6 and 21.9 months; reconstructed HR 1.05 versus published 1.00.
    The gefitinib reconstruction contains 102 versus 104 published deaths;
    all 95 chemotherapy deaths are reproduced.

  }
  \item{\code{EGFR_Lux-Lung3}}{
    \strong{Population (P):} Previously untreated advanced EGFR-mutated NSCLC
    in the overall LUX-Lung 3 population.

    \strong{Intervention (I):} Afatinib.

    \strong{Comparator (C):} Cisplatin plus pemetrexed.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=230,
    deaths=139 (60.4 percent); control: N=115, deaths=72 (62.6 percent).

    \strong{Validation and limitations:}
    First-line afatinib versus cisplatin/pemetrexed in the overall EGFR-
    mutated trial population. Reconstructed medians are 29.44 and 28.31
    months versus 28.2 months in both published arms. The reconstructed
    HR 0.883 closely matches the published 0.88 and the Simons reconstruction.
    One death is missing in each arm (139 versus 140 afatinib, 72 versus
    73 chemotherapy). Treatment effects differ by EGFR mutation subtype;
    individual subtype information is unavailable here.

  }
  \item{\code{EGFR_Lux-Lung6}}{
    \strong{Population (P):} Previously untreated advanced EGFR-mutated NSCLC
    in LUX-Lung 6.

    \strong{Intervention (I):} Afatinib.

    \strong{Comparator (C):} Gemcitabine plus cisplatin.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=242,
    deaths=161 (66.5 percent); control: N=122, deaths=83 (68.0 percent).

    \strong{Validation and limitations:}
    First-line afatinib versus gemcitabine/cisplatin. Reconstructed medians
    are 23.21 and 23.80 months versus 23.1 and 23.5 published. Reconstructed
    HR 0.89 matches Simons' recreated 0.89, compared with a published HR
    of 0.93. One death is missing in each arm (161 versus 162 afatinib,
    83 versus 84 chemotherapy). Published stratified and reconstructed
    unstratified p-values differ.

  }
  \item{\code{EGFR_Lux-Lung7}}{
    \strong{Population (P):} Previously untreated advanced NSCLC with EGFR exon
    19 deletion or L858R mutation.

    \strong{Intervention (I):} Afatinib.

    \strong{Comparator (C):} Gefitinib.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=160,
    deaths=109 (68.1 percent); control: N=159, deaths=118 (74.2 percent).

    \strong{Validation and limitations:}
    First-line afatinib versus gefitinib among patients with exon 19
    deletion or L858R EGFR mutations. Reconstructed medians (27.39 and
    24.34 months) and HR (0.852) closely match published values (27.9,
    24.5, and 0.86). The reconstruction has one additional gefitinib
    death (118 versus 117); all 109 afatinib deaths are reproduced.
    Late KM steps are based on small risk sets.

  }
  \item{\code{EGFR_NCT02959749}}{
    \strong{Population (P):} Advanced EGFR T790M-positive NSCLC previously
    treated with an EGFR-TKI and platinum-based chemotherapy.

    \strong{Intervention (I):} Osimertinib.

    \strong{Comparator (C):} Docetaxel plus bevacizumab.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=74, deaths=14
    (18.9 percent); control: N=73, deaths=15 (20.5 percent).

    \strong{Validation and limitations:}
    Osimertinib versus docetaxel/bevacizumab after prior EGFR-directed and
    platinum-based therapy in EGFR T790M-positive NSCLC. Both KM medians
    are NR and reconstructed HR 0.80 agrees with published 0.79.
    The original report contains internally inconsistent median OS values:
    its text and curves imply NR, whereas a table reports numerical
    medians with incompatible confidence intervals. The reconstruction
    has 14 and 15 deaths; the disposition diagram uses a different
    cutoff and does not unambiguously verify final OS event totals.

  }
  \item{\code{EGFR_OPTIMAL}}{
    \strong{Population (P):} Previously untreated advanced EGFR-mutated NSCLC
    in OPTIMAL.

    \strong{Intervention (I):} Erlotinib.

    \strong{Comparator (C):} Gemcitabine plus carboplatin.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=82, deaths=68
    (82.9 percent); control: N=72, deaths=53 (73.6 percent).

    \strong{Validation and limitations:}
    Erlotinib versus gemcitabine/carboplatin. Reconstructed medians
    (22.86 versus 27.22 months) and HR (1.20) closely reproduce published
    values (22.8 versus 27.2 months; HR 1.19). The reconstruction contains
    one extra chemotherapy death (53 versus 52), while all 68 erlotinib
    deaths are reproduced. The terminal numbers-at-risk interval was
    omitted during reconstruction because it was not adequately supported
    by available digitized KM coordinates; retain this provenance caveat.

  }
  \item{\code{EGFR_TITAN}}{
    \strong{Population (P):} Activating EGFR mutation-positive subgroup of
    previously treated advanced NSCLC in TITAN (11 patients).

    \strong{Intervention (I):} Erlotinib.

    \strong{Comparator (C):} Docetaxel or pemetrexed.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=7, deaths=3
    (42.9 percent); control: N=4, deaths=1 (25.0 percent).

    \strong{Validation and limitations:}
    EGFR mutation-positive subgroup of TITAN, with only seven erlotinib
    and four chemotherapy patients. The reconstruction contains three and
    one deaths, respectively; published subgroup death totals are not
    explicitly reported. The published chemotherapy KM plateau is about
    67 percent versus 75 percent reconstructed. Reconstructed HR 1.39
    versus 1.19 published has an extremely wide confidence interval.
    This subgroup is unsuitable for reliable comparative-effect estimation
    and particularly sensitive to individual reconstructed observations.

  }
  \item{\code{EGFR_TORCH}}{
    \strong{Population (P):} EGFR mutation-positive subgroup of TORCH;
    randomized comparison of treatment sequences.

    \strong{Intervention (I):} Erlotinib first, followed by cisplatin plus
    gemcitabine at progression.

    \strong{Comparator (C):} Cisplatin plus gemcitabine first, followed by
    erlotinib at progression.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=19, deaths=14
    (73.7 percent); control: N=20, deaths=10 (50.0 percent).

    \strong{Validation and limitations:}
    EGFR-mutated subgroup comparing \emph{treatment sequences}:
    erlotinib followed by cisplatin/gemcitabine versus the reverse
    sequence. All published deaths are reproduced (14 versus 10), and
    reconstructed HR 1.51 closely matches Simons' 1.51 and published
    subgroup HR 1.58. The chemotherapy-first median is 26.0 versus
    32.5 months published, although both curves lie very near the 50
    percent KM threshold; the discrepancy remains documented rather
    than corrected. Subgroup numbers at risk are not provided in the figure.

  }
  \item{\code{EGFR_WJOG-5108L}}{
    \strong{Population (P):} Previously treated EGFR mutation-positive lung
    adenocarcinoma in WJOG-5108L, including uncommon mutations.

    \strong{Intervention (I):} Erlotinib.

    \strong{Comparator (C):} Gefitinib.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=198,
    deaths=96 (48.5 percent); control: N=203, deaths=107 (52.7 percent).

    \strong{Validation and limitations:}
    EGFR mutation-positive subgroup comparing erlotinib
    (\code{intervention}) with gefitinib (\code{control}). The published
    HR is expressed in the opposite direction (gefitinib versus erlotinib):
    published 1.189 versus inverted reconstructed approximately 1.181.
    Reconstructed medians are 30.36 months for erlotinib and 26.60 for
    gefitinib, versus published 31.4 and 26.5. Subgroup-specific published
    death totals are unavailable; counts for the larger trial population
    are not interchangeable.

  }
  \item{\code{Immunotherapy_CheckMate-026_PDL1_5}}{
    \strong{Population (P):} Previously untreated advanced NSCLC with PD-L1
    expression of at least 5 percent.

    \strong{Intervention (I):} Nivolumab.

    \strong{Comparator (C):} Platinum-doublet chemotherapy.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=211,
    deaths=130 (61.6 percent); control: N=212, deaths=122 (57.5 percent).

    \strong{Validation and limitations:}
    First-line nivolumab versus platinum-doublet chemotherapy among patients
    with PD-L1 expression >=5 percent. Reconstructed medians are 14.57
    and 13.40 months versus 14.4 and 13.2 published; reconstructed HR
    1.04 versus published 1.02. Reconstructed deaths are 130 and 122;
    the figure alone does not establish arm-specific published OS totals.
    The small late risk sets affect the KM tails.

  }
  \item{\code{Immunotherapy_CheckMate_017_PDL1_unselected}}{
    \strong{Population (P):} Previously treated advanced squamous NSCLC,
    without a PD-L1 expression eligibility threshold.

    \strong{Intervention (I):} Nivolumab.

    \strong{Comparator (C):} Docetaxel.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=135,
    deaths=86 (63.7 percent); control: N=137, deaths=113 (82.5 percent).

    \strong{Validation and limitations:}
    Previously treated \emph{squamous} NSCLC, nivolumab versus docetaxel;
    PD-L1 expression was not an entry requirement. All published deaths
    are reproduced (86 nivolumab; 113 docetaxel). Reconstructed medians
    are 9.27 and 6.04 months versus 9.2 and 6.0 published; HR 0.591
    versus 0.59 published. The reconstructed upper confidence limit for
    nivolumab median OS is 15.39 versus 13.3 months published.

  }
  \item{\code{Immunotherapy_CheckMate_057_PDL1_unselected}}{
    \strong{Population (P):} Previously treated advanced nonsquamous NSCLC,
    without a PD-L1 expression eligibility threshold.

    \strong{Intervention (I):} Nivolumab.

    \strong{Comparator (C):} Docetaxel.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=292,
    deaths=192 (65.8 percent); control: N=290, deaths=221 (76.2 percent).

    \strong{Validation and limitations:}
    Previously treated \emph{nonsquamous} NSCLC, nivolumab versus
    docetaxel; PD-L1 unselected. The total 413 published deaths are
    reproduced, but assigned differently: 192 versus 190 published
    nivolumab deaths and 221 versus 223 docetaxel deaths. Reconstructed
    medians (12.16 and 9.45 months) and HR 0.75 agree closely with
    published values. The published primary HR has a 96 percent confidence
    interval, whereas the reconstruction reports 95 percent; curves show
    early crossing and a later nivolumab advantage.

  }
  \item{\code{Immunotherapy_IMpower150_PDL1_unselected}}{
    \strong{Population (P):} Previously untreated metastatic nonsquamous NSCLC,
    EGFR/ALK wild-type ITT population, irrespective of PD-L1; two of three
    randomized arms.

    \strong{Intervention (I):} Atezolizumab plus bevacizumab, carboplatin, and
    paclitaxel (ABCP).

    \strong{Comparator (C):} Bevacizumab plus carboplatin and paclitaxel (BCP).

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=359,
    deaths=179 (49.9 percent); control: N=337, deaths=196 (58.2 percent).

    \strong{Validation and limitations:}
    Atezolizumab/bevacizumab/carboplatin/paclitaxel (ABCP) versus
    bevacizumab/carboplatin/paclitaxel (BCP). This dataset includes only
    these two arms of the \emph{EGFR/ALK wild-type} ITT population, not
    the third randomized arm or all trial patients. Reconstructed medians
    (19.46 versus 14.72 months) and HR 0.784 closely match published
    estimates. ABCP deaths match (179); BCP has 196 versus 197 published.
    Data are from an interim OS analysis.

  }
  \item{\code{Immunotherapy_KEYNOTE-010_PDL1_1}}{
    \strong{Population (P):} Previously treated advanced NSCLC with PD-L1 TPS
    of at least 1 percent; includes the at least 50 percent subgroup.

    \strong{Intervention (I):} Pembrolizumab 2 mg/kg every three weeks; the 10
    mg/kg arm is excluded.

    \strong{Comparator (C):} Docetaxel.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=344,
    deaths=174 (50.6 percent); control: N=343, deaths=197 (57.4 percent).

    \strong{Validation and limitations:}
    PD-L1 TPS >=1 percent, pembrolizumab \emph{2 mg/kg} versus docetaxel.
    The pembrolizumab 10 mg/kg randomized arm is not included. Medians
    (10.42 and 8.57 months) almost exactly match published values;
    reconstructed HR 0.729 versus published 0.71. There are two extra
    pembrolizumab deaths (174 versus 172) and four extra docetaxel
    deaths (197 versus 193). This dataset overlaps the >=50 percent
    KEYNOTE-010 subgroup below.

  }
  \item{\code{Immunotherapy_KEYNOTE-010_PDL1_50}}{
    \strong{Population (P):} Previously treated advanced NSCLC with PD-L1 TPS
    of at least 50 percent; nested in the at least 1 percent population.

    \strong{Intervention (I):} Pembrolizumab 2 mg/kg every three weeks; the 10
    mg/kg arm is excluded.

    \strong{Comparator (C):} Docetaxel.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=139,
    deaths=59 (42.4 percent); control: N=152, deaths=88 (57.9 percent).

    \strong{Validation and limitations:}
    PD-L1 TPS >=50 percent subset of KEYNOTE-010, again comparing
    pembrolizumab \emph{2 mg/kg} with docetaxel. Published and reconstructed
    medians agree almost exactly (14.9 and 8.2 months); reconstructed
    HR 0.600 versus published 0.54. There is one extra pembrolizumab
    death (59 versus 58) and two extra docetaxel deaths (88 versus 86).
    This subset is not independent of the >=1 percent dataset.

  }
  \item{\code{Immunotherapy_KEYNOTE-024_PDL1_50}}{
    \strong{Population (P):} Previously untreated advanced NSCLC with PD-L1 TPS
    of at least 50 percent, without sensitizing EGFR/ALK alterations.

    \strong{Intervention (I):} Pembrolizumab 200 mg every three weeks.

    \strong{Comparator (C):} Platinum-doublet chemotherapy.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=154,
    deaths=44 (28.6 percent); control: N=151, deaths=64 (42.4 percent).

    \strong{Validation and limitations:}
    First-line pembrolizumab versus platinum-doublet chemotherapy,
    PD-L1 TPS >=50 percent without sensitizing EGFR/ALK alterations.
    Both median OS values were \emph{not reached} in the original 2016
    publication and in the current reconstruction. All 44 pembrolizumab
    and 64 chemotherapy deaths are reproduced; reconstructed HR 0.615
    versus published 0.60. Simons et al. (2020) instead report a 12.6-month
    chemotherapy median in both comparison columns, inconsistent with the
    original analysis. The chemotherapy KM curve is near 50 percent,
    making the median classification sensitive to reconstruction details.

  }
  \item{\code{Immunotherapy_KEYNOTE-042_PDL1_50}}{
    \strong{Population (P):} Previously untreated advanced NSCLC with PD-L1 TPS
    of at least 50 percent, without sensitizing EGFR/ALK alterations.

    \strong{Intervention (I):} Pembrolizumab 200 mg every three weeks.

    \strong{Comparator (C):} Platinum-doublet chemotherapy.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=299,
    deaths=159 (53.2 percent); control: N=300, deaths=201 (67.0 percent).

    \strong{Validation and limitations:}
    First-line pembrolizumab versus chemotherapy, PD-L1 TPS >=50 percent.
    Reconstructed medians are 20.02 and 12.18 months versus 20.0 and
    12.2 published; reconstructed HR 0.703 versus published 0.69.
    Published deaths inferred from the final cumulative censoring table
    are 157 pembrolizumab and 199 chemotherapy; reconstructed totals are
    159 and 201, respectively. Published curves cross in early follow-up.
    The corresponding event references must not be reversed between arms.

  }
  \item{\code{Immunotherapy_KEYNOTE-189_PDL1_50}}{
    \strong{Population (P):} Previously untreated metastatic nonsquamous NSCLC
    without sensitizing EGFR/ALK alterations, PD-L1 TPS at least 50 percent.

    \strong{Intervention (I):} Pembrolizumab plus platinum/pemetrexed, followed
    by pembrolizumab/pemetrexed maintenance.

    \strong{Comparator (C):} Placebo plus platinum/pemetrexed, followed by
    placebo/pemetrexed maintenance.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=132,
    deaths=33 (25.0 percent); control: N=70, deaths=35 (50.0 percent).

    \strong{Validation and limitations:}
    PD-L1 TPS >=50 percent subgroup, pembrolizumab plus platinum/pemetrexed
    versus placebo plus the same chemotherapy (132 versus 70 patients).
    The published Figure 2C reports an HR of 0.42 versus reconstructed
    0.408; the survival curves agree closely. Reconstructed deaths total
    68 versus 70 reported for the subgroup. An ancillary source reports
    34 intervention and 36 control deaths, versus 33 and 35 reconstructed;
    these arm-specific numbers are not in the main KM figure. The numerical
    subgroup medians are not reported in that figure. The subgroup overlaps
    the whole-trial dataset below.

  }
  \item{\code{Immunotherapy_KEYNOTE-189_PDL1_unselected}}{
    \strong{Population (P):} Previously untreated metastatic nonsquamous NSCLC
    without sensitizing EGFR/ALK alterations, irrespective of PD-L1.

    \strong{Intervention (I):} Pembrolizumab plus platinum/pemetrexed, followed
    by pembrolizumab/pemetrexed maintenance.

    \strong{Comparator (C):} Placebo plus platinum/pemetrexed, followed by
    placebo/pemetrexed maintenance.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=410,
    deaths=125 (30.5 percent); control: N=206, deaths=107 (51.9 percent).

    \strong{Validation and limitations:}
    Full ITT population with metastatic \emph{nonsquamous} NSCLC without
    sensitizing EGFR/ALK alterations, receiving pembrolizumab plus
    platinum/pemetrexed or placebo plus the same regimen, including
    maintenance treatment. Reconstructed median OS is NR versus 11.25
    months (published NR versus 11.3), and reconstructed HR 0.492 matches
    published 0.49. Reconstructed events are 125 and 107, compared with
    127 and 108 in regulatory efficacy information. The reconstructed
    chemotherapy median upper confidence limit is NR versus published
    15.1 months. The PD-L1 >=50 percent subgroup is nested in this cohort.

  }
  \item{\code{Immunotherapy_KEYNOTE-407_PDL1_unselected}}{
    \strong{Population (P):} Previously untreated metastatic squamous NSCLC,
    irrespective of PD-L1 expression.

    \strong{Intervention (I):} Pembrolizumab with carboplatin and
    paclitaxel/nab-paclitaxel, followed by pembrolizumab maintenance.

    \strong{Comparator (C):} Placebo with carboplatin and
    paclitaxel/nab-paclitaxel, followed by placebo maintenance.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=278,
    deaths=87 (31.3 percent); control: N=281, deaths=114 (40.6 percent).

    \strong{Validation and limitations:}
    Full first-line \emph{squamous} NSCLC population, pembrolizumab plus
    carboplatin and paclitaxel/nab-paclitaxel versus placebo plus the same
    chemotherapy, with subsequent study-treatment maintenance. Reconstructed
    medians (15.89 versus 11.44 months) agree with published 15.9 versus
    11.3; reconstructed HR 0.677 is close to Simons' 0.68 and published
    0.64. Death counts differ: 87 versus 85 intervention and 114 versus
    120 control. The source is an immature interim OS analysis.

  }
  \item{\code{Immunotherapy_OAK_PDL1_unselected}}{
    \strong{Population (P):} Previously treated advanced squamous or
    nonsquamous NSCLC, irrespective of PD-L1; first 850 randomized patients in
    the primary ITT analysis.

    \strong{Intervention (I):} Atezolizumab.

    \strong{Comparator (C):} Docetaxel.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=425,
    deaths=271 (63.8 percent); control: N=425, deaths=298 (70.1 percent).

    \strong{Validation and limitations:}
    Atezolizumab versus docetaxel in previously treated, PD-L1-unselected
    squamous or nonsquamous NSCLC. The reconstruction covers the first
    \emph{850 randomized patients} used in the primary efficacy analysis,
    not all 1,225 ultimately randomized. Published event totals agree
    exactly (271 atezolizumab, 298 docetaxel). Medians (13.68 and 9.63
    months) and reconstructed HR 0.731 closely match published values.
    The trial-flow diagram omits one subsequently ascertained death in
    each arm; the reconstruction correctly reflects the efficacy totals.

  }
  \item{\code{Immunotherapy_POPLAR_PDL1_unselected}}{
    \strong{Population (P):} Previously treated advanced squamous or
    nonsquamous NSCLC, irrespective of PD-L1; full randomized ITT population.

    \strong{Intervention (I):} Atezolizumab.

    \strong{Comparator (C):} Docetaxel.

    \strong{Outcome (O):} Overall survival (death from any cause).

    \strong{Reconstructed arm-level information:} Intervention: N=144,
    deaths=77 (53.5 percent); control: N=143, deaths=95 (66.4 percent).

    \strong{Validation and limitations:}
    Atezolizumab versus docetaxel in previously treated, PD-L1-unselected
    squamous or nonsquamous NSCLC. Reconstructed medians (12.63 versus
    9.72 months) and HR 0.724 are very close to published 12.6 versus
    9.7 months and HR 0.73. The reconstruction contains 77 versus 78
    published atezolizumab deaths and correctly reproduces all 95 docetaxel
    deaths. The trial-flow diagram lists 93 docetaxel deaths, but its
    footnote adds two deaths identified from public records. The upper
    reconstructed atezolizumab median confidence limit is NR versus
    published 16.4 months.

  }
}
}

\section{Appropriate and inappropriate uses}{
These data are intended for education, reproducible methodological
illustration, Kaplan-Meier estimation, proportional hazards diagnostics,
parametric and flexible survival modelling, and survival extrapolation.
They may also be used to study the sensitivity of survival analyses to
imperfect reconstruction, incomplete follow-up, small risk sets, and
alternative censoring allocations.

They should not be presented as verified original patient-level trial
data. Patient-specific covariates or causal effects of subsequent therapy
cannot be inferred from reconstructed records. Comparisons from nested or
overlapping subgroups should not be treated as statistically independent;
differences in analysis cutoffs should be respected. The particularly
large event-count mismatch in CTONG-0901 and the extremely small TITAN
subgroup warrant caution in analyses that depend on accurate event or
censoring histories. More broadly, late-tail predictions should account
for reconstruction uncertainty as well as sampling and extrapolation
uncertainty.

Survival times are expressed in months. For applications using years,
including \code{f_PERSUADE()}, divide \code{time_months} by 12.
}

\source{
The digitized Kaplan-Meier coordinates and associated numbers-at-risk
source files originate from:

Simons, M., et al. (2020). Digitized patient level time to event data of
overall survival. Zenodo.
\href{https://doi.org/10.5281/zenodo.4288605}{doi:10.5281/zenodo.4288605}

The original digitized files were released under the Creative Commons
Attribution 4.0 International (CC BY 4.0) licence. The accompanying
reconstructed pseudo-IPD in this package were assembled from these
sources. The companion \code{\link{nsclc_trials}} dataset supplies the
individual source publications and their URLs.
}

\references{
Simons, M., Ramaekers, B., Peeters, A., Mankor, J., Paats, M., Aerts, J.,
van Harten, W., Retel, V., and Joore, M. (2020). Observed versus modelled
lifetime overall survival of targeted therapies and immunotherapies for
advanced non-small cell lung cancer patients - A systematic review.
Critical Reviews in Oncology/Hematology, 153, 103035.
\href{https://doi.org/10.1016/j.critrevonc.2020.103035}{doi:10.1016/j.critrevonc.2020.103035}

Guyot, P., Ades, A. E., Ouwens, M. J. N. M., and Welton, N. J. (2012).
Enhanced secondary analysis of survival data: reconstructing the data
from published Kaplan-Meier survival curves.
BMC Medical Research Methodology, 12, 9.
\href{https://doi.org/10.1186/1471-2288-12-9}{doi:10.1186/1471-2288-12-9}

Trial-specific primary publications are identified by DOI and URL in
\code{\link{nsclc_trials}}.
}

\seealso{
\code{\link{nsclc_trials}}
}

\examples{
data("nsclc_ipd")
data("nsclc_trials")

# Inspect reconstructed patient histories and available comparisons
head(nsclc_ipd)
unique(nsclc_ipd$dataset_id)

# Arm-level reconstructed sample sizes, deaths, and event percentages
arm_summary <- aggregate(
  data.frame(
    N = rep(1L, nrow(nsclc_ipd)),
    deaths = nsclc_ipd$event
  ),
  by = list(
    dataset_id = nsclc_ipd$dataset_id,
    arm = as.character(nsclc_ipd$arm)
  ),
  FUN = sum
)
arm_summary$event_pct <- round(
  100 * arm_summary$deaths / arm_summary$N, 1
)

# Compare arms in a specific trial population
subset(
  arm_summary,
  dataset_id == "Immunotherapy_KEYNOTE-042_PDL1_50"
)

# Optionally attach the population and treatment descriptions
ix <- match(arm_summary$dataset_id, nsclc_trials$dataset_id)
arm_summary$population <- nsclc_trials$population[ix]
arm_summary$control_treatment <- nsclc_trials$control_treatment[ix]
arm_summary$intervention_treatment <- nsclc_trials$intervention_treatment[ix]

# Select LUX-Lung 3 and inspect deaths and censoring by treatment
luxlung3 <- subset(nsclc_ipd, dataset_id == "EGFR_Lux-Lung3")
with(luxlung3, table(arm, event))

# Recreate Kaplan-Meier estimates from the pseudo-IPD
km <- survival::survfit(
  survival::Surv(time_months, event) ~ arm,
  data = luxlung3
)
print(km)

# Trial metadata, citation, and original publication URL
subset(nsclc_trials, dataset_id == "EGFR_Lux-Lung3")

# Open the publication containing the original KM curves
# Row 28 is KEYNOTE-042 (PD-L1 >=50 percent) in the supplied order.
if (interactive()) {
  browseURL(nsclc_trials$trial_url[28])
}

# For reproducibility, locate the DOI URL using dataset_id instead
idx <- match("EGFR_Lux-Lung3", nsclc_trials$dataset_id)
if (interactive() && !is.na(idx)) {
  browseURL(nsclc_trials$trial_url[idx])
}
}

\keyword{datasets}
