% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_charge_states.R
\name{qc_charge_states}
\alias{qc_charge_states}
\title{Check charge state distribution}
\usage{
qc_charge_states(
  data,
  sample,
  grouping,
  charge_states,
  intensity = NULL,
  remove_na_intensities = TRUE,
  method = "count",
  plot = FALSE,
  interactive = FALSE
)
}
\arguments{
\item{data}{A data frame containing at least sample names, peptide or precursor identifiers and missed cleavage counts for each peptide or precursor.}

\item{sample}{the column in the data data frame containing the sample name.}

\item{grouping}{the column in the data data frame containing either precursor or peptide identifiers.}

\item{charge_states}{the column in the data data frame containing the different charge states assigned to the precursor or peptide.}

\item{intensity}{the name of the column containing the corresponding raw or normalised intensity values (not log2) for each peptide or precursor. Required when "intensity" is chosen as the method.}

\item{remove_na_intensities}{logical specifying if sample/grouping combinations with intensities that are NA (not quantified IDs) should
be dropped from the data frame for analysis of missed cleavages. Default is TRUE since we are usually
interested in quantifiable peptides. This is only relevant for method = "count".}

\item{method}{character vector indicating the method used for evaluation. "count" calculates the charge state distribution based on counts of the corresponding peptides or precursors in the charge state group, "intensity" calculates the percentage of precursors or peptides in each charge state group based on the corresponding intensity values.}

\item{plot}{logical indicating whether the result should be plotted.}

\item{interactive}{argument specifying whether the plot should be interactive (default is FALSE).}
}
\value{
A data frame that contains the calculated percentage made up by the sum of either all counts or intensities of peptides or precursors of the corresponding charge state (depending on which method is chosen).
}
\description{
Calculates the charge state distribution for each sample (by count or intensity).
}
\examples{
\dontrun{
qc_charge_states(
  data,
  sample = r_file_name,
  grouping = pep_stripped_sequence,
  charge_states = fg_charge,
  method = "count",
  plot = TRUE
)
}
}
