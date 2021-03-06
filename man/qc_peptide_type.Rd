% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_peptide_type.R
\name{qc_peptide_type}
\alias{qc_peptide_type}
\title{Check peptide type percentage share}
\usage{
qc_peptide_type(
  data,
  sample,
  peptide,
  pep_type,
  intensity,
  remove_na_intensities = TRUE,
  method = "count",
  plot = FALSE,
  interactive = FALSE
)
}
\arguments{
\item{data}{A data frame containing at least the input columns.}

\item{sample}{the name of the column containing the sample names.}

\item{peptide}{the name of the column containing the peptide sequence.}

\item{pep_type}{the name of the column containing the peptide type. Can be obtained using the \code{find_peptide} and \code{peptide_type} function together.}

\item{intensity}{the name of the column containing the corresponding raw or normalised intensity values (not log2) for each peptide or precursor. Required when "intensity" is chosen as the method.}

\item{remove_na_intensities}{Logical specifying if sample/peptide combinations with intensities that are NA (not quantified IDs) should
be dropped from the data frame for analysis of peptide type distributions. Default is TRUE since we are usually
interested in the peptide type distribution of quantifiable IDs. This is only relevant for method = "count".}

\item{method}{character vector indicating the method used for evaluation. \code{method = "intensity"} calculates the peptide type percentage by intensity, whereas \code{method = "count"} calculates the percentage by peptide ID count. Default is \code{method = count}.}

\item{plot}{a logical indicating whether the result should be plotted.}

\item{interactive}{a logical indicating whether the plot should be interactive.}
}
\value{
A data frame that contains the calculated percentage shares of each peptide type per sample. The \code{count} column contains the number of peptides with a specific type. The \code{peptide_type_percent} column contains the percentage share of a specific peptide type.
}
\description{
Calculates the percentage share of each peptide types (fully-tryptic, semi-tryptic, non-tryptic) for each sample.
}
\examples{
\dontrun{
qc_peptide_type(
  data,
  sample = r_file_name,
  peptide = pep_stripped_sequence,
  pep_type = pep_type,
  intensity = fg_quantity,
  method = "intensity",
  plot = TRUE
)
}
}
