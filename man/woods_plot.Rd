% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/woods_plot.R
\name{woods_plot}
\alias{woods_plot}
\title{Woods' plot}
\usage{
woods_plot(
  data,
  fold_change,
  start_position,
  end_position,
  protein_length,
  coverage = NULL,
  protein_id,
  targets = "all",
  facet = TRUE,
  colouring = NULL,
  fold_change_cutoff = 1,
  highlight = NULL,
  export = FALSE,
  export_name = "woods_plots"
)
}
\arguments{
\item{data}{a data frame containing differential abundance, start and end peptide or precursor positions, protein
length and optionally a variable based on which peptides or precursors should be coloured.}

\item{fold_change}{a column in the data frame containing log2 fold changes.}

\item{start_position}{a column in the data frame containing the start positions for each peptide or precursor.}

\item{end_position}{a column in the data frame containing the end positions for each peptide or precursor.}

\item{protein_length}{a column in the data frame containing the length of the protein.}

\item{coverage}{optional, column in the data frame containing coverage in percent. Will appear in the title of the Woods' plot 
if provided.}

\item{protein_id}{a column in the data frame containing protein identifiers.}

\item{targets}{a character vector that specifies the identifiers of the proteins (depending on \code{protein_id})
that should be plotted. This can also be \code{"all"} if plots for all proteins should be created. Default is \code{"all"}.}

\item{facet}{a logical indicating if plots should be summarised into facets of 20 plots. This is recommended for many plots.
Default is \code{facet = TRUE}.}

\item{colouring}{optional, column in the data frame containing information by which peptide or precursors should
be coloured.}

\item{fold_change_cutoff}{optional, numeric argument specifying the log2 fold change cutoff used in the plot. The default value is 2.}

\item{highlight}{optional logical column containing logicals, specifying whether specific peptides or precursors should be highlighted with an asterisk.}

\item{export}{a logical indicating if plots should be exported as PDF. The output directory will be the current working directory. The
name of the file can be chosen using the \code{export_name} argument. Default is \code{export = FALSE}.}

\item{export_name}{a character vector providing the name of the exported file if \code{export = TRUE}. Default is \code{export_name = "woods_plots"}}
}
\value{
A Woods' plot is returned. Plotting peptide or precursor fold changes across protein sequence.
}
\description{
Creates a Woods' plot that plots log2 fold change of peptides or precursors along the protein sequence.
}
\examples{
\dontrun{
woods_plot(test,
  fold_change = diff,
  start_position = start,
  end_position = end,
  protein_length = length,
  colouring = pep_type,
  protein_id, pg_protein_accessions,
  facet = TRUE,
  highlight = is_significant
)
}
}
