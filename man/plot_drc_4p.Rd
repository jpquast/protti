% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_drc_4p.R
\name{plot_drc_4p}
\alias{plot_drc_4p}
\title{Plotting of four-parameter dose response curves}
\usage{
plot_drc_4p(
  data,
  grouping,
  response,
  dose,
  targets,
  unit = "uM",
  y_axis_name = "Response",
  facet = TRUE,
  scales = "free",
  x_axis_scale_log10 = TRUE,
  export = FALSE,
  export_name = "dose-response_curves"
)
}
\arguments{
\item{data}{A data frame that is obtained by calling the \code{fit_drc_4p} function.}

\item{grouping}{The name of the column containing precursor, peptide or protein identifiers.}

\item{response}{The name of the column containing response values, eg. log2 transformed intensities.}

\item{dose}{The name of the column containing dose values, eg. the treatment concentrations.}

\item{targets}{A character vector that specifies the names of the precursors, peptides or proteins (depending on \code{grouping})
that should be plotted. This can also be \code{"all"} if plots for all curve fits should be created.}

\item{unit}{A character vector specifying the unit of the concentration.}

\item{y_axis_name}{A character vector specifying the name of the y-axis of the plot.}

\item{facet}{A logical indicating if plots should be summarised into facets of 20 plots. This is recommended for many plots.}

\item{scales}{A character vector that specifies if the scales in faceted plots (if more than one target was provided) should be \code{"free"} or \code{"fixed"}.}

\item{x_axis_scale_log10}{A logical indicating if the x-axis scale should be log10 transformed.}

\item{export}{A logical indicating if plots should be exported as PDF. The output directory will be the current working directory. The
name of the file can be chosen using the \code{export_name} argument. If only one target is selected and \code{export = TRUE},
the plot is exported and in addition returned in R.}

\item{export_name}{A character vector providing the name of the exported file if \code{export = TRUE}.}
}
\value{
If \code{targets = "all"} a list containing plots for every unique identifier in the \code{grouping} variable is created. Otherwise a plot for the specified targets is created with maximally 20 facets.
}
\description{
Function for plotting four-parameter dose response curves for each group (precursor, peptide or protein), based on output from \code{fit_drc_4p} function.
}
\examples{
\dontrun{
plot_drc_4p(
  data,
  grouping = eg_precursor_id,
  response = intensity,
  dose = concentration,
  targets = c("ABCDEFK")
)
}
}
