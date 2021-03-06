% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcano_protti.R
\name{volcano_protti}
\alias{volcano_protti}
\title{Volcano plot}
\usage{
volcano_protti(
  data,
  grouping,
  log2FC,
  significance,
  method,
  target_column = NULL,
  target = NULL,
  facet_by = NULL,
  title = "Volcano plot",
  x_axis_label = "log2(fold change)",
  y_axis_label = "-log10(q-value)",
  legend_label = "Target",
  log2FC_cutoff = 1,
  significance_cutoff = 0.01,
  interactive = FALSE
)
}
\arguments{
\item{data}{a data frame containing at least the input variables.}

\item{grouping}{the column in the data data frame containing either precursor or peptide identifiers.}

\item{log2FC}{the column in the data frame containing the log2 transfromed fold changes between two conditions.}

\item{significance}{the column containing the p-value or adjusted p-value for the corresponding fold changes. P-value is ideally adjusted using e.g. Benjamini-Hochberg correction.}

\item{method}{character verctor with the method used for the plot. \code{method = "target"} highlights your protein, proteins or any other entities of
interest (specified in the `target` argument) in the volcano plot. \code{method = "significant"} highlights all significantly changing entities.}

\item{target_column}{optional column required for \code{method = "target"}, can contain for example protein identifiers or a logical that marks
certain proteins such as proteins that are known to interact with the treatment. Can also be provided if \code{method = "significant"}
to label data points in an interactive plot.}

\item{target}{optional character vector argument required for \code{method = "target"}. It can contain one or more specific entities of the
column provided in \code{target_column}. This can be for example a protein ID if \code{target_column} contains protein IDs or TRUE or FALSE for a logical column.}

\item{facet_by}{optional argument specifying a column that contains information by which the data should be faceted into multiple plots.}

\item{title}{optional argument specifying the title of the volcano plot. Default is "Volcano plot".}

\item{x_axis_label}{optional argument specifying the x-axis label. Default is "log2(fold change)".}

\item{y_axis_label}{optional argument specifying the y-axis label. Default is "-log10(q-value)".}

\item{legend_label}{optional argument specifying the legend label. Default is "Target".}

\item{log2FC_cutoff}{optional argument specifying the log2 transformed fold change cutoff used for assessing whether changes are significant. Default value is 1.}

\item{significance_cutoff}{optional argument specifying the p-value cutoff used for assessing significance of changes. Default is 0.01.}

\item{interactive}{logical, indicating whether the plot should be interactive or not. Default is \code{interactive = FALSE}.}
}
\value{
Depending on the method used a volcano plot with either highlighted targets (\code{method = "target"}) or highlighted significant proteins (\code{method = "significant"}) is returned.
}
\description{
Plots a volcano plot for the given input.
}
\examples{
\dontrun{
volcano_protti(
  data,
  grouping = pep_stripped_sequence,
  log2FC = log2FC,
  significance = p_value,
  method = "target",
  target_column = uniprot_id,
  target = "Q9Y6K9",
  facet_by = comparison,
  title = "Finding Nemo",
  x_axis_label = "log2(fold change) treated vs untreated",
  y_axis_label = "-log10(p-value)",
  legend_label = "Target Protein",
  log2FC_cutoff = 2,
  significance_cutoff = 0.05,
  interactive = TRUE
)
}
}
