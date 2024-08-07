% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_uniprot.R
\name{fetch_uniprot}
\alias{fetch_uniprot}
\title{Fetch protein data from UniProt}
\usage{
fetch_uniprot(
  uniprot_ids,
  columns = c("protein_name", "length", "sequence", "gene_names", "xref_geneid",
    "xref_string", "go_f", "go_p", "go_c", "cc_interaction", "ft_act_site", "ft_binding",
    "cc_cofactor", "cc_catalytic_activity", "xref_pdb"),
  batchsize = 200,
  max_tries = 10,
  timeout = 20,
  show_progress = TRUE
)
}
\arguments{
\item{uniprot_ids}{a character vector of UniProt accession numbers.}

\item{columns}{a character vector of metadata columns that should be imported from UniProt (all
possible columns can be found \href{https://www.uniprot.org/help/return_fields}{here}. For
cross-referenced database provide the database name with the prefix "xref_", e.g. \code{"xref_pdb"})}

\item{batchsize}{a numeric value that specifies the number of proteins processed in a single
single query. Default and max value is 200.}

\item{max_tries}{a numeric value that specifies the number of times the function tries to download
the data in case an error occurs.}

\item{timeout}{a numeric value that specifies the maximum request time per try. Default is 20 seconds.}

\item{show_progress}{a logical value that determines if a progress bar will be shown. Default
is TRUE.}
}
\value{
A data frame that contains all protein metadata specified in \code{columns} for the
proteins provided. The \code{input_id} column contains the provided UniProt IDs. If an invalid ID
was provided that contains a valid UniProt ID, the valid portion of the ID is still fetched and
present in the \code{accession} column, while the \code{input_id} column contains the original not completely
valid ID.
}
\description{
Fetches protein metadata from UniProt.
}
\examples{
\donttest{
fetch_uniprot(c("P36578", "O43324", "Q00796"))

# Not completely valid ID
fetch_uniprot(c("P02545", "P02545;P20700"))
}
}
