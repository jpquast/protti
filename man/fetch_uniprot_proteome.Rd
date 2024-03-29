% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_uniprot_proteome.R
\name{fetch_uniprot_proteome}
\alias{fetch_uniprot_proteome}
\title{Fetch proteome data from UniProt}
\usage{
fetch_uniprot_proteome(organism_id, columns = c("accession"), reviewed = TRUE)
}
\arguments{
\item{organism_id}{a numeric value that specifies the NCBI taxonomy identifier (TaxId) for an
organism.}

\item{columns}{a character vector of metadata columns that should be imported from UniProt (all
possible columns can be found \href{https://www.uniprot.org/help/return_fields}{here}. For
cross-referenced database provide the database name with the prefix "xref_", e.g. \code{"xref_pdb"}).
Note: Not more than one or two columns should be selected otherwise the function will not be
able to efficiently retrieve the information. If more information is needed, \code{fetch_uniprot()}
can be used with the IDs retrieved by this function.}

\item{reviewed}{a logical value that determines if only reviewed protein entries will be retrieved.}
}
\value{
A data frame that contains all protein metadata specified in \code{columns} for the
organism of choice.
}
\description{
Fetches proteome data from UniProt for the provided organism ID.
}
\examples{
\donttest{
head(fetch_uniprot_proteome(9606))
}
}
