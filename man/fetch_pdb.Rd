% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_pdb.R
\name{fetch_pdb}
\alias{fetch_pdb}
\title{Fetch structure information from RCSB}
\usage{
fetch_pdb(pdb_ids, batchsize = 200, show_progress = TRUE)
}
\arguments{
\item{pdb_ids}{a character vector of PDB identifiers.}

\item{batchsize}{numeric, specifying the number of structures to be processed in a single query. Default is 2000.}

\item{show_progress}{logical, if true, a progress bar will be shown. Default is TRUE.}
}
\value{
A data frame that contains all structure metadata for the PDB IDs provided.
}
\description{
Fetches structure metadata from RCSB. If you want to retrieve atom data such as positions, use the function \code{fetch_pdb_structure()}.
}
\examples{
\donttest{
head(fetch_pdb(c("6HG1", "1E9I", "6D3Q", "4JHW")))
}
}
