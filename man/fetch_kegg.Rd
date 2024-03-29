% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_kegg.R
\name{fetch_kegg}
\alias{fetch_kegg}
\title{Fetch KEGG pathway data from KEGG}
\usage{
fetch_kegg(species)
}
\arguments{
\item{species}{a character value providing an abreviated species name. "hsa" for human, "eco"
for E. coli and "sce" for S. cerevisiae. Additional possible names can be found for
\href{https://www.genome.jp/kegg-bin/show_organism?category=Eukaryotes}{eukaryotes} and for
\href{https://www.genome.jp/kegg-bin/show_organism?category=Prokaryotes}{prokaryotes}.}
}
\value{
A data frame that contains gene IDs with corresponding pathway IDs and names for a
selected organism.
}
\description{
Fetches gene IDs and corresponding pathway IDs and names for the provided organism.
}
\examples{
\donttest{
kegg <- fetch_kegg(species = "hsa")

head(kegg)
}
}
