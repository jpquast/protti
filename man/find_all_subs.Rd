% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_all_subs.R
\name{find_all_subs}
\alias{find_all_subs}
\title{Find all ChEBI sub IDs of an ID}
\usage{
find_all_subs(data, id, type = "is_a")
}
\arguments{
\item{data}{A data frame that contains information on ChEBI IDs (id), their sub IDs (incoming) and their relationship (type). This data frame can be obtained by calling \code{fetch_chebi(relation = TRUE)}.}

\item{id}{A character vector of ChEBI IDs for which sub IDs should be retreived.}

\item{type}{A character vector containing the type of relationship that should be considered for retreival. It is possible to use "all" relationships. The default type is "is_a". A list of possible relationships can be found \href{https://docs.google.com/document/d/1_w-DwBdCCOh1gMeeP6yqGzcnkpbHYOa3AGSODe5epcg/edit#heading=h.hnsqoqu978s5}{here}.}
}
\value{
A list of character vector containing the provided ID and all of its sub IDs. It contains one element per input ID.
}
\description{
For a given ChEBI ID, find all ChEBI sub IDs (incoming IDs) and their sub IDs. The type of relationship can be selected too. This is a helper function for other functions.
}
