% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomise_queue.R
\name{randomise_queue}
\alias{randomise_queue}
\title{Randomise samples in MS queue}
\usage{
randomise_queue(data = NULL, rows = NULL, export = FALSE)
}
\arguments{
\item{data}{optional, a data frame containing a queue. If not provided a queue file can be chosen interactively.}

\item{rows}{optional, a range of rows in for which samples should be randomized.}

\item{export}{logical, if TRUE a \code{"randomised_queue.csv"} file will be saved in the working directory. If FALSE a data frame will be returned.}
}
\value{
If \code{export = TRUE} a \code{"randomised_queue.csv"} file will be saved in the working directory. If \code{export = FALSE} a data frame that contains the randomised queue is returned.
}
\description{
This function randomises the order of samples in an MS queue. QC and Blank samples are left in place. It is also possible to randomise only parts of the queue. The randomisation is every time the same for the same input due to \code{set.seed(123)}.
}
\examples{
\dontrun{
randomise_queue(rows = 195:235, export = TRUE)
}
}