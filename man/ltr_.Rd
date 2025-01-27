% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_operators.R
\name{ltr_}
\alias{ltr_}
\title{Extract values of matrix below the diagonal}
\usage{
ltr_(M)
}
\arguments{
\item{W}{A matrix}
}
\value{
A vector of values under the diagonal sorted as \code{\link[=lower.tri]{lower.tri()}} function does
}
\description{
Extract values of matrix below the diagonal
}
\examples{
data(leukemia)
M<- matrix(1:9, nrow=3)
ltr_(M)
}
