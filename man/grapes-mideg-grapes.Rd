% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_operators.R
\name{\%mideg\%}
\alias{\%mideg\%}
\title{Calculate mean inner degree of the cluster in network based ow weight matrix}
\usage{
W \%mideg\% mask
}
\arguments{
\item{W}{Numeric matrix}

\item{mask}{logical mask}
}
\value{
\code{sum(W[mask, mask])/sum(mask)}
}
\description{
Calculate mean inner degree of the cluster in network based ow weight matrix
}
