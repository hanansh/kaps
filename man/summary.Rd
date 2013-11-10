\name{summary}
\docType{methods}
\alias{summary}
\alias{summary,kaps-method}
\title{Methods for Function \code{summary} for the S4 classes "kaps".}
\usage{
\S4method{summary}{kaps}(object, K)
}
\description{ This function provides tabloid information such as survival median, 1-, 3-, and 5 years actual suvival time for each terminal group as a table format.}
\arguments{
  \item{object}{an object of S4 class \code{kaps}.}
  \item{K}{a scalar object to determine the number of subgroups K. If missing, the estimated subgroup K is selected.}
}
\seealso{
  \code{\link{kaps}} and \code{\linkS4class{kaps}} class \cr
}
\keyword{methods}
