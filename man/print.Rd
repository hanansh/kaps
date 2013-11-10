\name{print}
\docType{methods}
\alias{print}
\alias{print,kaps-method}
\title{Methods for Function \code{print} for the "kaps" class. }
\usage{
	\S4method{print}{kaps}(x, K)
}
\description{Print simple summary for \code{kaps} calss with the specified number of subgroups K}
\arguments{
  \item{x}{an object of S4 class \code{kaps}.}
  \item{K}{a scalar object to determine the number of subgroups K. If missing, the estimated subgroup K is selected.}
}
\seealso{
  \code{\link{kaps}} and \code{\linkS4class{kaps}} class
}
\keyword{methods}
