\name{plot}
\docType{methods}
\alias{plot}
\alias{plot,kaps-method}
\title{Methods for Function \code{plot} of an S4 class object "kaps".}
\description{
  Plot an elbow polot for subgroups of the \code{"kaps"} object.
}
\usage{
	\S4method{plot}{kaps}(x, K, ...)
}
\arguments{
 \item{x}{an object of the S4 class \code{"kaps"}. }
 \item{K}{a scalar object that plots the kaplan meire survival curves for the proposed number of K. If missing, it works with selected K in the fitting model.}
 \item{\dots}{other arguments to the \code{\link[=graphics]{plot.default}} function can be passed here.}
}
\details{ 
	This function provides elbow plot and Kaplan-Meire survival curves for S4 formulation \code{"kaps"}.
}
\seealso{
    \code{\link{kaps}} and \code{\linkS4class{kaps}} class
}
\keyword{methods}
