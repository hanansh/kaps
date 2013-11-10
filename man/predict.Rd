\name{predict}
\docType{methods}
\alias{predict}
\alias{predict,kaps-method}
\title{Predict new data witgh the fitted model \code{kaps} .}
\usage{
\S4method{predict}{kaps}(object, newdata, type = c("predict","kaps"))
}
\description{
This function produces predicted values, obtained by evaluationg the kaps in the frame newdata predict. 
}
\arguments{
  \item{object}{object of S4 class \code{kaps}.}
  \item{newdata}{An optinal data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}
  \item{type}{A logical argument. If "predict", predict an kaps class. The type of "kaps" is used for the internal useness.}
}
\seealso{
  \code{\link{kaps}}
}
\keyword{methods}
