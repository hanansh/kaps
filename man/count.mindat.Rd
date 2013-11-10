\name{count.mindat}
\alias{count.mindat}
\title{Caculate the sample size given the number of partitions or subgroups}
\description{
   This function caculates the sample size given the number of partitions or subgroups.
}
\usage{
	count.mindat(formula, data, part = 10)
}

\arguments{
  \item{formula}{an object with the type of Formula. }
  \item{data}{ a dataset used in the Formula object.}
  \item{part}{a numeric object to determine the number of partitions or subgroups we want.}
}
\seealso{
	\code{\link{kaps}}
}
\examples{
	data(toy)
	count.mindat(Surv(time,staus) ~ meta, data = toy, part =5)
	count.mindat(Surv(time,staus) ~ meta, data = toy, part =10)
	count.mindat(Surv(time,staus) ~ meta, data = toy)
}
\keyword{kaps}
