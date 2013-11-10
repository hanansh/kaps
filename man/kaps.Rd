\name{kaps}
\alias{kaps}
\title{K-adaptive partitioing for survival data}
\description{
    Conduct $K$-adaptive partitioning algorithm for censored survival data, especially with an application to SEER (Surveillance Epidemiology and End Results).
}
\usage{
	kaps(formula, data, K = 2:5, mindat, 
    type = c("perm", "Boot", "test", "NULL"), ...)
}
\arguments{
  \item{formula}{a Formula object with a response on the left hand side of the '~' operator, and the covariate terms on the right side. The response has to be a survival object with survival time and censoring status in the \link[survival:Surv]{Surv} function. For more details, see \link[Formula:Formula]{Formula} page. }
  \item{data}{a data frame with variables used in formula. It needs at least three variables including survival time, censoring status, and a covariate. Multivariate covariates can be supported with "+" sign.}
  \item{K}{
    the number of subgroups. The default value is 2:4. 
  }
  \item{type}{
    the type of algorithms for finding optimal number of subgroups K. 
  }
  \item{mindat}{
    the minimum number of observations at each subgroup. The default value is 5\% of data.
  }
  \item{\ldots}{
	a list of minor parameter with the type of S4 class, "kapsOptions". For more details, see \link{kaps.control} page.
  }
}
\details{
  The function \code{kaps()} provides some routines for K-adaptive partitioning algorithm which is desinged to classify cut-off vaules by the minimax rule and Bootstrap approaches. 
}
\value{
	the algorithm KAPS returns an object with class "kaps" with the below slots. 
    \item{\code{call}:}{evaluated function call with the the type of \code{"language"}  }
    \item{\code{formula}:}{formula to be used with the type of \code{"Formula"}  }
    \item{\code{data}:}{data to be used in the fitted model with the type of \code{"data.frame"}  }
    \item{\code{groupID}:}{subgroup classified with the type of \code{"vector"}  }
    \item{\code{index}:}{optimal K used in the fitted model with the type of \code{"integer"}  }
    \item{\code{X}:}{test statistic with the worst pair of subgroups for split s with the type of \code{"numeric"}  }
    \item{\code{Z}:}{log rank statistic with the type of \code{"numeric"} }
    \item{\code{pvalue}:}{p-value of the \code{X} with the type of \code{"numeric"} }
    \item{\code{WH}:}{Wilson-Hilferty approximation statistic with the type of \code{"numeric"} }
    \item{\code{t}:}{standardized value of WH statistic with the type of \code{"numeric"} }
    \item{\code{pair}:}{selected pair of subgroups for split s with the type of \code{"numeric"} }
    \item{\code{split.var}:}{selected covariate in the model formula with the type of \code{"character"} }
    \item{\code{split.pt}:}{set of cut-off points with the type of \code{"numeric"} }
    \item{\code{mindat}:}{minimum number of observations to be used (\code{"numeric"}) }
    \item{\code{elbow}:}{CV statistic matrix with the type of \code{"matrix"} }
    \item{\code{groups}:}{candidate K to be used with the type of \code{"vector"} }
    \item{\code{results}:}{results for each K with the type of \code{"list"} }
    \item{\code{Options}:}{minor parameters to be used with the type of \code{"apssOptions"} }
}
\author{
    Soo-Heang Eo \email{hanansh@korea.ac.kr} \cr
    HyungJun Cho \email{hj4cho@korea.ac.kr} \cr
}
\references{
Eo, SH, Hong, S-M and Cho, H (2013). K-adaptive partitioning for survival data with an application to SEER: The kaps add-on package for R, arXiv.
}
\seealso{
	 \code{\link{show}}, \code{\link{plot}}, \code{\link{predict}}, \code{\link{print}} and \code{\link{summary}} for the convenient use of kaps() \cr
	 \code{\link{kaps.control}} to control kaps() more detail \cr
}
\examples{
  \dontrun{
    data(toy)
    f <- Surv(time, status) ~ meta
    # Fit kaps algorithm without cross-validation.
    # It means the step to finding optimal K is not entered.
    fit1 <- kaps(f, data = toy, K = 3)

    # show the object of kaps (it contains apss S4 class)
    fit1

    # plot Kaplan-Meire estimates
    plot(fit1)

    # Fit kaps algorithm with cross-validation using V input argument.
    fit2 <- kaps(f, data = toy, K= 2:4, V=5) # V = 5 means that the number of fold is 5.
    fit2

    plot(fit2)
    print(fit2,K=2)
    summary(fit2)
    summary(fit2,K=2)

    # require(party)
    # fit4 <- ctree(f, data = toy)
  }
}
\keyword{kaps}
