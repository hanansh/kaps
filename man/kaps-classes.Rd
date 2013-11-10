\name{kaps-class}
\Rdversion{1.1}
\docType{class}
\alias{kaps-class}
\title{ S4 Class \code{"kaps"}}
\description{
	S4 class of the object "kaps".
}
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("kaps", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{evaluated function call with the type of \code{"language"}  }
    \item{\code{formula}:}{formula to be used with the type of \code{"Formula"}  }
    \item{\code{data}:}{data to be used in the fitted model with the type of \code{"data.frame"}  }
    \item{\code{groupID}:}{subgroup classified with the type of \code{"vector"}  }
    \item{\code{index}:}{optimal K used in the fitted model with the type of \code{"integer"}  }
    \item{\code{X}:}{test statistic with the worst pair of subgroups for split s with the type of \code{"numeric"}  }
    \item{\code{Z}:}{log rank statistic with the type of \code{"numeric"} }
    \item{\code{pvalue}:}{p-value from permutation test with the type of \code{"numeric"} }
    \item{\code{WH}:}{Wilson-Hilferty approximation statistic with the type of \code{"numeric"} }
    \item{\code{t}:}{standardized value of WH statistic with the type of \code{"numeric"} }
    \item{\code{pair}:}{selected pair of subgroups for split s with the type of \code{"numeric"} }
    \item{\code{split.var}:}{selected covariate in the model formula with the type of \code{"character"} }
    \item{\code{split.pt}:}{set of cut-off points with the type of \code{"numeric"} }
    \item{\code{mindat}:}{minimum number of observations to be used (\code{"numeric"}) }
    \item{\code{elbow}:}{p-value matrix with the type of \code{"matrix"}. The first and third columns mean overall and the worst pairs's p-values, respectively. The second and fourth columns denotes their standard deviations.}
    \item{\code{over.stat.sample}:}{overall test statistics calculated by Bootstrap}
    \item{\code{pair.stat.sample}:}{worst pairwise test statistics calculated by Bootstrap}
    \item{\code{groups}:}{candidate K to be used with the type of \code{"vector"} }
    \item{\code{results}:}{results for each K with the type of \code{"list"} }
    \item{\code{Options}:}{minor parameters to be used with the type of \code{"kapsOptions"} }
	}
}
\section{Methods}{
  \describe{
	\item{show}{\code{signature(object = "kaps")}: Same as the show method without the optional arguments }
	\item{print}{\code{signature(object = "kaps", K)}: Same as the print method with the specified number of subgroups K.}
	\item{plot}{\code{signature(object = "kaps")}: Plot an object with kaps class on the current graphical device}
	\item{predict}{\code{signature(object = "kaps")}: Predict an object with kaps class by the estimated cut-off points}
	\item{summary}{\code{signature(object = "kaps")}: Summarize an object with kaps class by survival times for each subgroup }
	}
}
\author{
Soo-Heang Eo <hanansh@korea.ac.kr>
}
\examples{
showClass("kaps")
}
\keyword{classes}
