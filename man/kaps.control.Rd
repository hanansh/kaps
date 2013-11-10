\name{kaps.control}
\alias{kaps.control}
\title{Control minor parameters for "apss" object}
\description{
    Minor parameters for kaps()}
\usage{
	kaps.control(pre.pt = list(), scope = list(), 
  rho = 0, ncl = 1, lower.limit = 0, upper.limit = 100,
  shortcut = TRUE, N.perm = 1999, N.boot = 200, alpha = 0.05,
  splits = c("logrank", "exact"),
  boot.sel = c("mean","median","trim","test"),
  p.adjust.methods = c("none", "holm", "hochberg", 
    "hommel", "bonferroni",  "BH", "BY", "fdr"))
}
\arguments{
  \item{pre.pt}{a parameter with the type of list in order to treat pre-specified split candidates. Use the option as list(var = split points), i.e., x = 1:100}
  \item{scope}{a parameter with the type of list in order to treat pre-determined split range. Use the option as list(var = ranges), i.e., x = c(1,100)}
  \item{rho}{a scalar parameter that controls the type of logrank test. See \link[=survdiff]{survdiff}.}
  \item{ncl}{a scalar parameter with the type of integer to determine the number of cores used in parallel computing}
  \item{lower.limit}{a scalar parameter with the type of numeric to treat pre-determined overall lower bound. Default is 0.}
  \item{upper.limit}{a scalar parameter with the type of numeric to treat pre-determined overall upper bound. Default is 12.}
  \item{shortcut}{a logical value. If TRUE, we skip the off-diagonal matrix in pairwise-comparison to reduce computational cost. Default is TRUE.}
  \item{N.perm}{a scalar parameter with the type of numeric to give the number of permutation samples used in the kaps with single permutation, and bootstrap and permuting kaps algorithms. The default value is 1999.}
  \item{N.boot}{a scalar parameter with the type of numeric to give the number of Bootstrap samples used in the bootstrap and permuting kaps algorithm. The default value is 200.}
  \item{alpha}{a scalar parameter with the type of numeric to proviede significant in the process of Bootstrap and permuting}
  \item{splits}{a character value to determine the kind of pairwise test. Default is logrank test.}
  \item{boot.sel}{a character value to determine the kind of selection methods for the pairwise test stestatic in Bootstrap algorithm. The \code{mean} and \code{median} denote the Bootstrap meand and median, \code{trim} describes censoring-related trimmed mean, and \code{test} is the method using proportional test for the Bootsrap samples. In \code{test}, the null hypothesis is p = 0.05.}
  \item{p.adjust.methods}{a character parameter to select the criteria for the multiple comparison.}
  }
\seealso{
	\code{\link{kaps}}
}
\keyword{kaps}
