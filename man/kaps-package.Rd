\name{kaps-package}
\alias{kaps-package}
\docType{package}
\title{
	Functions for K-adaptive partitioning and recursive partitioning algorithms for censored survival data.
}
\description{
This package provides some routines for K-adaptive partitioning and recursive partitioning algorithms for censored survival data.
}
\details{
\tabular{ll}{
Package: \tab kaps\cr
Type: \tab Package\cr
Version: \tab 0.9.15\cr
Date: \tab 2013-10-20\cr
License: \tab GPL-3\cr
LazyLoad: \tab no\cr
}
This package contains some routines to conduct the K-adaptive parititioning (kaps) and recursive partitioning (lrtree) models for survival data, with the application to SEER (Surveillance Epidemiology and End Results). A function kaps is an implementation version of our algorithm kaps which provides minimax based partitioing rule based on the log-rank test. 
}
\author{
Soo-Heang Eo <hanansh@korea.ac.kr> and HyungJun Cho <hj4cho@korea.ac.kr>
}

\references{
SH Eo, S-M Hong and H Cho (2013). K-adaptive partitioning for survival data with an application to SEER: The kaps add-on package for R, \emph{arXiv}, 1306.\cr

S-M Hong, H Cho, C Moskaluk , E Yu (2007). Measurement of the Invasion Depth of Extrahepatic Bile Duct Carcinoma: An Alternative Method Overcoming the Current T Classiffication Problems of the AJCC Staging System, \emph{American Journal of Surgical Pathology}, 31, 199-206.

MR Segal (1988). Regression trees for censored survival data, \emph{Biometrics}, 44, 35-47.
}
\seealso{
    \code{\link{kaps}} 
}
\keyword{ package }