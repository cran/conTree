% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{air_quality}
\alias{air_quality}
\title{Air Quality Data from UC Irvine Machine Learning Repository}
\format{
\subsection{\code{air_quality}}{

A list with 4 items.
\describe{
\item{xco}{data frame of 9357 observations on 13 variables}
\item{yco}{hourly averaged CO concentration}
\item{zco}{sample membership indicator}
\item{pr2}{probability propensity score}
}
}
}
\source{
\url{https://archive.ics.uci.edu/ml/datasets/air+quality}
}
\usage{
air_quality
}
\description{
The data set consists of hourly averaged measurements from an array
of 5 metal oxide chemical sensors embedded in an air quality
chemical multisensor device. The outcome variable \code{y} is the
corresponding true hourly averaged concentration CO taken from a
reference analyzer. The input variables \code{x} are taken to be the
corresponding hourly averaged measurements of the 13 other
quantities.
}
\keyword{datasets}
