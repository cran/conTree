% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{census}
\alias{census}
\title{Census Data Example from UC Irvine Machine Learning Repository}
\format{
\subsection{\code{census}}{

A list of 10 items.
\describe{
\item{x}{training data frame of 32561 observations on 12 predictor variables}
\item{y}{training binary response whether salary is above $50K or not}
\item{xt}{test data frame of 16281 observations predictor variables}
\item{yt}{test binary response whether salary is above $50K or not}
\item{gbl}{training GBL response variable}
\item{gblt}{test GBL response variable}
\item{gbp}{training GBP response variable}
\item{gbpt}{test GBP response variable}
\item{rf}{training RF response probabilities}
\item{rft}{test GBP response probabilities}
}
}
}
\source{
\url{https://archive.ics.uci.edu/ml/datasets/census+income}
}
\usage{
census
}
\description{
Includes a data frame of 1994 US census income from 48,842 people
divided into a training set of 32,561 and an independent test set
of 16,281. The training outcome variable \code{y} (\code{yt} for test) is
binary and indicates whether or not a person’s income is greater
than $50,000 per year. There are 12 predictor variables \code{x} (\code{xt}
for test) consisting of various demographic and financial
properties associated with each person. It also included estimates
of \eqn{Pr(y=1|x)} obtained by several machine learning methods:
gradient boosting on logistic scale using maximum likelihood (GBL),
random forest (RF), and gradient boosting on the probability scale
(GBP) using least–squares.
}
\keyword{datasets}
