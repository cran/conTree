% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contrast.R
\name{ydist}
\alias{ydist}
\title{Transform z-values t(z) such that the distribution of \eqn{p(t(z) | x)} approximates \eqn{p(t(y | x)} for type = 'dist' only}
\usage{
ydist(model, x, z, num = model$niter)
}
\arguments{
\item{model}{model object output from modtrast()}

\item{x}{vector of predictor variable values for a (single) observation}

\item{z}{sample of z-values drawn from \eqn{p(z | x)}}

\item{num}{number of trees used to compute model values}
}
\value{
vector of \code{length(z)} containing transformed values t(z) approximating \eqn{p(y | x)}
}
\description{
Transform z-values t(z) such that the distribution of \eqn{p(t(z) | x)} approximates \eqn{p(t(y | x)} for type = 'dist' only
}
\seealso{
\code{\link[=contrast]{contrast()}}
}
