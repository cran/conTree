% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contrast.R
\name{getnodes}
\alias{getnodes}
\title{Get terminal node observation assignments}
\usage{
getnodes(tree, x)
}
\arguments{
\item{tree}{model object output from contrast() or prune()}

\item{x}{training input predictor data matrix or data frame in same format as in contrast()}
}
\value{
vector of tree terminal node identifiers (numbers) corresponding to each observation (row of x)
}
\description{
Get terminal node observation assignments
}
\seealso{
\code{\link[=contrast]{contrast()}}
}
