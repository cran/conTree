#' Census Data Example from UC Irvine Machine Learning Repository
#'
#' Includes a data frame of 1994 US census income from 48,842 people
#' divided into a training set of 32,561 and an independent test set
#' of 16,281. The training outcome variable `y` (`yt` for test) is
#' binary and indicates whether or not a person’s income is greater
#' than $50,000 per year. There are 12 predictor variables `x` (`xt`
#' for test) consisting of various demographic and financial
#' properties associated with each person. It also included estimates
#' of \eqn{Pr(y=1|x)} obtained by several machine learning methods:
#' gradient boosting on logistic scale using maximum likelihood (GBL),
#' random forest (RF), and gradient boosting on the probability scale
#' (GBP) using least–squares.
#'
#' @format ## `census`
#' A list of 10 items.
#' \describe{
#'    \item{x}{training data frame of 32561 observations on 12 predictor variables}
#'    \item{y}{training binary response whether salary is above $50K or not}
#'    \item{xt}{test data frame of 16281 observations predictor variables}
#'    \item{yt}{test binary response whether salary is above $50K or not}
#'    \item{gbl}{training GBL response variable}
#'    \item{gblt}{test GBL response variable}
#'    \item{gbp}{training GBP response variable}
#'    \item{gbpt}{test GBP response variable}
#'    \item{rf}{training RF response probabilities}
#'    \item{rft}{test GBP response probabilities}
#' }
#' @source <https://archive.ics.uci.edu/ml/datasets/census+income>
"census"

#' Age and Demographics data
#'
#' The data come from 9243 questionnaires filled out by shopping mall
#' customers in the San Francisco Bay Area (Impact Resources, Inc.,
#' Columbus, OH). Here we attempt to estimate a persons age as a
#' function of the other 13 demographic variables. For this data set
#' age value is reported as being in one of seven intervals `{13-17,
#' 18-24, 25-34, 35-44, 45-54, 55-64, >= 65}`. Each persons age is
#' randomly generated uniformly within its corresponding reported
#' interval. For the last interval an exponential distribution was
#' used with mean corresponding to life expectancy after reaching age
#' 65.
#'
#' @format ## `age_data`
#' A list of 3 items.
#' \describe{
#'    \item{xage}{data frame of 8856 observations on 13 variables}
#'    \item{yage}{Randomly generated age in the range above}
#'    \item{gbage}{gradient boosting model for median age given x}
#' }
#' @source The Elements of Statistical Learning, Data Mining, Second Edition, by Hastie,
#' Tibshirani, and Friedman.
"age_data"

#' Air Quality Data from UC Irvine Machine Learning Repository
#'
#' The data set consists of hourly averaged measurements from an array
#' of 5 metal oxide chemical sensors embedded in an air quality
#' chemical multisensor device. The outcome variable `y` is the
#' corresponding true hourly averaged concentration CO taken from a
#' reference analyzer. The input variables `x` are taken to be the
#' corresponding hourly averaged measurements of the 13 other
#' quantities.
#'
#' @format ## `air_quality`
#' A list with 4 items.
#' \describe{
#'    \item{xco}{data frame of 9357 observations on 13 variables}
#'    \item{yco}{hourly averaged CO concentration}
#'    \item{zco}{sample membership indicator}
#'    \item{pr2}{probability propensity score}
#' }
#' @source <https://archive.ics.uci.edu/ml/datasets/air+quality>
"air_quality"
