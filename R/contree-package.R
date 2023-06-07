#' Contrast and Boosted Trees
#'
#' Contrast trees represent a new approach for assessing the accuracy
#' of many types of machine learning estimates that are not amenable
#' to standard (cross) validation methods. In situations where
#' inaccuracies are detected, boosted contrast trees can often improve
#' performance. Functions are provided to to build such trees in
#' addition to a special case, distribution boosting, an assumption
#' free method for estimating the full probability distribution of an
#' outcome variable given any set of joint input predictor variable
#' values.
#'
#' @name conTree-package
#' @docType package
#' @useDynLib conTree
#' @author Original code (C) by Jerome H. Friedman, minor modifications,
#'     formatting, and packaging by Balasubramanian Narasimhan
#' @references Jerome Friedman (2019). _Contrast Trees and
#'     Distribution Boosting_ \url{https://arxiv.org/abs/1912.03785}
NULL

