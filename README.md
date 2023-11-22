# Contrast Trees and Distribution Boosting

<!-- badges: start -->
[![R-CMD-check](https://github.com/bnaras/conTree/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bnaras/conTree/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/conTree)](https://cran.r-project.org/package=conTree)
[![](https://cranlogs.r-pkg.org/badges/conTree)](https://CRAN.R-project.org/package=conTree)
<!-- badges: end -->

Contrast trees are used to assess the accuracy of many types of
machine learning estimates that are not amenable to standard
validation techniques. These include properties of the conditional
distribution $p_{y}(y\,|\,\mathbf{x})$ (means, quantiles, complete
distribution) as functions of $\mathbf{x}$. Given a set of predictor
variables $\mathbf{x}=(x_{1},x_{2},$$,x_{p})$ and two outcome
variables $y$ and $z$ associated with each $\mathbf{x}$, a contrast
tree attempts to partition the space of $\mathbf{x}$ values into local
regions within which the respective distributions of
$y\,|\,\mathbf{x}$ and $z\,|\,\mathbf{x}$, or selected properties of
those distributions such as means or quantiles, are most different.

For more details, please see the
[tutorial](https://jhfhub.github.io/conTree_tutorial/).


