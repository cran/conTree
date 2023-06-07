## Builtin discrepancy codes for one sample
onesamp_types <- c(
    dist = 1L,
    diff = 3L,
    class = 4L,
    quant = 5L,
    prob = 7L,
    maxmean = 8L,
    diffmean = 7L
)

## Builtin discrepancy codes for two sample
twosamp_types <- c(
    dist = 10L,
    prob = 12L,
    diffmean = 12L,
    maxmean = 14L
)

## Parameters used in Fortran for replication/debugging purposes

#' Return the one sample parameters used in fortran discrepancy
#' functions
#' @return a named list for each of the types.
#' @description These functions are mostly useful when one wants to
#'     test one's own discrepancy function in R `f(y, z, w)` to
#'     determine if the results are correct. So a natural test
#'     is to experiment by programming one of the already implemented
#'     discrepancy functions in R. However, the Fortran
#'     implementations of such discrepancy measures use some
#'     parameters in the computations and therefore the returned
#'     results from a simple R implementation may not exactly
#'     match. Using these parameters, one can ensure that they
#'     do. These are to be interpreted as follows.  For one sample,
#'     the `type = "dist"` implementation in the package returns 0 if
#'     the length of `y` is less than `nmin` which is (100L). The `eps
#'     = 1.0e-5` parameter is used to ensure that the denominator in
#'     the formula for the Anderson-Darling statistic is at least
#'     `eps`. Next, for `type = "prob"`, if the length of the vector
#'     is less than `nmin = 20` the discrepancy is computed to be
#'     0. And so on. Refer to the R and Fortran source for further
#'     details as this is an advanced topic.
#' @rdname Fortran_details
#' @export
onesample_parameters  <- function() {
    list(
        dist = list(eps = 1.0e-5, nmin = 100),
        diff = list(),
        class = list(maxclass2 = 10000, nmin = 100, idum = 2),
        quant = list(nmin = 50),
        prob = list(nmin = 20),
        maxmean = list(nmin = 20, sml = -1.0e20),
        diffmean = list(nmin = 20)
    )
}

#' Return the two sample parameters used in fortran discrepancy
#' functions
#' @rdname Fortran_details
#' @export
twosample_parameters  <- function() {
    list(
        dist = list(eps = 1.0e-5, nmin = 100),
        prob = list(fmin = 20),
        diffmean = list(fmin = 20),
        maxmean = list(fmin = 20, sml = -1.0e20)
    )
}


