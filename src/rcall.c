#include <R.h>
#include <Rinternals.h>

/* The coding process is always: */
/* 1. Modify `rfcall` function below with appropriate arguments and matching code in delineated section. */

/* Then calling workflow has to always be: */
/* 1. Call `store_rfun` in R to save the user-specified function `rfun` before calling `.Fortran`. */
/* 2. Make the `.Fortran` fortran function call `rfcall` with appropriate arguments. */

/* The R function we save to avoid passing anything besides double and int items in fortran! */
static SEXP rfunc;

/* Store the R function */
SEXP store_rfun(SEXP rfun) {
  rfunc = rfun;
  return(R_NilValue);
}

void F77_NAME(rfcall)(int *n, double *y, double *z, double *w, double *result) {
  
  /* Start of section that customized for arguments passed to `rfcall` */
  /* Ensure you modify the arguments above, first */
  /* Then modify stuff below to match that. */

  int num_protected = 0;  
  SEXP ry, rz, rw;
  
  PROTECT(ry = allocVector(REALSXP, *n));
  num_protected++;
  PROTECT(rz = allocVector(REALSXP, *n));
  num_protected++;
  PROTECT(rw = allocVector(REALSXP, *n));
  num_protected++;
  
  double *yvec = REAL(ry); double *zvec = REAL(rz); double *wvec = REAL(rw);
  for (int i = 0; i < *n; i++) {
    yvec[i] = y[i]; zvec[i] = z[i]; wvec[i] = w[i];
  }
  /* End of section that customized for arguments passed to `rfcall` */

  SEXP rho = R_GetCurrentEnv();
  SEXP call = PROTECT(LCONS(rfunc, LCONS(ry, LCONS(rz, LCONS(rw, R_NilValue)))));
  num_protected++;
  
  SEXP r_result = R_forceAndCall(call, 3, rho);

  /* Start of further post-processing section. Modify as needed to check result */

  R_len_t len = length(r_result);
  if (len > 1) {
    error("R discrepancy function result length > 1");
  }

  *result = REAL(r_result)[0];

  /* End of further post-processing section. */  

  UNPROTECT(num_protected);
}

