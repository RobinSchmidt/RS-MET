#pragma once

namespace LaPackCPP {


//=================================================================================================
// Blas

// explicit template instantiations - todo: move them into separate files, one file per
// datatype - maybe the particular type for which we instantiate the templates can be #defined
// the i just need to write the instantiation code once and copy/rename the files and just change 
// the #define - easier to maintain

// Level 1:

template int axpy(long int* n, double *da, double *dx, long int *incx,
  double *dy, long int *incy);

template double asum(integer *n, double *dx, integer *incx);

template int copy(long *n, double *dx, long *incx, double *dy, long *incy);

template double dot(integer *n, double *dx, integer *incx, double *dy, integer *incy);


template long iamax(long *n, double *dx, long *incx);

template int scal(long *n, double *da, double *dx, long *incx);

template int swap(long *n, double *dx, long *incx, double *dy, long *incy);

// Level 2:

//template int ger(integer *m, integer *n, doublereal *alpha, doublereal *x, integer *incx, 
//  doublereal *y, integer *incy, doublereal *a, integer *lda);

template int ger(long *m, long *n, double *alpha, double *x, long *incx, double *y, long *incy,
  double *a, long *lda);

template int gbmv(char *trans, integer *m, integer *n, integer *kl, integer *ku, double *alpha,
  double *a, integer *lda, double *x, integer *incx, double *beta, double *y, integer *incy,
  ftnlen trans_len);

template int gemv(char *trans, integer *m, integer *n, double *alpha, double *a, integer *lda,
  double *x, integer *incx, double *beta, double *y, integer *incy, ftnlen trans_len);

template int tbsv(char *uplo, char *trans, char *diag, integer *n, integer *k, double *a,
  integer *lda, double *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);


// Level 3:

template int gemm(char *transa, char *transb, integer *m, integer *n, integer *k, double *alpha,
  double *a, integer *lda, double *b, integer *ldb, double *beta, double *c__, integer *ldc,
  ftnlen transa_len, ftnlen transb_len);

template int trsm(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n,
  double *alpha, double *a, integer *lda, double *b, integer *ldb, ftnlen side_len,
  ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len);


//=================================================================================================
// LaPack


// Driver routines:

template int gbsv(long int *n, long int *kl, long int *ku, long int *nrhs, double *ab, 
  long int *ldab, long int *ipiv, double *b, long int *ldb, long int *info);

template int gbsvx(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, 
  doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, integer *ipiv, char *equed, 
  doublereal *r__, doublereal *c__, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
  doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
  integer *info, ftnlen fact_len, ftnlen trans_len, ftnlen equed_len);

template int gbsvxx(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, 
  doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, integer *ipiv, char *equed, 
  doublereal *r__, doublereal *c__, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
  doublereal *rcond, doublereal *rpvgrw, doublereal *berr, integer *n_err_bnds__, 
  doublereal *err_bnds_norm__, doublereal *err_bnds_comp__, integer *nparams, doublereal *params, 
  doublereal *work, integer *iwork, integer *info, ftnlen fact_len, ftnlen trans_len, 
  ftnlen equed_len);


// Computational routines:

template int gbequb(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab, 
  double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, integer *info);

template int gbrfsx(char *trans, char *equed, integer *n, integer *kl, integer *ku, integer *nrhs, 
  double *ab, integer *ldab, double *afb, integer *ldafb, integer *ipiv, double *r__, double *c__, 
  double *b, integer *ldb, double *x, integer *ldx, double *rcond, double *berr, 
  integer *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, integer *nparams, 
  double *params, double *work, integer *iwork, integer *info, ftnlen trans_len, ftnlen equed_len);

template int gbtrf(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab, 
  integer *ipiv, integer *info);

template int gbtrs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, double *ab,
  integer *ldab, integer *ipiv, double *b, integer *ldb, integer *info, ftnlen trans_len);


// Auxiliary routines:

template int latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, 
  double *ab, integer *ldab, double *x, double *scale, double *cnorm, integer *info, 
  ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len, ftnlen normin_len);

//template int rscl(integer *n, double *sa, double *sx, integer *incx);
// this instantiation gives a linker error with gcc: "no definition available" but leaving the 
// instantiation out gives another linker error: "undefined reference"
// ...why is no definition available? -> figure out, where we include the definitions
// maybe try this tool: https://include-what-you-use.org/

}

