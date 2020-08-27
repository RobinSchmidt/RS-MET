
namespace LaPackCPP {

// some fiddling to make it compile and link - those functions are defined outside the LaPackCPP
// namespace, so in order to not have to enter a :: each time a function is called, we use this
// wrappers/delegators here inside the namspace - todo: clean this up! try to get rid and if
// impossible, at least inline them and maybe move to some other file

// some silly headers define these as macros - leading to utterly inscrutable compiler errors of
// the forms
// 'char' should be preceded by ';'
// 'double' followed by 'char' is illegal
// etc.
#undef min
#undef max
#undef small

//template<class T> inline T min(T x, T y) {  return std::min(x, y); }
//template<class T> inline T max(T x, T y) {  return std::max(x, y); }


// maybe move to Blas.cpp (done for min/max)
//inline long min(long x, long y) {  return std::min(x, y); }
//inline long max(long x, long y) {  return std::max(x, y); }
//
//inline double min(double x, double y) {  return std::min(x, y); }
//inline double max(double x, double y) {  return std::max(x, y); }

inline double sqrt(double x)  { return ::sqrt(x); }

template<class T> inline T log(T x)  {  return ::log((T)x); }
// do we need this?

/*
integer i_len(char *s, ftnlen n)
{
  //return ::i_len(s, n);
  return LibF2C::i_len(s, n);
}
integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
{
  return LibF2C::s_cmp(a0, b0, la, lb);
}
int s_copy(register char *a, register char *b, ftnlen la, ftnlen lb)
{
  LibF2C::s_copy(a, b, la, lb);
  return 0;
}
integer i_nint(f2c_real *x)
{
  return LibF2C::i_nint(x);
}
*/


blas_trans_type toTransType(integer* value)
{
  if( *value == 'T' ) return blas_trans;
  if( *value == 'C' ) return blas_conj_trans;
  else                return blas_no_trans;       // 'N'
}

blas_prec_type toPrecType(integer* value)
{
  if( *value == 'S' ) return blas_prec_single;
  if( *value == 'D' ) return blas_prec_double;
  if( *value == 'I' ) return blas_prec_indigenous;
  else                return blas_prec_extra;       // 'X'
}
// maybe raise assertion when the value is not one of the predefined ones

//=================================================================================================
// DRIVER routines

// translated from dgbsv, LAPACK driver routine (version 3.7.0) --
template<class T>
int gbsv(long int *n, long int *kl, long int *ku, long int *nrhs, T *ab, long int *ldab,
  long int *ipiv, T *b, long int *ldb, long int *info)
{
  // System generated locals
  long int ab_dim1, ab_offset, b_dim1, b_offset, i__1;

  // Local variables
  // Subroutine
  //extern  int
  //  gbtrf(long int *, long int *, long int *,
  //    long int *, double *, long int *, long int *, long int *),
  //  xerbla_(char *, long int *, ftnlen),
  //  gbtrs(char *, long int *,
  //    long int *, long int *, long int *, double *, long int *, long int *, double *,
  //    long int *, long int *, ftnlen);

  // Parameter adjustments
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --ipiv;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;

  // Function Body
  *info = 0;
  if(*n < 0) {
    *info = -1;
  }
  else if(*kl < 0) {
    *info = -2;
  }
  else if(*ku < 0) {
    *info = -3;
  }
  else if(*nrhs < 0) {
    *info = -4;
  }
  else if(*ldab < (*kl << 1) + *ku + 1) {
    *info = -6;
  }
  else if(*ldb < max(*n, 1)) {
    *info = -9;
  }
  if(*info != 0) {
    i__1 = -(*info);
    xerbla("DGBSV ", &i__1, (ftnlen)6);
    return 0;
  }

  // Compute the LU factorization of the band matrix A:
  gbtrf(n, n, kl, ku, &ab[ab_offset], ldab, &ipiv[1], info);
  if(*info == 0) {
    // Solve the system A*X = B, overwriting B with X:
    gbtrs("No transpose", n, kl, ku, nrhs, &ab[ab_offset], ldab, &ipiv[1], &b[b_offset], ldb,
      info, (ftnlen)12);
  }
  return 0;
}

//-------------------------------------------------------------------------------------------------

// translated from dgbsvx - LAPACK driver routine (version 3.7.0)
template<class T>
int gbsvx(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs,
  T *ab, integer *ldab, T *afb, integer *ldafb, integer *ipiv, char *equed, T *r__, T *c__, T *b,
  integer *ldb, T *x, integer *ldx, T *rcond, T *ferr, T *berr, T *work, integer *iwork,
  integer *info, ftnlen fact_len, ftnlen trans_len, ftnlen equed_len)
{
  /* Table of constant values */
  static integer c__1 = 1;

  /* System generated locals */
  integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset,
    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
  T d__1, d__2, d__3;

  /* Local variables */
  static integer i__, j, j1, j2;
  static T amax;
  static char norm[1];
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static T rcmin, rcmax, anorm;
  //extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *,
  //  doublereal *, integer *);
  static logical equil;
  //extern doublereal dlangb_(char *, integer *, integer *, integer *,
  //  doublereal *, integer *, doublereal *, ftnlen), dlamch_(char *,
  //    ftnlen);
  //extern /* Subroutine */ int dlaqgb_(integer *, integer *, integer *,
  //  integer *, doublereal *, integer *, doublereal *, doublereal *,
  //  doublereal *, doublereal *, doublereal *, char *, ftnlen),
  //  dgbcon_(char *, integer *, integer *, integer *, doublereal *,
  //    integer *, integer *, doublereal *, doublereal *, doublereal *,
  //    integer *, integer *, ftnlen);
  static T colcnd;
  //extern doublereal dlantb_(char *, char *, char *, integer *, integer *,
  //  doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
  //extern /* Subroutine */ int dgbequ_(integer *, integer *, integer *,
  //  integer *, doublereal *, integer *, doublereal *, doublereal *,
  //  doublereal *, doublereal *, doublereal *, integer *), dgbrfs_(
  //    char *, integer *, integer *, integer *, integer *, doublereal *,
  //    integer *, doublereal *, integer *, integer *, doublereal *,
  //    integer *, doublereal *, integer *, doublereal *, doublereal *,
  //    doublereal *, integer *, integer *, ftnlen), dgbtrf_(integer *,
  //      integer *, integer *, integer *, doublereal *, integer *, integer
  //      *, integer *);
  static logical nofact;
  //extern /* Subroutine */ int dlacpy_(char *, integer *, integer *,
  //  doublereal *, integer *, doublereal *, integer *, ftnlen),
  //  xerbla_(char *, integer *, ftnlen);
  static T bignum;
  //extern /* Subroutine */ int dgbtrs_(char *, integer *, integer *, integer
  //  *, integer *, doublereal *, integer *, integer *, doublereal *,
  //  integer *, integer *, ftnlen);
  static integer infequ;
  static logical colequ;
  static T rowcnd;
  static logical notran;
  static T smlnum;
  static logical rowequ;
  static T rpvgrw;

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  afb_dim1 = *ldafb;
  afb_offset = 1 + afb_dim1;
  afb -= afb_offset;
  --ipiv;
  --r__;
  --c__;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;
  x_dim1 = *ldx;
  x_offset = 1 + x_dim1;
  x -= x_offset;
  --ferr;
  --berr;
  --work;
  --iwork;

  /* Function Body */
  *info = 0;
  nofact = lsame(fact, "N", (ftnlen)1, (ftnlen)1);
  equil = lsame(fact, "E", (ftnlen)1, (ftnlen)1);
  notran = lsame(trans, "N", (ftnlen)1, (ftnlen)1);
  if (nofact || equil) {
    *(unsigned char *)equed = 'N';
    rowequ = FALSE_;
    colequ = FALSE_;
  } else {
    rowequ = lsame(equed, "R", (ftnlen)1, (ftnlen)1) || lsame(equed,
      "B", (ftnlen)1, (ftnlen)1);
    colequ = lsame(equed, "C", (ftnlen)1, (ftnlen)1) || lsame(equed,
      "B", (ftnlen)1, (ftnlen)1);
    smlnum = lamch("Safe minimum", (ftnlen)12);
    bignum = 1. / smlnum;
  }

  /*     Test the input parameters. */

  if (! nofact && ! equil && ! lsame(fact, "F", (ftnlen)1, (ftnlen)1)) {
    *info = -1;
  } else if (! notran && ! lsame(trans, "T", (ftnlen)1, (ftnlen)1) && !
    lsame(trans, "C", (ftnlen)1, (ftnlen)1)) {
    *info = -2;
  } else if (*n < 0) {
    *info = -3;
  } else if (*kl < 0) {
    *info = -4;
  } else if (*ku < 0) {
    *info = -5;
  } else if (*nrhs < 0) {
    *info = -6;
  } else if (*ldab < *kl + *ku + 1) {
    *info = -8;
  } else if (*ldafb < (*kl << 1) + *ku + 1) {
    *info = -10;
  } else if (lsame(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ
    || lsame(equed, "N", (ftnlen)1, (ftnlen)1))) {
    *info = -12;
  } else {
    if (rowequ) {
      rcmin = bignum;
      rcmax = 0.;
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        /* Computing MIN */
        d__1 = rcmin, d__2 = r__[j];
        rcmin = min(d__1,d__2);
        /* Computing MAX */
        d__1 = rcmax, d__2 = r__[j];
        rcmax = max(d__1,d__2);
        /* L10: */
      }
      if (rcmin <= 0.) {
        *info = -13;
      } else if (*n > 0) {
        rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
      } else {
        rowcnd = 1.;
      }
    }
    if (colequ && *info == 0) {
      rcmin = bignum;
      rcmax = 0.;
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        /* Computing MIN */
        d__1 = rcmin, d__2 = c__[j];
        rcmin = min(d__1,d__2);
        /* Computing MAX */
        d__1 = rcmax, d__2 = c__[j];
        rcmax = max(d__1,d__2);
        /* L20: */
      }
      if (rcmin <= 0.) {
        *info = -14;
      } else if (*n > 0) {
        colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
      } else {
        colcnd = 1.;
      }
    }
    if (*info == 0) {
      if (*ldb < max(1,*n)) {
        *info = -16;
      } else if (*ldx < max(1,*n)) {
        *info = -18;
      }
    }
  }

  if (*info != 0) {
    i__1 = -(*info);
    xerbla("DGBSVX", &i__1, (ftnlen)6);
    return 0;
  }

  if (equil) {

    /*        Compute row and column scalings to equilibrate the matrix A. */

    gbequ(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &rowcnd,
      &colcnd, &amax, &infequ);
    if (infequ == 0) {

      /*           Equilibrate the matrix. */

      laqgb(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &
        rowcnd, &colcnd, &amax, equed, (ftnlen)1);
      rowequ = lsame(equed, "R", (ftnlen)1, (ftnlen)1) || lsame(equed,
        "B", (ftnlen)1, (ftnlen)1);
      colequ = lsame(equed, "C", (ftnlen)1, (ftnlen)1) || lsame(equed,
        "B", (ftnlen)1, (ftnlen)1);
    }
  }

  /*     Scale the right hand side. */

  if (notran) {
    if (rowequ) {
      i__1 = *nrhs;
      for (j = 1; j <= i__1; ++j) {
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          b[i__ + j * b_dim1] = r__[i__] * b[i__ + j * b_dim1];
          /* L30: */
        }
        /* L40: */
      }
    }
  } else if (colequ) {
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        b[i__ + j * b_dim1] = c__[i__] * b[i__ + j * b_dim1];
        /* L50: */
      }
      /* L60: */
    }
  }

  if (nofact || equil) {

    /*        Compute the LU factorization of the band matrix A. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      /* Computing MAX */
      i__2 = j - *ku;
      j1 = max(i__2,1);
      /* Computing MIN */
      i__2 = j + *kl;
      j2 = min(i__2,*n);
      i__2 = j2 - j1 + 1;
      copy(&i__2, &ab[*ku + 1 - j + j1 + j * ab_dim1], &c__1, &afb[*
        kl + *ku + 1 - j + j1 + j * afb_dim1], &c__1);
      /* L70: */
    }

    gbtrf(n, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], info);

    /*        Return if INFO is non-zero. */

    if (*info > 0) {

      /*           Compute the reciprocal pivot growth factor of the */
      /*           leading rank-deficient INFO columns of A. */

      anorm = 0.;
      i__1 = *info;
      for (j = 1; j <= i__1; ++j) {
        /* Computing MAX */
        i__2 = *ku + 2 - j;
        /* Computing MIN */
        i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
        i__3 = min(i__4,i__5);
        for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
          /* Computing MAX */
          d__2 = anorm, d__3 = (d__1 = ab[i__ + j * ab_dim1], abs(
            d__1));
          anorm = max(d__2,d__3);
          /* L80: */
        }
        /* L90: */
      }
      /* Computing MIN */
      i__3 = *info - 1, i__2 = *kl + *ku;
      i__1 = min(i__3,i__2);
      /* Computing MAX */
      i__4 = 1, i__5 = *kl + *ku + 2 - *info;
      rpvgrw = lantb("M", "U", "N", info, &i__1, &afb[max(i__4,i__5)
        + afb_dim1], ldafb, &work[1], (ftnlen)1, (ftnlen)1, (
          ftnlen)1);
      if (rpvgrw == 0.) {
        rpvgrw = 1.;
      } else {
        rpvgrw = anorm / rpvgrw;
      }
      work[1] = rpvgrw;
      *rcond = 0.;
      return 0;
    }
  }

  /*     Compute the norm of the matrix A and the */
  /*     reciprocal pivot growth factor RPVGRW. */

  if (notran) {
    *(unsigned char *)norm = '1';
  } else {
    *(unsigned char *)norm = 'I';
  }
  anorm = langb(norm, n, kl, ku, &ab[ab_offset], ldab, &work[1], (ftnlen)
    1);
  i__1 = *kl + *ku;
  rpvgrw = lantb("M", "U", "N", n, &i__1, &afb[afb_offset], ldafb, &work[
    1], (ftnlen)1, (ftnlen)1, (ftnlen)1);
  if (rpvgrw == 0.) {
    rpvgrw = 1.;
  } else {
    rpvgrw = langb("M", n, kl, ku, &ab[ab_offset], ldab, &work[1], (
      ftnlen)1) / rpvgrw;
  }

  /*     Compute the reciprocal of the condition number of A. */

  gbcon(norm, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], &anorm, rcond,
    &work[1], &iwork[1], info, (ftnlen)1);

  /*     Compute the solution matrix X. */

  lacpy("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
  gbtrs(trans, n, kl, ku, nrhs, &afb[afb_offset], ldafb, &ipiv[1], &x[
    x_offset], ldx, info, (ftnlen)1);

  /*     Use iterative refinement to improve the computed solution and */
  /*     compute error bounds and backward error estimates for it. */

  gbrfs(trans, n, kl, ku, nrhs, &ab[ab_offset], ldab, &afb[afb_offset],
    ldafb, &ipiv[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &
    berr[1], &work[1], &iwork[1], info, (ftnlen)1);

  /*     Transform the solution matrix X to a solution of the original */
  /*     system. */

  if (notran) {
    if (colequ) {
      i__1 = *nrhs;
      for (j = 1; j <= i__1; ++j) {
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__) {
          x[i__ + j * x_dim1] = c__[i__] * x[i__ + j * x_dim1];
          /* L100: */
        }
        /* L110: */
      }
      i__1 = *nrhs;
      for (j = 1; j <= i__1; ++j) {
        ferr[j] /= colcnd;
        /* L120: */
      }
    }
  } else if (rowequ) {
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
      i__3 = *n;
      for (i__ = 1; i__ <= i__3; ++i__) {
        x[i__ + j * x_dim1] = r__[i__] * x[i__ + j * x_dim1];
        /* L130: */
      }
      /* L140: */
    }
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
      ferr[j] /= rowcnd;
      /* L150: */
    }
  }

  /*     Set INFO = N+1 if the matrix is singular to working precision. */

  if (*rcond < lamch("Epsilon", (ftnlen)7)) {
    *info = *n + 1;
  }

  work[1] = rpvgrw;
  return 0;

} /* gbsvx */


//-------------------------------------------------------------------------------------------------

// from  dgbsvxx -- LAPACK driver routine (version 3.7.0)
template<class T>
int gbsvxx(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab,
  integer *ldab, T *afb, integer *ldafb, integer *ipiv, char *equed, T *r__, T *c__, T *b,
  integer *ldb, T *x, integer *ldx, T *rcond, T *rpvgrw, T *berr, integer *n_err_bnds__,
  T *err_bnds_norm__, T *err_bnds_comp__, integer *nparams, T *params, T *work, integer *iwork,
  integer *info, ftnlen fact_len, ftnlen trans_len, ftnlen equed_len)
{
  /* System generated locals */
  integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset,
    x_dim1, x_offset, err_bnds_norm_dim1, err_bnds_norm_offset,
    err_bnds_comp_dim1, err_bnds_comp_offset, i__1, i__2;
  doublereal d__1, d__2;

  /* Local variables */
  static integer i__, j;
  static T amax;
  //extern doublereal dla_gbrpvgrw__(integer *, integer *, integer *, integer
  //  *, doublereal *, integer *, doublereal *, integer *);
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static T rcmin, rcmax;
  static logical equil;
  //extern doublereal dlamch_(char *, ftnlen);
  //extern /* Subroutine */ int dlaqgb_(integer *, integer *, integer *,
  //  integer *, doublereal *, integer *, doublereal *, doublereal *,
  //  doublereal *, doublereal *, doublereal *, char *, ftnlen);
  static T colcnd;
  //extern /* Subroutine */ int dgbtrf_(integer *, integer *, integer *,
  //  integer *, doublereal *, integer *, integer *, integer *);
  static logical nofact;
  //extern /* Subroutine */ int dlacpy_(char *, integer *, integer *,
  //  doublereal *, integer *, doublereal *, integer *, ftnlen),
  //  xerbla_(char *, integer *, ftnlen);
  static T bignum;
  //extern /* Subroutine */ int dgbtrs_(char *, integer *, integer *, integer
  //  *, integer *, doublereal *, integer *, integer *, doublereal *,
  //  integer *, integer *, ftnlen);
  static integer infequ;
  static logical colequ;
  static T rowcnd;
  static logical notran;
  static T smlnum;
  static logical rowequ;
  //extern /* Subroutine */ int dlascl2_(integer *, integer *, doublereal *,
  //  doublereal *, integer *), dgbequb_(integer *, integer *, integer *
  //    , integer *, doublereal *, integer *, doublereal *, doublereal *,
  //    doublereal *, doublereal *, doublereal *, integer *), dgbrfsx_(
  //      char *, char *, integer *, integer *, integer *, integer *,
  //      doublereal *, integer *, doublereal *, integer *, integer *,
  //      doublereal *, doublereal *, doublereal *, integer *, doublereal *,
  //      integer *, doublereal *, doublereal *, integer *, doublereal *,
  //      doublereal *, integer *, doublereal *, doublereal *, integer *,
  //      integer *, ftnlen, ftnlen);



  /* Parameter adjustments */
  err_bnds_comp_dim1 = *nrhs;
  err_bnds_comp_offset = 1 + err_bnds_comp_dim1;
  err_bnds_comp__ -= err_bnds_comp_offset;
  err_bnds_norm_dim1 = *nrhs;
  err_bnds_norm_offset = 1 + err_bnds_norm_dim1;
  err_bnds_norm__ -= err_bnds_norm_offset;
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  afb_dim1 = *ldafb;
  afb_offset = 1 + afb_dim1;
  afb -= afb_offset;
  --ipiv;
  --r__;
  --c__;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;
  x_dim1 = *ldx;
  x_offset = 1 + x_dim1;
  x -= x_offset;
  --berr;
  --params;
  --work;
  --iwork;

  /* Function Body */
  *info = 0;
  nofact = lsame(fact, "N", (ftnlen)1, (ftnlen)1);
  equil = lsame(fact, "E", (ftnlen)1, (ftnlen)1);
  notran = lsame(trans, "N", (ftnlen)1, (ftnlen)1);
  smlnum = lamch("Safe minimum", (ftnlen)12);
  bignum = 1. / smlnum;
  if (nofact || equil) {
    *(unsigned char *)equed = 'N';
    rowequ = FALSE_;
    colequ = FALSE_;
  } else {
    rowequ = lsame(equed, "R", (ftnlen)1, (ftnlen)1) || lsame(equed,
      "B", (ftnlen)1, (ftnlen)1);
    colequ = lsame(equed, "C", (ftnlen)1, (ftnlen)1) || lsame(equed,
      "B", (ftnlen)1, (ftnlen)1);
  }

  /*     Default is failure.  If an input parameter is wrong or */
  /*     factorization fails, make everything look horrible.  Only the */
  /*     pivot growth is set here, the rest is initialized in DGBRFSX. */

  *rpvgrw = 0.;

  /*     Test the input parameters.  PARAMS is not tested until DGBRFSX. */

  if (! nofact && ! equil && ! lsame(fact, "F", (ftnlen)1, (ftnlen)1)) {
    *info = -1;
  } else if (! notran && ! lsame(trans, "T", (ftnlen)1, (ftnlen)1) && !
    lsame(trans, "C", (ftnlen)1, (ftnlen)1)) {
    *info = -2;
  } else if (*n < 0) {
    *info = -3;
  } else if (*kl < 0) {
    *info = -4;
  } else if (*ku < 0) {
    *info = -5;
  } else if (*nrhs < 0) {
    *info = -6;
  } else if (*ldab < *kl + *ku + 1) {
    *info = -8;
  } else if (*ldafb < (*kl << 1) + *ku + 1) {
    *info = -10;
  } else if (lsame(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ
    || lsame(equed, "N", (ftnlen)1, (ftnlen)1))) {
    *info = -12;
  } else {
    if (rowequ) {
      rcmin = bignum;
      rcmax = 0.;
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        /* Computing MIN */
        d__1 = rcmin, d__2 = r__[j];
        rcmin = min(d__1,d__2);
        /* Computing MAX */
        d__1 = rcmax, d__2 = r__[j];
        rcmax = max(d__1,d__2);
        /* L10: */
      }
      if (rcmin <= 0.) {
        *info = -13;
      } else if (*n > 0) {
        rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
      } else {
        rowcnd = 1.;
      }
    }
    if (colequ && *info == 0) {
      rcmin = bignum;
      rcmax = 0.;
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        /* Computing MIN */
        d__1 = rcmin, d__2 = c__[j];
        rcmin = min(d__1,d__2);
        /* Computing MAX */
        d__1 = rcmax, d__2 = c__[j];
        rcmax = max(d__1,d__2);
        /* L20: */
      }
      if (rcmin <= 0.) {
        *info = -14;
      } else if (*n > 0) {
        colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
      } else {
        colcnd = 1.;
      }
    }
    if (*info == 0) {
      if (*ldb < max(1,*n)) {
        *info = -15;
      } else if (*ldx < max(1,*n)) {
        *info = -16;
      }
    }
  }

  if (*info != 0) {
    i__1 = -(*info);
    xerbla("DGBSVXX", &i__1, (ftnlen)7);
    return 0;
  }

  if (equil) {

    /*     Compute row and column scalings to equilibrate the matrix A. */

    gbequb(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &
      rowcnd, &colcnd, &amax, &infequ);
    if (infequ == 0) {

      /*     Equilibrate the matrix. */

      laqgb(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &
        rowcnd, &colcnd, &amax, equed, (ftnlen)1);
      rowequ = lsame(equed, "R", (ftnlen)1, (ftnlen)1) || lsame(equed,
        "B", (ftnlen)1, (ftnlen)1);
      colequ = lsame(equed, "C", (ftnlen)1, (ftnlen)1) || lsame(equed,
        "B", (ftnlen)1, (ftnlen)1);
    }

    /*     If the scaling factors are not applied, set them to 1.0. */

    if (! rowequ) {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        r__[j] = 1.;
      }
    }
    if (! colequ) {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        c__[j] = 1.;
      }
    }
  }

  /*     Scale the right hand side. */

  if (notran) {
    if (rowequ) {
      lascl2(n, nrhs, &r__[1], &b[b_offset], ldb);
    }
  } else {
    if (colequ) {
      lascl2(n, nrhs, &c__[1], &b[b_offset], ldb);
    }
  }

  if (nofact || equil) {

    /*        Compute the LU factorization of A. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      i__2 = (*kl << 1) + *ku + 1;
      for (i__ = *kl + 1; i__ <= i__2; ++i__) {
        afb[i__ + j * afb_dim1] = ab[i__ - *kl + j * ab_dim1];
        /* L30: */
      }
      /* L40: */
    }
    gbtrf(n, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], info);

    /*        Return if INFO is non-zero. */

    if (*info > 0) {

      /*           Pivot in column INFO is exactly 0 */
      /*           Compute the reciprocal pivot growth factor of the */
      /*           leading rank-deficient INFO columns of A. */

      *rpvgrw = la_gbrpvgrw(n, kl, ku, info, &ab[ab_offset], ldab, &
        afb[afb_offset], ldafb);
      return 0;
    }
  }

  /*     Compute the reciprocal pivot growth factor RPVGRW. */

  *rpvgrw = la_gbrpvgrw(n, kl, ku, n, &ab[ab_offset], ldab, &afb[
    afb_offset], ldafb);

  /*     Compute the solution matrix X. */

  lacpy("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
  gbtrs(trans, n, kl, ku, nrhs, &afb[afb_offset], ldafb, &ipiv[1], &x[
    x_offset], ldx, info, (ftnlen)1);

  /*     Use iterative refinement to improve the computed solution and */
  /*     compute error bounds and backward error estimates for it. */

  gbrfsx(trans, equed, n, kl, ku, nrhs, &ab[ab_offset], ldab, &afb[
    afb_offset], ldafb, &ipiv[1], &r__[1], &c__[1], &b[b_offset], ldb,
      &x[x_offset], ldx, rcond, &berr[1], n_err_bnds__, &
      err_bnds_norm__[err_bnds_norm_offset], &err_bnds_comp__[
        err_bnds_comp_offset], nparams, &params[1], &work[1], &iwork[1],
          info, (ftnlen)1, (ftnlen)1);

  /*     Scale solutions. */

  if (colequ && notran) {
    lascl2(n, nrhs, &c__[1], &x[x_offset], ldx);
  } else if (rowequ && ! notran) {
    lascl2(n, nrhs, &r__[1], &x[x_offset], ldx);
  }

  return 0;

  /*     End of DGBSVXX */

} /* dgbsvxx_ */

//=================================================================================================
// COMPUTATIONAL routines

// from dgbcon - LAPACK computational routine (version 3.7.0)
template<class T>
int gbcon(char *norm, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, integer *ipiv,
  T *anorm, T *rcond, T *work, integer *iwork, integer *info, ftnlen norm_len)
{
  /* Table of constant values */
  static integer c__1 = 1;

  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3;
  T d__1;

  /* Local variables */
  static integer j;
  static T t;
  static integer kd, lm, jp, ix, kase;
  //extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *,
  //  integer *);
  static integer kase1;
  static T scale;
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static integer isave[3];
  //extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *,
  //  integer *);
  static logical lnoti;
  //extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *,
  //  integer *, doublereal *, integer *), dlacn2_(integer *,
  //    doublereal *, doublereal *, integer *, doublereal *, integer *,
  //    integer *);
  //extern doublereal dlamch_(char *, ftnlen);
  //extern integer idamax_(integer *, doublereal *, integer *);
  //extern /* Subroutine */ int dlatbs_(char *, char *, char *, char *,
  //  integer *, integer *, doublereal *, integer *, doublereal *,
  //  doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
  //  ftnlen), xerbla_(char *, integer *, ftnlen);
  static T ainvnm;
  static logical onenrm;
  static char normin[1];
  static T smlnum;

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --ipiv;
  --work;
  --iwork;

  /* Function Body */
  *info = 0;
  onenrm = *(unsigned char *)norm == '1' || lsame(norm, "O", (ftnlen)1, (
    ftnlen)1);
  if (! onenrm && ! lsame(norm, "I", (ftnlen)1, (ftnlen)1)) {
    *info = -1;
  } else if (*n < 0) {
    *info = -2;
  } else if (*kl < 0) {
    *info = -3;
  } else if (*ku < 0) {
    *info = -4;
  } else if (*ldab < (*kl << 1) + *ku + 1) {
    *info = -6;
  } else if (*anorm < 0.) {
    *info = -8;
  }
  if (*info != 0) {
    i__1 = -(*info);
    xerbla("DGBCON", &i__1, (ftnlen)6);
    return 0;
  }

  /*     Quick return if possible */

  *rcond = 0.;
  if (*n == 0) {
    *rcond = 1.;
    return 0;
  } else if (*anorm == 0.) {
    return 0;
  }

  smlnum = lamch("Safe minimum", (ftnlen)12);

  /*     Estimate the norm of inv(A). */

  ainvnm = 0.;
  *(unsigned char *)normin = 'N';
  if (onenrm) {
    kase1 = 1;
  } else {
    kase1 = 2;
  }
  kd = *kl + *ku + 1;
  lnoti = *kl > 0;
  kase = 0;
L10:
  lacn2(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
  if (kase != 0) {
    if (kase == kase1) {

      /*           Multiply by inv(L). */

      if (lnoti) {
        i__1 = *n - 1;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MIN */
          i__2 = *kl, i__3 = *n - j;
          lm = min(i__2,i__3);
          jp = ipiv[j];
          t = work[jp];
          if (jp != j) {
            work[jp] = work[j];
            work[j] = t;
          }
          d__1 = -t;
          axpy(&lm, &d__1, &ab[kd + 1 + j * ab_dim1], &c__1, &
            work[j + 1], &c__1);
          /* L20: */
        }
      }

      /*           Multiply by inv(U). */

      i__1 = *kl + *ku;
      latbs("Upper", "No transpose", "Non-unit", normin, n, &i__1, &
        ab[ab_offset], ldab, &work[1], &scale, &work[(*n << 1) +
        1], info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
    } else {

      /*           Multiply by inv(U**T). */

      i__1 = *kl + *ku;
      latbs("Upper", "Transpose", "Non-unit", normin, n, &i__1, &ab[
        ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 1],
          info, (ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);

      /*           Multiply by inv(L**T). */

      if (lnoti) {
        for (j = *n - 1; j >= 1; --j) {
          /* Computing MIN */
          i__1 = *kl, i__2 = *n - j;
          lm = min(i__1,i__2);
          work[j] -= dot(&lm, &ab[kd + 1 + j * ab_dim1], &c__1, &
            work[j + 1], &c__1);
          jp = ipiv[j];
          if (jp != j) {
            t = work[jp];
            work[jp] = work[j];
            work[j] = t;
          }
          /* L30: */
        }
      }
    }

    /*        Divide X by 1/SCALE if doing so will not cause overflow. */

    *(unsigned char *)normin = 'Y';
    if (scale != 1.) {
      ix = iamax(n, &work[1], &c__1);
      if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.)
      {
        goto L40;
      }
      rscl(n, &scale, &work[1], &c__1);
    }
    goto L10;
  }

  /*     Compute the estimate of the reciprocal condition number. */

  if (ainvnm != 0.) {
    *rcond = 1. / ainvnm / *anorm;
  }

L40:
  return 0;

  /*     End of DGBCON */

} /* dgbcon_ */

//-------------------------------------------------------------------------------------------------

// from gbequ - LAPACK computational routine (version 3.7.0)
template<class T>
int gbequ(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *r__, T *c__,
  T *rowcnd, T *colcnd, T *amax, integer *info)
{
  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
  T d__1, d__2, d__3;

  /* Local variables */
  static integer i__, j, kd;
  static T rcmin, rcmax;
  //extern doublereal dlamch_(char *, ftnlen);
  //extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
  static T bignum, smlnum;

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --r__;
  --c__;

  /* Function Body */
  *info = 0;
  if (*m < 0) {
    *info = -1;
  } else if (*n < 0) {
    *info = -2;
  } else if (*kl < 0) {
    *info = -3;
  } else if (*ku < 0) {
    *info = -4;
  } else if (*ldab < *kl + *ku + 1) {
    *info = -6;
  }
  if (*info != 0) {
    i__1 = -(*info);
    xerbla("DGBEQU", &i__1, (ftnlen)6);
    return 0;
  }

  /*     Quick return if possible */

  if (*m == 0 || *n == 0) {
    *rowcnd = 1.;
    *colcnd = 1.;
    *amax = 0.;
    return 0;
  }

  /*     Get machine constants. */

  smlnum = lamch("S", (ftnlen)1);
  bignum = 1. / smlnum;

  /*     Compute row scale factors. */

  i__1 = *m;
  for (i__ = 1; i__ <= i__1; ++i__) {
    r__[i__] = 0.;
    /* L10: */
  }

  /*     Find the maximum element in each row. */

  kd = *ku + 1;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    /* Computing MAX */
    i__2 = j - *ku;
    /* Computing MIN */
    i__4 = j + *kl;
    i__3 = min(i__4,*m);
    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
      /* Computing MAX */
      d__2 = r__[i__], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1],
        abs(d__1));
      r__[i__] = max(d__2,d__3);
      /* L20: */
    }
    /* L30: */
  }

  /*     Find the maximum and minimum scale factors. */

  rcmin = bignum;
  rcmax = 0.;
  i__1 = *m;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /* Computing MAX */
    d__1 = rcmax, d__2 = r__[i__];
    rcmax = max(d__1,d__2);
    /* Computing MIN */
    d__1 = rcmin, d__2 = r__[i__];
    rcmin = min(d__1,d__2);
    /* L40: */
  }
  *amax = rcmax;

  if (rcmin == 0.) {

    /*        Find the first zero scale factor and return an error code. */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
      if (r__[i__] == 0.) {
        *info = i__;
        return 0;
      }
      /* L50: */
    }
  } else {

    /*        Invert the scale factors. */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
      /* Computing MIN */
      /* Computing MAX */
      d__2 = r__[i__];
      d__1 = max(d__2,smlnum);
      r__[i__] = 1. / min(d__1,bignum);
      /* L60: */
    }

    /*        Compute ROWCND = min(R(I)) / max(R(I)) */

    *rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
  }

  /*     Compute column scale factors */

  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    c__[j] = 0.;
    /* L70: */
  }

  /*     Find the maximum element in each column, */
  /*     assuming the row scaling computed above. */

  kd = *ku + 1;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    /* Computing MAX */
    i__3 = j - *ku;
    /* Computing MIN */
    i__4 = j + *kl;
    i__2 = min(i__4,*m);
    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
      /* Computing MAX */
      d__2 = c__[j], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
        d__1)) * r__[i__];
      c__[j] = max(d__2,d__3);
      /* L80: */
    }
    /* L90: */
  }

  /*     Find the maximum and minimum scale factors. */

  rcmin = bignum;
  rcmax = 0.;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    /* Computing MIN */
    d__1 = rcmin, d__2 = c__[j];
    rcmin = min(d__1,d__2);
    /* Computing MAX */
    d__1 = rcmax, d__2 = c__[j];
    rcmax = max(d__1,d__2);
    /* L100: */
  }

  if (rcmin == 0.) {

    /*        Find the first zero scale factor and return an error code. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      if (c__[j] == 0.) {
        *info = *m + j;
        return 0;
      }
      /* L110: */
    }
  } else {

    /*        Invert the scale factors. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      /* Computing MIN */
      /* Computing MAX */
      d__2 = c__[j];
      d__1 = max(d__2,smlnum);
      c__[j] = 1. / min(d__1,bignum);
      /* L120: */
    }

    /*        Compute COLCND = min(C(J)) / max(C(J)) */

    *colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
  }

  return 0;

  /*     End of DGBEQU */

} /* dgbequ_ */

//-------------------------------------------------------------------------------------------------

// dgbequb_ -- LAPACK computational routine (version 3.7.0) -
template<class T>
int gbequb(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *r__, T *c__,
  T *rowcnd, T *colcnd, T *amax, integer *info)
{
  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
  T d__1, d__2, d__3;

  /* Builtin functions */
  //double log(doublereal), pow_di(doublereal *, integer *);

  /* Local variables */
  static integer i__, j, kd;
  static T radix, rcmin, rcmax;
  //extern doublereal dlamch_(char *, ftnlen);
  //extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
  static T bignum, logrdx, smlnum;

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --r__;
  --c__;

  /* Function Body */
  *info = 0;
  if (*m < 0) {
    *info = -1;
  } else if (*n < 0) {
    *info = -2;
  } else if (*kl < 0) {
    *info = -3;
  } else if (*ku < 0) {
    *info = -4;
  } else if (*ldab < *kl + *ku + 1) {
    *info = -6;
  }
  if (*info != 0) {
    i__1 = -(*info);
    xerbla("DGBEQUB", &i__1, (ftnlen)7);
    return 0;
  }

  /*     Quick return if possible. */

  if (*m == 0 || *n == 0) {
    *rowcnd = 1.;
    *colcnd = 1.;
    *amax = 0.;
    return 0;
  }

  /*     Get machine constants.  Assume SMLNUM is a power of the radix. */

  smlnum = lamch("S", (ftnlen)1);
  bignum = 1. / smlnum;
  radix = lamch("B", (ftnlen)1);
  logrdx = log(radix);

  /*     Compute row scale factors. */

  i__1 = *m;
  for (i__ = 1; i__ <= i__1; ++i__) {
    r__[i__] = 0.;
    /* L10: */
  }

  /*     Find the maximum element in each row. */

  kd = *ku + 1;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    /* Computing MAX */
    i__2 = j - *ku;
    /* Computing MIN */
    i__4 = j + *kl;
    i__3 = min(i__4,*m);
    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
      /* Computing MAX */
      d__2 = r__[i__], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1],
        abs(d__1));
      r__[i__] = max(d__2,d__3);
      /* L20: */
    }
    /* L30: */
  }
  i__1 = *m;
  for (i__ = 1; i__ <= i__1; ++i__) {
    if (r__[i__] > 0.) {
      i__3 = (integer) (log(r__[i__]) / logrdx);
      r__[i__] = pow_di(&radix, &i__3);
    }
  }

  /*     Find the maximum and minimum scale factors. */

  rcmin = bignum;
  rcmax = 0.;
  i__1 = *m;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /* Computing MAX */
    d__1 = rcmax, d__2 = r__[i__];
    rcmax = max(d__1,d__2);
    /* Computing MIN */
    d__1 = rcmin, d__2 = r__[i__];
    rcmin = min(d__1,d__2);
    /* L40: */
  }
  *amax = rcmax;

  if (rcmin == 0.) {

    /*        Find the first zero scale factor and return an error code. */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
      if (r__[i__] == 0.) {
        *info = i__;
        return 0;
      }
      /* L50: */
    }
  } else {

    /*        Invert the scale factors. */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
      /* Computing MIN */
      /* Computing MAX */
      d__2 = r__[i__];
      d__1 = max(d__2,smlnum);
      r__[i__] = 1. / min(d__1,bignum);
      /* L60: */
    }

    /*        Compute ROWCND = min(R(I)) / max(R(I)). */

    *rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
  }

  /*     Compute column scale factors. */

  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    c__[j] = 0.;
    /* L70: */
  }

  /*     Find the maximum element in each column, */
  /*     assuming the row scaling computed above. */

  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    /* Computing MAX */
    i__3 = j - *ku;
    /* Computing MIN */
    i__4 = j + *kl;
    i__2 = min(i__4,*m);
    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
      /* Computing MAX */
      d__2 = c__[j], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
        d__1)) * r__[i__];
      c__[j] = max(d__2,d__3);
      /* L80: */
    }
    if (c__[j] > 0.) {
      i__2 = (integer) (log(c__[j]) / logrdx);
      c__[j] = pow_di(&radix, &i__2);
    }
    /* L90: */
  }

  /*     Find the maximum and minimum scale factors. */

  rcmin = bignum;
  rcmax = 0.;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    /* Computing MIN */
    d__1 = rcmin, d__2 = c__[j];
    rcmin = min(d__1,d__2);
    /* Computing MAX */
    d__1 = rcmax, d__2 = c__[j];
    rcmax = max(d__1,d__2);
    /* L100: */
  }

  if (rcmin == 0.) {

    /*        Find the first zero scale factor and return an error code. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      if (c__[j] == 0.) {
        *info = *m + j;
        return 0;
      }
      /* L110: */
    }
  } else {

    /*        Invert the scale factors. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      /* Computing MIN */
      /* Computing MAX */
      d__2 = c__[j];
      d__1 = max(d__2,smlnum);
      c__[j] = 1. / min(d__1,bignum);
      /* L120: */
    }

    /*        Compute COLCND = min(C(J)) / max(C(J)). */

    *colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
  }

  return 0;

  /*     End of DGBEQUB */

} /* dgbequb_ */


//-------------------------------------------------------------------------------------------------

// from dgbrfs - LAPACK computational routine (version 3.7.0)
template<class T>
int gbrfs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab, integer *ldab,
  T *afb, integer *ldafb, integer *ipiv, T *b, integer *ldb, T *x, integer *ldx, T *ferr, T *berr,
  T *work, integer *iwork, integer *info, ftnlen trans_len)
{
  /* Table of constant values */
  static integer c__1 = 1;
  static doublereal c_b15 = -1.;
  static doublereal c_b17 = 1.;


  /* System generated locals */
  integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset,
    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
  T d__1, d__2, d__3;

  /* Local variables */
  static integer i__, j, k;
  static T s;
  static integer kk;
  static T xk;
  static integer nz;
  static T eps;
  static integer kase;
  static T safe1, safe2;
  //extern /* Subroutine */ int dgbmv_(char *, integer *, integer *, integer *
  //  , integer *, doublereal *, doublereal *, integer *, doublereal *,
  //  integer *, doublereal *, doublereal *, integer *, ftnlen);
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static integer isave[3];
  //extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *,
  //  doublereal *, integer *), daxpy_(integer *, doublereal *,
  //    doublereal *, integer *, doublereal *, integer *);
  static integer count;
  //extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
  //  integer *, doublereal *, integer *, integer *);
  //extern doublereal dlamch_(char *, ftnlen);
  static T safmin;
  //extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dgbtrs_(
  //  char *, integer *, integer *, integer *, integer *, doublereal *,
  //  integer *, integer *, doublereal *, integer *, integer *, ftnlen);
  static logical notran;
  static char transt[1];
  static T lstres;

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  afb_dim1 = *ldafb;
  afb_offset = 1 + afb_dim1;
  afb -= afb_offset;
  --ipiv;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;
  x_dim1 = *ldx;
  x_offset = 1 + x_dim1;
  x -= x_offset;
  --ferr;
  --berr;
  --work;
  --iwork;

  /* Function Body */
  *info = 0;
  notran = lsame(trans, "N", (ftnlen)1, (ftnlen)1);
  if (! notran && ! lsame(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame(
    trans, "C", (ftnlen)1, (ftnlen)1)) {
    *info = -1;
  } else if (*n < 0) {
    *info = -2;
  } else if (*kl < 0) {
    *info = -3;
  } else if (*ku < 0) {
    *info = -4;
  } else if (*nrhs < 0) {
    *info = -5;
  } else if (*ldab < *kl + *ku + 1) {
    *info = -7;
  } else if (*ldafb < (*kl << 1) + *ku + 1) {
    *info = -9;
  } else if (*ldb < max(1,*n)) {
    *info = -12;
  } else if (*ldx < max(1,*n)) {
    *info = -14;
  }
  if (*info != 0) {
    i__1 = -(*info);
    xerbla("DGBRFS", &i__1, (ftnlen)6);
    return 0;
  }

  /*     Quick return if possible */

  if (*n == 0 || *nrhs == 0) {
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
      ferr[j] = 0.;
      berr[j] = 0.;
      /* L10: */
    }
    return 0;
  }

  if (notran) {
    *(unsigned char *)transt = 'T';
  } else {
    *(unsigned char *)transt = 'N';
  }

  /*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

  /* Computing MIN */
  i__1 = *kl + *ku + 2, i__2 = *n + 1;
  nz = min(i__1,i__2);
  eps = lamch("Epsilon", (ftnlen)7);
  safmin = lamch("Safe minimum", (ftnlen)12);
  safe1 = nz * safmin;
  safe2 = safe1 / eps;

  /*     Do for each right hand side */

  i__1 = *nrhs;
  for (j = 1; j <= i__1; ++j) {

    count = 1;
    lstres = 3.;
  L20:

    /*        Loop until stopping criterion is satisfied. */

    /*        Compute residual R = B - op(A) * X, */
    /*        where op(A) = A, A**T, or A**H, depending on TRANS. */

    copy(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
    gbmv(trans, n, n, kl, ku, &c_b15, &ab[ab_offset], ldab, &x[j *
      x_dim1 + 1], &c__1, &c_b17, &work[*n + 1], &c__1, (ftnlen)1);

    /*        Compute componentwise relative backward error from formula */

    /*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

    /*        where abs(Z) is the componentwise absolute value of the matrix */
    /*        or vector Z.  If the i-th component of the denominator is less */
    /*        than SAFE2, then SAFE1 is added to the i-th components of the */
    /*        numerator and denominator before dividing. */

    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
      /* L30: */
    }

    /*        Compute abs(op(A))*abs(X) + abs(B). */

    if (notran) {
      i__2 = *n;
      for (k = 1; k <= i__2; ++k) {
        kk = *ku + 1 - k;
        xk = (d__1 = x[k + j * x_dim1], abs(d__1));
        /* Computing MAX */
        i__3 = 1, i__4 = k - *ku;
        /* Computing MIN */
        i__6 = *n, i__7 = k + *kl;
        i__5 = min(i__6,i__7);
        for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
          work[i__] += (d__1 = ab[kk + i__ + k * ab_dim1], abs(d__1)
            ) * xk;
          /* L40: */
        }
        /* L50: */
      }
    } else {
      i__2 = *n;
      for (k = 1; k <= i__2; ++k) {
        s = 0.;
        kk = *ku + 1 - k;
        /* Computing MAX */
        i__5 = 1, i__3 = k - *ku;
        /* Computing MIN */
        i__6 = *n, i__7 = k + *kl;
        i__4 = min(i__6,i__7);
        for (i__ = max(i__5,i__3); i__ <= i__4; ++i__) {
          s += (d__1 = ab[kk + i__ + k * ab_dim1], abs(d__1)) * (
            d__2 = x[i__ + j * x_dim1], abs(d__2));
          /* L60: */
        }
        work[k] += s;
        /* L70: */
      }
    }
    s = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      if (work[i__] > safe2) {
        /* Computing MAX */
        d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
          i__];
        s = max(d__2,d__3);
      } else {
        /* Computing MAX */
        d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1)
          / (work[i__] + safe1);
        s = max(d__2,d__3);
      }
      /* L80: */
    }
    berr[j] = s;

    /*        Test stopping criterion. Continue iterating if */
    /*           1) The residual BERR(J) is larger than machine epsilon, and */
    /*           2) BERR(J) decreased by at least a factor of 2 during the */
    /*              last iteration, and */
    /*           3) At most ITMAX iterations tried. */

    if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

      /*           Update solution and try again. */

      gbtrs(trans, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1],
        &work[*n + 1], n, info, (ftnlen)1);
      axpy(n, &c_b17, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
      lstres = berr[j];
      ++count;
      goto L20;
    }

    /*        Bound error from formula */

    /*        norm(X - XTRUE) / norm(X) .le. FERR = */
    /*        norm( abs(inv(op(A)))* */
    /*           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X) */

    /*        where */
    /*          norm(Z) is the magnitude of the largest component of Z */
    /*          inv(op(A)) is the inverse of op(A) */
    /*          abs(Z) is the componentwise absolute value of the matrix or */
    /*             vector Z */
    /*          NZ is the maximum number of nonzeros in any row of A, plus 1 */
    /*          EPS is machine epsilon */

    /*        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B)) */
    /*        is incremented by SAFE1 if the i-th component of */
    /*        abs(op(A))*abs(X) + abs(B) is less than SAFE2. */

    /*        Use DLACN2 to estimate the infinity-norm of the matrix */
    /*           inv(op(A)) * diag(W), */
    /*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      if (work[i__] > safe2) {
        work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps *
          work[i__];
      } else {
        work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps *
          work[i__] + safe1;
      }
      /* L90: */
    }

    kase = 0;
  L100:
    lacn2(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &kase, isave);
    if (kase != 0) {
      if (kase == 1) {

        /*              Multiply by diag(W)*inv(op(A)**T). */

        gbtrs(transt, n, kl, ku, &c__1, &afb[afb_offset], ldafb,
          &ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          work[*n + i__] *= work[i__];
          /* L110: */
        }
      } else {

        /*              Multiply by inv(op(A))*diag(W). */

        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          work[*n + i__] *= work[i__];
          /* L120: */
        }
        gbtrs(trans, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &
          ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
      }
      goto L100;
    }

    /*        Normalize error. */

    lstres = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      /* Computing MAX */
      d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
      lstres = max(d__2,d__3);
      /* L130: */
    }
    if (lstres != 0.) {
      ferr[j] /= lstres;
    }

    /* L140: */
  }

  return 0;

  /*     End of DGBRFS */

} /* dgbrfs_ */

//-------------------------------------------------------------------------------------------------

// from dgbrfsx,  LAPACK computational routine (version 3.7.0) --
template<class T>
int gbrfsx(char *trans, char *equed, integer *n, integer *kl, integer *ku, integer *nrhs,
  T *ab, integer *ldab, T *afb, integer *ldafb, integer *ipiv, T *r__, T *c__, T *b, integer *ldb,
  T *x, integer *ldx, T *rcond, T *berr, integer *n_err_bnds__, T *err_bnds_norm__,
  T *err_bnds_comp__, integer *nparams, T *params, T *work, integer *iwork, integer *info,
  ftnlen trans_len, ftnlen equed_len)
{
  /* Table of constant values */
  static integer c_n1 = -1;
  static integer c__0 = 0;
  static integer c__1 = 1;


  /* System generated locals */
  integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset,
    x_dim1, x_offset, err_bnds_norm_dim1, err_bnds_norm_offset,
    err_bnds_comp_dim1, err_bnds_comp_offset, i__1;
  T d__1, d__2;

  /* Builtin functions */
  double sqrt(T);

  /* Local variables */
  static T illrcond_thresh__, unstable_thresh__, err_lbnd__;
  static integer ref_type__;
  //extern integer ilatrans_(char *, ftnlen);
  static integer j;
  static T rcond_tmp__;
  static integer prec_type__, trans_type__;
  //extern doublereal dla_gbrcond__(char *, integer *, integer *, integer *,
  //  doublereal *, integer *, doublereal *, integer *, integer *,
  //  integer *, doublereal *, integer *, doublereal *, integer *,
  //  ftnlen);
  static T cwise_wrong__;
  //extern /* Subroutine */ int dla_gbrfsx_extended__(integer *, integer *,
  //  integer *, integer *, integer *, integer *, doublereal *, integer
  //  *, doublereal *, integer *, integer *, logical *, doublereal *,
  //  doublereal *, integer *, doublereal *, integer *, doublereal *,
  //  integer *, doublereal *, doublereal *, doublereal *, doublereal *,
  //  doublereal *, doublereal *, doublereal *, integer *, doublereal *
  //  , doublereal *, logical *, integer *);
  static char norm[1];
  static logical ignore_cwise__;
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static T anorm;
  //extern doublereal dlangb_(char *, integer *, integer *, integer *,
  //  doublereal *, integer *, doublereal *, ftnlen), dlamch_(char *,
  //    ftnlen);
  //extern /* Subroutine */ int dgbcon_(char *, integer *, integer *, integer
  //  *, doublereal *, integer *, integer *, doublereal *, doublereal *,
  //  doublereal *, integer *, integer *, ftnlen), xerbla_(char *,
  //    integer *, ftnlen);
  static logical colequ, notran, rowequ;
  //extern integer ilaprec_(char *, ftnlen);
  static integer ithresh, n_norms__;
  static T rthresh;

  /* Parameter adjustments */
  err_bnds_comp_dim1 = *nrhs;
  err_bnds_comp_offset = 1 + err_bnds_comp_dim1;
  err_bnds_comp__ -= err_bnds_comp_offset;
  err_bnds_norm_dim1 = *nrhs;
  err_bnds_norm_offset = 1 + err_bnds_norm_dim1;
  err_bnds_norm__ -= err_bnds_norm_offset;
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  afb_dim1 = *ldafb;
  afb_offset = 1 + afb_dim1;
  afb -= afb_offset;
  --ipiv;
  --r__;
  --c__;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;
  x_dim1 = *ldx;
  x_offset = 1 + x_dim1;
  x -= x_offset;
  --berr;
  --params;
  --work;
  --iwork;

  /* Function Body */
  *info = 0;
  trans_type__ = ilatrans(trans, (ftnlen)1);
  ref_type__ = 1;
  if (*nparams >= 1) {
    if (params[1] < 0.) {
      params[1] = 1.;
    } else {
      ref_type__ = (integer) params[1];
    }
  }

  /*     Set default parameters. */

  illrcond_thresh__ = (doublereal) (*n) * lamch("Epsilon", (ftnlen)7);
  ithresh = 10;
  rthresh = .5;
  unstable_thresh__ = .25;
  ignore_cwise__ = FALSE_;

  if (*nparams >= 2) {
    if (params[2] < 0.) {
      params[2] = (doublereal) ithresh;
    } else {
      ithresh = (integer) params[2];
    }
  }
  if (*nparams >= 3) {
    if (params[3] < 0.) {
      if (ignore_cwise__) {
        params[3] = 0.;
      } else {
        params[3] = 1.;
      }
    } else {
      ignore_cwise__ = params[3] == 0.;
    }
  }
  if (ref_type__ == 0 || *n_err_bnds__ == 0) {
    n_norms__ = 0;
  } else if (ignore_cwise__) {
    n_norms__ = 1;
  } else {
    n_norms__ = 2;
  }

  notran = lsame(trans, "N", (ftnlen)1, (ftnlen)1);
  rowequ = lsame(equed, "R", (ftnlen)1, (ftnlen)1) || lsame(equed, "B", (
    ftnlen)1, (ftnlen)1);
  colequ = lsame(equed, "C", (ftnlen)1, (ftnlen)1) || lsame(equed, "B", (
    ftnlen)1, (ftnlen)1);

  /*     Test input parameters. */

  if (trans_type__ == -1) {
    *info = -1;
  } else if (! rowequ && ! colequ && ! lsame(equed, "N", (ftnlen)1, (
    ftnlen)1)) {
    *info = -2;
  } else if (*n < 0) {
    *info = -3;
  } else if (*kl < 0) {
    *info = -4;
  } else if (*ku < 0) {
    *info = -5;
  } else if (*nrhs < 0) {
    *info = -6;
  } else if (*ldab < *kl + *ku + 1) {
    *info = -8;
  } else if (*ldafb < (*kl << 1) + *ku + 1) {
    *info = -10;
  } else if (*ldb < max(1,*n)) {
    *info = -13;
  } else if (*ldx < max(1,*n)) {
    *info = -15;
  }
  if (*info != 0) {
    i__1 = -(*info);
    xerbla("DGBRFSX", &i__1, (ftnlen)7);
    return 0;
  }

  /*     Quick return if possible. */

  if (*n == 0 || *nrhs == 0) {
    *rcond = 1.;
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
      berr[j] = 0.;
      if (*n_err_bnds__ >= 1) {
        err_bnds_norm__[j + err_bnds_norm_dim1] = 1.;
        err_bnds_comp__[j + err_bnds_comp_dim1] = 1.;
      }
      if (*n_err_bnds__ >= 2) {
        err_bnds_norm__[j + (err_bnds_norm_dim1 << 1)] = 0.;
        err_bnds_comp__[j + (err_bnds_comp_dim1 << 1)] = 0.;
      }
      if (*n_err_bnds__ >= 3) {
        err_bnds_norm__[j + err_bnds_norm_dim1 * 3] = 1.;
        err_bnds_comp__[j + err_bnds_comp_dim1 * 3] = 1.;
      }
    }
    return 0;
  }

  /*     Default to failure. */

  *rcond = 0.;
  i__1 = *nrhs;
  for (j = 1; j <= i__1; ++j) {
    berr[j] = 1.;
    if (*n_err_bnds__ >= 1) {
      err_bnds_norm__[j + err_bnds_norm_dim1] = 1.;
      err_bnds_comp__[j + err_bnds_comp_dim1] = 1.;
    }
    if (*n_err_bnds__ >= 2) {
      err_bnds_norm__[j + (err_bnds_norm_dim1 << 1)] = 1.;
      err_bnds_comp__[j + (err_bnds_comp_dim1 << 1)] = 1.;
    }
    if (*n_err_bnds__ >= 3) {
      err_bnds_norm__[j + err_bnds_norm_dim1 * 3] = 0.;
      err_bnds_comp__[j + err_bnds_comp_dim1 * 3] = 0.;
    }
  }

  /*     Compute the norm of A and the reciprocal of the condition */
  /*     number of A. */

  if (notran) {
    *(unsigned char *)norm = 'I';
  } else {
    *(unsigned char *)norm = '1';
  }
  anorm = langb(norm, n, kl, ku, &ab[ab_offset], ldab, &work[1], (ftnlen)
    1);
  gbcon(norm, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], &anorm, rcond,
    &work[1], &iwork[1], info, (ftnlen)1);

  /*     Perform refinement on each right-hand side */

  if (ref_type__ != 0 && *info == 0) {
    prec_type__ = ilaprec("E", (ftnlen)1);
    if (notran) {
      la_gbrfsx_extended(&prec_type__, &trans_type__, n, kl, ku,
        nrhs, &ab[ab_offset], ldab, &afb[afb_offset], ldafb, &
        ipiv[1], &colequ, &c__[1], &b[b_offset], ldb, &x[x_offset]
        , ldx, &berr[1], &n_norms__, &err_bnds_norm__[
          err_bnds_norm_offset], &err_bnds_comp__[
            err_bnds_comp_offset], &work[*n + 1], &work[1], &work[(*n
              << 1) + 1], &work[1], rcond, &ithresh, &rthresh, &
              unstable_thresh__, &ignore_cwise__, info);
    } else {
      la_gbrfsx_extended(&prec_type__, &trans_type__, n, kl, ku,
        nrhs, &ab[ab_offset], ldab, &afb[afb_offset], ldafb, &
        ipiv[1], &rowequ, &r__[1], &b[b_offset], ldb, &x[x_offset]
        , ldx, &berr[1], &n_norms__, &err_bnds_norm__[
          err_bnds_norm_offset], &err_bnds_comp__[
            err_bnds_comp_offset], &work[*n + 1], &work[1], &work[(*n
              << 1) + 1], &work[1], rcond, &ithresh, &rthresh, &
              unstable_thresh__, &ignore_cwise__, info);
    }
  }
  /* Computing MAX */
  d__1 = 10., d__2 = sqrt((doublereal) (*n));
  err_lbnd__ = max(d__1,d__2) * lamch("Epsilon", (ftnlen)7);
  if (*n_err_bnds__ >= 1 && n_norms__ >= 1) {

    /*     Compute scaled normwise condition number cond(A*C). */

    if (colequ && notran) {
      rcond_tmp__ = la_gbrcond(trans, n, kl, ku, &ab[ab_offset],
        ldab, &afb[afb_offset], ldafb, &ipiv[1], &c_n1, &c__[1],
        info, &work[1], &iwork[1], (ftnlen)1);
    } else if (rowequ && ! notran) {
      rcond_tmp__ = la_gbrcond(trans, n, kl, ku, &ab[ab_offset],
        ldab, &afb[afb_offset], ldafb, &ipiv[1], &c_n1, &r__[1],
        info, &work[1], &iwork[1], (ftnlen)1);
    } else {
      rcond_tmp__ = la_gbrcond(trans, n, kl, ku, &ab[ab_offset],
        ldab, &afb[afb_offset], ldafb, &ipiv[1], &c__0, &r__[1],
        info, &work[1], &iwork[1], (ftnlen)1);
    }
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {

      /*     Cap the error at 1.0. */

      if (*n_err_bnds__ >= 2 && err_bnds_norm__[j + (err_bnds_norm_dim1
        << 1)] > 1.) {
        err_bnds_norm__[j + (err_bnds_norm_dim1 << 1)] = 1.;
      }

      /*     Threshold the error (see LAWN). */

      if (rcond_tmp__ < illrcond_thresh__) {
        err_bnds_norm__[j + (err_bnds_norm_dim1 << 1)] = 1.;
        err_bnds_norm__[j + err_bnds_norm_dim1] = 0.;
        if (*info <= *n) {
          *info = *n + j;
        }
      } else if (err_bnds_norm__[j + (err_bnds_norm_dim1 << 1)] <
        err_lbnd__) {
        err_bnds_norm__[j + (err_bnds_norm_dim1 << 1)] = err_lbnd__;
        err_bnds_norm__[j + err_bnds_norm_dim1] = 1.;
      }

      /*     Save the condition number. */

      if (*n_err_bnds__ >= 3) {
        err_bnds_norm__[j + err_bnds_norm_dim1 * 3] = rcond_tmp__;
      }
    }
  }
  if (*n_err_bnds__ >= 1 && n_norms__ >= 2) {

    /*     Compute componentwise condition number cond(A*diag(Y(:,J))) for */
    /*     each right-hand side using the current solution as an estimate of */
    /*     the true solution.  If the componentwise error estimate is too */
    /*     large, then the solution is a lousy estimate of truth and the */
    /*     estimated RCOND may be too optimistic.  To avoid misleading users, */
    /*     the inverse condition number is set to 0.0 when the estimated */
    /*     cwise error is at least CWISE_WRONG. */

    cwise_wrong__ = sqrt(lamch("Epsilon", (ftnlen)7));
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
      if (err_bnds_comp__[j + (err_bnds_comp_dim1 << 1)] <
        cwise_wrong__) {
        rcond_tmp__ = la_gbrcond(trans, n, kl, ku, &ab[ab_offset],
          ldab, &afb[afb_offset], ldafb, &ipiv[1], &c__1, &x[j *
          x_dim1 + 1], info, &work[1], &iwork[1], (ftnlen)1);
      } else {
        rcond_tmp__ = 0.;
      }

      /*     Cap the error at 1.0. */

      if (*n_err_bnds__ >= 2 && err_bnds_comp__[j + (err_bnds_comp_dim1
        << 1)] > 1.) {
        err_bnds_comp__[j + (err_bnds_comp_dim1 << 1)] = 1.;
      }

      /*     Threshold the error (see LAWN). */

      if (rcond_tmp__ < illrcond_thresh__) {
        err_bnds_comp__[j + (err_bnds_comp_dim1 << 1)] = 1.;
        err_bnds_comp__[j + err_bnds_comp_dim1] = 0.;
        if (params[3] == 1. && *info < *n + j) {
          *info = *n + j;
        }
      } else if (err_bnds_comp__[j + (err_bnds_comp_dim1 << 1)] <
        err_lbnd__) {
        err_bnds_comp__[j + (err_bnds_comp_dim1 << 1)] = err_lbnd__;
        err_bnds_comp__[j + err_bnds_comp_dim1] = 1.;
      }

      /*     Save the condition number. */

      if (*n_err_bnds__ >= 3) {
        err_bnds_comp__[j + err_bnds_comp_dim1 * 3] = rcond_tmp__;
      }
    }
  }

  return 0;

  /*     End of DGBRFSX */

} /* dgbrfsx_ */











//-------------------------------------------------------------------------------------------------

// translated from dgbtf2, LAPACK computational routine (version 3.7.0)
template<class T>
int gbtf2(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, integer *ipiv,
  integer *info)
{
  /* Table of constant values */
  static integer c__1 = 1;
  static doublereal c_b9 = -1.;

  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
  T d__1;

  /* Local variables */
  static integer i__, j, km, jp, ju, kv;
  //extern /* Subroutine */ int dger_(integer *, integer *, T *,
  //  T *, integer *, T *, integer *, T *,
  //  integer *), dscal_(integer *, T *, T *, integer
  //    *), dswap_(integer *, T *, integer *, T *,
  //      integer *);
  //extern integer idamax_(integer *, T *, integer *);
  //extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

  /*     KV is the number of superdiagonals in the factor U, allowing for */
  /*     fill-in. */

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --ipiv;

  /* Function Body */
  kv = *ku + *kl;

  /*     Test the input parameters. */

  *info = 0;
  if (*m < 0) {
    *info = -1;
  } else if (*n < 0) {
    *info = -2;
  } else if (*kl < 0) {
    *info = -3;
  } else if (*ku < 0) {
    *info = -4;
  } else if (*ldab < *kl + kv + 1) {
    *info = -6;
  }
  if (*info != 0) {
    i__1 = -(*info);
    xerbla("DGBTF2", &i__1, (ftnlen)6);
    return 0;
  }

  /*     Quick return if possible */

  if (*m == 0 || *n == 0) {
    return 0;
  }

  /*     Gaussian elimination with partial pivoting */

  /*     Set fill-in elements in columns KU+2 to KV to zero. */

  i__1 = min(kv,*n);
  for (j = *ku + 2; j <= i__1; ++j) {
    i__2 = *kl;
    for (i__ = kv - j + 2; i__ <= i__2; ++i__) {
      ab[i__ + j * ab_dim1] = 0.;
      /* L10: */
    }
    /* L20: */
  }

  /*     JU is the index of the last column affected by the current stage */
  /*     of the factorization. */

  ju = 1;

  i__1 = min(*m,*n);
  for (j = 1; j <= i__1; ++j) {

    /*        Set fill-in elements in column J+KV to zero. */

    if (j + kv <= *n) {
      i__2 = *kl;
      for (i__ = 1; i__ <= i__2; ++i__) {
        ab[i__ + (j + kv) * ab_dim1] = 0.;
        /* L30: */
      }
    }

    /*        Find pivot and test for singularity. KM is the number of */
    /*        subdiagonal elements in the current column. */

    /* Computing MIN */
    i__2 = *kl, i__3 = *m - j;
    km = min(i__2,i__3);
    i__2 = km + 1;
    jp = iamax(&i__2, &ab[kv + 1 + j * ab_dim1], &c__1);
    ipiv[j] = jp + j - 1;
    if (ab[kv + jp + j * ab_dim1] != 0.) {
      /* Computing MAX */
      /* Computing MIN */
      i__4 = j + *ku + jp - 1;
      i__2 = ju, i__3 = min(i__4,*n);
      ju = max(i__2,i__3);

      /*           Apply interchange to columns J to JU. */

      if (jp != 1) {
        i__2 = ju - j + 1;
        i__3 = *ldab - 1;
        i__4 = *ldab - 1;
        swap(&i__2, &ab[kv + jp + j * ab_dim1], &i__3, &ab[kv + 1 +
          j * ab_dim1], &i__4);
      }

      if (km > 0) {

        /*              Compute multipliers. */

        d__1 = 1. / ab[kv + 1 + j * ab_dim1];
        scal(&km, &d__1, &ab[kv + 2 + j * ab_dim1], &c__1);

        /*              Update trailing submatrix within the band. */

        if (ju > j) {
          i__2 = ju - j;
          i__3 = *ldab - 1;
          i__4 = *ldab - 1;
          ger(&km, &i__2, &c_b9, &ab[kv + 2 + j * ab_dim1], &c__1,
            &ab[kv + (j + 1) * ab_dim1], &i__3, &ab[kv + 1 +
            (j + 1) * ab_dim1], &i__4);
        }
      }
    } else {

      /*           If pivot is zero, set INFO to the index of the pivot */
      /*           unless a zero pivot has already been found. */

      if (*info == 0) {
        *info = j;
      }
    }
    /* L40: */
  }
  return 0;

  /*     End of DGBTF2 */

} /* gbtf2 */

//-------------------------------------------------------------------------------------------------

// translated from dgbtrf,LAPACK computational routine (version 3.7.0)
template<class T>
int gbtrf(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, integer *ipiv,
  integer *info)
{
  // Table of constant values
  static integer c__1 = 1;
  static integer c__65 = 65;
  static T c_b18 = -1.;
  static T c_b31 = 1.;

  // System generated locals
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
  T d__1;

  // Local variables
  static integer i__, j, i2, i3, j2, j3, k2, jb, nb, ii, jj, jm, ip, jp, km,
    ju, kv, nw;
  // we need to comment these declarations - otherwise, the linker tries to find those functions
  // in the LaPackCPP namespace (but they belong to the BlasCPP namespace)
  //extern int ger(integer *, integer *, T *,
  //  T *, integer *, T *, integer *, T *,
  //  integer *);
  static T temp;
  //extern int scal(integer *, T *, T *,
  //  integer *), gemm(char *, char *, integer *, integer *, integer *
  //    , T *, T *, integer *, T *, integer *,
  //    T *, T *, integer *, ftnlen, ftnlen), dcopy_(
  //      integer *, T *, integer *, T *, integer *),
  //  swap(integer *, T *, integer *, T *, integer *
  //    );
  static T work13[4160], work31[4160];
  //extern int trsm(char *, char *, char *, char *,
  //  integer *, integer *, T *, T *, integer *,
  //  T *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dgbtf2_(
  //    integer *, integer *, integer *, integer *, T *, integer
  //    *, integer *, integer *);
  //extern integer iamax(integer *, T *, integer *);
  //extern int xerbla_(char *, integer *, ftnlen);
  //extern integer ilaenv(integer *, char *, char *, integer *, integer *,
  //  integer *, integer *, ftnlen, ftnlen);
  //extern int laswp(integer *, T *, integer *,
  //  integer *, integer *, integer *, integer *);


  // KV is the number of superdiagonals in the factor U, allowing for fill-in

  // Parameter adjustments
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --ipiv;

  // Function Body
  kv = *ku + *kl;

  // Test the input parameters.
  *info = 0;
  if(*m < 0) {
    *info = -1;
  }
  else if(*n < 0) {
    *info = -2;
  }
  else if(*kl < 0) {
    *info = -3;
  }
  else if(*ku < 0) {
    *info = -4;
  }
  else if(*ldab < *kl + kv + 1) {
    *info = -6;
  }
  if(*info != 0) {
    i__1 = -(*info);
    xerbla("DGBTRF", &i__1, (ftnlen)6);
    return 0;
  }

  // Quick return if possible
  if(*m == 0 || *n == 0) {
    return 0;
  }

  // Determine the block size for this environment
  nb = ilaenv(&c__1, "DGBTRF", " ", m, n, kl, ku, (ftnlen)6, (ftnlen)1);
  // The block size must not exceed the limit set by the size of the
  // local arrays WORK13 and WORK31.

  nb = min(nb, 64);

  if(nb <= 1 || nb > *kl) {
    // Use unblocked code
    gbtf2(m, n, kl, ku, &ab[ab_offset], ldab, &ipiv[1], info);
  }
  else {

  // Use blocked code
  // Zero the superdiagonal elements of the work array WORK13
    i__1 = nb;
    for(j = 1; j <= i__1; ++j) {
      i__2 = j - 1;
      for(i__ = 1; i__ <= i__2; ++i__) {
        work13[i__ + j * 65 - 66] = 0.;
        // L10:
      }
      // L20:
    }

    // Zero the subdiagonal elements of the work array WORK31
    i__1 = nb;
    for(j = 1; j <= i__1; ++j) {
      i__2 = nb;
      for(i__ = j + 1; i__ <= i__2; ++i__) {
        work31[i__ + j * 65 - 66] = 0.;
        // L30:
      }
      // L40:
    }

    // Gaussian elimination with partial pivoting
    // Set fill-in elements in columns KU+2 to KV to zero
    i__1 = min(kv, *n);
    for(j = *ku + 2; j <= i__1; ++j) {
      i__2 = *kl;
      for(i__ = kv - j + 2; i__ <= i__2; ++i__) {
        ab[i__ + j * ab_dim1] = 0.;
        // L50:
      }
      // L60:
    }

    // JU is the index of the last column affected by the current
    // stage of the factorization
    ju = 1;

    i__1 = min(*m, *n);
    i__2 = nb;
    for(j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
      // Computing MIN
      i__3 = nb, i__4 = min(*m, *n) - j + 1;
      jb = min(i__3, i__4);

      // The active part of the matrix is partitioned
      // A11   A12   A13
      // A21   A22   A23
      // A31   A32   A33
      // Here A11, A21 and A31 denote the current block of JB columns
      // which is about to be factorized. The number of rows in the
      // partitioning are JB, I2, I3 respectively, and the numbers
      // of columns are JB, J2, J3. The superdiagonal elements of A13
      // and the subdiagonal elements of A31 lie outside the band.

      // Computing MIN
      i__3 = *kl - jb, i__4 = *m - j - jb + 1;
      i2 = min(i__3, i__4);
      // Computing MIN
      i__3 = jb, i__4 = *m - j - *kl + 1;
      i3 = min(i__3, i__4);

      // J2 and J3 are computed after JU has been updated.

      // Factorize the current block of JB columns
      i__3 = j + jb - 1;
      for(jj = j; jj <= i__3; ++jj) {

        // Set fill-in elements in column JJ+KV to zero

        if(jj + kv <= *n) {
          i__4 = *kl;
          for(i__ = 1; i__ <= i__4; ++i__) {
            ab[i__ + (jj + kv) * ab_dim1] = 0.;
            // L70:
          }
        }

        // Find pivot and test for singularity. KM is the number of
        // subdiagonal elements in the current column. */

        // Computing MIN
        i__4 = *kl, i__5 = *m - jj;
        km = min(i__4, i__5);
        i__4 = km + 1;
        jp = iamax(&i__4, &ab[kv + 1 + jj * ab_dim1], &c__1);
        ipiv[jj] = jp + jj - j;
        if(ab[kv + jp + jj * ab_dim1] != 0.) {
          // Computing MAX
          // Computing MIN
          i__6 = jj + *ku + jp - 1;
          i__4 = ju, i__5 = min(i__6, *n);
          ju = max(i__4, i__5);
          if(jp != 1) {

            // Apply interchange to columns J to J+JB-1
            if(jp + jj - 1 < j + *kl) {

              i__4 = *ldab - 1;
              i__5 = *ldab - 1;
              swap(&jb, &ab[kv + 1 + jj - j + j * ab_dim1], &
                i__4, &ab[kv + jp + jj - j + j * ab_dim1],
                &i__5);
            }
            else {

            // The interchange affects columns J to JJ-1 of A31
            // which are stored in the work array WORK31
              i__4 = jj - j;
              i__5 = *ldab - 1;
              swap(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1],
                &i__5, &work31[jp + jj - j - *kl - 1], &
                c__65);
              i__4 = j + jb - jj;
              i__5 = *ldab - 1;
              i__6 = *ldab - 1;
              swap(&i__4, &ab[kv + 1 + jj * ab_dim1], &i__5, &
                ab[kv + jp + jj * ab_dim1], &i__6);
            }
          }

          // Compute multipliers

          d__1 = 1. / ab[kv + 1 + jj * ab_dim1];
          scal(&km, &d__1, &ab[kv + 2 + jj * ab_dim1], &c__1);

          // Update trailing submatrix within the band and within
          // the current block. JM is the index of the last column
          // which needs to be updated.

          // Computing MIN
          i__4 = ju, i__5 = j + jb - 1;
          jm = min(i__4, i__5);
          if(jm > jj) {
            i__4 = jm - jj;
            i__5 = *ldab - 1;
            i__6 = *ldab - 1;
            ger(&km, &i__4, &c_b18, &ab[kv + 2 + jj * ab_dim1],
              &c__1, &ab[kv + (jj + 1) * ab_dim1], &i__5, &
              ab[kv + 1 + (jj + 1) * ab_dim1], &i__6);
          }
        }
        else {

       // If pivot is zero, set INFO to the index of the pivot
       // unless a zero pivot has already been found.

          if(*info == 0) {
            *info = jj;
          }
        }

        // Copy current column of A31 into the work array WORK31

        // Computing MIN
        i__4 = jj - j + 1;
        nw = min(i__4, i3);
        if(nw > 0) {
          copy(&nw, &ab[kv + *kl + 1 - jj + j + jj * ab_dim1], &
            c__1, &work31[(jj - j + 1) * 65 - 65], &c__1);
        }
        // L80:
      }
      if(j + jb <= *n) {

        // Apply the row interchanges to the other blocks.

        // Computing MIN
        i__3 = ju - j + 1;
        j2 = min(i__3, kv) - jb;
        // Computing MAX
        i__3 = 0, i__4 = ju - j - kv + 1;
        j3 = max(i__3, i__4);

        // Use DLASWP to apply the row interchanges to A12, A22, and
        // A32.
        i__3 = *ldab - 1;
        laswp(&j2, &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__3, &
          c__1, &jb, &ipiv[j], &c__1);

        // Adjust the pivot indices.
        i__3 = j + jb - 1;
        for(i__ = j; i__ <= i__3; ++i__) {
          ipiv[i__] = ipiv[i__] + j - 1;
          // L90:
        }

        // Apply the row interchanges to A13, A23, and A33
        // columnwise.
        k2 = j - 1 + jb + j2;
        i__3 = j3;
        for(i__ = 1; i__ <= i__3; ++i__) {
          jj = k2 + i__;
          i__4 = j + jb - 1;
          for(ii = j + i__ - 1; ii <= i__4; ++ii) {
            ip = ipiv[ii];
            if(ip != ii) {
              temp = ab[kv + 1 + ii - jj + jj * ab_dim1];
              ab[kv + 1 + ii - jj + jj * ab_dim1] = ab[kv + 1 +
                ip - jj + jj * ab_dim1];
              ab[kv + 1 + ip - jj + jj * ab_dim1] = temp;
            }
            // L100:
          }
          // L110:
        }

        // Update the relevant part of the trailing submatrix

        if(j2 > 0) {

          // Update A12
          i__3 = *ldab - 1;
          i__4 = *ldab - 1;
          trsm("Left", "Lower", "No transpose", "Unit", &jb, &j2,
            &c_b31, &ab[kv + 1 + j * ab_dim1], &i__3, &ab[kv
            + 1 - jb + (j + jb) * ab_dim1], &i__4, (ftnlen)4,
            (ftnlen)5, (ftnlen)12, (ftnlen)4);

          if(i2 > 0) {

            // Update A22
            i__3 = *ldab - 1;
            i__4 = *ldab - 1;
            i__5 = *ldab - 1;
            gemm("No transpose", "No transpose", &i2, &j2, &jb,
              &c_b18, &ab[kv + 1 + jb + j * ab_dim1], &i__3,
              &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__4,
              &c_b31, &ab[kv + 1 + (j + jb) * ab_dim1], &
              i__5, (ftnlen)12, (ftnlen)12);
          }

          if(i3 > 0) {

            // Update A32
            i__3 = *ldab - 1;
            i__4 = *ldab - 1;
            gemm("No transpose", "No transpose", &i3, &j2, &jb,
              &c_b18, work31, &c__65, &ab[kv + 1 - jb + (j
                + jb) * ab_dim1], &i__3, &c_b31, &ab[kv + *kl
              + 1 - jb + (j + jb) * ab_dim1], &i__4, (
                ftnlen)12, (ftnlen)12);
          }
        }

        if(j3 > 0) {

          // Copy the lower triangle of A13 into the work array
          // WORK13
          i__3 = j3;
          for(jj = 1; jj <= i__3; ++jj) {
            i__4 = jb;
            for(ii = jj; ii <= i__4; ++ii) {
              work13[ii + jj * 65 - 66] = ab[ii - jj + 1 + (jj
                + j + kv - 1) * ab_dim1];
              // L120:
            }
            // L130:
          }

          // Update A13 in the work array
          i__3 = *ldab - 1;
          trsm("Left", "Lower", "No transpose", "Unit", &jb, &j3,
            &c_b31, &ab[kv + 1 + j * ab_dim1], &i__3, work13,
            &c__65, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)
            4);

          if(i2 > 0) {
            // Update A23
            i__3 = *ldab - 1;
            i__4 = *ldab - 1;
            gemm("No transpose", "No transpose", &i2, &j3, &jb,
              &c_b18, &ab[kv + 1 + jb + j * ab_dim1], &i__3,
              work13, &c__65, &c_b31, &ab[jb + 1 + (j + kv)
              * ab_dim1], &i__4, (ftnlen)12, (ftnlen)12);
          }

          if(i3 > 0) {
            // Update A33
            i__3 = *ldab - 1;
            gemm("No transpose", "No transpose", &i3, &j3, &jb,
              &c_b18, work31, &c__65, work13, &c__65, &
              c_b31, &ab[*kl + 1 + (j + kv) * ab_dim1], &
              i__3, (ftnlen)12, (ftnlen)12);
          }

          // Copy the lower triangle of A13 back into place
          i__3 = j3;
          for(jj = 1; jj <= i__3; ++jj) {
            i__4 = jb;
            for(ii = jj; ii <= i__4; ++ii) {
              ab[ii - jj + 1 + (jj + j + kv - 1) * ab_dim1] =
                work13[ii + jj * 65 - 66];
              // L140:
            }
            // L150:
          }
        }
      }
      else {

      // Adjust the pivot indices.

        i__3 = j + jb - 1;
        for(i__ = j; i__ <= i__3; ++i__) {
          ipiv[i__] = ipiv[i__] + j - 1;
          // L160:
        }
      }

      // Partially undo the interchanges in the current block to
      // restore the upper triangular form of A31 and copy the upper
      // triangle of A31 back into place
      i__3 = j;
      for(jj = j + jb - 1; jj >= i__3; --jj) {
        jp = ipiv[jj] - jj + 1;
        if(jp != 1) {

          // Apply interchange to columns J to JJ-1

          if(jp + jj - 1 < j + *kl) {
            // The interchange does not affect A31
            i__4 = jj - j;
            i__5 = *ldab - 1;
            i__6 = *ldab - 1;
            swap(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
              i__5, &ab[kv + jp + jj - j + j * ab_dim1], &
              i__6);
          }
          else {
            // The interchange does affect A31
            i__4 = jj - j;
            i__5 = *ldab - 1;
            swap(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
              i__5, &work31[jp + jj - j - *kl - 1], &c__65);
          }
        }

        // Copy the current column of A31 back into place
        // Computing MIN
        i__4 = i3, i__5 = jj - j + 1;
        nw = min(i__4, i__5);
        if(nw > 0) {
          copy(&nw, &work31[(jj - j + 1) * 65 - 65], &c__1, &ab[
            kv + *kl + 1 - jj + j + jj * ab_dim1], &c__1);
        }
        // L170:
      }
      // L180:
    }
  }

  return 0;
} // gbtrf

//-------------------------------------------------------------------------------------------------

// translated from dgbtrs, LAPACK computational routine (version 3.7.0)
template<class T>
int gbtrs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab,
  integer *ldab, integer *ipiv, T *b, integer *ldb, integer *info, ftnlen trans_len)
{
  // Table of constant values
  static T c_b7 = -1.;
  static integer c__1 = 1;
  static T c_b23 = 1.;

  // System generated locals
  integer ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2, i__3;

  // Local variables
  static integer i__, j, l, kd, lm;
  //extern int dger_(integer *, integer *, doublereal *,
  //  doublereal *, integer *, doublereal *, integer *, doublereal *,
  //  integer *);
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  //extern int dgemv_(char *, integer *, integer *,
  //  doublereal *, doublereal *, integer *, doublereal *, integer *,
  //  doublereal *, doublereal *, integer *, ftnlen), dswap_(integer *,
  //    doublereal *, integer *, doublereal *, integer *), dtbsv_(char *,
  //      char *, char *, integer *, integer *, doublereal *, integer *,
  //      doublereal *, integer *, ftnlen, ftnlen, ftnlen);
  static logical lnoti;
  //extern int xerbla(char *, integer *, ftnlen);
  static logical notran;

  // Parameter adjustments
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --ipiv;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;

  // Function Body
  *info = 0;
  notran = lsame(trans, "N", (ftnlen)1, (ftnlen)1);
  if(!notran && !lsame(trans, "T", (ftnlen)1, (ftnlen)1) && !lsame(
    trans, "C", (ftnlen)1, (ftnlen)1)) {
    *info = -1;
  }
  else if(*n < 0) {
    *info = -2;
  }
  else if(*kl < 0) {
    *info = -3;
  }
  else if(*ku < 0) {
    *info = -4;
  }
  else if(*nrhs < 0) {
    *info = -5;
  }
  else if(*ldab < (*kl << 1) + *ku + 1) {
    *info = -7;
  }
  else if(*ldb < max(1, *n)) {
    *info = -10;
  }
  if(*info != 0) {
    i__1 = -(*info);
    xerbla("DGBTRS", &i__1, (ftnlen)6);
    return 0;
  }

  // Quick return if possible
  if(*n == 0 || *nrhs == 0) {
    return 0;
  }

  kd = *ku + *kl + 1;
  lnoti = *kl > 0;

  if(notran) {

    // Solve  A*X = B.
    // Solve L*X = B, overwriting B with X.
    // L is represented as a product of permutations and unit lower
    // triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
    // where each transformation L(i) is a rank-one modification of
    // the identity matrix.
    if(lnoti) {
      i__1 = *n - 1;
      for(j = 1; j <= i__1; ++j) {
        // Computing MIN
        i__2 = *kl, i__3 = *n - j;
        lm = min(i__2, i__3);
        l = ipiv[j];
        if(l != j) {
          swap(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
        }
        ger(&lm, nrhs, &c_b7, &ab[kd + 1 + j * ab_dim1], &c__1, &b[
          j + b_dim1], ldb, &b[j + 1 + b_dim1], ldb);
        // L10:
      }
    }

    i__1 = *nrhs;
    for(i__ = 1; i__ <= i__1; ++i__) {

      // Solve U*X = B, overwriting B with X.
      i__2 = *kl + *ku;
      tbsv("Upper", "No transpose", "Non-unit", n, &i__2, &ab[
        ab_offset], ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5,
          (ftnlen)12, (ftnlen)8);
      // L20:
    }

  }
  else {

 // Solve A**T*X = B.
    i__1 = *nrhs;
    for(i__ = 1; i__ <= i__1; ++i__) {

      // Solve U**T*X = B, overwriting B with X.
      i__2 = *kl + *ku;
      tbsv("Upper", "Transpose", "Non-unit", n, &i__2, &ab[ab_offset],
        ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)9,
        (ftnlen)8);
      // L30:
    }

    // Solve L**T * X = B, overwriting B with X.
    if(lnoti) {
      for(j = *n - 1; j >= 1; --j) {
        // Computing MIN
        i__1 = *kl, i__2 = *n - j;
        lm = min(i__1, i__2);
        gemv("Transpose", &lm, nrhs, &c_b7, &b[j + 1 + b_dim1], ldb,
          &ab[kd + 1 + j * ab_dim1], &c__1, &c_b23, &b[j +
          b_dim1], ldb, (ftnlen)9);
        l = ipiv[j];
        if(l != j) {
          swap(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
        }
        // L40:
      }
    }
  }
  return 0;

  // End of DGBTRS

} // gbtrs

//-------------------------------------------------------------------------------------------------

// dla_gbamv -- LAPACK computational routine (version 3.7.1)
template<class T>
int la_gbamv(integer *trans, integer *m, integer *n, integer *kl, integer *ku, T *alpha, T *ab,
  integer *ldab, T *x, integer *incx, T *beta, T *y, integer *incy)
{
  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
  T d__1;

  /* Builtin functions */
  //double d_sign(doublereal *, doublereal *);

  /* Local variables */
  //extern integer ilatrans_(char *, ftnlen);
  static integer i__, j;
  static logical symb_zero__;
  static integer kd, ke, iy, jx, kx, ky, info;
  static T temp;
  static integer lenx, leny;
  static T safe1;
  //extern doublereal dlamch_(char *, ftnlen);
  //extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --x;
  --y;

  /* Function Body */
  info = 0;
  if (! (*trans == ilatrans("N", (ftnlen)1) || *trans == ilatrans("T", (
    ftnlen)1) || *trans == ilatrans("C", (ftnlen)1))) {
    info = 1;
  } else if (*m < 0) {
    info = 2;
  } else if (*n < 0) {
    info = 3;
  } else if (*kl < 0 || *kl > *m - 1) {
    info = 4;
  } else if (*ku < 0 || *ku > *n - 1) {
    info = 5;
  } else if (*ldab < *kl + *ku + 1) {
    info = 6;
  } else if (*incx == 0) {
    info = 8;
  } else if (*incy == 0) {
    info = 11;
  }
  if (info != 0) {
    xerbla("DLA_GBAMV ", &info, (ftnlen)10);
    return 0;
  }

  /*     Quick return if possible. */

  if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
    return 0;
  }

  /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
  /*     up the start points in  X  and  Y. */

  if (*trans == ilatrans("N", (ftnlen)1)) {
    lenx = *n;
    leny = *m;
  } else {
    lenx = *m;
    leny = *n;
  }
  if (*incx > 0) {
    kx = 1;
  } else {
    kx = 1 - (lenx - 1) * *incx;
  }
  if (*incy > 0) {
    ky = 1;
  } else {
    ky = 1 - (leny - 1) * *incy;
  }

  /*     Set SAFE1 essentially to be the underflow threshold times the */
  /*     number of additions in each row. */

  safe1 = lamch("Safe minimum", (ftnlen)12);
  safe1 = (*n + 1) * safe1;

  /*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

  /*     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to */
  /*     the inexact flag.  Still doesn't help change the iteration order */
  /*     to per-column. */

  kd = *ku + 1;
  ke = *kl + 1;
  iy = ky;
  if (*incx == 1) {
    if (*trans == ilatrans("N", (ftnlen)1)) {
      i__1 = leny;
      for (i__ = 1; i__ <= i__1; ++i__) {
        if (*beta == 0.) {
          symb_zero__ = TRUE_;
          y[iy] = 0.;
        } else if (y[iy] == 0.) {
          symb_zero__ = TRUE_;
        } else {
          symb_zero__ = FALSE_;
          y[iy] = *beta * (d__1 = y[iy], abs(d__1));
        }
        if (*alpha != 0.) {
          /* Computing MAX */
          i__2 = i__ - *kl;
          /* Computing MIN */
          i__4 = i__ + *ku;
          i__3 = min(i__4,lenx);
          for (j = max(i__2,1); j <= i__3; ++j) {
            temp = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
              d__1));
            symb_zero__ = symb_zero__ && (x[j] == 0. || temp ==
              0.);
            y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
          }
        }
        if (! symb_zero__) {
          y[iy] += d_sign(&safe1, &y[iy]);
        }
        iy += *incy;
      }
    } else {
      i__1 = leny;
      for (i__ = 1; i__ <= i__1; ++i__) {
        if (*beta == 0.) {
          symb_zero__ = TRUE_;
          y[iy] = 0.;
        } else if (y[iy] == 0.) {
          symb_zero__ = TRUE_;
        } else {
          symb_zero__ = FALSE_;
          y[iy] = *beta * (d__1 = y[iy], abs(d__1));
        }
        if (*alpha != 0.) {
          /* Computing MAX */
          i__3 = i__ - *kl;
          /* Computing MIN */
          i__4 = i__ + *ku;
          i__2 = min(i__4,lenx);
          for (j = max(i__3,1); j <= i__2; ++j) {
            temp = (d__1 = ab[ke - i__ + j + i__ * ab_dim1], abs(
              d__1));
            symb_zero__ = symb_zero__ && (x[j] == 0. || temp ==
              0.);
            y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
          }
        }
        if (! symb_zero__) {
          y[iy] += d_sign(&safe1, &y[iy]);
        }
        iy += *incy;
      }
    }
  } else {
    if (*trans == ilatrans("N", (ftnlen)1)) {
      i__1 = leny;
      for (i__ = 1; i__ <= i__1; ++i__) {
        if (*beta == 0.) {
          symb_zero__ = TRUE_;
          y[iy] = 0.;
        } else if (y[iy] == 0.) {
          symb_zero__ = TRUE_;
        } else {
          symb_zero__ = FALSE_;
          y[iy] = *beta * (d__1 = y[iy], abs(d__1));
        }
        if (*alpha != 0.) {
          jx = kx;
          /* Computing MAX */
          i__2 = i__ - *kl;
          /* Computing MIN */
          i__4 = i__ + *ku;
          i__3 = min(i__4,lenx);
          for (j = max(i__2,1); j <= i__3; ++j) {
            temp = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
              d__1));
            symb_zero__ = symb_zero__ && (x[jx] == 0. || temp ==
              0.);
            y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
            jx += *incx;
          }
        }
        if (! symb_zero__) {
          y[iy] += d_sign(&safe1, &y[iy]);
        }
        iy += *incy;
      }
    } else {
      i__1 = leny;
      for (i__ = 1; i__ <= i__1; ++i__) {
        if (*beta == 0.) {
          symb_zero__ = TRUE_;
          y[iy] = 0.;
        } else if (y[iy] == 0.) {
          symb_zero__ = TRUE_;
        } else {
          symb_zero__ = FALSE_;
          y[iy] = *beta * (d__1 = y[iy], abs(d__1));
        }
        if (*alpha != 0.) {
          jx = kx;
          /* Computing MAX */
          i__3 = i__ - *kl;
          /* Computing MIN */
          i__4 = i__ + *ku;
          i__2 = min(i__4,lenx);
          for (j = max(i__3,1); j <= i__2; ++j) {
            temp = (d__1 = ab[ke - i__ + j + i__ * ab_dim1], abs(
              d__1));
            symb_zero__ = symb_zero__ && (x[jx] == 0. || temp ==
              0.);
            y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
            jx += *incx;
          }
        }
        if (! symb_zero__) {
          y[iy] += d_sign(&safe1, &y[iy]);
        }
        iy += *incy;
      }
    }
  }

  return 0;

  /*     End of DLA_GBAMV */

} /* dla_gbamv__ */

//-------------------------------------------------------------------------------------------------

// dla_gbrcond -- LAPACK computational routine (version 3.7.0) --
template<class T>
T la_gbrcond(char *trans, integer *n, integer *kl, integer *ku, T *ab, integer *ldab,
  T *afb, integer *ldafb, integer *ipiv, integer *cmode, T *c__, integer *info, T *work,
  integer *iwork, ftnlen trans_len)
{
  /* Table of constant values */
  static integer c__1 = 1;

  /* System generated locals */
  integer ab_dim1, ab_offset, afb_dim1, afb_offset, i__1, i__2, i__3, i__4;
  T ret_val, d__1;

  /* Local variables */
  static integer i__, j, kd, ke;
  static T tmp;
  static integer kase;
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static integer isave[3];
  //extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
  //  integer *, doublereal *, integer *, integer *), xerbla_(char *,
  //    integer *, ftnlen), dgbtrs_(char *, integer *, integer *, integer
  //      *, integer *, doublereal *, integer *, integer *, doublereal *,
  //      integer *, integer *, ftnlen);
  static T ainvnm;
  static logical notrans;

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  afb_dim1 = *ldafb;
  afb_offset = 1 + afb_dim1;
  afb -= afb_offset;
  --ipiv;
  --c__;
  --work;
  --iwork;

  /* Function Body */
  ret_val = 0.;

  *info = 0;
  notrans = lsame(trans, "N", (ftnlen)1, (ftnlen)1);
  if (! notrans && ! lsame(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame(
    trans, "C", (ftnlen)1, (ftnlen)1)) {
    *info = -1;
  } else if (*n < 0) {
    *info = -2;
  } else if (*kl < 0 || *kl > *n - 1) {
    *info = -3;
  } else if (*ku < 0 || *ku > *n - 1) {
    *info = -4;
  } else if (*ldab < *kl + *ku + 1) {
    *info = -6;
  } else if (*ldafb < (*kl << 1) + *ku + 1) {
    *info = -8;
  }
  if (*info != 0) {
    i__1 = -(*info);
    xerbla("DLA_GBRCOND", &i__1, (ftnlen)11);
    return ret_val;
  }
  if (*n == 0) {
    ret_val = 1.;
    return ret_val;
  }

  /*     Compute the equilibration matrix R such that */
  /*     inv(R)*A*C has unit 1-norm. */

  kd = *ku + 1;
  ke = *kl + 1;
  if (notrans) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      tmp = 0.;
      if (*cmode == 1) {
        /* Computing MAX */
        i__2 = i__ - *kl;
        /* Computing MIN */
        i__4 = i__ + *ku;
        i__3 = min(i__4,*n);
        for (j = max(i__2,1); j <= i__3; ++j) {
          tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1] * c__[j],
            abs(d__1));
        }
      } else if (*cmode == 0) {
        /* Computing MAX */
        i__3 = i__ - *kl;
        /* Computing MIN */
        i__4 = i__ + *ku;
        i__2 = min(i__4,*n);
        for (j = max(i__3,1); j <= i__2; ++j) {
          tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(d__1));
        }
      } else {
        /* Computing MAX */
        i__2 = i__ - *kl;
        /* Computing MIN */
        i__4 = i__ + *ku;
        i__3 = min(i__4,*n);
        for (j = max(i__2,1); j <= i__3; ++j) {
          tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1] / c__[j],
            abs(d__1));
        }
      }
      work[(*n << 1) + i__] = tmp;
    }
  } else {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      tmp = 0.;
      if (*cmode == 1) {
        /* Computing MAX */
        i__3 = i__ - *kl;
        /* Computing MIN */
        i__4 = i__ + *ku;
        i__2 = min(i__4,*n);
        for (j = max(i__3,1); j <= i__2; ++j) {
          tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1] * c__[j],
            abs(d__1));
        }
      } else if (*cmode == 0) {
        /* Computing MAX */
        i__2 = i__ - *kl;
        /* Computing MIN */
        i__4 = i__ + *ku;
        i__3 = min(i__4,*n);
        for (j = max(i__2,1); j <= i__3; ++j) {
          tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1], abs(d__1)
            );
        }
      } else {
        /* Computing MAX */
        i__3 = i__ - *kl;
        /* Computing MIN */
        i__4 = i__ + *ku;
        i__2 = min(i__4,*n);
        for (j = max(i__3,1); j <= i__2; ++j) {
          tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1] / c__[j],
            abs(d__1));
        }
      }
      work[(*n << 1) + i__] = tmp;
    }
  }

  /*     Estimate the norm of inv(op(A)). */

  ainvnm = 0.;
  kase = 0;
L10:
  lacn2(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
  if (kase != 0) {
    if (kase == 2) {

      /*           Multiply by R. */

      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        work[i__] *= work[(*n << 1) + i__];
      }
      if (notrans) {
        gbtrs("No transpose", n, kl, ku, &c__1, &afb[afb_offset],
          ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
      } else {
        gbtrs("Transpose", n, kl, ku, &c__1, &afb[afb_offset],
          ldafb, &ipiv[1], &work[1], n, info, (ftnlen)9);
      }

      /*           Multiply by inv(C). */

      if (*cmode == 1) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] /= c__[i__];
        }
      } else if (*cmode == -1) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] *= c__[i__];
        }
      }
    } else {

      /*           Multiply by inv(C**T). */

      if (*cmode == 1) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] /= c__[i__];
        }
      } else if (*cmode == -1) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] *= c__[i__];
        }
      }
      if (notrans) {
        gbtrs("Transpose", n, kl, ku, &c__1, &afb[afb_offset],
          ldafb, &ipiv[1], &work[1], n, info, (ftnlen)9);
      } else {
        gbtrs("No transpose", n, kl, ku, &c__1, &afb[afb_offset],
          ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
      }

      /*           Multiply by R. */

      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        work[i__] *= work[(*n << 1) + i__];
      }
    }
    goto L10;
  }

  /*     Compute the estimate of the reciprocal condition number. */

  if (ainvnm != 0.) {
    ret_val = 1. / ainvnm;
  }

  return ret_val;

} /* dla_gbrcond__ */


//-------------------------------------------------------------------------------------------------

// dla_gbrfsx_extended -- LAPACK computational routine (version 3.7.1) --
template<class T>
int la_gbrfsx_extended(integer *prec_type__, integer *trans_type__, integer *n, integer *kl,
  integer *ku, integer *nrhs, T *ab, integer *ldab, T *afb, integer *ldafb, integer *ipiv,
  logical *colequ, T *c__, T *b, integer *ldb, T *y, integer *ldy, T *berr_out__,
  integer *n_norms__, T *err_bnds_norm__, T *err_bnds_comp__, T *res, T *ayb, T *dy, T *y_tail__,
  T *rcond, integer *ithresh, T *rthresh, T *dz_ub__, logical *ignore_cwise__, integer *info)
{
  /* Table of constant values */
  static integer c__1 = 1;
  static T c_b6 = -1.;
  static T c_b8 = 1.;

  /* System generated locals */
  integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset,
    y_dim1, y_offset, err_bnds_norm_dim1, err_bnds_norm_offset,
    err_bnds_comp_dim1, err_bnds_comp_offset, i__1, i__2, i__3;
  T d__1, d__2;
  char ch__1[1];

  /* Local variables */
  static T dxratmax, dzratmax;
  static integer i__, j, m;
  //extern /* Subroutine */ int dla_gbamv__(integer *, integer *, integer *,
  //  integer *, integer *, doublereal *, doublereal *, integer *,
  //  doublereal *, integer *, doublereal *, doublereal *, integer *);
  static logical incr_prec__;
  static T prev_dz_z__, yk, final_dx_x__;
  //extern /* Subroutine */ int dla_wwaddw__(integer *, doublereal *,
  //  doublereal *, doublereal *);
  static T final_dz_z__, prevnormdx;
  static integer cnt;
  static T dyk, eps, incr_thresh__, dx_x__, dz_z__;
  //extern /* Subroutine */ int dla_lin_berr__(integer *, integer *, integer *
  //  , doublereal *, doublereal *, doublereal *);
  static T ymin;
  //extern /* Subroutine */ int blas_dgbmv_x__(integer *, integer *, integer *
  //  , integer *, integer *, doublereal *, doublereal *, integer *,
  //  doublereal *, integer *, doublereal *, doublereal *, integer *,
  //  integer *);
  static integer y_prec_state__;
  //extern /* Subroutine */ int blas_dgbmv2_x__(integer *, integer *, integer
  //  *, integer *, integer *, doublereal *, doublereal *, integer *,
  //  doublereal *, doublereal *, integer *, doublereal *, doublereal *,
  //  integer *, integer *), dgbmv_(char *, integer *, integer *,
  //    integer *, integer *, doublereal *, doublereal *, integer *,
  //    doublereal *, integer *, doublereal *, doublereal *, integer *,
  //    ftnlen), dcopy_(integer *, doublereal *, integer *, doublereal *,
  //      integer *);
  static T dxrat, dzrat;
  //extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *,
  //  integer *, doublereal *, integer *);
  static char trans[1];
  static T normx, normy;
  //extern doublereal dlamch_(char *, ftnlen);
  //extern /* Subroutine */ int dgbtrs_(char *, integer *, integer *, integer
  //  *, integer *, doublereal *, integer *, integer *, doublereal *,
  //  integer *, integer *, ftnlen);
  static T normdx;
  //extern /* Character */ VOID chla_transtype__(char *, ftnlen, integer *);
  static T hugeval;
  static integer x_state__, z_state__;

  /* Parameter adjustments */
  err_bnds_comp_dim1 = *nrhs;
  err_bnds_comp_offset = 1 + err_bnds_comp_dim1;
  err_bnds_comp__ -= err_bnds_comp_offset;
  err_bnds_norm_dim1 = *nrhs;
  err_bnds_norm_offset = 1 + err_bnds_norm_dim1;
  err_bnds_norm__ -= err_bnds_norm_offset;
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  afb_dim1 = *ldafb;
  afb_offset = 1 + afb_dim1;
  afb -= afb_offset;
  --ipiv;
  --c__;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;
  y_dim1 = *ldy;
  y_offset = 1 + y_dim1;
  y -= y_offset;
  --berr_out__;
  --res;
  --ayb;
  --dy;
  --y_tail__;

  /* Function Body */
  if (*info != 0) {
    return 0;
  }
  chla_transtype(ch__1, (ftnlen)1, trans_type__);
  *(unsigned char *)trans = *(unsigned char *)&ch__1[0];
  eps = lamch("Epsilon", (ftnlen)7);
  hugeval = lamch("Overflow", (ftnlen)8);
  /*     Force HUGEVAL to Inf */
  hugeval *= hugeval;
  /*     Using HUGEVAL may lead to spurious underflows. */
  incr_thresh__ = (doublereal) (*n) * eps;
  m = *kl + *ku + 1;
  i__1 = *nrhs;
  for (j = 1; j <= i__1; ++j) {
    y_prec_state__ = 1;
    if (y_prec_state__ == 2) {
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        y_tail__[i__] = 0.;
      }
    }
    dxrat = 0.;
    dxratmax = 0.;
    dzrat = 0.;
    dzratmax = 0.;
    final_dx_x__ = hugeval;
    final_dz_z__ = hugeval;
    prevnormdx = hugeval;
    prev_dz_z__ = hugeval;
    dz_z__ = hugeval;
    dx_x__ = hugeval;
    x_state__ = 1;
    z_state__ = 0;
    incr_prec__ = FALSE_;
    i__2 = *ithresh;
    for (cnt = 1; cnt <= i__2; ++cnt) {

      /*        Compute residual RES = B_s - op(A_s) * Y, */
      /*            op(A) = A, A**T, or A**H depending on TRANS (and type). */

      copy(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
      if (y_prec_state__ == 0) {
        gbmv(trans, &m, n, kl, ku, &c_b6, &ab[ab_offset], ldab, &y[
          j * y_dim1 + 1], &c__1, &c_b8, &res[1], &c__1, (
            ftnlen)1);
      } else if (y_prec_state__ == 1) {
        // call to blas_dgbmv_x was edited by Robin Schmidt to make it compatible to the xblas
        // routine - the original code is retained as comment below (todo: try to templatize)
        blas_gbmv_x(blas_colmajor,  // blas_colmajor added by Robin - guess!!!
          toTransType(trans_type__), *n, *n, *kl, *ku, c_b6, &ab[ab_offset], *ldab,
          &y[j * y_dim1 + 1], c__1, c_b8, &res[1], c__1, toPrecType(prec_type__));
        //blas_dgbmv_x(trans_type__, n, n, kl, ku, &c_b6, &ab[
        //  ab_offset], ldab, &y[j * y_dim1 + 1], &c__1, &c_b8, &
        //    res[1], &c__1, prec_type__);
      } else {
        // also edited by Robin, original below..
        blas_gbmv2_x(blas_colmajor,  // blas_colmajor added by Robin - guess!!!
          toTransType(trans_type__), *n, *n, *kl, *ku, c_b6, &ab[ab_offset], *ldab,
          &y[j * y_dim1 + 1], &y_tail__[1], c__1, c_b8, &res[1], c__1, toPrecType(prec_type__));
        //blas_dgbmv2_x(trans_type__, n, n, kl, ku, &c_b6, &ab[
        //  ab_offset], ldab, &y[j * y_dim1 + 1], &y_tail__[1], &
        //    c__1, &c_b8, &res[1], &c__1, prec_type__);
      }
      /*        XXX: RES is no longer needed. */
      copy(n, &res[1], &c__1, &dy[1], &c__1);
      gbtrs(trans, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1]
        , &dy[1], n, info, (ftnlen)1);

      /*         Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT. */

      normx = 0.;
      normy = 0.;
      normdx = 0.;
      dz_z__ = 0.;
      ymin = hugeval;
      i__3 = *n;
      for (i__ = 1; i__ <= i__3; ++i__) {
        yk = (d__1 = y[i__ + j * y_dim1], abs(d__1));
        dyk = (d__1 = dy[i__], abs(d__1));
        if (yk != 0.) {
          /* Computing MAX */
          d__1 = dz_z__, d__2 = dyk / yk;
          dz_z__ = max(d__1,d__2);
        } else if (dyk != 0.) {
          dz_z__ = hugeval;
        }
        ymin = min(ymin,yk);
        normy = max(normy,yk);
        if (*colequ) {
          /* Computing MAX */
          d__1 = normx, d__2 = yk * c__[i__];
          normx = max(d__1,d__2);
          /* Computing MAX */
          d__1 = normdx, d__2 = dyk * c__[i__];
          normdx = max(d__1,d__2);
        } else {
          normx = normy;
          normdx = max(normdx,dyk);
        }
      }
      if (normx != 0.) {
        dx_x__ = normdx / normx;
      } else if (normdx == 0.) {
        dx_x__ = 0.;
      } else {
        dx_x__ = hugeval;
      }
      dxrat = normdx / prevnormdx;
      dzrat = dz_z__ / prev_dz_z__;

      /*         Check termination criteria. */

      if (! (*ignore_cwise__) && ymin * *rcond < incr_thresh__ * normy
        && y_prec_state__ < 2) {
        incr_prec__ = TRUE_;
      }
      if (x_state__ == 3 && dxrat <= *rthresh) {
        x_state__ = 1;
      }
      if (x_state__ == 1) {
        if (dx_x__ <= eps) {
          x_state__ = 2;
        } else if (dxrat > *rthresh) {
          if (y_prec_state__ != 2) {
            incr_prec__ = TRUE_;
          } else {
            x_state__ = 3;
          }
        } else {
          if (dxrat > dxratmax) {
            dxratmax = dxrat;
          }
        }
        if (x_state__ > 1) {
          final_dx_x__ = dx_x__;
        }
      }
      if (z_state__ == 0 && dz_z__ <= *dz_ub__) {
        z_state__ = 1;
      }
      if (z_state__ == 3 && dzrat <= *rthresh) {
        z_state__ = 1;
      }
      if (z_state__ == 1) {
        if (dz_z__ <= eps) {
          z_state__ = 2;
        } else if (dz_z__ > *dz_ub__) {
          z_state__ = 0;
          dzratmax = 0.;
          final_dz_z__ = hugeval;
        } else if (dzrat > *rthresh) {
          if (y_prec_state__ != 2) {
            incr_prec__ = TRUE_;
          } else {
            z_state__ = 3;
          }
        } else {
          if (dzrat > dzratmax) {
            dzratmax = dzrat;
          }
        }
        if (z_state__ > 1) {
          final_dz_z__ = dz_z__;
        }
      }

      /*           Exit if both normwise and componentwise stopped working, */
      /*           but if componentwise is unstable, let it go at least two */
      /*           iterations. */

      if (x_state__ != 1) {
        if (*ignore_cwise__) {
          goto L666;
        }
        if (z_state__ == 3 || z_state__ == 2) {
          goto L666;
        }
        if (z_state__ == 0 && cnt > 1) {
          goto L666;
        }
      }
      if (incr_prec__) {
        incr_prec__ = FALSE_;
        ++y_prec_state__;
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__) {
          y_tail__[i__] = 0.;
        }
      }
      prevnormdx = normdx;
      prev_dz_z__ = dz_z__;

      /*           Update soluton. */

      if (y_prec_state__ < 2) {
        axpy(n, &c_b8, &dy[1], &c__1, &y[j * y_dim1 + 1], &c__1);
      } else {
        la_wwaddw(n, &y[j * y_dim1 + 1], &y_tail__[1], &dy[1]);
      }
    }
    /*        Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT. */
  L666:

    /*     Set final_* when cnt hits ithresh. */

    if (x_state__ == 1) {
      final_dx_x__ = dx_x__;
    }
    if (z_state__ == 1) {
      final_dz_z__ = dz_z__;
    }

    /*     Compute error bounds. */

    if (*n_norms__ >= 1) {
      err_bnds_norm__[j + (err_bnds_norm_dim1 << 1)] = final_dx_x__ / (
        1 - dxratmax);
    }
    if (*n_norms__ >= 2) {
      err_bnds_comp__[j + (err_bnds_comp_dim1 << 1)] = final_dz_z__ / (
        1 - dzratmax);
    }

    /*     Compute componentwise relative backward error from formula */
    /*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) ) */
    /*     where abs(Z) is the componentwise absolute value of the matrix */
    /*     or vector Z. */

    /*        Compute residual RES = B_s - op(A_s) * Y, */
    /*            op(A) = A, A**T, or A**H depending on TRANS (and type). */

    copy(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
    gbmv(trans, n, n, kl, ku, &c_b6, &ab[ab_offset], ldab, &y[j *
      y_dim1 + 1], &c__1, &c_b8, &res[1], &c__1, (ftnlen)1);
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      ayb[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
    }

    /*     Compute abs(op(A_s))*abs(Y) + abs(B_s). */

    la_gbamv(trans_type__, n, n, kl, ku, &c_b8, &ab[ab_offset], ldab, &
      y[j * y_dim1 + 1], &c__1, &c_b8, &ayb[1], &c__1);
    la_lin_berr(n, n, &c__1, &res[1], &ayb[1], &berr_out__[j]);

    /*     End of loop for each RHS */

  }

  return 0;
} /* dla_gbrfsx_extended__ */

//-------------------------------------------------------------------------------------------------

// from dla_gbrpvgrw -- LAPACK computational routine (version 3.7.0)
template<class T>
T la_gbrpvgrw(integer *n, integer *kl, integer *ku, integer *ncols, T *ab, integer *ldab,
  T *afb, integer *ldafb)
{
  /* System generated locals */
  integer ab_dim1, ab_offset, afb_dim1, afb_offset, i__1, i__2, i__3, i__4;
  T ret_val, d__1, d__2;

  /* Local variables */
  static integer i__, j, kd;
  static T amax, umax, rpvgrw;

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  afb_dim1 = *ldafb;
  afb_offset = 1 + afb_dim1;
  afb -= afb_offset;

  /* Function Body */
  rpvgrw = 1.;
  kd = *ku + 1;
  i__1 = *ncols;
  for (j = 1; j <= i__1; ++j) {
    amax = 0.;
    umax = 0.;
    /* Computing MAX */
    i__2 = j - *ku;
    /* Computing MIN */
    i__4 = j + *kl;
    i__3 = min(i__4,*n);
    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
      /* Computing MAX */
      d__2 = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(d__1));
      amax = max(d__2,amax);
    }
    /* Computing MAX */
    i__3 = j - *ku;
    i__2 = j;
    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
      /* Computing MAX */
      d__2 = (d__1 = afb[kd + i__ - j + j * afb_dim1], abs(d__1));
      umax = max(d__2,umax);
    }
    if (umax != 0.) {
      /* Computing MIN */
      d__1 = amax / umax;
      rpvgrw = min(d__1,rpvgrw);
    }
  }
  ret_val = rpvgrw;
  return ret_val;
} /* dla_gbrpvgrw__ */


//-------------------------------------------------------------------------------------------------

// dla_lin_berr__ -- LAPACK computational routine (version 3.7.0)
template<class T>
int la_lin_berr(integer *n, integer *nz, integer *nrhs, T *res, T *ayb, T *berr)
{
  /* System generated locals */
  integer ayb_dim1, ayb_offset, res_dim1, res_offset, i__1, i__2;
  T d__1;

  /* Local variables */
  static integer i__, j;
  static T tmp, safe1;
  //extern T dlamch_(char *, ftnlen);

  /*     Adding SAFE1 to the numerator guards against spuriously zero */
  /*     residuals.  A similar safeguard is in the SLA_yyAMV routine used */
  /*     to compute AYB. */

  /* Parameter adjustments */
  --berr;
  ayb_dim1 = *n;
  ayb_offset = 1 + ayb_dim1;
  ayb -= ayb_offset;
  res_dim1 = *n;
  res_offset = 1 + res_dim1;
  res -= res_offset;

  /* Function Body */
  safe1 = lamch("Safe minimum", (ftnlen)12);
  safe1 = (*nz + 1) * safe1;
  i__1 = *nrhs;
  for (j = 1; j <= i__1; ++j) {
    berr[j] = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      if (ayb[i__ + j * ayb_dim1] != 0.) {
        tmp = (safe1 + (d__1 = res[i__ + j * res_dim1], abs(d__1))) /
          ayb[i__ + j * ayb_dim1];
        /* Computing MAX */
        d__1 = berr[j];
        berr[j] = max(d__1,tmp);
      }
      /*     If AYB is exactly 0.0 (and if computed by SLA_yyAMV), then we know */
      /*     the true residual also must be exactly 0.0. */
    }
  }
  return 0;
} /* dla_lin_berr__ */


//-------------------------------------------------------------------------------------------------

// dlascl2 -- LAPACK computational routine (version 3.7.0) --
template<class T>
int lascl2(integer *m, integer *n, T *d__, T *x, integer *ldx)
{
  /* System generated locals */
  integer x_dim1, x_offset, i__1, i__2;

  /* Local variables */
  static integer i__, j;

  /* Parameter adjustments */
  --d__;
  x_dim1 = *ldx;
  x_offset = 1 + x_dim1;
  x -= x_offset;

  /* Function Body */
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    i__2 = *m;
    for (i__ = 1; i__ <= i__2; ++i__) {
      x[i__ + j * x_dim1] *= d__[i__];
    }
  }
  return 0;
} /* dlascl2_ */


//-------------------------------------------------------------------------------------------------

// LAPACK computational routine (version 3.7.0)
template<class T>
int la_wwaddw(integer *n, T *x, T *y, T *w)
{
  /* System generated locals */
  integer i__1;

  /* Local variables */
  static integer i__;
  static doublereal s;

  /* Parameter adjustments */
  --w;
  --y;
  --x;

  /* Function Body */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s = x[i__] + w[i__];
    s = s + s - s;
    y[i__] = x[i__] - s + w[i__] + y[i__];
    x[i__] = s;
    /* L10: */
  }
  return 0;
} /* dla_wwaddw__ */

//=================================================================================================
// Auxiliary routines:

//-------------------------------------------------------------------------------------------------

// -- LAPACK auxiliary routine (version 3.7.0)
template<class T>
int labad(T *small, T *large)
{
  // If it looks like we're on a Cray, take the square root of
  // SMALL and LARGE to avoid overflow and underflow problems.
  if(d_lg10(large) > 2e3) {
    *small = sqrt(*small);
    *large = sqrt(*large);
  }
  return 0;
}

//-------------------------------------------------------------------------------------------------

// dlacn2_ - LAPACK auxiliary routine (version 3.7.0) */
template<class T>
int lacn2(integer *n, T *v, T *x, integer *isgn, T *est, integer *kase, integer *isave)
{
  // Table of constant values
  static integer c__1 = 1;
  static T c_b11 = 1.;

  // System generated locals
  int i__1;
  T d__1;

  // Local variables
  static integer i__;
  static  T temp;
  static integer jlast;
  static T altsgn, estold;

  // Parameter adjustments
  --isave;
  --isgn;
  --x;
  --v;

  /* Function Body */
  if (*kase == 0) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      x[i__] = 1. / (doublereal) (*n);
      /* L10: */
    }
    *kase = 1;
    isave[1] = 1;
    return 0;
  }

  switch (isave[1]) {
  case 1:  goto L20;
  case 2:  goto L40;
  case 3:  goto L70;
  case 4:  goto L110;
  case 5:  goto L140;
  }

  /*     ................ ENTRY   (ISAVE( 1 ) = 1) */
  /*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

L20:
  if (*n == 1) {
    v[1] = x[1];
    *est = abs(v[1]);
    /*        ... QUIT */
    goto L150;
  }
  *est = asum(n, &x[1], &c__1);

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    x[i__] = d_sign(&c_b11, &x[i__]);
    isgn[i__] = i_dnnt(&x[i__]);
    /* L30: */
  }
  *kase = 2;
  isave[1] = 2;
  return 0;

  /*     ................ ENTRY   (ISAVE( 1 ) = 2) */
  /*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */

L40:
  isave[2] = iamax(n, &x[1], &c__1);
  isave[3] = 2;

  /*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

L50:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    x[i__] = 0.;
    /* L60: */
  }
  x[isave[2]] = 1.;
  *kase = 1;
  isave[1] = 3;
  return 0;

  /*     ................ ENTRY   (ISAVE( 1 ) = 3) */
  /*     X HAS BEEN OVERWRITTEN BY A*X. */

L70:
  copy(n, &x[1], &c__1, &v[1], &c__1);
  estold = *est;
  *est = asum(n, &v[1], &c__1);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    d__1 = d_sign(&c_b11, &x[i__]);
    if (i_dnnt(&d__1) != isgn[i__]) {
      goto L90;
    }
    /* L80: */
  }
  /*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
  goto L120;

L90:
  /*     TEST FOR CYCLING. */
  if (*est <= estold) {
    goto L120;
  }

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    x[i__] = d_sign(&c_b11, &x[i__]);
    isgn[i__] = i_dnnt(&x[i__]);
    /* L100: */
  }
  *kase = 2;
  isave[1] = 4;
  return 0;

  /*     ................ ENTRY   (ISAVE( 1 ) = 4) */
  /*     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */

L110:
  jlast = isave[2];
  isave[2] = iamax(n, &x[1], &c__1);
  if (x[jlast] != (d__1 = x[isave[2]], abs(d__1)) && isave[3] < 5) {
    ++isave[3];
    goto L50;
  }

  /*     ITERATION COMPLETE.  FINAL STAGE. */

L120:
  altsgn = 1.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    x[i__] = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) +
      1.);
    altsgn = -altsgn;
    /* L130: */
  }
  *kase = 1;
  isave[1] = 5;
  return 0;

  /*     ................ ENTRY   (ISAVE( 1 ) = 5) */
  /*     X HAS BEEN OVERWRITTEN BY A*X. */

L140:
  temp = asum(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
  if (temp > *est) {
    copy(n, &x[1], &c__1, &v[1], &c__1);
    *est = temp;
  }

L150:
  *kase = 0;
  return 0;

  /*     End of DLACN2 */

} /* dlacn2_ */

//-------------------------------------------------------------------------------------------------

// translated from dlangb - LAPACK auxiliary routine (version 3.7.0) */
template<class T>
T langb(char *norm, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *work,
  ftnlen norm_len)
{
  /* Table of constant values */
  static integer c__1 = 1;

  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
  T ret_val, d__1;

  /* Builtin functions */
  double sqrt(T);

  /* Local variables */
  static integer i__, j, k, l;
  static T sum, temp, scale;
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static T value;
  //extern logical disnan_(doublereal *);
  //extern /* Subroutine */ int dlassq_(integer *, T *, integer *, T *, T *);

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --work;

  /* Function Body */
  if (*n == 0) {
    value = 0.;
  } else if (lsame(norm, "M", (ftnlen)1, (ftnlen)1)) {

    /*        Find max(abs(A(i,j))). */

    value = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      /* Computing MAX */
      i__2 = *ku + 2 - j;
      /* Computing MIN */
      i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
      i__3 = min(i__4,i__5);
      for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
        temp = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
        if (value < temp || isnan(&temp)) {
          value = temp;
        }
        /* L10: */
      }
      /* L20: */
    }
  } else if (lsame(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
    norm == '1') {

    /*        Find norm1(A). */

    value = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      sum = 0.;
      /* Computing MAX */
      i__3 = *ku + 2 - j;
      /* Computing MIN */
      i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
      i__2 = min(i__4,i__5);
      for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
        sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
        /* L30: */
      }
      if (value < sum || isnan(&sum)) {
        value = sum;
      }
      /* L40: */
    }
  } else if (lsame(norm, "I", (ftnlen)1, (ftnlen)1)) {

    /*        Find normI(A). */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      work[i__] = 0.;
      /* L50: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      k = *ku + 1 - j;
      /* Computing MAX */
      i__2 = 1, i__3 = j - *ku;
      /* Computing MIN */
      i__5 = *n, i__6 = j + *kl;
      i__4 = min(i__5,i__6);
      for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
        work[i__] += (d__1 = ab[k + i__ + j * ab_dim1], abs(d__1));
        /* L60: */
      }
      /* L70: */
    }
    value = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = work[i__];
      if (value < temp || isnan(&temp)) {
        value = temp;
      }
      /* L80: */
    }
  } else if (lsame(norm, "F", (ftnlen)1, (ftnlen)1) || lsame(norm, "E", (
    ftnlen)1, (ftnlen)1)) {

    /*        Find normF(A). */

    scale = 0.;
    sum = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      /* Computing MAX */
      i__4 = 1, i__2 = j - *ku;
      l = max(i__4,i__2);
      k = *ku + 1 - j + l;
      /* Computing MIN */
      i__2 = *n, i__3 = j + *kl;
      i__4 = min(i__2,i__3) - l + 1;
      lassq(&i__4, &ab[k + j * ab_dim1], &c__1, &scale, &sum);
      /* L90: */
    }
    value = scale * sqrt(sum);
  }

  ret_val = value;
  return ret_val;

  /*     End of DLANGB */

} /* dlangb_ */

//-------------------------------------------------------------------------------------------------

// translated from dlacpy - LAPACK auxiliary routine (version 3.7.0)
template<class T>
int lacpy(char *uplo, integer *m, integer *n, T *a, integer *lda, T *b, integer *ldb,
  ftnlen uplo_len)
{
  /* System generated locals */
  integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

  /* Local variables */
  static integer i__, j;
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;

  /* Function Body */
  if (lsame(uplo, "U", (ftnlen)1, (ftnlen)1)) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      i__2 = min(j,*m);
      for (i__ = 1; i__ <= i__2; ++i__) {
        b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
        /* L10: */
      }
      /* L20: */
    }
  } else if (lsame(uplo, "L", (ftnlen)1, (ftnlen)1)) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      i__2 = *m;
      for (i__ = j; i__ <= i__2; ++i__) {
        b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
        /* L30: */
      }
      /* L40: */
    }
  } else {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      i__2 = *m;
      for (i__ = 1; i__ <= i__2; ++i__) {
        b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
        /* L50: */
      }
      /* L60: */
    }
  }
  return 0;

} /* dlacpy_ */

//-------------------------------------------------------------------------------------------------

// from dlantb - LAPACK auxiliary routine (version 3.7.0)
template<class T>
T lantb(char *norm, char *uplo, char *diag, integer *n, integer *k, T *ab, integer *ldab,
  T *work, ftnlen norm_len, ftnlen uplo_len, ftnlen diag_len)
{
  /* Table of constant values */
  static integer c__1 = 1;

  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5;
  T ret_val, d__1;

  /* Builtin functions */
  double sqrt(T);

  /* Local variables */
  static integer i__, j, l;
  static T sum, scale;
  static logical udiag;
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static T value;
  //extern logical disnan_(T *);
  //extern /* Subroutine */ int dlassq_(integer *, T *, integer *, T *, T *);

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --work;

  /* Function Body */
  if (*n == 0) {
    value = 0.;
  } else if (lsame(norm, "M", (ftnlen)1, (ftnlen)1)) {

    /*        Find max(abs(A(i,j))). */

    if (lsame(diag, "U", (ftnlen)1, (ftnlen)1)) {
      value = 1.;
      if (lsame(uplo, "U", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MAX */
          i__2 = *k + 2 - j;
          i__3 = *k;
          for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
            sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            if (value < sum || isnan(&sum)) {
              value = sum;
            }
            /* L10: */
          }
          /* L20: */
        }
      } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MIN */
          i__2 = *n + 1 - j, i__4 = *k + 1;
          i__3 = min(i__2,i__4);
          for (i__ = 2; i__ <= i__3; ++i__) {
            sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            if (value < sum || isnan(&sum)) {
              value = sum;
            }
            /* L30: */
          }
          /* L40: */
        }
      }
    } else {
      value = 0.;
      if (lsame(uplo, "U", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MAX */
          i__3 = *k + 2 - j;
          i__2 = *k + 1;
          for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
            sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            if (value < sum || isnan(&sum)) {
              value = sum;
            }
            /* L50: */
          }
          /* L60: */
        }
      } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MIN */
          i__3 = *n + 1 - j, i__4 = *k + 1;
          i__2 = min(i__3,i__4);
          for (i__ = 1; i__ <= i__2; ++i__) {
            sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            if (value < sum || isnan(&sum)) {
              value = sum;
            }
            /* L70: */
          }
          /* L80: */
        }
      }
    }
  } else if (lsame(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
    norm == '1') {

    /*        Find norm1(A). */

    value = 0.;
    udiag = lsame(diag, "U", (ftnlen)1, (ftnlen)1);
    if (lsame(uplo, "U", (ftnlen)1, (ftnlen)1)) {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        if (udiag) {
          sum = 1.;
          /* Computing MAX */
          i__2 = *k + 2 - j;
          i__3 = *k;
          for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
            sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            /* L90: */
          }
        } else {
          sum = 0.;
          /* Computing MAX */
          i__3 = *k + 2 - j;
          i__2 = *k + 1;
          for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
            sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            /* L100: */
          }
        }
        if (value < sum || isnan(&sum)) {
          value = sum;
        }
        /* L110: */
      }
    } else {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        if (udiag) {
          sum = 1.;
          /* Computing MIN */
          i__3 = *n + 1 - j, i__4 = *k + 1;
          i__2 = min(i__3,i__4);
          for (i__ = 2; i__ <= i__2; ++i__) {
            sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            /* L120: */
          }
        } else {
          sum = 0.;
          /* Computing MIN */
          i__3 = *n + 1 - j, i__4 = *k + 1;
          i__2 = min(i__3,i__4);
          for (i__ = 1; i__ <= i__2; ++i__) {
            sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            /* L130: */
          }
        }
        if (value < sum || isnan(&sum)) {
          value = sum;
        }
        /* L140: */
      }
    }
  } else if (lsame(norm, "I", (ftnlen)1, (ftnlen)1)) {

    /*        Find normI(A). */

    value = 0.;
    if (lsame(uplo, "U", (ftnlen)1, (ftnlen)1)) {
      if (lsame(diag, "U", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] = 1.;
          /* L150: */
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          l = *k + 1 - j;
          /* Computing MAX */
          i__2 = 1, i__3 = j - *k;
          i__4 = j - 1;
          for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
            work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
              d__1));
            /* L160: */
          }
          /* L170: */
        }
      } else {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] = 0.;
          /* L180: */
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          l = *k + 1 - j;
          /* Computing MAX */
          i__4 = 1, i__2 = j - *k;
          i__3 = j;
          for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
            work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
              d__1));
            /* L190: */
          }
          /* L200: */
        }
      }
    } else {
      if (lsame(diag, "U", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] = 1.;
          /* L210: */
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          l = 1 - j;
          /* Computing MIN */
          i__4 = *n, i__2 = j + *k;
          i__3 = min(i__4,i__2);
          for (i__ = j + 1; i__ <= i__3; ++i__) {
            work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
              d__1));
            /* L220: */
          }
          /* L230: */
        }
      } else {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] = 0.;
          /* L240: */
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          l = 1 - j;
          /* Computing MIN */
          i__4 = *n, i__2 = j + *k;
          i__3 = min(i__4,i__2);
          for (i__ = j; i__ <= i__3; ++i__) {
            work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
              d__1));
            /* L250: */
          }
          /* L260: */
        }
      }
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sum = work[i__];
      if (value < sum || isnan(&sum)) {
        value = sum;
      }
      /* L270: */
    }
  } else if (lsame(norm, "F", (ftnlen)1, (ftnlen)1) || lsame(norm, "E", (
    ftnlen)1, (ftnlen)1)) {

    /*        Find normF(A). */

    if (lsame(uplo, "U", (ftnlen)1, (ftnlen)1)) {
      if (lsame(diag, "U", (ftnlen)1, (ftnlen)1)) {
        scale = 1.;
        sum = (doublereal) (*n);
        if (*k > 0) {
          i__1 = *n;
          for (j = 2; j <= i__1; ++j) {
            /* Computing MIN */
            i__4 = j - 1;
            i__3 = min(i__4,*k);
            /* Computing MAX */
            i__2 = *k + 2 - j;
            lassq(&i__3, &ab[max(i__2,1) + j * ab_dim1], &c__1, &scale, &sum);
            /* L280: */
          }
        }
      } else {
        scale = 0.;
        sum = 1.;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MIN */
          i__4 = j, i__2 = *k + 1;
          i__3 = min(i__4,i__2);
          /* Computing MAX */
          i__5 = *k + 2 - j;
          lassq(&i__3, &ab[max(i__5,1) + j * ab_dim1], &c__1, &scale, &sum);
          /* L290: */
        }
      }
    } else {
      if (lsame(diag, "U", (ftnlen)1, (ftnlen)1)) {
        scale = 1.;
        sum = (doublereal) (*n);
        if (*k > 0) {
          i__1 = *n - 1;
          for (j = 1; j <= i__1; ++j) {
            /* Computing MIN */
            i__4 = *n - j;
            i__3 = min(i__4,*k);
            lassq(&i__3, &ab[j * ab_dim1 + 2], &c__1, &scale, &sum);
            /* L300: */
          }
        }
      } else {
        scale = 0.;
        sum = 1.;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MIN */
          i__4 = *n - j + 1, i__2 = *k + 1;
          i__3 = min(i__4,i__2);
          lassq(&i__3, &ab[j * ab_dim1 + 1], &c__1, &scale, &sum);
          /* L310: */
        }
      }
    }
    value = scale * sqrt(sum);
  }

  ret_val = value;
  return ret_val;
} /* dlantb_ */

//-------------------------------------------------------------------------------------------------

// from dlaqgb - LAPACK auxiliary routine (version 3.7.0)
template<class T>
int laqgb(integer *m, integer *n, integer *kl, integer *ku,
  T *ab, integer *ldab, T *r__, T *c__, T *rowcnd, T *colcnd, T *amax, char *equed,
  ftnlen equed_len)
{
  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;

  /* Local variables */
  static integer i__, j;
  static T cj, large, small;
  //extern T dlamch_(char *, ftnlen);

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --r__;
  --c__;

  /* Function Body */
  if (*m <= 0 || *n <= 0) {
    *(unsigned char *)equed = 'N';
    return 0;
  }

  /*     Initialize LARGE and SMALL. */

  small = lamch("Safe minimum", (ftnlen)12) / lamch("Precision", (ftnlen)9);
  large = 1. / small;

  if (*rowcnd >= .1 && *amax >= small && *amax <= large) {

    /*        No row scaling */

    if (*colcnd >= .1) {

      /*           No column scaling */

      *(unsigned char *)equed = 'N';
    } else {

      /*           Column scaling */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        cj = c__[j];
        /* Computing MAX */
        i__2 = 1, i__3 = j - *ku;
        /* Computing MIN */
        i__5 = *m, i__6 = j + *kl;
        i__4 = min(i__5,i__6);
        for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
          ab[*ku + 1 + i__ - j + j * ab_dim1] = cj * ab[*ku + 1 +
            i__ - j + j * ab_dim1];
          /* L10: */
        }
        /* L20: */
      }
      *(unsigned char *)equed = 'C';
    }
  } else if (*colcnd >= .1) {

    /*        Row scaling, no column scaling */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      /* Computing MAX */
      i__4 = 1, i__2 = j - *ku;
      /* Computing MIN */
      i__5 = *m, i__6 = j + *kl;
      i__3 = min(i__5,i__6);
      for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
        ab[*ku + 1 + i__ - j + j * ab_dim1] = r__[i__] * ab[*ku + 1 +
          i__ - j + j * ab_dim1];
        /* L30: */
      }
      /* L40: */
    }
    *(unsigned char *)equed = 'R';
  } else {

    /*        Row and column scaling */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      cj = c__[j];
      /* Computing MAX */
      i__3 = 1, i__4 = j - *ku;
      /* Computing MIN */
      i__5 = *m, i__6 = j + *kl;
      i__2 = min(i__5,i__6);
      for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
        ab[*ku + 1 + i__ - j + j * ab_dim1] = cj * r__[i__] * ab[*ku
          + 1 + i__ - j + j * ab_dim1];
        /* L50: */
      }
      /* L60: */
    }
    *(unsigned char *)equed = 'B';
  }

  return 0;

  /*     End of DLAQGB */

} /* dlaqgb_ */

//-------------------------------------------------------------------------------------------------

// from dlassq - LAPACK auxiliary routine (version 3.7.0)
template<class T>
int lassq(integer *n, T *x, integer *incx, T *scale, T *sumsq)
{
  /* System generated locals */
  integer i__1, i__2;
  T d__1;

  /* Local variables */
  static integer ix;
  static T absxi;
  //extern logical disnan_(doublereal *);

  /* Parameter adjustments */
  --x;

  /* Function Body */
  if (*n > 0) {
    i__1 = (*n - 1) * *incx + 1;
    i__2 = *incx;
    for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
      absxi = (d__1 = x[ix], abs(d__1));
      if (absxi > 0. || isnan(&absxi)) {
        if (*scale < absxi) {
          /* Computing 2nd power */
          d__1 = *scale / absxi;
          *sumsq = *sumsq * (d__1 * d__1) + 1;
          *scale = absxi;
        } else {
          /* Computing 2nd power */
          d__1 = absxi / *scale;
          *sumsq += d__1 * d__1;
        }
      }
      /* L10: */
    }
  }
  return 0;

  /*     End of DLASSQ */

} /* dlassq_ */

//-------------------------------------------------------------------------------------------------

// translated from dlaswp, LAPACK auxiliary routine (version 3.7.1)
template<class T>
int laswp(integer *n, T *a, integer *lda, integer *k1, integer *k2, integer *ipiv, integer *incx)
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

  /* Local variables */
  static integer i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
  static T temp;

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --ipiv;

  /* Function Body */
  if (*incx > 0) {
    ix0 = *k1;
    i1 = *k1;
    i2 = *k2;
    inc = 1;
  } else if (*incx < 0) {
    ix0 = *k1 + (*k1 - *k2) * *incx;
    i1 = *k2;
    i2 = *k1;
    inc = -1;
  } else {
    return 0;
  }

  n32 = *n / 32 << 5;
  if (n32 != 0) {
    i__1 = n32;
    for (j = 1; j <= i__1; j += 32) {
      ix = ix0;
      i__2 = i2;
      i__3 = inc;
      for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3)
      {
        ip = ipiv[ix];
        if (ip != i__) {
          i__4 = j + 31;
          for (k = j; k <= i__4; ++k) {
            temp = a[i__ + k * a_dim1];
            a[i__ + k * a_dim1] = a[ip + k * a_dim1];
            a[ip + k * a_dim1] = temp;
            /* L10: */
          }
        }
        ix += *incx;
        /* L20: */
      }
      /* L30: */
    }
  }
  if (n32 != *n) {
    ++n32;
    ix = ix0;
    i__1 = i2;
    i__3 = inc;
    for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
      ip = ipiv[ix];
      if (ip != i__) {
        i__2 = *n;
        for (k = n32; k <= i__2; ++k) {
          temp = a[i__ + k * a_dim1];
          a[i__ + k * a_dim1] = a[ip + k * a_dim1];
          a[ip + k * a_dim1] = temp;
          /* L40: */
        }
      }
      ix += *incx;
      /* L50: */
    }
  }

  return 0;

  /*     End of DLASWP */

} /* dlaswp_ */

//-------------------------------------------------------------------------------------------------

/* from dlatbs_ - LAPACK auxiliary routine (version 3.7.0) */
template<class T>
int latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, T *ab,
  integer *ldab, T *x, T *scale, T *cnorm, integer *info, ftnlen uplo_len, ftnlen trans_len,
  ftnlen diag_len, ftnlen normin_len)
{
  // Table of constant values
  static integer c__1 = 1;
  static doublereal c_b36 = .5;


  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
  doublereal d__1, d__2, d__3;

  /* Local variables */
  static integer i__, j;
  static T xj, rec, tjj;
  static integer jinc, jlen;
  //extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *,
  //  integer *);
  static T xbnd;
  static integer imax;
  static T tmax, tjjs, xmax, grow, sumj;
  //extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *,
  //  integer *);
  static integer maind;
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static T tscal, uscal;
  //extern doublereal dasum_(integer *, doublereal *, integer *);
  static integer jlast;
  //extern /* Subroutine */ int dtbsv_(char *, char *, char *, integer *,
  //  integer *, doublereal *, integer *, doublereal *, integer *,
  //  ftnlen, ftnlen, ftnlen), daxpy_(integer *, doublereal *,
  //    doublereal *, integer *, doublereal *, integer *);
  static logical upper;
  //extern doublereal dlamch_(char *, ftnlen);
  //extern integer idamax_(integer *, doublereal *, integer *);
  //extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
  static T bignum;
  static logical notran;
  static integer jfirst;
  static T smlnum;
  static logical nounit;

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --x;
  --cnorm;

  /* Function Body */
  *info = 0;
  upper = lsame(uplo, "U", (ftnlen)1, (ftnlen)1);
  notran = lsame(trans, "N", (ftnlen)1, (ftnlen)1);
  nounit = lsame(diag, "N", (ftnlen)1, (ftnlen)1);

  /*     Test the input parameters. */

  if (! upper && ! lsame(uplo, "L", (ftnlen)1, (ftnlen)1)) {
    *info = -1;
  } else if (! notran && ! lsame(trans, "T", (ftnlen)1, (ftnlen)1) && !
    lsame(trans, "C", (ftnlen)1, (ftnlen)1)) {
    *info = -2;
  } else if (! nounit && ! lsame(diag, "U", (ftnlen)1, (ftnlen)1)) {
    *info = -3;
  } else if (! lsame(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame(normin,
    "N", (ftnlen)1, (ftnlen)1)) {
    *info = -4;
  } else if (*n < 0) {
    *info = -5;
  } else if (*kd < 0) {
    *info = -6;
  } else if (*ldab < *kd + 1) {
    *info = -8;
  }
  if (*info != 0) {
    i__1 = -(*info);
    xerbla("DLATBS", &i__1, (ftnlen)6);
    return 0;
  }

  /*     Quick return if possible */

  if (*n == 0) {
    return 0;
  }

  /*     Determine machine dependent parameters to control overflow. */

  smlnum = lamch("Safe minimum", (ftnlen)12) / lamch("Precision", (
    ftnlen)9);
  bignum = 1. / smlnum;
  *scale = 1.;

  if (lsame(normin, "N", (ftnlen)1, (ftnlen)1)) {

    /*        Compute the 1-norm of each column, not including the diagonal. */

    if (upper) {

      /*           A is upper triangular. */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        /* Computing MIN */
        i__2 = *kd, i__3 = j - 1;
        jlen = min(i__2,i__3);
        cnorm[j] = asum(&jlen, &ab[*kd + 1 - jlen + j * ab_dim1], &
          c__1);
        /* L10: */
      }
    } else {

      /*           A is lower triangular. */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        /* Computing MIN */
        i__2 = *kd, i__3 = *n - j;
        jlen = min(i__2,i__3);
        if (jlen > 0) {
          cnorm[j] = asum(&jlen, &ab[j * ab_dim1 + 2], &c__1);
        } else {
          cnorm[j] = 0.;
        }
        /* L20: */
      }
    }
  }

  /*     Scale the column norms by TSCAL if the maximum element in CNORM is */
  /*     greater than BIGNUM. */

  imax = iamax(n, &cnorm[1], &c__1);
  tmax = cnorm[imax];
  if (tmax <= bignum) {
    tscal = 1.;
  } else {
    tscal = 1. / (smlnum * tmax);
    scal(n, &tscal, &cnorm[1], &c__1);
  }

  /*     Compute a bound on the computed solution vector to see if the */
  /*     Level 2 BLAS routine DTBSV can be used. */

  j = iamax(n, &x[1], &c__1);
  xmax = (d__1 = x[j], abs(d__1));
  xbnd = xmax;
  if (notran) {

    /*        Compute the growth in A * x = b. */

    if (upper) {
      jfirst = *n;
      jlast = 1;
      jinc = -1;
      maind = *kd + 1;
    } else {
      jfirst = 1;
      jlast = *n;
      jinc = 1;
      maind = 1;
    }

    if (tscal != 1.) {
      grow = 0.;
      goto L50;
    }

    if (nounit) {

      /*           A is non-unit triangular. */

      /*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
      /*           Initially, G(0) = max{x(i), i=1,...,n}. */

      grow = 1. / max(xbnd,smlnum);
      xbnd = grow;
      i__1 = jlast;
      i__2 = jinc;
      for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

        /*              Exit the loop if the growth factor is too small. */

        if (grow <= smlnum) {
          goto L50;
        }

        /*              M(j) = G(j-1) / abs(A(j,j)) */

        tjj = (d__1 = ab[maind + j * ab_dim1], abs(d__1));
        /* Computing MIN */
        d__1 = xbnd, d__2 = min(1.,tjj) * grow;
        xbnd = min(d__1,d__2);
        if (tjj + cnorm[j] >= smlnum) {

          /*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

          grow *= tjj / (tjj + cnorm[j]);
        } else {

          /*                 G(j) could overflow, set GROW to 0. */

          grow = 0.;
        }
        /* L30: */
      }
      grow = xbnd;
    } else {

      /*           A is unit triangular. */

      /*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

      /* Computing MIN */
      d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
      grow = min(d__1,d__2);
      i__2 = jlast;
      i__1 = jinc;
      for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

        /*              Exit the loop if the growth factor is too small. */

        if (grow <= smlnum) {
          goto L50;
        }

        /*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

        grow *= 1. / (cnorm[j] + 1.);
        /* L40: */
      }
    }
  L50:

    ;
  } else {

    /*        Compute the growth in A**T * x = b. */

    if (upper) {
      jfirst = 1;
      jlast = *n;
      jinc = 1;
      maind = *kd + 1;
    } else {
      jfirst = *n;
      jlast = 1;
      jinc = -1;
      maind = 1;
    }

    if (tscal != 1.) {
      grow = 0.;
      goto L80;
    }

    if (nounit) {

      /*           A is non-unit triangular. */

      /*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
      /*           Initially, M(0) = max{x(i), i=1,...,n}. */

      grow = 1. / max(xbnd,smlnum);
      xbnd = grow;
      i__1 = jlast;
      i__2 = jinc;
      for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

        /*              Exit the loop if the growth factor is too small. */

        if (grow <= smlnum) {
          goto L80;
        }

        /*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

        xj = cnorm[j] + 1.;
        /* Computing MIN */
        d__1 = grow, d__2 = xbnd / xj;
        grow = min(d__1,d__2);

        /*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

        tjj = (d__1 = ab[maind + j * ab_dim1], abs(d__1));
        if (xj > tjj) {
          xbnd *= tjj / xj;
        }
        /* L60: */
      }
      grow = min(grow,xbnd);
    } else {

      /*           A is unit triangular. */

      /*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

      /* Computing MIN */
      d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
      grow = min(d__1,d__2);
      i__2 = jlast;
      i__1 = jinc;
      for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

        /*              Exit the loop if the growth factor is too small. */

        if (grow <= smlnum) {
          goto L80;
        }

        /*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

        xj = cnorm[j] + 1.;
        grow /= xj;
        /* L70: */
      }
    }
  L80:
    ;
  }

  if (grow * tscal > smlnum) {

    /*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
    /*        elements of X is not too small. */

    tbsv(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &x[1], &c__1, (
      ftnlen)1, (ftnlen)1, (ftnlen)1);
  } else {

    /*        Use a Level 1 BLAS solve, scaling intermediate results. */

    if (xmax > bignum) {

      /*           Scale X so that its components are less than or equal to */
      /*           BIGNUM in absolute value. */

      *scale = bignum / xmax;
      scal(n, scale, &x[1], &c__1);
      xmax = bignum;
    }

    if (notran) {

      /*           Solve A * x = b */

      i__1 = jlast;
      i__2 = jinc;
      for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

        /*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

        xj = (d__1 = x[j], abs(d__1));
        if (nounit) {
          tjjs = ab[maind + j * ab_dim1] * tscal;
        } else {
          tjjs = tscal;
          if (tscal == 1.) {
            goto L100;
          }
        }
        tjj = abs(tjjs);
        if (tjj > smlnum) {

          /*                    abs(A(j,j)) > SMLNUM: */

          if (tjj < 1.) {
            if (xj > tjj * bignum) {

              /*                          Scale x by 1/b(j). */

              rec = 1. / xj;
              scal(n, &rec, &x[1], &c__1);
              *scale *= rec;
              xmax *= rec;
            }
          }
          x[j] /= tjjs;
          xj = (d__1 = x[j], abs(d__1));
        } else if (tjj > 0.) {

          /*                    0 < abs(A(j,j)) <= SMLNUM: */

          if (xj > tjj * bignum) {

            /*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
            /*                       to avoid overflow when dividing by A(j,j). */

            rec = tjj * bignum / xj;
            if (cnorm[j] > 1.) {

              /*                          Scale by 1/CNORM(j) to avoid overflow when */
              /*                          multiplying x(j) times column j. */

              rec /= cnorm[j];
            }
            scal(n, &rec, &x[1], &c__1);
            *scale *= rec;
            xmax *= rec;
          }
          x[j] /= tjjs;
          xj = (d__1 = x[j], abs(d__1));
        } else {

          /*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
          /*                    scale = 0, and compute a solution to A*x = 0. */

          i__3 = *n;
          for (i__ = 1; i__ <= i__3; ++i__) {
            x[i__] = 0.;
            /* L90: */
          }
          x[j] = 1.;
          xj = 1.;
          *scale = 0.;
          xmax = 0.;
        }
      L100:

        /*              Scale x if necessary to avoid overflow when adding a */
        /*              multiple of column j of A. */

        if (xj > 1.) {
          rec = 1. / xj;
          if (cnorm[j] > (bignum - xmax) * rec) {

            /*                    Scale x by 1/(2*abs(x(j))). */

            rec *= .5;
            scal(n, &rec, &x[1], &c__1);
            *scale *= rec;
          }
        } else if (xj * cnorm[j] > bignum - xmax) {

          /*                 Scale x by 1/2. */

          scal(n, &c_b36, &x[1], &c__1);
          *scale *= .5;
        }

        if (upper) {
          if (j > 1) {

            /*                    Compute the update */
            /*                       x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) - */
            /*                                             x(j)* A(max(1,j-kd):j-1,j) */

            /* Computing MIN */
            i__3 = *kd, i__4 = j - 1;
            jlen = min(i__3,i__4);
            d__1 = -x[j] * tscal;
            axpy(&jlen, &d__1, &ab[*kd + 1 - jlen + j * ab_dim1]
              , &c__1, &x[j - jlen], &c__1);
            i__3 = j - 1;
            i__ = iamax(&i__3, &x[1], &c__1);
            xmax = (d__1 = x[i__], abs(d__1));
          }
        } else if (j < *n) {

          /*                 Compute the update */
          /*                    x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) - */
          /*                                          x(j) * A(j+1:min(j+kd,n),j) */

          /* Computing MIN */
          i__3 = *kd, i__4 = *n - j;
          jlen = min(i__3,i__4);
          if (jlen > 0) {
            d__1 = -x[j] * tscal;
            axpy(&jlen, &d__1, &ab[j * ab_dim1 + 2], &c__1, &x[
              j + 1], &c__1);
          }
          i__3 = *n - j;
          i__ = j + iamax(&i__3, &x[j + 1], &c__1);
          xmax = (d__1 = x[i__], abs(d__1));
        }
        /* L110: */
      }

    } else {

      /*           Solve A**T * x = b */

      i__2 = jlast;
      i__1 = jinc;
      for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

        /*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
        /*                                    k<>j */

        xj = (d__1 = x[j], abs(d__1));
        uscal = tscal;
        rec = 1. / max(xmax,1.);
        if (cnorm[j] > (bignum - xj) * rec) {

          /*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

          rec *= .5;
          if (nounit) {
            tjjs = ab[maind + j * ab_dim1] * tscal;
          } else {
            tjjs = tscal;
          }
          tjj = abs(tjjs);
          if (tjj > 1.) {

            /*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

            /* Computing MIN */
            d__1 = 1., d__2 = rec * tjj;
            rec = min(d__1,d__2);
            uscal /= tjjs;
          }
          if (rec < 1.) {
            scal(n, &rec, &x[1], &c__1);
            *scale *= rec;
            xmax *= rec;
          }
        }

        sumj = 0.;
        if (uscal == 1.) {

          /*                 If the scaling needed for A in the dot product is 1, */
          /*                 call DDOT to perform the dot product. */

          if (upper) {
            /* Computing MIN */
            i__3 = *kd, i__4 = j - 1;
            jlen = min(i__3,i__4);
            sumj = dot(&jlen, &ab[*kd + 1 - jlen + j * ab_dim1],
              &c__1, &x[j - jlen], &c__1);
          } else {
            /* Computing MIN */
            i__3 = *kd, i__4 = *n - j;
            jlen = min(i__3,i__4);
            if (jlen > 0) {
              sumj = dot(&jlen, &ab[j * ab_dim1 + 2], &c__1, &x[j + 1], &c__1);
            }
          }
        } else {

          /*                 Otherwise, use in-line code for the dot product. */

          if (upper) {
            /* Computing MIN */
            i__3 = *kd, i__4 = j - 1;
            jlen = min(i__3,i__4);
            i__3 = jlen;
            for (i__ = 1; i__ <= i__3; ++i__) {
              sumj += ab[*kd + i__ - jlen + j * ab_dim1] *
                uscal * x[j - jlen - 1 + i__];
              /* L120: */
            }
          } else {
            /* Computing MIN */
            i__3 = *kd, i__4 = *n - j;
            jlen = min(i__3,i__4);
            i__3 = jlen;
            for (i__ = 1; i__ <= i__3; ++i__) {
              sumj += ab[i__ + 1 + j * ab_dim1] * uscal * x[j +
                i__];
              /* L130: */
            }
          }
        }

        if (uscal == tscal) {

          /*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j) */
          /*                 was not used to scale the dotproduct. */

          x[j] -= sumj;
          xj = (d__1 = x[j], abs(d__1));
          if (nounit) {

            /*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

            tjjs = ab[maind + j * ab_dim1] * tscal;
          } else {
            tjjs = tscal;
            if (tscal == 1.) {
              goto L150;
            }
          }
          tjj = abs(tjjs);
          if (tjj > smlnum) {

            /*                       abs(A(j,j)) > SMLNUM: */

            if (tjj < 1.) {
              if (xj > tjj * bignum) {

                /*                             Scale X by 1/abs(x(j)). */

                rec = 1. / xj;
                scal(n, &rec, &x[1], &c__1);
                *scale *= rec;
                xmax *= rec;
              }
            }
            x[j] /= tjjs;
          } else if (tjj > 0.) {

            /*                       0 < abs(A(j,j)) <= SMLNUM: */

            if (xj > tjj * bignum) {

              /*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

              rec = tjj * bignum / xj;
              scal(n, &rec, &x[1], &c__1);
              *scale *= rec;
              xmax *= rec;
            }
            x[j] /= tjjs;
          } else {

            /*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
            /*                       scale = 0, and compute a solution to A**T*x = 0. */

            i__3 = *n;
            for (i__ = 1; i__ <= i__3; ++i__) {
              x[i__] = 0.;
              /* L140: */
            }
            x[j] = 1.;
            *scale = 0.;
            xmax = 0.;
          }
        L150:
          ;
        } else {

          /*                 Compute x(j) := x(j) / A(j,j) - sumj if the dot */
          /*                 product has already been divided by 1/A(j,j). */

          x[j] = x[j] / tjjs - sumj;
        }
        /* Computing MAX */
        d__2 = xmax, d__3 = (d__1 = x[j], abs(d__1));
        xmax = max(d__2,d__3);
        /* L160: */
      }
    }
    *scale /= tscal;
  }

  /*     Scale the column norms by 1/TSCAL for return. */

  if (tscal != 1.) {
    d__1 = 1. / tscal;
    scal(n, &d__1, &cnorm[1], &c__1);
  }

  return 0;

  /*     End of DLATBS */

} /* dlatbs_ */

//-------------------------------------------------------------------------------------------------

//  drscl - LAPACK auxiliary routine (version 3.8.0)
template<class T>
int rscl(integer *n, T *sa, T *sx, integer *incx)
{
  static T mul, cden;
  static logical done;
  static T cnum, cden1, cnum1;
  //extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *,
  //  integer *), dlabad_(doublereal *, doublereal *);
  //extern doublereal dlamch_(char *, ftnlen);
  static T bignum, smlnum;

  /* Parameter adjustments */
  --sx;

  /* Function Body */
  if (*n <= 0) {
    return 0;
  }

  /*     Get machine parameters */

  smlnum = lamch("S", (ftnlen)1);
  bignum = 1. / smlnum;
  labad(&smlnum, &bignum);

  /*     Initialize the denominator to SA and the numerator to 1. */

  cden = *sa;
  cnum = 1.;

L10:
  cden1 = cden * smlnum;
  cnum1 = cnum / bignum;
  if (abs(cden1) > abs(cnum) && cnum != 0.) {

    /*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM. */

    mul = smlnum;
    done = FALSE_;
    cden = cden1;
  } else if (abs(cnum1) > abs(cden)) {

    /*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM. */

    mul = bignum;
    done = FALSE_;
    cnum = cnum1;
  } else {

    /*        Multiply X by CNUM / CDEN and return. */

    mul = cnum / cden;
    done = TRUE_;
  }

  /*     Scale the vector X by MUL */

  scal(n, &mul, &sx[1], incx);

  if (! done) {
    goto L10;
  }

  return 0;

  /*     End of DRSCL */

} /* drscl_ */

}
