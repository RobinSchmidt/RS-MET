
namespace LaPackCPP {

//=================================================================================================
// BLAS level 1 routines



//-------------------------------------------------------------------------------------------------

// from dasum - Reference BLAS level1 routine (version 3.8.0) 
template<class T>
T asum(integer *n, T *dx, integer *incx)
{
  /* System generated locals */
  integer i__1, i__2;
  T ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

  /* Local variables */
  static integer i__, m, mp1;
  static T dtemp;
  static integer nincx;

  /* Parameter adjustments */
  --dx;

  /* Function Body */
  ret_val = 0.;
  dtemp = 0.;
  if (*n <= 0 || *incx <= 0) {
    return ret_val;
  }
  if (*incx == 1) {
    /*        code for increment equal to 1 */
    /*        clean-up loop */
    m = *n % 6;
    if (m != 0) {
      i__1 = m;
      for (i__ = 1; i__ <= i__1; ++i__) {
        dtemp += (d__1 = dx[i__], abs(d__1));
      }
      if (*n < 6) {
        ret_val = dtemp;
        return ret_val;
      }
    }
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 6) {
      dtemp = dtemp + (d__1 = dx[i__], abs(d__1)) + (d__2 = dx[i__ + 1],
        abs(d__2)) + (d__3 = dx[i__ + 2], abs(d__3)) + (d__4 = 
          dx[i__ + 3], abs(d__4)) + (d__5 = dx[i__ + 4], abs(d__5)) 
        + (d__6 = dx[i__ + 5], abs(d__6));
    }
  } else {

    /*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
      dtemp += (d__1 = dx[i__], abs(d__1));
    }
  }
  ret_val = dtemp;
  return ret_val;
} /* dasum_ */

  //-------------------------------------------------------------------------------------------------

// translated from daxpy, Reference BLAS level1 routine (version 3.8.0)
template<class T>
int axpy(long int* n, T *da, T *dx, long int *incx, T *dy, long int *incy)
{
  // System generated locals
  long int i__1;

  // Local variables
  static long int i__, m, ix, iy, mp1;

  // Parameter adjustments
  --dy;
  --dx;

  // Function Body
  if(*n <= 0) {
    return 0;
  }
  if(*da == 0.) {
    return 0;
  }
  if(*incx == 1 && *incy == 1) {
  // code for both increments equal to 1
    m = *n % 4;
    if(m != 0) {
      i__1 = m;
      for(i__ = 1; i__ <= i__1; ++i__) {  // clean-up loop
        dy[i__] += *da * dx[i__];
      }
    }
    if(*n < 4) {
      return 0;
    }
    mp1 = m + 1;
    i__1 = *n;
    for(i__ = mp1; i__ <= i__1; i__ += 4) {
      dy[i__] += *da * dx[i__];
      dy[i__ + 1] += *da * dx[i__ + 1];
      dy[i__ + 2] += *da * dx[i__ + 2];
      dy[i__ + 3] += *da * dx[i__ + 3];
    }
  }
  else {
  // code for unequal increments or equal increments not equal to 1
    ix = 1;
    iy = 1;
    if(*incx < 0) {
      ix = (-(*n) + 1) * *incx + 1;
    }
    if(*incy < 0) {
      iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__) {
      dy[iy] += *da * dx[ix];
      ix += *incx;
      iy += *incy;
    }
  }
  return 0;
}

//-------------------------------------------------------------------------------------------------

// translated from dcopy,  Reference BLAS level1 routine (version 3.8.0)
template<class T>
int copy(integer *n, T *dx, integer *incx, T *dy, integer *incy)
{
  // System generated locals 
  integer i__1;

  // Local variables 
  static integer i__, m, ix, iy, mp1;

  // Parameter adjustments 
  --dy;
  --dx;

  // Function Body 
  if (*n <= 0) {
    return 0;
  }
  if (*incx == 1 && *incy == 1) {

    // code for both increments equal to 1
    // clean-up loop 

    m = *n % 7;
    if (m != 0) {
      i__1 = m;
      for (i__ = 1; i__ <= i__1; ++i__) {
        dy[i__] = dx[i__];
      }
      if (*n < 7) {
        return 0;
      }
    }
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
      dy[i__] = dx[i__];
      dy[i__ + 1] = dx[i__ + 1];
      dy[i__ + 2] = dx[i__ + 2];
      dy[i__ + 3] = dx[i__ + 3];
      dy[i__ + 4] = dx[i__ + 4];
      dy[i__ + 5] = dx[i__ + 5];
      dy[i__ + 6] = dx[i__ + 6];
    }
  } else {

    // code for unequal increments or equal increments
    // not equal to 1 

    ix = 1;
    iy = 1;
    if (*incx < 0) {
      ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
      iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      dy[iy] = dx[ix];
      ix += *incx;
      iy += *incy;
    }
  }
  return 0;
} // copy


//-------------------------------------------------------------------------------------------------

// ddot -- Reference BLAS level1 routine (version 3.8.0) -- 
template<class T>
T dot(integer *n, T *dx, integer *incx, T *dy, integer *incy)
{
  /* System generated locals */
  integer i__1;
  T ret_val;

  /* Local variables */
  static integer i__, m, ix, iy, mp1;
  static T dtemp;

  /* Parameter adjustments */
  --dy;
  --dx;

  /* Function Body */
  ret_val = 0.;
  dtemp = 0.;
  if (*n <= 0) {
    return ret_val;
  }
  if (*incx == 1 && *incy == 1) {

    /*        code for both increments equal to 1 */


    /*        clean-up loop */

    m = *n % 5;
    if (m != 0) {
      i__1 = m;
      for (i__ = 1; i__ <= i__1; ++i__) {
        dtemp += dx[i__] * dy[i__];
      }
      if (*n < 5) {
        ret_val = dtemp;
        return ret_val;
      }
    }
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
      dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + 
        dx[i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + 
        dx[i__ + 4] * dy[i__ + 4];
    }
  } else {

    /*        code for unequal increments or equal increments */
    /*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
      ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
      iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      dtemp += dx[ix] * dy[iy];
      ix += *incx;
      iy += *incy;
    }
  }
  ret_val = dtemp;
  return ret_val;
} /* ddot_ */

//-------------------------------------------------------------------------------------------------

// translated from idamax, Reference BLAS level1 routine (version 3.8.0)
template<class T>
integer iamax(integer *n, T *dx, integer *incx)
{
  /* System generated locals */
  integer ret_val, i__1;
  T d__1;

  /* Local variables */
  static integer i__, ix;
  static T dmax__;


  /* Parameter adjustments */
  --dx;

  /* Function Body */
  ret_val = 0;
  if (*n < 1 || *incx <= 0) {
    return ret_val;
  }
  ret_val = 1;
  if (*n == 1) {
    return ret_val;
  }
  if (*incx == 1) {

    /*        code for increment equal to 1 */

    dmax__ = abs(dx[1]);
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
      if ((d__1 = dx[i__], abs(d__1)) > dmax__) {
        ret_val = i__;
        dmax__ = (d__1 = dx[i__], abs(d__1));
      }
    }
  } else {

    /*        code for increment not equal to 1 */

    ix = 1;
    dmax__ = abs(dx[1]);
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
      if ((d__1 = dx[ix], abs(d__1)) > dmax__) {
        ret_val = i__;
        dmax__ = (d__1 = dx[ix], abs(d__1));
      }
      ix += *incx;
    }
  }
  return ret_val;
} /* idamax_ */

//-------------------------------------------------------------------------------------------------

// translated from lsame, Reference BLAS level1 routine (version 3.1)
logical lsame(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len)
{
  // System generated locals
  logical ret_val;

  // Local variables 
  static integer inta, intb, zcode;

  // Test if the characters are equal
  ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
  if(ret_val) {
    return ret_val;
  }

  // Now test for equivalence if both characters are alphabetic.
  zcode = 'Z';

  // Use 'Z' rather than 'A' so that ASCII can be detected on Prime 
  // machines, on which ICHAR returns a value with bit 8 set. 
  // ICHAR('A') on Prime machines returns 193 which is the same as 
  // ICHAR('A') on an EBCDIC machine.

  inta = *(unsigned char *)ca;
  intb = *(unsigned char *)cb;

  if(zcode == 90 || zcode == 122) {

    // ASCII is assumed - ZCODE is the ASCII code of either lower or upper case 'Z'. 
    if(inta >= 97 && inta <= 122) {
      inta += -32;
    }
    if(intb >= 97 && intb <= 122) {
      intb += -32;
    }

  }
  else if(zcode == 233 || zcode == 169) {

   // EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or 
   // upper case 'Z'.
    if(inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta
      >= 162 && inta <= 169) {
      inta += 64;
    }
    if(intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb
      >= 162 && intb <= 169) {
      intb += 64;
    }

  }
  else if(zcode == 218 || zcode == 250) {

   // ASCII is assumed, on Prime machines - ZCODE is the ASCII code
   // plus 128 of either lower or upper case 'Z'. 

    if(inta >= 225 && inta <= 250) {
      inta += -32;
    }
    if(intb >= 225 && intb <= 250) {
      intb += -32;
    }
  }
  ret_val = inta == intb;

  // End of LSAME

  return ret_val;
}

//-------------------------------------------------------------------------------------------------

// tranlated from dscal, Reference BLAS level1 routine (version 3.8.0)
template<class T>
int scal(integer *n, T *da, T *dx, integer *incx)
{
  /* System generated locals */
  integer i__1, i__2;

  /* Local variables */
  static integer i__, m, mp1, nincx;

  /* Parameter adjustments */
  --dx;

  /* Function Body */
  if (*n <= 0 || *incx <= 0) {
    return 0;
  }
  if (*incx == 1) {

    /*        code for increment equal to 1 */
    /*        clean-up loop */

    m = *n % 5;
    if (m != 0) {
      i__1 = m;
      for (i__ = 1; i__ <= i__1; ++i__) {
        dx[i__] = *da * dx[i__];
      }
      if (*n < 5) {
        return 0;
      }
    }
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
      dx[i__] = *da * dx[i__];
      dx[i__ + 1] = *da * dx[i__ + 1];
      dx[i__ + 2] = *da * dx[i__ + 2];
      dx[i__ + 3] = *da * dx[i__ + 3];
      dx[i__ + 4] = *da * dx[i__ + 4];
    }
  } else {

    /*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
      dx[i__] = *da * dx[i__];
    }
  }
  return 0;
} /* scal */

//-------------------------------------------------------------------------------------------------

// translated from dswap, Reference BLAS level1 routine (version 3.8.0)
template<class T>
int swap(integer *n, T *dx, integer *incx, T *dy, integer *incy)
{
  // System generated locals 
  integer i__1;

  // Local variables 
  static integer i__, m, ix, iy, mp1;
  static T dtemp;

  // Parameter adjustments 
  --dy;
  --dx;

  // Function Body 
  if (*n <= 0) {
    return 0;
  }
  if (*incx == 1 && *incy == 1) {

    // code for both increments equal to 1
    // clean-up loop 
    m = *n % 3;
    if (m != 0) {
      i__1 = m;
      for (i__ = 1; i__ <= i__1; ++i__) {
        dtemp = dx[i__];
        dx[i__] = dy[i__];
        dy[i__] = dtemp;
      }
      if (*n < 3) {
        return 0;
      }
    }
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 3) {
      dtemp = dx[i__];
      dx[i__] = dy[i__];
      dy[i__] = dtemp;
      dtemp = dx[i__ + 1];
      dx[i__ + 1] = dy[i__ + 1];
      dy[i__ + 1] = dtemp;
      dtemp = dx[i__ + 2];
      dx[i__ + 2] = dy[i__ + 2];
      dy[i__ + 2] = dtemp;
    }
  } else {

    // code for unequal increments or equal increments not equal to 1
    ix = 1;
    iy = 1;
    if (*incx < 0) {
      ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
      iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      dtemp = dx[ix];
      dx[ix] = dy[iy];
      dy[iy] = dtemp;
      ix += *incx;
      iy += *incy;
    }
  }
  return 0;
} // swap

//-------------------------------------------------------------------------------------------------

// todo: fix linker errors for s_wsfe, s_stop, len_trim__, do_fio - these seem to be functions from
// libF2C
// translated from xerbla, Reference BLAS level1 routine (version 3.7.0)
int xerbla(char *srname, integer *info, ftnlen srname_len)
{
  // Code of the function has been commented out by Robin Schmidt - at the moment, xerbla is just
  // a dummy function - i should probably set a debug-breakpoint here and re-implement it 
  // completely - it deals with some weird I/O functions from the f2c library which seems to be 
  // useless clutter...

  /*
  // Table of constant values
  static integer c__1 = 1;

  // Format strings 
  static char fmt_9999[] = "(\002 ** On entry to \002,a,\002 parameter num"
    "ber \002,i2,\002 had \002,\002an illegal value\002)";

  // Builtin functions 
  integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
  int s_stop(char *, ftnlen);

  // Local variables
  extern integer len_trim__(char *, ftnlen);
  // Note by Robin Schmidt:
  // this was declared as an intrinsic function in the original xerbla.f file, but it stumped the 
  // f2c translator, giving an error about an unknown intrinsic function, so i changed the 
  // "INTRINSIC" keyword to "EXTERNAL" - that allowed the translation, but i'm not sure, if it 
  // gives the correct behavior 

  // Fortran I/O blocks
  static cilist io___1 = { 0, 6, 0, fmt_9999, 0 };

  s_wsfe(&io___1);
  do_fio(&c__1, srname, len_trim__(srname, srname_len));
  do_fio(&c__1, (char *)&(*info), (ftnlen)sizeof(integer));
  e_wsfe();

  s_stop("", (ftnlen)0);

  // End of XERBLA
  */

  return 0;
} // xerbla

//=================================================================================================
// BLAS level 2 routines


//-------------------------------------------------------------------------------------------------

// translated from dger, Reference BLAS level2 routine (version 3.7.0)
template<class T>
int ger(integer *m, integer *n, T *alpha, T *x, integer *incx, T *y, integer *incy, T *a, 
  integer *lda)
{
  // System generated locals
  integer a_dim1, a_offset, i__1, i__2;

  // Local variables
  static integer i__, j, ix, jy, kx, info;
  static T temp;
  extern int xerbla(char *, integer *, ftnlen);

  // Parameter adjustments
  --x;
  --y;
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;

  /* Function Body */
  info = 0;
  if (*m < 0) {
    info = 1;
  } else if (*n < 0) {
    info = 2;
  } else if (*incx == 0) {
    info = 5;
  } else if (*incy == 0) {
    info = 7;
  } else if (*lda < max(1,*m)) {
    info = 9;
  }
  if (info != 0) {
    xerbla("DGER  ", &info, (ftnlen)6);
    return 0;
  }

  // Quick return if possible.

  if (*m == 0 || *n == 0 || *alpha == 0.) {
    return 0;
  }

  // Start the operations. In this version the elements of A are 
  // accessed sequentially with one pass through A.
  if (*incy > 0) {
    jy = 1;
  } else {
    jy = 1 - (*n - 1) * *incy;
  }
  if (*incx == 1) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      if (y[jy] != 0.) {
        temp = *alpha * y[jy];
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
          a[i__ + j * a_dim1] += x[i__] * temp;
          // L10: 
        }
      }
      jy += *incy;
      // L20:
    }
  } else {
    if (*incx > 0) {
      kx = 1;
    } else {
      kx = 1 - (*m - 1) * *incx;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      if (y[jy] != 0.) {
        temp = *alpha * y[jy];
        ix = kx;
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
          a[i__ + j * a_dim1] += x[ix] * temp;
          ix += *incx;
          // L30: 
        }
      }
      jy += *incy;
      // L40: 
    }
  }

  return 0;

  // End of DGER 

} // ger

//-------------------------------------------------------------------------------------------------

// translated from dgbmv, Reference BLAS level2 routine (version 3.7.0)
template<class T>
int gbmv(char *trans, integer *m, integer *n, integer *kl, integer *ku, T *alpha, T *a, 
  integer *lda, T *x, integer *incx, T *beta, T *y, integer *incy, ftnlen trans_len)
{
  // System generated locals
  integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;

  // Local variables
  static integer i__, j, k, ix, iy, jx, jy, kx, ky, kup1, info;
  static T temp;
  static integer lenx, leny;
  extern logical lsame(char *, char *, ftnlen, ftnlen);
  extern int xerbla(char *, integer *, ftnlen);

  // Parameter adjustments
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --x;
  --y;

  // Function Body
  info = 0;
  if (! lsame(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame(trans, "T", (
    ftnlen)1, (ftnlen)1) && ! lsame(trans, "C", (ftnlen)1, (ftnlen)1)
    ) {
    info = 1;
  } else if (*m < 0) {
    info = 2;
  } else if (*n < 0) {
    info = 3;
  } else if (*kl < 0) {
    info = 4;
  } else if (*ku < 0) {
    info = 5;
  } else if (*lda < *kl + *ku + 1) {
    info = 8;
  } else if (*incx == 0) {
    info = 10;
  } else if (*incy == 0) {
    info = 13;
  }
  if (info != 0) {
    xerbla("DGBMV ", &info, (ftnlen)6);
    return 0;
  }

  // Quick return if possible.
  if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
    return 0;
  }

  // Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  // up the start points in  X  and  Y.
  if (lsame(trans, "N", (ftnlen)1, (ftnlen)1)) {
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

  // Start the operations. In this version the elements of A are 
  // accessed sequentially with one pass through the band part of A. 
  // First form  y := beta*y.

  if (*beta != 1.) {
    if (*incy == 1) {
      if (*beta == 0.) {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[i__] = 0.;
          // L10: 
        }
      } else {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[i__] = *beta * y[i__];
          // L20: 
        }
      }
    } else {
      iy = ky;
      if (*beta == 0.) {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[iy] = 0.;
          iy += *incy;
          // L30: 
        }
      } else {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[iy] = *beta * y[iy];
          iy += *incy;
          // L40: 
        }
      }
    }
  }
  if (*alpha == 0.) {
    return 0;
  }
  kup1 = *ku + 1;
  if (lsame(trans, "N", (ftnlen)1, (ftnlen)1)) {

    // Form  y := alpha*A*x + y.
    jx = kx;
    if (*incy == 1) {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = *alpha * x[jx];
        k = kup1 - j;
        // Computing MAX 
        i__2 = 1, i__3 = j - *ku;
        // Computing MIN 
        i__5 = *m, i__6 = j + *kl;
        i__4 = min(i__5,i__6);
        for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
          y[i__] += temp * a[k + i__ + j * a_dim1];
          // L50: 
        }
        jx += *incx;
        // L60: 
      }
    } else {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = *alpha * x[jx];
        iy = ky;
        k = kup1 - j;
        // Computing MAX 
        i__4 = 1, i__2 = j - *ku;
        // Computing MIN 
        i__5 = *m, i__6 = j + *kl;
        i__3 = min(i__5,i__6);
        for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
          y[iy] += temp * a[k + i__ + j * a_dim1];
          iy += *incy;
          // L70: 
        }
        jx += *incx;
        if (j > *ku) {
          ky += *incy;
        }
        // L80: 
      }
    }
  } else {

    // Form  y := alpha*A**T*x + y.
    jy = ky;
    if (*incx == 1) {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = 0.;
        k = kup1 - j;
        // Computing MAX 
        i__3 = 1, i__4 = j - *ku;
        // Computing MIN 
        i__5 = *m, i__6 = j + *kl;
        i__2 = min(i__5,i__6);
        for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
          temp += a[k + i__ + j * a_dim1] * x[i__];
          // L90:
        }
        y[jy] += *alpha * temp;
        jy += *incy;
        // L100: 
      }
    } else {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = 0.;
        ix = kx;
        k = kup1 - j;
        // Computing MAX 
        i__2 = 1, i__3 = j - *ku;
        // Computing MIN 
        i__5 = *m, i__6 = j + *kl;
        i__4 = min(i__5,i__6);
        for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
          temp += a[k + i__ + j * a_dim1] * x[ix];
          ix += *incx;
          // L110:
        }
        y[jy] += *alpha * temp;
        jy += *incy;
        if (j > *ku) {
          kx += *incx;
        }
        // L120: 
      }
    }
  }

  return 0;

} // gbmv

//-------------------------------------------------------------------------------------------------

// translated from dgemv, Reference BLAS level2 routine (version 3.7.0)
template<class T>
int gemv(char *trans, integer *m, integer *n, T *alpha, T *a, integer *lda, T *x, integer *incx, 
  T *beta, T *y, integer *incy, ftnlen trans_len)
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;

  /* Local variables */
  static integer i__, j, ix, iy, jx, jy, kx, ky, info;
  static T temp;
  static integer lenx, leny;
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  //extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

  /*     Test the input parameters. */

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --x;
  --y;

  /* Function Body */
  info = 0;
  if (! lsame(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame(trans, "T", (
    ftnlen)1, (ftnlen)1) && ! lsame(trans, "C", (ftnlen)1, (ftnlen)1)
    ) {
    info = 1;
  } else if (*m < 0) {
    info = 2;
  } else if (*n < 0) {
    info = 3;
  } else if (*lda < max(1,*m)) {
    info = 6;
  } else if (*incx == 0) {
    info = 8;
  } else if (*incy == 0) {
    info = 11;
  }
  if (info != 0) {
    xerbla("DGEMV ", &info, (ftnlen)6);
    return 0;
  }

  /*     Quick return if possible. */

  if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
    return 0;
  }

  /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
  /*     up the start points in  X  and  Y. */

  if (lsame(trans, "N", (ftnlen)1, (ftnlen)1)) {
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

  /*     Start the operations. In this version the elements of A are */
  /*     accessed sequentially with one pass through A. */

  /*     First form  y := beta*y. */

  if (*beta != 1.) {
    if (*incy == 1) {
      if (*beta == 0.) {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[i__] = 0.;
          /* L10: */
        }
      } else {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[i__] = *beta * y[i__];
          /* L20: */
        }
      }
    } else {
      iy = ky;
      if (*beta == 0.) {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[iy] = 0.;
          iy += *incy;
          /* L30: */
        }
      } else {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[iy] = *beta * y[iy];
          iy += *incy;
          /* L40: */
        }
      }
    }
  }
  if (*alpha == 0.) {
    return 0;
  }
  if (lsame(trans, "N", (ftnlen)1, (ftnlen)1)) {

    /*        Form  y := alpha*A*x + y. */

    jx = kx;
    if (*incy == 1) {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = *alpha * x[jx];
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
          y[i__] += temp * a[i__ + j * a_dim1];
          /* L50: */
        }
        jx += *incx;
        /* L60: */
      }
    } else {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = *alpha * x[jx];
        iy = ky;
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
          y[iy] += temp * a[i__ + j * a_dim1];
          iy += *incy;
          /* L70: */
        }
        jx += *incx;
        /* L80: */
      }
    }
  } else {

    /*        Form  y := alpha*A**T*x + y. */

    jy = ky;
    if (*incx == 1) {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = 0.;
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
          temp += a[i__ + j * a_dim1] * x[i__];
          /* L90: */
        }
        y[jy] += *alpha * temp;
        jy += *incy;
        /* L100: */
      }
    } else {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = 0.;
        ix = kx;
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
          temp += a[i__ + j * a_dim1] * x[ix];
          ix += *incx;
          /* L110: */
        }
        y[jy] += *alpha * temp;
        jy += *incy;
        /* L120: */
      }
    }
  }

  return 0;

  /*     End of DGEMV . */

} /* dgemv_ */

//-------------------------------------------------------------------------------------------------

// translated from dtbsv, Reference BLAS level2 routine (version 3.7.0)
template<class T>
int tbsv(char *uplo, char *trans, char *diag, integer *n, integer *k, T *a, integer *lda, T *x, 
  integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

  /* Local variables */
  static integer i__, j, l, ix, jx, kx, info;
  static T temp;
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static integer kplus1;
  //extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
  static logical nounit;

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --x;

  /* Function Body */
  info = 0;
  if (! lsame(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame(uplo, "L", (
    ftnlen)1, (ftnlen)1)) {
    info = 1;
  } else if (! lsame(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame(trans, 
    "T", (ftnlen)1, (ftnlen)1) && ! lsame(trans, "C", (ftnlen)1, (
      ftnlen)1)) {
    info = 2;
  } else if (! lsame(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame(diag, 
    "N", (ftnlen)1, (ftnlen)1)) {
    info = 3;
  } else if (*n < 0) {
    info = 4;
  } else if (*k < 0) {
    info = 5;
  } else if (*lda < *k + 1) {
    info = 7;
  } else if (*incx == 0) {
    info = 9;
  }
  if (info != 0) {
    xerbla("DTBSV ", &info, (ftnlen)6);
    return 0;
  }

  /*     Quick return if possible. */

  if (*n == 0) {
    return 0;
  }

  nounit = lsame(diag, "N", (ftnlen)1, (ftnlen)1);

  /*     Set up the start point in X if the increment is not unity. This */
  /*     will be  ( N - 1 )*INCX  too small for descending loops. */

  if (*incx <= 0) {
    kx = 1 - (*n - 1) * *incx;
  } else if (*incx != 1) {
    kx = 1;
  }

  /*     Start the operations. In this version the elements of A are */
  /*     accessed by sequentially with one pass through A. */

  if (lsame(trans, "N", (ftnlen)1, (ftnlen)1)) {

    /*        Form  x := inv( A )*x. */

    if (lsame(uplo, "U", (ftnlen)1, (ftnlen)1)) {
      kplus1 = *k + 1;
      if (*incx == 1) {
        for (j = *n; j >= 1; --j) {
          if (x[j] != 0.) {
            l = kplus1 - j;
            if (nounit) {
              x[j] /= a[kplus1 + j * a_dim1];
            }
            temp = x[j];
            /* Computing MAX */
            i__2 = 1, i__3 = j - *k;
            i__1 = max(i__2,i__3);
            for (i__ = j - 1; i__ >= i__1; --i__) {
              x[i__] -= temp * a[l + i__ + j * a_dim1];
              /* L10: */
            }
          }
          /* L20: */
        }
      } else {
        kx += (*n - 1) * *incx;
        jx = kx;
        for (j = *n; j >= 1; --j) {
          kx -= *incx;
          if (x[jx] != 0.) {
            ix = kx;
            l = kplus1 - j;
            if (nounit) {
              x[jx] /= a[kplus1 + j * a_dim1];
            }
            temp = x[jx];
            /* Computing MAX */
            i__2 = 1, i__3 = j - *k;
            i__1 = max(i__2,i__3);
            for (i__ = j - 1; i__ >= i__1; --i__) {
              x[ix] -= temp * a[l + i__ + j * a_dim1];
              ix -= *incx;
              /* L30: */
            }
          }
          jx -= *incx;
          /* L40: */
        }
      }
    } else {
      if (*incx == 1) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          if (x[j] != 0.) {
            l = 1 - j;
            if (nounit) {
              x[j] /= a[j * a_dim1 + 1];
            }
            temp = x[j];
            /* Computing MIN */
            i__3 = *n, i__4 = j + *k;
            i__2 = min(i__3,i__4);
            for (i__ = j + 1; i__ <= i__2; ++i__) {
              x[i__] -= temp * a[l + i__ + j * a_dim1];
              /* L50: */
            }
          }
          /* L60: */
        }
      } else {
        jx = kx;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          kx += *incx;
          if (x[jx] != 0.) {
            ix = kx;
            l = 1 - j;
            if (nounit) {
              x[jx] /= a[j * a_dim1 + 1];
            }
            temp = x[jx];
            /* Computing MIN */
            i__3 = *n, i__4 = j + *k;
            i__2 = min(i__3,i__4);
            for (i__ = j + 1; i__ <= i__2; ++i__) {
              x[ix] -= temp * a[l + i__ + j * a_dim1];
              ix += *incx;
              /* L70: */
            }
          }
          jx += *incx;
          /* L80: */
        }
      }
    }
  } else {

    /*        Form  x := inv( A**T)*x. */

    if (lsame(uplo, "U", (ftnlen)1, (ftnlen)1)) {
      kplus1 = *k + 1;
      if (*incx == 1) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          temp = x[j];
          l = kplus1 - j;
          /* Computing MAX */
          i__2 = 1, i__3 = j - *k;
          i__4 = j - 1;
          for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
            temp -= a[l + i__ + j * a_dim1] * x[i__];
            /* L90: */
          }
          if (nounit) {
            temp /= a[kplus1 + j * a_dim1];
          }
          x[j] = temp;
          /* L100: */
        }
      } else {
        jx = kx;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          temp = x[jx];
          ix = kx;
          l = kplus1 - j;
          /* Computing MAX */
          i__4 = 1, i__2 = j - *k;
          i__3 = j - 1;
          for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
            temp -= a[l + i__ + j * a_dim1] * x[ix];
            ix += *incx;
            /* L110: */
          }
          if (nounit) {
            temp /= a[kplus1 + j * a_dim1];
          }
          x[jx] = temp;
          jx += *incx;
          if (j > *k) {
            kx += *incx;
          }
          /* L120: */
        }
      }
    } else {
      if (*incx == 1) {
        for (j = *n; j >= 1; --j) {
          temp = x[j];
          l = 1 - j;
          /* Computing MIN */
          i__1 = *n, i__3 = j + *k;
          i__4 = j + 1;
          for (i__ = min(i__1,i__3); i__ >= i__4; --i__) {
            temp -= a[l + i__ + j * a_dim1] * x[i__];
            /* L130: */
          }
          if (nounit) {
            temp /= a[j * a_dim1 + 1];
          }
          x[j] = temp;
          /* L140: */
        }
      } else {
        kx += (*n - 1) * *incx;
        jx = kx;
        for (j = *n; j >= 1; --j) {
          temp = x[jx];
          ix = kx;
          l = 1 - j;
          /* Computing MIN */
          i__4 = *n, i__1 = j + *k;
          i__3 = j + 1;
          for (i__ = min(i__4,i__1); i__ >= i__3; --i__) {
            temp -= a[l + i__ + j * a_dim1] * x[ix];
            ix -= *incx;
            /* L150: */
          }
          if (nounit) {
            temp /= a[j * a_dim1 + 1];
          }
          x[jx] = temp;
          jx -= *incx;
          if (*n - j >= *k) {
            kx -= *incx;
          }
          /* L160: */
        }
      }
    }
  }

  return 0;

  /*     End of DTBSV . */

} /* tbsv */

//=================================================================================================
// BLAS level 3 routines

// translated from dgemm, Reference BLAS level3 routine (version 3.7.0)
template<class T>
int gemm(char *transa, char *transb, integer *m, integer *n, integer *k, T *alpha, T *a, 
  integer *lda, T *b, integer *ldb, T *beta, T *c__, 
  integer *ldc, ftnlen transa_len, ftnlen transb_len)
{
  // System generated locals
  integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
    i__3;

  // Local variables
  static integer i__, j, l, info;
  static logical nota, notb;
  static T temp;
  static integer ncola;
  extern logical lsame(char *, char *, ftnlen, ftnlen);
  static integer nrowa, nrowb;
  extern int xerbla(char *, integer *, ftnlen);

  // Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not 
  // transposed and set  NROWA, NCOLA and  NROWB  as the number of rows 
  // and  columns of  A  and the  number of  rows  of  B  respectively. 

  // Parameter adjustments
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;
  c_dim1 = *ldc;
  c_offset = 1 + c_dim1;
  c__ -= c_offset;

  // Function Body
  nota = lsame(transa, "N", (ftnlen)1, (ftnlen)1);
  notb = lsame(transb, "N", (ftnlen)1, (ftnlen)1);
  if (nota) {
    nrowa = *m;
    ncola = *k;
  } else {
    nrowa = *k;
    ncola = *m;
  }
  if (notb) {
    nrowb = *k;
  } else {
    nrowb = *n;
  }

  // Test the input parameters.
  info = 0;
  if (! nota && ! lsame(transa, "C", (ftnlen)1, (ftnlen)1) && ! lsame(
    transa, "T", (ftnlen)1, (ftnlen)1)) {
    info = 1;
  } else if (! notb && ! lsame(transb, "C", (ftnlen)1, (ftnlen)1) && ! 
    lsame(transb, "T", (ftnlen)1, (ftnlen)1)) {
    info = 2;
  } else if (*m < 0) {
    info = 3;
  } else if (*n < 0) {
    info = 4;
  } else if (*k < 0) {
    info = 5;
  } else if (*lda < max(1,nrowa)) {
    info = 8;
  } else if (*ldb < max(1,nrowb)) {
    info = 10;
  } else if (*ldc < max(1,*m)) {
    info = 13;
  }
  if (info != 0) {
    xerbla("DGEMM ", &info, (ftnlen)6);
    return 0;
  }

  // Quick return if possible.

  if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
    return 0;
  }

  // And if  alpha.eq.zero.
  if (*alpha == 0.) {
    if (*beta == 0.) {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
          c__[i__ + j * c_dim1] = 0.;
          // L10: 
        }
        // L20: 
      }
    } else {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
          c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
          // L30: 
        }
        // L40: 
      }
    }
    return 0;
  }

  // Start the operations.

  if (notb) {
    if (nota) {

      // Form  C := alpha*A*B + beta*C. 

      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        if (*beta == 0.) {
          i__2 = *m;
          for (i__ = 1; i__ <= i__2; ++i__) {
            c__[i__ + j * c_dim1] = 0.;
            // L50: 
          }
        } else if (*beta != 1.) {
          i__2 = *m;
          for (i__ = 1; i__ <= i__2; ++i__) {
            c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
            // L60: 
          }
        }
        i__2 = *k;
        for (l = 1; l <= i__2; ++l) {
          temp = *alpha * b[l + j * b_dim1];
          i__3 = *m;
          for (i__ = 1; i__ <= i__3; ++i__) {
            c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
            // L70: 
          }
          // L80: 
        }
        // L90: 
      }
    } else {

      // Form  C := alpha*A**T*B + beta*C

      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
          temp = 0.;
          i__3 = *k;
          for (l = 1; l <= i__3; ++l) {
            temp += a[l + i__ * a_dim1] * b[l + j * b_dim1];
            // L100: 
          }
          if (*beta == 0.) {
            c__[i__ + j * c_dim1] = *alpha * temp;
          } else {
            c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
              i__ + j * c_dim1];
          }
          // L110: 
        }
        // L120: 
      }
    }
  } else {
    if (nota) {

      // Form  C := alpha*A*B**T + beta*C 

      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        if (*beta == 0.) {
          i__2 = *m;
          for (i__ = 1; i__ <= i__2; ++i__) {
            c__[i__ + j * c_dim1] = 0.;
            // L130: 
          }
        } else if (*beta != 1.) {
          i__2 = *m;
          for (i__ = 1; i__ <= i__2; ++i__) {
            c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
            // L140: 
          }
        }
        i__2 = *k;
        for (l = 1; l <= i__2; ++l) {
          temp = *alpha * b[j + l * b_dim1];
          i__3 = *m;
          for (i__ = 1; i__ <= i__3; ++i__) {
            c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
            // L150: 
          }
          // L160: 
        }
        // L170: 
      }
    } else {

      // Form  C := alpha*A**T*B**T + beta*C

      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
          temp = 0.;
          i__3 = *k;
          for (l = 1; l <= i__3; ++l) {
            temp += a[l + i__ * a_dim1] * b[j + l * b_dim1];
            // L180: 
          }
          if (*beta == 0.) {
            c__[i__ + j * c_dim1] = *alpha * temp;
          } else {
            c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
              i__ + j * c_dim1];
          }
          // L190: 
        }
        // L200: 
      }
    }
  }

  return 0;

} // gemm

//-------------------------------------------------------------------------------------------------

// translated from dtrsm, Reference BLAS level3 routine (version 3.7.0) */
template<class T>
int trsm(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, T *alpha, 
  T *a, integer * lda, T *b, integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, 
  ftnlen diag_len)
{
  // System generated locals 
  integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

  // Local variables
  static integer i__, j, k, info;
  static T temp;
  static logical lside;
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static integer nrowa;
  static logical upper;
  //extern int xerbla_(char *, integer *, ftnlen);
  static logical nounit;

  // Parameter adjustments
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;

  // Function Body 
  lside = lsame(side, "L", (ftnlen)1, (ftnlen)1);
  if (lside) {
    nrowa = *m;
  } else {
    nrowa = *n;
  }
  nounit = lsame(diag, "N", (ftnlen)1, (ftnlen)1);
  upper = lsame(uplo, "U", (ftnlen)1, (ftnlen)1);

  info = 0;
  if (! lside && ! lsame(side, "R", (ftnlen)1, (ftnlen)1)) {
    info = 1;
  } else if (! upper && ! lsame(uplo, "L", (ftnlen)1, (ftnlen)1)) {
    info = 2;
  } else if (! lsame(transa, "N", (ftnlen)1, (ftnlen)1) && ! lsame(transa,
    "T", (ftnlen)1, (ftnlen)1) && ! lsame(transa, "C", (ftnlen)1, (
      ftnlen)1)) {
    info = 3;
  } else if (! lsame(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame(diag, 
    "N", (ftnlen)1, (ftnlen)1)) {
    info = 4;
  } else if (*m < 0) {
    info = 5;
  } else if (*n < 0) {
    info = 6;
  } else if (*lda < max(1,nrowa)) {
    info = 9;
  } else if (*ldb < max(1,*m)) {
    info = 11;
  }
  if (info != 0) {
    xerbla("DTRSM ", &info, (ftnlen)6);
    return 0;
  }

  // Quick return if possible

  if (*m == 0 || *n == 0) {
    return 0;
  }

  // And when  alpha.eq.zero.

  if (*alpha == 0.) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      i__2 = *m;
      for (i__ = 1; i__ <= i__2; ++i__) {
        b[i__ + j * b_dim1] = 0.;
        // L10: 
      }
      // L20: 
    }
    return 0;
  }

  // Start the operations.

  if (lside) {
    if (lsame(transa, "N", (ftnlen)1, (ftnlen)1)) {

      // Form  B := alpha*inv( A )*B.

      if (upper) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          if (*alpha != 1.) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
              b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
                ;
              // L30:
            }
          }
          for (k = *m; k >= 1; --k) {
            if (b[k + j * b_dim1] != 0.) {
              if (nounit) {
                b[k + j * b_dim1] /= a[k + k * a_dim1];
              }
              i__2 = k - 1;
              for (i__ = 1; i__ <= i__2; ++i__) {
                b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
                  i__ + k * a_dim1];
                // L40: 
              }
            }
            // L50:
          }
          // L60:
        }
      } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          if (*alpha != 1.) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
              b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
                ;
              // L70: 
            }
          }
          i__2 = *m;
          for (k = 1; k <= i__2; ++k) {
            if (b[k + j * b_dim1] != 0.) {
              if (nounit) {
                b[k + j * b_dim1] /= a[k + k * a_dim1];
              }
              i__3 = *m;
              for (i__ = k + 1; i__ <= i__3; ++i__) {
                b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
                  i__ + k * a_dim1];
                // L80: 
              }
            }
            // L90: 
          }
          // L100: 
        }
      }
    } else {

      // Form  B := alpha*inv( A**T )*B.

      if (upper) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          i__2 = *m;
          for (i__ = 1; i__ <= i__2; ++i__) {
            temp = *alpha * b[i__ + j * b_dim1];
            i__3 = i__ - 1;
            for (k = 1; k <= i__3; ++k) {
              temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
              // L110: 
            }
            if (nounit) {
              temp /= a[i__ + i__ * a_dim1];
            }
            b[i__ + j * b_dim1] = temp;
            // L120: 
          }
          // L130: 
        }
      } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          for (i__ = *m; i__ >= 1; --i__) {
            temp = *alpha * b[i__ + j * b_dim1];
            i__2 = *m;
            for (k = i__ + 1; k <= i__2; ++k) {
              temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
              // L140: 
            }
            if (nounit) {
              temp /= a[i__ + i__ * a_dim1];
            }
            b[i__ + j * b_dim1] = temp;
            // L150: 
          }
          // L160: 
        }
      }
    }
  } else {
    if (lsame(transa, "N", (ftnlen)1, (ftnlen)1)) {

      // Form  B := alpha*B*inv( A ). 

      if (upper) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          if (*alpha != 1.) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
              b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
                ;
              // L170: 
            }
          }
          i__2 = j - 1;
          for (k = 1; k <= i__2; ++k) {
            if (a[k + j * a_dim1] != 0.) {
              i__3 = *m;
              for (i__ = 1; i__ <= i__3; ++i__) {
                b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
                  i__ + k * b_dim1];
                // L180: 
              }
            }
            // L190: 
          }
          if (nounit) {
            temp = 1. / a[j + j * a_dim1];
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
              b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
              // L200: 
            }
          }
          // L210: 
        }
      } else {
        for (j = *n; j >= 1; --j) {
          if (*alpha != 1.) {
            i__1 = *m;
            for (i__ = 1; i__ <= i__1; ++i__) {
              b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
                ;
              // L220: 
            }
          }
          i__1 = *n;
          for (k = j + 1; k <= i__1; ++k) {
            if (a[k + j * a_dim1] != 0.) {
              i__2 = *m;
              for (i__ = 1; i__ <= i__2; ++i__) {
                b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
                  i__ + k * b_dim1];
                // L230: 
              }
            }
            // L240: 
          }
          if (nounit) {
            temp = 1. / a[j + j * a_dim1];
            i__1 = *m;
            for (i__ = 1; i__ <= i__1; ++i__) {
              b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
              // L250: 
            }
          }
          // L260:
        }
      }
    } else {

      // Form  B := alpha*B*inv( A**T ).

      if (upper) {
        for (k = *n; k >= 1; --k) {
          if (nounit) {
            temp = 1. / a[k + k * a_dim1];
            i__1 = *m;
            for (i__ = 1; i__ <= i__1; ++i__) {
              b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
              // L270: 
            }
          }
          i__1 = k - 1;
          for (j = 1; j <= i__1; ++j) {
            if (a[j + k * a_dim1] != 0.) {
              temp = a[j + k * a_dim1];
              i__2 = *m;
              for (i__ = 1; i__ <= i__2; ++i__) {
                b[i__ + j * b_dim1] -= temp * b[i__ + k * 
                  b_dim1];
                // L280: 
              }
            }
            // L290: 
          }
          if (*alpha != 1.) {
            i__1 = *m;
            for (i__ = 1; i__ <= i__1; ++i__) {
              b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
                ;
              // L300: 
            }
          }
          // L310: 
        }
      } else {
        i__1 = *n;
        for (k = 1; k <= i__1; ++k) {
          if (nounit) {
            temp = 1. / a[k + k * a_dim1];
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
              b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
              // L320:
            }
          }
          i__2 = *n;
          for (j = k + 1; j <= i__2; ++j) {
            if (a[j + k * a_dim1] != 0.) {
              temp = a[j + k * a_dim1];
              i__3 = *m;
              for (i__ = 1; i__ <= i__3; ++i__) {
                b[i__ + j * b_dim1] -= temp * b[i__ + k * 
                  b_dim1];
                // L330: 
              }
            }
            // L340: 
          }
          if (*alpha != 1.) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
              b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
                ;
              // L350: 
            }
          }
          // L360: 
        }
      }
    }
  }

  return 0;

} // trsm 




//=================================================================================================

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
  ftnlen uplo_len, ftnlen transa_len,  ftnlen diag_len);




} // namespace BlasCPP