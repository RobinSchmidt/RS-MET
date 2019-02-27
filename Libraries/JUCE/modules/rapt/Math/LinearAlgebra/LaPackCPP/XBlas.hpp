#pragma once

namespace LaPackCPP {

// Contains code to set up the FPU control registers on x86 systems. The current double-double code
// requires that all arithmetic is done in double precision (as opposed to double-extended).

#ifdef x86
#ifdef _WIN32

#include <float.h>
#define FPU_FIX_DECL unsigned int __old_cw, __new_cw;
#define FPU_FIX_START \
  __old_cw = _control87(0, 0);  \
  __new_cw = _control87(0x00010000, 0x00030000);
#define FPU_FIX_STOP \
  _control87(*_old_cw, 0xFFFFFFFF);
#else  /* _WIN32 */

#if HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif

#ifndef _FPU_GETCW
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#endif

#ifndef _FPU_SETCW
#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#endif

#ifndef _FPU_EXTENDED
#define _FPU_EXTENDED 0x0300
#endif

#ifndef _FPU_DOUBLE
#define _FPU_DOUBLE 0x0200
#endif

#define FPU_FIX_DECL unsigned short __old_cw, __new_cw;
#define FPU_FIX_START \
  _FPU_GETCW(__old_cw); \
  __new_cw = (__old_cw & ~_FPU_EXTENDED) | _FPU_DOUBLE; \
  _FPU_SETCW(__new_cw); 
#define FPU_FIX_STOP \
  _FPU_SETCW(__old_cw);
#endif  /* else _WIN32 */

#else   /* x86 */
#define FPU_FIX_DECL
#define FPU_FIX_START
#define FPU_FIX_STOP
#endif  /* else x86 */


/* Split a double into 2 parts with at most 26 bits each. (2^27 + 1) */
//#define split 	(134217729.0)
static const double split = 134217729.0;

enum blas_order_type {
  blas_rowmajor = 101,
  blas_colmajor = 102 };

enum blas_trans_type {
  blas_no_trans   = 111,
  blas_trans      = 112,
  blas_conj_trans = 113 };

enum blas_prec_type {
  blas_prec_single     = 211,
  blas_prec_double     = 212,
  blas_prec_indigenous = 213,
  blas_prec_extra      = 214 };

// maybe wrap a namespace XBlasCPP around these functions, use capitalization (or not) of "BLAS" 
// in function names consistently

//-------------------------------------------------------------------------------------------------

/**
Purpose:
gbmv computes y = alpha * A * x + beta * y, where
A is a m x n banded matrix
x is a n x 1 vector
y is a m x 1 vector
alpha and beta are scalars


Arguments:
order  Order of AP; row or column major
trans  Transpose of AB; no trans, trans, or conjugate trans
m      Dimension of AB
n      Dimension of AB and the length of vector x
kl     Number of lower diagnols of AB
ku     Number of upper diagnols of AB
alpha
AB
lda    Leading dimension of AB lda >= ku + kl + 1
x
incx   The stride for vector x.
beta
y
incy   The stride for vector y.
prec   Specifies the internal precision to be used.
        = blas_prec_single: single precision.
        = blas_prec_double: double precision.
        = blas_prec_extra : anything at least 1.5 times as accurate
                            than double, and wider than 80-bits.
                            We use double-double in our implementation. */
template<class T>
void blas_gbmv_x(enum blas_order_type order, enum blas_trans_type trans, int m, int n, int kl,
  int ku, T alpha, const T *a, int lda, const T *x, int incx, T beta, T *y, int incy, 
  enum blas_prec_type prec);

//-------------------------------------------------------------------------------------------------

/**
Purpose:
This routines computes the matrix product:
  y  <-  alpha * op(A) * (x_head + x_tail) + beta * y
where
A is a m x n banded matrix
x is a n x 1 vector
y is a m x 1 vector
alpha and beta are scalars

Arguments:
order  Order of AB; row or column major*
trans  Transpose of AB; no trans, trans, or conjugate trans
m      Dimension of AB
n      Dimension of AB and the length of vector x and z
kl     Number of lower diagnols of AB
ku:    Number of upper diagnols of AB
alpha
AB
lda    Leading dimension of AB lda >= ku + kl + 1
head_x
tail_x
incx   The stride for vector x.
beta
y
incy   The stride for vector y.
prec   Specifies the internal precision to be used.
       = blas_prec_single: single precision.
       = blas_prec_double: double precision.
       = blas_prec_extra : anything at least 1.5 times as accurate
         than double, and wider than 80-bits. We use double-double in our implementation. */
template<class T>
void blas_gbmv2_x(enum blas_order_type order, enum blas_trans_type trans, int m, int n, int kl,
  int ku, T alpha, const T *a, int lda, const T *head_x, const T *tail_x,
  int incx, T beta, T *y, int incy, enum blas_prec_type prec);

//-------------------------------------------------------------------------------------------------

/**
Argument
rname     routine name
iflag     a negative value indicates that parameter number -IFLAG caused the error; a nonnegative
          value is an implementation-specific error code.
ival      the value of parameter number -IFLAG. */
void BLAS_error(const char *rname, int iflag, int ival, char *form, ...);

}