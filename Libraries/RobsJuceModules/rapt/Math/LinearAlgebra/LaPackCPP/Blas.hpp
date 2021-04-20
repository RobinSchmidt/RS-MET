#pragma once

namespace LaPackCPP {


// How can we make sure, that client code can use a drop-in replacement for these (unoptimized)
// reference implementations of the BLAS routines? Not allowing that would go totally against the
// design idea of LAPACK (which is to rely on lower-level linear algebra routines that can be
// optimized for particular target platform and then dropped in). Maybe the LAPACK routines call
// a general template function like axpy (not daxpy, saxpy, caxpy or zaxpy)...and then it must be
// somehow dispatched at compile-time which of the 4 (d,s,c,z) is chosen depending on the template
// parameter for axpy. The design goal of this C++ translation is to NOT have separate routines for
// different datatypes and instead just have one axpy function...maybe for some things that are 
// specific to complex numbers, that may not work out completely - we'll see - but for many of the 
// BLAS and LAPACK routines, it should be possible to write them in a type independent way as 
// templates. Maybe an explicit specialization of axpy for double should call daxpy and that daxpy
// may then be defined somewhere else. If we don't want to use such an explicit specialization, we 
// just request an explicit instantiation of axpy from the compiler - this can be done by letting
// the user decide which source files are compiled with the application - whether it's the ones 
// with explicit instantiations or the ones with explicit specializations...something along those
// lines...
// Maybe i should implement my own, very naive, BLAS routines (no loop unrolling, etc.). That stuff
// might better be left to the compiler anyway - maybe unrolling by 4 (as this reference 
// implementation does) is not optimal and an optimizing compiler would unroll by some other number 
// (8,16,..) instead? So maybe a more naive BLAS implementation, in addition to be more readable, 
// could indeed perfom better? -> try it! With such a "NaiveBlas" in place, we could also try the 
// replacement mechanism, once it's implemented, and swap between the reference implementation 
// maybe "RefBlas" and make some comparison benchmarks. ..also try https://www.netlib.org/atlas/

//=================================================================================================

/** \name BLAS level 1 routines (operations involving scalars and vectors) */

/** Takes the sum of the absolute values.
N:    number of elements in input vector(s)
DX:   array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
INCX: storage spacing between elements of DX  */
template<class T>
T asum(integer *n, T *dx, integer *incx);

/** Computes constant times a vector plus a vector: y = a*x + y
N:    Number of elements in input vector(s)
a:    On entry, a specifies the scalar alpha.
x:    array, dimension ( 1 + ( N - 1 )*abs( incX ) ), input
incX: storage spacing between elements of x
y:    array, dimension ( 1 + ( N - 1 )*abs( incY ) ), input/output
incY: storage spacing between elements of y  */
template<class T>
int axpy(long int* N, T* a, T* x, long int* incX, T* y, long int* incY);

/* Copies a vector, x, to a vector, y. uses unrolled loops for increments equal to 1.

Arguments:
N:    number of elements in input vector(s)
DX:   array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
INCX: storage spacing between elements of DX
DY:   array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
INCY: storage spacing between elements of DY    */
template<class T>
int copy(integer *n, T *dx, integer *incx, T *dy, integer *incy);

/** Forms the dot product of two vectors. uses unrolled loops for increments equal to one.

Arguments:
N:    number of elements in input vector(s)
DX:   array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
INCX: storage spacing between elements of DX
DY:   array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
INCY: storage spacing between elements of DY  */
template<class T>
T dot(integer *n, T *dx, integer *incx, T *dy, integer *incy);


/** iamax finds the index of the first element having maximum absolute value.

Arguments:
N:    number of elements in input vector(s)
DX:   array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
INCX: storage spacing between elements of DX  */
template<class T>
integer iamax(integer *n, T *dx, integer *incx);


/** LSAME returns .TRUE. if CA is the same letter as CB regardless of case. CA and CB specify the 
single characters to be compared. */
logical lsame(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len);
//inline logical lsame_(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len) 
//{ 
//  return lsame(ca, cb, ca_len, cb_len); 
//} // temporary - so code can use both lsame and lsame_ (less editing work)

/** scal scales a vector by a constant. uses unrolled loops for increment equal to 1.

Arguments:
N:    number of elements in input vector(s)
DA:   On entry, DA specifies the scalar alpha.
DX:   array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
INCX: storage spacing between elements of DX */
template<class T>
int scal(integer *n, T *da, T *dx, integer *incx);


/** Interchanges two vectors. uses unrolled loops for increments equal to 1. 

Arguments:
N:    number of elements in input vector(s)
DX:   array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
INCX: storage spacing between elements of DX
DY:   array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
INCY: storage spacing between elements of DY   */
template<class T>
int swap(integer *n, T *dx, integer *incx, T *dy, integer *incy);


/** xerbla is an error handler

Purpose:
XERBLA is an error handler for the LAPACK routines. It is called by an LAPACK routine if an input 
parameter has an invalid value. A message is printed and execution stops. Installers may consider 
modifying the STOP statement in order to call system-specific exception-handling facilities. 

Arguments:
SRNAME: The name of the routine which called XERBLA.
INFO:   The position of the invalid parameter in the parameter list of the calling routine. */
int xerbla(char *srname, integer *info, ftnlen srname_len);

int xerbla(const char *srname, integer *info, ftnlen srname_len);
// Needed, so we can call it like xerbla("blablabla") without getting compiler warnings from clang
// that ISO C++11 doesn't allow conversion from string literal to char*

//=================================================================================================

/** \name BLAS level 2 routines (operations involving matrices and vectors) */


/**
Purpose: tbsv solves one of the systems of equations
   A*x = b,   or   A**T*x = b,
where b and x are n element vectors and A is an n by n unit, or non-unit, upper or lower 
triangular band matrix, with ( k + 1 ) diagonals. No test for singularity or near-singularity is 
included in this routine. Such tests must be performed before calling this routine.

Arguments:
UPLO:  On entry, UPLO specifies whether the matrix is an upper or lower triangular matrix as follows:
       UPLO = 'U' or 'u'   A is an upper triangular matrix.
       UPLO = 'L' or 'l'   A is a lower triangular matrix.
TRANS: On entry, TRANS specifies the equations to be solved as follows:
       TRANS = 'N' or 'n'   A*x = b.
       TRANS = 'T' or 't'   A**T*x = b.
       TRANS = 'C' or 'c'   A**T*x = b.
DIAG:  On entry, DIAG specifies whether or not A is unit triangular as follows:
       DIAG = 'U' or 'u'   A is assumed to be unit triangular.
       DIAG = 'N' or 'n'   A is not assumed to be unit triangular.
N:     On entry, N specifies the order of the matrix A. N must be at least zero.
K:     On entry with UPLO = 'U' or 'u', K specifies the number of super-diagonals of the matrix A.
       On entry with UPLO = 'L' or 'l', K specifies the number of sub-diagonals of the matrix A.
       K must satisfy  0 .le. K.
A:     array, dimension ( LDA, N ). Before entry with UPLO = 'U' or 'u', the leading ( k + 1 ) by 
       n part of the array A must contain the upper triangular band part of the matrix of 
       coefficients, supplied column by column, with the leading diagonal of the matrix in row
       ( k + 1 ) of the array, the first super-diagonal starting at position 2 in row k, and so on.
       The top left k by k triangle of the array A is not referenced. The following program segment 
       will transfer an upper triangular band matrix from conventional full matrix storage to band 
       storage:
           DO 20, J = 1, N
               M = K + 1 - J
              DO 10, I = MAX( 1, J - K ), J
                 A( M + I, J ) = matrix( I, J )
        10    CONTINUE
        20 CONTINUE
       Before entry with UPLO = 'L' or 'l', the leading ( k + 1 ) by n part of the array A must 
       contain the lower triangular band part of the matrix of coefficients, supplied column by 
       column, with the leading diagonal of the matrix in row 1 of the array, the first 
       sub-diagonal starting at position 1 in row 2, and so on. The bottom right k by k triangle of
       the array A is not referenced. The following program segment will transfer a lower triangular
       band matrix from conventional full matrix storage  to band storage:
           DO 20, J = 1, N
              M = 1 - J
              DO 10, I = J, MIN( N, J + K )
                 A( M + I, J ) = matrix( I, J )
        10    CONTINUE
        20 CONTINUE
       Note that when DIAG = 'U' or 'u' the elements of the array A corresponding to the diagonal 
       elements of the matrix are not referenced, but are assumed to be unity.
LDA:   On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. 
       LDA must be at least ( k + 1 ).
X:     array, dimension at least ( 1 + ( n - 1 )*abs( INCX ) ). Before entry, the incremented array
       X must contain the n element right-hand side vector b. On exit, X is overwritten with the 
       solution vector x.
INCX:  On entry, INCX specifies the increment for the elements of X. INCX must not be zero.  */
template<class T>
int tbsv(char *uplo, char *trans, char *diag, integer *n, integer *k, T *a, integer *lda,
  T *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);
// move down (should be in alphabetical order)

//-------------------------------------------------------------------------------------------------

/** ger updates matrix A := alpha*x*y**T + A

Purpose: 
ger performs the rank 1 operation 
A := alpha*x*y**T + A
where alpha is a scalar, x is an m element vector, y is an n element vector and A is an m by n 
matrix.

Arguments:
M:     On entry, M specifies the number of rows of the matrix A. M must be at least zero.
N:     On entry, N specifies the number of columns of the matrix A. N must be at least zero.
ALPHA: On entry, ALPHA specifies the scalar alpha.
X:     array, dimension at least ( 1 + ( m - 1 )*abs( INCX ) ). Before entry, the incremented array
       X must contain the m element vector x.
INCX:  On entry, INCX specifies the increment for the elements of X. INCX must not be zero.
Y:     array, dimension at least ( 1 + ( n - 1 )*abs( INCY ) ). Before entry, the incremented array 
       Y must contain the n element vector y.
INCY:  On entry, INCY specifies the increment for the elements of Y. INCY must not be zero.
A:     array, dimension ( LDA, N ). Before entry, the leading m by n part of the array A must 
       contain the matrix of coefficients. On exit, A is overwritten by the updated matrix.
LDA:   On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. 
       LDA must be at least max( 1, m ). */
template<class T>
int ger(integer *m, integer *n, T *alpha, T *x, integer *incx, T *y, integer *incy, T *a, 
  integer *lda);

//-------------------------------------------------------------------------------------------------

/** gbmv performs general banded matrix-vector multiplication

Purpose:
gbmv performs one of the matrix-vector operations
y := alpha*A*x    + beta*y,   or
y := alpha*A**T*x + beta*y
where alpha and beta are scalars, x and y are vectors and A is an m by n band matrix, with kl 
sub-diagonals and ku super-diagonals.

Arguments: 
TRANS: On entry, TRANS specifies the operation to be performed as follows: 
       TRANS = 'N' or 'n'   y := alpha*A*x + beta*y
       TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y. 
       TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y. 
M:     On entry, M specifies the number of rows of the matrix A. M must be at least zero. 
N:     On entry, N specifies the number of columns of the matrix A. N must be at least zero. 
KL:    On entry, KL specifies the number of sub-diagonals of the matrix A. KL must satisfy 
       0 .le. KL. 
KU:    On entry, KU specifies the number of super-diagonals of the matrix A. KU must satisfy  
       0 .le. KU. 
ALPHA: On entry, ALPHA specifies the scalar alpha. 
A:     array, dimension ( LDA, N ). Before entry, the leading ( kl + ku + 1 ) by n part of the 
       array A must contain the matrix of coefficients, supplied column by column, with the leading
       diagonal of the matrix in row ( ku + 1 ) of the array, the first super-diagonal starting at 
       position 2 in row ku, the first sub-diagonal starting at position 1 in row ( ku + 2 ), and 
       so on. Elements in the array A that do not correspond to elements in the band matrix (such 
       as the top left ku by ku triangle) are not referenced. The following program segment will 
       transfer a band matrix from conventional full matrix storage to band storage: 
          DO 20, J = 1, N 
             K = KU + 1 - J 
             DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL ) 
                A( K + I, J ) = matrix( I, J ) 
       10    CONTINUE 
       20 CONTINUE 
LDA:   On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. LDA 
       must be at least ( kl + ku + 1 ). 
X:     array, dimension at least ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' and at least
       ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. Before entry, the incremented array X must contain the 
       vector x. 
INCX:  On entry, INCX specifies the increment for the elements of X. INCX must not be zero. 
BETA:  On entry, BETA specifies the scalar beta. When BETA is supplied as zero then Y need not be 
       set on input. 
Y:     array, dimension at least ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' and at least 
       ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. Before entry, the incremented array Y must contain the 
       vector y. On exit, Y is overwritten by the updated vector y.   
INCY:  On entry, INCY specifies the increment for the elements of Y. INCY must not be zero. 

Robin's notes: 
trans_len is not used in the routine and the returned value is always zero (so, it has no meaning). */
template<class T>
int gbmv(char *trans, integer *m, integer *n, integer *kl, integer *ku, T *alpha, T *a, 
  integer *lda, T *x, integer *incx, T *beta, T *y, integer *incy, ftnlen trans_len);

//-------------------------------------------------------------------------------------------------

/**
Purpose: gemv  performs one of the matrix-vector operations
  y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
where alpha and beta are scalars, x and y are vectors and A is an m by n matrix.

Arguments:
TRANS: On entry, TRANS specifies the operation to be performed as follows:
       TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
       TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
       TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
M:     On entry, M specifies the number of rows of the matrix A. M must be at least zero.
N:     On entry, N specifies the number of columns of the matrix A. N must be at least zero.
ALPHA: On entry, ALPHA specifies the scalar alpha.
A:     array, dimension ( LDA, N ). Before entry, the leading m by n part of the array A must 
       contain the matrix of coefficients.
LDA:   On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. 
       LDA must be at least max( 1, m ).
X:     array, dimension at least ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' and at 
       least ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. Before entry, the incremented array X must 
       contain the vector x.
INCX:  On entry, INCX specifies the increment for the elements of X. INCX must not be zero.
BETA:  On entry, BETA specifies the scalar beta. When BETA is supplied as zero then Y need not be 
       set on input.
Y:     array, dimension at least ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' and at 
       least ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. Before entry with BETA non-zero, the 
       incremented array Y must contain the vector y. On exit, Y is overwritten by the updated 
       vector y.
INCY:  On entry, INCY specifies the increment for the elements of Y. INCY must not be zero.

Further Details:
The vector and matrix arguments are not referenced when N = 0, or M = 0 */
template<class T>
int gemv(char *trans, integer *m, integer *n, T *alpha, T *a, integer *lda, T *x, integer *incx, 
  T *beta, T *y, integer *incy, ftnlen trans_len);

//-------------------------------------------------------------------------------------------------











//=================================================================================================

/** \name BLAS level 3 routines (operations involving two matrices) */

/**

Purpose:
GEMM  performs one of the matrix-matrix operations
 C := alpha*op( A )*op( B ) + beta*C
where  op( X ) is one of
 op( X ) = X   or   op( X ) = X**T,
alpha and beta are scalars, and A, B and C are matrices, with op( A ) an m by k matrix, op( B )  
a k by n matrix and  C an m by n matrix. 

Arguments:
TRANSA: On entry, TRANSA specifies the form of op( A ) to be used in the matrix multiplication as 
        follows: 
        TRANSA = 'N' or 'n',  op( A ) = A. 
        TRANSA = 'T' or 't',  op( A ) = A**T. 
        TRANSA = 'C' or 'c',  op( A ) = A**T. 
TRANSB: On entry, TRANSB specifies the form of op( B ) to be used in the matrix multiplication as 
        follows
        TRANSB = 'N' or 'n',  op( B ) = B.
        TRANSB = 'T' or 't',  op( B ) = B**T.
        TRANSB = 'C' or 'c',  op( B ) = B**T.
M:      On entry,  M  specifies  the number  of rows  of the  matrix op( A )  and of the  matrix C.  
        M  must  be at least  zero.
N:      On entry,  N  specifies the number  of columns of the matrix op( B ) and the number of 
        columns of the matrix C. N must be at least zero.
K:      On entry,  K  specifies  the number of columns of the matrix op( A ) and the number of rows 
        of the matrix op( B ). K must be at least  zero.
ALPHA:  On entry, ALPHA specifies the scalar alpha.
A:      array, dimension ( LDA, ka ), where ka is k  when  TRANSA = 'N' or 'n',  and is m otherwise.
        Before entry with  TRANSA = 'N' or 'n',  the leading  m by k part of the array A must contain 
        the matrix A, otherwise the leading  k by m  part of the array  A  must contain  the matrix A.
LDA:    On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. 
        When  TRANSA = 'N' or 'n' then LDA must be at least  max( 1, m ), otherwise  LDA must be at
        least  max( 1, k ).
B:      array, dimension ( LDB, kb ), where kb is n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
        Before entry with  TRANSB = 'N' or 'n',  the leading  k by n part of the array  B  must contain
        the matrix B, otherwise the leading n by k part of the array B must contain  the  matrix B.
LDB:    On entry, LDB specifies the first dimension of B as declared in the calling (sub) program. 
        When  TRANSB = 'N' or 'n' then LDB must be at least  max( 1, k ), otherwise  LDB must be at 
        least  max( 1, n ).
BETA:   On entry,  BETA  specifies the scalar  beta.  When  BETA  is supplied as zero then C need 
        not be set on input.
C:      array, dimension ( LDC, N ) Before entry, the leading  m by n  part of the array  C must 
        contain the matrix  C,  except when  beta  is zero, in which case C need not be set on 
        entry. On exit, the array  C  is overwritten by the  m by n  matrix 
        ( alpha*op( A )*op( B ) + beta*C ).
LDC:    On entry, LDC specifies the first dimension of C as declared in the calling (sub) program.   
        LDC must be at least max( 1, m ). */
template<class T>
int gemm(char *transa, char *transb, integer *m, integer *n, integer *k, T *alpha, T *a, 
  integer *lda, T *b, integer *ldb, T *beta, T *c__, integer *ldc, ftnlen transa_len, 
  ftnlen transb_len);

//-------------------------------------------------------------------------------------------------

/**
Purpose:
trsm solves one of the matrix equations
 op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
where alpha is a scalar, X and B are m by n matrices, A is a unit, or non-unit, upper or lower 
triangular matrix  and  op( A )  is one  of
 op( A ) = A   or   op( A ) = A**T.
The matrix X is overwritten on B.

Arguments:
SIDE:   On entry, SIDE specifies whether op( A ) appears on the left or right of X as follows:
        IDE = 'L' or 'l'   op( A )*X = alpha*B.
        SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
UPLO:   On entry, UPLO specifies whether the matrix A is an upper or lower triangular matrix as 
        follows:
        UPLO = 'U' or 'u'   A is an upper triangular matrix.
        UPLO = 'L' or 'l'   A is a lower triangular matrix.
TRANSA: On entry, TRANSA specifies the form of op( A ) to be used in the matrix multiplication as 
        follows:
        TRANSA = 'N' or 'n'   op( A ) = A.
        TRANSA = 'T' or 't'   op( A ) = A**T.
        TRANSA = 'C' or 'c'   op( A ) = A**T.
DIAG:   On entry, DIAG specifies whether or not A is unit triangular as follows:
        DIAG = 'U' or 'u'   A is assumed to be unit triangular.
        DIAG = 'N' or 'n'   A is not assumed to be unit triangular.
M:      On entry, M specifies the number of rows of B. M must be at least zero.
N:      On entry, N specifies the number of columns of B.  N must be at least zero.
ALPHA:  On entry,  ALPHA specifies the scalar  alpha. When  alpha is zero then  A is not 
        referenced and  B need not be set before entry.
A:      array, dimension ( LDA, k ), where k is m when SIDE = 'L' or 'l' and k is n when 
        SIDE = 'R' or 'r'. Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k upper 
        triangular part of the array  A must contain the upper triangular matrix  and the strictly 
        lower triangular part of A is not referenced. Before entry  with  UPLO = 'L' or 'l',  the 
        leading  k by k lower triangular part of the array  A must contain the lower triangular 
        matrix  and the strictly upper triangular part of A is not referenced. Note that when  
        DIAG = 'U' or 'u',  the diagonal elements of A are not referenced either, but are assumed 
        to be  unity.
LDA:    On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. 
        When  SIDE = 'L' or 'l'  then LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
        then LDA must be at least max( 1, n ).
B:      array, dimension ( LDB, N ) Before entry,  the leading  m by n part of the array  B must
        contain  the  right-hand  side  matrix  B,  and  on exit  is overwritten by the solution 
        matrix  X.
LDB:    On entry, LDB specifies the first dimension of B as declared in the calling (sub) program.
        LDB must be at least max( 1, m ).  */
template<class T>
int trsm(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, T *alpha, 
  T *a, integer *lda, T *b, integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, 
  ftnlen diag_len);
// todo: enable for complex case

}