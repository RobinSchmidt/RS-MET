#ifndef RAPT_LINEARALGEBRA_H_INCLUDED
#define RAPT_LINEARALGEBRA_H_INCLUDED

/** Collection of functions for linear algebra such as solving systems of linear equations, matrix 
inversion, etc.

ToDo: get rid of rs-prefix in the function names (they are now wrapped into a class which already 
has the prefix)
*/

class rsLinearAlgebra
{

public:


  /** Solves the 2x2 system A*x = y where A is a 2x2 matrix and x, y are 2-dimensional vectors. The
  result is returned in x, the other parameters will not be modified. */
  template<class T>
  static void rsSolveLinearSystem2x2(const T A[2][2], T x[2], const T y[2]);
    // rename to solve2x2

  /** Solves the 3x3 system A*x = y where A is a 3x3 matrix and x, y are 3-dimensional vectors. The
  result is returned in x, the other parameters will not be modified. */
  template<class T>
  static void rsSolveLinearSystem3x3(const T A[3][3], T x[3], const T y[3]);
    // todo: check, if this is the right form of const-ness
    // rename to solve3x3

  /** Computes the first eigenvalue of the matrix [[a,b],[c,d]] where "first" means the one with 
  the minus-sign in the term under the square-root in the formula. For real eigenvalues, this is
  the smaller of the two. If a matrix with real coefficients has complex eigenvalues ... it 
  currently doesn't work (encounters a sqrt of a negative number...todo: maybe return the real part
  of the result) */
  template<class T>
  static T eigenvalue2x2_1(T a, T b, T c, T d);

  /** Computes the second eigenvalue of the matrix [[a,b],[c,d]]. */
  template<class T>
  static T eigenvalue2x2_2(T a, T b, T c, T d);

  /** Computes the first eigenvector of the matrix [[a,b],[c,d]] and stores the result in vx, vy. 
  It may optionally normalize the the result to unit length. */
  template<class T>
  static void eigenvector2x2_1(T a, T b, T c, T d, T* vx, T* vy, bool normalize = true);

  /** Computes the second eigenvector of the matrix [[a,b],[c,d]] */
  template<class T>
  static void eigenvector2x2_2(T a, T b, T c, T d, T* vx, T* vy, bool normalize = true);

  /** Solves a*x + b*y = p subject to x^2 + y^2 = min. */
  template<class T>
  static void solveMinNorm(T a, T b, T p, T* x, T* y);

  /** Solves the linear system of equations:
  \f[ a_{00}x_0      + a_{01}x_1      + \ldots + a_{0 (N-1)}    x_{N-1} &= b_0      \f] 
  \f[ a_{10}x_0      + a_{11}x_1      + \ldots + a_{1 (N-1)}    x_{N-1} &= b_1      \f] 
  \f[ \vdots                                                                        \f] 
  \f[ a_{(N-1) 0}x_0 + a_{(N-1) 1}x_1 + \ldots + a_{(N-1)(N-1)} x_{N-1} &= b_{N-1}  \f] 
  which is respresented by the matrix equation:
  \f[ \mathbf{A x} = \mathbf{b} \f]  
  that is, it computes the solution vector x which satisfies the system of equations (if any).
  Note that the indexing scheme for the matrix and array entries above conforms with the 
  C-array-indexing (starting at index 0) rather than the math-textbook-indexing (starting at
  index 1). The matrix A is represented as a pointer to pointer to double - the first 
  dereferencing is therefore a pointer to double and must point to the beginning of a row. The 
  vectors x and b are simple pointers pointing to an array of doubles which must be of length N 
  as well. When accessing A as two-dimensional array, this means that the the first index in A 
  must indicate the row, the second indicates the column. The vector b represents the right hand 
  side of the equation and x will contain the solution vetor on return. The function uses 
  Gaussian elimination with partial pivoting and subsequent backsubstitution. The boolean 
  return value informs whether a solution could be computed (the algorithm will fail when the 
  matrix is singular) - when it returns false, it means that either there is no solution at all 
  or that there is not a unique solution.  */
  template<class T>
  //static bool rsSolveLinearSystem(const T **A, T *x, const T *b, int N); // compiler error
  static bool rsSolveLinearSystem(T** A, T* x, const T* b, int N);
  // deprecate! 

  //static bool rsSolveLinearSystem(T **A, T *x, T *b, int N);
    // maybe have possibly different types for the matrix elements and vector elements - some 
    // equations (for curves in Salomon's Computer Graphics... for example) are formulated in terms 
    // of matrices-of-numbers and vectors-of-points

  /** Solves the linear system just as solveLinearSystem() does - but doesn't allocate temporary 
  heap memory and instead destroys the coefficient matrix A and the target vector b during the 
  process because the computation is done in place. In fact, the function solveLinearSystem() just 
  makes temporary copies of the matrix A and target vector b and then calls this function with these 
  copies to do the actual computation. If you don't need the matrix or vector anymore after solving 
  the system, you can use this function directly to get rid of the copying overhead. The algorithm 
  is Gaussian elimination with partial pivoting (...i think -> verify this). */
  template<class T>
  static bool rsSolveLinearSystemInPlace(T **A, T *x, T *b, int N);
  // deprecate!..but keep around for the comments

  // todo: add functions to solve NxM systems with N != M (find minimum-norm solution for 
  // underdetermined systems and least-squares approximation for overdetermined systems...
  // maybe we should have two different functions solveUnderDeterminedSystem, 
  // solveOverDeterminedSystem...or something. the overdetermined case can then be used
  // inside polynomial curve-fitting and multiple linear regression routines 
  // (i think polynomial fits use a Vandermonde matrix (created from data vectors) and then works
  // the same a multiple linear regression -> look it up...)
  // (maybe) make a dispatcher function  rsSolveLinearSystem(T **A, T *x, T *b, int N, int M);
  // that dispatches the 3 cases to the 3 functions

  /** Inverts the matrix A via Gauss-Jordan elimination with partial pivoting. */
  template<class T>
  static bool rsInvertMatrix(T **A, int N);
  // deprecate!






  /** Solves the tridiagonal system of equations defined by a NxN matrix having the 3 nonzero 
  diagonals 'lowerDiagonal', 'mainDiagonal' and 'upperDiagonal' where the 'mainDiagonal' array 
  should have N elements and 'lowerDiagonal' and 'upperDiagonal' represent the the diagonals below 
  and above the main diagonal respectively and should have N-1 elements. The 'rightHandSide' 
  argument represents the right hand side of the equation (vector with N elements) and 'solution' 
  ís where the solution vector will be stored (N elements). */
  template<class T>
  static bool rsSolveTridiagonalSystem(T *lowerDiagonal, T *mainDiagonal, T *upperDiagonal, 
    T *rightHandSide, T *solution, int N);

  /** Solves a pentadiagonal linear system of equations with given diagonals and right-hand side 
  using a simple algorithm without pivot-search. lowerDiag1 is the one directly below the main 
  diagonal, lowerDiag2 the one below lowerDiag1 - and similarly for upperDiag1/upperDiag2. In the 
  process of the computations, the right hand side vector is destroyed. the same is true for 
  mainDiag and the two inner sub/superdiagonals lowerDiag1, upperDiag1. Note also that you can't 
  use the same array for lowerDiag1 and upperDiag1, even if your matrix is symmetric.
  ..What about lowerDiag2/upperDiag2? are these preserved and may these point to the same vector? 
  It's probably safest to assume that everything may get messed up and all arrays should be 
  distinct. */
  template<class T>
  static bool rsSolvePentaDiagonalSystem(T* lowerDiag2, T* lowerDiag1, T* mainDiag, T* upperDiag1, 
    T* upperDiag2, T* righHandSide, T *solution, int N);


  /** Given NxN matrices A and B whose columns are assumed to constitue two bases for R^N and the 
  coordinates va of a vector v with respect to basis A, this function returns the coordinates of 
  the same vector v with respect to basis B and return them in vb. The return value informs, if the
  solution of the linear system succeeded (it may not, if B is not a basis - but what if A is not a 
  basis? will it work anyway? quite possible, i think - because it just restricts the vector v to 
  lie in some subspace of R^N which is no problem). */
  template<class T>
  static bool rsChangeOfBasisColumnWise(T **A, T **B, T *va, T *vb, int N);

  /** Similar to rsChangeOfBasisColumnWise but here, the basis vectors are given by the rows of the 
  matrices A and B. */
  template<class T>
  static bool rsChangeOfBasisRowWise(T **A, T **B, T *va, T *vb, int N);

  /** Given NxN matrices A and B whose columns are assumed to constitue two bases for R^N, this 
  function computes the change-of-basis matrix C which converts coordinates with respect to basis A 
  into coordinates with respect to basis B such that for any vector v with coordinates va in A and 
  coordinates vb in B, we have vb = C * va. */
  template<class T>
  static bool rsChangeOfBasisMatrixColumnWise(T **A, T **B, T **C, int N);

  /** Similar to rsChangeOfBasisMatrixColumnWise but here, the basis vectors are given by the rows 
  of the matrices A and B. */
  template<class T>
  static bool rsChangeOfBasisMatrixRowWise(T **A, T **B, T **C, int N);



};

#endif
