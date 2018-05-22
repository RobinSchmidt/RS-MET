#ifndef RAPT_LINEARALGEBRA_H_INCLUDED
#define RAPT_LINEARALGEBRA_H_INCLUDED

class rsLinearAlgebra
{

public:

  /** Solves the 2x2 system A*x = y where A is a 2x2 matrix and x, y are 2-dimensional vectors. The
  result is returned in x, the other parameters will not be modified. */
  template<class T>
  static void rsSolveLinearSystem2x2(const T A[2][2], T x[2], const T y[2]);

  /** Solves the 3x3 system A*x = y where A is a 3x3 matrix and x, y are 3-dimensional vectors. The
  result is returned in x, the other parameters will not be modified. */
  template<class T>
  static void rsSolveLinearSystem3x3(const T A[3][3], T x[3], const T y[3]);
    // todo: check, if this is the right form of const-ness

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
  //bool rsSolveLinearSystem(const T **A, T *x, const T *b, int N);
  static bool rsSolveLinearSystem(T **A, T *x, T *b, int N);

  /** Solves the linear system just as solveLinearSystem() does - but doesn't allocate temporary 
  heap memory and instead destroys the coefficient matrix A and the target vector b during the 
  process because the computation is done in place. In fact, the function solveLinearSystem() just 
  makes temporary copies of the matrix A and target vector b and the calls this function with these 
  copies to do the actual computation. If you don't need the matrix or vector anymore after solving 
  the system, you can use this function directly to get rid of the copying overhead. The algorithm 
  is Gaussian elimination with partial pivoting (...i think -> verify this). */
  template<class T>
  static bool rsSolveLinearSystemInPlace(T **A, T *x, T *b, int N);

  // todo: add functions to solve NxM systems with N != M (find minimum-norm solution for 
  // underdetermined systems and least-squares approximation for overdetermined systems...
  // maybe we should have two different functions solveUnderDeterminedSystem, 
  // solveOverDeterminedSystem...or something. the overdetermined case can then be used
  // inside polynomial curve-fitting and multiple linear regression routines 
  // (i think polynomial fits use a Vandermonde matrix (created from data vectors) and then works
  // the same a multiple linear regression -> look it up...)

  /** Inverts the matrix A via Gauss-Jordan elimination with partial pivoting. */
  template<class T>
  static bool rsInvertMatrix(T **A, int N);

  /** Solves the tridiagonal system of equations defined by a NxN matrix having the 3 nonzero 
  diagonals 'lowerDiagonal', 'mainDiagonal' and 'upperDiagonal' where the 'mainDiagonal' array 
  should have N elements and 'lowerDiagonal' and 'upperDiagonal' represent the the diagonals below 
  and above the main diagonal respectively and should have N-1 elements. The 'rightHandSide' 
  argument represents the right hand side of the equation (vector with N elements) and 'solution' 
  ís where the solution vector will be stored (N elements). */
  template<class T>
  static bool rsSolveTridiagonalSystem(T *lowerDiagonal, T *mainDiagonal, T *upperDiagonal, 
    T *rightHandSide, T *solution, int N);

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
