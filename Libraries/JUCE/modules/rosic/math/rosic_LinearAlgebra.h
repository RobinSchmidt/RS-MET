#ifndef rosic_LinearAlgebra_h
#define rosic_LinearAlgebra_h

// rosic-indcludes:
#include "rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /*

  This is a collection of functions for linear algrebraic manipulation of matrices and vectors 
  represented in raw C-style pointer representation.

  \todo: (maybe) templatize the functions so as to work on matrices of types other than double

  */

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
  bool solveLinearSystem(double** A, double* x, double* b, int N);

  /** Solves the linear system just as solveLinearSystem() but destroys the coefficient matrix A 
  and the target vector b during the process because the computation is done in placce. In 
  fact, the function solveLinearSystem() just makes temporary copies of the matrix A and target 
  vector b and the calls this function with these copies to do the actual computation. If you 
  don't need the matrix or vector anymore after solving the system, you can use this function 
  directly to get rid of the copying overhead. */
  bool solveLinearSystemInPlace(double** A, double* x, double* b, int N);

  /** Inverts the matrix A via Gauss-Jordan elimination with partial pivoting. */
  bool invertMatrix(double** A, int N);

  /** Solves the tridiagonal system of equations defined by a NxN matrix having the 3 nonzero 
  diagonals 'lowerDiagonal', 'mainDiagonal' and 'upperDiagonal' where the 'mainDiagonal' array 
  should have N elements and 'lowerDiagonal' and 'upperDiagonal' represent the the diagonals below 
  and above the main diagonal respectively and should have N-1 elements. The 'rightHandSide' 
  argument represents the right hand side of the equation (vector with N elements) and 'solution' 
  ís where the solution vector will be stored (N elements). */
  bool solveTriDiagonalSystem(double *lowerDiagonal, double *mainDiagonal, double *upperDiagonal, 
    double *rightHandSide, double *solution, int N);

}  // end namespace rosic

#endif // rosic_LinearAlgebra_h
