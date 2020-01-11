#ifndef RAPT_LINEARALGEBRA_NEW_H_INCLUDED
#define RAPT_LINEARALGEBRA_NEW_H_INCLUDED


/** Collection of functions for linear algebra such as solving systems of linear equations, matrix 
inversion, etc. 

Many algorithms destroy the original input matrices and/or right-hand-sides in the course of their 
operation. If a parameter like a pointer to an array or a reference to an rsMatrixView is passed as
non-const, client code should always assume that the argument will contain either a desired result 
or garbage when the function returns. If you need to keep the inputs, pass a copy. */

// these implementations are based on the new class rsMatrixView which works with matrices in flat
// (row major) storage format - eventually, they should replace the old implementations (which 
// represent matrices as array-of-array/pointer-to-pointer) and the old ones shall be deprecated

class rsLinearAlgebraNew
{

public:


  // newer versions using rsMatrixView - document them (use documentation from old implementations, 
  // keep the old functions as deprecated legacy functions around as long as they are still needed)

  // maybe provide an API that takes the flat arrays of the matrices directly as inputs

  /** \name Convenience functions 
  These functions are convenient to use but probably not recommendable for realtime code (they
  allocate memory) */

  template<class T>
  static std::vector<T> solveLinearSystem(rsMatrixView<T>& A, std::vector<T>& b);
  // allocates
  // make inputs const
  // rename to solve

  template<class T>
  static rsMatrix<T> inverse(const RAPT::rsMatrixView<T>& A);
  // allocates


  /** \name Solvers */

  template<class T>
  static bool solveLinearSystem(rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B);
  // doesn't allocate, destroys A,B
  // rename to solve






  /** \name Subroutines */

  /** Transforms the augmented coefficient matrix A|B for a set of linear systems of equations of 
  the form A * X = B into an upper triangular form, a.k.a. "row echelon form". In this form, the 
  corresponding linear system is easy to solve. It is a *set* of linear system because B may 
  contain several right hand side vectors and the matrix X should contain just as many solutions. 
  The process used is Gaussian elimination with partial pivoting. */
  template<class T>
  static bool makeSystemUpperTriangular(rsMatrixView<T>& A, rsMatrixView<T>& B);
  // doesn't allocate
  // rename to makeUpperTriangular or makeTriangular or rowEchelon

  /** Simplified version that doesn't use pivoting - this may fail even for non-singular 
  matrices, so it's not recommended for general use. It has less requirements on the datatype T 
  (doesn't require comparison operators <,> or an rsAbs function defined for it - it does, 
  however, require an == operator). */
  template<class T>
  static bool makeSystemUpperTriangularNoPivot(rsMatrixView<T>& A, rsMatrixView<T>& B);
  // doesn't allocate

  /** This produces the so called *reduced* row echelon form of a linear system of equations. In 
  this form, the matrix is diagonal, so it's trivial to solve
  Note that this has *nothing* to do 
  with "diagonalization" (which refers to finding the eigenvalues and -vectors and expressing that 
  matrix as product: A = inv(M) * D * M where M is a matrix whose columns are the eigenvectors and 
  D is a diagonal matrix with the eigenvalues on the main diagonal - this is NOT what this function 
  does). */
  template<class T>
  static bool makeSystemDiagonal(rsMatrixView<T>& A, rsMatrixView<T>& B);
  // needs test, doesn't allocate, maybe remove from library - doesn't seem to be useful
  // rename to makeDiagonal or reducedRowEchelon


  template<class T>
  static void solveUpperTriangularSystem(
    rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B);
  // doesn't allocate, rename to solveTriangular - document, if it can be used in place, i.e. X==B

};



#endif