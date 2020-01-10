#ifndef RAPT_LINEARALGEBRA_NEW_H_INCLUDED
#define RAPT_LINEARALGEBRA_NEW_H_INCLUDED


/** Collection of functions for linear algebra such as solving systems of linear equations, matrix 
inversion, etc. 

Many algorithms destroy the original input matrices and/or right-hand-sides in the course of their 
operation. If a parameter like a pointer to an array or a reference to an rsMatrixView is passed as
non-const, always assume that the input will contain either a desired result or garbage when the 
function returns. If you need to keep the inputs, pass a copy. */

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

  template<class T>
  static rsMatrix<T> inverse(const RAPT::rsMatrixView<T>& A);
  // allocates


  /** \name Solvers */

  template<class T>
  static bool solveLinearSystem(rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B);
  // doesn't allocate, destroys A,B


  template<class T>
  static bool makeSystemDiagonal(rsMatrixView<T>& A, rsMatrixView<T>& B);
  // not yet finished, doesn't allocate



  /** \name Subroutines */

  template<class T>
  static bool makeSystemUpperTriangular(rsMatrixView<T>& A, rsMatrixView<T>& B);
  // doesn't allocate

  template<class T>
  static void solveUpperTriangularSystem(
    rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B);
  // doesn't allocate









};



#endif