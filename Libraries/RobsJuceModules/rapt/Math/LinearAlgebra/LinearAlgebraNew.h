#ifndef RAPT_LINEARALGEBRA_NEW_H_INCLUDED
#define RAPT_LINEARALGEBRA_NEW_H_INCLUDED


/** Collection of functions for linear algebra such as solving systems of linear equations, matrix 
inversion, etc. */

// these implementations are based on the new class rsMatrixView which works with matrices in flat
// (row major) storage format - eventually, they should replace the old implementations (which 
// represent matrices as array-of-array/pointer-to-pointer) and the old ones shall be deprecated

class rsLinearAlgebraNew
{

public:


  // newer versions using rsMatrixView - document them (use documentation from old implementations, 
  // keep the old functions as deprecated legacy functions around as long as they are still needed)

  // maybe provide an API that takes the flat arrays of the matrices directly as inputs

  template<class T>
  static bool makeSystemUpperTriangular(rsMatrixView<T>& A, rsMatrixView<T>& B);
  // doesn't allocate

  template<class T>
  static void solveUpperTriangularSystem(
    rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B);
  // doesn't allocate

  template<class T>
  static bool solveLinearSystem(rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B);
  // doesn't alloctae


  // convenience function:
  template<class T>
  static std::vector<T> solveLinearSystem(rsMatrixView<T>& A, std::vector<T>& b);
  // allocates
  // make inputs const

  template<class T>
  static rsMatrix<T> inverse(const RAPT::rsMatrixView<T>& A);
  // allocates


};



#endif