#ifndef RAPT_MATRIXTOOLS_H_INCLUDED
#define RAPT_MATRIXTOOLS_H_INCLUDED

/** A collection of functions for basic handling 2-dimensional arrays (i.e. matrices). Fancy linear 
algebra stuff is not included here - this is part of the math module. */

// note that this code is mostly only here for legacy reasons - the now preferred approach to 
// handle matrices ins embodied in clas rsMatrix (using flat arrays instead of an array-of-arrays)
// but some older code still does it the old way, so we still need to keep it around
// todo: move to legacy folder
// todo: rename to rsMatrixTools ..and get rid of the rs prefixes in the function names



class rsMatrixTools
{

public:

  template<class T>
  static void rsAllocateMatrix(T**& A, int N, int M);

  template<class T>
  static void rsDeAllocateMatrix(T**& A, int N, int M); 
  // rename to freeMatrix - ..and why is the M needed?

  template<class T>
  static void rsInitMatrix(T** A, int N, int M, T value = T(0));

  template<class T>
  static void rsCopyMatrix(T** source, T **destination, int N, int M);

  template<class T>
  static bool rsAreMatricesApproximatelyEqual(T **A, T **B, int N, int M, T tol);

  /* Pre-multiplies the M-dimensional vector x by NxM matrix A and stores the resulting 
  N-dimensional vector in y, such that y = A * x. */
  template<class T>
  static void rsMatrixVectorMultiply(T **A, T *x, T *y, int N, int M);

  /** Pre-multiplies the N-dimensional vector x by the transposed NxM matrix A (which is MxN) and 
  stores the resulting M-dimensional vector in y, such that y = A^T * x. */
  template<class T>
  static void rsTransposedMatrixVectorMultiply(T **A, T *x, T *y, int N, int M);

  /** Multiplies NxM matrix A by MxP matrix B and stores the result in NxP matrix C = A * B. */
  template<class T>
  static void rsMatrixMultiply(T **A, T **B, T **C, int N, int M, int P);

  /** Multiplies the transposed NxM matrix A (A^T is MxN) by NxP matrix B and stores the result in 
  MxP matrix C = A^T * B. */
  template<class T>
  static void rsMatrixMultiplyFirstTransposed(T **A, T **B, T **C, int N, int M, int P);

  /** Multiplies the NxM matrix A by the transposed PxM matrix B (B^T is MxP) and stores the result 
  in NxP matrix C = A * B^T. */
  template<class T>
  static void rsMatrixMultiplySecondTransposed(T **A, T **B, T **C, int N, int M, int P);

  /** Multiplies the transposed NxM matrix A (A^T is MxN) by the transposed PxN matrix B (B^T is 
  NxP) and stores the result in MxP matrix C = A^T * B^T. */
  template<class T>
  static void rsMatrixMultiplyBothTransposed(T **A, T **B, T **C, int N, int M, int P);

  /** Given an NxM matrix A and an MxM (square-)matrix B, this function computes the product and 
  stores it in A, such that A will be replaced by A * B. */
  template<class T>
  static void rsMatrixInPlaceMultiply(T **A, T **B, int N, int M);
    // rename to rsMatrixInPlacePostMultiply - a similar function can be written when the first
    // factor is a square matrix - then the result can replace the 2nd matrix and it could be 
    // called in-place-pre-multiply
    // allocates heap-memory - maybe make a version that takes a workspace parameter
};



#endif