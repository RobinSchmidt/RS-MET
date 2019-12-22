#ifndef RS_MATRIXFUNCTIONS_H
#define RS_MATRIXFUNCTIONS_H

namespace RSLib
{
  // Here, we have some functions for basic handling of 2D arrays (i.e. matrices). But fancy linear 
  // algebra stuff is not included here - this is part of the math module.

  template<class T>
  void allocateMatrix(T**& A, int N, int M);

  template<class T>
  void deallocateMatrix(T**& A, int N, int M);

  template<class T>
  void initMatrix(T** A, int N, int M, T value = T(0));

  template<class T>
  void copyMatrix(T** source, T **destination, int N, int M);

  template<class T>
  bool areMatricesApproximatelyEqual(T **A, T **B, int N, int M, T tol);

  /* Pre-multiplies the M-dimensional vector x by NxM matrix A and stores the resulting 
  N-dimensional vector in y, such that y = A * x. */
  template<class T>
  void rsMatrixVectorMultiply(T **A, T *x, T *y, int N, int M);

  /** Pre-multiplies the N-dimensional vector x by the transposed NxM matrix A (which is MxN) and 
  stores the resulting M-dimensional vector in y, such that y = A^T * x. */
  template<class T>
  void rsTransposedMatrixVectorMultiply(T **A, T *x, T *y, int N, int M);

  /** Multiplies NxM matrix A by MxP matrix B and stores the result in NxP matrix C = A * B. */
  template<class T>
  void rsMatrixMultiply(T **A, T **B, T **C, int N, int M, int P);

  /** Multiplies the transposed NxM matrix A (A^T is MxN) by NxP matrix B and stores the result in 
  MxP matrix C = A^T * B. */
  template<class T>
  void rsMatrixMultiplyFirstTransposed(T **A, T **B, T **C, int N, int M, int P);

  /** Multiplies the NxM matrix A by the transposed PxM matrix B (B^T is MxP) and stores the result 
  in NxP matrix C = A * B^T. */
  template<class T>
  void rsMatrixMultiplySecondTransposed(T **A, T **B, T **C, int N, int M, int P);

  /** Multiplies the transposed NxM matrix A (A^T is MxN) by the transposed PxN matrix B (B^T is 
  NxP) and stores the result in MxP matrix C = A^T * B^T. */
  template<class T>
  void rsMatrixMultiplyBothTransposed(T **A, T **B, T **C, int N, int M, int P);

  /** Given an NxM matrix A and an MxM (square-)matrix B, this function computes the product and 
  stores it in A, such that A will be replaced by A * B. */
  template<class T>
  void rsMatrixInPlaceMultiply(T **A, T **B, int N, int M);
    // rename to rsMatrixInPlacePostMultiply - a similar function can be written when the first
    // factor is a square matrix - then the result can replace the 2nd matrix and it could be 
    // called in-place-pre-multiply

}

#endif