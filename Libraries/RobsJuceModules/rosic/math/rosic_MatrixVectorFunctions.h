#ifndef rosic_MatrixVectorFunctions_h
#define rosic_MatrixVectorFunctions_h

//// rosic includes:
//#include "rosic_Vector.h" 
//#include "rosic_Matrix.h"

namespace rosic
{

/** A collection of functions that operate on matrices and/or vectors */

  /** Left-multiplies a vector by a matrix such that c = A*b with A being the input matrix, b being
  the input vector and c being the output vector. Thus, b is interpreted as a column-vector and the
  number of columns in the matrix A must equal the number of elements in the vector b. The result
  will be a vector with the number rows of A as its dimensionality. */
inline Vector operator*(const Matrix &A, const Vector &b)
{
  rassert(A.numColumns == b.dim);  // matrix and vector incompatible
  Vector c(A.numRows);
  for(int i=0; i<A.numRows; i++)
  {
    c.v[i] = 0.0;
    for(int j=0; j<A.numColumns; j++)
      c.v[i] += A.m[i][j] * b.v[j];
  }
  return c;
}
// check

/** Right-multiplies a transposed vector by a matrix such that c = b^T * A with A being the input
matrix, b being the input vector and c being the output vector. Thus, b is interpreted as a
row-vector and the number of rows in the matrix A must equal the number of elements in the
vector b. The result will be a matrix consisting of one row and a number of columns equal to the
dimenstionality of the vector b. */
inline Matrix operator*(const Vector &b, const Matrix &A)
{
  rassert(A.numRows == b.dim);  // matrix and vector incompatible
  Matrix C(1, b.dim);
  for(int i=0; i<C.numColumns; i++)
  {
    C.m[0][i] = 0.0;
    for(int j=0; j<b.dim; j++)
      C.m[0][i] += b.v[j] * A.m[j][i];
  }
  return C;
}
// check

/** Computes the outer product a * b^T between two vectors (which should have the same
dimensionality) - the result is a square matrix. Note that the outer product is not
commutative, specifically a * b^T == (b * a^T)^T. */
inline Matrix outerProduct(const Vector &a, const Vector &b)
{
  rassert(a.dim == b.dim);   // vectors incopatible
  Matrix result(a.dim, a.dim);
  for(int i=0; i<result.numRows; i++)
  {
    for(int j=0; j<result.numColumns; j++)
      result.m[i][j] = a.v[i] * b.v[j];
  }
  return result;
}
// check


// functions that operate on reference or pointer variables like 
// mul(const Matrix& A, const Vector& b, Vector& c)
// more efficient due to computation without copying in the assignment operators

} // end namespace rosic

#endif // #ifndef rosic_MatrixVectorFunctions_h