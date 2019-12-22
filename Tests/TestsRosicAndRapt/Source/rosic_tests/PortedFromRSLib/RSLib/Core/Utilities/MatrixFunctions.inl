#ifndef RS_MATRIXFUNCTIONS_INL
#define RS_MATRIXFUNCTIONS_INL

namespace RSLib
{

  template<class T>
  void allocateMatrix(T**& A, int N, int M)
  {
    A = new T*[N];
    for(int i = 0; i < N; i++)
      A[i] = new T[M];

    // can we rewrite this function in a way that ensures a contiguous chunk of memory to be
    // allocated?
    // maybe like so:
    //A = new T*[N];
    //T *tmp = new T[N*M];
    //for(int i = 0; i < N; i++)
    //  A[i] = &tmp[i*M];
    // in this case, the deallocation has to be modified also - but how? the "tmp" variable is local
    // can we access it as &A[0] or something?
    // make performance tests with both versions
  }

  template<class T>
  void deallocateMatrix(T**& A, int N, int M)
  {
    for(int i = 0; i < N; i++)
      delete[] A[i];
    delete[] A;
  }

  template<class T>
  void initMatrix(T** A, int N, int M, T value)
  {
    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < M; j++)
        A[i][j] = value;
    }
  }

  template<class T>
  void copyMatrix(T** source, T **destination, int N, int M)
  {
    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < M; j++)
        destination[i][j] = source[i][j];
    }
  }

  template<class T>
  bool areMatricesApproximatelyEqual(T **A, T **B, int N, int M, T tol)
  {
    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < M; j++)
      {
        T d = rsAbs( A[i][j] - B[i][j] );
        if( d > tol )
          return false;
      }
    }
    return true;
  }

  template<class T>
  void rsMatrixVectorMultiply(T **A, T *x, T *y, int N, int M)
  {
    for(int i = 0; i < N; i++)
    {
      y[i] = T(0);
      for(int j = 0; j < M; j++)
        y[i] += A[i][j] * x[j];
    }
  }

  template<class T>
  void rsTransposedMatrixVectorMultiply(T **A, T *x, T *y, int N, int M)
  {
    for(int i = 0; i < M; i++)
    {
      y[i] = T(0);
      for(int j = 0; j < N; j++)
        y[i] += A[j][i] * x[j];
    }
  }

  template<class T>
  void rsMatrixMultiply(T **A, T **B, T **C, int N, int M, int P)
  {
    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < P; j++)
      {
        C[i][j] = T(0);
        for(int k = 0; k < M; k++)
          C[i][j] += A[i][k] * B[k][j];
      }
    }
  }

  template<class T>
  void rsMatrixMultiplyFirstTransposed(T **A, T **B, T **C, int N, int M, int P)
  {
    for(int i = 0; i < M; i++)
    {
      for(int j = 0; j < P; j++)
      {
        C[i][j] = T(0);
        for(int k = 0; k < N; k++)
          C[i][j] += A[k][i] * B[k][j];
      }
    }
  }

  template<class T>
  void rsMatrixMultiplySecondTransposed(T **A, T **B, T **C, int N, int M, int P)
  {
    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < P; j++)
      {
        C[i][j] = T(0);
        for(int k = 0; k < M; k++)
          C[i][j] += A[i][k] * B[j][k];
      }
    }
  }

  template<class T>
  void rsMatrixMultiplyBothTransposed(T **A, T **B, T **C, int N, int M, int P)
  {
    for(int i = 0; i < M; i++)
    {
      for(int j = 0; j < P; j++)
      {
        C[i][j] = T(0);
        for(int k = 0; k < N; k++)
          C[i][j] += A[k][i] * B[j][k];
      }
    }
  }

  template<class T>
  void rsMatrixInPlaceMultiply(T **A, T **B, int N, int M)
  {
    T *Ai = new T[M];  // i-th row of A (for temporary storage)
    for(int i = 0; i < N; i++)
    {
      rsCopyBuffer(A[i], Ai, M);
      for(int j = 0; j < M; j++)
      {
        A[i][j] = T(0);
        for(int k = 0; k < M; k++)
          A[i][j] += Ai[k] * B[k][j];
      }
    }
    delete[] Ai;
  }

}

#endif