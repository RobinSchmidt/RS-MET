#ifndef RS_MATRIX_H
#define RS_MATRIX_H

namespace RSLib
{

  /**

  This is a class for representing matrices and doing mathematical operations with them. To make
  the class fast, no consistency checking is done in the mathematical operators. You must ensure
  yourself that the input arguments are compatible - for example don't try to add two matrices of
  different dimensionality or multiply an n times m matrix with a q times p matrix when m != q.

  \todo: optimize: Avoid excessive deep-copying of the 2D arrays by defining an rsMatrixData class
  and keeping a reference (pointer) to an object of that class in the actual rsMatrix class. 
  Asssignment operator and copy constructor only copy the reference pointer and the rsMatrixData 
  class keeps a reference counter. All non-const member functions will create a deep copy of the 
  actual data, iff they actually change some value and refCount > 1 (we probably should have a
  createDeepDataCopy member function for that purpose). The access operator should be made const
  and/or return a const reference - for actually setting values, we should use a separate 
  set-function.


  The more sophisticated matrix computations (like singular value decomposition and
  eigendecomposition) are powered by the Template Numerical Toolkit (TNT) ....not yet implemented


  References:
  -(1): Template Numerical Toolkit (TNT), http://math.nist.gov/tnt/

  */

  template<class T>
  class rsMatrix
  {

  public:

    /** \name Public Data Members \todo: make them protected */

    /** Number of rows. */
    int numRows;

    /** Number of colums. */
    int numColumns;

      // \todo maybe use rsUint32 instead of int

    /** The actual values of the matrix. The first index refers to the row, the second to the
    column. */
    T **m;

    /** The matrix as flat array - this is only used to facilitate contiguous memory allocation. */
    T *mFlat;


    /** \name Construction/Destruction */

    /** Default constructor - constructs a matrix with zero rows, zero colums and NULL pointers. */
    rsMatrix();
      // not tested

    /** Constructor. You must pass the number of rows and colums here. */
    rsMatrix(int numRows, int numColumns, bool initElementsWithZeros = false);
      // not tested

    /** Constructor. You must pass the number of rows and colums here and the values to intialize
    the matix. */
    rsMatrix(int numRows, int numColumns, T **values);
      // not tested

    /** Copy constructor. */
    rsMatrix(const rsMatrix& other);
      // not tested

    /** Returns a square matrix with diagonal elements given by 'scalarValue' and zero off-diagonal
    elements. */
    static rsMatrix scalarMatrix(T scalarValue, int dimension)
    {
      rsMatrix result(dimension, dimension);
      for(int r=0; r<dimension; r++)  // use rsInitMatrix instead
      {
        for(int c=0; c<dimension; c++)
          result.m[r][c] = 0.0;
      }
      for(int r=0; r<dimension; r++)
        result.m[r][r] = scalarValue;
      return result;
    }
      // not tested

    /** Returns a square matrix with diagonal elements given by the array 'diagonalValues' and zero
    off-diagonal elements. The array is assumed to contain a number of elements equal to
    'dimension'. */
    static rsMatrix diagonalMatrix(T *diagonalValues, int dimension)
    {
      rsMatrix result(dimension, dimension);
      for(int r=0; r<dimension; r++)
      {
        for(int c=0; c<dimension; c++)
          result.m[r][c] = 0.0;
      }
      for(int r=0; r<dimension; r++)
        result.m[r][r] = diagonalValues[r];
      return result;
    }
      // not tested

    /** Destructor. */
    ~rsMatrix();


    /** \name Manipulations */

    /** Sets up the size of the matrix (number of rows and columns). */
    void setSize(int newNumRows, int newNumColumns)
    {
      // re-allocate memory, if necessary:
      if( numRows != newNumRows || numColumns != newNumColumns )
      {
        delete[] m;
        delete[] mFlat;
        numRows    = newNumRows;
        numColumns = newNumColumns;
        mFlat      = new T[numRows*numColumns];
        m          = new T*[numRows];
        for(int r=0; r<numRows; r++)
          m[r] = &mFlat[r*numColumns];
      }
    }
      // not tested

    /** Sets one of the matrix elements. */
    void set(int row, int column, T newValue)
    {
      rsAssert( row < numRows && column < numColumns ); // invalid row or column index
      m[row][column] = newValue;
    }

    /** Transposes this matrix (exchanges the roles of rows and columns). */
    void transpose()
    {
      if( numRows == numColumns )
        rsTransposeSquareArray(m, numRows);
      else
      {
        // create a new 2D-array:
        int numRowsNew    = numColumns;
        int numColumnsNew = numRows;
        T *mFlatNew  = new T[numRowsNew*numColumnsNew];
        T **mNew     = new T*[numRowsNew];
        for(int r=0; r<numRowsNew; r++)
          mNew[r] = &mFlatNew[r*numColumnsNew];

        // copy the values:
        for(int r=0; r<numRowsNew; r++)
        {
          for(int c=0; c<numColumnsNew; c++)
            mNew[r][c] = m[c][r];
        }

        // free old memory and re-assign pointers
        delete[] m;
        delete[] mFlat;
        m          = mNew;
        mFlat      = mFlatNew;
        numRows    = numRowsNew;
        numColumns = numColumnsNew;
      }
    }

    /** Assigns random values between 'min' and 'max' to each element. */
    void randomizeElements(T min, T max)
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] = rsRandomUniform(min, max);
      }
    }
      // not tested

    /** Initializes all elements with zero values. */
    // \todo use rsInitMatrix
    void initWithZeros()
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] = 0.0;
      }
    }
      // not tested

    /** Initializes all elements with the sum of their row- and column indices. */
    // why is this useful?
    void initWithIndexSum()
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] = r+c;
      }
    }
      // not tested

    /** Applies the passed function to each element of the matrix. */
    void applyFunction( T (*f) (T) )
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] = f(m[r][c]);
      }
    }


    /** \name Inquiry */

    int getNumRows() const 
    {
      return numRows;
    }

    int getNumColumns() const 
    {
      return numColumns;
    }

    bool isSquare() const 
    { 
      return numRows == numColumns; 
    }

    // getRank, getConditionNumber, getTrace, getDeterminant, getFrobeniusNorm, getEuclideanNorm

    //rsMatrix getSubMatrix(int fromRow, int toRow, int fromColumn, int toColumn);
    //rsMatrix withRemovedRow(int rowToRemove);
    //rsMatrix withRemovedColumn(int columnToRemove);


    /** \name Operators */

    /** Accesses the element at given index-pair for reading. */
    const T& operator()(const int i, const int j) const
    {
      return m[i][j];
    }

    /** Assigns one matrix with another one. */
    rsMatrix& operator=(const rsMatrix& m2)
    {
      // re-allocate memory, if necessary:
      if( this->numRows != m2.numRows || this->numColumns != m2.numColumns )
      {
        delete[] m;
        delete[] mFlat;
        numRows    = m2.numRows;
        numColumns = m2.numColumns;
        mFlat      = new T[numRows*numColumns];
        m          = new T*[numRows];
        for(int r=0; r<numRows; r++)
          m[r] = &mFlat[r*numColumns];
      }

      // copy the values:
      memcpy(mFlat, m2.mFlat, numRows*numColumns*sizeof(T));

      return *this;
    }

    /** Compares two matrices for equality. */
    bool operator==(const rsMatrix& m2) const
    {
      if( numRows != m2.numRows || numColumns != m2.numColumns )
        return false;
      else
      {
        for(int r=0; r<numRows; r++)
        {
          for(int c=0; c<numColumns; c++)
          {
            if( m[r][c] != m2.m[r][c] )
              return false;
          }
        }
      }
      return true;
    }

    /** Compares two matrices for inequality. */
    bool operator!=(const rsMatrix& m2) const
    {
      return !(*this == m2);
    }

    /** Defines the negative of a matrix. */
    rsMatrix operator-()
    {
      rsMatrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = -m[r][c];
      }
      return result;
    }

    /** Adds another matrix to this matrix and returns the result. */
    rsMatrix& operator+=(const rsMatrix &m2)
    {
      rsAssert( numRows == m2.numRows && numColumns == m2.numColumns ); // matrices incompatible
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] += m2.m[r][c];
      }
      return *this;
    }

    /** Adds a scalar to this matrix and returns the result. */
    rsMatrix& operator+=(const T &x)
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] += x;
      }
      return *this;
    }

    /** Subtracts another matrix from this matrix and returns the result. */
    rsMatrix& operator-=(const rsMatrix &m2)
    {
      rsAssert( numRows == m2.numRows && numColumns == m2.numColumns ); // matrices incompatible
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] -= m2.m[r][c];
      }
      return *this;
    }

    /** Subtracts a scalar from this matrix and returns the result. */
    rsMatrix& operator-=(const T &x)
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] -= x;
      }
      return *this;
    }

    /** Right-multiplies this matrix with another matrix and returns the result. */
    rsMatrix& operator*=(const rsMatrix &m2)
    {
      rsAssert( this->numColumns == m2.numRows );  // matrices incompatible

      if( m2.numColumns == numColumns )
        rsMatrixInPlaceMultiply(m, m2.m, numRows, numColumns);
      else
      {
        // allocate memory for result:
        T *mFlatNew = new T[numRows*m2.numColumns];
        T **mNew    = new T*[numRows];
        for(int r = 0; r < numRows; r++)
          mNew[r] = &mFlatNew[r*m2.numColumns];

        // compute result:
        rsMatrixMultiply(m, m2.m, mNew, numRows, numColumns, m2.numColumns);

        // delete old data:
        delete[] m;
        delete[] mFlat;

        // update data/pointers:
        m          = mNew;
        mFlat      = mFlatNew;
        numColumns = m2.numColumns;
      }

      return *this;
    }

    /** Multiplies this matrix by a scalar and returns the result. */
    rsMatrix& operator*=(const T &x)
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] *= x;
      }
      return *this;
    }

    /** Divides this matrix by a scalar and returns the result. */
    rsMatrix& operator/=(const T &x)
    {
      T scale = 1.0 / x;
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] *= scale;
      }
      return *this;
    }

    /** Adds two matrices. */
    rsMatrix operator+(const rsMatrix &m2)
    {
      rsAssert( numRows == m2.numRows && numColumns == m2.numColumns ); // matrices incompatible
      rsMatrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] + m2.m[r][c];
      }
      return result;
    }

    /** Adds a matrix and a scalar. */
    rsMatrix operator+(const T &x)
    {
      rsMatrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] + x;
      }
      return result;
    }

    /** Subtracts two matrices. */
    rsMatrix operator-(const rsMatrix &m2)
    {
      rsAssert( numRows == m2.numRows && numColumns == m2.numColumns ); // matrices incompatible
      rsMatrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] - m2.m[r][c];
      }
      return result;
    }

    /** Subtracts a scalar from a matrix. */
    rsMatrix operator-(const T &x)
    {
      rsMatrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] - x;
      }
      return result;
    }

    /** Multiplies a matrix and a scalar. */
    rsMatrix operator*(const T &x)
    {
      rsMatrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] * x;
      }
      return result;
    }

    /** Multiplies two matrices. */
    rsMatrix operator*(const rsMatrix &m2)
    {
      rsAssert( this->numColumns == m2.numRows );  // matrices incompatible
      rsMatrix result(numRows, m2.numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<m2.numColumns; c++)
        {
          result.m[r][c] = 0.0;
          for(int k=0; k<numColumns; k++)
            result.m[r][c] += m[r][k] * m2.m[k][c];
        }
      }
      return result;
    }

    /** Divides a matrix by a scalar. */
    rsMatrix operator/(const T &x)
    {
      T scale = 1.0 / x;
      rsMatrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] * scale;
      }
      return result;
    }


    /** \name Decompositions */

    /** Returns the the singular value decomposition of this matrix (lets denote it by A) such that
    A = U * S * V' where U and V are orthogonal matrices (U' * U = E, V' * V = E) and S is a
    diagonal matrix with positive or zero elements (the singular values). The singular values are
    ordered such that S[0][0] >= S[1][1] >= ...  */
    void getSingularValueDecomposition(rsMatrix *U, rsMatrix *S, rsMatrix *V);

    // getEigenDecomposition, getLowerUpperDecomposition, ...


    //---------------------------------------------------------------------------------------------
    // others:

    /** Prints the values to the standard output - mainly for debugging purposes. */
    void print();

  }; // end of class rsMatrix

  // some binary operators are defined outside the class such that the left hand operand does
  // not necesarrily need to be of class rsMatrix
  // \todo: maybe get rid of the inlining

  /** Adds a scalar and a matrix. */
  template<class T>
  RS_INLINE rsMatrix<T> operator+(const T &x, const rsMatrix<T> &m)
  {
    rsMatrix<T> result(m.numRows, m.numColumns);
    for(int r=0; r<m.numRows; r++)
    {
      for(int c=0; c<m.numColumns; c++)
        result.m[r][c] = m.m[r][c] + x;
    }
    return result;
  }

  /** Subtracts a matrix from a scalar. */
  template<class T>
  RS_INLINE rsMatrix<T> operator-(const T &x, const rsMatrix<T> &m)
  {
    rsMatrix<T> result(m.numRows, m.numColumns);
    for(int r=0; r<m.numRows; r++)
    {
      for(int c=0; c<m.numColumns; c++)
        result.m[r][c] = x - m.m[r][c];
    }
    return result;
  }

  /** Multiplies a scalar and a matrix. */
  template<class T>
  RS_INLINE rsMatrix<T> operator*(const T &x, const rsMatrix<T> &m)
  {
    rsMatrix<T> result(m.numRows, m.numColumns);
    for(int r=0; r<m.numRows; r++)
    {
      for(int c=0; c<m.numColumns; c++)
        result.m[r][c] = x * m.m[r][c];
    }
    return result;
  }

  /** Returns matrix A transposed. */
  template<class T>
  RS_INLINE rsMatrix<T> trans(const rsMatrix<T>& A)
  {
    rsMatrix<T> B(A.numColumns, A.numRows);
    for(int r=0; r<B.numRows; r++)
    {
      for(int c=0; c<B.numColumns; c++)
        B.m[r][c] = A.m[c][r];
    }
    return B;
  }



  // matrix/vector functions:


  /** Left-multiplies a vector by a matrix such that c = A*b with A being the input matrix, b being
  the input vector and c being the output vector. Thus, b is interpreted as a column-vector and the
  number of columns in the matrix A must equal the number of elements in the vector b. The result
  will be a vector with the number rows of A as its dimensionality. */
  template<class T>
  RS_INLINE rsVector<T> operator*(const rsMatrix<T> &A, const rsVector<T> &b)
  {
    rsAssert( A.numColumns == b.dim );  // matrix and vector incompatible
    rsVector<T> c(A.numRows);
    for(int i=0; i<A.numRows; i++)
    {
      c.v[i] = 0.0;
      for(int j=0; j<A.numColumns; j++)
        c.v[i] += A.m[i][j] * b.v[j];
    }
    return c;
  }
      // not tested

  /** Right-multiplies a transposed vector by a matrix such that c = b^T * A with A being the input
  matrix, b being the input vector and c being the output vector. Thus, b is interpreted as a
  row-vector and the number of rows in the matrix A must equal the number of elements in the
  vector b. The result will be a matrix consisting of one row and a number of columns equal to the
  dimenstionality of the vector b. */
  template<class T>
  RS_INLINE rsMatrix<T> operator*(const rsVector<T> &b, const rsMatrix<T> &A)
  {
    rsAssert( A.numRows == b.dim );  // matrix and vector incompatible
    rsMatrix<T> C(1, b.dim);
    for(int i=0; i<C.numColumns; i++)
    {
      C.m[0][i] = 0.0;
      for(int j=0; j<b.dim; j++)
        C.m[0][i] += b.v[j] * A.m[j][i];
    }
    return C;
  }
      // not tested

  /** Computes the outer product a * b^T between two vectors (which should have the same
  dimensionality) - the result is a square matrix. Note that the outer product is not
  commutative, specifically a * b^T == (b * a^T)^T. */
  template<class T>
  RS_INLINE rsMatrix<T> outerProduct(const rsVector<T> &a, const rsVector<T> &b)
  {
    rsAssert( a.dim == b.dim );   // vectors incopatible
    rsMatrix<T> result(a.dim, a.dim);
    for(int i=0; i<result.numRows; i++)
    {
      for(int j=0; j<result.numColumns; j++)
        result.m[i][j] = a.v[i] * b.v[j];
    }
    return result;
  }
      // not tested

  // typedefs for explicit instantiations:
  typedef rsMatrix<double> RSLib_API rsMatrixDbl;

}

#endif
