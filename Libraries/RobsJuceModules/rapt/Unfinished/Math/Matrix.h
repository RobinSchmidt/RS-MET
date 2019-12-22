#ifndef RAPT_MATRIX_H_INCLUDED
#define RAPT_MATRIX_H_INCLUDED

  /** This is a class for representing matrices and doing mathematical operations with them. To 
  make the class fast, no consistency checking is done in the mathematical operators. You must 
  ensure yourself that the input arguments are compatible - for example don't try to add two 
  matrices of different dimensionality or multiply an n times m matrix with a q times p matrix when
  m != q.

  \todo: optimize: Avoid excessive deep-copying of the 2D arrays by defining an rsMatrixData class
  and keeping a reference (pointer) to an object of that class in the actual rsMatrixOld class.
  Asssignment operator and copy constructor only copy the reference pointer and the rsMatrixData
  class keeps a reference counter. All non-const member functions will create a deep copy of the
  actual data, iff they actually change some value and refCount > 1 (we probably should have a
  createDeepDataCopy member function for that purpose). The access operator should be made const
  and/or return a const reference - for actually setting values, we should use a separate
  set-function.

  \todo: done - but: we could go back to the old implementation and just implement move constructor
  and move assignment operator, see here:
  http://blog.smartbear.com/c-plus-plus/c11-tutorial-introducing-the-move-constructor-and-the-move-assignment-operator/



  The more sophisticated matrix computations (like singular value decomposition and
  eigendecomposition) are powered by the Template Numerical Toolkit (TNT) ....not yet implemented


  References:
  -(1): Template Numerical Toolkit (TNT), http://math.nist.gov/tnt/

  */

template<class T>
class rsMatrixOld  // rename to rsMatrixOld
{

public:



  /** \name Construction/Destruction */

  /** Default constructor - constructs a matrix with zero rows, zero colums and NULL pointers. */
  rsMatrixOld();
    // not tested

  /** Constructor. You must pass the number of rows and colums here. */
  rsMatrixOld(int numRows, int numColumns, bool initElementsWithZeros = false);
    // not tested

  /** Constructor. You must pass the number of rows and colums here and the values to intialize
  the matix. */
  rsMatrixOld(int numRows, int numColumns, T **values);
    // not tested

  /** Copy constructor. It only copies the value of the pointer to the data and increases its
  reference-count. No deep coping is taking place because that would lead to excessive overhead
  in function and operator calls.
  \todo: maybe make protected
  */
  rsMatrixOld(const rsMatrixOld& other);


  /** Returns a square matrix with diagonal elements given by 'scalarValue' and zero off-diagonal
  elements. */
  static rsMatrixOld scalarMatrix(T scalarValue, int dimension)
  {
    rsMatrixOld result(dimension, dimension);
    for(int r=0; r<dimension; r++)  // use rsInitMatrix instead
    {
      for(int c=0; c<dimension; c++)
        result(r,c) = 0.0;
    }
    for(int r=0; r<dimension; r++)
      result(r,r) = scalarValue;
    return result;
  }
    // not tested

  /** Returns a square matrix with diagonal elements given by the array 'diagonalValues' and zero
  off-diagonal elements. The array is assumed to contain a number of elements equal to
  'dimension'. */
  static rsMatrixOld diagonalMatrix(T *diagonalValues, int dimension)
  {
    rsMatrixOld result(dimension, dimension);
    for(int r=0; r<dimension; r++)
    {
      for(int c=0; c<dimension; c++)
        result(r,c) = 0.0;
    }
    for(int r=0; r<dimension; r++)
      result(r,r) = diagonalValues[r];
    return result;
  }
  // untested

  /** Destructor. */
  ~rsMatrixOld();


  /** \name Manipulations */

  /** Sets up the size of the matrix (number of rows and columns). */
  void setSize(int newNumRows, int newNumColumns)
  {
    if(newNumRows == getNumRows() && newNumColumns == getNumColumns())
      return;
    else
    {
      data->numReferences -= 1;
      if(data->numReferences == 0)
        delete data;
      data = new rsMatrixData<T>(newNumRows, newNumColumns);
      data->numReferences = 1;
    }
  }
  // untested

  /** Sets one of the matrix elements. */
  void set(int row, int column, T newValue)
  {
    rsAssert(row < data->numRows && column < data->numColumns); // invalid row or column index

    if(data->numReferences > 1)
      makeDeepCopy();

    data->m[row][column] = newValue;
  }

  /** Transposes this matrix (exchanges the roles of rows and columns). */
  void transpose()
  {
    rsMatrixData<T>* newData = new rsMatrixData<T>(data->numColumns, data->numRows);
    for(int i = 0; i < data->numRows; i++)
    {
      for(int j = 0; j < data->numColumns; j++)
        newData->m[j][i] = data->m[i][j];
    }
    updateDataPointer(newData);
  }

  /** Assigns random values between 'min' and 'max' to each element. */
  void randomizeElements(T min, T max)
  {
    //if( data->numReferences > 1 )
    //  makeDeepCopy();

    makeDeepCopyIfNecessary();// it's actually not necessary to copy the old values - get rid of it

    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < data->numColumns; c++)
        data->m[r][c] = rsRandomUniform(min, max);
    }
  }
  // untested


  void copyDataFrom(const rsMatrixOld<T>& mat)
  {
    rsAssert(getNumRows()    == mat.getNumRows());
    rsAssert(getNumColumns() == mat.getNumColumns());
    for(int r = 0; r < data->numRows; r++)
      for(int c = 0; c < data->numColumns; c++)
        data->m[r][c] = mat.data->m[r][c];
  }

  /** Initializes all elements with zero values. */
  // \todo use rsInitMatrix
  void initWithZeros()
  {
    if(data->numReferences > 1)
      makeDeepCopy();  // it's actually not necessary to copy the old values - get rid of it

    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < data->numColumns; c++)
        data->m[r][c] = 0.0;
    }
  }
  // untested

  /** Initializes all elements with the sum of their row- and column indices. */
  // why is this useful?
  void initWithIndexSum()
  {
    if(data->numReferences > 1)
      makeDeepCopy();  // it's actually not necessary to copy the old values - get rid of it

    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < data->numColumns; c++)
        data->m[r][c] = r+c;
    }
  }
  // untested

  /** Applies the passed function to each element of the matrix. */
  void applyFunction(T (*f) (T))
  {
    if(data->numReferences > 1)
      makeDeepCopy(); // maybe we can avoid the copying by directly setting the values of the new
                      // matrix to the function values

    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < data->numColumns; c++)
        data->m[r][c] = f(data->m[r][c]);
    }
  }


  /** \name Inquiry */

  int getNumRows() const
  {
    return data->numRows;
  }

  int getNumColumns() const
  {
    return data->numColumns;
  }

  int getNumElements() const
  {
    return data->numRows * data->numColumns;
  }

  bool isSquare() const
  {
    return data->numRows == data->numColumns;
  }


  /** Returns a pointer to the actual underlying data.
  \todo: it's kinda dangerous and bad style to let client code access (and possibly modify) this
  data - find a better solution. Returning a const pointer somehow doesn't work in conjunction
  with rsSolveLinearSystem - check why and fix it. */
  T** getDataPointer() const
  {
    return data->m;
  }

  /** Returns a pointer to the row with geiven index. */
  T* getRowPointer(int index) const
  {
    return data->m[index];
  }

  // getRank, getConditionNumber, getTrace, getDeterminant, getFrobeniusNorm, getEuclideanNorm

  //rsMatrixOld getSubMatrix(int fromRow, int toRow, int fromColumn, int toColumn);
  //rsMatrixOld withRemovedRow(int rowToRemove);
  //rsMatrixOld withRemovedColumn(int columnToRemove);


  /** \name Operators */

  /** Accesses the element at given index-pair for reading and writing. */
  T& operator()(const int i, const int j) const
  {
    return data->m[i][j];
  }



  ///** Accesses the element at given index-pair for reading. */
  //const T& operator()(const int i, const int j) const
  //{
  //  return data->m[i][j];
  //}



  /** Assigns one matrix with another one. */
  rsMatrixOld& operator=(const rsMatrixOld& m2)
  {
    updateDataPointer(m2.data);
    return *this;
  }

  /** Compares two matrices for equality. */
  bool operator==(const rsMatrixOld& m2) const
  {
    if(data->numRows != m2.data->numRows || data->numColumns != m2.data->numColumns)
      return false;
    else
    {
      for(int r = 0; r < data->numRows; r++)
      {
        for(int c = 0; c < data->numColumns; c++)
        {
          if(data->m[r][c] != m2.data->m[r][c])
            return false;
        }
      }
      // use rsAreBuffersEqual
    }
    return true;
  }

  /** Compares two matrices for inequality. */
  bool operator!=(const rsMatrixOld& m2) const
  {
    return !(*this == m2);
  }

  /** Defines the negative of a matrix. */
  rsMatrixOld operator-()
  {
    rsMatrixOld result(data->numRows, data->numColumns);
    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < data->numColumns; c++)
        result.data->m[r][c] = -data->m[r][c];
    } // use rsNegateBuffer
    return result;
  }

  /** Adds another matrix to this matrix and returns the result. */
  rsMatrixOld& operator+=(const rsMatrixOld &m2)
  {
    rsAssert(data->numRows == m2.data->numRows && data->numColumns == m2.data->numColumns); // matrices incompatible
    makeDeepCopyIfNecessary();
    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < data->numColumns; c++)
        data->m[r][c] += m2.data->m[r][c];
    } // use rsAddBuffers, make deep copy
    return *this;
  }

  /** Adds a scalar to this matrix and returns the result. */
  rsMatrixOld& operator+=(const T &x)
  {
    makeDeepCopyIfNecessary();
    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < data->numColumns; c++)
        data->m[r][c] += x;
    }
    return *this;
  }

  /** Subtracts another matrix from this matrix and returns the result. */
  rsMatrixOld& operator-=(const rsMatrixOld &m2)
  {
    rsAssert(data->numRows == m2.data->numRows && data->numColumns == m2.data->numColumns); // matrices incompatible
    makeDeepCopyIfNecessary();
    rsArrayTools::subtract(data->mFlat, m2.data->mFlat, data->mFlat, getNumElements());
    return *this;

    /*
    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < data->numColumns; c++)
        data->m[r][c] -= m2.data->m[r][c];
    }
    return *this;
    */
  }

  /** Subtracts a scalar from this matrix and returns the result. */
  rsMatrixOld& operator-=(const T &x)
  {
    /*
    makeDeepCopyIfNecessary();
    rsSubtract(data->mFlat, getNumElements(), x);
    return *this;
    */

    makeDeepCopyIfNecessary();
    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < data->numColumns; c++)
        data->m[r][c] -= x;
    }
    return *this;
  }

  /** Right-multiplies this matrix with another matrix and returns the result. */
  rsMatrixOld& operator*=(const rsMatrixOld &m2)
  {
    rsAssert(data->numColumns == m2.data->numRows);  // matrices incompatible

    if(m2.isSquare() && data->numReferences == 1)
      MatrixTools::rsMatrixInPlaceMultiply(data->m, m2.data->m, data->numRows, data->numColumns);
    else
    {
      rsMatrixData<T> *newData = new rsMatrixData<T>(getNumRows(), m2.getNumColumns());

      MatrixTools::rsMatrixMultiply(data->m, m2.data->m, newData->m,
        data->numRows, data->numColumns, m2.data->numColumns);

      updateDataPointer(newData);

      /*
      data->numReferences -= 1;
      if( data->numReferences == 0 )
        delete data;
      data = newData;
      data->numReferences = 1;
        // maybe these 5 lines can be put into a function updateDataPointer(oldData, newData)
        */
    }

    return *this;
  }

  /** Multiplies this matrix by a scalar and returns the result. */
  rsMatrixOld& operator*=(const T &x)
  {
    makeDeepCopyIfNecessary();
    rsArrayTools::scale(data->mFlat, getNumElements(), x);
    return *this;
  }

  /** Divides this matrix by a scalar and returns the result. */
  rsMatrixOld& operator/=(const T &x)
  {
    makeDeepCopyIfNecessary();
    rsArrayTools::scale(data->mFlat, getNumElements(), T(1)/x);
    return *this;
  }

  /** Adds two matrices. */
  rsMatrixOld operator+(const rsMatrixOld &m2)
  {
    rsAssert(data->numRows == m2.data->numRows && data->numColumns == m2.data->numColumns); // matrices incompatible
    rsMatrixOld result(data->numRows, data->numColumns);
    rsArrayTools::add(data->mFlat, m2.data->mFlat, result.data->mFlat, getNumElements());
    return result;
  }

  /** Adds a matrix and a scalar. */
  rsMatrixOld operator+(const T &x)
  {
    rsMatrixOld result(data->numRows, data->numColumns);
    rsArrayTools::add(data->mFlat, x, result.data->mFlat, getNumElements());
    return result;
  }

  /** Subtracts two matrices. */
  rsMatrixOld operator-(const rsMatrixOld &m2)
  {
    rsAssert(data->numRows == m2.data->numRows && data->numColumns == m2.data->numColumns); // matrices incompatible
    rsMatrixOld result(data->numRows, data->numColumns);
    rsArrayTools::subtract(data->mFlat, m2.data->mFlat, result.data->mFlat, getNumElements());
    return result;
  }

  /** Subtracts a scalar from a matrix. */
  rsMatrixOld operator-(const T &x)
  {
    rsMatrixOld result(data->numRows, data->numColumns);
    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < data->numColumns; c++)
        result.data->m[r][c] = data->m[r][c] - x;
    }
    return result;
  }

  /** Multiplies a matrix and a scalar. */
  rsMatrixOld operator*(const T &x)
  {
    rsMatrixOld result(data->numRows, data->numColumns);
    rsArrayTools::scale(data->mFlat, result.data->mFlat, getNumElements(), x);
    return result;
  }

  /** Multiplies two matrices. */
  rsMatrixOld operator*(const rsMatrixOld &m2)
  {
    rsAssert(data->numColumns == m2.data->numRows);  // matrices incompatible
    rsMatrixOld result(data->numRows, m2.data->numColumns);
    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < m2.data->numColumns; c++)
      {
        result.data->m[r][c] = 0.0;
        for(int k = 0; k < data->numColumns; k++)
          result.data->m[r][c] += data->m[r][k] * m2.data->m[k][c];
      }
    } // use rsMatrixMultiply
    return result;
  }

  /** Divides a matrix by a scalar. */
  rsMatrixOld operator/(const T &x)
  {
    T scale = 1.0 / x;
    rsMatrixOld result(data->numRows, data->numColumns);
    for(int r = 0; r < data->numRows; r++)
    {
      for(int c = 0; c < data->numColumns; c++)
        result.data->m[r][c] = data->m[r][c] * scale;
    }
    return result;
  }


  /** \name Decompositions */

  /** Returns the the singular value decomposition of this matrix (lets denote it by A) such that
  A = U * S * V' where U and V are orthogonal matrices (U' * U = E, V' * V = E) and S is a
  diagonal matrix with positive or zero elements (the singular values). The singular values are
  ordered such that S[0][0] >= S[1][1] >= ...  */
  //void getSingularValueDecomposition(rsMatrixOld *U, rsMatrixOld *S, rsMatrixOld *V);

  // getEigenDecomposition, getLowerUpperDecomposition, ...


  //---------------------------------------------------------------------------------------------
  // others:




  /** Prints the values to the standard output - mainly for debugging purposes. */
  void print();



protected:


  /** Class for holding the actual data for the rsMatrixOld class. This has been factored out to
  facilitate sharing data between rsMatrixOld objects and, most importantly, to avoid excessive
  deep copying of data in assignment-operations and matrix-valued function return values by
  (re)assigning only pointers.

  \todo - move the implementation into the .inl file

  */

  //template<class TT>
  template<class>
  class rsMatrixData
  {
    rsMatrixData(int numRows, int numColumns)
    {
      this->numRows    = numRows;
      this->numColumns = numColumns;
      numReferences    = 0;
      m                = nullptr;
      mFlat            = nullptr;
      allocateMemory();
    }

    ~rsMatrixData()
    {
      freeMemory();
    }

    void allocateMemory()
    {
      if(numRows > 0 && numColumns > 0)
      {
        mFlat = new T[numRows*numColumns];
        m     = new T*[numRows];
        for(int i = 0; i < numRows; i++)
          m[i] = &mFlat[i*numColumns];
      }
    }

    void freeMemory()
    {
      delete[] m;
      delete[] mFlat;
      m     = nullptr;
      mFlat = nullptr;
    }

    rsMatrixData<T>* getDeepCopy()
    {
      rsMatrixData<T>* copy = new rsMatrixData<T>(numRows, numColumns);
      memcpy(mFlat, copy->mFlat, numRows*numColumns*sizeof(T));
      for(int i = 0; i < numRows; i++)
        copy->m[i] = &copy->mFlat[i*numColumns];
      return copy;
    }

    //void addReference();
    //void removeReference

    int numRows;
    int numColumns;
    int numReferences;
     // \todo maybe use rsUint32 instead of int

    T **m;     // 2D pointer to the values: 1st index is the row, 2nd index the column.
    T *mFlat;  // matrix as flat array - to facilitate contiguous memory allocation

    friend class rsMatrixOld;
  };


  /** Creates a deep copy of the actual data, re-adjusts our pointer to the copy and decrements
  the reference counter of the original data. */
  void makeDeepCopy();

  RS_INLINE void makeDeepCopyIfNecessary();

  void updateDataPointer(rsMatrixData<T> *newData);


  rsMatrixData<T> *data;


  /** Copy constructor. It should not be used by client code but is needed internally for
  matrix-valued return values of functions and operators. It only copies the value of the pointer
  to the data and increases its reference-count. No deep coping is taking place because that
  would lead to excessive overhead in function and operator calls. */
  //rsMatrixOld(const rsMatrixOld& other);
    // untested

  /*
  // friend functions:
  friend RS_INLINE rsMatrixOld<T> trans(const rsMatrixOld<T>& A);
  friend RS_INLINE rsMatrixOld<T> outerProduct(const rsVector<T> &a, const rsVector<T> &b);

  // friend operators:
  friend RS_INLINE rsMatrixOld<T> operator+(const T &x, const rsMatrixOld<T> &m);
  friend RS_INLINE rsMatrixOld<T> operator-(const T &x, const rsMatrixOld<T> &m);
  friend RS_INLINE rsMatrixOld<T> operator*(const T &x, const rsMatrixOld<T> &m);
  friend RS_INLINE rsVector<T> operator*(const rsMatrixOld<T> &A, const rsVector<T> &b);
  friend RS_INLINE rsMatrixOld<T> operator*(const rsVector<T> &b, const rsMatrixOld<T> &A);
  */


}; // end of class rsMatrixOld

// some binary operators are defined outside the class such that the left hand operand does
// not necesarrily need to be of class rsMatrixOld
// \todo: maybe get rid of the inlining

/** Adds a scalar and a matrix. */
template<class T>
RS_INLINE rsMatrixOld<T> operator+(const T &x, const rsMatrixOld<T> &m)
{
  rsMatrixOld<T> result(m.getNumRows(), m.getNumColumns());
  for(int r = 0; r < m.getNumRows(); r++)
  {
    for(int c = 0; c < m.getNumColumns(); c++)
    {
      //result.m[r][c] = m.m[r][c] + x;
      result.set(r, c, m(r, c) + x); // maybe we need setUnsafe
    }
  }
  return result;
}

/** Subtracts a matrix from a scalar. */
template<class T>
RS_INLINE rsMatrixOld<T> operator-(const T &x, const rsMatrixOld<T> &m)
{
  rsMatrixOld<T> result(m.getNumRows(), m.getNumColumns());
  for(int r = 0; r < m.getNumRows(); r++)
  {
    for(int c = 0; c < m.getNumColumns(); c++)
    {
      //result.m[r][c] = x - m.m[r][c];
      result.set(r, c, x - m(r, c)); // maybe we need setUnsafe
    }
  }
  return result;
}

/** Multiplies a scalar and a matrix. */
template<class T>
RS_INLINE rsMatrixOld<T> operator*(const T &x, const rsMatrixOld<T> &m)
{
  rsMatrixOld<T> result(m.getNumRows(), m.getNumColumns());
  for(int r = 0; r < m.getNumRows(); r++)
  {
    for(int c = 0; c < m.getNumColumns(); c++)
    {
      //result.m[r][c] = x * m.m[r][c];
      result.set(r, c, x * m(r, c));  // use setUnsafe, rename m to A
    }
  }
  return result;

  /*
  rsMatrixOld<T> result(m.numRows, m.numColumns);
  for(int r=0; r<m.numRows; r++)
  {
    for(int c=0; c<m.numColumns; c++)
      result.m[r][c] = x * m.m[r][c];
  }
  return result;
  */
}

/** Returns matrix A transposed. */
template<class T>
RS_INLINE rsMatrixOld<T> trans(const rsMatrixOld<T>& A)
{
  rsMatrixOld<T> B(A.getNumColumns(), A.getNumRows());
  for(int r = 0; r < B.getNumRows(); r++)
  {
    for(int c = 0; c < B.getNumColumns(); c++)
    {
      B.set(r, c, A(c, r));  // use setUnsafe
      //B.m[r][c] = A.m[c][r];
    }
  }
  return B;


  /*
  rsMatrixOld<T> B(A.numColumns, A.numRows);
  for(int r=0; r<B.numRows; r++)
  {
    for(int c=0; c<B.numColumns; c++)
      B.m[r][c] = A.m[c][r];
  }
  return B;
  */
}



// matrix/vector functions:


/** Left-multiplies a vector by a matrix such that c = A*b with A being the input matrix, b being
the input vector and c being the output vector. Thus, b is interpreted as a column-vector and the
number of columns in the matrix A must equal the number of elements in the vector b. The result
will be a vector with the number rows of A as its dimensionality. */
template<class T>
RS_INLINE rsVector<T> operator*(const rsMatrixOld<T> &A, const rsVector<T> &b)
{
  rsAssert(A.getNumColumns() == b.dim);  // matrix and vector incompatible
  rsVector<T> c(A.getNumRows());
  for(int i = 0; i < A.getNumRows(); i++)
  {
    c.v[i] = 0.0;
    for(int j = 0; j < A.getNumColumns(); j++)
      c.v[i] += A(i, j) * b.v[j];
  }
  return c;


  /*
  rsAssert( A.data->numColumns == b.dim );  // matrix and vector incompatible
  rsVector<T> c(A.data->numRows);
  for(int i = 0; i < A.data->numRows; i++)
  {
    c.v[i] = 0.0;
    for(int j = 0; j < A.data->numColumns; j++)
      c.v[i] += A.data->m[i][j] * b.v[j];
  }
  return c;
  */
}
// untested

/** Right-multiplies a transposed vector by a matrix such that c = b^T * A with A being the input
matrix, b being the input vector and c being the output vector. Thus, b is interpreted as a
row-vector and the number of rows in the matrix A must equal the number of elements in the
vector b. The result will be a matrix consisting of one row and a number of columns equal to the
dimenstionality of the vector b. */
template<class T>
RS_INLINE rsMatrixOld<T> operator*(const rsVector<T> &b, const rsMatrixOld<T> &A)
{
  rsAssert(A.numRows == b.dim);  // matrix and vector incompatible
  rsMatrixOld<T> C(1, b.dim);
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
RS_INLINE rsMatrixOld<T> outerProduct(const rsVector<T> &a, const rsVector<T> &b)
{
  rsAssert(a.dim == b.dim);   // vectors incopatible
  rsMatrixOld<T> result(a.dim, a.dim);
  for(int i=0; i<result.numRows; i++)
  {
    for(int j=0; j<result.numColumns; j++)
      result.m[i][j] = a.v[i] * b.v[j];
  }
  return result;
}
    // not tested

/*
template<class T>
rsMatrix<T> matrixMagnitudes(const rsMatrix<std::complex<T>>& A)
{
  int N = A.getNumRows();
  int M = A.getNumColumns();
  rsMatrix<T> mags(N, M);
  for(int i = 0; i < N; i++)
    for(int j = 0; j < M; j++)
      mags(i, j) = abs(A(i, j));
  return mags;
}
template<class T>
rsMatrix<T> matrixPhases(const rsMatrix<std::complex<T>>& A)
{
  int N = A.getNumRows();
  int M = A.getNumColumns();
  rsMatrix<T> phases(N, M);
  for(int i = 0; i < N; i++)
    for(int j = 0; j < M; j++)
      phases(i, j) = arg(A(i, j));
  return phases;
}
// maybe factor out common code...maybe something like applyMatrixFunction with different
// input and output types for the template parameter
*/


// typedefs for explicit instantiations:
typedef rsMatrixOld<double> rsMatrixDbl;


/*
Note:
There was a design decision to be made with regard to when we want to invoke deep copying and
when only shallow copying. There are 3 possibilities worth to consider: (1) deep copying in
copy-constructor and assignment operator, (2) deep copying in assignment operator, shallow
copying in copy-constructor, (3) shallow copying in assignment operator and copy-constructor. In
scenario (1), we would be free to let client code manipulate matrix elements at will without
checking any reference counter, because each matrix object would have its own unique data. In
scenario (2), we could make the copy-constructor unaccessible from client code which ensures that
it will only be used for function/operator return values, so an operator call as in A+B would
avoid copying but an assigment afterwards as in C = A+B would still copy data. We could still
let client code manipulate matrix-data at will without checking a reference counter. In scenario
(3), we would avoid as much copying as possible, but after a statement like B = A, A and B would
reference the same data internally, so in order to let client code modify B but leave A
untouched, we would have to check a reference counter in any attempt to modify B (or A) and if
its > 1, we would need to make a deep copy at the moment of the manipluation. Since it is
supposed to be much more common to have more assignments involving matrix-objects than to
manipulate matrix elements directly, it has been opted for variant (3) - so operator/function
returns and assignments will be fast, but write access for single elements (via the set-function)
might be somewhat slow.

Maybe provide a setUnsafe function which avoids reference-counter checking (and possibly the
creation of a deep copy) and is therefore fast but unsafe. Client code which uses this function
should be sure, that the side effect of altering other matrix-objects which reference  the same
data may occur.

This strategy could be applied to other data-heavy classes as well. It could be considered a
pattern, maybe "Lazy Deep Copy" pattern or something. Maybe it already exists - look up the
literature (maybe "Lazy Initialization")

*/

#endif
