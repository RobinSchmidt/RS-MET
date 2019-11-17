#ifndef rosic_Matrix_h
#define rosic_Matrix_h

// deprecated - get rid of this class - use RAPT::rsMatrix<double> instead

namespace rosic
{

  /**

  This is a class for representing matrices and doing mathematical operations with them. To make
  the class fast, no consistency checking is done in the mathematical operators. You must ensure
  yourself that the input arguments are compatible - for example don't try to add two matrices of
  different dimensionality or multiply an n times m matrix with a q times p matrix when m != q.
  The more sophisticated matrix computations (like singular value decomposition and
  eigendecomposition) are powered by the Template Numerical Toolkit (TNT).

  References:
  -(1): Template Numerical Toolkit (TNT), http://math.nist.gov/tnt/

  */

  class Matrix
  {

  public:

    //---------------------------------------------------------------------------------------------
    // public member variables:

    /** Number of rows. */
    int numRows;

    /** Number of colums. */
    int numColumns;

    /** The actual values of the matrix. The first index refers to the row, the second to the
    column. */
    double **m;

    /** The matrix as flat array - this is only used to facilitate contiguous memory allocation. */
    double *mFlat;

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Default constructor - constructs a matrix with zero rows, zero colums and NULL pointers. */
    Matrix();
    // check

    /** Constructor. You must pass the number of rows and colums here. */
    Matrix(int numRows, int numColumns, bool initElementsWithZeros = false);
    // check

    /** Constructor. You must pass the number of rows and colums here and the values to intialize
    the matix. */
    Matrix(int numRows, int numColumns, double **values);

    /** Copy constructor. */
    Matrix(const Matrix& other);
    // check

    /** Returns a square matrix with diagonal elements given by 'scalarValue' and zero off-diagonal
    elements. */
    static Matrix scalarMatrix(double scalarValue, int dimension)
    {
      Matrix result(dimension, dimension);
      for(int r=0; r<dimension; r++)
      {
        for(int c=0; c<dimension; c++)
          result.m[r][c] = 0.0;
      }
      for(int r=0; r<dimension; r++)
        result.m[r][r] = scalarValue;
      return result;
    }
    // check

    /** Returns a square matrix with diagonal elements given by the array 'diagonalValues' and zero
    off-diagonal elements. The array is assumed to contain a number of elements equal to
    'dimension'. */
    static Matrix diagonalMatrix(double *diagonalValues, int dimension)
    {
      Matrix result(dimension, dimension);
      for(int r=0; r<dimension; r++)
      {
        for(int c=0; c<dimension; c++)
          result.m[r][c] = 0.0;
      }
      for(int r=0; r<dimension; r++)
        result.m[r][r] = diagonalValues[r];
      return result;
    }
    // check

    /** Destructor. */
    ~Matrix();

    //---------------------------------------------------------------------------------------------
    // manipulations:

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
        mFlat      = new double[numRows*numColumns];
        m          = new double*[numRows];
        for(int r=0; r<numRows; r++)
          m[r] = &mFlat[r*numColumns];
      }
    }

    /** Sets one of the matrix elements. */
    void setElement(int row, int column, double newValue)
    {
      rassert( row < numRows && column < numColumns ); // invalid row or column index
      m[row][column] = newValue;
    }
    // check

    /** Transposes this matrix (exchanges the roles of rows and columns). */
    void transpose()
    {
      // treat square matrices seperately using temporary stack-memory (via alloca) to avoid heap
      // fragmentation - nope - don't do it like that - todo: do it in place!
      if( numRows == numColumns )
      {
        //double *tmp = (double*) alloca(numRows*numColumns*sizeof(double));
        double *tmp = new double[numRows*numColumns]; // should not be needed at all!
        memcpy(tmp, mFlat, numRows*numColumns*sizeof(double));
        for(int r=0; r<numRows; r++)
        {
          for(int c=0; c<numColumns; c++)
            m[r][c] = tmp[c*numRows+r];   //m[r][c] = tmp[c][r];
        }
        delete[] tmp;
      }
      else
      {
        // create a new 2D-array:
        int numRowsNew    = numColumns;
        int numColumnsNew = numRows;
        double *mFlatNew  = new double[numRowsNew*numColumnsNew];
        double **mNew     = new double*[numRowsNew];
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
        m         = mNew;
        mFlat     = mFlatNew;
        numRows    = numRowsNew;
        numColumns = numColumnsNew;
      }
    }
    // check

    /** Assigns random values between 'min' and 'max' to each element. */
    void randomizeElements(double min, double max)
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] = RAPT::rsRandomUniform(min, max);
      }
    }

    /** Initializes all elements with zero values. */
    void initWithZeros()
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] = 0.0;
      }
    }

    /** Initializes all elements with the sum of their row- and column indices. */
    void initWithIndexSum()
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] = r+c;
      }
    }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    //Matrix getSubMatrix(int fromRow, int toRow, int fromColumn, int toColumn);

    //---------------------------------------------------------------------------------------------
    // overloaded operators:

    /** Assigns one matrix with another one. */
    Matrix& operator=(const Matrix& m2)
    {
      // re-allocate memory, if necessary:
      if( this->numRows != m2.numRows || this->numColumns != m2.numColumns )
      {
        delete[] m;
        delete[] mFlat;
        numRows    = m2.numRows;
        numColumns = m2.numColumns;
        mFlat      = new double[numRows*numColumns];
        m          = new double*[numRows];
        for(int r=0; r<numRows; r++)
          m[r] = &mFlat[r*numColumns];
      }

      // copy the values:
      memcpy(mFlat, m2.mFlat, numRows*numColumns*sizeof(double));

      return *this;
    }
    // check

    /** Compares two matrices for equality. */
    bool operator==(const Matrix& m2) const
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
    // check

    /** Compares two matrices for inequality. */
    bool operator!=(const Matrix& m2) const
    {
      return !(*this == m2);
    }
    // check

    /** Defines the negative of a matrix. */
    Matrix operator-()
    {
      Matrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = -m[r][c];
      }
      return result;
    }
    // check

    /** Adds another matrix to this matrix and returns the result. */
    Matrix& operator+=(const Matrix &m2)
    {
      rassert( numRows == m2.numRows && numColumns == m2.numColumns ); // matrices incompatible
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] += m2.m[r][c];
      }
      return *this;
    }
    // check

    /** Adds a scalar to this matrix and returns the result. */
    Matrix& operator+=(const double &x)
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] += x;
      }
      return *this;
    }
    // check

    /** Subtracts another matrix from this matrix and returns the result. */
    Matrix& operator-=(const Matrix &m2)
    {
      rassert( numRows == m2.numRows && numColumns == m2.numColumns ); // matrices incompatible
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] -= m2.m[r][c];
      }
      return *this;
    }
    // check

    /** Subtracts a scalar from this matrix and returns the result. */
    Matrix& operator-=(const double &x)
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] -= x;
      }
      return *this;
    }
    // check

    /** Multiplies this matrix by a scalar and returns the result. */
    Matrix& operator*=(const double &x)
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] *= x;
      }
      return *this;
    }
    // check

    /** Divides this matrix by a scalar and returns the result. */
    Matrix& operator/=(const double &x)
    {
      double scale = 1.0 / x;
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] *= scale;
      }
      return *this;
    }
    // check

    /** Adds two matrices. */
    Matrix operator+(const Matrix &m2)
    {
      rassert( numRows == m2.numRows && numColumns == m2.numColumns ); // matrices incompatible
      Matrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] + m2.m[r][c];
      }
      return result;
    }
    // check

    /** Adds a matrix and a scalar. */
    Matrix operator+(const double &x)
    {
      Matrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] + x;
      }
      return result;
    }
    // check

    /** Subtracts two matrices. */
    Matrix operator-(const Matrix &m2)
    {
      rassert( numRows == m2.numRows && numColumns == m2.numColumns ); // matrices incompatible
      Matrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] - m2.m[r][c];
      }
      return result;
    }
    // check

    /** Subtracts a scalar from a matrix. */
    Matrix operator-(const double &x)
    {
      Matrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] - x;
      }
      return result;
    }
    // check

    /** Multiplies a matrix and a scalar. */
    Matrix operator*(const double &x)
    {
      Matrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] * x;
      }
      return result;
    }
    // check

    /** Multiplies two matrices. */
    Matrix operator*(const Matrix &m2)
    {
      rassert( this->numColumns == m2.numRows )  // matrices incompatible
      Matrix result(numRows, m2.numColumns);
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
    // check

    /** Divides a matrix by a scalar. */
    Matrix operator/(const double &x)
    {
      double scale = 1.0 / x;
      Matrix result(numRows, numColumns);
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          result.m[r][c] = m[r][c] * scale;
      }
      return result;
    }
    // check

    //---------------------------------------------------------------------------------------------
    // setup:

    // removeRow, removeColumn, ...

    /** Applies the passed function to each element of the matrix. */
    void applyFunction( double (*f) (double) )
    {
      for(int r=0; r<numRows; r++)
      {
        for(int c=0; c<numColumns; c++)
          m[r][c] = f(m[r][c]);
      }
    }
    // check

    //---------------------------------------------------------------------------------------------
    // inquiry:

    bool isSquare() const { return numRows == numColumns; }
    // check


    //---------------------------------------------------------------------------------------------
    // matrix computations:

    ///** Returns the the singular value decomposition of this matrix (lets denote it by A) such that
    //A = U * S * V' where U and V are orthogonal matrices (U' * U = E, V' * V = E) and S is a
    //diagonal matrix with positive or zero elements (the singular values). The singular values are
    //ordered such that S[0][0] >= S[1][1] >= ...  */
    //void getSingularValueDecomposition(Matrix *U, Matrix *S, Matrix *V);


    // getDeterminant, ...

    //---------------------------------------------------------------------------------------------
    // others:

    /** Prints the values to the standard output - mainly for debugging purposes. */
    void print();

  }; // end of class Matrix

  // some binary operators are defined outside the class such that the left hand operand does
  // not necesarrily need to be of class Matrix

  /** Adds a scalar and a matrix. */
  inline Matrix operator+(const double &x, const Matrix &m)
  {
    Matrix result(m.numRows, m.numColumns);
    for(int r=0; r<m.numRows; r++)
    {
      for(int c=0; c<m.numColumns; c++)
        result.m[r][c] = m.m[r][c] + x;
    }
    return result;
  }
  // check

  /** Subtracts a matrix from a scalar. */
  inline Matrix operator-(const double &x, const Matrix &m)
  {
    Matrix result(m.numRows, m.numColumns);
    for(int r=0; r<m.numRows; r++)
    {
      for(int c=0; c<m.numColumns; c++)
        result.m[r][c] = x - m.m[r][c];
    }
    return result;
  }
  // check

  /** Multiplies a scalar and a matrix. */
  inline Matrix operator*(const double &x, const Matrix &m)
  {
    Matrix result(m.numRows, m.numColumns);
    for(int r=0; r<m.numRows; r++)
    {
      for(int c=0; c<m.numColumns; c++)
        result.m[r][c] = x * m.m[r][c];
    }
    return result;
  }
  // check

  /** Returns matrix A transposed. */
  inline Matrix trans(const Matrix& A)
  {
    Matrix B(A.numColumns, A.numRows);
    for(int r=0; r<B.numRows; r++)
    {
      for(int c=0; c<B.numColumns; c++)
        B.m[r][c] = A.m[c][r];
    }
    return B;
  }
  // check


}  // end namespace rosic

#endif // rosic_Matrix_h
