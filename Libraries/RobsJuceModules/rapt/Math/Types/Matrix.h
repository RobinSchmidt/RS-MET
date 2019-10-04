#ifndef RAPT_MATRIX_H
#define RAPT_MATRIX_H



template<class T>
class rsMatrix2x2
{

public:

  T a, b, c, d;
  // matrix coefficients |a b|
  //                     |c d|


  /** Stadard constructor. You can pass the matrix elements. If you pass nothing, an identity
  matrix will be created. */
  rsMatrix2x2(T a = T(1), T b = T(0), T c = T(0), T d = T(1)) { setValues(a, b, c, d); }
  // todo: maybe require arguments to be passed - or initialze teh matrix to the zero matrix
  // by default


  /** \name Setup */

  void setValues(T a, T b, T c, T d) { this->a = a; this->b = b; this->c = c; this->d = d; }




  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the determinant of this matrix. */
  T determinant() const { return a*d - b*c; }

  /** Returns the trace (sum of diagonal elements) of this matrix.  */
  T trace() const { return a+d; }

  /** Returns the first eigenvalue of this matrix. */
  T eigenvalue1() const { return rsLinearAlgebra::eigenvalue2x2_1(a, b, c, d); }

  /** Returns the second eigenvalue of this matrix. */
  T eigenvalue2() const { return rsLinearAlgebra::eigenvalue2x2_2(a, b, c, d); }

  /** Returns the first eigenvector of this matrix. */
  rsVector2D<T> eigenvector1() const
  { rsVector2D<T> v; rsLinearAlgebra::eigenvector2x2_1(a, b, c, d, &v.x, &v.y, true); return v; }

  /** Returns the second eigenvector of this matrix. */
  rsVector2D<T> eigenvector2() const
  { rsVector2D<T> v; rsLinearAlgebra::eigenvector2x2_2(a, b, c, d, &v.x, &v.y, true); return v; }

  /** Returns the inverse of this matrix. */
  rsMatrix2x2<T> inverse() const
  { T D = determinant(); T s = T(1) / D; return rsMatrix2x2<T>(s*d, -s*b, -s*c, s*a); }

  // maybe these functions should be named getDeterminant, etc. - more consistent with other
  // classes and states more explicitly what they do

  /** Tests, if another matrix B is close to this matrix within a given tolerance (all components
  of the difference must be <= tolerance). */
  bool isCloseTo(const rsMatrix2x2<T>& B, T tol)
  {
    if(rsAbs(a-B.a) <= tol && rsAbs(b-B.b) <= tol && rsAbs(c-B.c) <= tol && rsAbs(d-B.d) <= tol)
      return true;
    return false;
  }
  // doesn't work with complex matrices

  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Adds two matrices: C = A + B. */
  rsMatrix2x2<T> operator+(const rsMatrix2x2<T>& B) const
  { rsMatrix2x2<T> C; C.a = a + B.a; C.b = b + B.b; C.c = c + B.c; C.d = d + B.d; return C; }

  /** Subtracts two matrices: C = A - B. */
  rsMatrix2x2<T> operator-(const rsMatrix2x2<T>& B) const
  { rsMatrix2x2<T> C; C.a = a - B.a; C.b = b - B.b; C.c = c - B.c; C.d = d - B.d; return C; }

  /** Multiplies two matrices: C = A * B. */
  rsMatrix2x2<T> operator*(const rsMatrix2x2<T>& B) const
  { rsMatrix2x2<T> C; C.a = a*B.a + b*B.c; C.b = a*B.b + b*B.d;
    C.c = c*B.a + d*B.c; C.d = c*B.b + d*B.d; return C; }

  /** Multiplies the left matrix operand with the inverse of the right matrix operand. */
  rsMatrix2x2<T> operator/(const rsMatrix2x2<T>& B) const { return *this * B.inverse(); }

  /** Compares matrices for equality */
  bool operator==(const rsMatrix2x2<T>& B) const
  { return a == B.a && b == B.b && c == B.c && d == B.d; }

  /** Multiplies matrix by a vector: w = A*v */
  rsVector2D<T> operator*(const rsVector2D<T>& v) const
  {
    rsVector2D<T> w;
    w.x = a * v.x  +  b * v.y;
    w.y = c * v.x  +  d * v.y;
    return w;
  }
  // todo: left multiplication w = v^H * A

  // todo: operators that take a scalar as left or right argument


  //-----------------------------------------------------------------------------------------------
  /** \name Factory */

  static rsMatrix2x2<T> zero()     { return rsMatrix2x2<T>(T(0), T(0), T(0), T(0)); }

  static rsMatrix2x2<T> identity() { return rsMatrix2x2<T>(T(1), T(0), T(0), T(1)); }


  /** Returns the commutator of the two matrices A and B: C = A*B - B*A. In general, matrix
  multiplication is non-commutative, but for some special cases, it may be commutative nonetheless.
  The commutator captures, how non-commutative two matrices behave when being multiplied. If the
  two matrices commute (i.e. behave commutatively), their commutator is the zero matrix. */
  static rsMatrix2x2<T> commutator(const rsMatrix2x2<T>& A, const rsMatrix2x2<T>& B)
  {
    return A*B - B*A;
  }

};

/** Multiplies a scalar and a matrix. */
template<class T>
inline rsMatrix2x2<T> operator*(const T& s, const rsMatrix2x2<T>& A)
{
  return rsMatrix2x2<T>(s*A.a, s*A.b, s*A.c, s*A.d);
}


//=================================================================================================

/** This is a class for treating C-arrays as matrices. It does not store/own the actual matrix
data, it just acts as wrapper around an existing array for more conveniently accessing and
manipulating matrix elements. */

template<class T>
class rsMatrixView
{

public:

  /** \name Construction/Destruction */

  /**  */
  rsMatrixView(int numRows = 0, int numColumns = 0, T* data = nullptr)
  {
    this->numRows = numRows;
    this->numCols = numColumns;
    d = data;
  }

  /** \name Setup */

  inline void setAllValues(T value) { rsArray::fillWithValue(d, getSize(), value); }

  inline void scale(T factor) { rsArray::scale(d, getSize(), factor); }


  inline void reshape(int newNumRows, int newNumColumns)
  {
    rsAssert(newNumRows*newNumColumns == numRows*numCols);
    numRows = newNumRows;
    numCols = newNumColumns;
  }

  // void setToIdentityMatrix(T scaler = 1);



  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /**  */
  static bool areSameShape(const rsMatrixView<T>& A, const rsMatrixView<T>& B)
  {
    return A.numRows == B.numRows && A.numCols == B.numCols;
  }

  int getNumRows()    const { return numRows; }

  int getNumColumns() const { return numCols; }

  int getSize()       const { return numRows * numCols; }



  //-----------------------------------------------------------------------------------------------
  /** \name Arithmetic */

  /** Adds elements of A to corresponding elements in B and stores results in C. */
  static void add(const rsMatrixView<T>* A, const rsMatrixView<T>* B, rsMatrixView<T>* C)
  {
    rsAssert(areSameShape(*A, *B) && areSameShape(*A, *C), "arguments incompatible");
    for(int i = 0; i < A->numRows; i++)
      for(int j = 0; j < A->numCols; j++)
        (*C)(i, j) = A->at(i, j) + B->at(i, j);
  }

  /** Subtracts elements of B from corresponding elements A in and stores results in C. */
  static void sub(const rsMatrixView<T>* A, const rsMatrixView<T>* B, rsMatrixView<T>* C)
  {
    rsAssert(areSameShape(*A, *B) && areSameShape(*A, *C), "arguments incompatible");
    for(int i = 0; i < A->numRows; i++)
      for(int j = 0; j < A->numCols; j++)
        (*C)(i, j) = A->at(i, j) - B->at(i, j);
  }

  /** Computes the matrix product C = A*B. */
  static void mul(const rsMatrixView<T>* A, const rsMatrixView<T>* B, rsMatrixView<T>* C)
  {
    // maybe factor out a function: areMultiplicable(A, B)
    rsAssert(A->numCols == B->numRows);
    rsAssert(C->numCols == B->numCols);
    rsAssert(C->numRows == A->numRows);
    for(int i = 0; i < C->numRows; i++) {
      for(int j = 0; j < C->numCols; j++) {
        (*C)(i,j) = T(0);
        for(int k = 0; k < A->numCols; k++)
          (*C)(i,j) += A->at(i,k) * B->at(k,j); }}
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Read and write access to matrix elements with row-index i and column-index j. */
  inline T& operator()(const int i, const int j) { return d[flatIndex(i, j)]; }

  /** Read only accees - used mainly internally with const reference arguments (for example,
  in add). */
  inline const T& at(const int i, const int j) const { return d[flatIndex(i, j)]; }

  /** Converts a row index i and a column index j to a flat array index. */
  inline int flatIndex(const int i, const int j) const
  {
    return numCols*i + j;
    // todo:
    //  -be more general: colStride*i + rowStride*j. goal: allow row-major and column-major storage
    //   while the syntax of the operator is always row-major (as is conventional in math)
    //   regardless whatever the internal storage format is - column major storage is required for
    //   compatibility with lapack
    // -maybe be even more general: colOffset + colStride*i + (rowOffset + rowStride)*j
    //  -> may allow to access sub-matrices with the same syntax (todo: verify formula)
  }


protected:

  /** \name Data */

  //size_t N, M;    // number of rows and columns

  int numRows, numCols;
  T *d;           // data pointer

};

//=================================================================================================

/** This is a class for representing matrices and doing mathematical operations with them. */

template<class T>
class rsMatrixNew : public rsMatrixView<T>
{

public:

  /** \name Construction/Destruction */


  /** Standard constructor. You must pass the initial number of rows and columns */
  rsMatrixNew(int numRows = 0, int numColumns = 0)
  {
    setSize(numRows, numColumns);
    // todo: optionally init with zeros
  }

  /** Creates matrix from a std::vector - convenient to initialize elements.  */
  rsMatrixNew(int numRows, int numColumns, const std::vector<T>& newData)
    : data(newData)
  {
    rsAssert(numRows*numColumns == newData.size());
    this->numRows = numRows;
    this->numCols = numColumns;
    updateDataPointer();
  }

  /** Copy constructor. */
  rsMatrixNew(const rsMatrixNew& B)
  {
    setSize(B.numRows, B.numCols);
    rsArray::copyBuffer(B.d, this->d, this->getSize());
  }

  /** Move constructor. */
  rsMatrixNew(const rsMatrixNew&& B)
  {
    setSize(B.numRows, B.numCols);
    rsArray::copyBuffer(B.d, this->d, this->getSize());
  }

  rsMatrixNew<T>& operator=(const rsMatrixNew<T>& other) // copy assignment
  {
    if (this != &other) { // self-assignment check expected
      setSize(other.numRows, other.numCols);
      rsArray::copyBuffer(other.d, this->d, this->getSize());
    }
    return *this;
  }

  rsMatrixNew<T>& operator=(const rsMatrixNew<T>&& other) // move assignment
  {
    if (this != &other) { // self-assignment check expected
      setSize(other.numRows, other.numCols);
      rsArray::copyBuffer(other.d, this->d, this->getSize());
    }
    return *this;
  }
  // can we avoid the copying? ..i mean, that's the wohle point of move operators



  // todo: implement the various copy/move assigment operators and -constructors - this should
  // optimize returning values from functions and operators (avoid unnessary copying)
  // https://en.cppreference.com/w/cpp/language/operators#Assignment_operator
  // ..i tried - but no avail yet
  // but we have to implement them because the standard versions will copy the pointer variable
  // inherited from the baseclass - we must call updateDataPointer in the copy/move
  // construtors/assigners


  /** \name Setup */

  /** Sets the number of rows and columns, this matrix should have. ToDo: provide a way to retain
  the data (optionally) - what does std::vector's resize do? Does it retain data...but if it does,
  it would be useless anyway in case the number of columns changed. */
  void setSize(int numRows, int numColumns)
  {
    this->numRows = numRows;
    this->numCols = numColumns;
    data.resize(this->numRows * this->numCols);
    updateDataPointer();
    // optionally initialize with zeros
  }


  /** Computes the Kronecker product between matrices A and B. For a 3x2 matrix A, it looks like
  that:
              |a11*B a12*B|
  A (x) B  =  |a21*B a22*B|
              |a31*B a32*B|

  Where each entry aij*B is a submatrix of dimensions of B with the entries of b scaled by an
  appropriate element from A. */
  static rsMatrixNew<T> kroneckerProduct(const rsMatrixNew<T>& A, const rsMatrixNew<T>& B)
  {
    rsMatrixNew<T> C(A.numRows*B.numRows, A.numCols*B.numCols);
    for(int ia = 0; ia < A.numRows; ia++) {
      for(int ja = 0; ja < A.numCols; ja++) {
        int startRow = ia*B.numRows;
        int startCol = ja*B.numCols;
        for(int ib = 0; ib < B.numRows; ib++) {
          for(int jb = 0; jb < B.numCols; jb++) {
            C(startRow+ib, startCol+jb) = A.at(ia,ja) * B.at(ib, jb); }}}}
    return C;
  }
  // see https://rosettacode.org/wiki/Kronecker_product#C



  /** \name Manipulations */


  /** \name Inquiry */


  /** \name Decompositions */


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Compares matrices for equality */
  bool operator==(const rsMatrixNew<T>& B) const
  {
    if(this->numRows != B.numRows || this->numCols != B.numCols)
      return false;
    return rsArray::areBuffersEqual(this->d, B.d, this->getSize());
  }
  // move to rsMatrixView

  /** Compares matrices for inequality */
  bool operator!=(const rsMatrixNew<T>& B) const { return !(*this == B); }
  // move to rsMatrixView


  /** Defines the negative of a matrix. */
  rsMatrixNew<T> operator-()
  {
    rsMatrixNew<T> C(this->numRows, this->numCols);
    for(int i = 0; i < this->getSize(); i++)
      C.d[i] = -d[i]; // maybe factor out into "neg" function in baseclass
    return C;
  }


  /** Adds two matrices: C = A + B. */
  rsMatrixNew<T> operator+(const rsMatrixNew<T>& B) const
  { rsMatrixNew<T> C(this->numRows, this->numCols); this->add(this, &B, &C); return C; }

  /** Subtracts two matrices: C = A - B. */
  rsMatrixNew<T> operator-(const rsMatrixNew<T>& B) const
  { rsMatrixNew<T> C(this->numRows, this->numCols); this->sub(this, &B, &C); return C; }

  /** Multiplies two matrices: C = A * B. */
  rsMatrixNew<T> operator*(const rsMatrixNew<T>& B) const
  { rsMatrixNew<T> C(this->numRows, B.numCols); this->mul(this, &B, &C); return C; }

  // todo: /, ==,,+=,-=,*=,-


protected:



  /** Updates the data-pointer inherited from rsMatrixView to point to the begin of our std::vector
  that holds the actual data. */
  void updateDataPointer()
  {
    if(data.size() > 0)
      this->d = &data[0];
    else
      this->d = nullptr;
  }


  /** \name Data */

  std::vector<T> data;
  // maybe we should just work with our inherited d pointer - but then we can't look easily at the
  // data in the debugger anymore -> very bad! so, nope!

};

/** Multiplies a scalar and a matrix. */
template<class T>
inline rsMatrixNew<T> operator*(const T& s, const rsMatrixNew<T>& A)
{
  rsMatrixNew<T> B(A);
  B.scale(s);
  return B;
}



  // some binary operators are defined outside the class such that the left hand operand does
  // not necesarrily need to be of class rsMatrix

  // matrix/vector functions:

  // todo: make some special-case classes for 2x2, 3x3 matrices which can use simpler algorithms
  // for some of the computations

#endif
