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
manipulating matrix elements.

*/

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

  inline void setAllValues(T value)
  {
    rsArray::fillWithValue(d, int(numRows * numCols), value);
  }

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
  static void sub(const rsMatrixView<T>& A, const rsMatrixView<T>& B, rsMatrixView<T>& C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, C), "arguments incompatible");
    for(int i = 0; i < A.numRows; i++)
      for(int j = 0; j < A.numCols; j++)
        C(i, j) = A.at(i, j) - B.at(i, j);
  }

  /** Computes the matrix product C = A*B. */
  static void mul(const rsMatrixView<T>& A, const rsMatrixView<T>& B, rsMatrixView<T>& C)
  {
    int N = A.numRows;
    int M = B.numRows;
    int P = A.numCols;
    rsAssert(P == B.numCols);
    rsAssert(N == C.numRows);
    rsAssert(P == C.numCols);
    // A: NxM, B: MxP, C: NxP 
    // verify these conditions - factor them out into areMultiplicable(A, B, C)

    for(int i = 0; i < N; i++) {
      for(int j = 0; j < P; j++) {
        C(i,j) = T(0);
        for(int k = 0; k < M; k++)
          C(i,j) += A.at(i,k) * B.at(k,j); }}
  }



  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Read and write access to matrix elements with row-index i and column-index j. */
  inline T& operator()(const int i, const int j)
  {
    return d[flatIndex(i, j)];

  }

  /** Read only accees - used mainly internally with const reference arguments (for example, 
  in add). */
  inline const T& at(const int i, const int j) const
  {
    return d[flatIndex(i, j)];
  }

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
  rsMatrixNew(int numRows = 0, int numColumns = 0);
  //rsMatrix(size_t numRows = 1, size_t numColumns = 1); // leads to memory leaks

 
  rsMatrixNew(int numRows, int numColumns, const std::vector<T>& data);


  /** Copy constructor. */
  //rsMatrix(const rsMatrix& other);
                              
  /** Move constructor. */
  //rsMatrix(const rsMatrix&& other);


  // todo: implement the various copy/move assigment operators and -constructors - this should
  // optimize returning values from functions and operators (avoid unnessary copying)



  /** Destructor. */
  ~rsMatrixNew() {}


  /** \name Setup */

  /** Sets the number of rows and columns, this matrix should have. ToDo: provide a way to retain 
  the data (optionally) - what does std::vector's resize do? Does it retain data...but if it does,
  it would be useless anyway in case the number of columns changed. */
  void setSize(int numRows, int numColumns);

    
  /** \name Manipulations */


  /** \name Inquiry */


  /** \name Decompositions */


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Compares matrices for equality */
  bool operator==(const rsMatrixNew<T>& B) const 
  { 
    if(numRows != B.numRows || numCols != B.numCols)
      return false;
    return rsArray::areBuffersEqual(d, B.d, getSize());
  }

  /** Compares matrices for inequality */
  bool operator!=(const rsMatrixNew<T>& B) const { return !(*this == B); }

  /** Adds two matrices: C = A + B. */
  rsMatrixNew<T> operator+(const rsMatrixNew<T>& B) const
  { rsMatrixNew<T> C(numRows, numCols); add(this, &B, &C); return C; }

  /** Subtracts two matrices: C = A - B. */
  //rsMatrixNew<T> operator-(const rsMatrixNew<T>& B) const
  //{ rsMatrixNew<T> C(numRows, numCols); sub(this, &B, &C); return C; }

  /** Multiplies two matrices: C = A * B. */
  //rsMatrixNew<T> operator*(const rsMatrixNew<T>& B) const
  //{ rsMatrixNew<T> C(numRows, B.numCols); mul(this, &B, &C); return C; }

  // todo: /, ==,,+=,-=,*= 


protected:

  /** \name Data */

  std::vector<T> data;
  
}; 

  // some binary operators are defined outside the class such that the left hand operand does
  // not necesarrily need to be of class rsMatrix

  // matrix/vector functions:

  // todo: make some special-case classes for 2x2, 3x3 matrices which can use simpler algorithms
  // for some of the computations

#endif
