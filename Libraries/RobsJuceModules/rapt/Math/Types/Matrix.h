#ifndef RAPT_MATRIX_H
#define RAPT_MATRIX_H

//=================================================================================================

/** A class for representing 2x2 matrices. They are treated as a special case because a lot of
things which are impractical in the general case can be done for the 2x2 case. For example, it's
possible to compute eigenvalues and eigenvectors via closed form formulas where in the general
case numerical algorithms are needed. */

template<class T>
class rsMatrix2x2
{

public:

  T a, b, c, d;
  // matrix coefficients [a b]
  //                     [c d]


  /** Stadard constructor. You can pass the matrix elements. If you pass nothing, an identity
  matrix will be created. */
  //rsMatrix2x2(T a = T(1), T b = T(0), T c = T(0), T d = T(1)) { setValues(a, b, c, d); }
  // ToDo: Maybe require arguments to be passed or initialize the matrix to the zero matrix
  // by default. The zero matrix seems to make more sense than the identity for default 
  // construction

  /** Constructor. Initializes elements with  given values. */
  rsMatrix2x2(T a, T b, T c, T d) { setValues(a, b, c, d); }

  /** Standard constructor. Leaves elements uninitialized.
  (...try to avoid using it - prefer RAII) */
  rsMatrix2x2() {}

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the elements of the matrix. */
  void setValues(T a, T b, T c, T d) { this->a = a; this->b = b; this->c = c; this->d = d; }

  /** Sets all elements of the matrix to zero. */
  void setZero() { this->a = this->b = this->c = this->d = T(0); }

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the determinant of this matrix. */
  T getDeterminant() const { return a*d - b*c; }

  /** Returns the trace (sum of diagonal elements) of this matrix.  */
  T getTrace() const { return a+d; }

  /** Returns the first eigenvalue of this matrix. */
  T getEigenvalue1() const { return rsLinearAlgebra::eigenvalue2x2_1(a, b, c, d); }

  /** Returns the second eigenvalue of this matrix. */
  T getEigenvalue2() const { return rsLinearAlgebra::eigenvalue2x2_2(a, b, c, d); }

  /** Returns the first eigenvector of this matrix. */
  rsVector2D<T> getEigenvector1() const
  { rsVector2D<T> v; rsLinearAlgebra::eigenvector2x2_1(a, b, c, d, &v.x, &v.y, true); return v; }

  /** Returns the second eigenvector of this matrix. */
  rsVector2D<T> getEigenvector2() const
  { rsVector2D<T> v; rsLinearAlgebra::eigenvector2x2_2(a, b, c, d, &v.x, &v.y, true); return v; }

  // ToDo: 
  // -Document, what counts as 1 and 2...i think, the smaller eigenvalue counts as eigenvalue 1, due 
  //  to the way, we do the +- business for the square-root in the solution formula for the 
  //  quadratic equation. That seems to be the most sensible convention.
  // -Document what happens when the eigenvalue of a real matrix happens to be complex. In this 
  //  case, the functions above must fail -> document in which way and/or maybe implement functions
  //  to compute complex eigenvalues and -vectors
  // -Implement a getEigenSystem() function that computes both eigenvalues and -vectors at once, 
  //  avoiding redundant calculations. Maybe also getEigenValues, getEigenVectors

  /** Returns the inverse of this matrix. */
  rsMatrix2x2<T> getInverse() const
  { T D = getDeterminant(); T s = T(1) / D; return rsMatrix2x2<T>(s*d, -s*b, -s*c, s*a); }

  /** Tests, if another matrix B is close to this matrix within a given tolerance (all components
  of the difference must be <= tolerance). */
  bool isCloseTo(const rsMatrix2x2<T>& B, T tol)
  {
    if(rsAbs(a-B.a) <= tol && rsAbs(b-B.b) <= tol && rsAbs(c-B.c) <= tol && rsAbs(d-B.d) <= tol)
      return true;
    return false;
  }
  // doesn't work with complex matrices

  /** Computes n-th power of the matrix using the closed form formula from:
  https://people.math.carleton.ca/~williams/papers/pdf/175.pdf , Eq. 2
  https://distill.pub/2017/momentum/ (section "Dynamics of momentum")
  ...it has been tested with real matrices with real eigenvalues (distinct and equal) - but what
  about real matrices with complex eigenvalues? The function won't work for them, i guess - try
  it!  */
  rsMatrix2x2<T> getPower(int n)
  {
    using Mat = rsMatrix2x2<T>;
    T ev1 = getEigenvalue1();
    T ev2 = getEigenvalue2();
    Mat I = identity();
    if(ev1 != ev2) {       // ToDo: maybe a tolerance is needed
      T ab  = T(1) / (ev1 - ev2);
      Mat X = ab * ((*this) - ev2*I);
      Mat Y = ab * (ev1*I - (*this));
      return pow(ev1, n) * X + pow(ev2, n) * Y; }
    else
      return pow(ev1, n-1) * (T(n) * *this - T(n-1) * ev1*I);
  }
  // ToDo:
  // Needs more tests, especially for the case when the eigenvalues are almost equal. I think,
  // we should use a relative tolerance, maybe based on ev1/ev2. If the quotient is too close 
  // to +1, use the 2nd branch. Then, test it for double and float.
  // What about non-integer powers? Does the formula also work for those? Try it and compare 
  // results to matrix power computed by Taylor expansion, maybe using sqrt as example, i.e. 
  // n = 0.5. What about n = -1? Would that be the inverse matrix?
  // What about the matrix exponential? And maybe sine and cosine, too? Maybe even logarithm? See
  // the code for geometric algebra in the Research repo for possible algorithms.


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Adds two matrices: C = A + B. */
  rsMatrix2x2<T> operator+(const rsMatrix2x2<T>& B) const
  { return rsMatrix2x2<T>(a + B.a, b + B.b, c + B.c, d + B.d); }

  /** Subtracts two matrices: C = A - B. */
  rsMatrix2x2<T> operator-(const rsMatrix2x2<T>& B) const
  { return rsMatrix2x2<T>(a - B.a, b - B.b, c - B.c, d - B.d); }

  /** Multiplies two matrices: C = A * B. */
  rsMatrix2x2<T> operator*(const rsMatrix2x2<T>& B) const
  { return rsMatrix2x2<T>(a*B.a + b*B.c, a*B.b + b*B.d, c*B.a + d*B.c, c*B.b + d*B.d); }

  /** Multiplies the left matrix operand with the inverse of the right matrix operand. */
  rsMatrix2x2<T> operator/(const rsMatrix2x2<T>& B) const { return *this * B.getInverse(); }

  /** Compares matrices for equality */
  bool operator==(const rsMatrix2x2<T>& B) const
  { return a == B.a && b == B.b && c == B.c && d == B.d; }

  /** Multiplies matrix by a vector: w = A*v */
  rsVector2D<T> operator*(const rsVector2D<T>& v) const
  { return rsVector2D<T>(a*v.x + b*v.y, c*v.x + d*v.y); }

  /** Divides matrix by a scalar divisor. */
  rsMatrix2x2<T> operator/(const T divisor) const
  { T s = T(1) / divisor; return rsMatrix2x2<T> (s*a, s*b, s*c, s*d); }

  // ToDo: 
  // -left multiplication w = v^H * A
  // -operators that take a scalar as left or right argument, +=, -=, *=, /=
  // -allow element access via A(i,j) operator to resemble interface of the general rsMatrix - that
  //  may be inefficient but sometimes convenient


  //-----------------------------------------------------------------------------------------------
  /** \name Factory */

  static rsMatrix2x2<T> zero()     { return rsMatrix2x2<T>(T(0), T(0), T(0), T(0)); }

  static rsMatrix2x2<T> identity() { return rsMatrix2x2<T>(T(1), T(0), T(0), T(1)); }

  /** Returns the commutator C of the two matrices A and B: C = A*B - B*A. In general, matrix
  multiplication is non-commutative, but for some special cases, it may be commutative nonetheless.
  The commutator captures, how non-commutative two matrices behave when being multiplied. If the
  two matrices commute (i.e. behave commutatively), their commutator is the zero matrix. */
  static rsMatrix2x2<T> commutator(const rsMatrix2x2<T>& A, const rsMatrix2x2<T>& B)
  { return A*B - B*A; }
  // see: https://en.wikipedia.org/wiki/Commutator#Ring_theory
  // -Maybe implement also the anticommutator defined there as: {A,B} = A*B + B*A
  // -Maybe implement it as a free function rsCommutator for arbitrary types. The notion of a 
  //  commutator may make sense for other things as well (for example, multivectors)

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Solves A*x = b for x. The matrix A is assumed to be regular. If it's singular, we'll 
  encounter a division by zero and produce garbage. */
  static void solve(const rsMatrix2x2<T>& A, rsVector2D<T>& x, const rsVector2D<T>& b)
  {
    //T tol = 1000 * RS_EPS(T);
    //rsAssert(rsAbs(A.getDeterminant()) > tol, "Handling of singular matrices not implemented");
    // maybe we should divide the determinant by the maximum of the matrix elements or something
    // to get a relative measure - we hit this, when all elements are very small, wich should not 
    // be considered to be a problem

    T D = A.getDeterminant();
    T s = T(1) / D;
    x.x = s * (A.d * b.x - A.b * b.y);
    x.y = s * (A.a * b.y - A.c * b.x);

    // Old implementation for reference:
    //rsMatrix2x2<T> Ai = A.getInverse();
    //x.x = Ai.a * b.x + Ai.b * b.y;
    //x.y = Ai.c * b.x + Ai.d * b.y;
  }
  // ToDo: 
  // -Handle singluar matrices: in the overdetermined case, produce a least squares
  //  approximation, in the underdetermined case, produce a minimum norm solution, maybe that 
  //  should be in a function solveSingular
  // -Maybe implement it without using rsVector2D - just take references to coordinates instead
  // -Maybe turn into free function rsSolve, which in general may have an abstract signature
  //  TRhs rsSolve(const TOp& op, const TRhs& rhs) using templates for the operator type and 
  //  right-hand-side type and returning an object of same type as rhs or mabye use a non-const ref
  //  for the result.

  /** Like solve, but checks for division by zero and assigns zero to the result in this case. */
  static void solveSave(const rsMatrix2x2<T>& A, rsVector2D<T>& x, const rsVector2D<T>& b)
  {
    T D = A.getDeterminant();
    if(D == T(0)) { x.x = x.y = T(0); return; }
    T s = T(1) / D;
    x.x = s * (A.d * b.x - A.b * b.y);
    x.y = s * (A.a * b.y - A.c * b.x);
  }
  // ToDo: 
  // -Maybe use an absolute threshold that defaults to zero, like: if(abs(D) <= thresh) { ... } 
  // -Maybe don't return zero but instead a least-squares or minimum-norm solution. But maybe both
  //  variants should be available under different names. Maybe call the function above solveOrZero
  //  and the variant using least-squares/min-norm solveOrApprox or just solveApprox or 
  //  solveLeastSquares or solveBest

};

/** Multiplies a scalar and a matrix. */
template<class T>
inline rsMatrix2x2<T> operator*(const T& s, const rsMatrix2x2<T>& A)
{
  return rsMatrix2x2<T>(s*A.a, s*A.b, s*A.c, s*A.d);
}

// ToDo: 
// -Implement matrix exponential for 2x2 matrices as free function rsExp. If A is invertible, use
//  the formula based on diagonalization. If it isn't...dunno...maybe use Taylor expansion similar
//  to the geometric algebra code in the research codebase? Or maybe Jordan normal form? See also
//  https://www.youtube.com/watch?v=Iz7PSlTpjyI -> exp(tr(A)) = det(exp(A)) -> verify 
//  experimentally. also: tr(A) = sum of eigenvalues, det(A) = product of eigenvalues

//=================================================================================================

/** Under construction...it works but is still lacking functionality...

Class for representing 3x3 matrices...  */

template<class T>
class rsMatrix3x3
{

public:

  rsMatrix3x3() {}

  rsMatrix3x3(T a, T b, T c, T d, T e, T f, T g, T h, T i)
  { setValues(a, b, c, d, e, f, g, h, i); }


  void setValues(T a, T b, T c, T d, T e, T f, T g, T h, T i) 
  { 
    A[0] = a; A[1] = b; A[2] = c;
    A[3] = d; A[4] = e; A[5] = f;
    A[6] = g; A[7] = h; A[8] = i;
  }


  /** Element access for read and write. */
  T& operator()(const int i, const int j) { return A[flatIndex(i, j)]; }

  /** Element access for read only. */
  const T& operator()(const int i, const int j) const { return A[flatIndex(i, j)]; }

  // ToDo: arithmetic operators, 

protected:

  inline int flatIndex(const int i, const int j) const
  {
    rsAssert(i >= 0 && i < 3, "invalid row index");
    rsAssert(j >= 0 && j < 3, "invalid column index");
    return 3*i + j;
  }

  T A[9];  // array of coeffs, capitalized  because we want to use "a" for the 1st element

};


// ToDo: 
// -rsMatrix2x3 (useful for 3D -> 2D projections in graphics)
// -rsMatrix3x2 (maybe useful in differential geometry of 2D surfaces embedded in 3D)
// -maybe the small matrices should go into a dedicated SmallMatrices.h file. Maybe the same
//  should be done for the small vectors in Vector.h -> move 2D/3D to SmallVectors.h - wait:
//  there are only small vectors in Vector.h, so maybe just rename it to SmallVectors.h
// -These "small" versions should go up to 4 - 4D is sometimes needed when using homogeneous
//  coordinates in graphics (and also in relativity)


//=================================================================================================

/** This is a class for treating raw C-arrays as matrices. It does not store/own the actual matrix
data. It just acts as wrapper around an existing array for conveniently accessing and manipulating
matrix elements via row/column indicies using the () operator with two integers. The matrix data is 
treated as being stored in row-major format. 

Note: LaPack assumes column-major format, so the storage is not compatible. To mix rsMatrix with
LaPack routines, you will need to do transpositions in memory. I think, row-major is better for
computing matrix-vector products, access pattern wise and its also sometimes convenient to be able
to address a row of a matrix as vector/array. Column-major storage may be more convenient and 
efficient for certain other operations that i found to be less common in my typical audio related 
applications, so i opted for row-major. ....your mileage may vary. */

// ToDo:
// -create benchmarks for and optimize elementary row- and column operations and matrix-vector
//  multiplication (these are the most important linear algebra operations, so they should be fast)
//  -> implement benchmarks
// -maybe make the storage format more flexible by using strides (allows for column major storage
//  and also block submatrices without copying)...but measure the performance impact - if it's too 
//  costly, then maybe don't

template<class T>
class rsMatrixView
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Default constructor. */
  rsMatrixView() {}

  /** Creates a matrix view with the given number of rows and columns for the given raw array of
  values in "data". The view will *not* take ownership over the data. */
  rsMatrixView(int numRows, int numColumns, T* data)
  {
    rsAssert(numRows >= 1 && numColumns >= 1 && data != nullptr);
    this->numRows = numRows;
    this->numCols = numColumns;
    dataPointer = data;
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Re-interprets the arrangement of the underlying array as having the new given numbers of rows
  and columns. Their product must remain the same, though. For example, you can reshape a 3x4
  matrix into 4x3, 2x6, 6x2, 1x12, 12x1 but nothing else. ..todo: maybe lift that restriction? */
  void setShape(int newNumRows, int newNumColumns)
  {
    rsAssert(newNumRows*newNumColumns == numRows*numCols);
    numRows = newNumRows;
    numCols = newNumColumns;
  }
  // maybe rename or make an alias "reshape" to be consistent with NumPy, MatLab, etc.

  /** Resets the number of rows and columns to zero and the dataPointer to nullptr. Should be 
  called whenever you need to invalidate our pointer member. */
  void reset()
  {
    numRows = 0;
    numCols = 0;
    dataPointer = nullptr;
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */
  // maybe split out a Retrieval section (copyRow, getRowPointer, getDataPointer, etc.)

  /** Returns true, iff matrices A and B have the same number of rows and columns, such that
  their sum A+B and difference A-B can be formed. */
  static bool areSameShape(const rsMatrixView<T>& A, const rsMatrixView<T>& B)
  { return A.numRows == B.numRows && A.numCols == B.numCols; }

  /** Returns true, iff the matrices A and B have dimensions, such that the matrix-matrix product
  A*B can be formed. */
  static bool areMultiplicable(const rsMatrixView<T>& A, const rsMatrixView<T>& B)
  { return A.numCols == B.numRows; }

  /** Returns true, iff the rhs matrix is equal to this matrix with an optional tolerance. */
  bool equals(const rsMatrixView<T>& rhs, T tolerance = T(0)) const
  {
    return areSameShape(*this, rhs)
      && rsArrayTools::almostEqual(dataPointer, rhs.dataPointer, getSize(), tolerance);
  }


  /** Returns true, iff this matrix has the given shape. */
  bool hasShape(int numRows, int numCols) const
  { return this->numRows == numRows && this->numCols == numCols; }

  /** Returns true, iff B has the same shape as this matrix. */
  bool hasSameShapeAs(const rsMatrixView<T>& B) const
  { return areSameShape(*this, B); }

  /** Returns the number of rows. */
  int getNumRows() const { return numRows; }

  /** Returns the number of columns. */
  int getNumColumns() const { return numCols; }

  /** Returns the total size, i.e. the number of elements in the matrix. */
  int getSize() const { return numRows * numCols; }

  /** Returns true, iff this matrix is a row vector. */
  bool isRowVector() const { return numRows == 1; }

  /** Returns true, iff this matrix is a column vector. */
  bool isColumnVector() const { return numCols == 1; }

  /** Returns true, iff this matrix is either a row-vector or column-vector. */
  bool isVector() const { return isRowVector() || isColumnVector(); }

  /** Returns true, iff this matrix is a square matrix, i.e. has the same number of rows and
  columns. */
  bool isSquare() const { return numRows == numCols; }

  /** Returns true, iff given index i is a valid row index. */
  bool isValidRowIndex(int i) const { return i >= 0 && i < numRows; }

  /** Returns true, iff given index i is a valid column index. */
  bool isValidColumnIndex(int j) const { return j >= 0 && j < numCols; }

  /** Returns true, iff all elements are zero. */
  bool isZero() const { return rsArrayTools::isAllZeros(dataPointer, getSize()); }
  // todo: allow for a tolerance ..or maybe have a function isAlmostZero for that

  bool isZero(T tol) const
  { return rsArrayTools::isAllZeros(dataPointer, getSize(), tol); }

  bool isRowZero(int rowIndex, T tol = T(0)) const
  {
    for(int j = 0; j < getNumColumns(); ++j)
      //if( rsAbs(at(rowIndex, j)) > tol )
      if( rsGreaterAbs(at(rowIndex, j), tol) )
        return false;
    return true;
  }

  bool areRowsZero(int startRow, int endRow, T tol = T(0)) const
  {
    for(int i = startRow; i <= endRow; ++i)
      if(!isRowZero(i, tol))
        return false;
    return true;
  }

  bool isColumnZero(int columnIndex, T tol = T(0)) const
  {
    for(int i = 0; i < getNumRows(); ++i)
      if( rsAbs(at(i, columnIndex)) > tol )
        return false;
    return true;
  }


  /** Returns true, iff one of the columns of the matrix consists of all zeros. */
  bool containsZeroColumn(T tol = T(0)) const
  {
    for(int j = 0; j < numCols; j++)
      if( isColumnZero(j, tol) )
        return true;
    return false;
  }
  // needs test, todo: document for what this is good for...i think, it's a condition for which 
  // certain linear algebra operations don't work

  /** Returns true, iff this matrix (denoted as A) is symmetric (up to a given tolerance), i.e. 
  A(i,j) == A(j,i) for all i,j. Symmetry considerations usually apply only to square matrices. If 
  A isn't a square matrix, it will be considered non-symmetric, regardless of its content. */
  bool isSymmetric(T tol = T(0)) const
  {
    if(numRows != numCols) return false;  // non-square matrices are never considered symmetric
    for(int i = 1; i < numRows; i++) {
      for(int j = i; j < numCols; j++) {  
        T d = rsAbs(at(i,j) - at(j,i));
        if(d > tol)
          return false; }}
    return true;
  }
  // needs test, todo: implement test for antisymmetry - the structure is the same, just that we 
  // need to use at(i,j) + at(j,i) instead of at(i,j) - at(j,i)
  // Shouldn't the inner j-loop start at i+1? A(i,j) is always equal to A(j,i) when i==j, so it 
  // seems, we are checking one value too many. It doesn't make the result wrong - we are just 
  // doing (slightly) more work than strictly necessary. But that will not translate to the check
  // for antisymmetry because A(i,i) is only equal to -A(i,i) when it's zero.


  /** Returns a const pointer to the data for read access as a flat array. */
  const T* getDataPointerConst() const { return dataPointer; }

  /** Returns a pointer to the data for read and write access as a flat array. */
  T* getDataPointer() { return dataPointer; }

  /** Returns a pointer to the stored data. When using this, be sure that you know exactly what
  you are doing.... */
  //T* getData() { return dataPointer; }

  /** Returns a pointer to a given row. */
  T* getRowPointer(int rowIndex)  { return &dataPointer[rowIndex*numCols]; }
  // if we later support column-major storage, we should assert that the matrix in row-major
  // storage

  /** Returns a const pointer to a given row. */
  const T* getRowPointerConst(int rowIndex) const { return &dataPointer[rowIndex*numCols]; }

  /** Return the number of data entries that overlap in matrix-views A and B. For many
  computations, it is required, that A and B have non-overlapping data, so this function can be
  used to check, if those operations can be done on two particular matrices. */
  /*
  static size_t getDataOverlap(const rsMatrixView<T>& A, const rsMatrixView<T>& B)
  {
    if(B.dataPointer < A.dataPointer)  // A's data must start before or at the same memory
      return getDataOverlap(B, A);     // location as B's data - otherwise use swapped arguments

    size_t safeStartForB = A.dataPointer + A.getSize(); // 1 position after A's last slot
    if(B.dataPointer >= safeStartForB)
      return 0;  // no overlap
    else
      return rsMin(safeStartForB - B.dataPointer, B.getSize());
  }
  */
  // actually, we should move this to rsArrayTools::getOverlap(T* x, size_t Nx, T*y, size_t Ny)
  // needs unit test

  /** Returns the minimum value of all elements in the matrix. */
  T getMinimum() const { return rsArrayTools::minValue(dataPointer, getSize()); }

  /** Returns the maximum value of all elements in the matrix. */
  T getMaximum() const { return rsArrayTools::maxValue(dataPointer, getSize()); }

  /** Returns the maximum absolute value of all elements in the matrix. */
  T getAbsoluteMaximum() const { return rsArrayTools::maxAbs(dataPointer, getSize()); }

  /** Returns the number of nonzero elements in this matrix. */
  int getNumNonZeros() const { return rsArrayTools::numNonZeros(dataPointer, getSize()); }

  /** Returns the density of this matrix which is a number between 0 and 1 defined as the number 
  of nonzero elements divided by the total number of elements. Sparse matrices have a low density
  (near zero), fully populated matrices have a density of 1. */
  T getDensity() const { return T(getNumNonZeros()) / T(getSize()); }

  /** Computes the trace of the matrix which is the sum of the diagonal elements. */
  T getTrace() const
  {
    T t = T(0);
    for(int i = 0; i < rsMin(this->numRows, this->numCols); ++i)
      t = t + this->at(i, i);   // todo: use +=
    return t;
  }

  /** Computes the product of the diagonal elements (does this also have a special name?) */
  T getDiagonalProduct() const
  {
    T p = T(1);
    for(int i = 0; i < rsMin(this->numRows, this->numCols); ++i)
      p = p * this->at(i, i);
    return p;
  }
  // todo: use *= - (needs implementation of that operator in rsRationalFunction)

  /** Returns a copy of the i-th row as a std::vector. */
  std::vector<T> getRowAsVector(int i) const { return toVector(getRowPointerConst(i), numCols); }

  /** Copies data of i-th row into given array arr which should be of length numCols. */
  void copyRow(int i, T* arr) const
  {
    rsArrayTools::copy(getRowPointerConst(i), arr, numCols);
  }
  // needs test

  /** Copies data of j-th column into given array arr which should be of length numRows. */
  void copyColumn(int j, T* arr) const
  {
    for(int i = 0; i < numRows; i++)
      arr[i] = at(i, j);
  }
  // needs test

  //-----------------------------------------------------------------------------------------------
  /** \name Manipulation */

  /** Sets all matrix elements to zero. */
  void setToZero() { rsArrayTools::fillWithZeros(dataPointer, getSize()); }

  /** Sets the matrix elements to the identity matrix, i.e. fills the main diagonal with ones and
  the rest with zeros. If the matrix is not square, then the overhanging portion to the right or
  bottom will be all zeros as well (-> verify this). */
  void setToIdentity() { setToZero(); setDiagonalValues(T(1)); }
  // needs test

  /** Sets all elements in the matrix to the given value. */
  void setAllValues(T value) { rsArrayTools::fillWithValue(dataPointer, getSize(), value); }
  // maybe rename to fill, or make an alias

  /** Initializes all elements with given value. */
  //void init(T value = T(0)) { RAPT::rsArrayTools::fillWithValue(dataPointer, getSize(), value); }
  // maybe remove - is redundant with setAllValues

  /** Sets all elements on the main diagonal to the given value. If the matrix is not square, only
  the top- or left square submatrix will be affected. */
  void setDiagonalValues(T value)
  {
    for(int i = 0; i < rsMin(this->numRows, this->numCols); i++)
      dataPointer[i*numCols + i] = value;
  }
  // needs test

  /** Sets the elements on the main diagonal to the values in the given array which must be of
  length min(numRows, numColumns). */
  void setDiagonalValues(T* values)
  {
    for(int i = 0; i < rsMin(this->numRows, this->numCols); i++)
      dataPointer[i*numCols + i] = values[i];
  }
  // needs test

  void setRow(int i, const T* values) 
  { 
    rsAssert(i >= 0 && i < numRows);
    rsArrayTools::copy(values, getRowPointer(i), numCols); 
  }
  // needs test

  void setRow(int i, const std::vector<T>& values)
  {
    rsAssert((int)values.size() == numCols);
    setRow(i, &values[0]);
  }
  // convenience function

  void setColumn(int j, const T* values)
  {
    rsAssert(j >= 0 && j < numCols);
    for(int i = 0; i < numRows; i++)
      (*this)(i, j) = values[i];
  }
  // needs test..maybe rsArrayTools::copy should have a version with srcStride, dstStride
  // make convenience function fo std::vector like we have for setRow





  /** Scales all elements in the matrix by a given factor. */
  void scale(T factor) { rsArrayTools::scale(dataPointer, getSize(), factor); }

  /** Negates all values of the matrix, i.e. inverts their sign. */
  void negate() { rsArrayTools::negate(dataPointer, dataPointer, getSize()); }

  // todo: conjugate


  /** Scales the row with given index by the given scale factor. */
  void scaleRow(int rowIndex, T scaler)
  {
    rsAssert(isValidRowIndex(rowIndex), "row index out of range");
    for(int j = 0; j < numCols; ++j)
      (*this)(rowIndex, j) *= scaler;
  }
  // maybe use rsArrayTools::scale

  /** Swaps the two rows with given row indices i1 and i2. */
  void swapRows(int i1, int i2)
  {
    rsAssert(isValidRowIndex(i1) && isValidRowIndex(i1), "row index out of range");
    for(int j = 0; j < numCols; ++j)
      rsSwap((*this)(i1, j), (*this)(i2, j));
  }
  // may be optimized by using fixed base-pointers to each row and loop increment 1 - no
  // recompuation of the row-start in the iterations - the (i,j) operator always does a
  // multiplication

  /** Adds a multiple of the row with index iSrc to the row with index iDst. The multiplier is
  given by weight. */
  void addWeightedRowToOther(int iSrc, int iDst, T weight)
  {
    rsAssert(isValidRowIndex(iSrc) && isValidRowIndex(iDst), "row index out of range");
    for(int j = 0; j < numCols; ++j)
      (*this)(iDst, j) += weight * (*this)(iSrc, j);
  }
  // maybe rename to addScaledRowToOther - 2 chars shorter
  // optimize using base-pointers

  /** Like above but does it only for the column-indices in between minCol and maxCol (both
  ends inclusive). This is mainly for optimizing row operations when it is known that beyond
  minCol...maxCol only zeros occur in the source-row (this occurs in Gaussian elimination). */
  void addWeightedRowToOther(int iSrc, int iDst, T weight, int minCol, int maxCol)
  {
    rsAssert(isValidRowIndex(iSrc) && isValidRowIndex(iDst), "row index out of range");
    rsAssert(isValidColumnIndex(minCol) && isValidColumnIndex(maxCol), "column index out of range");

    //rsAssert(minCol >= 0 && maxCol < numCols, "column index out of range");
    // minCol > numCols is allowed - in this case, the loop is not entered (Gaussian elimination
    // may produce such values)

    for(int j = minCol; j <= maxCol; ++j)
      (*this)(iDst, j) += weight * (*this)(iSrc, j);
  }
  // optimize using base-pointers


  /** Scales the column with given index by the given scale factor. */
  void scaleColumn(int columnIndex, T scaler)
  {
    rsAssert(isValidColumnIndex(columnIndex), "column index out of range");
    for(int i = 0; i < numRows; ++i)
      (*this)(i, columnIndex) *= scaler;
  }
  // needs test, maybe use rsArrayTools::scale (needs to be generalized to accept an optional
  // stride parameter that defaults to 1 - but generalizing like that may slow it down for the
  // common case -> benchmark)

  /** Swaps the two columns with given indices j1 and j2. */
  void swapColumns(int j1, int j2)
  {
    rsAssert(isValidColumnIndex(j1) && isValidColumnIndex(j1), "column index out of range");
    for(int i = 0; i < numRows; ++i)
      rsSwap((*this)(i, j1), (*this)(i, j2));
  }

  /** Adds a multiple of the column with index jSrc to the column with index jDst. The multiplier
  is given by weight. */
  void addWeightedColumnToOther(int jSrc, int jDst, T weight)
  {
    rsAssert(isValidColumnIndex(jSrc) && isValidColumnIndex(jDst), "column index out of range");
    for(int i = 0; i < numRows; ++i)
      (*this)(i, jDst) += weight * (*this)(i, jSrc);
  }
  // can also be optimized with something like:
  // T* s = flatIndex(0, jSrc);
  // T* d = flatIndex(0, jDst);
  // int stride = numCols
  // for(int i = 0; i < numRows; ++i) {
  //   *d += weight * *s; s += stride; d += stride; }

  // optimize at least the elementary row-operations because these occur frequently in matrix
  // algorithms like Gaussian elimination


  //-----------------------------------------------------------------------------------------------
  /** \name Arithmetic */

  /** Adds elements of A to corresponding elements in B and stores results in C. */
  static void add(const rsMatrixView<T>& A, const rsMatrixView<T>& B, rsMatrixView<T>* C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, *C), "arguments incompatible");
    rsArrayTools::add(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }

  /** Subtracts elements of B from corresponding elements A in and stores results in C. */
  static void subtract(const rsMatrixView<T>& A, const rsMatrixView<T>& B, rsMatrixView<T>* C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, *C), "arguments incompatible");
    rsArrayTools::subtract(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }

  /** Weighted sum of the mA x nA matrix A and the mB x nB matrix B. The result matrix C needs 
  to have a shape max(mA, mB) x max(nA, nB). For rows or columns that are present in one matrix 
  but not in the other, the other matrix is zero padded appropriately. */
  static void weightedSum(const rsMatrixView<T>& A, T wA, const rsMatrixView<T>& B, T wB,
    rsMatrixView<T>& C)
  {
    int M = rsMax(A.getNumRows(),    B.getNumRows());
    int N = rsMax(A.getNumColumns(), B.getNumColumns());
    rsAssert(C.hasShape(M, N));
    for(int m = 0; m < M; m++)
      for(int n = 0; n < N; n++)
        C(m, n) = wA * A.getElementPadded(m, n) + wB * B.getElementPadded(m, n);
  }
  // needs test

  /** Multiplies the two matrices element-wise. */
  static void elementwiseMultiply(
    const rsMatrixView<T>& A, const rsMatrixView<T>& B, rsMatrixView<T>* C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, *C), "arguments incompatible");
    rsArrayTools::multiply(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }

  /** Divides the two matrices element-wise. */
  static void elementwiseDivide(
    const rsMatrixView<T>& A, const rsMatrixView<T>& B, rsMatrixView<T>* C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, *C), "arguments incompatible");
    rsArrayTools::divide(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }

  /** Computes the matrix product C = A * B. */
  static void matrixMultiply(
    const rsMatrixView<T>& A, const rsMatrixView<T>& B, rsMatrixView<T>* C)
  {
    rsAssert(A.numCols  == B.numRows);
    rsAssert(C->numCols == B.numCols);
    rsAssert(C->numRows == A.numRows);
    for(int i = 0; i < C->numRows; i++) {
      for(int j = 0; j < C->numCols; j++) {
        (*C)(i, j) = T(0);
        for(int k = 0; k < A.numCols; k++)
          (*C)(i, j) += A.at(i, k) * B.at(k, j); }}
  }
  // -maybe implement matrixDivide: A/B := A * inverse(B)


  /** Performs the multiply-and-accumulate operation: C = C + A * B. */
  static void matrixMultiplyAccumulate(
    const rsMatrixView<T>& A, const rsMatrixView<T>& B, rsMatrixView<T>* C)
  {
    rsAssert(A.numCols  == B.numRows);
    rsAssert(C->numCols == B.numCols);
    rsAssert(C->numRows == A.numRows);
    for(int i = 0; i < C->numRows; i++) {
      for(int j = 0; j < C->numCols; j++) {
        for(int k = 0; k < A.numCols; k++)
          (*C)(i, j) += A.at(i, k) * B.at(k, j); }}
  }







  /** Computes the matrix product C = A^T * B. */
  static void matrixMultiplyFirstTransposed(
    const rsMatrixView<T>& A, const rsMatrixView<T>& B, rsMatrixView<T>* C)
  {
    // get rid, add assertions
    int N = A.numRows;
    int M = A.numCols;
    int P = B.numCols;
    for(int i = 0; i < M; i++) {
      for(int j = 0; j < P; j++) {
        (*C)(i, j) = T(0);
        for(int k = 0; k < N; k++)
          (*C)(i, j) += A.at(k, i) * B.at(k, j); }}
  }
  // needs test
  // todo: matrixMultiplySecondTransposed

  // implement functions that compute A^T * A and A * A^T, a "sandwich" X^T A X, X A X^T
  // but maybe more generally A^T B where A and B may or may not be the same
  // matrixMultiplyFirstTransposed, matrixMultiplySecondTransposed matrixmultiplyBothTransposed,
  // adapt code from rsMatrixTools


  /** Fills the matrix B with the transpose of matrix A. Assumes that A and B have compatible
  shapes. Matrix A and B must point to non-overlapping arrays */
  static void transpose(const rsMatrixView<T>& A, rsMatrixView<T>* B)
  {
    //rsAssert(A.dataPointer != B->dataPointer, "can't be used in place");
    //rsAssert(getDataOverlap(A, *B) == 0, "can't be used with overlapping matrices");
    rsAssert(A.numRows == B->numCols);
    rsAssert(A.numCols == B->numRows);
    for(int i = 0; i < A.numRows; i++)
      for(int j = 0; j < A.numCols; j++)
        (*B)(j,i) = A.at(i,j);
  }
  // B is a pointer and not a reference because it's an output - adopt that idiom generally:
  // output variables are always passed as pointers, never as references to make it visible in 
  // client code when a function parameter is an output

  /** Transposes the square matrix A in place. */
  static void transposeSquare(rsMatrixView<T>* A)
  {
    rsAssert(A->isSquare());
    for(int i = 0; i < A->getNumRows(); i++)
      for(int j = i+1; j < A->getNumRows(); j++)
        rsSwap((*A)(i,j), (*A)(j,i));
  }

  /** Computes the Kronecker product between matrices A and B and stores the result in C. Assumes,
  that C has the right dimensions. For more info, see the documentation of 
  rsMatrix::getKroneckerProduct. */
  static void kroneckerProduct(
    const rsMatrixView<T>& A, const rsMatrixView<T>& B, rsMatrixView<T>* C)
  {
    rsAssert(C->numRows == A.numRows * B.numRows);
    rsAssert(C->numCols == A.numCols * B.numCols);
    for(int ia = 0; ia < A.numRows; ia++) {
      for(int ja = 0; ja < A.numCols; ja++) {
        int startRow = ia*B.numRows;
        int startCol = ja*B.numCols;
        for(int ib = 0; ib < B.numRows; ib++) {
          for(int jb = 0; jb < B.numCols; jb++) {
            (*C)(startRow+ib, startCol+jb) = A.at(ia,ja) * B.at(ib, jb); }}}}
  }

  /** Computes the matrix-vector product y = A*x of this matrix with the given vector x and stores
  the result in y where x and y are given as raw arrays. The length of x must match the number of
  columns and the length of y must match the number of rows. */
  void product(const T* x, T* y) const
  {
    rsAssert(x != y, "Can't be used in place");
    for(int i = 0; i < getNumRows(); i++) {
      y[i] = T(0);
      for(int j = 0; j < getNumColumns(); j++)
        y[i] += at(i, j) * x[j]; }
  }

  /** Computes the matrix-vector product y = A^T * x where ^T denotes taking the matrix transpose.
  @see product */
  void transProduct(const T* x, T* y) const
  {
    rsAssert(x != y, "Can't be used in place");
    for(int i = 0; i < getNumColumns(); i++) {
      y[i] = T(0);
      for(int j = 0; j < getNumRows(); j++)
        y[i] += at(j, i) * x[j]; }
  }
  // todo: maybe rename to transposeProduct or productTransposed and/or implement functions that 
  // use the  Hermitian transpose maybe productHermitian


  /** Convenience function to compute matrix-vector product y = A*x, taking a raw array for x as
  input and producing the result as a std::vector. */
  std::vector<T> productWith(const T* x) const
  { std::vector<T> y(getNumRows()); product(x, &y[0]); return y; }

  /** 2D convolution of the two matrices A and B into the result matrix C. The number of rows and 
  columns of C must equal the sum of the respective numbers of A and B minus 1. */
  static void convolve(
    const rsMatrixView<T>& A, const rsMatrixView<T>& B, rsMatrixView<T>* C)
  {
    int Ma = A.numRows; int Na = A.numCols;
    int Mb = B.numRows; int Nb = B.numCols;
    int Mc = Ma+Mb-1;   int Nc = Na+Nb-1;
    rsAssert(C->numRows == Mc && C->numCols == Nc, "Result matrix has wrong shape");
    for(int m = 0; m < Mc; m++) {
      for(int n = 0; n < Nc; n++) {
        T s = T(0);
        for(int i = rsMax(0, m-Ma+1); i <= rsMin(Mb-1, m); i++) {
          for(int j = rsMax(0, n-Na+1); j <= rsMin(Nb-1, n); j++) {
            s += B(i, j) * A(m-i, n-j);  }}
        (*C)(m, n) = s; }}
  }
  // ToDo: 
  // -make it possible to work in place -> reverse directions of outer loops, i.e. run m from
  //  Mc-1 down to 0 and likewise for n
  // -implement this algo - it is nice: it reduces the 2D convolution problem to a 1D 
  //  convolution: Polynomial multiplication and FFT
  //  https://cseweb.ucsd.edu/~slovett/teaching/SP15-CSE190/poly-mult-and-FFT.pdf

  //-----------------------------------------------------------------------------------------------
  /** \name Accessors */

  /** Read and write access to matrix elements with row-index i and column-index j. */
  T& operator()(const int i, const int j) { return dataPointer[flatIndex(i, j)]; }

  /** Read only access to matrix elements with row-index i and column-index j. */
  const T& operator()(const int i, const int j) const { return dataPointer[flatIndex(i, j)]; }

  /** Read only access - used mainly internally with const reference arguments (for example,
  in add). */
  const T& at(const int i, const int j) const { return dataPointer[flatIndex(i, j)]; }
  // maybe rename to get - do we actually need this? - if not, get rid! it's sometimes convenient
  // when we have a pointer to a matrix. The syntax A->at(i,j) clearer nicer than (*A)(i,j)

  // void set(i, j, val) ...to make it compatible with old implementation

  /** Returns the element at index pair i,j or padding, if i or j is out of range. Can be used to 
  conveniently access elements of a (zero-)padded matrix. */
  T getElementPadded(const int i, const int j, const T padding = T(0)) const
  {
    if(!isValidRowIndex(i) || !isValidColumnIndex(j))
      return padding;
    return dataPointer[flatIndex(i, j)];
  }

  /** Converts a row index i and a column index j to a flat array index. */
  int flatIndex(const int i, const int j) const
  {
    rsAssert(i >= 0 && i < numRows, "invalid row index");
    rsAssert(j >= 0 && j < numCols, "invalid column index");
    return numCols*i + j;
    // todo:
    //  -be more general: colStride*i + rowStride*j. goal: allow row-major and column-major storage
    //   while the syntax of the () operator is always row-major (as is conventional in math)
    //   regardless whatever the internal storage format is - column major storage is required for
    //   compatibility with lapack
    // -maybe be even more general: colOffset + colStride*i + (rowOffset + rowStride)*j
    //  -> may allow to access sub-matrices with the same syntax (todo: verify formula)
    //  ...but maybe that should be done in a class rsSubMatrixView
  }


protected:

  /** \name Data */

  int numRows = 0, numCols = 0;  // number of rows and columns
  T *dataPointer = nullptr;      // pointer to the actual data

};

//=================================================================================================

/** This is a class for representing matrices and doing mathematical operations with them. It's
implemented as subclass of rsMatrixView and by default stores the actual matrix data in a 
std::vector. Copy- and move constructors and -assignment operators have been implemented in order 
to avoid unnecessary heap allocations in arithmetic expressions with matrices (return value copy
elision). 

The second template parameter allows client code to select the underlying storage container for the
actual matrix data. It defaults to std::vector and if you want to use anything else, your 
replacement needs to partially conform to the interface and implementation of std::vector. It needs
to store the data in a contiguous memory area, it needs to have a resize() method...tbc... 

ToDo: 
-those member functions that use std::vector for parameters or return types should be changed 
 to use the template parameter V instead (the transition from hardcoded std::vector and V is not
 yet complete)

*/

template<class T, class V = std::vector<T>>
class rsMatrix : public rsMatrixView<T>
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Assignment/Destruction */

  /** Standard constructor. You must pass the initial number of rows and columns */
  rsMatrix(int numRows = 0, int numColumns = 0)
  {
    setShape(numRows, numColumns);
    // todo: optionally init with zeros
  }

  /** Constructor to create a matrix from a flat raw array. */
  rsMatrix(int numRows, int numColumns, const T* data)
  {
    setShape(numRows, numColumns);
    rsArrayTools::copy(data, this->getDataPointer(), this->getSize());
  }

  /** Constructor to create a matrix from an array-of-arrays - mostly for conveniently converting
  matrices in the old representation into the new one. */
  rsMatrix(int numRows, int numColumns, T** data)
  {
    setShape(numRows, numColumns);
    for(int i = 0; i < numRows; i++)
      for(int j = 0; j < numColumns; j++)
        (*this)(i,j) = data[i][j];
  }

  /** Destructor. */
  ~rsMatrix()
  {
    //int dummy = 0; // to figure out, when it gets called for debugging
  }

  /** Creates matrix from a std::vector.  */
  rsMatrix(int numRows, int numColumns, const std::vector<T>& newData) : data(newData)
  {
    rsAssert(numRows*numColumns == newData.size());
    this->numRows = numRows;
    this->numCols = numColumns;
    updateDataPointer();
  }

  /** Creates matrix from an unnamed/temporary/rvalue std::vector - convenient to initialize
  elements. You can initialize matrices like this:
    rsMatrix<double> A(2, 3, {1.,2.,3., 4.,5.,6.});   */
  rsMatrix(int numRows, int numColumns, std::vector<T>&& newData) : data(std::move(newData))
  {
    rsAssert(newData.size() == 0);
    rsAssert(numRows*numColumns == data.size());
    this->numRows = numRows;
    this->numCols = numColumns;
    updateDataPointer();
  }

  rsMatrix(int numRows, int numColumns, std::initializer_list<T> l) : data(l) 
  {
    rsAssert(numRows*numColumns == l.size());
    this->numRows = numRows;
    this->numCols = numColumns;
    updateDataPointer();
  }
  // needs tests

  /** Copy constructor. Copies data from B into this object.  */
  rsMatrix(const rsMatrix& B)
  {
    setShape(B.numRows, B.numCols);
    rsArrayTools::copy(B.dataPointer, this->dataPointer, this->getSize());
  }

  /** Move constructor. Takes over ownership of the data stored in B. */
  rsMatrix(rsMatrix&& B) : data(std::move(B.data))
  {
    rsAssert(B.data.size() == 0); // B's data has now become our data
    this->numRows = B.numRows;
    this->numCols = B.numCols;
    updateDataPointer();
    B.reset();                    // invalidates pointer in B
  }

  /** Copy assignment operator. Copies data from rhs into this object. */
  rsMatrix<T, V>& operator=(const rsMatrix<T, V>& rhs)
  {
    if(this != &rhs) { // self-assignment check expected
      setShape(rhs.numRows, rhs.numCols);
      rsArrayTools::copy(rhs.dataPointer, this->dataPointer, this->getSize());
    }
    return *this;
  }

  /** Move assignment operator. Takes over ownership of the data stored in rhs. */
  rsMatrix<T, V>& operator=(rsMatrix<T, V>&& rhs)
  {
    data = std::move(rhs.data);
    rsAssert(rhs.data.size() == 0);
    this->numRows = rhs.numRows;
    this->numCols = rhs.numCols;
    updateDataPointer();
    rhs.reset();
    return *this;
  }


  /** Needs test! */
  /*
  rsMatrix<T>& operator=(const rsMatrixView<T>& rhs)
  {
    if (this != &rhs) { // self-assignment check expected
      setSize(rhs.getNumRows(), rhs.getNumColumns());
      rsArrayTools::copy(rhs.getDataPointerConst(), this->dataPointer, this->getSize());
    }
    return *this;
  }
  */





  /** Creates a zero matrix with given number of rows and columns. */
  static rsMatrix<T, V> zero(int numRows, int numColumns)
  { rsMatrix<T, V> Z(numRows, numColumns); Z.setToZero(); return Z; }

  /** Creates an identity matrix of given size. */
  static rsMatrix<T, V> identity(int size)
  { rsMatrix<T, V> E(size, size); E.setToIdentity(); return E; }

  /** Creates a diagonal matrix of given size with all diagonal values set to d. */
  static rsMatrix<T, V> diag(int size, const T& d)
  { rsMatrix<T, V> D(size, size); D.setDiagonalValues(d); return D; }

  // todo: diag(int size, T* data), diag(int size, T value)

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the number of rows and columns, this matrix should have. ToDo: provide a way to retain
  the data (optionally) - what does std::vector's resize do? Does it retain data...but if it does,
  it would be useless anyway in case the number of columns changed. */
  //void setSize(int numRows, int numColumns)
  void setShape(int numRows, int numColumns)
  {
    if(numRows == this->numRows && numColumns == this->numCols)
      return;  // nothing to do

    this->numRows = numRows;
    this->numCols = numColumns;
    data.resize(this->numRows * this->numCols);
    updateDataPointer();
    // optionally initialize with zeros
  }
  // rename to setShape - shall override inherited setShape

  /** Copies the data from another matrix into this one, converting the datatype, if necessarry. */
  template<class T2>
  void copyDataFrom(const rsMatrixView<T2>& A)
  {
    setShape(A.getNumRows(), A.getNumColumns());
    rsArrayTools::convert(A.getDataPointerConst(), this->dataPointer, this->getSize());
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Manipulations */

  /** Transposes this matrix, i.e. the rows become columns and vice versa. Avoids reallocation in
  case of square-matrices and row- and column vectors. */
  void transpose()
  {
    if(this->isSquare()) { rsMatrixView<T>::transposeSquare(this); return; }
    if(this->isVector()) { rsSwap(this->numRows, this->numCols); return; }

    V v(this->getSize());     // the general case needs reallocation...
    rsMatrixView<T> B(this->numCols, this->numRows, &v[0]);
    rsMatrixView<T>::transpose(*this, &B);
    rsSwap(this->numRows, this->numCols);
    data = v;
    // don't we have to call updateDataPointer()? -> check, if there's a unit test
  }
  // maybe move to cpp file



  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns a const reference to the std::vector that stores the data. Mostly for unit tests. */
  const std::vector<T>& getDataVectorConst() const { return data; }

  // todo: getDeterminant, getInverse, getFrobeniusNorm, get[Other]Norm, isPositiveDefinite,
  // getEigenvalues, getTrace, isUpperLeftTriangular, getTransposed, getConjugateTransposed


  //-----------------------------------------------------------------------------------------------
  /** \name Computations */

  /** Computes the Kronecker product between matrices A and B. For a 3x2 matrix A, it looks like:

              |a11*B a12*B|
  A (x) B  =  |a21*B a22*B|
              |a31*B a32*B|

  Where each entry aij*B is a submatrix of dimensions of B with the entries of B scaled by the
  respective element from A. This product is also sometimes called tensor product, but i think,
  Kronecker product is more appropriate, as it explicitly deals with matrices rather than vector
  spaces. Also, in the context of tensor algebra, the tensor product of two rank-2 tensors
  (i.e. matrices) would give a rank-4 tensor (with 4 indicies), whereas this product here still has
  just 2 indices - it is again a matrix and not some 4-dimensional "block". See here:
  https://en.wikipedia.org/wiki/Kronecker_product
  https://en.wikipedia.org/wiki/Tensor_product     */
  static rsMatrix<T, V> getKroneckerProduct(const rsMatrix<T, V>& A, const rsMatrix<T, V>& B)
  {
    rsMatrix<T, V> C(A.numRows*B.numRows, A.numCols*B.numCols);
    rsMatrixView<T>::kroneckerProduct(A, B, &C);
    return C;
  }

  /** Returns a matrix that is the element-wise product of this matrix an the right operand. */
  rsMatrix<T, V> getElementwiseProduct(const rsMatrixView<T>& rightOperand) const
  {
    rsMatrix<T, V> result(this->numRows, this->numCols);
    rsMatrixView<T>::elementwiseMultiply(*this, rightOperand, &result);
    return result;
  }
  // it's intentional that the rsMatrixView method is called multiply and this here is called
  // product - the view method just does something (-> use a verb), the method here returns an
  // objcet (-> use a noun)

  // todo: getInverse, getTranspose, getConjugate, getConjugateTranspose

  /** Returns a transposed version of this matrix. */
  rsMatrix<T, V> getTranspose() const
  {
    rsMatrix<T, V> t(this->numCols, this->numRows);
    rsMatrixView<T>::transpose(*this, &t);
    return t;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Decompositions */

  // getLowerUpperDecomposition ...or decomposeLU, decomposeQR, decomposeSVD 
  // ...but i think, they should go into rsLinearAlgebra


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Compares matrices for equality */
  bool operator==(const rsMatrix<T, V>& rhs) const
  {
    if(this->numRows != rhs.numRows || this->numCols != rhs.numCols)
      return false;
    return rsArrayTools::equal(this->dataPointer, rhs.dataPointer, this->getSize());
  }
  // maybe move to rsMatrixView, if possible

  /** Compares matrices for inequality */
  bool operator!=(const rsMatrix<T, V>& rhs) const { return !(*this == rhs); }

  /** Returns the negative of this matrix. */
  rsMatrix<T, V> operator-()
  {
    rsMatrix<T, V> C(this->numRows, this->numCols);
    for(int i = 0; i < this->getSize(); i++)
      C.dataPointer[i] = -this->dataPointer[i]; // maybe factor out into "neg" function in baseclass
    return C;
  }

  /** Adds two matrices: C = A + B. */
  rsMatrix<T, V> operator+(const rsMatrix<T, V>& B) const
  { rsMatrix<T, V> C(this->numRows, this->numCols); this->add(*this, B, &C); return C; }

  /** Subtracts two matrices: C = A - B. */
  rsMatrix<T, V> operator-(const rsMatrix<T, V>& B) const
  { rsMatrix<T, V> C(this->numRows, this->numCols); this->subtract(*this, B, &C); return C; }

  /** Multiplies two matrices: C = A * B. */
  rsMatrix<T, V> operator*(const rsMatrix<T, V>& B) const
  {
    //if(!areMultiplicable(*this, B))
    //  return rsMatrix<T>(0, 0); // return empty matrix when attempting to multiply incompatible matrices
    rsMatrix<T, V> C(this->numRows, B.numCols);
    this->matrixMultiply(*this, B, &C);
    return C;
  }
  // maybe it should return an empty matrix, when attempting to multiply incompatible matrices
  // ...or maybe we should throw an exception in such cases?


  /** Adds another matrix to this matrix and returns the result. */
  rsMatrix<T, V>& operator+=(const rsMatrix<T, V>& B)
  { this->add(*this, B, this); return *this; }

  /** Subtracts another matrix from this matrix and returns the result. */
  rsMatrix<T, V>& operator-=(const rsMatrix<T, V>& B)
  { this->subtract(*this, B, this); return *this; }

  /** Multiplies this matrix by another and returns the result. This is not an in-place process, i.e. it
  will allocate temporary heap-memory. */
  rsMatrix<T, V>& operator*=(const rsMatrix<T, V>& B)
  { *this = *this * B; return *this; }

  /** Multiplies this matrix with a scalar s: B = A*s. The scalar is to the right of the matrix. */
  rsMatrix<T, V> operator*(const T& s) const
  { rsMatrix<T, V> B(*this); B.scale(s); return B; }

  /** Multiplies this matrix by a scalar and returns the result. */
  rsMatrix<T, V>& operator*=(const T& s)
  { this->scale(s); return *this; }

  /** Divides this matrix by a scalar and returns the result. */
  rsMatrix<T, V>& operator/=(const T& s)
  { this->scale(T(1)/s); return *this; }



  /** Multiplies a matrix with a std::vector to give another vector: y = A * x. */
  std::vector<T> operator*(const std::vector<T>& x) const
  {
    rsAssert((int) x.size() == this->numCols, "Vector incompatible for left multiply by matrix");
    std::vector<T> y(this->numRows);
    for(int i = 0; i < this->numRows; i++) {
      y[i] = T(0);
      for(int j = 0; j < this->numCols; j++)
        y[i] += this->at(i, j) * x[j]; }
    return y;
  }
  // needs test - maybe optimize inner loop by avoiding re-computation of row base index

  /** Convolves this matrix with matrix B and returns the result. */
  rsMatrix<T, V> getConvolutionWith(const rsMatrix<T, V>& B)
  {
    rsMatrix<T, V> C(this->numRows + B.numRows - 1, this->numCols + B.numCols - 1);
    rsMatrixView<T>::convolve(*this, B, &C);
    return C;
  }

protected:

  /** Updates the data-pointer inherited from rsMatrixView to point to the begin of our std::vector
  that holds the actual data. */
  void updateDataPointer()
  {
    if(data.size() > 0)
      this->dataPointer = &data[0];
    else
      this->dataPointer = nullptr;
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  V data;  // V must be a vector type similar to std::vector

};

/** Multiplies a scalar and a matrix. */
template<class T, class V>
inline rsMatrix<T, V> operator*(const T& s, const rsMatrix<T, V>& A)
{
  rsMatrix<T, V> B(A);
  B.scale(s);
  return B;
}

/** Multiplies a row vector x with a matrix: y = x * A. The result y is another row vector: 
(MxN)*(PxQ) = MxQ with M=1 */
template<class T, class V>
std::vector<T> operator*(const std::vector<T>& x, const rsMatrix<T, V>& A)
{
  rsAssert((int) x.size() == A.getNumRows(), "vector incompatible for right multiply by matrix");
  std::vector<T> y(A.getNumColumns());
  for(int i = 0; i < A.getNumColumns(); i++) {
    y[i] = T(0);
    for(int j = 0; j < A.getNumRows(); j++)
      y[i] += x[j] * A(j, i); }
  return y;
}

template<class T>
rsMatrix<T, std::vector<T>> matrixMagnitudes(const rsMatrix<std::complex<T>>& A)  
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
rsMatrix<T, std::vector<T>> matrixPhases(const rsMatrix<std::complex<T>>& A)
{
  int N = A.getNumRows();
  int M = A.getNumColumns();
  rsMatrix<T> phases(N, M);
  for(int i = 0; i < N; i++)
    for(int j = 0; j < M; j++)
      phases(i, j) = arg(A(i, j));
  return phases;
}
// why do we have to be specific in the type of output vector, i.e. have to return a 
// rsMatrix<T, std::vector<T>> instead of rsMatrix<T, V> for some 2nd template parameter V 
// specifying the storage for the output? ...maybe because it cannot be deduced? but shouldn't
// it then fall back to std::vector?

// maybe factor out common code (keeping the above as covenience functions)...maybe something like 
// applyMatrixFunction with different input and output types for the template parameter T:

/*
template<class TIn, class TOut, class F>
rsMatrix<TOut> matrixFunction(const rsMatrix<TIn>& A, F func)
{
  int N = A.getNumRows();
  int M = A.getNumColumns();
  rsMatrix<TOut> out(N, M);
  for(int i = 0; i < N; i++)
    for(int j = 0; j < M; j++)
      out(i, j) = func(A(i, j));
  return out;
}
*/

// ...but how is the compiler supposed to infer TOut? maybe it should be a member function "apply"
// of rsMatrix<TOut>, so TOut can be infered from that


//-------------------------------------------------------------------------------------------------
// Notes:
// -design goals:
//  -use std::vector to hold the data in a flat array (so we can inspect it in the debugger)
//  -storage format should be compatible with lapack routines (maybe not by default, but can be
//   made so) - that means to support column major storage
// -maybe we should just work with our inherited dataPointer and use new/delete
//  -saves a little bit of storage
//  -but then we can't look easily at the data in the debugger anymore -> very bad! so, nope!
//  -the little storage overhead of std::vector becomes negligible for all but the smallest
//   matrices - on the other hand, very small matrices may be common
//  -maybe use std::vector in debug builds and new/delete in release builds
//  -maybe the numHeapAllocations variable can also be used only in debug builds? i guess we may
//   need to define some function incrementAllocationCounter that reduces to no-op in release
//   builds
//  -in expressions like rsMatrix<float> C = B*(A + B) + B; we want to avoid copying the data
//   unnecessarily - i.e. avoid that the temporaries that occur inside this expression use heap
//   allocation only when absolutely necessarry
//   ...this especially means, we need to pass the return values of the arithmetic operators by
//   reference rather than by value - typically, in an implementation like
//  -ideally, i want to be able to use it in production code for realtime processing...but that
//   may not be possible...but maybe with rsMatrixView, it is?
//
//    rsMatrix<T> operator+(const rsMatrix<T>& B) const
//    { rsMatrix<T> C(this->numRows, this->numCols); this->add(this, &B, &C); return C; }
//
//   we have one constructor call (-> heap allocation) to create C - but C is a local variable - to
//   return it to the caller, there would be *another* (copy?)constructor call - ...right? and that
//   second call is what we want to get rid of

// ToDo:
//  -maybe rename to rsMatrixOnHeap and have a similar class rsMatrixOnStack (maybe using
//   std::dynarray instead of std::vector)...or by defining the size via template parameters at
//   compile time ...maybe the std::vector vs std::dynarrray distinction can determined by passing
//   the storage container as template argument like so:
//   template<class ElemType, class ContainerType>
//   class rsMatrix { ContainerType<ElemType> data; };
//  -maybe make subclasses rsRowVector, rsColumnVector with a simplified element access operator
//   that takes only one index


// maybe see here:
// https://www.youtube.com/watch?v=PNRju6_yn3o
// https://www.ibm.com/developerworks/community/blogs/5894415f-be62-4bc0-81c5-3956e82276f3/entry/RVO_V_S_std_move?lang=en
// https://en.cppreference.com/w/cpp/language/copy_elision

#endif
