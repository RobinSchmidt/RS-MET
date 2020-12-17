// contains some unfinished code "scratches" and code under construction which may eventually be
// moved to somewhere else once it works and does something useful

using namespace std;

template<class T>
void cleanUpIntegers(T* a, int N, T tol)
{
  for(int i = 0; i < N; ++i) {
    T rounded = round(a[i]);
    if( rsAbs(a[i] - rounded) < tol )
      a[i] = rounded; }
}
// move to rsArray, make a version for complex numbers that does the same thing for real and
// imaginary parts separately

//-------------------------------------------------------------------------------------------------
// Linear Algebra stuff:

/** Returns rank of a matrix assumed to be in row echelon form. This is the number of nonzero
rows. */
template<class T>
int getRankRowEchelon(const rsMatrixView<T>& A, T tol)
{
  // rsAssert(isRowEchelon(A));
  int i = 0;
  while(i < A.getNumRows()) {
    bool nonZeroElemFound = false;
    int j = i;
    while(j < A.getNumColumns()) {
      if( rsGreaterAbs(A(i, j), tol) )  {
        nonZeroElemFound = true;
        break; } // i-th row is not all-zeros
      j++; }
    if(!nonZeroElemFound)
      return i;  // no non-zero element was found - i-th row is all zeros
    i++; }
  return i;
}
// verify, if this is correct - maybe make unit test with weird matrices - actually this may be
// overly complicated - we could just use getNumNonZeroRows

/** Returns true, if the space spanned by the columns of x is within the span of the columns of B.
That means each column of x can be expressed as some linear combination of the columns of B. */
template<class T>
bool isInSpanOf(rsMatrix<T> B, rsMatrix<T> x, T tol)
{
  //rsAssert(B.hasSameShapeAs(x));
  rsAssert(B.getNumRows() == x.getNumRows());
  RAPT::rsLinearAlgebraNew::makeTriangular(B, x);  // use rowEchelon2
  int rankB = getRankRowEchelon(B, tol);
  return x.areRowsZero(rankB, x.getNumRows()-1, tol);
}

template<class T>
bool spanSameSpace(const rsMatrix<T>& A, const rsMatrix<T>& B, T tol)
{
  return isInSpanOf(A, B, tol) && isInSpanOf(B, A, tol);
}
// needs test

// we replicate the function from rsLinearAlgebraNew  here but taking full rsMatrix objects instead
// of views, so we can expand the content
template<class T>
bool makeTriangularNoPivot(rsMatrix<T>& A, rsMatrix<T>& B)
{
  rsAssert(A.isSquare()); // can we relax this?
  int N = A.getNumRows();
  for(int i = 0; i < N; i++) {
    if(A(i, i) == T(0)) {
      rsError("This matrix needs pivoting"); return false; }
    for(int j = i+1; j < N; j++) {
      T s = -A(j, i) / A(i, i);
      A.addWeightedRowToOther(i, j, s, i, A.getNumColumns()-1); // start at i: avoid adding zeros
      B.addWeightedRowToOther(i, j, s); }}
  return true;
}
// when using rational numbers or functions, maybe the pivoting could based on which element is the
// simplest (lowest degree or smallest) nonzero element - because then, there may be less complex
// computations


template<class T>
int getPivotRow(const rsMatrixView<T>& A, int row, int column)
{
  T biggest = T(0);
  int pivRow = row;
  for(int i = row; i < A.getNumRows(); i++)  {
    if( rsGreaterAbs(A(i, column), biggest) ) {
      biggest = A(i, column);
      pivRow = i; }}
  return pivRow;
}

template<class T>
int getLeadCoeffIndex(const rsMatrixView<T>& A, int row, T tol, int startColumn = 0)
{
  int j;
  for(j = startColumn; j < A.getNumColumns(); j++)
    if( rsGreaterAbs( A(row, j), tol ) )
      return j;
  return j;  // will return numCols (i.e. invalid index) when the row is all zeros
}


/** Puts the augmented coefficient matrix A|B into row echelon form via Gaussian elimination.
In this form, each row has its leading coefficient at least one position further to the right than
the previous row and rows of all zeros are at the bottom. */
template<class T>
void rowEchelon(rsMatrixView<T>& A, rsMatrixView<T>& B, T tol)
{
  //bool reduced = false; // make parameter - switch for producing *reduced* row echelon form

  int i = 0;
  int j = 0;
  while( i < A.getNumRows() && j < A.getNumColumns() )
  {
    int p = getPivotRow(A, i, j);
    if(i != p) {
      A.swapRows(i, p);
      B.swapRows(i, p); }
    if(!rsGreaterAbs(A(i, j), tol)) {   // column of all zeros encountered...
      j++; continue;    }               // ...try with next column but same row

    j = getLeadCoeffIndex(A, i, tol, j);
    if(j >= A.getNumColumns())
      break; // no pivot was found anymore - we are done

    // pivot row subtraction
    for(int k = i+1; k < A.getNumRows(); k++) {
      T w = -A(k, j) / A(i, j);     // A(i, j) is now known to be nonzero - right?
      A.addWeightedRowToOther(i, k, w, j, A.getNumColumns()-1);
      B.addWeightedRowToOther(i, k, w); }

    //// experimental - is that correct? - do tests!
    //if(reduced)
    //{
    //  for(int k = i-1; k >= 0; k--) {
    //    T w = -A(k, j) / A(i, j);
    //    A.addWeightedRowToOther(i, k, w, 0, j);
    //    B.addWeightedRowToOther(i, k, w); }
    //}

    i++;
  }

}
// we are finished when no pivots can found in the current row anymore
// needs test
// todo: 
// -maybe let the function getPivotRow be a functor that can be passed in by client code
// -we may want to use different pivot-search strategies dependign on the type T - if it's a 
//  floating point type, we want to use the element with largest absolute value to reduce roundoff
//  error, for rsFraction, we want to use a nonzero number with small denominator to avoid blowup
//  and overflow of num and den (that can also be solved by providing explicit specializations for 
//  different types T - however, using a functor is more flexible - it allows client code to change
//  the strategy)
// -maybe implement a version that does full pivoting...i think, that turns the algo form O(N^3) to
//  O(N^4)...right? or not? no! it makes the pivot search O(N^2) and that is called in a simple
//  loop - so we remain at O(N^3), just get a larger constant factor

template<class T>
void rowEchelon(rsMatrixView<T>& A, T tol)
{
  rsMatrix<T> dummy(A.getNumRows(), 1);
  rowEchelon(A, dummy, tol);
}
// allocates because of dummy - try to get rid of the allocation

// make a function reducedRowEchelon that also includes a backward/upward elimination pass - maybe
// that can be integrated as option into the function above - if reduced == true also eliminate
// upward - does that work



template<class T>
RAPT::rsPolynomial<T> getCharacteristicPolynomial(const rsMatrixView<T>& A, T tol)
{
  using RatFunc = RAPT::rsRationalFunction<T>;
  using Matrix  = RAPT::rsMatrix<RatFunc>;

  //T tol = 1.e-8;  // make parameter
  //T tol = 1.e-12;  // make parameter
  //T tol = T(0);  // make parameter

  // Create matrix B = A - x*I as matrix of rational functions:
  Matrix B(A.getNumRows(), A.getNumColumns());
  for(int i = 0; i < B.getNumRows(); ++i)
    for(int j = 0; j < B.getNumColumns(); ++j)
      if(i == j)
        B(i, j) = RatFunc({A(i, j), -1}, {1}, tol);  // function (A(i, j) - x) / 1 on the diagonal
      else
        B(i, j) = RatFunc({A(i, j)},     {1}, tol);  // constant function A(i,j) / 1 off the diagonal

  // Create a dummy right-hand-side (todo: allow function to be called without rhs) - maybe
  // make an LU decomposition function that fills an array of swaps at each index, the number says
  // with which row this row has to be swapped
  Matrix R(A.getNumRows(), 1);
  for(int i = 0; i < R.getNumRows(); ++i)
    R(i, 0) = RatFunc({ T(1) }, { T(1) }, tol);

  // compute row echelon form of B:
  //RAPT::rsLinearAlgebraNew::makeTriangularNoPivot(B, R);
  makeTriangularNoPivot(B, R);
  // i think, we really need pivoting for this here too - we may encounter the zero-function - but
  // maybe we should swap only when the zero-function is encountered - i.e. we don't search for the
  // largest element - instead, we check against zero and if we encounter a zero, we search for the
  // next nonzero element to swap with - the current code should not go into the library - it's
  // useless for production - see weitz book pg 310
  // or use rowEchelon2

  // Compute determinant. For a triangular matrix, this is the product of the diagonal elements.
  // The computed determinant is still a rational function but it should come out as a polynomial,
  // i.e. the denominator should have degree 0 (be a constant). I think, it should always be +1 or
  // -1 because the elementary row operations can only flip the determinant.
  RatFunc d = B.getDiagonalProduct();
  //d.reduce(1.e-8);   // test - doesn't seem to help
  rsAssert(d.getDenominatorDegree() == 0);
  return d.getNumerator() / d.getDenominator()[0];
}
// todo: make a version that uses Laplace expansion of the determinant with a matrix of polynomials
// -> avoids use of rsRationalFunction, needs only rsPolynomial -> should give the same result (but
// this is only for testing, not for production - Laplace expansion is ridiculously expensive)


/** Represents the eigenspace of a matrix with complex coefficients. Each eigenspace consists of
an eigenvalue and an associated set of eigenvectors represented as columns of a matrix. The columns
can be seen as basis vectors that span the eigenspace associated with the given eigenvalue. */
template<class T>
struct rsEigenSpace
{
  int getAlgebraicMultiplicity() const { return algebraicMultiplicity; }
  int getGeometricMultiplicity() const { return (int) eigenSpace.getNumColumns(); }
  //void orthonormalize();

  complex<T> eigenValue;
  int algebraicMultiplicity;
  rsMatrix<complex<T>> eigenSpace; // basis of nullspace of A - eigenvalue * I
};

/** Takes a matrix of real numbers and turns the elements into complex numbers. Needed for technical
reasons in eigenvector computations - mainly because the eigenvalues of a real matrix may
nevertheless be complex, so we need to lift some computations on real matrices into the complex
domain. */
template<class T>
RAPT::rsMatrix<complex<T>> complexify(const RAPT::rsMatrix<T>& A)
{
  RAPT::rsMatrix<complex<T>> Ac(A.getNumRows(), A.getNumColumns());
  Ac.copyDataFrom(A);
  return Ac;
}



/*
// meant to be used from the debugger:
template<class T>
void findEigenSpacesReal(const RAPT::rsMatrix<T>& A) // for real matrices - includes complexification
{
  //double tol = 1.e-12;  // make parameter
  complex<T> tol = 1.e-12;  // make parameter

  rsAssert(A.isSquare());  // try to relax later
  using Matrix  = RAPT::rsMatrix<T>;
  using MatrixC = RAPT::rsMatrix<complex<T>>;
  vector<complex<T>> eigenvalues = getEigenvalues(A);
  int N = A.getNumRows();
  MatrixC Ac = complexify(A);
  MatrixC I(N,N); I.setToIdentity();
  vector<MatrixC> eigenspaces(N);
  for(int i = 0; i < N; i++) {
    MatrixC Ai = Ac - eigenvalues[i] * I;

    eigenspaces[i] = getNullSpace(Ai, tol);
    //eigenspaces[i] = getNullSpaceTailParams(Ai, tol);

    // test, if A * v = eigenvalue * v for v in the eigenspace at once
    MatrixC test1 = Ac * eigenspaces[i];
    MatrixC test2 = eigenvalues[i] * eigenspaces[i];
    MatrixC error = test1 - test2;
    //rsAssert(error.isZero(tol));

    eigenspaces[i] = eigenspaces[i].getTranspose(); // for convenient inspection in the debugger
    int dummy = 0;
  }                                                 // todo: have a function that transpose in place

  // maybe getEigenValues should return an array of root objects in which each root is unique and
  // has a multiplicity - maybe we should have a functioon getUniqueRoots

  // This function is currently only good for inspecting things in the debugger - but eventually
  // it should return something - but what should the datastructure look like? and should we
  // perhaps take a complex matrix as input?

  int dummy = 0;
}
// move useful parts of the code eleswhere, then delete
*/


// get rid of the duplication:
template<class T>
vector<complex<T>> getPolynomialRoots(const RAPT::rsPolynomial<T>& p)
{
  vector<complex<T>> roots(p.getDegree());
  RAPT::rsPolynomial<T>::roots(p.getCoeffPointerConst(), p.getDegree(), &roots[0]);
  return roots;
}

template<class T>
vector<complex<T>> getPolynomialRoots(const RAPT::rsPolynomial<complex<T>>& p)
{
  vector<complex<T>> roots(p.getDegree());
  RAPT::rsPolynomial<T>::roots(p.getCoeffPointerConst(), p.getDegree(), &roots[0]);
  return roots;
}
// maybe make member getRoots - that would be convenient - but the problem is that even for real
// polynomials, the roots may be complex - so i don't know if it should return a vector<T> or a
// vector<complex<T>> - if the polynomial is itself complex, the former would be the way to go,
// if it's real, the latter ...maybe make two functions getRootsReal, getRootsComplex - or maybe
// in the case of real polynomials, it should return only the real roots unless one calls
// getComplexRoots - or maybe, if we implement it as template once for reals and once for
// complexes, the right one will be compiled into the class automatically? ..like
// vector<complex<T>> getRoots(const RAPT::rsPolynomial<T>& p); // T is real type
// vector<T> getRoots(const RAPT::rsPolynomial<complex<R>>& p); // T is complex, R is real


template<class R>   // R is a real-number datatype (float, double, etc.)
vector<complex<R>> getEigenvalues(const rsMatrixView<R>& A, R tol)
{
  RAPT::rsPolynomial<R> p = getCharacteristicPolynomial(A, tol);
  return getPolynomialRoots(p);
}

/** Takes a vector of complex numbers and returns a vector of their real parts. */
template<class R>
std::vector<R> getRealParts(const std::vector<std::complex<R>>& v)
{
  std::vector<R> r(v.size());
  for(size_t i = 0; i < v.size(); i++)
    r[i] = v[i].real();
  return r;
}

/** Returns the real parts of the eigenvalues of matrix A. This function may make sense when you
know that the eigenvalues are real anyway, because - for example - A is symmetric. */
template<class R>   // R is a real-number datatype (float, double, etc.)
std::vector<R> getEigenvaluesReal(const rsMatrixView<R>& A, R tol)
{
  std::vector<std::complex<R>> evc = getEigenvalues(A, tol); // complex eigenvalues
  // maybe we should assert that the imainary parts of evc are all zero
  return getRealParts(evc);
}
// todo: maybe a special algorithm to find real roots of a real polynomial can be used - avoid
// the intermediate complexification of the roots (-> optimization)





  // see also:
  // https://lpsa.swarthmore.edu/MtrxVibe/EigMat/MatrixEigen.html
  // https://www.scss.tcd.ie/Rozenn.Dahyot/CS1BA1/SolutionEigen.pdf
  // http://www.sosmath.com/matrix/eigen2/eigen2.html
  // http://wwwf.imperial.ac.uk/metric/metric_public/matrices/eigenvalues_and_eigenvectors/eigenvalues2.html

// http://doc.sagemath.org/html/en/constructions/linear_algebra.html


/** This changes the matrices A,B in random ways but without changing the solution set to
A * X = B. Can be used for investigations on numerical precision issues in matrix
computations. We can construct a matrix from a diagonal matrix by shuffling - for the diagonal
matrix, the exact solutions are easy to compute..... */
template<class T>
void shuffle(rsMatrix<T>& A, rsMatrix<T>& B, int range, int seed = 0)
{
  rsNoiseGenerator<T> prng;

  T w; int i;

  // downward sweep:
  for(i = 0; i < A.getNumRows()-1; ++i) {
    //plotMatrix(A, true);

    if(prng.getSample() >= T(0))
      w = T(+1);
    else
      w = T(-1);
    //w = round(range * prng.getSample());
    A.addWeightedRowToOther(i, i+1, w);
    B.addWeightedRowToOther(i, i+1, w);
  }

  // upward sweep:
  for(i = A.getNumRows()-1; i > 0; --i)
  {
    //plotMatrix(A, true);

    if(prng.getSample() >= T(0))
      w = T(+1);
    else
      w = T(-1);
    //w = prng.getSample();
    //w = round(range * prng.getSample());
    A.addWeightedRowToOther(i, i-1, w);
    B.addWeightedRowToOther(i, i-1, w);
  }

  // maybe make also a rightward and leftwar sweep
}

template<class T>
rsMatrix<T> getSubMatrix(
  const rsMatrixView<T>& A, int startRow, int startCol, int numRows, int numCols)
{
  // assert that it all makes sense

  rsMatrix<T> S(numRows, numCols);
  for(int i = 0; i < numRows; ++i)
    for(int j = 0; j < numCols; ++j)
      S(i, j) = A(startRow+i, startCol+j);
  return S;
}
// make member of rsMatrix

/** If A is an MxN matrix, this function returns the (M-1)x(N-1) matrix that results from removing
th i-th row and j-th column. See:
https://en.wikipedia.org/wiki/Adjugate_matrix  */
template<class T>
rsMatrix<T> getAdjugate(const rsMatrix<T>& A, int i, int j)
{
  // assert that it all makes sense
  int M = A.getNumRows();
  int N = A.getNumColumns();
  rsMatrix<T> Aij(M-1, N-1);
  int ii, jj;                                                    // indices into A matrix
  for(ii = 0; ii < i; ++ii) {                                    // loop over top rows
    for(jj = 0;   jj < j; ++jj)  Aij(ii,   jj  ) = A(ii, jj);    //   loop over left columns
    for(jj = j+1; jj < N; ++jj)  Aij(ii,   jj-1) = A(ii, jj); }  //   loop over right columns
  for(ii = i+1; ii < N; ++ii) {                                  // loop over bottom rows
    for(jj = 0;   jj < j; ++jj)  Aij(ii-1, jj  ) = A(ii, jj);    //   loop over left columns
    for(jj = j+1; jj < N; ++jj)  Aij(ii-1, jj-1) = A(ii, jj); }  //   loop over right columns
  return Aij;
}
// make member of rsMatrix, the ii,jj are ugly - replace with i,j - maybe rename current i,j into
// row, col


// for computing nullspaces, we will need a function that removes zero columns - it sould take an
// array of th indices of the column to remove
// rsMatrix<T> getWithRemovedColumns(const rsMatrix<T>& A, const std::vector<int>& colIndices)

/** Expands the determinant column-wise along the j-th column via the Laplace expansion method.
This method has extremely bad scaling of the complexity. The function calls itself recursively in a
loop (!!!). I think, the complexity may scale with the factorial function (todo: verify) - which
would be super-exponential - so it's definitely not meant for use in production code. For
production, use Gaussian elimination (we need to keep track of whether we have an odd or even
number of swaps - in the former case det = -1 * product(diagonal-elements of upper triangular
form), in the later case det = +1 * product(...) */
template<class T>
T getDeterminantColumnWise(const rsMatrix<T>& A, int j = 0)
{
  rsAssert(A.isSquare());
  int N = A.getNumRows();
  if(N == 1)
    return A(0, 0);
  T det = T(0);
  for(int i = 0; i < N; i++) {
    rsMatrix<T> Aij = getAdjugate(A, i, j);
    det += pow(-1, i+j) * A(i, j) * getDeterminantColumnWise(Aij, 0); }
  return det;

  // Notes:
  // -in the recursive call we always expand along the 0-th column to make it work also in the
  //  base-case which we will eventually reach
  // -the function could be optimized by implementing special rules for 2x2 and 3x3 matrices:
  //  https://en.wikipedia.org/wiki/Determinant#2_%C3%97_2_matrices to avoid the lowest level of
  //  recursion - this was not done because it's not meant for production use anyway - it's
  //  insanely inefficient anyway
  // -the formulas can be found for example here
  //  https://en.wikipedia.org/wiki/Determinant#Laplace's_formula_and_the_adjugate_matrix
  // -the function could be used to find the characteristic polynomial without resorting to
  //  rational functions - only class rsPolynomial itself is needed. For this, we must promote the
  //  given matrix of numbers to a matrix of polynomials instead of rational functions.
}
// rename to getDeterminantViaLaplaceExpansion

template<class T>
T getDeterminantRowWise(const rsMatrix<T>& A, int i = 0)
{
  rsAssert(A.isSquare());
  int N = A.getNumRows();
  if(N == 1) return A(0, 0);
  T det = T(0);
  for(int j = 0; j < N; j++) {
    rsMatrix<T> Aij = getAdjugate(A, i, j);
    det += pow(-1, i+j) * A(i, j) * getDeterminantRowWise(Aij, 0); }
  return det;
}


template<class T>
int getNumNonZeroRows(const rsMatrixView<T>& A, T tol)
{
  int i = A.getNumRows() - 1;
  while(i > 0) {
    if(!A.isRowZero(i, tol))
      return i+1;
    i--; }
  return i;
}
// function assumes that zero rows can occur only at the bottom - the name should reflect that
// maybe make a function getNumBottomZeroRows instead

template<class T>
int getNumBottomZeroRows(const rsMatrix<T>& A, T tol)
{
  for(int i = A.getNumRows()-1; i >= 0; i--)
    if(!A.isRowZero(i, tol))
      return A.getNumRows() - i;
  return 0;
}



template<class T>
rsMatrix<T> getWithoutBottomZeroRows(const rsMatrix<T>& A, T tol)
{
  int rank = getRankRowEchelon(A, tol); // A is assumed to be in row-echelon form
  return getSubMatrix(A, 0, 0, rank, A.getNumColumns());
}
// make member of rsMatrix

/** Returns the space spanned by the rows of matrix A... see Karpf. pg 140 */
template<class T>
rsMatrix<T> getRowSpace(rsMatrix<T> A, T tol)
{
  rsMatrix<T> z(A.getNumRows(), 1);    // dummy
  RAPT::rsLinearAlgebraNew::makeTriangular(A, z);
  return getWithoutBottomZeroRows(A, tol);
}
// needs more tests

template<class T>
rsMatrix<T> getColumnSpace(rsMatrix<T> A, T tol)
{
  return getRowSpace(A.getTranspose(), tol).getTranspose();
}

/** Returns a matrix whose columns are a basis of the nullspace (a.k.a. kernel) of the matrix A.
The basis is not orthogonal or normalized. If the nullspace contains only the zero vector, an
empty matrix is returned. This function returns wrong results when there are leading columns of
zeros in the row-echelon form of A - this can be tested with matrices that are already in this
form (in which case LA::makeTriangularLA::makeTriangular will do nothing and return 0).

...soo i think, this means that this function can be used safely for regular matrices but for
singular ones, it may or may not fail? ...or will it always fail for singluar matrices? */
template<class T>
rsMatrix<T> getNullSpaceTailParams(rsMatrix<T> A, T tol)
{
  //T tol = T(1.e-12);  // make parameter

  using Matrix = RAPT::rsMatrix<T>;
  using LA     = RAPT::rsLinearAlgebraNew;

  int numRows  = A.getNumRows();          // dimensionality of input space
  int numCols  = A.getNumColumns();       // dimensionality of output space
  Matrix z(numRows, 1);                   // dummy - needed by function
  LA::makeTriangular(A, z);               // reduces A to row echelon form
  int rank = getRankRowEchelon(A, tol);   // rank, dimensionality of image
  int nullity = numCols - rank;           // dimensionality of nullspace (see karpf. 142)

  // extract rank x rank system with nullity rhs vectors
  Matrix S = getSubMatrix(A, 0, 0, rank, rank);
  Matrix R = Matrix(rank, nullity);
  for(int j = 0; j < nullity; ++j)   // loop over the rhs
    for(int i = 0; i < rank; ++i)    // loop over rows in current rhs
      R(i, j) = -A(i, rank+j);

  // compute first nullity elements of basis vectors:
  Matrix b = Matrix(rank, nullity);
  LA::solveTriangular(S, b, R);

  // compute filled up basis vectors:
  Matrix B = Matrix(A.getNumColumns(), nullity);
  B.setToZero();
  for(int j = 0; j < nullity; ++j) {
    for(int i = 0; i < rank; i++)
      B(i, j) = b(i, j);
    B(rank+j, j) = 1; }

  return B;
}






// http://www.cfm.brown.edu/people/dobrush/am34/sage/kernel.html


template<class T>
bool isRowEchelon(const rsMatrixView<T>& A, T tol)
{
  int col = -1;
  for(int i = 0; i < A.getNumRows(); i++) {
    int j = getLeadCoeffIndex(A, i, tol);
    if(j <= col && j < A.getNumColumns()) // there was no right-step in this row, the 2nd condition
      return false;                       // is needed for dealing with zero-rows at the bottom
    col = j;
  }
  return true;
}
// for each step down, we must make at least one step right and we never make steps back left
// not yet tested



// https://www.wikihow.com/Find-the-Null-Space-of-a-Matrix
// it says:
// the pivots - the leading coefficients - rest in columns 1 and 3. That means that x1,x3 have
// their identifying equations. The result is that x2,x4,x5 are all free variables.
// -> we cannot just freely choose, which of the variables we use as free parameters - that seems
// to be the key: make a first pass to identify the pivots - these are the columns of the leading
// coefficients - the free parameters are then all the other variables - set up a system for these
// other variables in terms of the free ones and solve the system as many times as there are free
// variables
// http://www.eng.fsu.edu/~dommelen/aim/style_a/GEspc.html
// https://yutsumura.com/how-to-find-a-basis-for-the-nullspace-row-space-and-range-of-a-matrix/

// assumes A to be in row echelon form - maybe write a function that checks that and make an assert
template<class T>
std::vector<int> getPivots(const rsMatrix<T>& A, T tol)
{
  rsAssert(isRowEchelon(A, tol), "makes sense only for row echelon matrices");
  std::vector<int> pivots; // maybe pre-allocate - the number is known before
  int row = 0;
  while(true) {
    if( row >= A.getNumRows() )
      break;
    int leadCoeffIndex = getLeadCoeffIndex(A, row, tol);
    if( leadCoeffIndex >= A.getNumColumns() )
      break;
    pivots.push_back(leadCoeffIndex);
    row++; }
  return pivots;
}
// needs more tests

template<class T>
std::vector<int> getNonPivots(const rsMatrix<T>& A, T tol)
{
  rsAssert(isRowEchelon(A, tol), "makes sense only for row echelon matrices");
  std::vector<int> nonPivots;
  int i = 0, j = 0;
  while(i < A.getNumRows()) {
    while(j < A.getNumColumns() && !rsGreaterAbs(A(i, j), tol)) {
      nonPivots.push_back(j);
      j++;  }
    i++;
    j++; }
  return nonPivots;
}
// needs more tests


/** Returns true iff array A contains element x exactly once. */
template<class T>
bool containsOnce(const T* A, int N, T x)
{
  bool found = false;
  for(int i = 0; i < N; ++i) {
    if(found == true && A[i] == x)
      return false;       // x was found twice
    if(A[i] == x)
      found = true; }
  return found;
}
// move to rsArrayTools

/** Returns true, if each of the indices from 0 to numIndices-1 (both inclusive) are contained
exactly once in the indices array. That means the array is a permutation of the numbers from
0 to numIndices-1. */
bool isIndexPermutation(const int* indices, int numIndices)
{
  for(int i = 0; i < numIndices; ++i)
    if( !containsOnce(indices, numIndices, i) )
      return false;
  return true;
}

/** Returns true, iff every number between 0 and N-1 (both included) is contained exactly once in
either subset1 or subset2. This means the two subsets split the set of indices from 0 to N-1 into
two disjoint sets that together combine to the whole set. */
bool isIndexSplit(int numIndices, const int* subset1, int size1, const int* subset2, int size2)
{
  if(size1 + size2 != numIndices)
    return false;
  for(int i = 0; i < numIndices; ++i)
    if( !(containsOnce(subset1, size1, i) || containsOnce(subset2, size2, i)) )
      return false;
  return true;
}
// maybe make a more general function that instead of checking for a split of indices checks for a
// split of arbitrary things - it would have to take another array for the full set as parameter
// and we would have to do containsOnce(subset1, size1, fullSet[i]), etc. - maybe call it isSplit
// or isDisjointSplit

bool isIndexSplit(int numIndices, const std::vector<int> subset1, const std::vector<int> subset2)
{
  if(subset1.size() == 0)
    return isIndexPermutation(&subset2[0], numIndices);
  if(subset2.size() == 0)
    return isIndexPermutation(&subset1[0], numIndices);
  return isIndexSplit(numIndices,
    &subset1[0], (int)subset1.size(),
    &subset2[0], (int)subset2.size());
}

/** If the coefficient matrix A has rows full of zeros at the bottom, the system is singular. If
the augment B has also a zero row for each of the zero rows in A, the singular system is
consistent, so there are infinitely many solutions, so we get to choose some free parameters. We
assume that to be the case here - the matrices A and B are assumed to both have rankA nonzero rows.
We make the choice that the bottom elements in the solution vectors in X should be zero. This
amounts to solving the sub-system with the top-left section of the original matrix (setting the
bottom variables zero renders the right section of the matrix ineffective - for other choices, we
would have to adapt the right-hand side). So, this function returns a particular solution from the
infinietly many. The full set of solutions is given by that (or another) particular solution plus
any linear combination of the basis vectors of the matrix A's nullspace. */
template<class T>
bool solveUnderDeterminedRowEchelon(
  rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B, int rankA, T tol)
{
  rsAssert(isRowEchelon(A, tol), "A is assumed to be in row echelon form - but isn't");
  rsMatrix<T> a = getSubMatrix(A, 0, 0, rankA, rankA);
  rsMatrix<T> b = getSubMatrix(B, 0, 0, rankA, B.getNumColumns());
  rsMatrix<T> x(b.getNumRows(), b.getNumColumns());
  RAPT::rsLinearAlgebraNew::solveTriangular(a, x, b); // A is in echelon form - and so is a
  X.setToZero();
  for(int i = 0; i < x.getNumRows(); ++i)
    for(int j = 0; j < x.getNumColumns(); ++j)
      X(i, j) = x(i, j);
  // maybe factor out into X.copySubMatrixFrom(x, 0, 0, rankA, X.getNumColumns())
  return true;
}
// todo: make another function that computes a minimum-norm solution

template<class T>
bool solve2(rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B, T tol)
{
  rsAssert(A.getNumColumns() == X.getNumRows()); // A*X = B or A*x = b must make sense
  rsAssert(X.hasSameShapeAs(B));                 // num of solutions == num of rhs-vectors
  rsAssert(A.isSquare());                        // do we really need this? maybe not!
  rowEchelon(A, B, tol);
  int rankA  = getRankRowEchelon(A, tol); // number of nonzero rows, rank of the coeff matrix
  int rankAB = getNumNonZeroRows(B, tol); // same for the augmented coeff matrix A|B
  if(rankA == A.getNumColumns()) {                           // system was regular
    RAPT::rsLinearAlgebraNew::solveTriangular(A, X, B);      //   -> unique solution
    return true; }
  else  {                                                    // system was singular...
    if(rankAB > rankA)                                       // ...and inconsistent/overdetermined
      return false;                                          //        -> no solution possible
    else {
      solveUnderDeterminedRowEchelon(A, X, B, rankA, tol);   // ...and consistent/underdetermined
      return true;  }}                                       //        -> infinitely many solutions
}
// -maybe return an int: 0: no solution, 1: unique solution, 2: many solutions
// -maybe have a function getSolutionSet



/** Computes the set of vectors v which solve the homogenous linear system of equations A * v = 0
where v is some vector and 0 is the zero vector. We assume v to be M-dimensional, so the matrix A
must have M columns. The set of vectors v that solve this equation will in general span a subspace
of R^M. This subspace is called the nullspace of the matrix A. This function returns a basis for
this subspace represented as matrix. The columns of the matrix are the basis vectors. Note that
this may be the empty matrix which indicates that the nullspace of A consists only of the
zero-vector (todo: maybe we should return the Mx1 zero-vector in this case? ...decide later by
which convention is more convenient when dealing with eigenspaces */
template<class T>
rsMatrix<T> getNullSpace(rsMatrix<T> A, T tol)  // maybe to needs to be a separate type for complex matrices
{
  // Algorithm:
  // We bring the matrix into row-echelon form and figure out its rank and nullity. The nullity
  // gives the number of basis vectors that we must produce. Then we split (conceptually, not
  // literally) the M columns of A into those which correspond to our free parameters and those
  // which correspond to the dependent variables (they depend on the choice we make for our free
  // parameters). The number of dependent variables is equal to the rank R of the matrix and the
  // number of free parameters N gives the dimensionality of the nullspace - this is also called
  // the nullity of the matrix. By the rank-nullity theorem, they must sum up to M: M = R + N,
  // which says that the dimensionality of the embedding R^M vector space equals the dimensionality
  // of the A's nullspace plus the dimensionality of A's column space (or is it the row-space? the
  // space spanned by the rows makes more sense because the rows live in R^M while the columns may
  // not -> FIGURE OUT). We set up an RxR linear system and solve it for N different right hand
  // sides corresponding to N different choices for assigning the free parameters. The most natural
  // choice is to set one to 1 and all others to 0 in each assignment and select a different
  // one to set to 1 in each of the N cases. The solution of the linear system gives us
  // R elements for each of the N basis vectors. The remaining N elements must the be filled up
  // with ones and zeros according to our choices for the parameter assignments. To set up the
  // system and to combine the solution, we use our pivots and params arrays to gather and scatter
  // the numbers. See:
  // https://www.wikihow.com/Find-the-Null-Space-of-a-Matrix
  // http://www.eng.fsu.edu/~dommelen/aim/style_a/GEspc.html
  // https://en.wikipedia.org/wiki/Rank%E2%80%93nullity_theorem

  using Matrix = RAPT::rsMatrix<T>;
  rowEchelon(A, tol);

  // maybe factor out a function that takes a row-echelon matrix getNullSpaceEchelon(Matrix&)

  // find out which dimensions are free and which dependent:
  std::vector<int> pivots = getPivots(   A, tol);  // indices of dependent variables
  std::vector<int> params = getNonPivots(A, tol);  // indices of free variables (parameters)
  int nEqn = (int) pivots.size();                  // number of equations R (#dependents, rank)
  int nRhs = (int) params.size();                  // number of right-hand sides N (#parameters)

  // sanity checks for debug:
  rsAssert(isIndexSplit(A.getNumColumns(), pivots, params));
  rsAssert(getRankRowEchelon(A, tol) == (int) pivots.size());

  // set up the linear system ("gather") and solve it:
  Matrix M(nEqn, nEqn);                   // coefficient matrix of the NxN system
  Matrix R(nEqn, nRhs);                   // right hand side matrix
  Matrix b(nEqn, nRhs);                   // solution
  int i, j;
  for(i = 0; i < nEqn; i++) {
    for(j = 0; j < nEqn; j++)
      M(i, j) =  A(i, pivots[j]);         // copy relevant coeffs
    for(j = 0; j < nRhs; j++)
      R(i, j) = -A(i, params[j]); }       // resulting rhs from setting j-th param to 1
  bool success = solve2(M, b, R, tol);
  rsAssert(success);                      // this should never fail - right?

  // write solutions into output ("scatter") and fill up with ones and zeros:
  Matrix B(A.getNumColumns(), nRhs);               // final result
  B.setToZero();
  for(i = 0; i < nEqn; i++)
    for(j = 0; j < nRhs; j++)
      B(pivots[i], j) = b(i, j);
  for(i = 0; i < nRhs; i++)
    for(j = 0; j < nRhs; j++)
      //B(params[i], j) = 1;
      B(params[i], i) = 1;  // should it be params[i], j ? ...no!

  return B;
}
// move to library, rename to nullSpace, also move rowSpace and columnSpace over
// it sometimes produces bases that contain the zero vector

/** Structure for representing the occurence of a value together with its multiplicity (i.e. number
of times, it occurred). Useful, for example, for polynomial roots */
template<class T>
struct rsOccurrence
{
  rsOccurrence(T _value, int _multiplicity) : value(_value), multiplicity(_multiplicity) {}
  T value;
  int multiplicity;
};

template<class TItem, class TTol>
int rsFindOccurrence(const rsOccurrence<TItem>* items, int numItems, const TItem& item, TTol tol)
{
  for(int i = 0; i < numItems; i++)
    if( !rsGreaterAbs(items[i].value - item, TItem(tol)) )
      return i;
  return -1;
}
// try this with tol == 0 - done - seems to work

template<class TItem, class TTol>
int rsFindOccurrence(const std::vector<rsOccurrence<TItem>>& items, const TItem& item, TTol tol)
{
  if(items.size() == 0)
    return -1;
  return rsFindOccurrence(&items[0], (int) items.size(), item, tol);
}

template<class TItem, class TTol>
std::vector<rsOccurrence<TItem>> collectOccurrences(const std::vector<TItem>& items, TTol tol)
{
  std::vector<rsOccurrence<TItem>> occurrences;
  for(size_t i = 0; i < items.size(); i++) {
    int j = rsFindOccurrence(occurrences, items[i], tol);
    if(j != -1)
      occurrences[j].multiplicity++;
    else
      occurrences.push_back(rsOccurrence<TItem>(items[i], 1)); }
  return occurrences;
}



// this function should probably take a complex polynomial as input, so we can use it with complex
// matrices too

template<class T>
std::vector<rsOccurrence<std::complex<T>>>
getRootsWithMultiplicities(const rsPolynomial<std::complex<T>> p, T tol)
{
  return collectOccurrences(getPolynomialRoots(p), tol);
}

template<class T>
bool isValidEigenSpaceSet(const std::vector<rsEigenSpace<T>>& ess)
{
  int algMulSum = 0; // sum of algebraic multiplicities of eigenvalues

  // for each eigenvalue, we must have: 1 <= geoMul <= algMul
  for(size_t i = 0; i < ess.size(); i++) {
    int algMul = ess[i].getAlgebraicMultiplicity();
    int geoMul = ess[i].getGeometricMultiplicity();
    algMulSum += algMul;
    if(geoMul < 1 || geoMul > algMul)
      return false; }

  // sum of algebraic multiplicities must equal dimensionality of the embedding space:
  for(size_t i = 0; i < ess.size(); i++)
    if(ess[i].eigenSpace.getNumRows() != algMulSum)
      return false;

  return true;
}
// are there any more conditions for a sane result?
// sum of eigenvalues equals the trace of the matrix A and the product of the eigenvalues equals
// the determinant (Karpf. pg 412). - each eigenvalue appears as often in the sum or product as
// its multiplicity says - to check that, we would have to take the matrix A (and/or its
// row-echelon from) as additional inputs. however - the row-echelon form may have a determinant
// with the sign swapped (if the number of swaps in gaussian elimination was odd) - maybe we should
// keep track of that in the elimination process...

/** Returns the eigenspaces of the matrix A as an array of rsEigenSpace objects. Each such object
contains the eigenspace represented as matrix whose columns form a basis of the eigenspace. It also
contains the associated eigenvalue together with its algebraic multiplicity. The geometric
multiplicity is given by the number of columns of the matrix of basis-vectors. */
template<class T>
std::vector<rsEigenSpace<T>> getEigenSpaces(rsMatrix<std::complex<T>> A, T tol)
{
  using Complex = std::complex<T>;
  rsPolynomial<Complex> p = getCharacteristicPolynomial(A, Complex(tol));
  std::vector<rsOccurrence<Complex>> eigenValues = getRootsWithMultiplicities(p, tol);
  int numRoots = (int) eigenValues.size();
  std::vector<rsEigenSpace<T>> eigenSpaces(numRoots);
  rsMatrix<Complex> Ai = A;
  for(int i = 0; i < numRoots; i++) {
    eigenSpaces[i].eigenValue = eigenValues[i].value;
    eigenSpaces[i].algebraicMultiplicity = eigenValues[i].multiplicity;
    for(int j = 0; j < rsMin(A.getNumRows(), A.getNumColumns()); j++)
      Ai(j, j) = A(j, j) - eigenValues[i].value;                // Ai = A - eigenValue[i] * I
    eigenSpaces[i].eigenSpace = getNullSpace(Ai, Complex(tol)); // complexifying tol is unelegant!
  }
  rsAssert(isValidEigenSpaceSet(eigenSpaces));  // sanity check
  return eigenSpaces;
}
// this is the textbook method - using the characteristic polynomial and finding its roots is not
// the right way to do it numerically - this function is for proof of concept and should not be
// used in production

// convenience function for matrices of real numbers:
template<class T>
std::vector<rsEigenSpace<T>> getEigenSpaces(rsMatrix<T> A, T tol)
{
  return getEigenSpaces(complexify(A), tol);
}

// Here: https://en.wikipedia.org/wiki/Eigenvalue_algorithm it says:
// "...the columns of the matrix prod_{i != j} (A - x_i * I)^m_i must be either 0 or generalized
// eigenvectors of the eigenvalue x_j..." (i've adapted the notation a bit). The index i runs over
// the distinct eigenvalues. So, as an alternative to finding nullspaces, we may do repeated matrix
// multiplication? It also says "generalized eigenvectors" - so maybe this algo works even if the
// matrix is not diagonalizable? ...or will we get the zero-columns case in this case?

// todo:
// -write a function to compute generalized eigenspaces and a Jordan normal form - see:
// https://en.wikipedia.org/wiki/Generalized_eigenvector


template<class T>
rsMatrix<T> getOrthogonalComplement(rsMatrix<T> A, T tol)
{
  return getNullSpace(A.getTranspose(), tol);  // verify if that's correct
}
// needs test

// todo: getProjection(rsMatrix A, Vector v) - should project the vector v onto the basis spanned
// by the columns of A. this can be done by forming a linear combination of the columns of A with
// coeffs given by the scalar products of the repsective column with the target vector v -  i think

/** Computes scalar product of columns i and j of matrix A. */
template<class T>
T getColumnScalarProduct(const rsMatrix<T>& A, int i, int j)
{
  T sum = T(0);
  for(int k = 0; k < A.getNumRows(); k++)
    sum += A(k, i) * A(k, j);
  return sum;
}
// todo: maybe make a complex version - it should conjugate one of the inputs in the product
// (which?)

template<class T>
T getColumnSquaredNorm(const rsMatrix<T>& A, int j)
{
  return getColumnScalarProduct(A, j, j);
}

template<class T>
T getColumnNorm(const rsMatrix<T>& A, int j)
{
  return sqrt(getColumnSquaredNorm(A, j));
}

/** Normalizes the Euclidean of column j in matrix A (seen as column-vector) to unity. */
template<class T>
void normalizeColumn(rsMatrix<T>& A, int j)
{
  T norm = getColumnNorm(A, j);
  A.scaleColumn(j, T(1)/norm);
}

template<class T>
void normalizeColumns(rsMatrix<T>& A)
{
  for(int j = 0; j < A.getNumColumns(); j++)
    normalizeColumn(A, j);
}


/** Copies the j-th column of matrix A into vector v. */
template<class T>
void copyColumn(const rsMatrix<T>& A, int j, std::vector<T>& v)
{
  v.resize(A.getNumRows());
  for(int i = 0; i < A.getNumRows(); i++)
    v[i] = A(i, j);
}

template<class T>
void pasteColumn(rsMatrix<T>& A, int j, const std::vector<T>& v)
{
  rsAssert(A.getNumRows() == (int) v.size());
  for(int i = 0; i < A.getNumRows(); i++)
    A(i, j) = v[i];
}

/** Naive implementation of Gram-Schmidt orthonormalization of the columns of A. This algorithm is
numerically very bad and not recommended for use in production. */
template<class T>
void orthonormalizeColumns1(rsMatrix<T>& A)
{
  normalizeColumn(A, 0);
  std::vector<T> tmp(A.getNumRows());
  for(int i = 1; i < A.getNumColumns(); i++) {
    copyColumn(A, i, tmp);
    for(int j = 0; j < i; j++) {
      T w = getColumnScalarProduct(A, i, j);
      for(int k = 0; k < A.getNumRows(); k++)
        tmp[k] -= w * A(k, j); }
    pasteColumn(A, i, tmp);
    normalizeColumn(A, i); }
}


/** Implements orthonormalization of the columns of A via Gaussian eliminations. See:
https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process#Via_Gaussian_elimination */
template<class T>
void orthonormalizeColumns2(rsMatrix<T>& A, T tol)
{
  A = A.getTranspose();
  rsMatrix<T> A2 = A * A.getTranspose();
  rowEchelon(A2, A, tol);
  A = A.getTranspose();
  normalizeColumns(A);
}
// try to get rid of some transpositions
// ...how does this algorithm fare numerically? Gaussian elimination itself is supposed to be
// numerically well behvaed, right?


// todo: implement numerically stabilized version, version based on Gaussian elimination (done) and
// based on Householder transformations - see:
// https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process#Numerical_stability
// https://en.wikipedia.org/wiki/Householder_transformation


template<class T>
bool areColumnsNormalized(rsMatrix<T>& A, T tol)
{
  for(int i = 0; i < A.getNumColumns(); i++)
    if(rsAbs(getColumnNorm(A, i) - T(1)) > tol)
      return false;
  return true;
}


/** Returns true, iff all the columns of the matrix A (seen as vectors) are mutually orthogonal.
This means, the scalar product of any pair of distinct columns must be zero. Note that we do
not require, that the scalar product of a column with itself is unity. ...i think, this would be
the definition of an orthogonal matrix - math terminology seems a bit inconsistent here: a matrix
being orthogonal seems a stronger requirement than its columns beings mutually orthogonal - the
columns must be orthoNORMal for the matrix being orthoGONal....-> look up
See: https://en.wikipedia.org/wiki/Orthogonality  */
template<class T>
bool areColumnsOrthogonal(rsMatrix<T>& A, T tol)
{
  for(int i = 0; i < A.getNumColumns(); i++) {
    for(int j = i+1; j < A.getNumColumns(); j++) {
      T sp = getColumnScalarProduct(A, i, j);
      if(rsAbs(sp) > tol)
        return false; }}
  return true;
}

template<class T>
bool areColumnsOrthonormal(rsMatrix<T>& A, T tol)
{
  return areColumnsOrthogonal(A, tol) && areColumnsNormalized(A, tol);
}
// This means that A^T * A = I where I is the identity matrix

template<class T>
bool isOrthogonal(rsMatrix<T>& A, T tol)
{
  return areColumnsOrthonormal(A, tol);
}
// It seems, math terminology is inconsistent here: for a matrix to count as orthoGONal, its
// columns must be orthoNORMal. the definition of matrix orthogonality is that the linear map
// doesn't change the scalar product, i.e <x, y> = <A*x, A*y> which implies A^T * A = I
// if a matrix A is orthogonal, det(A) = +-1 (what if A is complex?)


template<class T>
rsMatrix<T> getHouseholderReflection(rsMatrix<T>& a)
{
  rsAssert(a.isColumnVector());
  T w = T(2) / getColumnSquaredNorm(a, 0);
  return rsMatrix<T>::identity(a.getNumRows()) - w * a * a.getTranspose();
}
// needs test, maybe optimize
// H_a = I - (2/(a^T * a)) * a * a^T - reflection along vector a (i think, this means reflection
// about a plane whose normal is a - figure out - see Karpf, pg 156)
// if possible, make a function applyHousholderReflection(Matrix& A, const Matrix& a) to apply the
// reflection in place - maybe it should apply the reflection simultaneously to two matrices - in
// algorithms such as the QR-algo, we need to apply it to matrix R but also keep track of what we
// have done in matrix Q - matrix Q accumulates the steps taken

// make a function getGivensRotation
// Householder reflections and Givens rotations are also called elementary orthogonal
// transformations -  they may be usd to produce zeros somewhere (see section 1.3 of:
// https://www.researchgate.net/publication/277069471_Numerical_Linear_Algebra)


template<class T>
rsMatrix<T> getGivensRotation(int N, int i, int j, T c, T s)
{
  rsAssert(i >= 0 && i < N);
  rsAssert(j >= 0 && j < N);
  rsMatrix<T> G(N, N);
  G.setToIdentity();
  G(i, i) =  c;
  G(j, j) =  c;
  G(i, j) =  s;
  G(j, i) = -s;
  return G;
}
// needs test
// c = cos(a), s = sin(a) - but we don't take the angle a as paremeter because the c,s values may
// actually be computed by different formulas, for example using:
// d = sqrt(xi^2 + xj^2), c = xi/d, s = xj/d
// will transform a vector x = (...,xi,...,xj,...) with nonzero xi,xj into y = G*x where
// yk = xk for k != i,j and yi = c*xi + s*xj, yj = c*xj - s*xi, see pg. 7 of:
// https://www.researchgate.net/publication/277069471_Numerical_Linear_Algebra
// so we leave the computation of c,s to client code. what if xi == 0 or xj == 0 - will the formula
// still work - just may be not as useful?
// needs test - also, we need a function to apply a Givens rotation in place - in practice, it's
// silly to actually create the whole matrix - this is just for prototyping



template<class T>
void decomposeQR(const rsMatrix<T>& A, rsMatrix<T>& Q, rsMatrix<T>& R)
{
  int n = A.getNumRows();
  int r = A.getNumColumns();
  rsAssert(n >= r);    // Karpf. pg.181 - do we need this?
  Q.setShape(n, n);
  Q.setToIdentity();
  R = A;
  rsMatrix<T> s(n, 1), a(n, 1), H(n, n);
  for(int j = 0; j < r; j++)
  {
    // copy column j into s, but with zeros up to j-1:
    int k;
    for(k = 0; k < j; k++) s(k, 0) = T(0);
    for(k = j; k < n; k++) s(k, 0) = R(k, j);

    // compute weight alpha:
    T sNorm = getColumnNorm(s, 0);
    T alpha = s(j, 0) < T(0) ? sNorm : -sNorm;

    // compute reflection vector a = s - alpha * e_j:
    a = s; a(j, 0) -= alpha;  // can be streamlined - copy can be avoided

    // compute Householder reflection matrix H_a about vector a:
    H = getHouseholderReflection(a);

    // update Q and R:
    Q = Q*H;
    R = H*R;
  }
}
// QR-decomposition based on Householder reflections (see Karpf. pg. 184)
// prototype - can be streamlined/optimized - lots of copying and allocation can be avoided
// the Householder reflection matrices may not have to be constructed explicitly - maybe they can
// be applied directly pre/postMultiplyByHoudeholderReflection(A, a)
// todo: implement recipies pg. 187,188 - using the QR decomposition to solve a linear system of
// equations and an overdetermined system (needs reduced QR decomposition)
// how would it be done when A is complex? will the transpositions in the Householder reflection be
// replaced with conjugate transposes?
// maybe make another version based on Givens rotations


/** Pastes submatrix S into matrix A starting at row-index iStart and column-index jStart. */
template<class T>
void pasteSubMatrix(rsMatrixView<T>& A, const rsMatrixView<T>& S, int iStart, int jStart)
{
  rsAssert(A.getNumRows()    >= S.getNumRows()    + iStart);
  rsAssert(A.getNumColumns() >= S.getNumColumns() + jStart);
  for(int i = 0; i < S.getNumRows(); i++)
    for(int j = 0; j < S.getNumColumns(); j++)
      A(iStart + i, jStart + j) = S(i, j);
}
// make member function of rsMatrixView so it may be called like A.pasteSubMatrix(S,..)


/** Computes the singular value decomposition of a real-valued MxN matrix A. This expresses A as a
product of 3 matrices: A = U * S * V^T where U is an MxM orthogonal matrix, V is an NxN orthogonal
matrix and S is and MxN diaogonal matrix (filled up with zeros at the bottom or right if M != N).
If M = N = 2, such a decomposition can be visualized as breaking up any linear map into a
(possibly improper) rotation followed by a scaling along the coordinate axes followed by another
(possibly improper) rotation - where "improper" means that there could be a reflection involved as
well. N is the dimensionality of the input space, M is the dimensionality of the output space and
the rank R of A is the dimensionality of the image of A within the M-dimensional output space.


...write more about interpretation and applications */
template<class R> // R is a real-number datatype
void decomposeRealUSV(const rsMatrix<R>& A, rsMatrix<R>& U, rsMatrix<R>& S, rsMatrix<R>& V, R tol)
{
  // Algorithm:
  // -Compute the eigenvalues of A^T * A and sort them in descencding order. Because A^T *A is a
  //  symmetric matrix, all of these eigenvalues are real and nonnegative.
  //   -The square-roots of these eigenvalues are called the singular values of A.
  //   -The number R of nonzero singular values equals the rank of A. We have R <= min(M,N). (verify!)
  // -Construct an orthonormal basis for N-dimensional space from the eigenvectors of A^T * A.
  //  These basis vectors v_i are written as columns into the matrix V. This is always possible
  //  because symmetric matrices are always diagonalizable (verify!). We also have that
  //  eigenvectors to different eigenvalues are orthogonal. If an eigenvalue has an algebraic
  //  multiplicity > 1, it's geometric multiplicity will be the same (really? verify!) - so we get
  //  indeed a full set of basis vectors
  // -Construct the matrix S by writing the singular values sigma_i on its diagonal and zeros
  //  everywhere else: S(i,i) = sigma_i
  // -Construct the matrix U by taking a column u_i as u_i = (1/sigma_i) * A * v_i for
  //  i = 0,...,R-1. If R < M, use for the remaining columns a basis for the orthogonal complement
  //  of the vectors u_i constructed so far.
  //  -The u_i are also eigenvectors of A * A^T (verify) and A^T * A and A * A^T have the same
  //   eigenvalues (i think, more generally A*B has the same eigenvalues as B*A, if the dimensions
  //   are such that the products make sense - verify!)


  // A is an m-by-n matrix:
  int m = A.getNumRows();      // m is dimensionality of output space
  int n = A.getNumColumns();   // n is dimensionality of input space

  // Find eigenvalues of A^T * A (they are all non-negative), sort them in descending order and
  // figure out r, the number of nonzero eigenvalues (which is also the rank of A):
  rsMatrix<R>    ATA    = A.getTranspose() * A;          // A^T * A is an n-by-n matrix
  std::vector<R> lambda = getEigenvaluesReal(ATA, tol);  // eigenvalues of A^T * A, lambda_i >= 0
  rsAssert((int) lambda.size() == n);                    // sanity check for debug
  rsHeapSort(&lambda[0], n, rsGreater);                  // sort descending
  int r = 0;                                             // rank of A, dimensionality of image under A
  while(r < n && lambda[r] > tol)                        // figure out rank r of A
    r++;

  // For each eigenvalue lambda_i, compute the eigenspace. From those eigenspaces, construct the
  // matrix V: (v_1,...,v_n) such that (A^T * A) * v_i = lambda_i * v_i. The v_i become the columns
  // of V:
  int i = 0, j, k;
  V.setShape(n, n);
  rsMatrix<R> v_i;
  rsMatrix<R> tmp = ATA;
  while(i < n) {
    for(k = 0; k < n; k++)
      tmp(k, k) = ATA(k, k) - lambda[i];  // form matrix A^T * A - lambda_i * I
    v_i = getNullSpace(tmp, tol);         // its nullspace is the eigenspace to lambda_i
    orthonormalizeColumns1(v_i);          // make basis an ONB (use better algo later)
    int d_i = v_i.getNumColumns();        // dimensionality of i-th eigenspace
    pasteSubMatrix(V, v_i, 0, i);         // paste d_i columns into V
    i += d_i;                             // we filled d_i columns of V in this iteration
  }
  //orthonormalizeColumns1(V);                // ...we need this! (use better algo later!)
  // maye it's enough to call orthonormalizeColumns1(vi) inside the loop? this should be sufficient
  // iff the eigenspaces belonging to different eigenvalues are all mutually orthogonal - is this
  // the case? ...at least, with our tests so far, it seems to be - make more tests and research -
  // doing it inside the loop is preferable (less work, less error-accumulation) - but we need to
  // ensure that it is actually valid to do it like this! if not, call orthonormalizeColumns1(V)
  // after the loop and remove orthonormalizeColumns1(v_i) from the loop
  // yes - this holds indeed true for symmetric matrices:
  // https://www.mathelounge.de/607249/eigenvektoren-verschiedenen-eigenwerten-sind-orthogonal
  // https://resources.mpi-inf.mpg.de/departments/d1/teaching/ss10/MFI2/kap46.pdf
  // https://matheplanet.com/default3.html?call=viewtopic.php?topic=129218&ref=https%3A%2F%2Fwww.google.de%2F
  // and A^T * A is symmetric, so we should be fine

  // Construct the diagonal matrix S from the singular values sigma_i, which are the square-roots
  // of the eigenvalues lambda_i:
  S.setShape(m, n);
  S.setToZero();
  //for(i = 0; i < rsMin(m, n); i++)
  for(i = 0; i < rsMin(r, n); i++)    // for i >= r, sqrt(lambda[i]) == 0 - no need to compute it
    S(i, i) = sqrt(lambda[i]);
  // todo: avoid taking square-roots of zero - just loop up to i < rsMin(m, r)


  // Construct matrix U = (u_1,...,u_m) where u_1,..,u_r are computed from the nonzero singular
  // values sigma_i and corresponding basis-vectors v_i as: u_i = (1/sigma_i) * A * v_i and the
  // remaining u_{r+1},...,u_m (if any) are a basis of the orthogonal complement of u_1,..,u_r:
  U.setShape(m, m);
  U.setToZero();
  for(i = 0; i < r; i++) {             // i: col-index into U
    for(j = 0; j < m; j++) {           // j: row-index into U
      for(k = 0; k < n; k++)           // k: col-index into A, row-index into V (...is n correct?)
        U(j, i) += A(j, k) * V(k, i);  // check indices
      U(j, i) /= S(i, i); }}           // S(i, i) is sigma_i
  if(r < m) {
    rsMatrix<R> Uo = getOrthogonalComplement(U, tol); // Uo is not necessarily orthonormal
    orthonormalizeColumns1(Uo);                       // naive Gram-Schmidt - use better algo later
    pasteSubMatrix(U, Uo, 0, r); }


  // in this video here, he says, that U can also be computed as the matrix of eigenvectors of
  // A * A^T:
  // https://www.youtube.com/watch?v=mBcLRGuAFUk
  // maybe we should do that? and/or maybe compare to the results of the algo above? do the
  // eigenvalues of this matrix also have any meaning? oh - he also says, it has the same
  // eigenvalues because A*B has the same eigenvalues as B*A - maybe to compute the eigenvalues, we
  // should select the smaller of the two matrices A^T * A, A * A^T -> less work, less error
  // -> try it


  // what about computing the whole U matrix at once like:
  // U = A * V;  // is that correct? the book says, we only obtain r vectors u_i this way
  // for(i = 0; i < n; i++)  // verify - book says, we only obtain r vectors u_i this way
  //   for(j = 0; j < m; j++)
  //     U(j, i) /= S(i, i);
  // will we get an m-by-m matrix this way? the book says, we only obtain r vectors u_i by
  // multiplying A with the columns v_i (and dividing by the i-th singular value)


  // maybe we should use:
  //std::vector<rsEigenSpace<R>> es = getEigenSpaces(ATA, tol);
  // ...and forget the stuff above...at least in a prototype - sorting the eigenspace involves more
  // data moving that just sorting the eigenvalues ...but we should really use an array of
  // rsOccurence - we don't really want to compute the same eigenspace twice - but maybe we don't
  // have to - we may just skip an eigenvalue, if it already occurred before - we may actually also
  // detect this from the dimensionality of the eigenspace - it it's d, we increment our
  // array-index into ev by d
  // we may factor out a function getEigenSpace(const rsMatrix<T>& A, T ev)
}
// singular value decomposition (see Karpf. pg 447)
// if A is real, A^T * A is symmetric and this in turn implies that all eigenvalues are real (and
// nonnegative) and A^T * A is diagonalizable - so we don't need to worry about having to consider
// complex eigenvalues and/or defective eigenspaces (wher the geometric multiplicity is less than
// the algebraic)
// Karpf: pg 448: because A^T * A is positive semidefinite, it's eigenvalues are >= 0

// see also:
// https://blog.statsbot.co/singular-value-decomposition-tutorial-52c695315254
// https://machinelearningmastery.com/singular-value-decomposition-for-machine-learning/
// https://towardsdatascience.com/singular-value-decomposition-example-in-python-dab2507d85a0

// https://math.stackexchange.com/questions/158219/is-a-matrix-multiplied-with-its-transpose-something-special
// https://en.wikipedia.org/wiki/Spectral_theorem
// how does it generalize to the complex case? would we form A^H * A isntead of A^T * A? (A^H means Hermitian
// transpose aka conjugate transpose)
// https://en.wikipedia.org/wiki/Hermitian_matrix

// implement recipies: Karpf., pg.138,140,153,154,159(done),166(done?),172,174,176,184,187,188,
// 443
// formulas: 156

// Constraints for eigenvalues:
// Gerschgorin circles:
// -for an NxN matrix A over the complex numbers, all eigenvalues are within the union of the
//  circles with centers given by A(i, i) and radii given by sum_j(abs(A(i,j))) where j runs from
//  0 to N-1 but j=i is left out in the summation (Karpf. pg. 419)
//  -not every one of these circles must contain an eigenvalue - but it must contain one if it's
//   disjoint from all other circles
//  -if none of the circles contains 0, A is invertible, because none of its eigenvalues can be
//   zero, so the determinant is nonzero (because it's the product of the eigenvalues)
// -make a function vector<rsCircle<T>> getGerschgorinCircles(rsMatrix A) and use it to plot them
// Also:
// -for every matrix-norm that is induced by a vector-norm, we must have that the absolute value of
//  each eigenvalue is <= norm(A), so we have
//  |ev| <= max_i(sum_j(abs(A(i,j)))) and |ev| <= max_j(sum_i(abs(A(i,j)))) - Karpf. pg. 482
//  where the first inequality comes from the L^inf norm and the 2nd from the L^1 norm (verify!),
//  so we may use the smaller of these two norms as our limit (what about the L^2 norm? is it
//  always in between these two or can it be used to further constrain the eigenvalues? what about
//  yet other norms?)


// make a class rsSubSpace that defines arithmetic operations:
// -subspaces of a R^M are represented by MxN matrices whose columns define a basis of R^M
// -equals(A, B): spanSameSpace(A, B), A,B must have the same embedding space, i.e. their number
//  of rows must match
// -add(A, B): set union: throw A,B together, then remove linearly dependent vectors via
//  transpose -> row-echelon -> remove bottom zero lines -> transpose
// -mul(A, B): set intersection: use those vectors of A which are in the span of B
// -negate(A): orthogonal complement of A
// -sub(A,B): set difference remove those elements from A which are in the span of B..or no - that
//  doesn't seem right, maybe just use negate and add

// -let V be a vector space, W a subspace and W^C it's orthogonal complement
// -we want the set operations union and intersection and maybe difference
// -union should be like addition, intersection like multiplication
// -the full space V = W + W^C is the neutral element with respect to intersection/multiplication,
//  W^C is the multiplicative inverse to W (right?)
// -the zero vector 0 is the neutral element with respect to union/addition
// -we can define subtraction W-U as keeping only those vectors in W which are not also in U
//  ...is this any useful?
// -there are no inverse elements with respect to addition - subtraction is a fundamentally
//  different operation (it can not be expressed as union with some sort of inverse element)
// -what sort of algebraic structure do we get from this? is this a topological space?
// -subspace intersection should probably work as follows: project the basis vectors
//  of A onto those of B (by taking the scalar product and using the result as coefficient?)


// Compute solution sets of inhomogeneous systems of linear equations - for example, consider the
// plane in R^3 defined by: 2*x + 3*y + 5*z = 7. This plane can be represented by the linear system
//   2 3 5 | 7
//   0 0 0 | 0
//   0 0 0 | 0
// we may mangle the system a bit by applying row transformations - it will still always define the
// same plane. To compute the solution set, we need to compute one particluar solution x_0 of the
// inhomogenous system and the nullspace of the matrix. The full solution is then given by
// { x € R^3 : x = x_0 + a * v_1 + b * v_2 } where a,b are scalars and v_1,v_2 are the columns of
// nullspace matrix (it should come out as a two column matrix). ...i think

// todo: can we also compute a basis for the image in a similar way?
// If N is numRows and K is the rank, in order to find the nullspace, we have solve a KxK system
// fo (N-K)x(N-K) different right-hand sides. This gives

// here:
// https://www.mathwizurd.com/linalg/2018/12/10/orthogonal-complement
// https://textbooks.math.gatech.edu/ila/orthogonal-complements.html
// are some relations bewtween row-spaces, column-spaces and nullspaces. Let's denote row-,
// column- and nullspace as rsp,csp,nsp and the orthogonal complement as comp - then we have
// comp(rsp(A)) == nsp(A), csp(A) = rsp(A^T) = comp(nsp(A^T)), rsp(A) = comp(nsp(A)),
// comp(csp(A)) = nsp(A) - do we also have comp(A) = nsp(A^T)?
// checking equality for spaces means to check if the bases span the same space - so we could
// write checks like r &= spanSameSpace(orthogonalComplement(rowSpace(A)), nullSpace(A)) etc.


// https://www.wikihow.com/Find-the-Null-Space-of-a-Matrix

// https://en.wikipedia.org/wiki/Row_and_column_spaces
// says: the null space of A is the orthogonal complement to the row space - so maybe it's easier
// to compute the row-space and from that its orthogonal complement?
// https://en.wikipedia.org/wiki/Orthogonal_complement
// https://www.mathwizurd.com/linalg/2018/12/10/orthogonal-complement

// this has a nice worked through example:
// https://textbooks.math.gatech.edu/ila/orthogonal-complements.html

// https://www.youtube.com/watch?v=Kpc5ELrOt5E

//=================================================================================================

// Newton iteration with numric derivatives - todo: make a version that takes a second function to
// compute the analytic derivative
// x is initial estimate, y is target value for y
template<class T>
T newton(const std::function<T(T)>& f, T x, T y = T(0))
{
  static const int maxNumIterations = 100;
  T tol = std::numeric_limits<T>::epsilon(); // maybe that's too strict - maybe 10*epsilon is better
  T h   = 1.e-8;  // make parameter - maybe take sqrt(epsilon)
  for(int i = 1; i <= maxNumIterations; i++) {
    T err = f(x) - y;                // current error
    if(rsAbs(err) <= tol) break;     // converged
    T fp = (f(x+h)-f(x-h)) / (2*h);  // estimate of f'(x) by central difference
    x -= err/fp;                     // Newton update step
  }
  rsError("rsRootFinder::newton failed to converge");
  return x;
}
// include a check, if the error decreased from one iteration to the next - if this happens, report
// maybe use a while(true) loop - break if converged or error has increased
// failure - it means that Newton iteration did not converge
// move to rsRootFinder
// maybe make a version that uses a one-sided etsimated for f'(x) - this avoids one function
// evaluation per iteration (-> 2 instead of 3) - maybe have 3 versions: newtonLeft, newtonRight,
// newtonCentral that uses left-sidedn, right-sided and central estimates - what about convergence?
// the centered difference is 2nd order accurate while the one-sided is only 1st order accurate -
// do we loose the quadratic convergence when making such estimates - with the one-sided, most
// probably yes - but what about the two-sided?
// see the "damped-newton" method in "Python Hacking for Math Junkies", pg 306 -maybe it can be
// further improved by also taking an interval-halving step in cases of slow convergence, i.e.
// when the error does decrease but not fast enough

//
// try to make a function f(x,y) that has exponential spirals as contour lines
// the parametric equation for the exponential spiral is:
//   f(t) = exp(a*t)*cos(t), g(t) = exp(a*t)*sin(t)
// define the function d(x,y) as the distance of the point (x,y) to the nearest point on the
// spiral:
//   d(x,y) = min_t [ (x-f(t))^2 + (y-g(t))^2 ]
// where min_t means, we need to find the minimum with respect to t of the term in the brackets,
// then plug that value of t into the parametric spiral equation and them compute the distance of
// the point x,y to the resulting point on the spiral. to find the minimum, we have to take the
// derivative of the term in the brackets with respect to t:
//  d/dt [ (x-exp(a*t)*cos(t))^2 + (y-exp(a*t)*sin(t)) ]
//   = 2*(a*cos(t)*e^(a*t) - e^(a*t)*sin(t))*(cos(t)*e^(a*t) - x) + 2*(a*e^(a*t)*sin(t) + cos(t)*e^(a*t))*(e^(a*t)*sin(t) - y)
//
// because the distances increase exponentially with the radius of (x,y), define weighted distance
//   D(x,y) = d(x,y) / exp( sqrt(x^2 + y^2) )
// the function R(x,y) = 1 / (1 + D^2) is a sort of ridge that has the shape of the exponential
// spiral, so we may use R(x,y) = c it as our implicit equation that should have exponential
// spirals as level lines
//   ...(check, if this weighting is good) - we wnat something that goes exponentially to infinity
//   as x^2+y^2 goes to zero and exponentially to zero as x^2 + y^2 goes to infinity maby
//   exp(1 / (x^2+y^2))
//
// sage:
// var("x y a t")
// f(t) = exp(a*t) * cos(t)
// g(t) = exp(a*t) * sin(t)
// d(t,x,y) = (x-f(t))^2 + (y-g(t))^2
// d_dt = diff(d(t,x,y), t)
// solve(d_dt == 0, t)
//
// leads to unwieldy expressions that are not explicitly solved for t (only for sin(t) and in the
// rhs t also appears in cos(t) terms - sooo, it seems we need to solve that equation numerically
// for t - Netwon iteration or something
//
// strategy to find the distance:
// -input: the point x,y
// -find values xl < x, xr >= x such that f(t) = exp(a*t)*cos(t) = xl or xr
// -find a vlaue for t such that exp(a*t)*cos(t) = x, or maybe two values - incre
//
// rsImageGenerator::spiralRidge does implement a different function - but it works for the
// intended purpose of drawing spiral ridges just as well (maybe even better) - but the ideas here
// may be applicable to other, similar problems


// or maybe use a simpler linear spiral:
//   f(t) = t*cos(t), g(t) = t*sin(t)

//=================================================================================================
// Differential Geometry


/** Numeric approximation of the first derivative of function f at the value x with approximation
step-size h. Uses a central difference and is 2nd order accurate in h. */
template<class T, class F>
T derivative(F f, T x, T h)
{
  return (f(x+h) - f(x-h)) / (T(2)*h);
}

/** Numeric approximation of the first derivative. 2nd order accurate in h. */
template<class T, class F>
T secondDerivative(F f, T x, T h)
{
  return (f(x-h) - T(2)*f(x) + f(x+h)) / (h*h);  // verify formula
}
// https://en.wikipedia.org/wiki/Finite_difference#Higher-order_differences

template<class T, class F>
T thirdDerivative(F f, T x, T h)
{
  return (-f(x-2*h) + 2*f(x-h) - 2*f(x+h) + 1*f(x+2*h)) / (2*h*h*h);
}
// coeffs found by:
// http://web.media.mit.edu/~crtaylor/calculator.html
// f_xxx = (-1*f[i-2]+2*f[i-1]+0*f[i+0]-2*f[i+1]+1*f[i+2])/(2*1.0*h**3)

// under construction:
/** Computes derivatives 0..3 with a 5-point stencil. */
template<class T, class F>
void derivativesUpTo3(F f, T x, T h, T* f0, T* f1, T* f2, T* f3)
{
  // evaluate function at stencil points:
  T fm2 = f(x-2*h);
  T fm1 = f(x - h);
  T fc  = f(x    );
  T fp1 = f(x + h);
  T fp2 = f(x+2*h);

  // form linear combinations to approximate derivatives:
  *f0 =                     fc;
  *f1 = ( fm2 -  8*fm1         +  8*fp1 - fp2) / (12*h);
  *f2 = (-fm2 + 16*fm1 - 30*fc + 16*fp1 - fp2) / (12*h*h);
  *f3 = (-fm2 +  2*fm1         -  2*fp1 + fp2) / (2*h*h*h);
}
// needs test
// stencil: -2,-1,0,1,2
// f_x = (1*f[i-2]-8*f[i-1]+0*f[i+0]+8*f[i+1]-1*f[i+2])/(12*1.0*h**1)
// f_xx = (-1*f[i-2]+16*f[i-1]-30*f[i+0]+16*f[i+1]-1*f[i+2])/(12*1.0*h**2)
// f_xxx = (-1*f[i-2]+2*f[i-1]+0*f[i+0]-2*f[i+1]+1*f[i+2])/(2*1.0*h**3)


// move these into rsNumericDifferentiator, make functions that compute all derivatives up to a
// given order - avoid re-evaluation of f(x+h) etc - evaluate them once and use different linear
// combinations to form approximations of the derivatives

// third derivative needs a 4-point-stencil at least

// maybe make a function that computes derivatives up to order 3 with a 5-point stencil using
// x-2*h,x-h,x,x+h,x+2*h - this should be good enough for differential geometry for graphical
// purposes




/** A class to represent n-dimensional parametric curves. */

template<class TScl, class TVec>  // scalar (for the parameter) and vector (for output) types
class rsParametricCurve
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the function that computes a position vector from a parameter value. You must pass a
  function object that takes a scalar of type TScl as input and returns a vector of type TVec. The
  function object must be assignable to a std::function...(lambda-functions, function-pointers or
  std::function objects will work, i think) */
  template<class F>
  void setPositionFunction(const F& newFunc) { f = newFunc; }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the position vector on the curve for the given value of the parameter t. The
  parameter may be interpreted as time. */
  TVec getPosition(TScl t) const { return f(t);  }

  /** Returns the velocity vector on the curve for the given value of the parameter t. The velocity
  is given by the derivative of the position vector with respect to the time parameter t. The
  velocity vector is tangent to the curve. This function uses a central difference approximation
  with stepsize h to numerically compute the derivative.  */
  TVec getVelocity(TScl t, TScl h) const { return NumDiff::derivative(f, t, h); }

  /** Returns the acceleration vector at the given value of parameter t. This is the second
  derivative of the position vector with respect to t. Iff the curve is parametrized by arc-length,
  i.e. the speed is always unity, the acceleration vector will be orthogonal to the velocity
  vector. That means, the acceleration will only change the direction but not the speed of the
  imagined moving point. But this really only holds form parametrizations by arc-length, aka
  natural parametrizations. */
  TVec getAcceleration(TScl t, TScl h) const { return NumDiff::secondDerivative(f, t, h); }



  TVec getJerk(TScl t, TScl h) const { return NumDiff::thirdDerivative(f, t, h); }
  // not yet tested
  // the 3rd derivative is used in some formulas for 3D curves - do we need any higher derivatives?
  // do they have names, too?
  // yes 4th derivative is calle jounce or snap:
  // https://en.wikipedia.org/wiki/Jounce/
  // jerk is also called jolt:
  // https://en.wikipedia.org/wiki/Jerk_(physics)
  // after snap comes crackle and pop
  // https://en.wikipedia.org/wiki/Crackle_(physics)
  // https://en.wikipedia.org/wiki/Pop_(physics)

  /** Joint computation of position, velocity and acceleration to save some function evaluations
  compared to evaluating them all separately (reduces the number of evaluations from 6 to 3). */
  void getPosVelAcc(TScl t, TScl h, TVec* p, TVec* v, TVec* a)
  { NumDiff::derivativesUpTo2(this->f, t, h, p, v, a); }
  //{ NumDiff::derivativesUpTo2(this->f, t, h, f0, f1, f2); }



  /** Returns an array with values of the arc-lengths corresponding to the parameter values given
  in t. */
  std::vector<TScl> getArcLengthFunction(const std::vector<TScl>& t)
  {
    // create and fill s-array with values of length-elements ds:
    size_t N = t.size();
    std::vector<TScl> s(N);
    s[0] = TScl(0);
    TVec v0 = getPosition(t[0]);     // tail of the current vector
    for(size_t n = 1; n < N; n++) {  // loop over length elements
      TVec v1 = getPosition(t[n]);   // tip of the current vector
      TVec dv = v1 - v0;             // connection between tip and tail
      s[n] = rsNorm(dv);             // distance ds between tip and tail (ds, integrand)
      v0 = v1; }                     // old tip becomes new tail

    // numerically integrate the ds values to obtain s itself:
    rsArrayTools::cumulativeSum(&s[0], &s[0], (int) N);
    // todo: use trapezoidal rule later - or maybe even better methods
    // don't we need the t-array here?

    return s;
  }
  // needs test
  // -maybe to improve accuracy, we could approximate the lenght element ds not by a straight line
  //  but by a cubic spline segment? for this, we would need (numeric approximations of) the
  //  derivative (dx/dt, dy/dt, ...) - then we could create Hermite spline segments or a natural
  //  cubic spline through the points
  // maybe have a function that computes the arc-length for a particular value of t - but document
  // that it should not be used in a loop to compute the arc-length function for multiple values of
  // ascending t-values (because that would result in a "Shlemiel-the-painter" algo)
  // this should use a generic integrate(Functor f, T a, T b, int numPoints) function

  /** Returns an array of values of t from t0 to t1, such that the arc-length between successive
  points of the curve for successive array-values of the parameter t is equal */
  std::vector<TScl> getArcLengthParametrization(TScl t0, TScl t1, int N)
  {
    using Vec = std::vector<TScl>;
    Vec  t = rsRangeLinear(t0, t1, N);      // equally spaced values of raw parameter t
    Vec  s = getArcLengthFunction(t);       // arc-length as function of raw parameter t
    TScl L = s[N-1];                        // total length of curve between t0 and t1
    Vec  r = rsRangeLinear(TScl(0), L, N);  // equally spaced values from 0..L
    Vec  t2(N);                             // new, modified array of time-stamps
    resampleNonUniformLinear(&s[0], &t[0], N, &r[0], &t2[0], N);
    return t2;
  }
  // maybe rename to getNaturalParameterMap or getArcLengthToParamMap, getNaturalReParametrization
  // getReparametrizationMap
  // needs test

  /** Returns the total length of the curve between parameter values t0 and t1 (using stepsize h
  for the numeric approximation of the derivative and N sample points for the numeric approximation
  of the integral). */
  TScl getTotalLength(TScl t0, TScl t1, TScl h, int N)
  {
    using Vec = std::vector<TScl>;
    Vec t = rsRangeLinear(t0, t1, N);
    Vec s = getArcLengthFunction(t);
    return s[N-1];
  }
  // needs test
  // todo: implement it in a way that doesn't need temporary arrays


  //std::vector<TScl> getArcLengthFunction(TScl t0, TScl t1, int numPoints);
  //std::vector<TScl> getArcLengthParameterMapping(TScl t0, TScl t1, int numPoints);
  // should use a generic rsNorm() function that returns the Euclidean norm of any vector -
  // implement it for rsVector2D, rsVector3D, std::vector

  // see plotParametricCurve, arcLengthFunction - but there, we use a Riemann sum - use trapezoidal
  // integration instead for better accurace - but the skeleton of the algo can be used



protected:

  using Func    = std::function<TVec(TScl)>;
  //using NumDiff = rsNumericDifferentiator<TScl, TVec, Func>; // old: <Tx,Ty,F>
  using NumDiff = rsNumericDifferentiator<TVec>;               // new: <Ty>

  Func f;

  // have optional members for (analytic computation of or specialized numerical algos) velocity,
  // acceleration
  // make it possible to pass functions for derivatives (velocity, acceleration) and fall back to
  // numeric derivatives only if nothing is passed

};



/** Partial specialization of rsParametricCurve for 2D curves. Still templatized on the scalar
type, i.e. the type to represent real numbers (the vector type is the chosen as just a vector of
the used scalar type).

References:
  (1) Edmund Weitz: Elementare Differentialgeometrie (nicht nur) für Informatiker

*/

template<class T>
class rsParametricCurve2D : public rsParametricCurve<T, rsVector2D<T>>
{

public:

  using Vec2 = rsVector2D<T>;


  //-----------------------------------------------------------------------------------------------
  // \name Setup


  //-----------------------------------------------------------------------------------------------
  // \name Local Features
  // for all these functions, you need to pass a parameter value t and a numeric approximation
  // stepsize h

  /** Computes the (signed) curvature. */
  T getCurvature(T t, T h)
  {
    Vec2 v = this->getVelocity(t, h);
    Vec2 a = this->getAcceleration(t, h);
    T d = rsDet(v, a);
    T s = v.getEuclideanNorm();  // speed
    return d / (s*s*s);          // (1), Eq. 4.2
  }
  // can be optimized: velocity and acceleration can be computed in a combined function - avoids
  // 2 of otherwise 5 evaluations of f
  // Note: for a curve of the form c(t) = (t,f(t)), i.e. a function, the curvature comes out as
  // f'' / (1 + f' * f')^(3/2)

  /** Computes the (signed) radius of curvature which is the reciprocal of the curvature itself. */
  T getRadiusOfCurvature(T t, T h) // is this the correct name?
  {
    return T(1) / getCurvature(t, h);
  }

  /** Returns the normal to the curve. The normal is at right angle to the velocity vector and
  points to the left, as seen from a point that traverses the curve. The normal vector returned
  here is not normalized to unit length (it will have a length given by the instantaneous
  speed). */
  Vec2 getNormal(T t, T h)
  {
    Vec2 v = this->getVelocity(t, h);
    return Vec2(-v.y, v.x);
  }
  // maybe getNormal should return a normalized vector already - it should be consistent with
  // getNormal of 3D curves and maybe even of surfaces (whether it's normalized or not and/or
  // if there's an extra function for normalized ones)..maybe have a boolean parameter "normalized"
  // that defaults to true

  /** Results the unit-length normal to the curve. */
  Vec2 getUnitNormal(T t, T h)
  {
    Vec2 n = getNormal(t, h);
    n.normalize();
    return n;
  }

  /** Computes the center of the osculating circle. This is the circle that best approximates the
  curve at any given point - it will be tangent to the curve at that point and have the same
  curvature. As a point traces out a curve, the center of the osculating circle traces out another
  curve which is called the evolute to the curve. */
  Vec2 getOsculatingCircleCenter(T t, T h)
  {
    Vec2 p = this->getPosition(t);
    Vec2 n = getUnitNormal(t, h);
    T r = getRadiusOfCurvature(t, h);
    return p + r*n;  // is this correct?
  }
  // can be optimized for production code
  // maybe rename to getEvolute

  // todo: implement the angle function

  //-----------------------------------------------------------------------------------------------
  // \name Global Features

  /*
  T getTotalCurvature(T t0, T t1, T h, int N)
  {

  }
  */


  // global features: total length (in an interval t = a..b) winding number around a given point
  // (applied to closed curves, needs the angle function), isClosed() or isInjectiveOn(a, b) -
  // create a list of points and compute the minimum of all the mutual distances - if it is below
  // a threshold, the points are considered the same and the curve is not injective...maybe call it
  // isSelfIntersecting, or maybe better: getNumSelfIntersections(a, b, threshold, numPoints) -
  // maybe for this, it makes sense to have a natural parameterization - otherwise, values for
  // t0,t1 close to each other may be below the threshold because the curve is very slow between
  // t0,t1
  // getGlobalCurvature (Totalkrümmung), Umlaufzahl (what's the english term?), winding number



  //-----------------------------------------------------------------------------------------------
  // \name Static functions (maybe get rid - or maybe not - maybe for optimized code, it's better
  // to directly operate on pointers rather than returning stack-created vector objects)

  /** Numerically computes the tangent vector to a 2D curve given by r(t) = (x(t),y(t)) at the given
  parameter t. */
  template<class F>
  static void velocity(F fx, F fy, T t, T* vx, T* vy, T h)
  {
    *vx = derivative(fx, t, h);
    *vy = derivative(fy, t, h);
  }

  template<class F>
  static void acceleration(F fx, F fy, T t, T* ax, T* ay, T h)
  {
    *ax = secondDerivative(fx, t, h);
    *ay = secondDerivative(fy, t, h);
  }


protected:

};


template<class T>
class rsParametricCurve3D : public rsParametricCurve<T, rsVector3D<T>>
{

public:

  using Vec3 = rsVector3D<T>;


  /*
  Vec3 getUnitTangent(T t, T h)
  {
    Vec3 v = getVelocity(t, h);
    return v / v.getEuclideanNorm();
  }
  */

  Vec3 getNormal(T t, T h, bool normalized = true)
  {
    Vec3 a = getAcceleration(t, h);
    if(normalized)
      return a / a.getEuclideanNorm();
    else
      return a;
  }
  // maybe rename to getNormal - does the non-normalized version have any significance? i mean, yes
  // it's the acceleration, but geometrically?

  Vec3 getBinormal(T t, T h)
  {
    Vec3 v = getVelocity(t, h);
    Vec3 n = getNormal(t, h, true);
    Vec3 b = cross(v, n);    // (1) pg. 96
    return b;                // the binormal vector is not normalized because v is not normalized!
  }


  void getFrenetFrame(T t, T h, Vec3* v, Vec3* n, Vec3* b)
  {
    getPosVelAcc(t, h, b, v, n);  // b is used as dummy here, t is done - it's the velocity
    n.normalize();                // n is normalized acceleration
    b = cross(v, n);              // binormal is cross-product between tangent and normal
  }
  // maybe it's useful, if this function would also return the position p
  // todo: figure out the usual conventions, which of these vectors are supposed to be normalized
  // and which don't
  // maybe we should not normalize anything - in an animation, it would perhaps make sense to have
  // longer vectors for stronger "effects" like acceleration



  T getCurvature(T t, T h)
  {
    Vec3 v = getVelocity(t, h);
    Vec3 a = getAcceleration(t, h);
    T d = rsNorm(cross(v, a));        // in 2D, it's det(v, a) - is this the same?
    T s = v.getEuclideanNorm();       // speed
    return d / (s*s*s);               // (1), Eq. 8.5
  }
  // needs test
  // check, how the norm of the cross-product relates to the determinant that is used in the
  // formula for 2D curves - if the norm of the cross-prodcut equals the 3x3 determinant (up to
  // sign), we may actually also define a signed curvature here, too via the determinant - would
  // that make any sense?

  T getTorsion(T t, T h)
  {
    // todo: compute them jointly:
    Vec3 v  = getVelocity(t, h);
    Vec3 a  = getAcceleration(t, h);
    Vec3 j  = getJerk(t, h);
    Vec3 va = cross(v, a);
    return dot(va, j) / va.getSquaredEuclideanNorm();  // (1), Eq. 8.5
  }
  // how does this relate to other formulas for the torsion like Eq 8.3
  // tau = -dotProduct(timeDerivative(binormal), normal)

  // todo: getTorsionAndCurvature - joint computation - avoid recomputations


  // todo getNormal, getBinormal, getFrenetFrame/FrenetTrihedron (begleitendes Dreibein)
  // updateFrenetSerret(curvature, torsion, *x, *y, *z), getCurvature(t, h), getTorsion(t, h),




  /*
  template<class T, class F>
  static void velocity(F fx, F fy, F fz, T t, T* vx, T* vy, T* vz)
  {
    T h = 1.e-8;  // use something better
    *vx = derivative(fx, t, h);
    *vy = derivative(fy, t, h);
    *vz = derivative(fz, t, h);
  }
  */

};


template<class T>
class rsParametricSurface
{

public:

  // todo: getPosition(T u, T v),
  // getCurveConstantU(T u), getCurveConstantV(T v)
  // returns std::vector<rsVector3D<T>> and/or rsParametricCurve3D objects

};





//=================================================================================================
// coordinate transformations:

// y = A*x
template<class T>
void applyMatrix4D(const T A[4][4], const T x[4], T y[4])
{
  rsAssert(x != y, "Cannot be used in place" );
  y[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2] + A[0][3]*x[3];
  y[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2] + A[1][3]*x[3];
  y[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2] + A[2][3]*x[3];
  y[3] = A[3][0]*x[0] + A[3][1]*x[1] + A[3][2]*x[2] + A[3][3]*x[3];
}

// C = A*B
template<class T>
void multiplyMatrices4D(const T A[4][4], const T B[4][4], T C[4][4])
{
  rsAssert(A != C && B != C, "Cannot be used in place" );
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      C[i][j] = T(0);
      for(int k = 0; k < 4; k++)
        C[i][j] += A[i][k] * B[k][j]; }}
}

template<class T>
void zeroMatrix4D(T A[4][4])
{
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      A[j][i] = T(0);
}

template<class T>
void identityMatrix4D(T A[4][4])
{
  zeroMatrix4D(A);
  for(int i = 0; i < 4; i++)
    A[i][i] = T(1);
}

// translation
template<class T>
void translationMatrix4D(T A[4][4], T dx, T dy, T dz)
{
  identityMatrix4D(A);
  A[0][3] = dx;
  A[1][3] = dy;
  A[2][3] = dz;
  // 4th column is used for translations (4th row is for perspective stuff)
}

// scaling
template<class T>
void scalingMatrix4D(T A[4][4], T sx, T sy, T sz)
{
  identityMatrix4D(A);
  A[0][0] = sx;
  A[1][1] = sy;
  A[2][2] = sz;
}

// rotation
template<class T>
void rotationMatrix4D(T A[4][4], T rx, T ry, T rz)
{
  identityMatrix4D(A);

  // sines/cosines:
  T sx = sin(rx); T cx = cos(rx);
  T sy = sin(ry); T cy = cos(ry);
  T sz = sin(rz); T cz = cos(rz);

  // rotation matrix coeffs:
  A[0][0] =  cz*cy;
  A[0][1] = -sz*cx + cz*sy*sx;
  A[0][2] =  sz*sx + cz*sy*cx;
  A[1][0] =  sz*cy;
  A[1][1] =  cz*cx + sz*sy*sx;
  A[1][2] = -cz*sx + sz*sy*cx;
  A[2][0] = -sy;
  A[2][1] =  cy*sx;
  A[2][2] =  cy*cx;

}
// code copied from rsRotationXYZ<T>::updateCoeffs - dirty! don't copy - re-use

template<class T>
void lookAtMatrix4D(T A[4][4],
  const rsVector3D<T>& eye, const rsVector3D<T>& center, const rsVector3D<T>& up, T zoom = T(1))
{
  rsVector3D<T> forward, side, up2;

  // forward: direction we look into
  forward = center - eye;
  forward.normalize();

  // Side = forward x up
  side = cross(forward, up);
  side.normalize();

  // Recompute up as: up = side x forward  (...why?)
  up2 = cross(side, forward);


  T m[4][4];
  identityMatrix4D(m);

  m[0][0] = zoom*side.x;
  m[0][1] = zoom*side.y;
  m[0][2] = zoom*side.z;

  m[1][0] = zoom*up2.x;
  m[1][1] = zoom*up2.y;
  m[1][2] = zoom*up2.z;

  m[2][0] = -zoom*forward.x;
  m[2][1] = -zoom*forward.y;
  m[2][2] = -zoom*forward.z;

  T t[4][4];
  translationMatrix4D(t, -eye.x, -eye.y, -eye.z);

  multiplyMatrices4D(m, t , A);  // A = m*t
  int dummy = 0;
}
// code adapted from vmath.h
// -compare results to what OpenGL's vmath::lookAt produces
// -why do we need the zoom - should the length of the forward vector, i.e. the distance between
//  eye and center determine the "closeness" - eye = (0,0,1) or eye = (0,0,2) produces the exact
//  same result. without zoom, we'll always have the same size of the frustum and some objects are
//  not visible unless we zoom out...





template<class T>
rsVector2D<T> project(const rsVector3D<T> v, const T A[4][4])
{
  T t1[4], t2[4];  // temporary homogeneous 4D vectors

  // copy 3D input vector into homogeneous 4D vector:
  t1[0] = v.x;   // x
  t1[1] = v.y;   // y
  t1[2] = v.z;   // z
  t1[3] = T(1);  // w = 1

  // apply homogeneous 4x4 matrix:
  applyMatrix4D(A, t1, t2);

  // extract x,y components and put into rsVector2D:
  return rsVector2D<T>(t2[0], t2[1]);
}
// this is not optimized code!


//=================================================================================================

/** A new experimental colorspace that is similar to HSL or HSV but with a twist that hopefully
overcomes the disadvantages of these... */

template<class T>
class rsColorBHS
{

public:

  /** Converts an (r,g,b) triple into a (B,H,S) triple. */
  void rgb2bhs(T r, T g, T b, T* B, T* H, T* S)
  {
    // compute brightness:
    *B = wr*r + wg*g + wb*b;

    // compute saturation:
    T max = rsMax(r, g, b);
    T min = rsMin(r, g, b);
    if(min == max) { // a gray-level (including black and white)...maybe have a tolerance
      *S = T(0);
      *H = T(0);
      return;   }
    *S = (max-min) / max;

    // compute hue:
    r -= min; g -= min; b -= min;
    if(     r == T(0)) *H = (g * T(1./3.) + b * T(2./3.)) / (g+b); // between green and blue
    else if(g == T(0)) *H = (b * T(2./3.) + r * T(3./3.)) / (b+r); // between blue and red
    else               *H = (r * T(0./3.) + g * T(1./3.)) / (r+g); // between red and gren, b == 0
  }

  /** Converts a (B,H,S) triple into an (r,g,b) triple. */
  void bhs2rgb(T B, T H, T S, T* r, T* g, T* b)
  {
    T k;
    if(H < T(1./3.)) {                                 // between red and green, so b = min
      if(H < T(1./6.)) {                               //   r = max
        k  = (3*H - 1);
        *r = -B*k/((k*wb + (6*H - 1)*wg)*S - k);
        *b = *r * (1-S);
        *g = (B - wb * *b - wr * *r) / wg;  }
      else {                                           //   g = max
        k  =  3*H;
        *g = -B*k/((k*wb + (6*H - 1)*wr)*S - k);
        *b = *g * (1-S);
        *r = (B - wb * *b - wg * *g) / wr; }}
    else if(H < T(2./3.))  {                           // between green and blue, so r = min
      if(T(2./3.) - H < H - T(1./3.)) {                //   b = max
        k  = (3*H - 1);
        *b = -B*k / ((3*(2*H - 1)*wg + k*wr)*S - k);
        *r = *b * (1-S);
        *g = (B - wb * *b - wr * *r) / wg;  }
      else  {                                          //   g = max
        k  = (3*H - 2);
        *g = -B*k/((3*(2*H - 1)*wb + k*wr)*S - k);
        *r = *g * (1-S);
        *b = (B - wg * *g - wr * *r) / wb; }}
    else {                                             // between blue and red, so g = min
      if(T(1)-H < H - T(2./3.)) {                      //   r = max
        k  = (3*H - 2);
        *r = -B*k/(((6*H - 5)*wb + k*wg)*S - k);
        *g = *r * (1-S);
        *b = (B - wg * *g - wr * *r) / wb; }
      else {                                           //   b = max
        k  = 3*(H - 1);
        *b = -B*k/((k*wg + (6*H - 5)*wr)*S - k);
        *g = *b * (1-S);
        *r = (B - wb * *b - wg * *g) / wr;  }}
  }


  /** Sets the weights for the 3 color channels in the brightness computation formula:
  brightness =  wr*r + wg*g + wb*b. The weights must sum up to unity. This weighting can be used to
  account for different sensitivity of the red, green and blue receptors such that a green with a
  given strength conributes more to the computed brightness than a blue of the same strength.  */
  void setWeights(T redWeight, T greenWeight, T blueWeight)
  {
    wr = redWeight;
    wg = greenWeight;
    wb = blueWeight;
    // todo: assert that the weights sum up to 1
  }

protected:

  // weights for the r,g,b components in the brightness formula (they should sum up to 1)
  T wr = T(1)/T(3);
  T wg = T(1)/T(3);
  T wb = T(1)/T(3);

};

// sage code for solving case blue..green closer to blue
// var("r,g,b,wr,wg,wb,B,H,S")
// eq1 = S == (b-r)/b
// eq2 = B == wr*r + wg*g + wb*b
// eq3 = H == ((g-r)*(1/3) + (b-r)*(2/3)) / ((g-r)+(b-r))
// solve([eq1,eq2,eq3],[r,g,b])
//
// r == (B*(3*H - 1)*S - B*(3*H - 1))/((3*(2*H - 1)*wg + (3*H - 1)*wr)*S - (3*H - 1)*wb - (3*H - 1)*wg - (3*H - 1)*wr),
// g == (3*B*(2*H - 1)*S - B*(3*H - 1))/((3*(2*H - 1)*wg + (3*H - 1)*wr)*S - (3*H - 1)*wb - (3*H - 1)*wg - (3*H - 1)*wr),
// b == -B*(3*H - 1)/((3*(2*H - 1)*wg + (3*H - 1)*wr)*S - (3*H - 1)*wb - (3*H - 1)*wg - (3*H - 1)*wr)
//
// -> use only formula for b, compute r,g via back-substitution when b is known
// for the branch blue..green, closer to blue, the 1st equation is replaced with S == (g-r)/g, the
// rest stays the same and we get:
//
// r == (B*(3*H - 2)*S - B*(3*H - 2))/((3*(2*H - 1)*wb + (3*H - 2)*wr)*S - (3*H - 2)*wb - (3*H - 2)*wg - (3*H - 2)*wr),
// g == -B*(3*H - 2)/((3*(2*H - 1)*wb + (3*H - 2)*wr)*S - (3*H - 2)*wb - (3*H - 2)*wg - (3*H - 2)*wr),
// b == (3*B*(2*H - 1)*S - B*(3*H - 2))/((3*(2*H - 1)*wb + (3*H - 2)*wr)*S - (3*H - 2)*wb - (3*H - 2)*wg - (3*H - 2)*wr)
//
// this time, we use the result for green, then compute red and blue by back-substitution
// for the 1st branch, we have S == (r-b)/r, H == ((r-b)*(0/3) + (g-b)*(1/3)) / ((r-b)+(g-b))
// and get:
// r == -B*(3*H - 1)/(((3*H - 1)*wb + (6*H - 1)*wg)*S - (3*H - 1)*wb - (3*H - 1)*wg - (3*H - 1)*wr),
// g == (B*(6*H - 1)*S - B*(3*H - 1))/(((3*H - 1)*wb + (6*H - 1)*wg)*S - (3*H - 1)*wb - (3*H - 1)*wg - (3*H - 1)*wr),
// b == (B*(3*H - 1)*S - B*(3*H - 1))/(((3*H - 1)*wb + (6*H - 1)*wg)*S - (3*H - 1)*wb - (3*H - 1)*wg - (3*H - 1)*wr)
// 2nd branch: use S == (g-b)/g, leave rest as is
// r == (B*(6*H - 1)*S - 3*B*H)/((3*H*wb + (6*H - 1)*wr)*S - 3*H*wb - 3*H*wg - 3*H*wr),
// g == -3*B*H/((3*H*wb + (6*H - 1)*wr)*S - 3*H*wb - 3*H*wg - 3*H*wr),
// b == 3*(B*H*S - B*H)/((3*H*wb + (6*H - 1)*wr)*S - 3*H*wb - 3*H*wg - 3*H*wr)
// between red and blue:
// H == ((r-g)*(3/3) + (b-g)*(2/3)) / ((r-g)+(b-g)),  S == (r-g)/r
// r == -B*(3*H - 2)/(((6*H - 5)*wb + (3*H - 2)*wg)*S - (3*H - 2)*wb - (3*H - 2)*wg - (3*H - 2)*wr),
// g == (B*(3*H - 2)*S - B*(3*H - 2))/(((6*H - 5)*wb + (3*H - 2)*wg)*S - (3*H - 2)*wb - (3*H - 2)*wg - (3*H - 2)*wr),
// b == (B*(6*H - 5)*S - B*(3*H - 2))/(((6*H - 5)*wb + (3*H - 2)*wg)*S - (3*H - 2)*wb - (3*H - 2)*wg - (3*H - 2)*wr)
// or  S = (b-g)/b:
// r == (B*(6*H - 5)*S - 3*B*(H - 1))/((3*(H - 1)*wg + (6*H - 5)*wr)*S - 3*(H - 1)*wb - 3*(H - 1)*wg - 3*(H - 1)*wr),
// g == 3*(B*(H - 1)*S - B*(H - 1))/((3*(H - 1)*wg + (6*H - 5)*wr)*S - 3*(H - 1)*wb - 3*(H - 1)*wg - 3*(H - 1)*wr),
// b == -3*B*(H - 1)/((3*(H - 1)*wg + (6*H - 5)*wr)*S - 3*(H - 1)*wb - 3*(H - 1)*wg - 3*(H - 1)*wr)

//=================================================================================================

template<class T, class F>
int minimizePartialParabolic(const F& f, T* v, int N, const T* h, T tol = 1.e-8)
{
  // Under construction

  // Algorithm:
  // Notation: x: current position vector, f(x): error funcion, f0n,f1n,f2n: value and 1st and 2nd
  // partial derivatives with respect to n-th coordinate
  // -at each step until convergence:
  //  -loop through the coordinates (n = 0..N-1):
  //   -compute value and 1st and 2nd partial derivatives f0n,f1n,f2n
  //   -if f2n > 0 (parabola along n-th coordinate has minimum):
  //    -jump into minimum of parabola along n-th coordinate
  //   -else:
  //    -jump an equal distance away from the maximum of the parabola (is this a good idea? or maybe
  //     we should use a smaller distance?)
  //  -compute function value at new location, if less than previous, accept step else reject and
  //   continue with next coordinate (or maybe try a half-step, then quarter, etc...before
  //   continuing)

  // Give the algorithm a name - maybe minimizePartialParabolic - maybe make a similar
  // minimizeGradientDescent function as baseline algo to compare against - however, it's probably
  // not advisable to use numeric gradients in gradient descent for efficiency reasons (each
  // gradient computation needs 3*N evaluations of f)
  // move into a class rsNumericMinimizer
  // it's similar to:
  // https://en.wikipedia.org/wiki/Coordinate_descent
  // just that we don't do line-search but instead use the numerical gradient info to jump into
  // the minimum of a parabola

  using NumDiff = rsNumericDifferentiator<T>;

  bool converged = false;
  int evals = 0;
  int iterations = 0;
  while(!converged)
  {

    // maybe factor out into a function minimizeStep1
    for(int n = 0; n < N; n++)
    {
      T x = v[n];      // variable of our 1D parabola
      T f0, f1, f2;
      NumDiff::partialDerivativesUpTo2(f, v, N, n, h[n], &f0, &f1, &f2); // 3 evals of f
      evals += 3;

      // compute parameters of parabola f(x) = a*x^2 + b*x + c
      T a, b, c;
      c = f0;             // not needed in the formulas
      a = f2/2;
      b = f1 - 2*a*x;     // 2*a = f2

      T fNew;
      T xEx = -b/(2*a);   // extremum (minimum or maximum) of the 1D parabola
      T dx  = xEx - x;    // update vector "delta-x"


      // this should probably be done in an acceptance loop - if the value of f has decreased,
      // accept the step, otherwise, reduce the stepsize by a factor (of 2?) and try again.
      if(f2 > 0)  // parabola has minimum
      {
        x += dx;          // todo: use x += step*dx;
        int dummy = 0;
      }
      else        // parabola has maximum
      {
        x -= dx;          // todo: use x -= step*dx;
        int dummy = 0;
        // this branch needs tests - to test it, we need a function with a local maximum somewhere
        // maybe use something like sin(x+y) or sin(x*y) - make contour-plots to get a feel for the
        // function
      }
      v[n] = x;
      fNew = f(v);
      evals++;



      if(rsAbs(f0-fNew) < tol)
        converged = true;
      else
        converged = false;
        // do we need this? can it happen, that in one inner iteration it gets set to true and in a
        // later one back to false? and if so - is this desirable?

      int dummy = 0;
    }

    iterations++;
    int dummy = 0;
  }

  return evals;  // return the number of evaluations of f
}
// maybe provide a means to keep track of the trajectory - it may make sense to plot that
// trajectory in a contour plot to see, what the algo does

// Other idea, for each step until convergence, do:
// -compute gradient
// -compute hessian * gradient
// -this determines a 1D parabola (right?) in the plane that intersects the function and contains
//  the gradient
// -jump into the minimum of the resulting 1D parabola - this is similar to jumping into the
//  minimum of the parabola above, but it uses the gradient direction instead of a coordinate
//  direction
// -see minimizeViaConjugateGradient: stepsize = - (d*g) / (dtH*d);
//  optimal stepsize - d is current direction, g is gradient, dtH is approximation of Hessian times
//  d (the t is for transposed, i think) - so in our case, d == g

template<class T, class F>
int minimizeGradientDescent(const F& f, T* v, int N, const T* h, T stepSize, T tol = 1.e-8)
{
  bool converged = false;
  int evals = 0;
  int iterations = 0;

  using AT     = rsArrayTools;
  using NumDif = rsNumericDifferentiator<T>;
  using Vec    = std::vector<T>;
  Vec g(N);  // gradient

  T fOld = f(v);
  evals++;
  T fNew = fOld;

  while(!converged)
  {
    NumDif::gradient(f, v, N, &g[0], h);       // g = numerical gradient
    evals += 2*N;                              // gradient estimation costs 2N evaluations of f
    AT::addWithWeight(v, N, &g[0], -stepSize); // v -= stepSize * g


    // maybe convergence can be assumed when all elements of the gradient have less absolute value
    // than their corresponding h-value? does that make sense?


    Vec dbg = toVector(v, N);

    fNew = f(v);
    evals++;

    // maybe, if the error has increased, reject the step, reduce the stepsize and try again

    if(rsAbs(fOld-fNew) < tol)
      converged = true;

    fOld = fNew;



    int dummy = 0;
  }
  return evals;
}
// https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
// https://docs.scipy.org/doc/scipy/reference/optimize.html


template<class T, class F>
int minimizeGradientAutoStep(const F& f, T* v, int N, const T* h, T tol = 1.e-8)
{
  // uses the gradient and automatic choice of the stepsize based on the product
  // Hessian times gradient...
  using AT     = rsArrayTools;
  using NumDif = rsNumericDifferentiator<T>;
  using Vec    = std::vector<T>;
  Vec g(N);     // gradient: g
  Vec Hg(N);    // Hessian times gradient: H*g
  Vec wrk(2*N); // workspace


  bool converged = false;
  int evals = 0;

  T fOld = f(v);
  evals++;
  T fNew = fOld;

  T k = pow(2, -14);  // make parameter


  while(!converged)
  {
    NumDif::gradient(f, v, N, &g[0], h);       // g = numerical gradient
    evals += 2*N;                              // gradient estimation costs 2N evaluations of f


    // compute optimal step:
    NumDif::hessianTimesVector(f, v, &g[0], N, &Hg[0], h, k, &wrk[0]);
    evals += 4*N;
    T step = AT::sumOfSquares(&g[0], N) / AT::sumOfProducts(&g[0], &Hg[0], N); // (g*g) / (g*H*g)
    // verify, if that formual is correct

    // do step:
    AT::addWithWeight(v, N, &g[0], -step); // or -step?



    //Vec dbg = toVector(v, N);


    fNew = f(v);
    evals++;

    // maybe, if the error has increased, reject the step, reduce the stepsize and try again

    if(rsAbs(fOld-fNew) < tol)
      converged = true;
    fOld = fNew;
  }


  return evals;
}

// does Newton steps into the minimum of a parabolic approximation of f where the gradient and
// Hessian are computed numerically, involving a lot of function evaluations - so this is probably
// not an efficient method
template<class T, class F>
int minimizeNewton(const F& f, T* x, int N, const T* h, T tol)
{
  using LinAlg = rsLinearAlgebraNew;
  using NumDif = rsNumericDifferentiator<T>;
  using Vec    = std::vector<T>;
  using AT     = rsArrayTools;
  rsMatrix<T> H(N, N); T* pH = H.getDataPointer(); // Hessian
  rsMatrix<T> g(N, 1); T* pg = g.getDataPointer(); // gradient
  rsMatrix<T> d(N, 1); T* pd = d.getDataPointer(); // update vector "delta-x"
  rsMatrix<T> X(N, 1); T *pX = X.getDataPointer(); // temporary vector for tentative new x
  T fNew; T fOld = f(x); int evals = 1;            // old and new value of f(x), # evaluations
  bool converged = false;
  while(!converged)
  {
    // compute gradient and Hessian at current position x:
    NumDif::gradient(f, x, N, g, h); evals += 2*N;
    NumDif::hessian( f, x, N, H, h); evals += 2*N*N + 1;
    // todo: make a function that computes both - this may save a couple of evals


    AT::negate(pg, pg, N);   // g becomes -g
    LinAlg::solve(H, d, g);  // solve H*d = -g for d which is required difference vector
    AT::add(x, pd, pX, N);   // add d to x and store result in X, our tentative new X

    // if the quadratic approximation to f has a minimum, the solution vector X is the minimum
    // location - but it can also be a maximum or saddle - we can figure this out, by evaluating
    // f at X and comparing the value to the old f at x....
    fNew = f(pX); evals += 1;

    if(fNew < fOld)
      rsArrayTools::copy(pX, x, N);
    else
    {
      AT::addWithWeight(x, N, pd, T(-1)); // jump away from the maximum - needs test
      // perhaps, we should also use a tentative X vector here and do an acceptance loop
      // maybe we should have a fallback to go a step along the negative gradient
    }
    // maybe instead of just checking if fNew < fOld, we should predict the change in f with our
    // quadratic approximation and compare that to the actual change - if it's too far off, it
    // means out approximation is not good and we should probably do a line search into the
    // gradient direction instead - the change can be up or down, though - depending on whether
    // we are near a minimum or maximum

    if(rsAbs(fOld-fNew) < tol)
      converged = true;
    fOld = fNew;
  }
  return evals;
}
// it looks like it's about time to make a class rsMinimizer

// Ideas:
// Gradient descent makes larger steps at steep cliffs and smaller steps in shallow regions because
// the length of the step is proportional to the norm of the gradient. That seems intuitively
// undesirable: into a direction where there is a steep decay, we may actually want to take a
// smaller step than into a direction where there's only a shallow decrease. Shouldn't we desire
// the step-sizes along the various directions to be inversely proportional to the steepnesses
// along these directions? But near a minimum where the steepness goes to zero, the stepsize
// would become infinite....hmmm...
// take f(x,y) = x^2 + 100*y^2 -> fx(x,y) = 2*x, fy(x,y) = 200*y
// at (x,y) = (2,2) - the gradient is (2*2,200*2) = (4,400) and we certainly don't want the y-step
// to 100 times larger than the x-step - we actually want them to be the same - the ideal step
// would be given by the vector (-2,-2) which would let us jump into the minimum at (0,0)
// ..maybe we should use the element-wise sign of the gradient, multiplied by a scalar stepsize?
// this would make all the steps have the same length...maybe each direction should have its
// own stepsize-scaler which increases when two successive steps have the same sign for that
// direction and decreases when they have opposite sign? i think, the effect would be similar to
// using a smoothed gradient as in the "rolling ball" algorithm

//=================================================================================================







//=================================================================================================

/** Class for representing floating point numbers as fractional and integer part to avoid
precision loss for greate numbers.. */

/*
template<class T>
struct rsFractionalIndex2
{
  rsFractionalIndex2 operator+(const rsFractionalIndex2& b)
  {
    rsFractionalIndex2 c;
    c.i = this->i + b.i;
    c.f = this->f + b.f;
    if(c.f >= 1) { c.f -= 1; c.i += 1; }
    if(c.f <  0) { c.f += 1; c.i -= 1; }
    return c;
  }

  int i;  // integer part of the number
  T   f;  // fractional part of the number
};
*/
