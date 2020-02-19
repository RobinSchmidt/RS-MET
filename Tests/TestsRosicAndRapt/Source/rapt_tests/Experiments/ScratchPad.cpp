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
// move to rsArray, make a cersion for complex numbers that does the same thing for real and 
// imaginary parts separately

//-------------------------------------------------------------------------------------------------
// Linear Algebra stuff:

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
  //bool reduced = false; // make parameter

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
// make member of rsMatrix, the ii,jj are ugly - replace with i,j - maybe renqame current i,j into
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

/** Returns true, if each of the indices form 0 to numIndices-1 (both inclusive) are contained 
exactly once in the idices array. That means the array is a permutation of the numbers from
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
  Q.setSize(n, n);
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
  V.setSize(n, n);
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
  S.setSize(m, n);
  S.setToZero();
  //for(i = 0; i < rsMin(m, n); i++)
  for(i = 0; i < rsMin(r, n); i++)    // for i >= r, sqrt(lambda[i]) == 0 - no need to compute it
    S(i, i) = sqrt(lambda[i]);
  // todo: avoid taking square-roots of zero - just loop up to i < rsMin(m, r)


  // Construct matrix U = (u_1,...,u_m) where u_1,..,u_r are computed from the nonzero singular 
  // values sigma_i and corresponding basis-vectors v_i as: u_i = (1/sigma_i) * A * v_i and the 
  // remaining u_{r+1},...,u_m (if any) are a basis of the orthogonal complement of u_1,..,u_r:
  U.setSize(m, m);
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
// x is initial etsimate, y is target value for y
template<class T>
T newton(const std::function<T(T)>& f, T x, T y = T(0))
{
  static const int maxNumIterations = 100;
  T tol = std::numeric_limits<T>::epsilon();
  T h   = 1.e-8;  // make parameter
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



template<class T>
T squaredDistance(T x1, T y1, T x2, T y2)
{
  T dx = x2 - x1;
  T dy = y2 - y1;
  return dx*dx + dy*dy;
}

template<class T>
T distance(T x1, T y1, T x2, T y2)
{
  return sqrt(squaredDistance(x1, y1, x2, y2));
}
// Euclidean distance

/** Computes the minimum value of the squared distances from (x0,y0) to the points in the x,y 
arrays. */
template<class T>
T minSquaredDistance(T x0, T y0, T* x, T* y, int N)
{
  T d2min = RS_INF(T);
  for(int i = 0; i < N; i++) {
    T d2 = squaredDistance(x0, y0, x[i], y[i]);
    if(d2 < d2min)
      d2min = d2; }
  return d2min;
}

template<class T>
T minDistance(T x0, T y0, T* x, T* y, int N)
{
  return sqrt(minSquaredDistance(x0, y0, x, y, N));
}


/** Given plane coordinates x,y, this function computes a height above the plane that has the shape
of ridge (of height 1) that spirals around in a logarithmic spiral with the parametric equation:
  x(t) = exp(a*t) * cos(sign * t + p); y(t) = exp(a*t) * sin(sign * t + p)
For points that are on the spiral, the function will return zero and for points that are "halfway" 
in between two "arcs", it will return 1. Halfway is to be understood in the logarithmic sense - for 
example, if (1,0) and (2,0) are points on the spiral, the point (sqrt(2),0) would be considered 
halfway between them. If the exponential growth parameter "a" is equal to log(2)/(2*pi), the spiral 
will grow by a factor of 2 in each revolution. The "sign" parameter should be +1 or -1 and 
determines the direction of the rotation. */
/*
double spiralRidge(double x, double y, double a = 1.0, double p = 0.0, double sign = 1.0, 
  int profile = 0, double exponent = 1.0)
{
  // sanity check inputs:
  rsAssert(sign == +1 || sign == -1,     "sign must be +-1");
  rsAssert(profile >= 0 && profile <= 3, "invalid profile");

  // compute raw height:
  double r = sqrt(x*x + y*y);
  if(r == 0.0) return 0.0;             // avoid log-of-zero
  double t  = log(r) / a;              // parameter t for point on the spiral with radius r
  double xs = r * cos(sign * t + p);   // x on the spiral for the given t
  double ys = r * sin(sign * t + p);   // y on the spiral for the given t
  double d  = distance(xs, ys, x, y);  // distance of input point to point on the spiral
  double h  = pow(0.5*d/r, exponent);  // height

  // apply shaping of the height profile:
  if(profile == 2) return h;                        // 2: rectified sine (comes out by raw formula)
  if(profile == 3) return 1-h;                      // 3: inverted rectified sine
  h = asin(h) / (0.5*PI);                           // convert to triangular
  if(profile == 0) return h;                        // 0: triangular
  if(profile == 1) return 0.5*(sin(PI*(h-0.5))+1);  // 1: sinusoidal
  return 0;                                         // unknown profile
}
// moved to rsImageGenerator
*/

// optimize: the sqrt in the distance computation can be avoided: compute the distance-squared and 
// then use 0.5*exponent in the subsequent pow call
//
// The algo computes the distance of (x,y) to a point on the spiral that has the same radius as 
// (x,y). It happens that the height-profile (as function of radius for a given angle) comes out as
// a rectified sine shape (when the radius is used a x-axis and the x-axis is logarithmically 
// scaled)
// what if sign is not +-1? what if we use different factors for x- and y: 
//   xs = r * cos(wx * t + px); ys = r * sin(wy * t + py);
// ..in this case, the profile computations will very likely become invalid because the raw profile 
// is not a rectified sine anymore - yes - using, for example cos(2*..),sin(3*..) gives a nice 
// effect - one should increase the shrink-factor accordingly because it get denser otherwise
// for the original profile and the profile converted to a full sine, it makes visually not 
// qualitative difference, when we invert all color channels - for the triangular profile, it does.
// maybe to create audio-signals, we could use dx = (xs-x)/r; dy = (ys-y)/r; as left and right 
// channel signal - but what should the input be? we don't have x,y pixel coordinates as inputs but
// time instants
// 
// when we use d / r, the birghtness of the white ridges is independent for the distance to the 
// center - using a power with exponent < 1, we get a darkening effect towrd the center - but mybe 
// such an effect can be applied as post-processing: 
// circularDarkening(img, x, y, amount)
//   img(i,j) /= pow(r, amount)
//
// try tL = atan2(y,x) + 2*k*pi where k = floor(t0/(2*pi)), tR = tL + 2*pi, compute (xL,yL),(xR,yR)
// by the parametric spiral equations, compute distances dL,dR and use minimum
//
// -these are not the actual distances to the nearest points on the spiral but rather the distances 
//  to two concentric circles that approximate the spiral at the given angle - but they can be used 
//  as an initial estimate for computing the actual distance via netwon iteration - maybe this 
//  refinement can be made optional, controlled by a boolean parameter
//
// see:
// https://en.wikipedia.org/wiki/Logarithmic_spiral
// https://en.wikipedia.org/wiki/Archimedean_spiral



// obsolete:
double spiralRidge2(double x, double y, double a = 1.0, double p = 0.0, double sign = 1.0)
{
  double r = sqrt(x*x + y*y);
  if(r == 0.0) return 0.0;            // avoid log-of-zero
  double t0  = log(r) / a;
  double k   = floor(t0/(2*PI));
  double phi = rsWrapToInterval(atan2(y,x), 0, 2*PI);

  double tL  = phi + k*2*PI;

  // should not be needed (?):
  while( tL > t0       )  tL -= 2*PI;
  while( tL + 2*PI < t0)  tL += 2*PI;


  // test:
  tL = t0;
  double inc = 0.1;
  int maxNumIterations =  (int) ceil(2*PI / inc);
  int i = 0;
  while(true)
  {
    double xL = exp(a*tL)*cos(tL);
    double yL = exp(a*tL)*sin(tL);
    double pL = rsWrapToInterval(atan2(yL, xL), 0, 2*PI);

    if(i > maxNumIterations)     break;

    
    if(y >= 0)  {
      if(pL < p + 2*PI)
        break; }
    else {
      break;  // test - i expect the whole y < 0 halfplane to be colored black - but it gets colored correctly - why?
      if(pL < p)   
        break;  
    }
    

    tL -= inc;
    i++;
  }

  double tR = tL + 2*PI;
  rsAssert(tL <= t0 && tR > t0);
  // sanity check - this does sometimes trigger - maybe the formula is imprecise? maybe we should 
  // do: if(tL > t0) tL -= 2*PI

  double xL = exp(a*tL) * cos(sign * tL + p);
  double yL = exp(a*tL) * sin(sign * tL + p);
  double xR = exp(a*tR) * cos(sign * tR + p);
  double yR = exp(a*tR) * sin(sign * tR + p);
  double dL = distance(xL, yL, x, y);
  double dR = distance(xR, yR, x, y);

  // test:
  double phiL = rsWrapToInterval(atan2(yL,xL), 0, 2*PI); // should be equal to phi
  double phiR = rsWrapToInterval(atan2(yR,xR), 0, 2*PI); // dito
  double rL   = sqrt(xL*xL + yL*yL);                     // shold be <= r
  double rR   = sqrt(xR*xR + yR*yR);                     // should be > r
  // phiL and phiR are equal to each other but different from phi - maybe we indeed need the loop
  // form the function below


  double d  = rsMin(dL, dR);
  //double d = (dL + dR) / 2;  // test

  //return dL / r;

  return d / r;  // use pow(r, distanceWeight)
}
// currently jumps discontinuously from black to white
// todo: echk, if (xL,yL),(xR,yR) have the same angle as (x,y) - the idea is that these two points
// should be *on* spiral arms at the same angle ats x,y and one should be on the inner and one on
// the outer arm, with respect to point x,y - but apparently, this doesn't work yet
// maybe instaed of L,R use I,O for inner/outer

// obsolete:
double spiralRidgeOld(double x, double y, double a = 1.0)
{
  // under construction


  //if(y < 0)
  //  return 0;  // test

  function<double(double)> fx =  [=](double t) -> double { return exp(a*t)*cos(t); };
  function<double(double)> fy =  [=](double t) -> double { return exp(a*t)*sin(t); };
  //double t = 1.5;
  //double x1 = fx(t);
  //double y1 = fy(t);
  //GNUPlotter plt;
  //plt.plotCurve2D(fx, fy, 1000, 0.0, 50.0);

  // try an input point (x,y) = (50,30)
  //x = 50;
  //y = 30;


  //// test:
  //x = fx(t);
  //y = fy(t);


  // convert x,y to polar coordinates
  double r = sqrt(x*x + y*y);

  if(rsAbs(r) < 0.00001)  // ad hoc
    return 0;

  double p = rsWrapToInterval(atan2(y, x), 0, 2*PI);

  // find a value for the that corresponds to the given radius r:
  double t0 = log(r) / a;
  // -this value of t, when plugged into the parametric spiral equation, will lead a point on the 
  //  spiral that has the same radius as our given input point - but it will in general not have 
  //  the same angle p
  // -the task is now to find tL <= t0, tr >= t0 that give points *on* the spiral with the same p
  //  tR = tL + 2*pi - so we just need to find tL
  // -then compute the distances of the input point x,y to both of these points on the spiral and
  //  select the smaller of them as our distance

  // very crude algorithm to find tL
  double tL = t0;
  double inc = 0.1;

  
  int maxNumIterations =  (int) ceil(2*PI / inc);
  int i = 0;
  /*
  while(true)
  {
    double xL = exp(a*tL)*cos(tL);
    double yL = exp(a*tL)*sin(tL);
    double pL = rsWrapToInterval(atan2(yL, xL), 0, 2*PI);

    if(i > maxNumIterations)     break;

    if(y >= 0)  {
      if(pL < p)
        break; }
    else {
      break;  // test - i expect the whole y < 0 halfplane to be colored black - but it gets colored correctly - why?
      if(pL < p + 2*PI)   
        break;  
      // maybe this condition is alway true, so it doesn't make a difference - but it would mean, 
      // that
    }



    //if(pL < p + 2*PI && y <  0)  break;

    tL -= inc;
    i++;
  }
  */
  

  /*
  // other algo - but doesn't work:
  tL = 2*PI * floor(t0 / (2*PI));  // tL is multiple of 2*PI and tL <= t0
  while(true)
  {
    double xL = exp(a*tL)*cos(tL);
    double yL = exp(a*tL)*sin(tL);
    double pL = rsWrapToInterval(atan2(yL, xL), 0, 2*PI);
    if(pL >= p)
      break;
    tL += inc;
  }
  */


  double tR = tL + 2*PI;

  double xL = fx(tL);
  double yL = fy(tL);
  double xR = fx(tR);
  double yR = fy(tR);
  double dL = distance(xL, yL, x, y);
  double dR = distance(xR, yR, x, y);
  double d  = rsMin(dL, dR);

  //return dR;  // test

  return dL/r;  // test

  return d;

  // it works!!! but why?! it does not compute the distance i wanted - but the other distance that
  // it actually computes makes a nice spiral-ridge function, too. we may need to pass the output
  // through some nonlinare function that expands the center values - the extreme black and white 
  // values are a bit overrepresented and the middle gray values underrepresented - maybe something
  // based on the logistic function - or maybe one of these:
  // https://en.wikipedia.org/wiki/Generalised_logistic_function
  // https://en.wikipedia.org/wiki/Gompertz_function


  // i think, when the point (x,y) is in the lower half-plane (i.e. the angle is > pi), it may not 
  // work correctly because the the computed angles pL in the loop are *less* than p for points
  // on the spiral further outside - i think, instead of comparing pL < p, we must compare 
  // pL < p+2*pi when y < 0 
  // ..we still seem to jump out of the loop early under certain conditions - figure out what the
  // breaking conditions must be depending on the angle - we probably need different conditions
  // for when (x,y) is in the 4 different quadrants - in the top-right quadrant, checking pL < p
  // should be fine





  /*
  // find a point on the unit circle with the same angle as (x,y)
  double xn = x / r;  // normalized x
  double yn = y / r;  // normalized y

  double ac = acos(xn);  // should be equal to asin(yn)?
  */


  // estimate - try to come up with something better - maybe involving acos(x) and/or asin(y)?


  // -try 
  // -to do so, compute



  return 0; // preliminary





  /*
  // ...use newton iteration to refine t - we want to find a t such that
  // u(t) :=   2*(a*cos(t)*e^(a*t) - e^(a*t)*sin(t))*(cos(t)*e^(a*t) - x) 
  //        + 2*(a*e^(a*t)*sin(t) + cos(t)*e^(a*t))*(e^(a*t)*sin(t) - y) = 0

  // simplify: s = sin(t), c = cos(t), r = e^(a*t)
  // u(t) := 2*(a*c*r-r*s)*(c*r-x) + 2*(a*r*s+c*r)*(r*s-y) = 0

  function<double(double)> u;
  u = [=](double t)     // needs to capture "a"
  { 
    double s = sin(t);
    double c = cos(t);
    double r = exp(a*t);
    return 2*(a*c*r-r*s)*(c*r-x) + 2*(a*r*s+c*r)*(r*s-y); 
  };

  t = newton(u, t);
  */
}
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

// because the distances increase exponentially with the radius of (x,y), define weighted distance
//   D(x,y) = d(x,y) / exp( sqrt(x^2 + y^2) ) 
// the function R(x,y) = 1 / (1 + D^2) is a sort of ridge that has the shape of the exponential 
// spiral, so we may use R(x,y) = c it as our implicit equation that should have exponential 
// spirals as level lines
//   ...(check, if this weighting is good) - we wnat something that goes exponentially to infinity
//   as x^2+y^2 goes to zero and exponentially to zero as x^2 + y^2 goes to infinity maby
//   exp(1 / (x^2+y^2))

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

// strategy to find the distance:
// -input: the point x,y
// -find values xl < x, xr >= x such that f(t) = exp(a*t)*cos(t) = xl or xr
// -find a vlaue for t such that exp(a*t)*cos(t) = x, or maybe two values - incre

// or maybe use a simpler linear spiral:
//   f(t) = t*cos(t), g(t) = t*sin(t)