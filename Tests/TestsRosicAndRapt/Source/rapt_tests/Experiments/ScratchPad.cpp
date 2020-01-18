// contains some unfinished code "scratches" and code under construction which may eventually be 
// moved to somewhere else once it works and does something useful

using namespace std;

template<class T>
void cleanUpIntegers(T* a, int N, T tol)
{
  for(int i = 0; i < N; ++i) {
    T rounded = round(a[i]);
    if( rsAbs(a[i] - rounded) <= tol )
      a[i] = rounded; }
}
// move to rsArray

//-------------------------------------------------------------------------------------------------
// Linear Algebra stuff:



/** Returns true, if the space spanned by the columns of x is within the span of the columns of B.
That means each column of x can be expressed as some linear combination of the columns of B. */
template<class T>
bool isInSpanOf(rsMatrix<T> B, rsMatrix<T> x, T tol)
{
  //rsAssert(B.hasSameShapeAs(x));
  rsAssert(B.getNumRows() == x.getNumRows());
  RAPT::rsLinearAlgebraNew::makeTriangular(B, x);
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
// rename to getLeadCoeffIndex


/** Puts the augmented coefficient matrix A|B into row echelon form via Gaussian elimination. 
In this form, each row has its leading coefficient at least one position further to the right than
the previous row and rows of all zeros are at the bottom. */
template<class T>
void rowEchelon2(rsMatrixView<T>& A, rsMatrixView<T>& B, T tol)  // change to rsMatrixView for production
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
void rowEchelon2(rsMatrixView<T>& A, T tol)
{
  rsMatrix<T> dummy(A.getNumRows(), 1);
  rowEchelon2(A, dummy, tol);
}

// make a function reducedRowEchelon that also includes a backward/upward elimination pass - maybe 
// that can be integrated as option into the function above - if reduced == true also eliminate 
// upward - does that work

/*
// obsolete - doesn't work when matrix has columns of zeros - use rowEchelon2:
template<class T>
void rowEchelonOld(rsMatrixView<T>& A, rsMatrixView<T>& B, T tol)
{
  //T tooSmall = T(1000) * RS_EPS(T) * A.getAbsoluteMaximum();    // ad hoc -> todo: research

  int numRows = A.getNumRows();
  int numCols = A.getNumColumns();
  for(int i = 0; i < rsMin(numRows, numCols); i++) {
    //rsMatrix<T> dbg; dbg.copyDataFrom(A);  // uncomment for debugging

    int p = i;

    // search pivot row (maybe factor out int inte findPivot(int column)
    T biggest = T(0);
    for(int j = i; j < numRows; j++) {                          
      if( rsGreaterAbs(A(j, i), biggest) ) { 
        biggest = A(j, i); 
        p = j; }}

    if(p != i) {
      A.swapRows(i, p);
      B.swapRows(i, p); }

    for(int j = i+1; j < numRows; j++) {              // pivot row subtraction
      if( !rsGreaterAbs(A(i, i), tol) )
        break;  // this is an all-zeros column
      T w = -A(j, i) / A(i, i);                       // weight
      A.addWeightedRowToOther(i, j, w, i, numCols-1); // start at i: avoid adding zeros
      B.addWeightedRowToOther(i, j, w);
    }
  }
}
template<class T>
void rowEchelonOld(rsMatrixView<T>& A, T tol)
{
  rsMatrix<T> dummy(A.getNumRows(), 1);
  rowEchelon(A, dummy);
}
*/





template<class T>
RAPT::rsPolynomial<T> getCharacteristicPolynomial(const rsMatrixView<T>& A)
{
  using RatFunc = RAPT::rsRationalFunction<T>;
  using Matrix  = RAPT::rsMatrix<RatFunc>;

  // Create matrix B = A - x*I as matrix of rational functions:
  Matrix B(A.getNumRows(), A.getNumColumns());
  for(int i = 0; i < B.getNumRows(); ++i)
    for(int j = 0; j < B.getNumColumns(); ++j)
      if(i == j)
        B(i, j) = RatFunc({A(i, j), -1}, {1});
      else
        B(i, j) = RatFunc({A(i, j)},     {1});

  // Create a dummy right-hand-side (todo: allow function to be called without rhs) - maybe
  // make an LU decomposition function that fills an array of swaps at each index, the number says
  // with which row this row has to be swapped
  Matrix R(A.getNumRows(), 1);
  for(int i = 0; i < R.getNumRows(); ++i)
    R(i, 0) = RatFunc({ 1 }, { 1 });

  // compute row echelon form of B:
  //RAPT::rsLinearAlgebraNew::makeTriangularNoPivot(B, R);
  makeTriangularNoPivot(B, R);
  // i think, we really need pivoting for this here too - we may encounter the zero-function - but 
  // maybe we should swap only when the zero-function is encountered - i.e. we don't search for the
  // largest element - instead, we check against zero and if we encounter a zero, we search for the 
  // next nonzero element to swap with - the current code should not go into the library - it's
  // useless for production - see weitz book pg 310

  // Compute determinant. For a triangular matrix, this is the product of the diagonal elements. 
  // The computed determinant is still a rational function but it should come out as a polynomial, 
  // i.e. the denominator should have degree 0 (be a constant). I think, it should always be +1 or
  // -1 because the elementary row operations can only flip the determinant.
  RatFunc d = B.getDiagonalProduct();
  rsAssert(d.getDenominatorDegree() == 0);
  return d.getNumerator() / d.getDenominator()[0]; 
}
// todo: make a version that uses Laplace expansion of the determinant with a matrix of polynomials
// -> should give the same result 

/** Represents the root of a polynomial along with its multiplicity. The datatype of the 
coefficients is assumed to be a complex number type. */
/*
template<class T>
struct rsPolynomialRoot
{
  T value;   // T is assumed to be a complex type
  int multiplicity;
};
*/

/** Represents the eigenspace of a matrix with complex coefficients. Each eigenspace consists of 
an eigenvalue and an associated set of eigenvectors. */
template<class T>
struct rsEigenSpace
{

  //int getAlgebraicMultiplicity() const { return eigenvalue.multiplicity; }

  int getAlgebraicMultiplicity() const { return algebraicMultiplicity; }

  //int getGeometricMultiplicity() const { return (int) eigenspace.getNumColumns(); }

  //void orthonormalize();

  //rsPolynomialRoot<T> eigenvalue;


  complex<T> eigenvalue;
  rsMatrix<complex<T>> eigenspace; 
  //rsMatrix<T> eigenspace;  // basis of nullspace of A - eigenvalue * I
  // should be a complex matrix

               
  int algebraicMultiplicity; // not yet assigned
};

/** Represents the set of eigenspaces of a matrix with complex coefficients. This set consists of a
set of (eigenvalues, each with an algebraic multiplicity which is the order of the
polynomial root. Associated with each eigenvalue is a set of eigenvectors...  */
/*
template<class T>  // T should be a complex type
class rsEigenSpaceSet
{

public:

  rsPolynomial<T> characteristicPolynomial;
  std::vector<rsEigenSpace<T>> eigenspaces;

protected:

};
*/

// for matrices with real coefficients, we must promote them to complex numbers because even real
// matrices may have complex eigenvalues and eigenvectors
/*
template<class T>
class rsEigenStuffReal : public rsEigenSpaceSet<std::complex<T>>
{

};
// 
*/

template<class T> 
RAPT::rsMatrix<complex<T>> complexify(const RAPT::rsMatrix<T>& A)
{
  RAPT::rsMatrix<complex<T>> Ac(A.getNumRows(), A.getNumColumns());
  Ac.copyDataFrom(A);
  return Ac;
}
// needed for technical reasons

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

template<class T> 
vector<complex<T>> getPolynomialRoots(const RAPT::rsPolynomial<T>& p)
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

// this function should probably take a complex polynomial as input, so we can use it with complex
// matrices too


template<class T> 
vector<complex<T>> getEigenvalues(const rsMatrixView<T>& A)
{
  RAPT::rsPolynomial<T> p = getCharacteristicPolynomial(A);
  return getPolynomialRoots(p);
}


/*
// i think, this is still very wrong:
template<class T>
vector<vector<complex<T>>> getEigenvectors(const rsMatrixView<T>& A)
{
  // for each eigenvalue x, we must form the matrix A - x*I and solve the system
  // A - x*I = 0 where 0 is the zero vector in this case - is that correct?
  // ...but the eigenvalues are complex - does that mean the eigenvectors may also be complex?
  // probably

  rsAssert(A.isSquare());

  using LA = RAPT::rsLinearAlgebraNew;
  int N = A.getNumRows();

  vector<complex<T>> eigenValues = getEigenvalues(A);
  rsAssert((int)eigenValues.size() == N);
  vector<complex<T>> rhs(N);
  rsMatrix<complex<T>> I(N, N), B(N, N);
  B.copyDataFrom(A);
  vector<vector<complex<T>>> eigenVectors(N); 

  //rhs[0] = 1;
  // i think, we can make a choice - or maybe even several choices - we are trying to solve a
  // singular system of equations when solving (A - x_i * I) * x = 0, so our regular Gaussian
  // elimination procedure is not suitable - but what else?

  // i think, this is wrong:
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
      B(j, j) = A(j, j) - eigenValues[i];   // set diagonal elements of B = (A - x_i * I)
    eigenVectors[i] = LA::solve(B, rhs);    // solve B * x == 0
  }
  // the system we are trying to solve is singular and actually the rhs is zero anyway - what we
  // need instead is to find the nullspace of A - x_i * I = B_i for each i

  // triggers assertion that matrix is singular -> try elimination by hand and compare to what
  // the algo does
  // let B_i be given by A - x_i * I, then B_0 = [[-3,6],[-3,6]], B_1 = [[-6,6],[-3,3]]
  // so yes, the shifted matrix is indeed singular - seems we can't solve for the eigenvectors by 
  // regular linear system solving...hmmm - how else can we find them? maybe we need to find 
  // coefficients for which the linear-combination of the rows gives the zero-row?

  // also - in general, there is no one-to-one correspondence between eigenvalues and eigenvectors
  // there may be zero or multiple eigenvectors to any given eigenvalue - the set of eigenvectors
  // belonging to eigenvalue x_i is the eigenspace of x_i and the dimensionality of this space
  // is the geometric multiplicity of x_i

  // ah - i think the problem is that an eigenvector represents infinitely many solutions to a 
  // linear system because we can scale it by an arbitrary factor and still have a solution - so we
  // have a situation where there's a parametric continuum of solutions - the number of free 
  // parameters is the dimenstion of the eigenspace? and that value is probably related to
  // rank(A - x_i * I)

  // ----INFO STILL RELEVANT - DON'T DELETE----:
  // It seems, we may not have to solve any system of equations as all - here it says:
  // https://en.wikipedia.org/wiki/Eigenvalue_algorithm
  // 
  // the columns of the matrix prod_{i != j} (A - x_i * I)^m_i must be either 0 or generalized 
  // eigenvectors of the eigenvalue x_j - i've adapted the notation a bit - but note that the i
  // runs over the *distinct* eigenvalues, not just over all eigenvalues ..soo, it seems to find 
  // the eigenspace to an eigenvalue x_i, we must do repeated matrix multiplication rather than
  // solving a linear system - or is there a better way? ...solving a linear system is actually 
  // wrong - we need to find a nullspace

  // see also:
  // https://lpsa.swarthmore.edu/MtrxVibe/EigMat/MatrixEigen.html
  // https://www.scss.tcd.ie/Rozenn.Dahyot/CS1BA1/SolutionEigen.pdf
  // http://www.sosmath.com/matrix/eigen2/eigen2.html
  // http://wwwf.imperial.ac.uk/metric/metric_public/matrices/eigenvalues_and_eigenvectors/eigenvalues2.html

  return eigenVectors;
}
// http://doc.sagemath.org/html/en/constructions/linear_algebra.html
*/


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
// make member of rsMatrix

// for computing nullspaces, we will need a function that removes zero columns - it sould take an
// array of th indices of the column to remove
// rsMatrix<T> getWithRemovedColumns(const rsMatrix<T>& A, const std::vector<int>& colIndices)

/** Expands the determinant column-wise along the j-th column. This is the textbook method and has
extremely bad scaling of the complexity. The function calls itself recursively in a loop (!!!). I 
think, the complexity may scale with the factorial function (todo: verify) - which would be 
super-exponential - so it's definitely not meant for use in production code. For production, use
Gaussian elimination (we need to keep track of whether we have an odd or even number of swaps - in 
the former case det = -1 * product(diagonal-elements of upper triangular form), in the later case
det = +1 * product(...) */
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
      if( rsGreaterAbs(A(i, j), tol) )  {  // we need a tolerance
        nonZeroElemFound = true;
        break; } // i-th row is not all-zeros
      j++; }
    if(!nonZeroElemFound)
      return i;  // no non-zero element was found - i-th row is all zeros
    i++; }
  return i;
}
// 
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
form (in which case LA::makeTriangularLA::makeTriangular will do nothing and return 0) */
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
infinietly many. The full set of solutions if given by that (or another) particular solution plus 
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
  rowEchelon2(A, B, tol);
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



/** Computes the set of vectors v which solve the homogenous linear system of equations A * v = 0 
where v is some vector and 0 is the zero vector. We assume v to be M-dimensional, so the matrix A 
must have M columns. The set of vectors v that solve this equation will in general span a subspace
of R^M. This subspace is called the nullspace of the matrix A. This function returns a basis for 
this subspace represented as matrix. The columns of the matrix are the basis vectors. Note that 
this may be the empty matrix which indicates that the nullspace of A consists only of the 
zero-vector (todo: maybe we should return the Mx1 zero vector in this case? ...decide later by 
which convention is more convenient when dealing with eigenspaces */
template<class T>
rsMatrix<T> getNullSpace(rsMatrix<T> A, T tol)
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
  // not -> FIGURE OUT). We set up an RxR linear system solve it for N different right hand sides
  // corresponding to N different choices for assigning the free parameters. The most natural 
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
  //rowEchelon(A);  // old
  rowEchelon2(A, tol);  // new

  // maybe factor out a function that takes a triangular matrix getNullSpaceEchelon(Matrix&)

  // find out which dimensions are free and which dependent:
  std::vector<int> pivots = getPivots(   A, tol);  // indices of dependent variables
  std::vector<int> params = getNonPivots(A, tol);  // indices of free variables (parameters)
  int nEqn = (int) pivots.size();                  // number of equations R (#dependents, rank)
  int nRhs = (int) params.size();                  // number of right-hand sides N (#parameters)

  // sanity checks for debug:
  rsAssert(isIndexSplit(A.getNumColumns(), pivots, params));
  rsAssert(getRankRowEchelon(A, tol) == (int) pivots.size());

  // set up the linear system ("gather") and solve it:
  Matrix M(nEqn, nEqn);                            // coefficient matrix of the NxN system
  Matrix R(nEqn, nRhs);                            // right hand side matrix
  Matrix b(nEqn, nRhs);                            // solution
  int i, j;
  for(i = 0; i < nEqn; i++) {
    for(j = 0; j < nEqn; j++)
      M(i, j) =  A(i, pivots[j]);                  // copy relevant coeffs
    for(j = 0; j < nRhs; j++)
      R(i, j) = -A(i, params[j]); }                // resulting rhs from setting j-th param to 1
  bool success = solve2(M, b, R, tol);
  rsAssert(success);

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



/*
template<class T>
rsEigenSpaceSet<complex<T>> getEigenSpaces(rsMatrix<T> A, T tol)
{
  rsAssert(A.isSquare());  // try to relax later

  rsEigenSpaceSet<complex<T>> eigenSet;
  using Matrix  = RAPT::rsMatrix<T>;
  using MatrixC = RAPT::rsMatrix<complex<T>>;

  // get characteristci polynomial and eigenvalues:
  RAPT::rsPolynomial<T> p = getCharacteristicPolynomial(A);

  //eigenSet.characteristicPolynomial = p; 
  // still incomatible with respect to what is complex and what is real - perhaps is easiest to
  // just make everything complex

  vector<complex<T>> eigenvalues = getPolynomialRoots(p);

  // convert roots array to an array that has unique roots with a multiplicity currently, we compute
  // each eigenspace multiple times - once for each occurence of its eigenvalue - that's wasteful

  // compute eigenspaces to each eigenvalue:
  int N = A.getNumRows();              // or columns?
  MatrixC Ac = complexify(A);
  MatrixC I(N,N); I.setToIdentity();


  for(int i = 0; i < N; i++) {
    MatrixC Ai = Ac - eigenvalues[i] * I;        // (A - x_i * I) * v = 0

    MatrixC eigenSpace = getNullSpace(Ai, complex<T>(tol));  //  get rid

    rsEigenSpace<complex<T>> es;
    es.eigenspace = eigenSpace;
    es.eigenvalue = eigenvalues[i];
    eigenSet.eigenspaces.push_back(es);

    // test, if A * v = eigenvalue * v for v in the eigenspace at once
    MatrixC test1 = Ac * eigenSpace;
    MatrixC test2 = eigenvalues[i] * eigenSpace;
    MatrixC error = test1 - test2;
    rsAssert(error.isZero(tol));

    //eigenspaces[i] = eigenspaces[i].getTranspose(); // for convenient inspection in the debugger
  }           

  return eigenSet;
}
// i get linker errors due to rsPolynomial<complex<double>>
*/


template<class T>
rsMatrix<T> getOrthogonalComplement(rsMatrix<T> A, T tol)
{
  return getNullSpace(A.getTranspose(), tol);  // verify if that's correct
}
// needs test






















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

