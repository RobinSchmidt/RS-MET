template<class T>
std::vector<T> rsLinearAlgebraNew::solve(const RAPT::rsMatrixView<T>& A_, const std::vector<T>& b_)
{
  int N = (int) b_.size();
  std::vector<T> x(N), b = b_; // temporaries
  rsMatrix<T> A(N, N, A_.getDataPointerConst());
  rsMatrixView<T> vx(N, 1, &x[0]), vb(N, 1, &b[0]);
  solve(A, vx, vb);
  return x;
}

template<class T>
RAPT::rsMatrix<T> rsLinearAlgebraNew::inverse(const RAPT::rsMatrixView<T>& A)
{
  rsAssert(A.isSquare());
  int N = A.getNumRows();
  RAPT::rsMatrix<T> tmp(N, N, A.getDataPointerConst()), E(N, N);
  E.setToIdentity();
  solve(tmp, E, E);
  return E; 
}

template<class T>
bool rsLinearAlgebraNew::solve(rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B)
{
  rsAssert(A.getNumColumns() == X.getNumRows()); // A*X = B or A*x = b must make sense
  rsAssert(X.hasSameShapeAs(B));                 // num of solutions == num of rhs-vectors
  rsAssert(A.isSquare()); 
  if(makeTriangular(A, B) != A.getNumRows())
    return false;                                // matrix A was singular -> report failure
  solveTriangular(A, X, B);
  return true;                                   // matrix A was regular -> report success
}

template<class T>
T rsLinearAlgebraNew::determinant(const rsMatrixView<T>& A)
{
  int N = A.getNumRows();
  rsAssert(A.getNumColumns() == N, "Matrix A must be square");
  rsMatrix<T> A2(N, N, A.getDataPointerConst()), B(N, 1);
  int numSwaps;
  int rank = makeTriangular(A2, B, &numSwaps);
  if(rank < N)
    return T(0);
  T det = A2.getDiagonalProduct();
  if(rsIsOdd(numSwaps))
    return -det;
  else
    return det;
}


// subspaces:

/*
template<class T>
rsMatrix<T> rsLinearAlgebraNew::rowSpace(rsMatrix<T>, T tol)
{
  rsMatrix<T> z(A.getNumRows(), 1);    // dummy
  RAPT::rsLinearAlgebraNew::makeTriangular(A, z);
  return getWithoutBottomZeroRows(A, tol);
}
*/




// subroutines:

template<class T>
int rsLinearAlgebraNew::makeDiagonal(rsMatrixView<T>& A, rsMatrixView<T>& B)
{
  int N = makeTriangular(A, B);
  for(int i = N-1; i >= 0; i--) {
    T s = T(1) / A(i, i);
    A.scaleRow(i, s);
    B.scaleRow(i, s);
    for(int j = i-1; j >= 0; j--) {
      s = -A(j, i);
      A.addWeightedRowToOther(i, j, s);
      B.addWeightedRowToOther(i, j, s); }} 
  return N;
  // here, also restrict the column-range as optimization (see makeSystemUpperTriangular)
  // A.addWeightedRowToOther(i, j, s, i, N-1)? or (i, j, s, j, N-1)? or (i, j, s, j+1, N-1)?
  // implement unit test and figure out
  // what about pivoting? what if A(i,i) == 0? makeTriangular already makes sure that there are 
  // nonzero diagonal elements down to the row N-1, where it stopped
}
// the reduced row echelon form is unique:
// https://en.wikipedia.org/wiki/Row_echelon_form
// todo: scale the rows by their leading coeff

/*
template<class T>
inline bool operator>(const std::complex<T>& L, const std::complex<T>& R)
{
  if(L.real() > R.real()) return true;
  if(L.imag() > R.imag()) return true;
  return false;
}
*/

template<class T>
int rsLinearAlgebraNew::makeTriangular(rsMatrixView<T>& A, rsMatrixView<T>& B)
{
  int numSwaps;  // dummy
  return makeTriangular(A, B, &numSwaps);
}

template<class T>
int rsLinearAlgebraNew::makeTriangular(rsMatrixView<T>& A, rsMatrixView<T>& B, int* numSwaps)
{
  rsAssert(A.getNumRows() == B.getNumRows());
  *numSwaps = 0;
  T tooSmall = T(1000) * RS_EPS(T) * A.getAbsoluteMaximum();    // ad hoc -> todo: research
  int i, numRows = A.getNumRows();
  for(i = 0; i < rsMin(numRows, A.getNumColumns()); i++) {
    //rsMatrix<T> dbg; dbg.copyDataFrom(A);  // uncomment for debugging
    int p = i; 
    T biggest = T(0);
    for(int j = i; j < numRows; j++) {                          // search pivot row
      if( rsGreaterAbs(A(j, i), biggest) ){ 
        biggest = A(j, i); 
        p = j; }}
    if(rsIsCloseTo(biggest, T(0), tooSmall))                    // no pivot found - return early
      return i;
    if(p != i) {                                                // turn pivot row into current row
      A.swapRows(i, p); 
      B.swapRows(i, p);
      (*numSwaps)++;     }                                      // keep track of number of swaps
    for(int j = i+1; j < numRows; j++) {                        // pivot row subtraction
      T w = -A(j, i) / A(i, i);                                 // weight
      A.addWeightedRowToOther(i, j, w, i, A.getNumColumns()-1); // start at i: avoid adding zeros
      B.addWeightedRowToOther(i, j, w); }}
  return i;
}
// maybe in if(rsIsCloseTo... we should not return early, if at the same time A(i,i) is zero - in 
// this case the i-th column is already zero from i downward - this is ok - or wait - no - this
// check is already includes in the for(int j=i ...loop
// pass a tol

// Maybe allow the function to be called without an rhs B. It may make sense to use it with a 
// single input in order to compute determinants - when the function returns, the determinant is
// the product of diagonal elements - up to a sign flip, which occurs when we had an odd number of
// row-swaps. Maybe keep track of whether the number of swaps was even or odd by doing
// oddSwaps *= -1 in if(p !=i) and return +1 or -1 - or better: return the determinant! if we
// run into the error branch, immediately return zero - but no - often, we don't need the 
// determinant, so this extra computation should be avoided - but maybe return the rank which is
// the iteration number i - callers may look at it and if it's less than N, they conclude that
// the matrix was singular
// instead of actually writing the zeros into the rows below, write the weights w - 
// addWeightedRow should then start at i+1 instead of at i - doing it this way produces the LU
// decomposition of a permutation of A ...but in order to be useful for later solving other 
// systems with ethe same matrix but other right-hand-sides, we would need to keep track of the
// permutations - here it is no problem, because we immediately apply the swaps to the RHS as well

// i think, this function is useful also for singular matrices - in this case, it should stop as
// soon as it encounters a situation where there are only zeros in th i-th column below the 
// diagonla element A(i,i) such that no pivot can be found - ith should then return i - it should
// always return the number of successful elimination steps - how does this number relate to the 
// rank - it can't be the rank itself because when the matix already is triangular, we take no step
// at all but it may still have full rank -  i think, the rank is given by the greatest index i for
// which A(i,i) is nonzero after the function returns

/*
template<class T>
void rsLinearAlgebraNew::makeTriangularNoPivot(rsMatrixView<T>& A, rsMatrixView<T>& B)
{
  rsAssert(A.isSquare()); // can we relax this?
  for(int i = 0; i < A.getNumRows(); i++) {  
    rsAssert(!(A(i, i) == T(0)), "This matrix needs pivoting");
    for(int j = i+1; j < A.getNumRows(); j++) {
      T w = -A(j, i) / A(i, i);
      A.addWeightedRowToOther(i, j, w, i, A.getNumColumns()-1); // start at i: avoid adding zeros
      B.addWeightedRowToOther(i, j, w); }}
}
*/
// maybe we should have a version that swaps when encountering a zero, maybe make an even simpler
// version without the test for zero - if client code knows that the matrix is in a form such that
// a zero pivot is never encountered, it may use this version (for optimization purposes)
// ...or replace the test by an assertion and make the function void


template<class T>
void rsLinearAlgebraNew::solveTriangular(
  rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B)
{
  T tol = 1.e-12;
  int M = X.getNumColumns();  // number of required solution vectors
  int N = A.getNumRows();     // number of elements in each solution vector
  for(int k = 0; k < M; k++) {
    for(int i = N-1; i >= 0; i--) {
      if( rsIsCloseTo(A(i,i), T(0), tol) ) {
        X(i, k) = T(0);       
        continue;      }
      T tmp = T(0);
      for(int j = i+1; j < N; j++)
        tmp += A(i, j) * X(j, k);
      X(i, k) = (B(i, k) - tmp) / A(i, i); }}
}
// what if A(i,i) == 0? it means that x_i doesn't occur in that equation for b_i, so it can have 
// any value - maybe we should simply set it zero? or one?

//-------------------------------------------------------------------------------------------------

template<class T>
std::vector<T> rsLinearAlgebraNew::solveOld(rsMatrix<T> A, std::vector<T> b)
{
  std::vector<T> x(b.size());
  int M = A.getNumRows();
  int N = A.getNumColumns();
  T** B;
  B = new T*[M];
  for(int i = 0; i < M; i++)
    B[i] = A.getRowPointer(i);
  RAPT::rsLinearAlgebra::rsSolveLinearSystemInPlace(B, &x[0], &b[0], N);
  delete[] B;
  return x;
}
// This is a temporary solution using the old Gaussian elimination code - todo: adapt that code 
// for the new rsMatrix class - use the new elementary row-operations - try to use as little extra 
// memory as possible - and if some is needed, use workspace parameters.
// It's a bit of a frankensteinization of new code (using rsMatrix as input) and old code (using 
// the old implementation of Gaussian elimination)


/*

ToDo:
-maybe provide an API that takes the flat arrays of the matrices directly as inputs
-implement solveOver/UnderDeterminedSystem: compute a least-squares solution or minimum-norm 
 solution respectively
 -based on that, implement pseudoInverse
-implement makeSystemDiagonal/solveDiagonalSystem. should call makeSystemUpperTriangular and then
 further reduce the triangular system to a diagonal one - if possible, this should also produce
 the diagonalizer - this can be used for diagonalization later..or wait - no - maybe not
-functions that ares still valid in the new framework (the band-diagonal and 2x2 stuff) shall
 eventually be moved over to here
-implement Gram-Schmidt orthogonalization of a set of basis vectors
-figure out and document the conditions under which we may use "solve" in place, i.e. the 3rd
 argument may equal the 2nd (E,E here) - i think, E must be upper triangular such that during the
 elimination only zeros get subtracted - but the algo may swap rows - so perhaps it works only if
 E is diagonal? in this case, whatever gets swapped, we will only ever add zeros to the rows


 -why is it that Gaussian elimination doesn't need to keep track of the swaps - or is this just 
  accidentally with our particular test examples? -> make unit tests with random matrices and 
  vectors. i think, it's because when you swap rows of the coefficient matrix and do the same
  swap in the rhs vector, we don't need to swap anything in the x-vector

maybe make a member function getInverse of rsMatrixView...hmm - but maybe not - this introduces 
too much circular dependencies between rsMatrixView/rsMatrix/rsLinearAlgebraNew - at the moment, 
the dependency is rsMatrixView < rsMatrix and rsMatrixView < rsLinearAlgebraNew (where < means
"right depends on left"). I actually want rsMatrix to remain independent of rsLinearAlgebraNew
because linear algebra is something *on top* of matrices...or maybe move basic, *simple* LA 
algorithms into the class rsMatrixView and use this class here only for more sophisticated stuff
...hmm - but no - i think, rsMatrix/View should be responsibel only for basic arithmetic, otoh
matrix "division" is an arithmetic operation - we'll see
...oh - actually LA depends on rsMatrix - not only the view - because we use it for temporaries
-to decouple it, we may use new/delete to allocate temporary memory and wrap rsMatrixView around
 the arrays - oh - but inverse returns a matrix and we really want it to be that way because 
 it's the most convenient way - so the dependency is perhaps inavoidable

*/