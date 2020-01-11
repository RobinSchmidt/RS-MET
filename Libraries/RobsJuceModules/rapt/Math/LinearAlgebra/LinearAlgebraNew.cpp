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
  bool invertible = makeTriangular(A, B);
  if(!invertible)
    return false;                                // matrix was found to be singular
  solveTriangular(A, X, B);
  return true;                                   // matrix was regular -> success!
}

template<class T>
bool rsLinearAlgebraNew::makeDiagonal(rsMatrixView<T>& A, rsMatrixView<T>& B)
{
  bool success = makeTriangular(A, B);
  if(!success)
    return false;
  int N = A.getNumRows();
  for(int i = N-1; i >= 0; i--) {
    for(int j = i-1; j >= 0; j--) {
      T s = -A(j, i) / A(i, i);
      A.addWeightedRowToOther(i, j, s);
      B.addWeightedRowToOther(i, j, s); }
  } 
  // here, also restrict the column-range as optimization (see makeSystemUpperTriangular)
  // A.addWeightedRowToOther(i, j, s, i, N-1)? or (i, j, s, j, N-1)? or (i, j, s, j+1, N-1)?
  // implement unit test and figure out

  return true;
}

template<class T>
bool rsLinearAlgebraNew::makeTriangular(rsMatrixView<T>& A, rsMatrixView<T>& B)
{
  rsAssert(A.isSquare()); // can we relax this?

  T tooSmall = 1.e-12;  // if pivot is less than that, the matrix is singular
                        // todo: use something based on RS_EPS(T)
  int N = A.getNumRows();
  for(int i = 0; i < N; i++) {
    int p = i; 
    T maxAbs = T(0);
    for(int j = i; j < N; j++) {       // search pivot row
      if(rsAbs(A(j, i)) > maxAbs) { 
        maxAbs = rsAbs(A(j, i)); 
        p = j; }}
    if(rsIsCloseTo(maxAbs, 0.0, tooSmall)) {
      rsError("Matrix (numerically) singular");
      return false; }
    if(p != i) {                       // turn pivot row into current row
      A.swapRows(i, p); 
      B.swapRows(i, p); } 
    for(int j = i+1; j < N; j++) {     // pivot row subtraction
      T w = -A(j, i) / A(i, i);        // weight
      A.addWeightedRowToOther(i, j, w, i, A.getNumColumns()-1); // start at i: avoid adding zeros
      B.addWeightedRowToOther(i, j, w); }}
  return true;
}

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
// maybe we should have a version that swaps when encountering a zero, maybe make an even simpler
// version without the test for zero - if client code knows that the matrix is in a form such that
// a zero pivot is never encountered, it may use this version (for optimization purposes)
// ...or replace the test by an assertion and make the function void


template<class T>
void rsLinearAlgebraNew::solveTriangular(
  rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B)
{
  int M = X.getNumColumns();  // number of required solution vectors
  int N = A.getNumRows();     // number of elements in each solution vector
  for(int k = 0; k < M; k++) {
    for(int i = N-1; i >= 0; i--) {
      T tmp = T(0);
      for(int j = i+1; j < N; j++)
        tmp += A(i, j) * X(j, k);
      X(i, k) = (B(i, k) - tmp) / A(i, i); }}
}


/*

ToDo:
-maybe provide an API that takes the flat arrays of the matrices directly as inputs
-implement solveOver/UnderDeterminedSystem: compute a least-squares solution or minimum-norm 
 solution respectively
 -based on that, implement pseudoInverse
-implement makeSystemDiagonal/solveDiagonalSystem. should call makeSystemUpperTriangular and then
 further reduce the triangular system to a diagonal one - if possible, this should also produce
 the diagonalizer - this can be used for diagonalization later
-functions that ares still valid in the new framework (the band-diagonal and 2x2 stuff) shall
 eventually be moved over to here
-implement Gram-Schmidt orthogonalization of a set of basis vectors
-figure out and document the conditions under which we may use "solve" in place, i.e. the 3rd
 argument may equal the 2nd (E,E here) - i think, E must be upper triangular such that during the
 elimination only zeros get subtracted - but the algo may swap rows - so perhaps it works only if
 E is diagonal? in this case, whatever gets swapped, we will only ever add zeros to the rows


 -why is it that Gaussian elimination doesn't need to keep track of the swaps - or is this just 
  accidentally he with out particular test examples? -> make unit tests with random matrices and 
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