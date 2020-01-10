
template<class T>
std::vector<T> rsLinearAlgebraNew::solveLinearSystem(RAPT::rsMatrixView<T>& A, std::vector<T>& b)
{
  int N = (int) b.size();
  std::vector<T> x(N);
  rsMatrixView<T> vx(N, 1, &x[0]), vb(N, 1, &b[0]);
  solveLinearSystem(A, vx, vb);
  return x;
}
// make a version that operates on raw arrays

template<class T>
RAPT::rsMatrix<T> rsLinearAlgebraNew::inverse(const RAPT::rsMatrixView<T>& A)
{
  rsAssert(A.isSquare()); // relax later - compute pseudoinverse in non-square case
  int N = A.getNumRows();
  RAPT::rsMatrix<T> tmp(N, N, A.getDataPointerConst()), E(N, N);
  E.setToIdentity();
  solveLinearSystem(tmp, E, E); // why does it work to use E for both - because E is the identity?
  return E; 
}
// maybe make a member function of rsMatrixView...hmm - but maybe not - this introduces too much
// circular dependencies between rsMatrixView/rsMatrix/rsLinearAlgebraNew - at the moment, the
// dependency is rsMatrixView < rsMatrix and rsMatrixView < rsLinearAlgebraNew (where < means
// "right depends on left")


template<class T>
bool rsLinearAlgebraNew::solveLinearSystem(
  rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B)
{
  // check, if everything makes sense:
  rsAssert(X.getNumRows()    == A.getNumColumns());
  rsAssert(B.getNumRows()    == A.getNumColumns());
  rsAssert(X.getNumColumns() == B.getNumColumns());
  rsAssert(A.isSquare()); 
  // relax last requirement later - if not square compute approximate solution in overdetermined 
  // cases and minimum-norm solution in underdetermined cases

  bool invertible = makeSystemUpperTriangular(A, B);
  if(!invertible)
    return false;  // matrix was found to be singular
  solveUpperTriangularSystem(A, X, B);
  return true; 
}

template<class T>
bool rsLinearAlgebraNew::makeSystemDiagonal(rsMatrixView<T>& A, rsMatrixView<T>& B)
{
  //rsError("not yet implemented");
  bool success = makeSystemUpperTriangular(A, B);
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

  return true;
}

template<class T>
bool rsLinearAlgebraNew::makeSystemUpperTriangular(rsMatrixView<T>& A, rsMatrixView<T>& B)
{
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
      T s = -A(j, i) / A(i, i);        // scaler
      A.addWeightedRowToOther(i, j, s, 0, A.getNumColumns()-1);
      B.addWeightedRowToOther(i, j, s, 0, B.getNumColumns()-1); }}
  return true;
}
// why is it that the algo doesn't need to keep track of the swaps - or is this just accidentally 
// the with out particular test examples? -> make unit tests with random matrices and vectors

// addWeightedRowToOther should allow to specify minimum and maximum column-index - we don't need
// to subtract zeros from zeros in the columns that are already zero from a previous step - see old
// code - the k-loop starts at i, not at 0 - we should pass i and N-1:
// A.addWeightedRowToOther(i, j, s, i, N-1), B.add..

template<class T>
void rsLinearAlgebraNew::solveUpperTriangularSystem(
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
-implement solveOver/UnderDeterminedSystem: compute a least-squares solution or minimum-norm 
 solution respectively
 -based on that, implement pseudoInverse
-implement makeSystemDiagonal/solveDiagonalSystem. should call makeSystemUpperTriangular and then
 further reduce the triangular system to a diagonal one - if possible, this should also produce
 the diagonalizer - this can be used for diagonalization later
-functions that ares still valid in the new framework (the band-diagonal and 2x2 stuff) shall
  eventually be moved over to here
-Gram-Schmidt orthogonalization of a set of basis vectors

*/