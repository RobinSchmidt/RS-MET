// contains some unfinished code "scratches" and code under construction which may eventually be 
// moved to somewhere else once it works and does something useful

template<class T>
bool isRowZero(const rsMatrix<T> x, int rowIndex, T tol)
{
  for(int j = 0; j < x.getNumColumns(); ++j)
    if( rsAbs(x(rowIndex, j)) > tol )
      return false;
  return true;
}

template<class T>
bool areRowsZero(const rsMatrix<T> x, int startRow, int endRow, T tol)
{
  for(int i = startRow; i <= endRow; ++i)
    if(!isRowZero(x, i, tol))
      return false;
  return true;
}
// make members of rsMatrixView

/** Returns true, if the space spanned by the columns of x is within the span of the columns of B.
That means each column of x can be expressed as some linear combination of the columns of B. */
template<class T>
bool isInSpanOf(rsMatrix<T> B, rsMatrix<T> x, T tol)
{
  RAPT::rsLinearAlgebraNew::makeTriangular(B, x);
  int rankB = getRowEchelonRank(B, tol);
  return areRowsZero(x, rankB, x.getNumRows()-1, tol);
}

template<class T>
bool spanSameSpace(const rsMatrix<T>& A, const rsMatrix<T>& B, T tol)
{
  return isInSpanOf(A, B) && isInSpanOf(B, A);
}
// needs test