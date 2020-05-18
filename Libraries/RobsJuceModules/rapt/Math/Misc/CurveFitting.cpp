template<class T>
bool rsCurveFitter::fitExponentialSum(T* y, int numValues, T* A, T* a, int numExponentials)
{
  rsAssert( numExponentials == numValues/2 );
    // later, we may generalize this to allow numExponentials <= numValues/2
    // for a least squares fit

  int N = numValues;
  int k = numExponentials;
  bool result = true;
  T *C   = new T[k+1];                  // characteristic polynomial coefficients
  T *rhs = new T[k];                    // right hand sides
  T **M; rsArrayTools::allocateSquareArray2D(M, k); // matrices

  // find coefficients and roots of the characteristic polynomial:
  int i, j;
  for(i = 0; i < k; i++)
  {
    rhs[i] = -y[k+i];
    for(j = 0; j < k; j++)
      M[i][j] = y[i+j];
  }
  rsLinearAlgebra::rsSolveLinearSystem(M, C, rhs, k);
  C[k] = 1.0;
  //rsComplexDbl *roots = new rsComplexDbl[k];
  std::complex<T> *roots = new std::complex<T>[k];
  rsPolynomial<T>::roots(C, k, roots);

  // if all roots are real, compute exponents and weights:
  if( !rsAreAllValuesReal(roots, k, RS_EPS(T)) )
    result = false;
  else
  {
    for(i = 0; i < k; i++)
      a[i] = log(roots[i].real());
    for(i = 0; i < k; i++)
    {
      rhs[i] = y[i];
      for(j = 0; j < k; j++)
        M[i][j] = exp(a[j]*i); // == pow(roots[j].re, i)
    }
    rsLinearAlgebra::rsSolveLinearSystem(M, A, rhs, k);
  }

  // clean up and return result:
  delete[] C;
  delete[] rhs;
  rsArrayTools::deAllocateSquareArray2D(M, k);
  delete[] roots;
  return result;

  // Can the algorithm be generalized to fit complex exponentials to complex datapoints? Or maybe
  // it already can do this when being instantiated with a complex datatype? If yes, does that mean
  // we can estimate frequencies of sinusoids with it?

  // see Numerical Methods for Scientists and Engineers (Hamming), Ch.39, page 617f
}

template<class T>
std::vector<T> rsCurveFitter::multipleLinearRegression(const rsMatrix<T>& X, const T* y)
{
  rsMatrix<T>    XX  = X * X.getTranspose();    // XX  = X * X^T
  std::vector<T> rhs = X.productWith(y);        // rhs = X * y
  return rsLinearAlgebraNew::solveOld(XX, rhs); // solve (X * X^T) * b = X * y for b
}
// in the assignment: rsMatrix<T> XX = X * X.getTranspose(); 
// the product of X * X^T can be allocated and filled without explicitly creating the transposed
// matrix as temporary object - use some matrixMultiplySecondTransposed function (to be written)

template<class T>
rsPolynomial<T> rsCurveFitter::fitPolynomial(T* x, T* y, int numDataPoints, int degree)
{
  return rsPolynomial<T>(fitPolynomialStdVec(x, y, numDataPoints, degree));
}
// maybe it can be avoided to first create a vector and then wrapping it into a polynomial...but 
// maybe a move-constructor can be used to move the vector into the polynomial without additional
// memory allocation? that means, if a polynomial is created from a std::vector, we want to avoid
// memory allocation in the constructor of rsPolynomial and instead take over ownership of the 
// passed vector - if possible...

template<class T>
std::vector<T> rsCurveFitter::fitPolynomialStdVec(T* x, T* y, int numDataPoints, int degree)
{
  rsMatrix<T> X = polyFitDataMatrix(numDataPoints, x, degree);  // MxN data matrix X...
  return multipleLinearRegression(X, y);                        // ...used for the regressors
}

template<class T>
rsMatrix<T> rsCurveFitter::polyFitDataMatrix(int numDataPoints, T* x, int degree)
{
  int M = degree+1;       // # rows
  int N = numDataPoints;  // # cols
  rsMatrix<T> X(M, N);
  typedef rsArrayTools AT;
  AT::fillWithValue(X.getRowPointer(0), N, T(1));   // 1st row is all ones
  for(int i = 1; i < M; i++)                        // i-th row is (i-1)th row times x
    AT::multiply(X.getRowPointer(i-1), x, X.getRowPointer(i), N);
  return X;
}