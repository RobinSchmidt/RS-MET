//using namespace RSLib;

template<class T>
bool rsFitSumOfExponentials(T* y, int numValues, T* A, T* a, int numExponentials)
{
  rsAssert( numExponentials == numValues/2 );
    // later, we may generalize this to allow numExponentials <= numValues/2
    // for a least squares fit

  int N = numValues;
  int k = numExponentials;
  bool result = true;
  T *C   = new T[k+1];                  // characteristic polynomial coefficients
  T *rhs = new T[k];                    // right hand sides
  T **M; rsArray::allocateSquareArray2D(M, k); // matrices

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
  rsPolynomial<T>::findPolynomialRoots(C, k, roots);

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
  rsArray::deAllocateSquareArray2D(M, k);
  delete[] roots;
  return result;
}
