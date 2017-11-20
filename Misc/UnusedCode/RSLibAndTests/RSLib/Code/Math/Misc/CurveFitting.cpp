using namespace RSLib;

bool RSLib::rsFitSumOfExponentials(double *y, int numValues, double *A, double *a,
                                   int numExponentials)
{
  rsAssert( numExponentials == numValues/2 );
    // later, we may generalize this to allow numExponentials <= numValues/2
    // for a least squares fit

  int N = numValues;
  int k = numExponentials;
  bool result = true;
  double *C   = new double[k+1];             // characteristic polynomial coefficients
  double *rhs = new double[k];               // right hand sides
  double **M; rsAllocateSquareArray2D(M, k); // matrices

  // find coefficients and roots of the characteristic polynomial:
  int i, j;
  for(i = 0; i < k; i++)
  {
    rhs[i] = -y[k+i];
    for(j = 0; j < k; j++)
      M[i][j] = y[i+j];
  }
  rsSolveLinearSystem(M, C, rhs, k);
  C[k] = 1.0;
  rsComplexDbl *roots = new rsComplexDbl[k];
  findPolynomialRoots(C, k, roots);

  // if all roots are real, compute exponents and weights:
  if( !rsAreAllValuesReal(roots, k, RS_EPS(double)) )
    result = false;
  else
  {
    for(i = 0; i < k; i++)
      a[i] = log(roots[i].re);
    for(i = 0; i < k; i++)
    {
      rhs[i] = y[i];
      for(j = 0; j < k; j++)
        M[i][j] = exp(a[j]*i); // == pow(roots[j].re, i)
    }
    rsSolveLinearSystem(M, A, rhs, k);
  }

  // clean up and return result:
  delete[] C;
  delete[] rhs;
  rsDeAllocateSquareArray2D(M, k);
  delete[] roots;
  return result;
}
