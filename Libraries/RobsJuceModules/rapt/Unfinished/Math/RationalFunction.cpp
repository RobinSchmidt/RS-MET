template<class T>
bool rsRationalFunction<T>::reduce(T tol)
{
  std::vector<T> gcd = polyGCD(num.coeffs, den.coeffs, tol);
  if(gcd.size() == 1)
    return false;  // was already in reduced form
  num.coeffs = polyDiv(num.coeffs, gcd, tol);
  den.coeffs = polyDiv(den.coeffs, gcd, tol);
  return true;
}

template<class T>
void rsRationalFunction<T>::partialFractionExpansion(
  std::complex<T> *num, int numDeg, std::complex<T> *den, int denDeg,
  std::complex<T> *poles, int *multiplicities, int numDistinctPoles,
  std::complex<T> *pfeCoeffs, std::complex<T>* polyCoeffs)
{
  // todo: if(numDeg >= denDeg), do a polynomial division to obtain the polynomial part

  // sanity check for inputs:
  rsAssert(numDeg < denDeg);
  rsAssert(rsArray::sum(multiplicities, numDistinctPoles) == denDeg);
  // todo: make a function that handles the case numeratorDegree >= denominatorDegree
  // -> this must do a polynomial division to obtain a purely polynomial part and a stricly proper
  // rational function, the latter of which can then be fed into this routine

  // make denominator monic:
  rsArray::scale(num, numDeg+1, T(1)/den[denDeg]); //!! wrong? should use numerator?
  rsArray::scale(den, denDeg+1, T(1)/den[denDeg]);
  // hmm - modifying the input arrays is no good idea - maybe use temporary memory - or require the
  // inputs to be monic - client code should deal with making it monic

  // establish coefficient matrix:
  std::complex<T> **A; rsArray::allocateSquareArray2D(A, denDeg);
  std::complex<T> *tmp = new std::complex<T>[denDeg+1]; // deflated denominator
  std::complex<T> remainder;                            // always zero
  for(int i = 0, k = 0; i < numDistinctPoles; i++) {
    rsArray::copyBuffer(den, tmp, denDeg+1);
    for(int m = 0; m < multiplicities[i]; m++) {
      rsPolynomial<T>::divideByMonomialInPlace(tmp, denDeg-m, poles[i], &remainder);
      for(int j = 0; j < denDeg; j++)
        A[j][k] = tmp[j];
      k++;
    }
  }

  // solve the linear system using an appropriately zero-padded numerator as RHS:
  rsArray::copyBuffer(num, tmp, numDeg+1);
  rsArray::fillWithZeros(&tmp[numDeg+1], denDeg-(numDeg+1));
  rsLinearAlgebra::rsSolveLinearSystem(A, pfeCoeffs, tmp, denDeg);

  // clean up:
  rsArray::deAllocateSquareArray2D(A, denDeg);
  delete[] tmp;
}
// use a more efficent algorithm if all poles are simple, (see also Experiments - there's something
// said about that)
// see also here:
// https://en.wikipedia.org/wiki/Partial_fraction_decomposition


/*

maybe implement a conversion to a Taylor series:
https://en.wikipedia.org/wiki/Rational_function#Taylor_series
that may be useful to approximate an IIR filter by an FIR filter or vice versa

*/