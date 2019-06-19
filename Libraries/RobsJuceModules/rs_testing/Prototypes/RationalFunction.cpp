

template<class T>
void rsRationalFunction<T>::partialFractionExpansion(
  std::complex<T> *numerator,   int numeratorDegree,
  std::complex<T> *denominator, int denominatorDegree,
  std::complex<T> *poles, int *multiplicities, int numDistinctPoles,
  std::complex<T> *pfeCoeffs)
{
  // sanity check for inputs:
  rsAssert(numeratorDegree < denominatorDegree);
  rsAssert(rsArray::sum(multiplicities, numDistinctPoles) == denominatorDegree);
  // todo: make a function that handles the case numeratorDegree >= denominatorDegree
  // -> this must do a polynomial division to obtain a purely polynomial part and a stricly proper
  // rational function, the latter of which can then be fed into this routine

  // make denominator monic:
  rsArray::scale(numerator,   numeratorDegree+1,   T(1)/denominator[denominatorDegree]);
  rsArray::scale(denominator, denominatorDegree+1, T(1)/denominator[denominatorDegree]);
  // hmm - modifying the input arrays is no good idea - maybe use temporary memory - or require the
  // inputs to be monic - client code should deal with making it monic

  // establish coefficient matrix:
  std::complex<T> **A; rsArray::allocateSquareArray2D(A, denominatorDegree);
  std::complex<T> *tmp = new std::complex<T>[denominatorDegree+1]; // deflated denominator
  std::complex<T> remainder;                                   // always zero
  for(int i = 0, k = 0; i < numDistinctPoles; i++) {
    rsArray::copyBuffer(denominator, tmp, denominatorDegree+1);
    for(int m = 0; m < multiplicities[i]; m++) {
      rsPolynomial<T>::divideByMonomialInPlace(tmp, denominatorDegree-m, poles[i], &remainder);
      for(int j = 0; j < denominatorDegree; j++)
        A[j][k] = tmp[j];
      k++;
    }
  }

  // solve the linear system using an appropriately zero-padded numerator as RHS:
  rsArray::copyBuffer(numerator, tmp, numeratorDegree+1);
  rsArray::fillWithZeros(&tmp[numeratorDegree+1], denominatorDegree-(numeratorDegree+1));
  rsLinearAlgebra::rsSolveLinearSystem(A, pfeCoeffs, tmp, denominatorDegree);

  // clean up:
  rsArray::deAllocateSquareArray2D(A, denominatorDegree);
  delete[] tmp;
}
// use a more efficent algorithm if all poles are simple, (see also Experiments - there's something
// said about that)
// see also here:
// https://en.wikipedia.org/wiki/Partial_fraction_decomposition

template class rsRationalFunction<double>; // remove when moving code to library

/*

maybe implement a conversion to a Taylor series:
https://en.wikipedia.org/wiki/Rational_function#Taylor_series
that may be useful to approximate an IIR filter by an FIR filter or vice versa

*/