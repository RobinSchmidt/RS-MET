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
int actualDegree(std::complex<T>* p, int maxDegree, T tol)
{
  int i = maxDegree;
  while(rsAbs(p[i]) < tol && i > 0)
    i--;
  return i;
}
// maybe move to rsPolynomial


template<class T>
void rsRationalFunction<T>::partialFractionExpansionDistinctPoles(
  std::complex<T>* num, int numDeg, std::complex<T>* den, int denDeg,
  const std::complex<T>* poles, std::complex<T>* pfeCoeffs)
{
  typedef RAPT::rsPolynomial<T> PolyR;
  typedef RAPT::rsPolynomial<std::complex<T>> PolyC;
  std::complex<T> numVal, denVal;
  for(int i = 0; i < denDeg; i++) {  // denDeg == # poles == # pfeCoeffs
    numVal = PolyC::evaluate(poles[i], num, numDeg);
    denVal = PolyR::evaluateFromRootsOneLeftOut(poles[i], den, denDeg, i);
    pfeCoeffs[i] = numVal/denVal;
  }
}
// as an alternative to evaluateFromRootsOneLeftOut, we could compute denVal as the derivative of 
// denominator - maybe try it and compare numerical precision of both ways....
// ...maybe get rid of the local variables numVal, denVal


template<class T>
void rsRationalFunction<T>::partialFractionExpansionMultiplePoles(
  const std::complex<T>* num, int numDeg, const std::complex<T>* den, int denDeg,
  const std::complex<T>* poles, const int* multiplicities, int numDistinctPoles,
  std::complex<T>* pfeCoeffs)
{
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
// todo: try to figure out an extended version of the cover-up method that is used for distinct 
// poles and use it as alternative algorithm... and/or use the residue method


template<class T>
void rsRationalFunction<T>::partialFractionExpansion(
  std::complex<T> *num, int numDeg, std::complex<T> *den, int denDeg,
  const std::complex<T> *poles, const int *multiplicities, int numDistinctPoles,
  std::complex<T> *pfeCoeffs, std::complex<T>* polyCoeffs)
{
  // make denominator monic:
  std::complex<T> s = T(1)/den[denDeg];
  rsArray::scale(num, numDeg+1, s);
  rsArray::scale(den, denDeg+1, s);
  // hmm - modifying the input arrays is no good idea - maybe use temporary memory - or require the
  // inputs to be monic - client code should deal with making it monic
  // ..or well, actually, we potentially destroy the numerator array anyway due to the in-place
  // polynomial division above - i should probably write into the documentation that this function
  // may destroy the original content of the input arrays
  // maybe this should be done before the polynomial division?

  // obtain polynomial ("FIR") part by polynomial division, see:
  // https://ccrma.stanford.edu/~jos/filters/FIR_Part_PFE.html
  T tol = 1.e-12; // ad hoc - use something based on numeric_limits::epsilon
  if(numDeg >= denDeg) {
    rsAssert(polyCoeffs != nullptr, "function has a polynomial part"); 
    rsPolynomial<std::complex<T>>::divide(num, numDeg, den, denDeg, polyCoeffs, num);
    numDeg = actualDegree(num, numDeg, tol); // new degree of numerator
    // todo: maybe zero out the higher coeffs that are close to zero totally

  }
  else if(polyCoeffs != nullptr)
    rsArray::fillWithZeros(polyCoeffs, denDeg+1);  // or should it be numDeg+1, does it matter?
  // maybe factor out into function polynomialPart?


  // sanity check for inputs:
  rsAssert(numDeg < denDeg);
  rsAssert(rsArray::sum(multiplicities, numDistinctPoles) == denDeg);

  // if all poles are distinct, we can use a simpler (more efficient and hopefully also more
  // numerically precise) algorithm:
  if(denDeg == numDistinctPoles){

    partialFractionExpansionMultiplePoles(num, numDeg, den, denDeg,
      poles, multiplicities, numDistinctPoles, pfeCoeffs);
    // ...preliminary - until function below is tested...

    //partialFractionExpansionDistinctPoles(num, numDeg, den, denDeg, poles, pfeCoeffs);
  }
  else {
    partialFractionExpansionMultiplePoles(num, numDeg, den, denDeg,
      poles, multiplicities, numDistinctPoles, pfeCoeffs);
  }

  // https://ccrma.stanford.edu/~jos/filters/Partial_Fraction_Expansion.html
}




// use a more efficent algorithm if all poles are simple, (see also Experiments - there's something
// said about that)
// see also here:
// https://en.wikipedia.org/wiki/Partial_fraction_decomposition


/*

maybe implement a conversion to a Taylor series:
https://en.wikipedia.org/wiki/Rational_function#Taylor_series
that may be useful to approximate an IIR filter by an FIR filter or vice versa

for integrating, see here: ftp://ftp.cs.wisc.edu/pub/techreports/1970/TR91.pdf
page 9 in particluar - we could let the function return another rational function for the rational 
part and the alpha_i, b_i coeffs for the transcendental part

*/