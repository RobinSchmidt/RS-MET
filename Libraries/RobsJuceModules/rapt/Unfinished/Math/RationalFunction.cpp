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
  //while(rsAbs(p[i]) < tol && i > 0)
  while( rsLessAbs(p[i], tol) && i > 0 )
    i--;
  return i;
}
// maybe move to rsPolynomial


template<class T>
template<class R>
void rsRationalFunction<T>::partialFractionExpansionDistinctPoles(
  std::complex<R>* num, int numDeg, std::complex<R>* den, int denDeg,
  const std::complex<R>* poles, std::complex<R>* pfeCoeffs)
{
  typedef RAPT::rsPolynomial<R> PolyR;
  typedef RAPT::rsPolynomial<std::complex<R>> PolyC;
  std::complex<R> numVal, denVal;
  for(int i = 0; i < denDeg; i++) {  // denDeg == # poles == # pfeCoeffs
    numVal = PolyC::evaluate(poles[i], num, numDeg);
    denVal = PolyR::evaluateFromRootsOneLeftOut(poles[i], poles, denDeg, i);
    pfeCoeffs[i] = numVal/denVal;
  }
}
// as an alternative to evaluateFromRootsOneLeftOut, we could compute denVal as the derivative of 
// denominator - maybe try it and compare numerical precision of both ways....
// ...maybe get rid of the local variables numVal, denVal


template<class T>
template<class R>
void rsRationalFunction<T>::partialFractionExpansionMultiplePoles(
  const std::complex<R>* num, int numDeg, const std::complex<R>* den, int denDeg,
  const std::complex<R>* poles, const int* multiplicities, int numDistinctPoles,
  std::complex<R>* pfeCoeffs)
{
  // establish coefficient matrix:
  std::complex<R> **A; rsArrayTools::allocateSquareArray2D(A, denDeg);
  std::complex<R> *tmp = new std::complex<R>[denDeg+1]; // deflated denominator
  std::complex<R> remainder;                            // always zero
  for(int i = 0, k = 0; i < numDistinctPoles; i++) {
    rsArrayTools::copy(den, tmp, denDeg+1);
    for(int m = 0; m < multiplicities[i]; m++) {
      rsPolynomial<T>::divideByMonomialInPlace(tmp, denDeg-m, poles[i], &remainder);
      for(int j = 0; j < denDeg; j++)
        A[j][k] = tmp[j];
      k++;
    }
  }

  // solve the linear system using an appropriately zero-padded numerator as RHS:
  rsArrayTools::copy(num, tmp, numDeg+1);
  rsArrayTools::fillWithZeros(&tmp[numDeg+1], denDeg-(numDeg+1));
  rsLinearAlgebra::rsSolveLinearSystem(A, pfeCoeffs, tmp, denDeg);

  // clean up:
  rsArrayTools::deAllocateSquareArray2D(A, denDeg);
  delete[] tmp;
}
// todo: try to figure out an extended version of the cover-up method that is used for distinct 
// poles and use it as alternative algorithm... and/or use the residue method


template<class T>
template<class R>
void rsRationalFunction<T>::partialFractionExpansion(
  std::complex<R> *num, int numDeg, std::complex<R> *den, int denDeg,
  const std::complex<R> *poles, const int *multiplicities, int numDistinctPoles,
  std::complex<R> *pfeCoeffs, std::complex<R>* polyCoeffs)
{
  // make denominator monic:
  std::complex<R> s = T(1)/den[denDeg];
  rsArrayTools::scale(num, numDeg+1, s);
  rsArrayTools::scale(den, denDeg+1, s);

  // obtain polynomial ("FIR") part by polynomial division (maybe factor out):
  T tol = 1.e-12; // ad hoc - use something based on numeric_limits::epsilon
  if(numDeg >= denDeg) {
    rsAssert(polyCoeffs != nullptr, "function has a polynomial part"); 
    rsPolynomial<std::complex<T>>::divide(num, numDeg, den, denDeg, polyCoeffs, num);
    numDeg = actualDegree(num, numDeg, tol); // new degree of numerator
    // todo: maybe zero out the higher coeffs that are close to zero totally
  }
  else if(polyCoeffs != nullptr)
    rsArrayTools::fillWithZeros(polyCoeffs, denDeg+1);  // or should it be numDeg+1, does it matter?

  // sanity checks:
  rsAssert(numDeg < denDeg);
  rsAssert(rsArrayTools::sum(multiplicities, numDistinctPoles) == denDeg);

  // dispatch between all-poles-distinct or poles-with-multiplicities algorithm: 
  if(denDeg == numDistinctPoles)
    partialFractionExpansionDistinctPoles(num, numDeg, den, denDeg, poles, pfeCoeffs);
  else
    partialFractionExpansionMultiplePoles(
      num, numDeg, den, denDeg, poles, multiplicities, numDistinctPoles, pfeCoeffs);
}

template<class T>
template<class R>
std::vector<std::complex<R>> rsRationalFunction<T>::partialFractions(
  const std::vector<std::complex<R>>& numerator,
  const std::vector<std::complex<R>>& denominator,
  const std::vector<std::complex<R>>& poles)
{
  typedef std::vector<std::complex<R>> Vec;
  Vec num = numerator;   // local copies
  Vec den = denominator; 
  rsAssert(num.size() < den.size()); // function must be strictly proper
  Vec pfeCoeffs(den.size()-1);
  partialFractionExpansionDistinctPoles(
    &num[0], (int) num.size()-1, &den[0], (int) den.size()-1, &poles[0], &pfeCoeffs[0]);
  return pfeCoeffs;
}

template<class T>
template<class R>
std::vector<std::complex<R>> rsRationalFunction<T>::partialFractions(
  const std::vector<std::complex<R>>& numerator,
  const std::vector<std::complex<R>>& denominator,
  const std::vector<std::complex<R>>& poles,
  const std::vector<int>& muls)
{
  typedef std::vector<std::complex<R>> Vec;
  Vec num = numerator;   // local copies
  Vec den = denominator; 
  rsAssert(num.size() < den.size()); // function must be strictly proper
  Vec pfeCoeffs(den.size()-1);
  partialFractionExpansionMultiplePoles(
    &num[0], (int) num.size()-1, &den[0], (int) den.size()-1, 
    &poles[0], &muls[0], (int) poles.size(), &pfeCoeffs[0]);
  return pfeCoeffs;
}


// resources:
// https://en.wikipedia.org/wiki/Partial_fraction_decomposition
// https://ccrma.stanford.edu/~jos/filters/Partial_Fraction_Expansion.html
// https://ccrma.stanford.edu/~jos/filters/FIR_Part_PFE.html


/*

maybe implement a conversion to a Taylor series:
https://en.wikipedia.org/wiki/Rational_function#Taylor_series
that may be useful to approximate an IIR filter by an FIR filter or vice versa

for integrating, see here: ftp://ftp.cs.wisc.edu/pub/techreports/1970/TR91.pdf
page 9 in particluar - we could let the function return another rational function for the rational 
part and the alpha_i, b_i coeffs for the transcendental part

*/