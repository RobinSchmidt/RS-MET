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
void rsRationalFunction<T>::valueAndSlopeAt(const T& x, T* y, T* yp) const
{
  T n, np, d, dp;
  num.valueAndSlopeAt(x, &n, &np); // compute value n and derivative np of numerator
  den.valueAndSlopeAt(x, &d, &dp); // compute value d and derivative dp of denominator
  *y  = n / d;                     // compute and assign function value 
  *yp = (np*d-dp*n) / (d*d);       // compute and assign derivative via quotient rule

  // ToDo: 
  // -maybe precompute 1/d and replace the divisions by multiplications
}

//-------------------------------------------------------------------------------------------------
// Computations on std::vector
// Functions that operate on std::vectors to perform polynomial coefficient array manipulations,
// translated from my python implementation. They should be moved into rsRationalFunction as static
// meber functions. They are sort of low-level, although they use std::vector...maybe mid-level, but
// they may need to resize the vectors. maybe factor out true low-level functions (operating on raw 
// arrays) and if the need to resize, they don't actually resize anything but just inform the caller
// about the new size by a return value. ...maybe some of them should be moved into rsPolynomial

template<class T>
T rsRationalFunction<T>::polyEval(std::vector<T>& p, T x)
{
  int k = (int)p.size()-1;  // last valid index
  if(k < 0)
    return 0;
  T y = p[k];
  while(k > 0) {
    k -= 1;
    y = y*x + p[k]; }
  return y;
}

template<class T>
void rsRationalFunction<T>::polyTrunc(std::vector<T>& p, T tol)
{
  int i = (int)p.size();
  while(i > 1) {  // a polynomial should have at least 1 coeff
                  //if(fabs(p[i-1]) > tol)
    if( rsGreaterAbs(p[i-1], tol) )
      break;
    i -= 1; }
  p.resize(i);
}

template<class T>
T rsRationalFunction<T>::makeMonic(std::vector<T>& p)
{
  T lc = rsLast(p);
  for(size_t i = 0; i < p.size(); i++)
    p[i] /= lc;
  return lc;
}

template<class T>
std::vector<T> rsRationalFunction<T>::polyAdd(
  const std::vector<T>& p, const std::vector<T>& q, 
  T tol, T wp, T wq)
{
  int np = (int) p.size();
  int nq = (int) q.size();
  std::vector<T> r;
  if(np >= nq) {
    r = wp*p;
    for(int i = 0; i < nq; i++) {
      r[i] += wq*q[i];
      polyTrunc(r, tol); }}
  else {
    r = wq*q;
    for(int i = 0; i < np; i++) {
      r[i] += wp*p[i];
      polyTrunc(r, tol); }}
  return r;
}

template<class T>
std::vector<T> rsRationalFunction<T>::polySub(const std::vector<T>& p, const std::vector<T>& q,
  T tol)
{
  return polyAdd(p, q, 1, -1, tol);
}

template<class T>
std::vector<T> rsRationalFunction<T>::polyMul(const std::vector<T>& x, const std::vector<T>& h,
  T tol)
{
  int L = (int)x.size() + (int)h.size() - 1;  // length of result
  std::vector<T> y(L);
  for(int n = 0; n < L; n++) {
    y[n] = 0;  
    for(int k = std::max(0, n-(int)x.size()+1); k < std::min(n+1, (int)h.size()); k++)
      y[n] += h[k] * x[n-k]; }
  polyTrunc(y, tol);
  return y;
}

template<class T>
void rsRationalFunction<T>::polyDivMod(std::vector<T> p, std::vector<T> d, 
  std::vector<T>& q, std::vector<T>& r, T tol)
{ 
  q.resize(p.size());
  r = p;                 // init remainder with copy of product
  rsFill(q, T(0));       // init quotient to all zeros
  int k = (int)p.size() - (int)d.size();
  while(k >= 0) {
    q[k] = r[(int)d.size()-1+k] / d[(int)d.size()-1];
    int j = (int)d.size()+k-2; 
    while(j >= k) {
      r[j] -= q[k] * d[j-k];
      j -= 1; }
    k -= 1; }
  for(int i = (int) d.size()-1; i < (int) p.size(); i++)
    r[i] = 0;
  polyTrunc(q, tol);
  polyTrunc(r, tol);
}

template<class T>
std::vector<T> rsRationalFunction<T>::polyDiv(std::vector<T> p, std::vector<T> d, T tol)
{
  std::vector<T> q, r;
  polyDivMod(p, d, q, r, tol);
  return q;
}

template<class T>
std::vector<T> rsRationalFunction<T>::polyMod(std::vector<T> p, std::vector<T> d, T tol)
{
  std::vector<T> q, r;
  polyDivMod(p, d, q, r, tol);
  return r;
}

template<class T>
bool rsRationalFunction<T>::isAllZeros(const std::vector<T>& v, T tol)
{
  for(size_t i = 0; i < v.size(); i++)
    //if(fabs(v[i]) > tol)
    if( rsGreaterAbs(v[i], tol) )
      return false;
  return true;
}

template<class T>
std::vector<T> rsRationalFunction<T>::polyGCD(
  const std::vector<T>& p, const std::vector<T>& q, T tol, bool monic)
{
  std::vector<T> a = p, b = q, t;
  while(!isAllZeros(b, tol)) {
    t = b;
    b = polyMod(a, b, tol);
    a = t; }
  if(monic)
    makeMonic(a);
  return a;
}
// See: https://cp-algorithms.com/algebra/polynomial.html. It has an algorithm that is potentially
// faster - the "half-GCD-Algorithm"


template<class T>
std::vector<T> rsRationalFunction<T>::polyNest(const std::vector<T>& a, const std::vector<T>& b)
{
  int aN = (int)a.size()-1;               // degree of a
  int bN = (int)b.size()-1;               // degree of b
  int cN = aN*bN;                         // degree of result c
  std::vector<T> an(cN+1), c(cN+1); 
  rsFill(an, T(0)); an[0] = 1;            // powers of a, i.e. a^n - initially a^0 = [1 0 0...]
  rsFill(c,  T(0));                       // accumulator for result
  int K = 1;
  for(int n = 1; n <= bN; n++) {
    an = polyMul(an, a, 0.0);
    K += aN;
    for(int k = 0; k < K; k++)  
      c[k] += b[n] * an[k]; }
  return c;
}

template<class T>
void rsRationalFunction<T>::ratReduce(const std::vector<T>& pIn, const std::vector<T>& qIn,
  std::vector<T>& pOut, std::vector<T>& qOut, T tol)
{
  std::vector<T> gcd = polyGCD(pIn, qIn, tol);
  pOut = polyDiv(pIn, gcd, tol);
  qOut = polyDiv(qIn, gcd, tol);
}

template<class T>
void rsRationalFunction<T>::rsRationalFunction<T>::ratMul(
  const std::vector<T>& p, const std::vector<T>& q,
  const std::vector<T>& r, const std::vector<T>& s,
  std::vector<T>& u, std::vector<T>& v, T tol, bool reduced)
{
  u = polyMul(p, r, tol);
  v = polyMul(q, s, tol);
  if(reduced)
    ratReduce(u, v, u, v, tol);
}

template<class T>
void rsRationalFunction<T>::ratDiv(
  const std::vector<T>& p, const std::vector<T>& q,
  const std::vector<T>& r, const std::vector<T>& s,
  std::vector<T>& u, std::vector<T>& v, T tol, bool reduced)
{
  ratMul(p, q, s, r, u, v, tol, reduced); // r and s are swapped
}

template<class T>
void rsRationalFunction<T>::ratAdd(
  const std::vector<T>& n1, const std::vector<T>& d1,
  const std::vector<T>& n2, const std::vector<T>& d2,
  std::vector<T>& nr, std::vector<T>& dr, 
  T tol, T w1, T w2)
{
  std::vector<T> gcd, f1, f2, s1, s2;
  gcd = polyGCD(d1, d2, tol);
  f1 = polyDiv(d2, gcd, tol);
  f2 = polyDiv(d1, gcd, tol);
  dr = polyMul(f1, d1, tol);
  s1 = polyMul(f1, n1, tol);         // 1st summand in numerator of result
  s2 = polyMul(f2, n2, tol);         // 2nd summand
  nr = polyAdd(s1, s2, tol, w1, w2); // numerator of result
}

template<class T>
void rsRationalFunction<T>::ratPolyNest(
  const std::vector<T>& ni, const std::vector<T>& di,
  const std::vector<T>& po,
  std::vector<T>& nr, std::vector<T>& dr, T tol)
{
  std::vector<T> nt;
  nr = { po[0] };   // numerator of result
  dr = { T(1)  };   // denominator of result
  nt = ni;          // temporary numerator (for convolutive accumulation)
  for(size_t k = 1; k < po.size(); k++) {
    dr = polyMul(dr, di, tol);  
    nr = polyMul(nr, di, tol);
    nr = polyAdd(nr, nt, tol, T(1), po[k]);
    nt = polyMul(nt, ni, tol); }
}

template<class T>
void rsRationalFunction<T>::ratNest(
  const std::vector<T>& nI, const std::vector<T>& dI,
  const std::vector<T>& nO, const std::vector<T>& dO,
  std::vector<T>& nR, std::vector<T>& dR, T tol)
{
  std::vector<T> nU, dU, nL, dL;
  ratPolyNest(nI, dI, nO, nU,	dU, tol);  // compute upper num and den
  ratPolyNest(nI, dI, dO, nL,	dL, tol);  // compute lower num and den
  ratDiv(nU, dU, nL, dL, nR, dR, tol);
}
// It's not optimal to call ratPolyNest two times - inside this function, there are values that are
// calculated just the same in both calls, namely the successive powers of ni - but this is not 
// meant to be optimized, high performance code. Maybe in production code, this optimization should
// be done.























//-------------------------------------------------------------------------------------------------
// Computations on raw coefficient arrays


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
// then move this implementation into prototypes


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