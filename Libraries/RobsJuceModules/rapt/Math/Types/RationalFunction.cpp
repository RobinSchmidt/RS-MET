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


  std::complex<R> tol = 1.e-12;
  // TODO: use something based on numeric_limits::epsilon or let the user pass it in


  rsLinearAlgebra::rsSolveLinearSystem(A, pfeCoeffs, tmp, denDeg, tol);

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
  // Make denominator monic:
  std::complex<R> s = T(1)/den[denDeg];
  rsArrayTools::scale(num, numDeg+1, s);
  rsArrayTools::scale(den, denDeg+1, s);
  // TODO: factor out - or maybe try to get rid and declare num and den const

  // Obtain polynomial ("FIR") part by polynomial division (maybe factor out):

  T tol = 1.e-12; 
  // TODO: use something based on numeric_limits::epsilon or let the user pass it in

  if(numDeg >= denDeg) {
    rsAssert(polyCoeffs != nullptr, "function has a polynomial part"); 
    rsPolynomial<std::complex<T>>::divide(num, numDeg, den, denDeg, polyCoeffs, num);
    numDeg = actualDegree(num, numDeg, tol); // new degree of numerator
    // todo: maybe zero out the higher coeffs that are close to zero totally
  }
  else if(polyCoeffs != nullptr)
    rsArrayTools::fillWithZeros(polyCoeffs, denDeg+1);  // or should it be numDeg+1, does it matter?

  // Sanity checks:
  rsAssert(numDeg < denDeg);
  rsAssert(rsArrayTools::sum(multiplicities, numDistinctPoles) == denDeg);

  // Dispatch between all-poles-distinct or poles-with-multiplicities algorithm: 
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


//=================================================================================================
/*

ToDo:

- Implement evaluation with derivative. Use quotient rule  (u/v)' = (u'*v - v'*u) / v^2. The 
  polynomial class has a function to evaluate the polynomial with derivative. Use that. Oh - we 
  already have valueAndSlopeAt. 

- Derive formula for (u/v)'' in terms of u,u',u'',v,v',v'' and implement that, too. Use quotient 
  rule on ((u'*v - v'*u) / v^2)'   v^2' = 2*v*v', so I think, we get:
     ((u'*v - v'*u)' * v^2 - (u'*v - v'*u)*2*v*v) / (v^4) 
  =  ( ( (u''*v + u'*v')  - (v''*u + v'*u') )' * v^2 - (u'*v - v'*u)*2*v*v) / (v^4) 
  ...Verify this! Check the formulas numerically by comparing to numerical differentiation result.
  Add a function valueSlopeAndCurvature or valueAndDerivatives2. Look into rsPolynomial (it has 
  such functions) and use matching naming conventions and API. I think, there are generalizations 
  of the quotient rule for higher derivatives. Look up the .tex files for the math book. There may
  be something about that in there, IIRC.

- For integrating, see here: ftp://ftp.cs.wisc.edu/pub/techreports/1970/TR91.pdf  page 9 in 
  particluar. We could let the function return another rational function for the rational part and 
  the alpha_i, b_i coeffs for the transcendental part

- Add a function evaluateInverse(y, xGuess). It should use Newton or Halley iteration (using the 
  evaluation with derivative(s) function) to find an x such that r(x) = y. The rational function 
  r(x) may not be invertible, so the user should provide an initial guess for x and the function 
  will produce a value near the given x

- Maybe implement a conversion to a Taylor series:
  https://en.wikipedia.org/wiki/Rational_function#Taylor_series
  that may be useful to approximate an IIR filter by an FIR filter or vice versa. See also the
  functions in the Prototypes.cpp to convert from Taylor to Pade approximation.

- Implement a "less-than" comparison function according to page 16 in "Counterexamples in 
  Analysis". The book defines an "ordered field" to a be field in which there exists a subset P of
  the underlying set F, called the "positive elements", for which the following holds:
  (1) x in P and y in P  implies  x + y in P, (2) x in P and y in P  implies  x * y in P, (3) for 
  any x in F exactly one of the 3 statements is true: x in P, x = 0, -x in P. The less-than 
  relation is then defined as: x < y  iff  (y - x) in P. So: x >= y  iff  (x - y) in P  or  x = y.
  For rational functions, the set P of positive elements is defined to be those, whose leading 
  coefficients of numerator and denominator have the same sign. So, what we would have to do to 
  implement the < relation of two rational functions R,S is  (1) compute D = S - R,  (2) check if 
  leading coeffs of numerator and denominator of D have the same sign. Beware of roundoff error 
  issues, though. Also, the idea works only for rational functions with integer, rational or 
  real coefficients but not for rational functions with complex coefficients, for example. I think,
  the arguments may still be complex though - but the coeffs must be real (verify!). Maybe the 
  underlying set of coefficients (which can be a field or ring, I guess?) must itself be ordered 
  for the idea to work? Figure out! The same idea can be used for polynomials, too because they are
  just a special case of rational functions (those with denominator 1). Based on the < relation, we
  may also implement an absloute value function: |f| = f  iff  f >= 0  and  f = -f  iff  f < 0  
  where 0 is the zero rational function (represented as  0*x^0 / 1*x^0  here, i.e. the zero 
  polynomial over the one polynomial...I think - verify!). All these things may be needed when we
  want to do linear algebra with matrices of rational functions because the pivoting steps in 
  Gaussian elemination need these operations. And we do indeed need matrix inversion of matrices of
  rational functions when we want to compute the transfer function of a state space filter 
  symbolically - so this stuff may actually be relevant in practice and not just an academic 
  excercise. ...well - I don't think, pivot selection should be based on that definition of an
  absolute value. Instead, we should use an appropriate rsIsBetterPivot() function. This is 
  currently in the works anyway.

- Try to implement an algorithm that decomposes rational functions into continued fraction 
  expansions. See:
  https://www.fim.uni-passau.de/fileadmin/dokumente/fakultaeten/fim/lehrstuhl/sauer/geyer/Kettenbrueche.pdf
  Maybe that could be useful to express filters in a nested structure of low order filters (e.g.
  biquads?). Maybe the z^-1 factors (i.e. the delays) would then be replaced by more general 
  filters?

- Implement an algorithm that decomposes an arbitrary filter transfer function H(z) into an
  FIR part plus (or maybe times?) a series connection of a minimum phase IIR filter and an allpass.
  Decomposing a proper (i.e. #zeros <= #poles) IIR transfer function into a minimum phase part and
  allpass could perhaps be achieved by factoring numerator and denominator, then finding the zeros
  outside of the unit circle and factor them out. Maybe they are already combined with 
  corresponding poles inside the unit circle - in this case, we would factor out a pair of zeros 
  and poles at a time. For zeros outside the unit circle that do not have a partner pole (i.e. are
  not part of an allpass), maybe we need to introduce one "ghost-pole" in the allpass (to keep it 
  allpass) and then perhaps compensate for that by introducing a cancelling zero into the min-phase 
  IIR part or into the FIR part. Probably the latter - otherwise the min-phase IIR would not be a 
  proper rational function (it may then have more zeros than poles - which isn't necessarily a 
  problem though) ..not sure about that, though.

- I think, rational functions R with rational coefficients can always be re-expressed as rational
  functions with integer coeffs by doing the following: Find the lcm of all the denominators of
  the coeffs (both numerator and denominator coeffs), then multiply numerator and denominator of R
  by that lcm. Now every coeff should be an integer, i.e. the denominator of each coeff should be 
  1. It now could still be the case that these integer coeffs have common factors, so it may make 
  sense to compute the gcd of all coeffs (again considering both numerator and denominator of the 
  function R) and then divide numerator and denominator by that gcd.

*/