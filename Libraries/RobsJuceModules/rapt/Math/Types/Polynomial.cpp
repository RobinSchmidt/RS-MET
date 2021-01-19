
// construction:

template<class T>
rsPolynomial<T>::rsPolynomial(int degree, bool initWithZeros)
{
  coeffs.resize(degree+1);
  if(initWithZeros)
    rsArrayTools::fillWithZeros(&coeffs[0], degree+1);
}


// operators:
/*
template<class T>
rsPolynomial<T> operator+(const rsPolynomial<T>& p, const rsPolynomial<T>& q)
{
  rsPolynomial<T> r(rsMax(getDegree(), q.getDegree()), false);
  weightedSum(coeffs.data(), getDegree(), T(1),
    q.coeffs.data(), q.getDegree(), T(1),
    r.coeffs.data());
  return r;
}
*/


// setup:

template<class T>
void rsPolynomial<T>::setCoeffs(const T* newCoeffs, int newDegree)
{
  coeffs.resize(newDegree+1);
  rsArrayTools::copy(newCoeffs, &coeffs[0], newDegree+1);
}

template<class T>
void rsPolynomial<T>::setRoots(const T* newRoots, int numRoots, T scaler)
{
  coeffs.resize(numRoots+1);
  rootsToCoeffs(newRoots, &coeffs[0], numRoots, scaler);
}

template<class T>
void rsPolynomial<T>::truncateTrailingZeros(const T& thresh)
{
  int i = (int) coeffs.size()-1;
  while(i > 0) {
    //if( rsAbs(coeffs[i]) > thresh ) // must be strictly greater to work correctly with thresh == 0
    if( rsGreaterAbs(coeffs[i], thresh) ) // must be strictly greater to work correctly with thresh == 0
      break;
    i--;
  }
  coeffs.resize(i+1);
}







//=================================================================================================
// raw number crunching routines operating on coefficient arrays:

//-------------------------------------------------------------------------------------------------
// evaluation:

template <class T>
T rsPolynomial<T>::evaluate(const T& x, const T *a, int degree)
{
  if(degree < 0)
    return T(0);
  T y = a[degree];
  for(int i = degree-1; i >= 0; i--)
    y = y*x + a[i];
  return y;
}

template<class T>
std::complex<T> rsPolynomial<T>::evaluateFromRoots(const std::complex<T>& s,
  const std::complex<T>* r, int N)
{
  std::complex<T> result = T(1);
  for(int i = 0; i < N; i++) {
    if(!isInfinite(r[i]))
      result *= (s - r[i]);
  }
  return result;
}
// actually, the roots don't have to be complex - the function should work with real roots just
// as well

template<class T>
std::complex<T> rsPolynomial<T>::evaluateFromRootsOneLeftOut(
  const std::complex<T>& x, const std::complex<T>* r, int nr, int i)
{
  return evaluateFromRoots(x, r, i) * evaluateFromRoots(x, &r[i+1], nr-i-1);
}

template<class T>
T rsPolynomial<T>::evaluateDerivative(const T& x, const T* a, int N)
{
  T y = T(N) * a[N];
  for(int i = N-1; i >= 1; i--)
    y = y*x + T(i) * a[i];
  return y;
}

template<class T>
T rsPolynomial<T>::evaluateDerivative(const T& x, const T* a, int N, int n)
{
  rsAssert(n >= 0, "derivative order must be non-negative");
  if(n > N)
    return T(0); // avoid evaluating products including zero (starting at negative indices)
  T y = T(rsProduct(N-n+1, N)) * a[N];
  for(int i = N-1; i >= n; i--)
    y = y*x + T(rsProduct(i-n+1, i)) * a[i];
  return y;
}
// runtime: (N-n)*n? ..the number of terms in each product is n, the i-loop runs N-n times
// we could probably also init y to 0 and run the loop from N to n - maybe the code would be neater
// but we would add a 0*x term - so we would have have one extra multiplication and one extra
// addition that doesn't really do anything useful

template <class T>
void rsPolynomial<T>::evaluateWithDerivative(const T& x, const T *a, int degree, T *y, T *yd)
{
  *y  = a[degree];
  *yd = 0.0;
  for(int i = degree-1; i >= 0; i--) {
    *yd = *yd * x + *y;
    *y  = *y  * x + a[i];
  }
}

template <class T>
void rsPolynomial<T>::evaluateWithDerivatives(const T& x, const T *a, int degree, T *results,
  int numDerivatives)
{
  rsAssert(numDerivatives < 32, "numDerivatives must be < 32"); // rsFactorials has 32 entries
  results[0] = a[degree];
  rsArrayTools::fillWithZeros(&results[1], numDerivatives);
  for(int i = degree-1; i >= 0; i--) {
    int n = rsMin(numDerivatives, degree-1);
    for(int j = n; j >= 1; j--)
      results[j] = results[j]*x + results[j-1];
    results[0] = results[0]*x + a[i];
  }
  rsArrayTools::multiply(&results[2], &rsFactorials[2], &results[2], numDerivatives-1);
  // todo: maybe lift the restriction to < 32 derivatives by computing additional factorials on 
  // the fly, if needed - but we need to be careful about overflow - i think 31! will already 
  // overflow int64...yep...easily
}

template<class T>
template<class R>
R rsPolynomial<T>::evaluateWithTwoDerivativesAndError(
  const std::complex<R>* a, int degree, std::complex<R> z, std::complex<R>* P)
{
  P[0] = a[degree];                    // P(z)
  P[1] = std::complex<R>(R(0), R(0));  // P'(z)
  P[2] = std::complex<R>(R(0), R(0));  // P''(z)
  R err = rsAbs(P[0]);                 // estimated roundoff error in evaluation of the polynomial
  R zA  = rsAbs(z);                    // absolute value of z
  for(int j = degree-1; j >= 0; j--) {
    P[2] = z * P[2] + P[1];
    P[1] = z * P[1] + P[0];
    P[0] = z * P[0] + a[j];
    err  = abs(P[0]) + zA*err; }
  P[2] *= R(2);
  return err;
}

template<class T>
T rsPolynomial<T>::evaluateIntegral(const T& x, const T* a, int N, T c)
{
  T y = a[N] / T(N+1);
  for(int i = N-1; i >= 0; i--)
    y = y*x + a[i] / T(i+1);
  return y*x + c;
}
// maybe make a variant that takes lower and upper integration limits a,b - this can be optimized:
// the final + c can be thrown away and more importantly, the evaluation at a and b can be done 
// within a single loop where a[i] / T(i+1) needs to be computed only once - in the loop, do:
//   Ai = a[i] / T(i+1); ya = ya*x + Ai; yb = yb*x + Ai; 
// where ya, yb are both initialized like y above. then return yb - ya

template<class T>
T rsPolynomial<T>::evaluateHermite(const T& x, int n)
{
  if(n == 0) return T(1);
  if(n == 1) return T(2)*x;
  T h0 = T(1);
  T h1 = T(2)*x;
  for(int i = 1; i < n; i++) {
    T tmp = T(2) * (x*h1 - T(i)*h0); // H[n+1](x) = 2*x*H[n](x) - 2*n*H[n-1](x)
    h0 = h1;
    h1 = tmp;
  }
  return h1;
}
// this is the physicist's version of Hermite polynomials, the probabilist's version would use the
// recursion: tmp = x*h1 - i*h0. ...without the factor two - maybe make the recursion factor 
// another parameter

template<class T>
T rsPolynomial<T>::evaluateNewton(const T& x, const T* c, const T* r, int N)
{
  T p = T(1);   // accumulator for Newton basis polynomials
  T y = c[0];   // accumulator for result
  for(int i = 0; i < N; i++) {
    p *= x - r[i];
    y += c[i+1] * p; }
  return y;
}

/*
template<class T>
int rsPolynomial<T>::actualDegree(T* p, int maxDegree, T tol)
{
  int i = maxDegree;
  while(rsAbs(p[i]) <= tol && i > 0)
    i--;
  return i;
}
*/
//-------------------------------------------------------------------------------------------------
// arithmetic:

template<class T>
void rsPolynomial<T>::weightedSum(
  const T* p, int pN, const T& wp, const T* q, int qN, const T& wq, T* r)
{
  int i;
  if(pN >= qN) {
    for(i = 0; i <= qN; i++)
      r[i] = wp*p[i] + wq*q[i];
    for(i = qN+1; i <= pN; i++)
      r[i] = wp*p[i];
  }
  else {
    for(i = 0; i <= pN; i++)
      r[i] = wp*p[i] + wq*q[i];
    for(i = pN+1; i <= qN; i++)
      r[i] = wq*q[i];
  }
}

template <class T>
void rsPolynomial<T>::divide(const T *p, int pDegree, const T *d, int dDegree, T *q, T *r)
{
  rsArrayTools::copy(p, r, pDegree+1); // init remainder with p
  rsArrayTools::fillWithZeros(q, pDegree+1); // init quotient with zeros
  for(int k = pDegree-dDegree; k >= 0; k--) {
    q[k] = r[dDegree+k] / d[dDegree];
    for(int j = dDegree+k-1; j >= k; j--)
      r[j] -= q[k] * d[j-k];
  }
  rsArrayTools::fillWithZeros(&r[dDegree], pDegree-dDegree+1);
  // maybe return the degree of the quotient, the degree of the remainder is then
  // pDegree-qDegree - ...what if dDegree > pDegree? the most sensibe thing in this case would be,
  // if the remainder is equal to p ...and the quotient should be 0 ...i think, the model is:
  // p(x) = q(x)*d(x) + r(x)
}

template<class T>
template<class S>
void rsPolynomial<T>::divideByMonomialInPlace(S *dividendAndResult, int dividendDegree,
  S x0, S *remainder)
{
  *remainder = dividendAndResult[dividendDegree];
  dividendAndResult[dividendDegree] = T(0);
  for(int i=dividendDegree-1; i>=0; i--) {
    S swap               = dividendAndResult[i];
    dividendAndResult[i] = *remainder;
    *remainder           = swap + *remainder*x0;
  }
}

/*
template <class T>
void rsPolynomial<T>::dividePolynomialByMonomialInPlace(T *dividendAndResult, int dividendDegree,
  T x0, T *remainder)
{
  *remainder = dividendAndResult[dividendDegree];
  dividendAndResult[dividendDegree] = T(0);
  for(int i=dividendDegree-1; i>=0; i--)
  {
    T swap               = dividendAndResult[i];
    dividendAndResult[i] = *remainder;
    *remainder           = swap + *remainder*x0;
  }
}
*/

template <class T>
void rsPolynomial<T>::greatestCommonDivisor(
  const T* p, int pDeg, const T* q, int qDeg, T* gcd, int* gcdDeg, T tol)
{
  rsError("Not yet implemented");
}

template <class T>
void rsPolynomial<T>::powers(const T* a, int N, T** aPowers, int highestPower)
{
  aPowers[0][0] = 1;
  if(highestPower < 1)
    return;
  rsArrayTools::copy(a, aPowers[1], N+1);
  for(int k = 2; k <= highestPower; k++)
    rsArrayTools::convolve(aPowers[k-1], (k-1)*N+1, a, N+1, aPowers[k]);
}

template <class T>
void rsPolynomial<T>::powers(const T* a, int N, T* aPowers, int highestPower, int stride)
{
  //rsError("Not yet tested");
  rsAssert(stride >= (N+1)*highestPower-1); // is this correct? or N*highestPower + 1

  aPowers[0] = 1;
  if(highestPower < 1)
    return;
  rsArrayTools::copy(a, &aPowers[stride], N+1);
  for(int k = 2; k <= highestPower; k++)
    rsArrayTools::convolve(&aPowers[(k-1)*stride], (k-1)*N+1, a, N+1, &aPowers[k*stride]);
}
// maybe make a version that assumes a different memory layout with "rows" having lengths:
// 1, N+1, 2*N+1, 3*N+1, ... i.e. k*N+1 such that each row es exactly long enough for the number
// of coeffs it holds - the way it's implemented now is convenient, especially when used together
// with rsMatrix, but wastes memory

template <class T>
void rsPolynomial<T>::compose(const T* a, int aN, const T* b, int bN, T* c, T* workspace)
{
  rsAssert(c != a && c != b, "Does not work in place");
  int cN = aN*bN;     // degree of c
  T*  an = workspace; // array for the successive powers of a[]
  an[0]  = T(1);      // initialize to a[]^0

  // accumulation:
  rsArrayTools::fillWithZeros(c, cN+1);
  c[0] = b[0];
  int K = 1;
  for(int n = 1; n <= bN; n++) {
    rsArrayTools::convolveInPlace(an, K, a, aN+1);
    K += aN;
    for(int k = 0; k < K; k++)
      c[k] += b[n] * an[k]; }
}

template <class T>
void rsPolynomial<T>::compose(const T* a, int aN, const T* b, int bN, T* c)
{
  int cN = aN*bN;               // degree of c
  T* an = new T[cN+1];          // allocate array for the successive powers of a[]
  compose(a, aN, b, bN, c, an); // call workspace based function
  delete[] an;                  // clean up
}

template<class T>
void rsPolynomial<T>::composeLinearWithCubic(T* a, T* c, T b0, T b1)
{
  T b02 = b0*b0;
  T b12 = b1*b1;
  c[0]  = a[3]*b0*b02 + a[2]*b02 + a[1]*b0 + a[0];
  c[1]  = T(3)*a[3]*b02*b1 + T(2)*a[2]*b0*b1 + a[1]*b1;
  c[2]  = T(3)*a[3]*b0*b12 + a[2]*b12;
  c[3]  = a[3]*b1*b12;
}
// We can compute the coeffs of the nested polynomial easily with sage:
//   var("a0 a1 a2 a3 b0 b1 c0 c1 c2 c3")
//   a(x) = a0 + a1*x + a2*x^2 + a3*x^3   # outer polynomial
//   b(x) = b0 + b1*x                     # inner polynomial
//   c(x) = a(b(x))                       # composed polynomial
//   (expand(c)).collect(x)
// which gives:
//   a3*b1^3*x^3 + a3*b0^3 + a2*b0^2 + (3*a3*b0*b1^2 + a2*b1^2)*x^2 + a1*b0 
//   + (3*a3*b0^2*b1 + 2*a2*b0*b1 + a1*b1)*x + a0
// so:
//   c0 = a3*b0^3 + a2*b0^2 + a1*b0 + a0
//   c1 = 3*a3*b0^2*b1 + 2*a2*b0*b1 + a1*b1
//   c2 = 3*a3*b0*b1^2 + a2*b1^2
//   c3 = a3*b1^3

template <class T>
void rsPolynomial<T>::negateArgument(const T *a, T *am, int N)
{
  scaleArgument(a, am, N, T(-1));
}

template <class T>
void rsPolynomial<T>::scaleArgument(const T* a, T* as, int N, T scaler)
{
  T s = T(1);
  for(int n = 0; n <= N; n++) {
    as[n] = s*a[n]; 
    s    *= scaler;  }
}

template <class T>
void rsPolynomial<T>::shiftArgument(const T *a, T *as, int N, T x0)
{
  T r[2] = { -x0, T(1) };
  compose(r, 1, a, N, as);
  return;
  // todo: provide workspace based version that uses the non-allocating version of compose


  // inefficient old implementation - todo: move to prototypes, it's worth to keep the code 
  // somewhere because it shows the algorithm derived from the binomial theorem:
  rsUint32 Nu = rsUint32(N); // used to fix warnings
  rsUint32 numLines = N+1;
  rsUint32 length   = (numLines*(numLines+1))/2;
  rsUint32 *pt = new rsUint32[length];
  rsCreatePascalTriangle(pt, numLines);
  T *x0n = new T[N+1];  // +- x0^n, the alternation comes from the minus in (x-x0)^n
  x0n[0] = 1.0;
  for(rsUint32 n = 1; n <= Nu; n++)
    x0n[n] = -x0 * x0n[n-1];
  for(rsUint32 n = 0; n <= Nu; n++) {
    as[n] = 0.0;
    for(rsUint32 k = n; k <= Nu; k++)
      as[n] += T(rsPascalTriangle(pt, k, k-n)) * x0n[k-n] * a[k]; }
  delete[] pt;
  delete[] x0n;
}

//-------------------------------------------------------------------------------------------------
// calculus:

template <class T>
void rsPolynomial<T>::derivative(const T *a, T *ad, int N)
{
  for(int n = 1; n <= N; n++)
    ad[n-1] = T(n) * a[n];
}

template <class T>
void rsPolynomial<T>::integral(const T *a, T *ai, int N, T c)
{
  for(int n = N+1; n >= 1; n--)
    ai[n] = a[n-1] / T(n);
  ai[0] = c;
}

template<class T>
void rsPolynomial<T>::integrateWithPolynomialLimits(
  const T* p, int pN, const T* a, int aN, const T* b, int bN, T* q)
{
  int PN = pN+1;  T* P = new T[PN+1];
  int AN = aN*PN; T* A = new T[AN+1];
  int BN = bN*PN; T* B = new T[BN+1];

  integral(p, P, pN);        // P(x) is the antiderivative of p(x)
  compose(a, aN, P, PN, A);  // A(x) = P(a(x))
  compose(b, bN, P, PN, B);  // B(x) = P(b(x))
  subtract(B, BN, A, AN, q); // q(x) = B(x) - A(x)

  delete[] P;
  delete[] A;
  delete[] B;
}
// todo: make a workspace based version, keep this as convenience function (or move it to protoypes
// and provide a convenience function that needs only 1 allocation)

template <class T>
void rsPolynomial<T>::finiteDifference(const T *a, T *ad, int N, int direction, T h)
{
  // (possibly alternating) powers of the stepsize h:
  T *hk = new T[N+1];
  T hs  = T(direction)*h;
  hk[0] = T(1);
  for(int k = 1; k <= N; k++)
    hk[k] = hk[k-1] * hs;

  // binomial coefficients:
  rsUint32 numCoeffs    = N+1;   // maybe use int
  rsUint32 triangleSize = (numCoeffs*(numCoeffs+1))/2;
  rsUint32 *binomCoeffs = new rsUint32[triangleSize];
  rsCreatePascalTriangle(binomCoeffs, numCoeffs);

  // actual coefficient computation for ad:
  rsArrayTools::fillWithZeros(ad, N);
  for(unsigned int n = 0; n <= (rsUint32)N; n++)
  {
    for(unsigned int k = 1; k <= n; k++)
      ad[n-k] += a[n] * T(rsPascalTriangle(binomCoeffs, n, k)) * hk[k];
  }
  if(direction == -1)
    rsArrayTools::scale(ad, N, -1);

  delete[] hk;
  delete[] binomCoeffs;
}

// maybe implement tha "antidifference" operator:
// http://faculty.cs.tamu.edu/klappi/csce411-s15/csce411-setFDCalculus.pdf

// see also here for more info about finite difference calculus:
//https://web.archive.org/web/20090419132601/http://www.stanford.edu/~dgleich/publications/finite-calculus.pdf

//-------------------------------------------------------------------------------------------------
// roots:

template<class T>
template<class R>
void rsPolynomial<T>::roots(const std::complex<R>* a, int degree, std::complex<R>* roots)
{
  const R eps = R(2.0e-14); // for float, it was 2.0e-6 - use template numeric_limit<T>
                            // maybe something like 100*epsilon?

  // allocate memory for the coefficients of the deflated polynomial and initialize it as
  // non-deflated polynomial:
  std::complex<R>* ad = new std::complex<R>[degree+1];
  rsArrayTools::copy(a, ad, degree+1);

  // loop over the roots:
  for(int j = degree; j >= 1; j--)
  {
    // find a root of the deflated polynomial using the Laguerre-algorithm with 0 as initial guess:
    std::complex<R> r = convergeToRootViaLaguerre(ad, j, std::complex<R>(0.0, 0.0));

    // polish the root by using the Laguerre method with the undeflated polynomial and the
    // non-polished root as initial guess:
    r = convergeToRootViaLaguerre(a, degree, r);

    // maybe move into a member function Complex::zeroNegligibleImaginaryPart(T ratio);
    // -> ratio = 2*eps (maybe leave this to client code):
    if(fabs(r.imag()) <= 2.0*eps*fabs(r.real()))
      r.imag(0);

    // store root in output array:
    roots[j-1] = r;

    // deflate the deflated polynomial again by the monomial that corresponds to our most recently
    // found root:
    std::complex<R> rem = ad[j];  // remainder - not used, needed for passing a dummy pointer
    rsPolynomial<std::complex<R>>::divideByMonomialInPlace(ad, j, r, &rem);
  }

  rsSortComplexArrayByReIm(roots, degree);  // maybe leave this to client code
  delete[] ad;
}

template<class T>
template<class R>
void rsPolynomial<T>::roots(const R* a, int degree, std::complex<R>* r)
{
  std::complex<R>* ac = new std::complex<R>[degree+1];
  rsArrayTools::convert(a, ac, degree+1);               // complexify real coeffs
  roots(ac, degree, r);
  delete[] ac;
}

template<class T>
template<class R>
std::complex<R> rsPolynomial<T>::convergeToRootViaLaguerre(
  const std::complex<R>* a, int degree, std::complex<R> initialGuess)
{
  const R eps = std::numeric_limits<R>::epsilon();

  static const int numFractions = 8; // number of fractions minus 1 (for breaking limit cycles)
  static const int itsBeforeFracStep = 10;  // number of iterations after which a fractional step
                                            // is taken (to break limit cycles)
  static const int maxNumIterations = itsBeforeFracStep*numFractions;

  // fractions for taking fractional update steps to break a limit cycles:
  static R fractions[numFractions+1] =
    { R(0.0),  R(0.5),  R(0.25), R(0.75), R(0.13), R(0.38), R(0.62), R(0.88), R(1.0) };

  std::complex<R> r = initialGuess; // the current estimate for the root
  for(int i = 1; i <= maxNumIterations; i++)
  {
    std::complex<R> P[3];    // holds P, P', P''
    R  err = eps * evaluateWithTwoDerivativesAndError(a, degree, r, P);

    if(abs(P[0]) <= err)
      return r;
    // "simplified stopping criterion due to Adams", referred to on page 373 (?)
    // can we get rid of this? if so, we might also replace the above loop by
    // evaluatePolynomialAndDerivativesAt

  // Laguerre's formulas:
    std::complex<R> G  = P[1]/P[0];       // Eq. 9.5.6
    std::complex<R> H  = G*G - P[2]/P[0]; // Eq. 9.5.7
    std::complex<R> sq = sqrt(R(degree-1)*(R(degree)*H-G*G)); // square-root Eq. 9.5.11
    std::complex<R> Gp = G + sq;  // denominator in 9.5.11 with positive sign for square-root
    std::complex<R> Gm = G - sq;  // denominator in 9.5.11 with negative sign for square-root

    // choose Gp or Gm according to which has larger magnitude (page 372, bottom), re-use Gp for
    // the result:
    R GpA = abs(Gp);
    R GmA = abs(Gm);
    if(GpA < GmA)
    {
      Gp  = Gm;
      GpA = GmA;
    }

    // compute difference between old and new estimate for the root r (the 'a' variable in
    // Eq. 9.5.8)
    std::complex<R> dr;
    if(GpA > 0.0)
      dr = std::complex<R>(R(degree), 0.0) / Gp;  // Eq. 9.5.11
    else
      dr = exp(log(R(1)+abs(r))) * std::complex<R>(cos((R)i), sin((R)i));
    // \todo use sinCos() or std::polar

    // compute new estimate for the root:
    std::complex<R> rNew = r - dr;
    if(r == rNew)
      return r;  // converged

    // update our r-variable to the new estimate:
    if(i % itsBeforeFracStep != 0)
      r = rNew;
    else
      r = r - fractions[i/itsBeforeFracStep]*dr; // fractional step to break limit cycle
  }

  rsError("Too many iterations taken, algorithm did not converge.");
  return 0.0;
}

template<class T>
T rsPolynomial<T>::rootLinear(const T& a, const T& b)
{
  if(a == T(0)) {  // hmm...maybe returning (+-)inf as root would actually be appropriate
    RS_DEBUG_BREAK;
    return T(0);
  }
  else
    return -b/a;
}

template<class T>
template<class R>
std::vector<std::complex<R>> rsPolynomial<T>::rootsQuadratic(const R& a, const R& b, const R& c)
{
  // catch degenerate case with zero leading coefficient:
  if(a == 0.0) {
    std::vector<std::complex<R>> roots(1);
    roots[0] = rootLinear(b, c);
    return roots;
  }

  std::vector<std::complex<R>> roots(2); // array to be returned
  R D      = b*b - R(4)*a*c;   // discriminant ...use discriminant-function
  R factor = R(1) / (R(2)*a);  // common factor that appears everywhere
  if(D > 0.0) {
    // D > 0: two distinct real roots:
    R rsSqrt_D = rsSqrt(D);
    roots[0]   = factor * (-b+rsSqrt_D);
    roots[1]   = factor * (-b-rsSqrt_D);
  }
  else if(D == 0.0) {
    // D == 0: a real root with multiplicity 2:
    roots[1] = roots[0] = std::complex<R>(-b * factor);
  }
  else {
    // D < 0: two complex conjugate roots:
    R imag   = rsSqrt(-D) * factor;
    R real   = -b       * factor;
    roots[0] = std::complex<R>(real, imag);
    roots[1] = std::complex<R>(real, -imag);
  }

  return roots;
}

template<class T>
template<class R>
void rsPolynomial<T>::rootsQuadraticReal(const R& c, const R& b, const R& a, R* x1, R* x2)
{
  // Equation:  a*x^2 + b*x + c = 0
  // Solutions: x1,x2 = (-b +- sqrt(b^2-4*a*c)) / (2*a):
  R s = T(1) / (2*a);         // scaler
  R d = b*b - 4*a*c;          // discriminant
  d   = sqrt(rsMax(d, R(0))); // we return the real part of the complex conjugate pair in case...
  *x1 = (-b-d) * s;           // ...of a negative discriminant, we return the roots in ascending...
  *x2 = (-b+d) * s;           // order, so the one with minus in the formula comes first
}
// is the formula (numerically) the same as the pq-formula? if not, which one is better -> test
// what about the degenerate case a=0?

template<class T>
template<class R>
void rsPolynomial<T>::rootsQuadraticComplex(
  const std::complex<R>& c, const std::complex<R>& b, const std::complex<R>& a,
  std::complex<R>* x1, std::complex<R>* x2)
{
  std::complex<R> s = R(1) / (R(2)*a);
  std::complex<R> d = sqrt(b*b - R(4)*a*c); // sqrt of discriminant
  *x1 = (-b-d) * s;
  *x2 = (-b+d) * s;
}

template<class T>
template<class R>
std::vector<std::complex<R>> rsPolynomial<T>::rootsCubic(
  const R& a, const R& b, const R& c, const R& d)
{
  // catch degenerate cases where the leading coefficient is zero:
  if(a == 0.0)
    return rootsQuadratic(b, c, d);

  std::vector<std::complex<R>> y(3);
  std::vector<std::complex<R>> roots(3);

  // compute p,q as in the Bronstein page 40, Eq. 1.154c and the offset for the substitution
  // y = x + b/(3*a):
  R p = (R(3)*a*c-b*b)/(R(9)*a*a);
  R q = (b*b*b)/(R(27)*a*a*a) - (b*c)/(R(6)*a*a) + d/(R(2)*a);

  R u, r, D, phi, ch, sh, re, im, tmp;

  if(p == 0.0 && q == 0.0)
  {
    y[0] = y[1] = y[2] = 0.0;         // a triple real root at y=0
    // checked with y = (x-1)^3 = x^3-3*x^2+3*x-1
  }
  else if(p != 0.0 && q == 0.0)
  {
    y[2] = 0.0;                       // a real root at y=0 and ...
    u    = -R(3)*p;
    if(u > 0.0)
    {
      tmp  =  rsSqrt(u);
      y[0] =  tmp;
      y[1] = -tmp;                    // ... two additional real roots or ...
      // checked with y = (x-1)*(x-2)*(x-3) = x^3-6*x^2+11*x-6
    }
    else // u < 0.0
    {
      tmp  =  rsSqrt(-u);
      y[0] =  std::complex<R>(0.0, tmp);
      y[1] =  std::complex<R>(0.0, -tmp);     // ... two imaginary roots
      // checked with y = (x-4i)*(x+4i)*(x-0) = x^3+16*x
    }
  }
  else if(p == 0.0 && q != 0.0)
  {
    u = -R(2)*q;
    if(u > 0.0)
    {
      tmp  = pow(u, R(1.0/3.0));
      y[2] = tmp;                     // a real root at a positive y or ...
      phi  = R((2.0/3.0)*PI);
      // checked with x^3+3*x^2+3*x
    }
    else // u < 0.0
    {
      tmp  = pow(-u, R(1.0/3.0));
      y[2] = -tmp;                    // ... a real root at a negative y and ...
      phi  = R(PI/3.0);
      // checked with x^3+3*x^2+3*x+10
    }
    rsSinCos(phi, &im, &re);
    re  *= tmp;
    im  *= tmp;
    y[0] = std::complex<R>(re, im);
    y[1] = std::complex<R>(re, -im);           // ... two complex conjugate roots
  }
  else // both p and q are nonzero
  {
    r = rsSign(q) * rsSqrt(fabs(p));
    if(p > 0.0)
    {
      phi = (R)rsAsinh(q/(r*r*r));
      rsSinhCosh(phi/R(3), &sh, &ch);
      y[0] = std::complex<R>(r*sh,  rsSqrt(R(3))*r*ch);
      y[1] = std::complex<R>(r*sh, -rsSqrt(R(3))*r*ch);
      y[2] = -R(2)*r*sh;
      // checked with y = (x-i)*(x+i)*(x-1) = x^3-x^2+x-1
    }
    else // p < 0.0
    {
      D = q*q + p*p*p;
      if(D > 0.0)
      {
        phi  = (R)rsAcosh(q/(r*r*r));
        rsSinhCosh(phi/R(3), &sh, &ch);
        y[0] = std::complex<R>(r*ch,  rsSqrt(R(3))*r*sh);
        y[1] = std::complex<R>(r*ch, -rsSqrt(R(3))*r*sh);
        y[2] = -R(2)*r*ch;
        // checked with y = (x-i)*(x+i)*(x-3) = x^3-3*x^2+x-3
      }
      else // D <= 0.0
      {
        phi  = acos(q/(r*r*r));
        y[0] =  R(2)*r*cos(R(PI/3.0 + phi/3.0));
        y[1] =  R(2)*r*cos(R(PI/3.0 - phi/3.0));
        y[2] = -R(2)*r*cos(R(phi/3.0));       // three distinct real roots
        // checked
      }
    }
  }

  // obtain the results for the original equation (back-substitution):
  R s = b/(R(3)*a);
  roots[0] = y[0]-s;
  roots[1] = y[1]-s;
  roots[2] = y[2]-s;

  return roots;
}
// this is a mess - use formulas from
// http://mathworld.wolfram.com/CubicFormula.html
// eq. 54..56


template<class T>
template<class R>
R rsPolynomial<T>::cubicDiscriminant(const R& d, const R& c, const R& b, const R& a)
{
  return b*b*c*c - R(4)*(a*c*c*c + b*b*b*d) - R(27)*a*a*d*d + R(18)*a*b*c*d;
}

template<class T>
T rsCubeRoot(T x)
{
  //return cbrt(x);
  return pow(x, T(1)/T(3));
}

template<class T>
template<class R>
void rsPolynomial<T>::rootsCubicComplex(
  std::complex<R> a0, std::complex<R> a1,
  std::complex<R> a2, std::complex<R> a3,
  std::complex<R>* r1, std::complex<R>* r2, std::complex<R>* r3)
{
  //rsAssert(false); // does not yet work - produces wrong results when roots are not real
  // formulas from http://mathworld.wolfram.com/CubicFormula.html

  // intermediate variables:
  std::complex<R> q, r, d, s, t, u, v, w; // maybe can use less intermediate variables by re-using
  q = R(1) / a3; a0 *= q; a1 *= q; a2 *= q;  // make monic (such that a3 == 1)
  q = (R(3)*a1 - a2*a2) * R(1.0/9.0);
  r = (R(9)*a2*a1 - R(27)*a0 - R(2)*a2*a2*a2) * R(1.0/54.0);
  d = q*q*q + r*r;     // discriminant
  d = sqrt(d);         // we actually need the square root of it
  s = rsCubeRoot(r+d);
  t = rsCubeRoot(r-d);
  u = -a2*R(1.0/3.0);
  v = (s+t) * R(0.5);
  w = (s-t) * (R(0.5*sqrt(3.0)) * std::complex<R>(0, 1)); // factor is constant - optimize

  // roots:
  *r1 = u - v + w;
  *r2 = u - v - w;
  *r3 = u + v + v;
}

template<class T>
template<class R>
R rsPolynomial<T>::cubicRootNear(R x, const R& a, const R& b, const R& c, const R& d,
  const R& min, const R& max, int maxIterations)
{
  R f    = ((a*x+b)*x+c)*x+d;
  R df   = (R(3)*a*x+R(2)*b)*x+c;
  R xNew = x - f/df;
  int i = 1;
  while(xNew != x && i < maxIterations) {
    x    = xNew;
    f    = ((a*x+b)*x+c)*x+d;
    df   = (R(3)*a*x+R(2)*b)*x+c;
    xNew = x - f/df;
    i++;
  }
  return rsClip(xNew, min, max);
}
// todo: re-write the algorithm such that the formulas only appear once (reduce code size)

template<class T>
template<class R>
R rsPolynomial<T>::rootNear(R x, const R* a, int degree, const R& min, const R& max,
  int maxIterations)
{
  // Newton/Raphson iteration:
  R f, df, xNew;
  evaluateWithDerivative(x, a, degree, &f, &df);
  xNew  = x - f/df;
  int i = 1;
  while(xNew != x && i < maxIterations) {  // maybe needs tolerance
    x    = xNew;
    evaluateWithDerivative(x, a, degree, &f, &df);
    xNew = x - f/df;
    i++;
  }
  return rsClip(xNew, min, max);
}
// drage the 1st iteration into the loop

//-------------------------------------------------------------------------------------------------
// conversions:

template<class T>
std::vector<std::complex<T>> rsPolynomial<T>::rootsToCoeffs(
  const std::vector<std::complex<T>>& roots)
{
  std::vector<std::complex<T>> coeffs(roots.size()+1);
  if(roots.size() < 1)
    return coeffs;
  rootsToCoeffs(&roots[0], &coeffs[0], (int) roots.size());
  return coeffs;
  /*
  // old:
  std::vector<std::complex<T>> coeffs;
  coeffs.reserve(roots.size()+1);
  coeffs.push_back(1.0); // init with one for convolutional accumulation - but in the special case
                         // roots.size() < 1, this would be wrong, right? ..hmm
  if(roots.size() < 1)
    return coeffs;
  for(int i = 0; i < roots.size(); i++) {
    std::complex<T> z = roots[i];
    coeffs.push_back(coeffs[i]);
    for(int j=i; j>=1; j--)
      coeffs[j] = coeffs[j-1] - z * coeffs[j];
    coeffs[0] = -z * coeffs[0];
  }
  return coeffs;
  */
}

template<class T>
void rsPolynomial<T>::rootsToCoeffs(const std::complex<T>* r, std::complex<T>* a, int N)
{
  std::complex<T>* rF = new std::complex<T>[N]; // only the finite roots
  int nF = rsCopyFiniteValues(r, rF, N);
  rsArrayTools::fillWithZeros(a, N+1);
  if(nF == 0)
    a[0] = 1.0;
  else {
    a[0] = -rF[0];
    a[1] = 1.0;
    for(int M = 2; M <= nF; M++) {
      a[M] = a[M-1];
      std::complex<T> rM = rF[M-1];
      for(int n = M-1; n >= 1; n--)
        a[n] = a[n-1] - rM*a[n];
      a[0] = -rM*a[0];
    }
  }
  delete[] rF;
  // todo: avoid memory allocation - check against infinity on the fly and skip the root, if it is
  // infinite, also - the function should not care whether the template type is real or complex -
  // the algo is the same in both cases
}
// why is this so complicated anyway? can'T we just use the version
//   rootsToCoeffs(const T* r, T* a, int N, T scaler)
// ...maybe make a second version that checks the root against infinity

template<class T>
template<class R>
void rsPolynomial<T>::complexRootsToRealCoeffs(const std::complex<R>* r, R* a, int N)
{
  std::complex<R>* ac = new std::complex<R>[N+1];
  rootsToCoeffs(r, ac, N);
  for(int n = 0; n <= N; n++)
  {
    constexpr R tol = R(1000) * RS_EPS(R);
    rsAssert(rsIsCloseTo(ac[n].imag(), R(0), tol),
      "roots do not occur in complex conjugate pairs");
    a[n] = ac[n].real();
  }
  delete[] ac;
}
// rename to complexRootsToRealCoeffs, assert that the imaginary parts are numerically close to
// zero

template<class T>
void rsPolynomial<T>::rootsToCoeffs(const T* r, T* a, int N, T scaler)
{
  rsArrayTools::fillWithZeros(a, N+1);
  a[0] = scaler;
  for(int n = 1; n <= N; n++)
    rsArrayTools::convolveWithTwoElems(a, n, -r[n-1], T(1), a);
}

template<class T>
void rsPolynomial<T>::newtonToMonomialCoeffs(T* x, T* a, int N)
{
  T x0 = x[0]; 
  x[0] = T(1);     // x is re-used as our convolutive accumulator (and destroyed in the process)
  for(int i = 1; i < N; i++) {
    T x1 = x[i];                                            // save x[i] because the next line..
    rsArrayTools::convolveWithTwoElems(x, i, -x0, T(1), x); // ..overwrites x up to x[i] but we..
    x0 = x1;                                                // ..still need it in next iteration
    for(int j = 0; j < i; j++)
      a[j] += a[i] * x[j]; }                                // accumulation of final coeffs
}


//-------------------------------------------------------------------------------------------------
// fitting:

template<class T>
void rsPolynomial<T>::cubicCoeffsTwoPointsAndDerivatives(T *a, const T *x, const T *y, const T *dy)
{
  // compute intermediate variables:
  T x0_2 = x[0]*x[0]; // x[0]^2
  T x0_3 = x0_2*x[0]; // x[0]^3
  T x1_2 = x[1]*x[1]; // x[1]^2
  T x1_3 = x1_2*x[1]; // x[1]^3
  T k1   = T(3)*x[0]*x1_2;
  T k2   = -T(3)*x[1]*y[1];
  T k3   = dy[1]-dy[0];
  T s    = T(1)/(-x1_3+k1-T(3)*x0_2*x[1]+x0_3);  // scaler

  a[0] =  s*(x0_2*(x1_2*k3+k2) + x0_3*(y[1]-x[1]*dy[1]) + x[0]*x1_3*dy[0] + y[0]*(-x1_3+k1));
  a[1] = -s*(x[0]*(x1_2*(T(2)*dy[1]+dy[0])-T(6)*x[1]*y[1]) - x0_3*dy[1]
             + x0_2*x[1]*(-dy[1]-T(2)*dy[0]) + x1_3*dy[0] + T(6)*x[0]*x[1]*y[0]);
  a[2] =  s*(x[0]*(x[1]*k3-T(3)*y[1]) + x1_2*(dy[1]+T(2)*dy[0]) + x0_2*(-dy[0]-T(2)*dy[1]) + k2
             + y[0]*(T(3)*x[1]+T(3)*x[0]));
  a[3] = -s*(x[1]*(dy[1]+dy[0]) + x[0]*(-dy[1]-dy[0]) - T(2)*y[1] + T(2)*y[0]);
}

template<class T>
void rsPolynomial<T>::cubicCoeffsTwoPointsAndDerivatives(T *a, const T *y, const T *dy)
{
  a[0] = y[0];
  a[1] = dy[0];
  a[2] = T(3)*(y[1]-a[1]-a[0])-dy[1]+a[1];
  a[3] = T(1.0/3.0)*(dy[1]-T(2)*a[2]-a[1]);
}

template<class T>
void rsPolynomial<T>::cubicCoeffsFourPoints(T *a, const T *y)
{
  a[0] = y[0];
  a[2] = T(0.5)*(y[-1]+y[1]-T(2)*a[0]);
  a[3] = T(1.0/6.0)*(y[2]-y[1]+y[-1]-a[0]-T(4)*a[2]);
  a[1] = T(0.5)*(y[1]-y[-1]-T(2)*a[3]);
}

// move to a function rsMatrixView:vandermonde(const T* x, int N, T* V)
template<class T>
T** rsPolynomial<T>::vandermondeMatrix(const T *x, int N)
{
  T **A; rsArrayTools::allocateSquareArray2D(A, N);
  for(int i = 0; i < N; i++) {
    T xi  = x[i];
    T xij = 1.0;  // xi^j
    for(int j = 0; j < N; j++) {
      A[i][j] = xij;
      xij *= xi; }}
  return A;
}

/*
// Old implementation using Gaussian elimination with the Vandermonde matrix:
template<class T>
void rsPolynomial<T>::interpolant(T *a, const T *x, const T *y, int N)
{
  T **A = vandermondeMatrix(x, N);
  rsLinearAlgebra::rsSolveLinearSystem(A, a, y, N); // use rsSolveLinearSystemInPlace
  rsArrayTools::deAllocateSquareArray2D(A, N);

  // For higher degree polynomials, this simple and direct approach may become numerically ill
  // conditioned. In this case, we could first normalize the data, such than xMin = yMin = -1 and
  // xMax = yMax = +1, find coefficients for Chebychev polynomials, such that our normalized
  // interpolant is given by P(x) = b0*T0(x) + b1*T1(x) + ... + (bN-1)*(TN-1)(x) where Tn is the
  // n-th degree Chebychev polynomial. Then, the b-coefficients could be converted back to
  // a-coefficients for powers of x and finally, we could denormalize using rsShiftPolynomial,
  // rsStretchPolynomial (for x-denormalization) and scaling the coeffs and adding an offset to
  // a[0] for y-denormalization
}
*/

// New implementation, using Lagrange's idea - uses only O(N) memory ..i think, the runtime is
// still O(N^3) - we have a loop nesting level of 2 here and one of the inner loops calls an O(N)
// convolution (with only two elements, that's why it's only O(N)). I think, we can reduce this by
// getting rid of the convolution inside the loop - just build up the product of all factors up
// front once (without leaving out the n-th linear factor) - then, inside the loop, just divide
// out the n-th factor again - this one call to dividePolynomialByMonomial which has O(N) instead
// of O(N^2) for the multiplicative accumulation loop - so the overall runtime would be O(N^2)...
// but roundoff error might be higher due to first multiplying in and then later dividing out
// linear factors

template<class T>
void rsPolynomial<T>::interpolant(T* a, const T* x, const T* y, int N)
{
  T* wrk = new T[N+1];  // do we really need N+1? isn't N enough?
  interpolant(a, x, y, N, &wrk[0]);
  delete[] wrk;
}


template<class T>
void rsPolynomial<T>::interpolant(T* a, const T* x, const T* y, int N, T* wrk)
{
  using AT = rsArrayTools;
  AT::fillWithZeros(a, N);
  T* num = &wrk[0];
  for(int n = 0; n < N; n++)
  {
    // init num and den to 1:
    AT::fillWithZeros(num, N);
    num[0] = 1;
    T den = 1;

    // convolutive and multiplicative accumulation of num and den:
    for(int k = 0; k < N; k++) {
      if(k != n) {
        AT::convolveWithTwoElems(num, k+1, -x[k], T(1), num);
        den *= x[n] - x[k];  }}

    // accumulate this result additively into coeff-array:
    T s =  y[n]/den;
    for(int k = 0; k < N; k++)
      a[k] += num[k] * s;
  }
}

// ToDo: implement Newton's formula using divided differences. I think, it needs a workspace of 
// size 2N (twice as much as Lagrange, but it's still O(N)) and is O(N^2) in time (compared to 
// O(N^3) for the Lagrange algorithm

template<class T>
void rsPolynomial<T>::interpolant(T *a, const T& x0, const T& dx, const T *y, int N)
{
  T *x = new T[N];
  for(int n = 0; n < N; n++)
    x[n] = x0 + T(n)*dx;
  interpolant(a, x, y, N);
  delete[] x;
}

template<class T>
void rsPolynomial<T>::interpolantViaNewton(T* a, const T* x, const T* y, int N, T* w)
{
  rsArrayTools::copy(x, w, N);
  rsArrayTools::copy(y, a, N);
  interpolantViaNewtonInPlace(w, a, N);
}

template<class T>
void rsPolynomial<T>::interpolantViaNewtonInPlace(T* x, T* y, int N)
{
  coeffsNewton(          x, y, N);   // overwrites y
  newtonToMonomialCoeffs(x, y, N);   // overwrites y again and also x
}

template<class T>
void rsPolynomial<T>::fitQuadraticDirect(T *a, const T *x, const T *y)
{
  T k1 = y[1] - y[0];
  T k2 = x[0]*x[0] - x[1]*x[1];
  T k3 = x[1] - x[0];
  T k4 = k1/k3;
  T k5 = k2/k3;

  a[2] = (y[2]-y[0]+k4*(x[0]-x[2])) / (x[2]*(x[2]+k5)-x[0]*(x[0]+k5));
  a[1] = (k1+k2*a[2])/k3;
  a[0] = y[0]-a[2]*x[0]*x[0]-a[1]*x[0];

  // operations: add: 4, sub: 8, mul: 9, div: 4, ass: 8, tmp: 5
}
// todo: benchmark against fitQuadraticLagrange and also compare numeric accuracy of both

template<class T>
void rsPolynomial<T>::fitQuadraticLagrange(T* a, const T* x, const T* y)
{
  T k1 = y[0] / ((x[0]-x[1])*(x[0]-x[2]));
  T k2 = y[1] / ((x[1]-x[0])*(x[1]-x[2]));
  T k3 = y[2] / ((x[2]-x[0])*(x[2]-x[1]));
  T b1 = -k1*(x[1]+x[2]);
  T b2 = -k2*(x[0]+x[2]);
  T b3 = -k3*(x[0]+x[1]);
  T c1 = k1*x[1]*x[2];
  T c2 = k2*x[0]*x[2];
  T c3 = k3*x[0]*x[1];
  a[2] = k1 + k2 + k3;
  a[1] = b1 + b2 + b3;
  a[0] = c1 + c2 + c3;

  // operations: add: 9, sub: 6, mul: 12, div: 3, neg: 3, ass: 12, tmp: 9
}
// todo: optimize out the b and c variables so we need only 3 temporaries and 6 assignments, but
// i guess, the compiler can do this, too - check at godbolt.org
// maybe derive and implement simplified formulas for the common special case
// x1 = -1, x2 = 0, x3 = +1...hmm - but i think, we will not get formulas any better than what we
// already have below:

template<class T>
void rsPolynomial<T>::fitQuadratic_0_1_2(T *a, const T *y)
{
  a[2] = T(0.5)*(y[0]+y[2])-y[1];
  a[1] = y[1]-y[0]-a[2];
  a[0] = y[0];
}

template<class T>
void rsPolynomial<T>::fitQuadratic_m1_0_1(T *a, const T *y)
{
  a[0] = y[1];
  a[1] = T(0.5)*(y[2]-y[0]);
  a[2] = y[2] - a[0] - a[1];
}

template<class T>
T rsPolynomial<T>::quadraticExtremumPosition(const T *a)
{
  return T(-0.5) * a[1]/a[2]; // it's the client's responsibility to ensure that a[2] is nonzero
}

template<class T>
void rsPolynomial<T>::fitQuarticWithDerivatives(T *a, const T *y, const T& s0, const T& s2)
{
  a[0] = y[0];
  a[1] = s0;
  a[2] = -(T(5)*y[2]-T(16)*y[1]+T(11)*y[0]-T(2)*s2+T(8)*s0)/T(4);
  a[3] =  (T(7)*y[2]-T(16)*y[1]+T( 9)*y[0]-T(3)*s2+T(5)*s0)/T(4);
  a[4] = -(T(2)*y[2]-T( 4)*y[1]+T( 2)*y[0]-  s2  +s0)/T(4);
}

template<class T>
template<class R>
bool rsPolynomial<T>::areRootsOnOrInsideUnitCircle(const R& a0, const R& a1, const R& a2)
{
  // p and q values for p-q formula
  R p = a1/a2;
  R q = a0/a2;
  R d = p*p/R(4) - q; // value under the square-root
  R rr, ri;
  if( d < 0 )
  {
    // complex conjugate roots
    rr = -p/R(2);     // real part
    ri = rsSqrt(-d); // imaginary part
    return rsSqrt(rr*rr + ri*ri) <= 1.0;
  }
  else
  {
    // 2 real roots
    d = rsSqrt(d);
    rr = -p/2 + d;
    if( rsAbs(rr) > 1.0 )
      return false;
    rr = -p/2 - d;
    if( rsAbs(rr) > 1.0 )
      return false;
    return true;
  }
}

//-------------------------------------------------------------------------------------------------
// coefficient generation for special polynomials:

template <class T>
void rsPolynomial<T>::threeTermRecursion(T* a, const T& w0, int degree, const T* a1, const T& w1,
  const T& w1x, const T* a2, const T& w2)
{
  rsAssert(degree >= 2);
  int n = degree;
  a[n] = (w1x*a1[n-1]) / w0;
  n--;
  a[n] = (w1*a1[n] + w1x*a1[n-1]) / w0;
  for(n = n-1; n > 0; n--)
    a[n] = (w1*a1[n] + w1x*a1[n-1] + w2*a2[n]) / w0;
  a[0] = (w1*a1[0] + w2*a2[0]) / w0;
  // optimize: replace divisions by w0 by multiplications
}

template<class T>
void rsPolynomial<T>::besselPolynomial(T *a, int degree)
{
  int m, n;
  for(n=0; n<=degree; n++)
    a[n] = 0.0;

  if( degree == 0 ) {
    a[0] = 1.0;
    return; }
  else if( degree == 1 )  {
    a[0] = 1.0;
    a[1] = 1.0;
    return; }

  // the general case is treated by recursion:
  a[0]  = 1.0;
  T *b1 = new T[degree+1]; // why +1?
  T *b2 = new T[degree+1];
  b2[0] = 1.0; b2[1] = 0.0;
  b1[0] = 1.0; b1[1] = 1.0;
  for(n = 2; n <= degree; n++) {
    T c = (T) (2*n-1);
    for(m = n; m >  0; m--)   a[m]  = c*b1[m-1];
    a[0] = 0;
    for(m = 0; m <  n-1; m++) a[m] += b2[m];
    for(m = 0; m <  n;   m++) b2[m] = b1[m];
    for(m = 0; m <= n;   m++) b1[m] = a[m];  }
  delete[] b1;
  delete[] b2;
}
// can this be optimized? i think so

template<class T>
template<class R>
void rsPolynomial<T>::legendrePolynomial(R *a, int degree)
{
  if(degree == 0) { a[0] = 1.0;             return; }
  if(degree == 1) { a[0] = 0.0; a[1] = 1.0; return; }

  a[0] = -0.5;
  a[1] =  0.0;
  a[2] =  1.5;
  if(degree == 2)
    return;

  int i, j;
  R *b1 = new R [degree+1];
  R *b2 = new R [degree+1];

  for(i = 0; i <= degree; i++) {
    b1[i] = b2[i] = 0.0;
  }
  b2[1] = 1.0;

  for(i = 3; i <= degree; i++) {
    for(j = 0; j <= i; j++) {
      b1[j] = b2[j];
      b2[j] = a[j];
      a[j]  = 0.0;
    }
    for(j = i-2; j >= 0; j-=2) {
      a[j] -= (i-1)*b1[j]/i;
    }
    for(j = i-1; j >= 0; j-=2) {
      a[j+1] += (2*i-1)*b2[j]/i;
    }
  }
  delete [] b1;
  delete [] b2;
}

template<class T>
void rsPolynomial<T>::jacobiRecursionCoeffs(int n, T a, T b, T *w0, T *w1, T *w1x, T *w2)
{
  T k  = T(2*n)+a+b;
  *w0  = T(2*n)*(T(n)+a+b)*(k-T(2));
  *w1  = (k-T(1))*(a*a-b*b);
  *w1x = k*(k-T(1))*(k-T(2));
  *w2  = -T(2)*k*(T(n)+a-T(1))*(T(n)+b-T(1));
  // from: https://en.wikipedia.org/wiki/Jacobi_polynomials#Recurrence_relation
}

template<class T>
void rsPolynomial<T>::jacobiRecursion(T *c, int n, T *c1, T *c2, T a, T b)
{
  // initialization:
  if( n == 0 ) {
    c[0] = 1.0;
    return; }
  if( n == 1 ) {
    c[0] = T(0.5)*(a-b);
    c[1] = T(0.5)*(a+b+T(2));
    return; }

  // recursion:
  T w0, w1, w1x, w2;
  jacobiRecursionCoeffs(n, a, b, &w0, &w1, &w1x, &w2);
  threeTermRecursion(c, w0, n, c1, w1, w1x, c2, w2);
}

template<class T>
void rsPolynomial<T>::jacobiPolynomials(T **c, T a, T b, int maxDegree)
{
  for(int n = 0; n <= maxDegree; n++)
    jacobiRecursion(c[n], n, c[n-1], c[n-2], a, b);
}

template<class T>
void rsPolynomial<T>::legendreRecursion(T *a, int n, T *a1, T *a2)
{
  if( n == 0 ) { a[0] = 1.0;             return; }
  if( n == 1 ) { a[0] = 0.0; a[1] = 1.0; return; }
  threeTermRecursion(a, T(n), n, a1, 0.0, T(2*n)-T(1), a2, -(T(n)-T(1))); // soem Ts can be moved out
  // Legendre polynomials are a special case of Jacobi polynomials, so this would also work:
  // jacobiRecursion(a, n, a1, a2, 0.0, 0.0);
}

template<class T>
template<class R>
void rsPolynomial<T>::maxSlopeMonotonic(R *w, int n)
{
  R *a,*p,*s,*v,c0,c1;
  int i,j,k;

  a = new R[n+1];
  p = new R[2*n+1];
  s = new R[2*n+1];
  v = new R[2*n+4];

  k = (n-1)/2;
  //  n = 2k + 1 for odd 'n' and n = 2k + 2 for even 'n':
  //  n: 1 2 3 4 5 6 ...
  //  k: 0 0 1 1 2 2 ...

  //  form vector of 'a' constants:
  if(n & 1)                   // odd
  {
    for(i = 0; i <= k; i++)
    {
      //a[i] = (2.0*i+1.0)/(M_SQRT2*(k+1.0));
      a[i] = R(2*i+1) / R(SQRT2*(k+1));
    }
  }                           // even
  else
  {
    for(i = 0; i < k+1; i++)
    {
      a[i] = 0.0;
    }
    if(k & 1)
    {
      for(i = 1; i <= k; i+=2)
      {
        a[i] = (2*i+1)/sqrt((R) ((k+1)*(k+2)));
      }
    }
    else
    {
      for(i = 0; i <= k; i+=2)
      {
        a[i] = (2*i+1)/sqrt((R) ((k+1)*(k+2)));
      }
    }
  }
  for(i = 0; i <= n; i++)
  {
    s[i] = 0.0;
    w[i] = 0.0;
  }

  // form s[] = sum of a[i]*P[i]
  s[0] = a[0];
  s[1] = a[1];
  for(i = 2; i <= k; i++)
  {
    legendrePolynomial(p, i);
    for(j = 0; j <= i; j++)
    {
      s[j] += a[i]*p[j];
    }
  }

  //  form v[] = square of s[]:
  for(i = 0; i <= 2*k+2; i++)
  {
    v[i] = 0.0;
  }
  for(i = 0; i <= k; i++)
  {
    for(j = 0; j <= k;j++)
    {
      v[i+j] += s[i]*s[j];
    }
  }

  // modify integrand for even 'n':
  v[2*k+1] = 0.0;
  if((n & 1) == 0)
  {
    for(i = n; i >= 0; i--)
    {
      v[i+1] += v[i];
    }
  }

  // form integral of v[]:
  for(i = n+1; i >= 0; i--)
  {
    v[i+1] = v[i]/(T)(i+1.0);
  }
  v[0] = 0.0;

  // clear s[] for use in computing definite integral:
  for(i = 0; i < (n+2); i++)
  {
    s[i] = 0.0;
  }
  s[0] = -1.0;
  s[1] =  2.0;

  // calculate definite integral:
  for(i = 1; i <= n; i++)
  {
    if(i > 1)
    {
      c0 = -s[0];
      for(j = 1; j < i+1; j++)
      {
        c1 = -s[j] + 2*s[j-1];
        s[j-1] = c0;
        c0 = c1;
      }
      c1 = 2*s[i];
      s[i] = c0;
      s[i+1] = c1;
    }
    for(j = i; j > 0; j--)
    {
      w[j] += (v[i]*s[j]);
    }
  }
  if((n & 1) == 0)
    w[1] = 0.0;

  delete [] v;
  delete [] p;
  delete [] s;
  delete [] a;
}
// optimize this - less allocations

template<class T>
void rsPolynomial<T>::coeffsNewton(const T* x, T* y, int N)
{
  for(int i = 0; i < N; i++)
    for(int j = i+1; j < N; j++)
      y[j] = (y[j] - y[i]) / (x[j] - x[i]);
}
// https://en.wikipedia.org/wiki/Polynomial_interpolation#Non-Vandermonde_solutions
// https://en.wikipedia.org/wiki/Neville%27s_algorithm

/*
ToDo:
-check the rootsToCoeffs functions - these can (should!) be improved

-for those static functions that explicitly expect real or complex parameters, use a different
 template parameter - not T - but rather R for real and complex<R> for complex values
 -this prepares the class to be instantiated for real and complex coefficient types
 -see rootsQuadraticReal for how this works
 -...maybe it would be more elegant, if these functions are factored out into a separate class
  like rsPolynomialAlgorithms which has only static functions - so we don't have carry along the
  superfluous T template parameter
 -todo: evaluateFromRoots


*/


/*
// moved to rsRationalFunction:
template<class T>
void rsPolynomial<T>::rsPartialFractionExpansion(
  std::complex<T> *numerator,   int numeratorDegree,
  std::complex<T> *denominator, int denominatorDegree,
  std::complex<T> *poles, int *multiplicities, int numDistinctPoles,
  std::complex<T> *pfeCoeffs)
{
  // sanity check for inputs:
  rsAssert(numeratorDegree < denominatorDegree);
  rsAssert(rsArrayTools::sum(multiplicities, numDistinctPoles) == denominatorDegree);

  // make denominator monic:
  rsArrayTools::scale(numerator,   numeratorDegree+1,   T(1)/denominator[denominatorDegree]);
  rsArrayTools::scale(denominator, denominatorDegree+1, T(1)/denominator[denominatorDegree]);

  // establish coefficient matrix:
  std::complex<T> **A; rsArrayTools::allocateSquareArray2D(A, denominatorDegree);
  std::complex<T> *tmp = new std::complex<T>[denominatorDegree+1]; // deflated denominator
  std::complex<T> remainder;                                   // always zero
  for(int i = 0, k = 0; i < numDistinctPoles; i++)
  {
    rsArrayTools::copy(denominator, tmp, denominatorDegree+1);
    for(int m = 0; m < multiplicities[i]; m++)
    {
      divideByMonomialInPlace(tmp, denominatorDegree-m, poles[i], &remainder);
      for(int j = 0; j < denominatorDegree; j++)
        A[j][k] = tmp[j];
      k++;
    }
  }

  // solve the linear system using an appropriately zero-padded numerator as RHS:
  rsArrayTools::copy(numerator, tmp, numeratorDegree+1);
  rsArrayTools::fillWithZeros(&tmp[numeratorDegree+1], denominatorDegree-(numeratorDegree+1));
  rsLinearAlgebra::rsSolveLinearSystem(A, pfeCoeffs, tmp, denominatorDegree);

  // clean up:
  rsArrayTools::deAllocateSquareArray2D(A, denominatorDegree);
  delete[] tmp;
}

*/



/*
todo: define the set of rational functions R(x) = P(x)/Q(x) where P, Q are polynomials
-find out, if the set is closed under addition, subtraction, multiplication, division, composition,
 (maybe) decomposition, differentiation, integration (indefinite and definite with other rational
 functions as limits)
-implement all these operations
-note: when there's a direct part added to a proper rational function, the whole thing can be
 expressed as improper rational function - so we don't need extra provisions for the direct part
-maybe it's more convenient to work with pole/zero/gain representations instead of coefficients for
 some operations

composition: let y = R(x) = P(x)/Q(x), z = U(y) = S(y)/T(y) such that
z = S( P(x)/Q(x) ) / T( P(x)/Q(x) ) is the composition of U with R
define W(x) = ((P/Q)°S)(x) = S( P(x)/Q(x) ) and
       V(x) = ((P/Q)°T)(x) = T( P(x)/Q(x) )
we see, that we first need the composition of an (inner) rational function with an (outer)
polynomial
example: let P(x) = x^2 + 2x, Q(x) = 2x^3 - x^2, y = P(x)/Q(x), S(y) = 2y^2 - y + 3, so:
S(x) = 2((x^2+2x)/(2x^3-x^2))^2 - ((x^2+2x)/(2x^3-x^2))^1 + 3((x^2+2x)/(2x^3-x^2))^0
     = 2((x^2+2x)/(2x^3-x^2))^2 - (x^2+2x)/(2x^3-x^2) + 3
S = 2*(P/Q)^2 - 1*(P/Q)^1 + 3*(P/Q)^0
  = (2*P^2*Q^0 - 1*P^1*Q^1 + 3*P^0*Q^2) / Q^2
-> denominator and numerator polynomial can be computed
   degree of new denominator = degree(S)*degree(Q)
   degree of new numerator = ?
-> make a class for representing polynomials, so we can implement the algorithms conveniently, like
   S.num = 2*P^2*Q^0 - 1*P^1*Q^1 + 3*P^0*Q^2
   S.den = Q^2
   when this works, optimize to avoid excessive memory copying and shlemiel algorithms

differentiation: use the quotient rule

integration: make use of the partial fraction expansion (and later re-assemble the partial integrals
 into a canonical representation)
 -maybe provide a function for summing an array of rational functions

division: can be reduced to multiplication: R = P/Q, U = S/T
 -> W = R/U = (P/Q) / (S/T) = (P*T)/(Q*S)
 W.num = P*T; W.den = Q*S
 -maybe it can then be reduced to lowest terms - divide out common roots in numerator and denominator
  (we should have a function "reduce" for that that implements a polynomial gcd algorithm - or it
  can be based on a linear-factor decomposition)
 -or maybe it should first find the roots

maybe for testing, it will be convenient to use polynomials and rational functions with rational
coefficients -> make a rational number class (using unsigned integers for numerator and denominator
and a bool for the negaive sign)..the (unsigned) integer class could also be a template parameter,
so we may use rsBigInteger
...rational numbers could be constructed from floating point numbers by computing their continued
fraction expansion

*/

