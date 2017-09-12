
//-----------------------------------------------------------------------------------------------
// template function definitions - move to .cpp-file:

template <class T>
T rsPolynomial<T>::evaluatePolynomialAt(T x, T *a, int order)
{
  if(order < 0)
    return T(0);
  T y = a[order];
  for(int i = order-1; i >= 0; i--)
    y = y*x + a[i];
  return y;
}

template <class T>
void rsPolynomial<T>::evaluatePolynomialAndDerivativeAt(T x, T *a, int order, T *y, T *yd)
{
  *y  = a[order];
  *yd = 0.0;
  for(int i = order-1; i >= 0; i--)
  {
    *yd = *yd * x + *y;
    *y  = *y  * x + a[i];
  }
}

template <class T>
void rsPolynomial<T>::evaluatePolynomialAndDerivativesAt(T x, T *a, int order, T *results, 
  int numDerivatives)
{
  results[0] = a[order];
  rsFillWithZeros(&results[1], numDerivatives);
  for(int i = order-1; i >= 0; i--)
  {
    int n = rsMin(numDerivatives, order-1);
    for(int j = n; j >= 1; j--)
      results[j] = results[j]*x + results[j-1];
    results[0] = results[0]*x + a[i];
  }
  rsMultiply(&results[2], &rsFactorials[2], &results[2], numDerivatives-1);
}

template <class T>
void rsPolynomial<T>::multiplyPolynomials(T *a, int aOrder, T *b, int bOrder, T *result)
{
  rsConvolve(a, aOrder+1, b, bOrder+1, result);
}

template <class T>
void rsPolynomial<T>::dividePolynomials(T *p, int pOrder, T *d, int dOrder, T *q, T *r)
{
  rsCopyBuffer(p, r, pOrder+1); // init remainder with p
  rsFillWithZeros(q, pOrder+1); // init quotient with zeros
  for(int k = pOrder-dOrder; k >= 0; k--)
  {
    q[k] = r[dOrder+k] / d[dOrder];
    for(int j = dOrder+k-1; j >= k; j--)
      r[j] -= q[k] * d[j-k];
  }
  rsFillWithZeros(&r[dOrder], pOrder-dOrder+1);
}

template <class T>
void rsPolynomial<T>::dividePolynomialByMonomialInPlace(T *dividendAndResult, int dividendOrder, 
  T x0, T *remainder)
{
  *remainder = dividendAndResult[dividendOrder];
  dividendAndResult[dividendOrder] = T(0);
  for(int i=dividendOrder-1; i>=0; i--)
  {
    T swap               = dividendAndResult[i];
    dividendAndResult[i] = *remainder;
    *remainder           = swap + *remainder*x0;
  }
}

template <class T>
void rsPolynomial<T>::polyCoeffsForNegativeArgument(T *a, T *am, int N)
{
  double s = 1.0;
  for(int n = 0; n <= N; n++)
  {
    am[n]  = s*a[n];
    s     *= -1.0;
  }
}

// todo: polyCoeffsForScaledArgument: aScaled[n] = a[n] * scaler^n - when the scaler equals -1, 
// it reduces to polyCoeffsForNegativeArgument - this function is superfluous then

template <class T>
void rsPolynomial<T>::polyCoeffsForShiftedArgument(T *a, T *as, int N, T x0)
{
  int numLines = N+1;
  int length   = (numLines*(numLines+1))/2;
  rsUint32 *pt = new rsUint32[length];
  rsCreatePascalTriangle(pt, numLines);
  T *x0n = new T[N+1];  // +- x0^n
  x0n[0] = 1.0;
  for(int n = 1; n <= N; n++)
    x0n[n] = -x0 * x0n[n-1];
  for(int n = 0; n <= N; n++)
  {
    as[n] = 0.0;
    for(int k = n; k <= N; k++)
      as[n] += rsPascalTriangle(pt, k, k-n) * x0n[k-n] * a[k];
  }
  delete[] pt;
  delete[] x0n;
}

template <class T>
void rsPolynomial<T>::polyDerivative(T *a, T *ad, int N)
{
  for(int n = 1; n <= N; n++)
    ad[n-1] = n * a[n];
}

template <class T>
void rsPolynomial<T>::polyFiniteDifference(T *a, T *ad, int N, int direction, T h)
{
  // (possibly alternating) powers of the stepsize h:
  T *hk = new T[N+1];
  T hs  = direction*h;
  hk[0] = T(1);
  for(int k = 1; k <= N; k++)
    hk[k] = hk[k-1] * hs;

  // binomial coefficients:
  int numCoeffs    = N+1;
  int triangleSize = (numCoeffs*(numCoeffs+1))/2;
  rsUint32 *binomCoeffs = new rsUint32[triangleSize];
  rsCreatePascalTriangle(binomCoeffs, numCoeffs);

  // actual coefficient computation for ad:
  rsFillWithZeros(ad, N);
  for(int n = 0; n <= N; n++)
  {
    for(int k = 1; k <= n; k++)
      ad[n-k] += a[n] * rsPascalTriangle(binomCoeffs, n, k) * hk[k];
  }
  if(direction == -1)
    rsScale(ad, N, -1);

  delete[] hk;
  delete[] binomCoeffs;
}

template <class T>
void rsPolynomial<T>::polyIntegral(T *a, T *ai, int N, T c)
{
  for(int n = N+1; n >= 1; n--)
    ai[n] = a[n-1] / n;
  ai[0] = c;
}

template <class T>
void rsPolynomial<T>::createPolynomialPowers(T *a, int N, T **aPowers, int highestPower)
{
  aPowers[0][0] = 1;
  if(highestPower < 1)
    return;
  rsCopyBuffer(a, aPowers[1], N+1);
  for(int k = 2; k <= highestPower; k++)
    rsConvolve(aPowers[k-1], (k-1)*N+1, a, N+1, aPowers[k]);
}

template <class T>
void rsPolynomial<T>::composePolynomials(T *a, int aN, T *b, int bN, T *c)
{
  int cN = aN*bN;
  T *an  = new T[cN+1];  // array for the successive powers of a[]
  an[0]  = T(1);         // initialize to a[]^0

                         // accumulation:
  rsFillWithZeros(c, cN+1);
  c[0] = b[0];
  int K = 1;
  for(int n = 1; n <= bN; n++)
  {
    rsConvolveInPlace(an, K, a, aN+1);
    K += aN;
    for(int k = 0; k < K; k++)
      c[k] += b[n] * an[k];
  }

  delete[] an;
}

template <class T>
void rsPolynomial<T>::rsPolynomialRecursion(T *a, T w0, int order, T *a1, T w1, T w1x, T *a2, T w2)
{
  rsAssert(order >= 2);
  int n = order;
  a[n] = (w1x*a1[n-1]) / w0;
  n--;
  a[n] = (w1*a1[n] + w1x*a1[n-1]) / w0;
  for(n = n-1; n > 0; n--)
    a[n] = (w1*a1[n] + w1x*a1[n-1] + w2*a2[n]) / w0;
  a[0] = (w1*a1[0] + w2*a2[0]) / w0;
}

template<class T>
void rsPolynomial<T>::weightedSumOfPolynomials(T *p, int pN, T wp, T *q, int qN, T wq, T *r)
{
  int i;
  if(pN >= qN)
  {
    for(i = 0; i <= qN; i++)
      r[i] = wp*p[i] + wq*q[i];
    for(i = qN+1; i <= pN; i++)
      r[i] = wp*p[i];
  }
  else
  {
    for(i = 0; i <= pN; i++)
      r[i] = wp*p[i] + wq*q[i];
    for(i = pN+1; i <= qN; i++)
      r[i] = wq*q[i];
  }
}

template<class T>
void rsPolynomial<T>::subtractPolynomials(T *p, int pN, T *q, int qN, T *r)
{
  weightedSumOfPolynomials(p, pN, T(1), q, qN, T(-1), r);
}

template<class T>
void rsPolynomial<T>::integratePolynomialWithPolynomialLimits(T *p, int pN, T *a, int aN, T *b, 
  int bN, T *q)
{
  int PN = pN+1;
  int AN = aN*PN;
  int BN = bN*PN;

  T *P = new T[PN+1];
  T *A = new T[AN+1];
  T *B = new T[BN+1];

  polyIntegral(p, P, pN);               // P(x) is the antiderivative of p(x)
  composePolynomials(a, aN, P, PN, A);  // A(x) = P(a(x))
  composePolynomials(b, bN, P, PN, B);  // B(x) = P(b(x)) 
  subtractPolynomials(B, BN, A, AN, q); // q(x) = B(x) - A(x)

  delete[] P;
  delete[] A;
  delete[] B;
}

template<class T>
bool rsPolynomial<T>::rsPolynomialBaseChange(T **Q, T *a, T **R, T *b, int order)
{
  return rsChangeOfBasisRowWise(Q, R, a, b, order+1);
}

// end of being moved from .h file
//-----------------------------------------------------------------------------------------------

template<class T>
std::complex<T> rsPolynomial<T>::evaluatePolynomialWithRoots(std::complex<T> s, 
  std::complex<T> *r, int N)
{
  std::complex<T> result = 1.0;
  for(int i = 0; i < N; i++)
  {
    if( !r[i].isInfinite() )
      result *= (s - r[i]);
  }
  return result;
}

// used in convergeToRootViaLaguerre - maybe turn into member function:
template<class T>
double evaluatePolynomialWithTwoDerivativesAndError(std::complex<T> *a, int order, 
  std::complex<T> z, std::complex<T> *P)
{
  P[0] = a[order];                // P(z)
  P[1] = std::complex<T>(0.0, 0.0);  // P'(z)
  P[2] = std::complex<T>(0.0, 0.0);  // P''(z)
  double err = P[0].getRadius();  // estimated roundoff error in evaluation of the polynomial
  double zA  = z.getRadius();     // absolute value of z
  for(int j = order-1; j >= 0; j--) 
  {
    P[2] = z * P[2] + P[1];
    P[1] = z * P[1] + P[0];
    P[0] = z * P[0] + a[j];
    err  = P[0].getRadius() + zA*err;
  }
  P[2] *= 2.0;
  return err;
}

template<class T>
std::complex<T> rsPolynomial<T>::convergeToRootViaLaguerre(std::complex<T> *a, int order, 
                                               std::complex<T> initialGuess)
{
  const double eps = std::numeric_limits<double>::epsilon(); 

  static const int numFractions = 8; // number of fractions minus 1 (for breaking limit cycles)
  static const int itsBeforeFracStep = 10;  // number of iterations after which a fractional step 
                                            // is taken (to break limit cycles)
  static const int maxNumIterations = itsBeforeFracStep*numFractions;  

  // fractions for taking fractional update steps to break a limit cycles:
  static double fractions[numFractions+1] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0}; 
       
  std::complex<T> r = initialGuess; // the current estimate for the root
  for(int i = 1; i <= maxNumIterations; i++) 
  {
    std::complex<T> P[3];    // holds P, P', P'' 
    double  err = eps * evaluatePolynomialWithTwoDerivativesAndError(a, order, r, P);

    if( P[0].getRadius() <= err ) 
      return r;   
      // "simplified stopping criterion due to Adams", referred to on page 373 (?)
      // can we get rid of this? if so, we might also replace the above loop by 
      // evaluatePolynomialAndDerivativesAt

    // Laguerre's formulas:
    std::complex<T> G  = P[1]/P[0];       // Eq. 9.5.6
    std::complex<T> H  = G*G - P[2]/P[0]; // Eq. 9.5.7

    std::complex<T> sq = rsSqrtC((order-1.0)*((double)order*H-G*G)); // square-root Eq. 9.5.11
    std::complex<T> Gp = G + sq;  // denominator in 9.5.11 with positive sign for square-root
    std::complex<T> Gm = G - sq;  // denominator in 9.5.11 with negative sign for square-root

    // choose Gp or Gm according to which has larger magnitude (page 372, bottom), re-use Gp for 
    // the result:
    double GpA = Gp.getRadius();
    double GmA = Gm.getRadius();
    if( GpA < GmA ) 
    {
      Gp  = Gm;
      GpA = GmA;
    }

    // compute difference between old and new estimate for the root r (the 'a' variable in 
    // Eq. 9.5.8)
    std::complex<T> dr;                       
    if( GpA > 0.0 ) 
      dr = std::complex<T>(order, 0.0) / Gp;  // Eq. 9.5.11
    else
      dr = exp(log(1.0+r.getRadius())) * std::complex<T>(cos((double)i), sin((double)i)); 
        // \todo use sinCos()

    // compute new estimate for the root:
    std::complex<T> rNew = r - dr;   
    if( r == rNew )
      return r;  // converged

    // update our r-variable to the new estimate:
    if( i % itsBeforeFracStep != 0 ) 
      r = rNew;
    else 
      r = r - fractions[i/itsBeforeFracStep]*dr; // fractional step to break limit cycle
  }

  //rsAssert(false);  // error - too many iterations taken, algorithm did not converge
  rsError("Too many iterations taken, algorithm did not converge.");

  return 0.0;
}

template<class T>
void rsPolynomial<T>::findPolynomialRoots(std::complex<T> *a, int order, std::complex<T> *roots)
{
  const double eps = 2.0e-14; // for float, it was 2.0e-6 - use template numeric_limit<T>

  // allocate memory for the coefficients of the deflated polynomial and initialize it as 
  // non-deflated polynomial:
  std::complex<T> *ad = new std::complex<T>[order+1];
  rsCopyBuffer(a, ad, order+1);

  // loop over the roots:
  for(int j = order; j >= 1; j--) 
  {
    // find a root of the deflated polynomial using the Laguerre-algorithm with 0 as initial guess:
    std::complex<T> r = convergeToRootViaLaguerre(ad, j, std::complex<T>(0.0, 0.0));    

    // polish the root by using the Laguerre method with the undeflated polynomial and the 
    // non-polished root as initial guess: 
    r = convergeToRootViaLaguerre(a, order, r);

    // maybe move into a member function Complex::zeroNegligibleImaginaryPart(double ratio); 
    // -> ratio = 2*eps (maybe leave this to client code):
    if( fabs(r.im) <= 2.0*eps*fabs(r.re) ) 
      r.im = 0.0;

    // store root in output array:
    roots[j-1] = r; 

    // deflate the deflated polynomial again by the monomial that corresponds to our most recently 
    // found root:
    std::complex<T> rem = ad[j];  // remainder - not used, needed for passing a dummy pointer
    dividePolynomialByMonomialInPlace(ad, j, r, &rem);
  }

  rsSortComplexArrayByReIm(roots, order);  // maybe leave this to client code
  delete[] ad;
}

template<class T>
void rsPolynomial<T>::findPolynomialRoots(double *a, int order, std::complex<T> *roots)
{
  std::complex<T> *ac = new std::complex<T>[order+1];
  rsConvertBuffer(a, ac, order+1);
  findPolynomialRoots(ac, order, roots);
  delete[] ac;
}

template<class T>
std::vector<std::complex<T>> rsPolynomial<T>::getPolynomialCoefficientsFromRoots(
  std::vector<std::complex<T>> roots)
{
  std::vector<std::complex<T>> coeffs;

  coeffs.ensureAllocatedSize(roots.getNumElements()+1);
  coeffs.appendElement(1.0);

  if( roots.getNumElements() < 1 ) 
    return coeffs;

  for(int i=0; i<roots.getNumElements(); i++)
  {
    std::complex<T> z = roots[i];
    coeffs.appendElement(coeffs[i]);
    for(int j=i; j>=1; j--)
      coeffs[j] = coeffs[j-1] - z * coeffs[j];
    coeffs[0] = -z * coeffs[0];
  }

  return coeffs;
}

template<class T>
void rsPolynomial<T>::rootsToCoeffs(std::complex<T> *r, std::complex<T> *a, int N)
{
  std::complex<T> *rF = new std::complex<T>[N]; // only the finite roots
  int nF = rsCopyFiniteValues(r, rF, N);
  rsFillWithZeros(a, N+1);
  if( nF == 0 )
    a[0] = 1.0;
  else
  {
    a[0] = -rF[0];
    a[1] = 1.0;
    for(int M = 2; M <= nF; M++)
    {
      a[M] = a[M-1];
      std::complex<T> rM = rF[M-1];
      for(int n = M-1; n >= 1; n--)
        a[n] = a[n-1] - rM*a[n];
      a[0] = -rM*a[0];
    }
  }
  delete[] rF;
}

template<class T>
void rsPolynomial<T>::rootsToCoeffs(std::complex<T> *r, double *a, int N)
{
  std::complex<T> *ac = new std::complex<T>[N+1];
  rootsToCoeffs(r, ac, N);
  for(int n = 0; n <= N; n++)
    a[n] = ac[n].re;
  delete[] ac;
}

template<class T>
double rsPolynomial<T>::getRootOfLinearEquation(double a, double b)
{
  if( a == 0.0 )
  {
    RS_DEBUG_BREAK;
    return 0.0;
  }
  else
    return -b/a;
}

template<class T>
std::vector<std::complex<T>> rsPolynomial<T>::getRootsOfQuadraticEquation(
  double a, double b, double c)
{
  // catch degenerate case with zero leading coefficient:
  if( a == 0.0 )
  {
    std::vector<std::complex<T>> roots(1);
    roots[0] = getRootOfLinearEquation(b, c);
    return roots;
  }

  std::vector<std::complex<T>> roots(2); // array to be returned
  double D      = b*b - 4.0*a*c;  // discriminant
  double factor = 1.0 / (2.0*a);  // common factor that appears everywhere
  if( D > 0.0 )
  {
    // D > 0: two distinct real roots:
    double rsSqrt_D = rsSqrt(D);
    roots[0]      = factor * (-b+rsSqrt_D);
    roots[1]      = factor * (-b-rsSqrt_D);
  }
  else if( D == 0.0 )
  {
    // D == 0: a real root with multiplicity 2:
    roots[1] = roots[0] = std::complex<T>( -b * factor );
  }
  else
  {
    // D < 0: two complex conjugate roots:
    double imag = rsSqrt(-D) * factor;
    double real = -b       * factor;
    roots[0]    = std::complex<T>(real,  imag);
    roots[1]    = std::complex<T>(real, -imag);
  }

  return roots;
}

template<class T>
std::vector<std::complex<T>> rsPolynomial<T>::getRootsOfCubicEquation(
  double a, double b, double c, double d)
{
  // catch degenerate cases where the leading coefficient is zero:
  if( a == 0.0 )
    return getRootsOfQuadraticEquation(b, c, d);

  std::vector<std::complex<T>> y(3);
  std::vector<std::complex<T>> roots(3);

  // compute p,q as in the Bronstein page 40, Eq. 1.154c and the offset for the substitution
  // y = x + b/(3*a):
  double p = (3.0*a*c-b*b)/(9.0*a*a);
  double q = (b*b*b)/(27.0*a*a*a) - (b*c)/(6.0*a*a) + d/(2.0*a);

  double u, r, D, phi, ch, sh, re, im, tmp;

  if( p == 0.0 && q == 0.0 )
  {
    y[0] = y[1] = y[2] = 0.0;         // a triple real root at y=0
    // checked with y = (x-1)^3 = x^3-3*x^2+3*x-1
  }
  else if( p != 0.0 && q == 0.0 )
  {
    y[2] = 0.0;                       // a real root at y=0 and ...
    u    = -3.0*p;
    if( u > 0.0 )
    {
      tmp  =  rsSqrt(u);
      y[0] =  tmp;
      y[1] = -tmp;                    // ... two additional real roots or ...
      // checked with y = (x-1)*(x-2)*(x-3) = x^3-6*x^2+11*x-6
    }
    else // u < 0.0
    {
      tmp  =  rsSqrt(-u);
      y[0] =  std::complex<T>(0.0,  tmp);
      y[1] =  std::complex<T>(0.0, -tmp);     // ... two imaginary roots
      // checked with y = (x-4i)*(x+4i)*(x-0) = x^3+16*x
    }
  }
  else if( p == 0.0 && q != 0.0 )
  {
    u = -2.0*q;
    if( u > 0.0 )
    {
      tmp  = pow(u, 1.0/3.0);
      y[2] = tmp;                     // a real root at a positive y or ...
      phi  = (2.0/3.0)*PI;
      // checked with x^3+3*x^2+3*x
    }
    else // u < 0.0
    {
      tmp  = pow(-u, 1.0/3.0);
      y[2] = -tmp;                    // ... a real root at a negative y and ...
      phi  = PI/3.0;
      // checked with x^3+3*x^2+3*x+10
    }
    rsSinCos(phi, &im, &re);
    re  *= tmp;
    im  *= tmp;
    y[0] = std::complex<T>(re,  im);
    y[1] = std::complex<T>(re, -im);           // ... two complex conjugate roots
  }
  else // both p and q are nonzero
  {
    r = rsSign(q) * rsSqrt(fabs(p));
    if( p > 0.0 )
    {
      phi       = rsAsinh( q/(r*r*r) );
      rsSinhCosh(phi/3.0, &sh, &ch);
      y[0] = std::complex<T>(r*sh,  rsSqrt(3.0)*r*ch);
      y[1] = std::complex<T>(r*sh, -rsSqrt(3.0)*r*ch);
      y[2] = -2.0*r*sh;
      // checked with y = (x-i)*(x+i)*(x-1) = x^3-x^2+x-1
    }
    else // p < 0.0
    {
      D = q*q + p*p*p;
      if( D > 0.0 )
      {
        phi  = rsAcosh( q/(r*r*r) );
        rsSinhCosh(phi/3.0, &sh, &ch);
        y[0] = std::complex<T>(r*ch,  rsSqrt(3.0)*r*sh);
        y[1] = std::complex<T>(r*ch, -rsSqrt(3.0)*r*sh);
        y[2] = -2.0*r*ch;
        // checked with y = (x-i)*(x+i)*(x-3) = x^3-3*x^2+x-3
      }
      else // D <= 0.0
      {
        phi  = acos( q/(r*r*r) );
        y[0] =  2.0*r*cos(PI/3.0 + phi/3.0);
        y[1] =  2.0*r*cos(PI/3.0 - phi/3.0);
        y[2] = -2.0*r*cos(         phi/3.0);       // three distinct real roots
        // checked
      }
    }
  }

  // obtain the results for the original equation (back-substitution):
  double s = b/(3.0*a);
  roots[0] = y[0]-s;
  roots[1] = y[1]-s;
  roots[2] = y[2]-s;

  return roots;
}

template<class T>
double rsPolynomial<T>::getCubicRootNear(double x, double a, double b, double c, double d,
                                double min, double max, int maxIterations)
{
  double f, df, xNew;

  f    = ((a*x+b)*x+c)*x+d;
  df   = (3.0*a*x+2.0*b)*x+c;
  xNew = x - f/df;
  int i = 1;
  while( xNew != x && i < maxIterations )
  {
    x    = xNew;
    f    = ((a*x+b)*x+c)*x+d;
    df   = (3.0*a*x+2.0*b)*x+c;
    xNew = x - f/df;
    i++;
  }

  return rsClipToRange(xNew, min, max);
}

template<class T>
double rsPolynomial<T>::getRootNear(double x, double *a, int order, double min, double max,
                          int maxIterations)
{
  // Newton/Raphson iteration:
  double f, df, xNew;
  evaluatePolynomialAndDerivativeAt(x, a, order, &f, &df);
  xNew  = x - f/df;
  int i = 1;
  while( xNew != x && i < maxIterations )
  {
    x    = xNew;
    evaluatePolynomialAndDerivativeAt(x, a, order, &f, &df);
    xNew = x - f/df;
    i++;
  }
  return rsClipToRange(xNew, min, max);
}

template<class T>
void rsPolynomial<T>::cubicCoeffsTwoPointsAndDerivatives(double *a, double *x, double *y, double *dy)
{
  // compute intermediate variables:
  double x0_2 = x[0]*x[0]; // x[0]^2
  double x0_3 = x0_2*x[0]; // x[0]^3
  double x1_2 = x[1]*x[1]; // x[1]^2
  double x1_3 = x1_2*x[1]; // x[1]^3
  double k1   = 3*x[0]*x1_2;
  double k2   = -3*x[1]*y[1];
  double k3   = dy[1]-dy[0];
  double s    = 1/(-x1_3+k1-3*x0_2*x[1]+x0_3);  // scaler

  a[0] =  s*(x0_2*(x1_2*k3+k2) + x0_3*(y[1]-x[1]*dy[1]) + x[0]*x1_3*dy[0] + y[0]*(-x1_3+k1));
  a[1] = -s*(x[0]*(x1_2*(2*dy[1]+dy[0])-6*x[1]*y[1]) - x0_3*dy[1] + x0_2*x[1]*(-dy[1]-2*dy[0])
             + x1_3*dy[0] + 6*x[0]*x[1]*y[0]);
  a[2] =  s*(x[0]*(x[1]*k3-3*y[1]) + x1_2*(dy[1]+2*dy[0]) + x0_2*(-dy[0]-2*dy[1]) + k2 
             + y[0]*(3*x[1]+3*x[0]));
  a[3] = -s*(x[1]*(dy[1]+dy[0]) + x[0]*(-dy[1]-dy[0]) - 2*y[1] + 2*y[0]);
}

template<class T>
void rsPolynomial<T>::rsCubicCoeffsTwoPointsAndDerivatives(double *a, double *y, double *dy)
{
  a[0] = y[0];
  a[1] = dy[0];
  a[2] = 3*(y[1]-a[1]-a[0])-dy[1]+a[1];
  a[3] = (1.0/3.0)*(dy[1]-2*a[2]-a[1]);
}

template<class T>
void rsPolynomial<T>::rsCubicCoeffsFourPoints(double *a, double *y)
{
  a[0] = y[0];
  a[2] = 0.5*(y[-1]+y[1]-2*a[0]);
  a[3] = (1.0/6.0)*(y[2]-y[1]+y[-1]-a[0]-4*a[2]);
  a[1] = 0.5*(y[1]-y[-1]-2*a[3]);
}

template<class T>
double** rsPolynomial<T>::rsVandermondeMatrix(double *x, int N)
{
  double **A; 
  rsAllocateSquareArray2D(A, N);
  for(int i = 0; i < N; i++)
  {
    double xi  = x[i];
    double xij = 1.0;  // xi^j
    for(int j = 0; j < N; j++)
    {
      A[i][j] = xij;
      xij *= xi;
    }
  }
  return A;
}

template<class T>
void rsPolynomial<T>::rsInterpolatingPolynomial(double *a, double *x, double *y, int N)
{
  double **A = rsVandermondeMatrix(x, N);
  rsSolveLinearSystem(A, a, y, N);
  rsDeAllocateSquareArray2D(A, N);  

  // For higher order polynomials, this simple and direct approach may become numerically ill 
  // conditioned. In this case, we could first normalize the data, such than xMin = yMin = -1 and
  // xMax = yMax = +1, find coefficients for Chebychev polynomials, such that our normalized
  // interpolant is given by P(x) = b0*T0(x) + b1*T1(x) + ... + (bN-1)*(TN-1)(x) where Tn is the
  // n-th order Chebychev polynomial. Then, the b-coefficients could be converted back to 
  // a-coefficients for powers of x and finally, we could denormalize using rsShiftPolynomial, 
  // rsStretchPolynomial (for x-denormalization) and scaling the coeffs and adding an offset to
  // a[0] for y-denormalization
}

template<class T>
void rsPolynomial<T>::rsInterpolatingPolynomial(double *a, double x0, double dx, double *y, int N)
{
  double *x = new double[N];
  for(int n = 0; n < N; n++)
    x[n] = x0 + n*dx;
  rsInterpolatingPolynomial(a, x, y, N);
  delete[] x;
}

template<class T>
void rsPolynomial<T>::fitQuadratic(double *a, double *x, double *y)
{
  double k1 = y[1] - y[0];
  double k2 = x[0]*x[0] - x[1]*x[1];
  double k3 = x[1] - x[0];
  double k4 = k1/k3;
  double k5 = k2/k3;

  a[2] = (y[2]-y[0]+k4*(x[0]-x[2])) / (x[2]*(x[2]+k5)-x[0]*(x[0]+k5));
  a[1] = (k1+k2*a[2])/k3;
  a[0] = y[0]-a[2]*x[0]*x[0]-a[1]*x[0];
}

template<class T>
void rsPolynomial<T>::fitQuadratic_0_1_2(double *a, double *y)
{
  a[2] = 0.5*(y[0]+y[2])-y[1];
  a[1] = y[1]-y[0]-a[2];
  a[0] = y[0];
}

template<class T>
void rsPolynomial<T>::fitQuarticWithDerivatives(double *a, double *y, double s0, double s2)
{
  a[0] = y[0];
  a[1] = s0;
  a[2] = -(5*y[2]-16*y[1]+11*y[0]-2*s2+8*s0)/4;
  a[3] =  (7*y[2]-16*y[1]+ 9*y[0]-3*s2+5*s0)/4;
  a[4] = -(2*y[2]- 4*y[1]+ 2*y[0]-  s2  +s0)/4;
}

template<class T>
bool rsPolynomial<T>::areRootsOnOrInsideUnitCircle(double a0, double a1, double a2)
{
  // p and q values for p-q formula
  double p = a1/a2;
  double q = a0/a2;
  double d = p*p/4 - q; // value under the square-root
  double rr, ri;
  if( d < 0 )
  {
    // complex conjugate roots
    rr = -p/2;     // real part
    ri = rsSqrt(-d); // imaginary part
    return rsSqrt(rr*rr + ri*ri) <= 1.0;
  }
  else
  {
    // 2 real roots
    d = rsSqrt(d); 
    rr = -p/2 + d;
    if( fabs(rr) > 1.0 )
      return false;
    rr = -p/2 - d;
    if( fabs(rr) > 1.0 )
      return false;
    return true;
  }
}

template<class T>
void rsPolynomial<T>::besselPolynomial(double *a, int order)
{
  int m, n;
  for(n=0; n<=order; n++)
    a[n] = 0.0;

  if( order == 0 )
  {
    a[0] = 1.0;
    return;
  }
  else if( order == 1 )
  {
    a[0] = 1.0;
    a[1] = 1.0;
    return;
  }
   
  // the general case is treated by recursion:
  a[0]       = 1.0;
  double *b1 = new double[order+1]; 
  double *b2 = new double[order+1]; 
  b2[0]      = 1.0;
  b2[1]      = 0.0;
  b1[0]      = 1.0;
  b1[1]      = 1.0;
  for(n=2; n<=order; n++)
  {
    double c = (double) (2*n-1);
    for(m=n; m>0; m--)
      a[m] = c*b1[m-1];
    a[0] = 0;
    for(m=0; m<n-1; m++)
      a[m] += b2[m];
    for(m=0; m<n; m++)
      b2[m] = b1[m];
    for(m=0; m<=n; m++)
      b1[m] = a[m];
  }
  delete[] b1; 
  delete[] b2; 
}

template<class T>
void rsPolynomial<T>::legendrePolynomial(double *a, int order)
{
  if(order == 0) 
  {
    a[0] = 1.0;
    return;
  }
  if(order == 1) 
  {
    a[0] = 0.0;
    a[1] = 1.0;
    return;
  }

  a[0] = -0.5;
  a[1] =  0.0;
  a[2] =  1.5;
  if(order == 2) 
    return;

  int i, j;
  double *b1 = new double [order+1];
  double *b2 = new double [order+1];

  for(i = 0; i <= order; i++) 
  {
    b1[i] = b2[i] = 0.0;
  }
  b2[1] = 1.0;

  for(i = 3; i <= order; i++) 
  {
    for(j = 0; j <= i; j++) 
    {
      b1[j] = b2[j];
      b2[j] = a[j];
      a[j]  = 0.0;
    }
    for(j = i-2; j >= 0; j-=2) 
    {
      a[j] -= (i-1)*b1[j]/i;
    }
    for(j = i-1; j >= 0; j-=2) 
    {
      a[j+1] += (2*i-1)*b2[j]/i;
    }
  }
  delete [] b1;
  delete [] b2;
}

template<class T>
void rsPolynomial<T>::rsJacobiPolynomialRecursionCoeffs(int n, double a, double b, double *w0, double *w1,
  double *w1x, double *w2)
{
  double k = 2*n+a+b;
  *w0  = 2*n*(n+a+b)*(k-2);
  *w1  = (k-1)*(a*a-b*b);
  *w1x = k*(k-1)*(k-2);
  *w2  = -2*k*(n+a-1)*(n+b-1);
  // from: https://en.wikipedia.org/wiki/Jacobi_polynomials#Recurrence_relation
}

template<class T>
void rsPolynomial<T>::rsJacobiPolynomialRecursion(double *c, int n, double *c1, double *c2, double a, 
  double b)
{
  // initialization:
  if( n == 0 )
  {
    c[0] = 1.0;
    return;
  }
  if( n == 1 )
  {
    c[0] = 0.5*(a-b);
    c[1] = 0.5*(a+b+2);
    return;
  }

  // recursion:
  double w0, w1, w1x, w2;
  rsJacobiPolynomialRecursionCoeffs(n, a, b, &w0, &w1, &w1x, &w2);
  rsPolynomialRecursion(c, w0, n, c1, w1, w1x, c2, w2);
}

template<class T>
void rsPolynomial<T>::rsJacobiPolynomials(double **c, double a, double b, int maxOrder)
{
  for(int n = 0; n <= maxOrder; n++)
    rsJacobiPolynomialRecursion(c[n], n, c[n-1], c[n-2], a, b);
}

template<class T>
void rsPolynomial<T>::rsLegendrePolynomialRecursion(double *a, int n, double *a1, double *a2)
{
  if( n == 0 )
  {
    a[0] = 1.0;
    return;
  }
  if( n == 1 )
  {
    a[0] = 0.0;
    a[1] = 1.0;
    return;
  }
  rsPolynomialRecursion(a, double(n), n, a1, 0.0, 2.0*n-1.0, a2, -(n-1.0));

  // Legendre polynomials are a special case of Jacobi polynomials, so this would also work:
  // rsJacobiPolynomialRecursion(a, n, a1, a2, 0.0, 0.0); 
}

template<class T>
void rsPolynomial<T>::rsPartialFractionExpansion(
  std::complex<T> *numerator,   int numeratorOrder,
  std::complex<T> *denominator, int denominatorOrder, 
  std::complex<T> *poles, int *multiplicities, int numDistinctPoles,
  std::complex<T> *pfeCoeffs)
{
  // sanity check for inputs:
  rsAssert(numeratorOrder < denominatorOrder);
  rsAssert(rsSum(multiplicities, numDistinctPoles) == denominatorOrder);

  // make denominator monic:
  rsScale(numerator,   numeratorOrder+1,   1.0/denominator[denominatorOrder]);
  rsScale(denominator, denominatorOrder+1, 1.0/denominator[denominatorOrder]);

  // establish coefficient matrix:
  std::complex<T> **A; rsAllocateSquareArray2D(A, denominatorOrder);
  std::complex<T> *tmp = new std::complex<T>[denominatorOrder+1]; // deflated denominator
  std::complex<T> remainder;                                   // always zero
  for(int i = 0, k = 0; i < numDistinctPoles; i++)
  {
    rsCopyBuffer(denominator, tmp, denominatorOrder+1);
    for(int m = 0; m < multiplicities[i]; m++)
    {
      dividePolynomialByMonomialInPlace(tmp, denominatorOrder-m, poles[i], &remainder);
      for(int j = 0; j < denominatorOrder; j++)
        A[j][k] = tmp[j];
      k++;
    }
  }

  // solve the linear system using an appropriately zero-padded numerator as RHS:
  rsCopyBuffer(numerator, tmp, numeratorOrder+1);
  rsFillWithZeros(&tmp[numeratorOrder+1], denominatorOrder-(numeratorOrder+1));
  rsSolveLinearSystem(A, pfeCoeffs, tmp, denominatorOrder);

  // clean up:
  rsDeAllocateSquareArray2D(A, denominatorOrder);  
  delete[] tmp;
}


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
   order of new denominator = order(S)*order(Q)
   order of new numerator = ?
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
  (we should have a function "reduce" for that)
 -or maybe it should first find the roots

*/

