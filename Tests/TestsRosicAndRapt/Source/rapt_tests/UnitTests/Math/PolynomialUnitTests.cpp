
using namespace RAPT;

/*
  template <class T>
  void rsDeConvolve(T *y, int yLength, T *h, int hLength, T *x)
  {
    int m = rsFirstIndexWithNonZeroValue(h, hLength);
    if( m == -1 )
    {
      // h is all zeros - return an all-zero x-signal:
      rsFillWithZeros(x, yLength-hLength+1);
      return;
    }
    T scaler = T(1) / h[m];
    x[0]     = scaler * y[m];
    for(int n = 1; n < yLength-hLength+1; n++)
    {
      x[n] = y[n+m];
      for(int k = m+1; k <= rsMin(hLength-1, n+m); k++)
        x[n] -= h[k] * x[n-k+m];
      x[n] *= scaler;
    }
  }
*/


  /** Given the sequence y of length yLength, this function returns a sequence x which, when
  convolved with itself, gives y. yLength is assumed to be odd, the first nonzero value in y is
  assumed to be positive and the index of first nonzero value is assumed to be even (because this
  is what will happen, when covolving a sequence with itself). To disambiguate the square-root, the
  function will return a sequence with its 1st nonzero value being positive. If the original
  sequence x (before it was convolved with itself to give y) started with a negative value, the
  result of taking the square-root of the squared sequence y will have all signs reversed with
  respect to the original sequence x. The length of x will be (yLength+1)/2. */
/*
  template <class T>
  void rsSequenceSqrt(T *y, int yLength, T *x)
  {
    int m2 = rsFirstIndexWithNonZeroValue(y, yLength);
    if( m2 == -1 )
    {
      // y is all zeros - return an all-zero x-sequence:
      rsFillWithZeros(x, (yLength+1)/2);
      return;
    }
    int m = m2/2;              // 1st index of nonzero value in x
    rsFillWithZeros(x, m);
    x[m] = rsSqrt(y[m2]);
    double scaler = 1.0 / (2*x[m]);
    for(int n = m+1; n < (yLength+1)/2; n++)
    {
      x[n] = y[m+n];
      for(int k = 1; k <= n-m-1; k++)
        x[n] -= x[m+k] * x[n-k];
      x[n] *= scaler;
    }
  }
  */

  /*
  template <class T>
  void rsSequenceSqrt(T *y, int yLength, T *x)
  {
    x[0] = rsSqrt(y[0]);
    double scaler = 1.0 / (2*x[0]);
    for(int n = 1; n < (yLength+1)/2; n++)
    {
      x[n] = y[n];
      for(int k = 1; k <= n-1; k++)
        x[n] -= x[k]*x[n-k];
      x[n] *= scaler;
    }
    // \todo - include a means to deal with sequences that have leading zeros
  }
  */


bool testConvolution(std::string &reportString)
{
  std::string testName = "(De)Convolution";
  bool testResult = true;

  static const int xN = 10;
  static const int hN = 5;
  static const int yN = xN + hN - 1;

  double x[xN]  = {1,2,-2,3,-1,-3,-5,4,-3,-1};              // input sequence
  double h[hN]  = {2,-3,1,2,-1};                            // impulse response
  double yt[yN] = {2,1,-9,16,-10,-6,6,15,-28,4,13,-11,1,1}; // target output sequence
  double y[yN];                                             // output sequence

  // test algorithm when all pointers are distinct:
  rsArrayTools::convolve(x, xN, h, hN, y);
  testResult &= rsArrayTools::equal(y, yt, yN);

  // test in-place convolution where x == y:
  rsArrayTools::fillWithZeros(y, yN);
  rsArrayTools::copy(x, y, xN);
  rsArrayTools::convolve(y, xN, h, hN, y);
  testResult &= rsArrayTools::equal(y, yt, yN);

  // test in-place convolution where h == y:
  rsArrayTools::fillWithZeros(y, yN);
  rsArrayTools::copy(h, y, hN);
  rsArrayTools::convolve(x, xN, y, hN, y);
  testResult &= rsArrayTools::equal(y, yt, yN);

  // test in-place convolution where x == h == y:
  rsArrayTools::fillWithZeros(y, yN);
  rsArrayTools::copy(h, y, hN);
  rsArrayTools::convolve(y, xN, y, hN, y);
  testResult &= y[0]  ==   4;
  testResult &= y[1]  == -12;
  testResult &= y[2]  ==  13;
  testResult &= y[3]  ==   2;
  testResult &= y[4]  == -15;
  testResult &= y[5]  ==  10;
  testResult &= y[6]  ==   2;
  testResult &= y[7]  ==  -4;
  testResult &= y[8]  ==   1;
  testResult &= y[9]  ==   0;
  testResult &= y[10] ==   0;
  testResult &= y[11] ==   0;
  testResult &= y[12] ==   0;
  testResult &= y[13] ==   0;

  // test deconvolution - recover the signal x:
  rsArrayTools::convolve(x, xN, h, hN, y);
  double xx[xN];
  rsArrayTools::deConvolve(y, yN, h, hN, xx);
  testResult &= rsArrayTools::almostEqual(x, xx, xN, 1.e-13);

  // convolve and deconvolve with an impulse response with leading zeros:
  h[0] = 0.0;
  h[1] = 0.0;
  rsArrayTools::convolve(x, xN, h, hN, y);
  rsArrayTools::deConvolve(y, yN, h, hN, xx);
  testResult &= rsArrayTools::almostEqual(x, xx, xN, 1.e-13);

  // recover the impulse response h:
  double hh[hN];
  rsArrayTools::deConvolve(y, yN, x, xN, hh);
  testResult &= rsArrayTools::almostEqual(h, hh, hN, 1.e-13);

  // test (de)convolution with all-zero impulse response:
  rsArrayTools::fillWithZeros(h, hN);
  rsArrayTools::convolve(x, xN, h, hN, y);
  rsArrayTools::deConvolve(y, yN, h, hN, xx);
  testResult &= rsArrayTools::isAllZeros(xx, xN);

  // test "square-root" of a sequence - convolve h with itself and recover h from the convolved
  // sequence:
  h[0]=2; h[1]=-3; h[2]=1; h[3]=2; h[4]=-1; // because we messed with it
  rsArrayTools::convolve(h, hN, h, hN, y);       // h convolved with itself ("h^2")
  int h2N = 2*hN-1;                         // length of h^2
  rsArrayTools::fillWithZeros(hh, hN);
  rsArrayTools::sequenceSqrt(y, h2N, hh);
  testResult &= rsArrayTools::almostEqual(h, hh, hN, 1.e-13);

  // test sequence square-root, when the sequence has leading zeros:
  h[0]=0; h[1]=0; h[2]=4; h[3]=-8; h[4]=2;
  rsArrayTools::convolve(h, hN, h, hN, y);
  rsArrayTools::sequenceSqrt(y, h2N, hh);
  testResult &= rsArrayTools::almostEqual(h, hh, hN, 1.e-13);

  // if we try to take the square-root x of an arbitrary sequence y (which was not constructed
  // by squaring some given sequence), and convolve the computed square-root with itself again,
  // in an attempt to reconstruct y, the result matches y only up to the n-th term where n is the
  // length of the square-root sequence. this is not surprising since in computing the square-root
  // sequence, we will not use any y[k] for k > n:
  y[0]=+0.2; y[1]=-0.3; y[2]=-0.1; y[3]=+0.4; y[4]=+0.2;
  y[5]=-0.3; y[6]=-0.2; y[7]=-0.4; y[8]=-0.2; y[9]=+0.3; y[10]=0.2;
  rsArrayTools::sequenceSqrt(y, 11, x);
  double yy[11];
  rsArrayTools::convolve(x, 6, x, 6, yy);
  testResult &= rsArrayTools::almostEqual(y, yy, 6, 1.e-13);

  return testResult;
}


bool testCubicCoeffsFourPoints(std::string &reportString)
{
  std::string testName = "CubicCoeffsFourPoints";
  bool testResult = true;

  double  y[4] = {3, -2, 5, 1};
  double  a[4];

  rsPolynomial<double>::cubicCoeffsFourPoints(a, &y[1]);

  double yc;            // computed value
  double tol = 1.e-14;  // tolerance

  yc = rsPolynomial<double>::evaluate(-1.0, a, 3);
  testResult &= rsIsCloseTo(yc,  y[0], tol);
  yc = rsPolynomial<double>::evaluate(0.0, a, 3);
  testResult &= rsIsCloseTo(yc,  y[1], tol);
  yc = rsPolynomial<double>::evaluate(1.0, a, 3);
  testResult &= rsIsCloseTo(yc,  y[2], tol);
  yc = rsPolynomial<double>::evaluate(2.0, a, 3);
  testResult &= rsIsCloseTo(yc,  y[3], tol);

  return testResult;
}

bool testCubicCoeffsTwoPointsAndDerivatives(std::string &reportString)
{
  std::string testName = "CubicCoeffsTwoPointsAndDerivatives";
  bool testResult = true;

  double  x[2] = {-3, 2};
  double  y[2] = {-2, 5};
  double dy[2] = { 7, 3};
  double  a[4];

  rsPolynomial<double>::cubicCoeffsTwoPointsAndDerivatives(a, x, y, dy);

  // check results:
  double yc, dyc;       // computed values
  double tol = 1.e-14;  // tolerance

  rsPolynomial<double>::evaluateWithDerivative(x[0], a, 3, &yc, &dyc);
  testResult &= rsIsCloseTo( yc,  y[0], tol);
  testResult &= rsIsCloseTo(dyc, dy[0], tol);

  rsPolynomial<double>::evaluateWithDerivative(x[1], a, 3, &yc, &dyc);
  testResult &= rsIsCloseTo( yc,  y[1], tol);
  testResult &= rsIsCloseTo(dyc, dy[1], tol);

  // test, if in the special case where x1=0, x2=1:
  x[0] = 0.0;
  x[1] = 1.0;
  rsPolynomial<double>::cubicCoeffsTwoPointsAndDerivatives(a, x, y, dy);
  rsPolynomial<double>::evaluateWithDerivative(x[0], a, 3, &yc, &dyc);
  testResult &= rsIsCloseTo( yc,  y[0], tol);
  testResult &= rsIsCloseTo(dyc, dy[0], tol);
  rsPolynomial<double>::evaluateWithDerivative(x[1], a, 3, &yc, &dyc);
  testResult &= rsIsCloseTo( yc,  y[1], tol);
  testResult &= rsIsCloseTo(dyc, dy[1], tol);

  // test, if the simplified algorithm for the special case returns the same coeffs:
  double b[4];
  rsPolynomial<double>::cubicCoeffsTwoPointsAndDerivatives(b, y, dy);
  testResult &= rsIsCloseTo(a[0], b[0], tol);
  testResult &= rsIsCloseTo(a[1], b[1], tol);
  testResult &= rsIsCloseTo(a[2], b[2], tol);
  testResult &= rsIsCloseTo(a[3], b[3], tol);

  return testResult;
}


bool testPolynomialEvaluation(std::string &reportString)
{
  std::string testName = "PolynomialEvaluation";
  bool testResult = true;

  // establish 5th order polynomial and its 1st 3 derivatives:
  double a[6]  = {2,-1,5,7,-3,2};   // a(x) = 2*x^5 - 3*x^4 + 7*x^3 + 5*x^2 - 1*x^1 + 2*x^0
  double a1[5] = {-1,10,21,-12,10}; // a'(x) = 10*x^4 - 12*x^3 + 21*x^2 + 10*x^1 - 1*x^0
  double a2[4] = {10,42,-36,40};    // a''(x) = 40*x^3 - 36*x^2 42*x^1 + 10*x^0
  double a3[3] = {42,-72,120};      // a'''(x) = 120*x^2 - 72*x^1 + 42*x^0

  double x0 = 2.0;     // point, where to evaluate
  double y[4], yt[4];  // evaluation results and target values

  yt[0] = rsPolynomial<double>::evaluate(x0, a,  5);
  yt[1] = rsPolynomial<double>::evaluate(x0, a1, 4);
  yt[2] = rsPolynomial<double>::evaluate(x0, a2, 3);
  yt[3] = rsPolynomial<double>::evaluate(x0, a3, 2);

  // test evaluation of polynomial and 1st derivative:
  rsPolynomial<double>::evaluateWithDerivative(x0, a, 5, &y[0], &y[1]);
  testResult &= yt[0] == y[0];
  testResult &= yt[1] == y[1];

  // test evaluation of polynomial and 1st 3 derivatives:
  rsPolynomial<double>::evaluateWithDerivatives(x0, a, 5, y, 3);
  testResult &= yt[0] == y[0];
  testResult &= yt[1] == y[1];
  testResult &= yt[2] == y[2];
  testResult &= yt[3] == y[3];

  return testResult;
}


bool testPolynomialDivision(std::string &reportString)
{
  std::string testName = "PolynomialDivision";
  bool testResult = true;

  // given a polynomial p and a divisor polynomial d, polynomial division finds two polynomials
  // q (quotient) and r (remainder) such that p = d*q + r

  // establish denominator d, quotient q and remainder r:
  double d[6] = {2, -1, 5,  7, -3, 2};  // d(x) = 2*x^5 - 3*x^4 + 7*x^3 + 5*x^2 - 1*x^1 + 2*x^0
  double q[4] = {2, -3, 6,  2};         // q(x) =                 2*x^3 + 6*x^2 - 3*x^1 + 2*x^0
  double r[5] = {3,  1, 4, -5, 3};      // r(x) =         3*x^4 - 5*x^3 + 4*x^2 + 1*x^1 - 3*x^0

  // establish polynomial p(x) = d(x)*q(x) + r(x):
  double p[9];                                          // 8th degree, 9 coeffs
  rsPolynomial<double>::multiply(d, 5, q, 3, p);               // p(x) = d(x)*q(x)
  rsPolynomial<double>::weightedSum(p, 8, 1.0, r, 4, 1.0, p);  // p(x) = d(x)*q(x) + r(x);

  // retrieve q(x) and r(x):
  double qq[9], rr[9];
  rsPolynomial<double>::divide(p, 8, d, 5, qq, rr);

  // p(x)/d(x) = q(x) + r(x)/d(x)

  testResult &= rsArrayTools::equal(q, qq, 4);
  testResult &= rsArrayTools::equal(r, rr, 5);

  return testResult;
}

bool testPolynomialArgumentShift(std::string &reportString)
{
  std::string testName = "PolynomialArgumentShift";
  bool testResult = true;

  // define polynomial p(x) and the shift value x0:
  static const int order = 6;
  double p[order+1]  = {2,1,-5,7,-3,2,-2}; // p(x) = -2x^6+2x^5-3x^4+7x^3-5x^2+1x^1+2x^0
  double x0          = 2.0;                // shift value

  // establish coeffs of q(x) = p(x-x0):
  double q[order+1];
  rsPolynomial<double>::shiftArgument(p, q, order, x0);

  // check, if q-coeffs have correct values:
  testResult &= q[0] == -316;
  testResult &= q[1] ==  745;
  testResult &= q[2] == -759;
  testResult &= q[3] ==  431;
  testResult &= q[4] == -143;
  testResult &= q[5] ==   26;
  testResult &= q[6] ==   -2;

  return testResult;
}

bool testPolynomialDiffAndInt(std::string &reportString)
{
  std::string testName = "PolynomialDiffAndInt";
  bool testResult = true;

  double a[6]  = {2, -1, 5, 7, -3, 2};
  double ad[5];
  double ai[7];

  using Poly = rsPolynomial<double>;

  Poly::derivative(a, ad, 5);
  testResult &= (ad[0] == -1);
  testResult &= (ad[1] == 10);
  testResult &= (ad[2] == 21);
  testResult &= (ad[3] == -12);
  testResult &= (ad[4] == 10);

  Poly::integral(a, ai, 5, 2.0);
  testResult &= (ai[0] ==  2.0);
  testResult &= (ai[1] ==  2.0/1.0);
  testResult &= (ai[2] == -1.0/2.0);
  testResult &= (ai[3] ==  5.0/3.0);
  testResult &= (ai[4] ==  7.0/4.0);
  testResult &= (ai[5] == -3.0/5.0);
  testResult &= (ai[6] ==  2.0/6.0);

  Poly p;           // should be the zero polynomial and have 1 coeff which is zero
  p.shiftY(1);      // now the coeff should be one
  p.integrate(1.0);
  testResult &= p == Poly({1,1});
  p.integrate(1.0);
  testResult &= p == Poly({1,1,1./2});
  p.integrate(1.0);
  testResult &= p == Poly({1,1,1./2,1./6});
  p.integrate(1.0);
  testResult &= p == Poly({1,1,1./2,1./6,1./24});

  // Test storing polynomials in arrays and integrating the array elements:
  std::vector<Poly> polys;
  polys.push_back(p);
  polys.push_back(p);
  polys.push_back(p);
  polys[1].integrate(1.0);
  testResult &= polys[1] == Poly({1,1,1./2,1./6,1./24,1./120});



  return testResult;
}


  template <class T>
  T evaluateFallingFactorialAt(T x, int n, T h = 1)
  {
    T y = T(0);
    for(int k = 0; k < n; k++)
      y *= (x - k*h);
    return y;
  }
  template <class T>
  T evaluateFallingFactorialPolynomialAt(T x, T *a, int order, T h = 1)
  {
    T y = T(0);
    for(int n = 0; n <= order; n++)
      y += a[n] * evaluateFallingFactorialAt(x, n, h);

    // the above is naive and preliminary - it can be optimized like so (converts O(N^2) to O(N)):
    /*
    T y  = a[0];
    T xn = x;
    for(int n = 1; n <= order; n++)
    {
      y  += a[n] * xn;
      xn *= (x-n*h);
    }
    */

    return y;

  }
  // \todo: continue here - implement evaluation of factorial polynomials (falling and rising)
  // implement conversion of polynomial coefficients from to falling/rising factorial polynomial
  // representations via striling numbers - maybe we can generalize the stirling numbers to include
  // a stepsize (the standard stirling numbers assume h = 1)
  // write unit tests for all this stuff

bool testPolynomialFiniteDifference(std::string &reportString)
{
  std::string testName = "PolynomialFiniteDifference";
  bool testResult = true;

  static const int order = 5;
  //static const int numCoeffs = order+1;
  double a[6]  = {2, -1, 5, 7, -3, 2};
  double ad[5];  // 1st difference polynomial

  static const int numValues = 20;
  double h    = 0.25; // stepsize;
  double xMin = -2.0;
  double x[numValues], y[numValues];
  double yf[numValues], yfc[numValues];  // forward differences (true and computed)
  double yb[numValues], ybc[numValues];  // backward differences (true and computed)

  // create the reference data to match:
  int n;
  for(n = 0; n < numValues; n++)
  {
    x[n] = xMin + n*h;
    y[n] = rsPolynomial<double>::evaluate(x[n], a, order);
  }
  for(n = 0; n < numValues-1; n++)
    yf[n] = y[n+1] - y[n];
  yf[numValues-1] = 0;      // not existent, actually

  for(n = 1; n < numValues; n++)
    yb[n] = y[n] - y[n-1];
  yb[0] = 0;               // not existent, actually

  // check forward difference:
  rsPolynomial<double>::finiteDifference(a, ad, order, 1, h);
  for(n = 0; n < numValues; n++)
    yfc[n] = rsPolynomial<double>::evaluate(x[n], ad, order-1);
  testResult &= rsArrayTools::equal(yf, yfc, numValues-1);

  // check backward difference:
  rsPolynomial<double>::finiteDifference(a, ad, order, -1, h);
  for(n = 0; n < numValues; n++)
    ybc[n] = rsPolynomial<double>::evaluate(x[n], ad, order-1);
  testResult &= rsArrayTools::equal(&yb[1], &ybc[1], numValues-1);

  return testResult;
}

bool testPolynomialComposition(std::string &reportString)
{
  std::string testName = "PolynomialComposition";
  bool testResult = true;

  using Poly = rsPolynomial<double>;

  static const int na = 5;
  static const int nb = 4;
  static const int nc = na*nb;
  double a[na+1] = {2, -1, 5,  7, -3, 2}; // 2*x^5 - 3*x^4 + 7*x^3 + 5*x^2 - 1*x^1 + 2*x^0
  double b[nb+1] = {3,  1, 4, -5,  3};    //         3*x^4 - 5*x^3 + 4*x^2 + 1*x^1 - 3*x^0
  double c[nc+1];
  Poly::compose(a, na, b, nb, c);

  // check, if the composed c-polynomial returns the same result as applying the 2nd b-polynomial
  // to the result of the 1st a-polynomial:
  double x = -3.0; // input value
  double y1, y2;
  y1 = Poly::evaluate(x,  a, na);
  y1 = Poly::evaluate(y1, b, nb);
  y2 = Poly::evaluate(x,  c, nc);
  testResult &= (y1 == y2);

  Poly p({2, -1, 5,  7, -3, 2});
  Poly q({3,  1, 4, -5,  3});
  Poly pq = p(q);
  // Sage:
  // p(x) = 2 - 1*x + 5*x^2 + 7*x^3 - 3*x^4 + 2*x^5
  // q(x) = 3 + 1*x + 4*x^2 - 5*x^3 + 3*x^4
  // p(q).expand()
  testResult &= pq == Poly({476,704,3262,199,6627,-9747,16849,-26397,46908,-61886,82425,-106821,
                            124440,-126114,120574,-106210,78177,-43290,16740,-4050,486});

  return testResult;
}

bool testPolynomialWeightedSum(std::string &reportString)
{
  std::string testName = "PolynomialWeightedSum";
  bool testResult = true;

  static const int pN = 5;
  static const int qN = 3;
  static const int rN = pN;              // == max(pN, qN);
  double p[pN+1] = {3, -1, 5, 7, -3, 2}; // p(x) = 2*x^5 - 3*x^4 + 7*x^3 + 5*x^2 - 1*x^1 + 3*x^0
  double q[qN+1] = {2, -3, 6, 2};        // q(x) =                 2*x^3 + 6*x^2 - 3*x^1 + 2*x^0
  double r[rN+1];

  // r(x) = 2*p(x) + 3*q(x) = 4*x^5 - 6*x^4 + 20*x^3 + 28*x^2 - 11*x^1 + 12*x^0
  rsPolynomial<double>::weightedSum(p, pN, 2.0, q, qN, 3.0, r);
  testResult &= (r[0] ==  12);
  testResult &= (r[1] == -11);
  testResult &= (r[2] ==  28);
  testResult &= (r[3] ==  20);
  testResult &= (r[4] == - 6);
  testResult &= (r[5] ==   4);

  // exchange roles of function parameters (function takes the other branch, result should be the
  // same):
  rsPolynomial<double>::weightedSum(q, qN, 3.0, p, pN, 2.0, r);
  testResult &= (r[0] ==  12);
  testResult &= (r[1] == -11);
  testResult &= (r[2] ==  28);
  testResult &= (r[3] ==  20);
  testResult &= (r[4] == - 6);
  testResult &= (r[5] ==   4);

  // use a truncated polynomial for p (such that p and q are of the same order):
  rsArrayTools::fillWithZeros(r, rN+1);
  rsPolynomial<double>::weightedSum(p, qN, 2.0, q, qN, 3.0, r);
  testResult &= (r[0] ==  12);
  testResult &= (r[1] == -11);
  testResult &= (r[2] ==  28);
  testResult &= (r[3] ==  20);
  testResult &= (r[4] ==   0);
  testResult &= (r[5] ==   0);

  return testResult;
}

bool testPolynomialIntegrationWithPolynomialLimits(std::string &reportString)
{
  std::string testName = "PolynomialIntegrationWithPolynomialLimits";
  bool testResult = true;

  static const int np = 5;
  static const int na = 3;
  static const int nb = 4;
  static const int nP = np+1;
  static const int nq = nP*nb;  // == nP*rmax(na, nb);

  double p[np+1] = {2, -1, 5,  7, -3, 2};  // 2*x^5 - 3*x^4 + 7*x^3 + 5*x^2 - 1*x^1 + 2*x^0
  double a[na+1] = {2, -3, 6,  2};         //                 2*x^3 + 6*x^2 - 3*x^1 + 2*x^0
  double b[nb+1] = {3,  1, 4, -5, 3};      //         3*x^4 - 5*x^3 + 4*x^2 + 1*x^1 - 3*x^0
  double q[nq+1];

  double x = 1.5; // input value
  double y1, y2;

  // obtain coefficients of indefinite integral:
  double P[nP+1];
  rsPolynomial<double>::integral(p, P, np);

  // compute integration limits for definite integral:
  double lowerLimit = rsPolynomial<double>::evaluate(x, a, na);
  double upperLimit = rsPolynomial<double>::evaluate(x, b, nb);

  // evaluate definite integral:
  y1 = rsPolynomial<double>::evaluate(upperLimit, P, nP) - rsPolynomial<double>::evaluate(lowerLimit, P, nP);


  rsPolynomial<double>::integrateWithPolynomialLimits(p, np, a, na, b, nb, q);
  y2 = rsPolynomial<double>::evaluate(x, q, nq);

  testResult &= rsIsCloseTo(y2, y1, 1.e-13 * fabs(y1));

  return testResult;
}

bool testPolynomialInterpolation(std::string &reportString)
{
  std::string testName = "PolynomialInterpolation";
  bool testResult = true;

  double tol = 1.e-13; // error tolerance

  using Poly = rsPolynomial<double>;

  // establish dataset to interpolate:
  static const int N = 7;  // N is number of data points, not polynomial degree!
  double x[N] = {-0.5, 1.0, 0.7, 1.5, -2.0, -1.3, 2.2};
  double y[N] = { 1.2, 1.4, 0.2, 2.5, -1.7, -2.3, 1.3};

  // get polynomial coefficients:
  double a[N];
  Poly::interpolant(a, x, y, N);

  // check, if the polynomial really matches the data:
  double yc[N];
  int n;
  for(n = 0; n < N; n++) {
    yc[n] = Poly::evaluate(x[n], a, N-1);
    testResult &= rsIsCloseTo(yc[n], y[n], tol); }

  // test computing coeffs for Newton basis polynomials and evaluation of Newton polynomial:
  double c[N];
  rsArrayTools::copy(y, c, N);
  Poly::coeffsNewton(x, c, N);
  for(n = 0; n < N; n++) {
    yc[n] = Poly::evaluateNewton(x[n], c, x, N-1);
    testResult &= rsIsCloseTo(yc[n], y[n], tol); }

  // test in-place implementation of finding polynomial coeffs via Newton expansion:
  Poly::interpolantViaNewton(a, x, y, N, c); // after the call, a will contain the desired coeffs 
  for(n = 0; n < N; n++) {
    yc[n] = Poly::evaluate(x[n], a, N-1);
    testResult &= rsIsCloseTo(yc[n], y[n], tol); }

  // check, if the workspace c contains coeffs of polynomial with roots at the original x[n] - but
  // only the values up to x[N-2] are roots - the last original x is not a root (this is weird - 
  // but obviously, the polynomial can only have N-1 roots because it's of degree N-1):
  for(n = 0; n < N-1; n++) {   // loop only up to N-1 because x[N-1] is not a root of c
    yc[n] = Poly::evaluate(x[n], c, N-1);
    testResult &= rsIsCloseTo(yc[n], 0.0, tol); }
  // what are the implications of the fact that the last x is not a root - does it mean, the order
  // of stored x-values (and associated y-values) can have an impact on the roundoff error? i 
  // think so...it probably has anyway

  // test function for equidistant abscissa values:
  double x0 = -3.2;
  double dx =  1.1;
  Poly::interpolant(a, x0, dx, y, N);
  for(n = 0; n < N; n++) {
    yc[n] = Poly::evaluate(x0+n*dx, a, N-1);
    testResult &= rsIsCloseTo(yc[n], y[n], tol); }

  return testResult;
}
// todo: test the different algorithms (Vandermonde-matrix, Lagrange, Newton) to see, if they give 
// the same results

bool testPolynomialRootFinder(std::string &reportString)
{
  std::string testName = "PolynomialRootFinder";
  bool testResult = true;

  // we use the polynomial p(x) = x^4 - 7x^3 + 21*x^2 - 23*x - 52 with roots at 2+3i, 2-3i, -1, 4 
  // as test function:
  double a1[5] = {-52, -23, 21, -7, 1};
  std::complex<double> r1[4];
  rsPolynomial<double>::roots(a1, 4, r1);

  // now we 
  static const int maxN     = 20;
  static const int numTests = 1000;
  double range = 10.0;                // range for the real and imaginary parts of the roots
  double tol   = 5.e-8;               // tolerance
  std::complex<double> a[maxN+1];     // polynomial coefficients
  std::complex<double> rTrue[maxN];   // true roots
  std::complex<double> rFound[maxN];  // roots that were found
  rsRandomUniform(-range, range, 0);  // set seed
  int i, j, k;
  for(i = 1; i <= numTests; i++) {
    // polynomial degree for this test:
    int N = (int) rsRandomUniform(1.0, maxN);

    // generate a bunch of random roots:
    for(k = 0; k < N; k++) {
      rTrue[k].real(rsRandomUniform(-range, range));
      rTrue[k].imag(rsRandomUniform(-range, range)); }

    // obtain polynomial coeffs:
    rsPolynomial<double>::rootsToCoeffs(rTrue, a, N);

    // find the roots:
    rsPolynomial<double>::roots(a, N, rFound);

    // try to find a matching root in the found roots for each of the true roots:
    for(j = 0; j < N; j++) {
      bool matchFound = false;
      for(k = 0; k < N; k++) {
        if( abs(rFound[j]-rTrue[k]) < tol ) {
          matchFound = true; break; }}
      rsAssert(matchFound);
      testResult &= matchFound; 
    }
  }

  // we need a rather high tolerance - the precision of the root finding algorithm seems to be not 
  // very good - can it be improved? Setting tol to 1.e-9 triggers the assert - and it
  // maybe by performing one or two steps of newton iteration...but
  // actually, the algo already does this "polishing" thing


  // try with p(x) = 27 + 9*x - 3*x^2 - 1*x^3 = -(x+3)^2 * (x-3). this has simple root at +3 and a 
  // double root at -3
  a[0] = 27, a[1] = 9, a[2] = -3, a[3] = -1;
  rsPolynomial<double>::roots(a, 3, rFound);
  // the double-root at -3 is found as: -3.0000000215539959, -2.9999999784460041 - this seems like
  // a very bad precision - can it be improved by more/better polishing? what happens, if we try
  // Netwon iteration? it will probably not work, because this is a double-root:
  // https://www.wolframalpha.com/input/?i=27+%2B+9*x+-+3*x%5E2+-+1*x%5E3
  // maybe in case of double roots, we could polish by finding a root of the derivative? could it 
  // be that multiple roots generally produce imprecise results? But in our random polynomials 
  // above, we probably don't get any multiple roots - the results seem to be imprecise 
  // nonetheless. Maybe, within the Newton-iteration, we should detect, whether we have a double 
  // (or triple, whatever multiple) root and if this condition is detected - find the root of the 
  // derivative. 
  // using numpy (can be done in sagecell):
  //
  // import numpy as np
  // r = np.roots([-1,-3,9,27])
  // r, format(r[1], '.16g')
  //
  // we get results in which the real parts are more precise but which have some false imaginary
  // part of the order of 1.e-8 (maybe we don't get any imaginary parts here because the flush it 
  // to zero based on a threshold? -> figure out). the numpy doc 
  // https://docs.scipy.org/doc/numpy/reference/generated/numpy.roots.html says:
  // "The algorithm relies on computing the eigenvalues of the companion matrix" which could mean 
  // that numpy uses the Jenkins/Traub algorithm (which seems to be the most popular polynomial
  // root findning algorithm anyway). one possible improvement could be to post-polish all roots 
  // via newton-iteration with multiple-root-detection and using the derivative in such a detected
  // case. detection could be based on the absolute value of the derivative (or the ratio f/f')
  // or maybe the ratio of the absolute value of the step-size and the absolute value of the 
  // current estimate - steps are supposed to be not too large - but no - roots near zero do not
  // seem to warrant smaller thresholds






  return testResult;
}

bool testPartialFractionExpansion(std::string &reportString)
{
  std::string testName = "PartialFractionExpansion";
  bool testResult = true;

  typedef std::complex<double> Complex;
  typedef rsRationalFunction<double> RF;

  // local variables:
  Complex j(0.0, 1.0);      // imaginary unit
  Complex n[9], d[9], p[9]; // numerator, denominator, poles
  Complex a[9];             // partial fraction expansion coefficients
  int numPoles;             // number of distinc poles
  int m[9];                 // pole multiplicities
  int numDeg, denDeg;       // numerator and denominator degrees
  double tol = 1.e-14;      // tolerance for equality checks

  // 1st test - contains a double root:
  // f(x) = 3/(x+5) - 4/(x+3) + 2/(x-1) + 5/(x-1)^2 - 3/(x-5)
  //      = n(x)/d(x) = numerator(x) / denominator(x)
  //      = (-2*x^4-13*x^3+25*x^2-275*x-215)/(x^5+x^4-30*x^3-22*x^2+125*x-75)
  numDeg   = 4;  // numerator degree
  denDeg   = 5; 
  numPoles = 4;  // only 4, because the 3rd one is a double-root
  n[0] = -215; n[1] = -275; n[2] = 25;  n[3] = -13; n[4] = -2;
  d[0] = -75;  d[1] = 125;  d[2] = -22; d[3] = -30; d[4] =  1; d[5] = 1;
  p[0] = -5; p[1] = -3; p[2] = +1; p[3] =  5;
  m[0] =  1; m[1] =  1; m[2] =  2; m[3] =  1;
  RF::partialFractionExpansion(n, numDeg, d, denDeg, p, m, numPoles, a);
  testResult &= abs(a[0]-3.0) < tol;
  testResult &= abs(a[1]+4.0) < tol;
  testResult &= abs(a[2]-2.0) < tol;
  testResult &= abs(a[3]-5.0) < tol;
  testResult &= abs(a[4]+3.0) < tol;

  // 2nd test, has numeratorOrder < denominatorOrder-1 (needs zero padding of RHS), has also a
  // double root:
  // f(x) = 2/(x+2) - 2/(x+1) + 1/(x+1)^2
  //      = -x / (x^3+4*x^2+5*x+2)
  numDeg   = 1;
  denDeg   = 3;
  numPoles = 2;
  n[0] =  0; n[1] = -1;
  d[0] =  2; d[1] =  5; d[2] = 4; d[3] = 1;
  p[0] = -2; p[1] = -1;
  m[0] =  1; m[1] =  2;
  RF::partialFractionExpansion(n, numDeg, d, denDeg, p, m, numPoles, a);
  testResult &= abs(a[0]-2.0) < tol;
  testResult &= abs(a[1]+2.0) < tol;
  testResult &= abs(a[2]-1.0) < tol;

  // \todo: check a couple of other cases, one where the denominator is not monic is still missing
  // maybe move these tests into function below with more convenient code - merge all tests into 
  // one function

  return testResult;
}

bool testPartialFractionExpansion2(std::string& reportString)
{
  std::string testName = "PartialFractionExpansion2";
  bool testResult = true;

  typedef std::vector<std::complex<double>> Vec;
  typedef rsRationalFunction<double> RF;
  double tol = 1.e-14;    // tolerance for our equality checks
  Vec num, den, poles;    // numerator, denominator and distinct poles
  std::vector<int> muls;  // pole multiplicities
  Vec pfeCofs, polyCofs;  // coeffs of partial fraction expansion and polynomial part

  // f(x) = 2/(x+1) - 3/(x-2) + 5/(x-5) 
  //      = (4x^2 - 7x + 25) / ((x+1)*(x-2)*(x-5))
  //      = (4x^2 - 7x + 25) / (x^3 - 6x^2 + 3x + 10)
  num   = { 25,-7,  4    };  // p(x) =       4x^2 - 7x + 25
  den   = { 10, 3, -6, 1 };  // q(x) = x^3 - 6x^2 + 3x + 10
  poles = { -1, 2,  5    };  // q(x) = (x+1) * (x-2) * (x-5)
  muls  = {  1, 1,  1    };  // all multiplicities are 1
  pfeCofs = RF::partialFractions(num, den, poles, muls); // function for multiple poles
  testResult &= pfeCofs == Vec({ 2,-3, 5 });
  pfeCofs = RF::partialFractions(num, den, poles);       // function for distinct poles
  testResult &= pfeCofs == Vec({ 2,-3, 5 });

  // f(x) = 2/(x-1) + 3/(x-1)^2 + 1/(x-1)^3 + 1/(x-2)
  //      = (3x^3 - 8x^2 + 5x - 1) / (x^4 - 5x^3 + 9x^2 - 7x + 2)
  // 2 poles at 1 and 2, first is triple pole, second simple pole
  num     = { -1, 5,-8, 3    };
  den     = {  2,-7, 9,-5, 1 };
  poles   = {  1, 2 };
  muls    = {  3, 1 };
  pfeCofs = RF::partialFractions(num, den, poles, muls);
  testResult &= rsAreVectorsEqual(pfeCofs, Vec({2,3,1,1}), tol);

  // f(x) = 3/(x+5) - 4/(x+3) + 2/(x-1) + 5/(x-1)^2 - 3/(x-5)
  //      = (-2x^4 - 13x^3 + 25x^2 - 275x - 215) / (x^5 + x^4 - 30x^3 - 22x^2 + 125x-75)
  num   = {-215,-275,25,-13,-2};
  den   = {-75,125,-22,-30,1,1};
  poles = {-5,-3,+1,+5};
  muls  = { 1, 1, 2, 1}; 
  pfeCofs = RF::partialFractions(num, den, poles, muls);
  testResult &= rsAreVectorsEqual(pfeCofs, Vec({3,-4,2,5,-3}), tol);


  // todo: take function from partialFractionExpansion3 - this has a polynomial part, so we need a
  // convenience function that may also return a polynomial part...


  // try it also with functions where den is not monic


  return testResult;
}

bool testPolynomialBaseChange(std::string &reportString)
{
  std::string testName = "PolynomialBaseChange";
  bool testResult = true;

  static const int N = 7; // polynomial order
  double **Q; rsArrayTools::allocateSquareArray2D(Q, N+1);
  double **R; rsArrayTools::allocateSquareArray2D(R, N+1);
  double  *a = new double[N+1];
  double  *b = new double[N+1];

  // create two polynomial bases Q and R and the expansion coeffs in terms of the Q-polynomials
  // randomly:
  int i, j;
  rsRandomUniform(-1.0, 1.0, 0);
  for(i = 0; i <= N; i++)
  {
    a[i] = rsRandomUniform(-9.0, 9.0);
    for(j = 0; j <= N; j++)
    {
      Q[i][j] = rsRandomUniform(-9.0, 9.0);
      R[i][j] = rsRandomUniform(-9.0, 9.0);
    }
  }

  // get the expansion coeffs in terms of R-polynomials:
  rsPolynomial<double>::baseChange(Q, a, R, b, N);

  // select a value for the argument:
  double x = 2.0;

  // compute y in terms of Q-polynomials:
  double yQ = 0.0;
  for(i = 0; i <= N; i++)
    yQ += a[i] * rsPolynomial<double>::evaluate(x, Q[i], N);

  // compute y in terms of R-polynomials:
  double yR = 0.0;
  for(i = 0; i <= N; i++)
    yR += b[i] * rsPolynomial<double>::evaluate(x, R[i], N);

  testResult &= rsIsCloseTo(yQ, yR, 1.e-11);

  rsArrayTools::deAllocateSquareArray2D(Q, N+1);
  rsArrayTools::deAllocateSquareArray2D(R, N+1);
  delete[] a;
  delete[] b;

  return testResult;
}

void rsPowersToChebychev(double *a, double *b, int N)
{
  rsArrayTools::fillWithZeros(b, N+1);
  double tmp, tmp2;         // temporary values
  int k, i;                 // loop indices
  int s = 0;                // recursion stage
  b[0] = a[N];
  b[1] = a[N-1];
  for(k = N-2; k >= 0; k--)
  {
    s++;
    tmp  = b[0];
    b[0] = a[k] + 0.5*b[1];
    tmp2 = b[1];
    b[1] = tmp  + 0.5*b[2];
    tmp  = tmp2;
    for(i = 2; i <= s-1; i++)
    {
      tmp2 = b[i];
      b[i] = 0.5*(tmp + b[i+1]);
      tmp  = tmp2;
    }
    tmp2   = b[i];     // i == max(s, 2) here - this is what we use in the backwards algo
    b[i]   = 0.5*tmp;
    b[i+1] = 0.5*tmp2;
  }
}

void rsChebychevToPowers(double *b, double *a, int N)
{
  double tmp, tmp2;
  double *bb = new double[N+1]; // use a tmp-buffer, because it will be modified
  rsArrayTools::copy(b, bb, N+1);
  int k, i;

  // this is basically rsPowersToChebychev run backwards:
  int s = N-1;
  for(k = 0; k <= N-2; k++)
  {
    i    = rsMax(s, 2);
    tmp2 = 2*bb[i+1];
    tmp  = 2*bb[i];
    for(i = s-1; i >= 2; i--)
    {
      tmp2 = tmp;
      tmp  = 2*bb[i] - bb[i+1];
      bb[i] = tmp2;
    }
    tmp2 = tmp;
    tmp  = bb[1] - 0.5*bb[2];
    bb[1] = tmp2;
    a[k] = bb[0] - 0.5*bb[1];
    bb[0] = tmp;
    s--;
  }
  a[N-1] = bb[1];
  //a[N]   = bb[0];                  // wrong - why?
  a[N]   = bb[N] * rsPowInt(2, N-1); // works

  delete[] bb;
}

double rsEvaluateChebychevPolynomial(double x, int n)
{
  double t0 = 1.0;
  double t1 = x;
  double tn = 1.0;
  for(int i = 0; i < n; i++)
  {
    tn = 2*x*t1 - t0;
    t0 = t1;
    t1 = tn;
  }
  return t0;
}
// move to rsPolynomial - done - use it everywhere and delete this here

double rsEvaluateChebychevExpansion(double x, double *a, int N)
{
  double y = 0.0;
  for(int i = 0; i <= N; i++)
    y += a[i] * rsEvaluateChebychevPolynomial(x, i);
     // optimize this - resuse evaluation results from previous iterations, maybe lookup cleshaw
     // algorithm for a generalization
  return y;
}



bool testPowersChebychevExpansionConversion(std::string &reportString)
{
  std::string testName = "PowersChebychevExpansionConversion";
  bool testResult = true;

  // we my express any polynomial P(x) as linear combination of powers of x:
  // P(x) = a0 x^0 + a1 x^1 + a2 x^2 + ... + aN x^N
  // but we may also express it as linear cobination of Chebychev polynomials:
  // P(x) = b0 T0(x) + b1 T1(x) + b2 T2(x) + ... + bN TN(x)

  static const int N = 5;
  //double a[N+1] = {4, 1, -8, -8, 16, 16}; // polynomial coeffs, gives the Chebychev expansion
  //                                         // 6*T0 + 5*T1 + 4*T2 + 3*T3 + 2*T4 + 1*T5
  double a[N+1] = {9, 6, -10, -20, 24, 32}; // polynomial coeffs, gives the Chebychev expansion
                                            // 13*T0 + 11*T1 + 7*T2 + 5*T3 + 3*T4 + 2*T5
  double b[N+1]; // Chebychev expansion coeffs
  //double c[N+1]; // for reconstructed a-coeffs

  int k;
  int i;
  int s;      // recursion stage

  rsPowersToChebychev(a, b, N);  // nope - this function is still wrong

  // for test, we implement the algorithm in a way that stores all the intermediate arrays
  // in a 2D array:

  double B[N][N+1];
  memset(B, 0, N*(N+1)*sizeof(double));
  s = 0;      // recursion stage
  B[s][0] = a[N];
  B[s][1] = a[N-1]; // intitialization
  for(k = N-2; k >= 0; k--)
  {
    s++;
    B[s][0] = a[k]      + 0.5*B[s-1][1];
    B[s][1] = B[s-1][0] + 0.5*B[s-1][2];
    for(i = 2; i <= s+1; i++)
    {
      if( i < s )
        B[s][i] = 0.5*(B[s-1][i-1] + B[s-1][i+1]);
      else
        B[s][i] = 0.5*B[s-1][i-1];
    }


    /*
    for(i = 2; i < s; i++)
      B[s][i] = 0.5*(B[s-1][i-1] + B[s-1][i+1]);
    B[s][i] = 0.5*B[s-1][i-1];  // i == s here
    i++;
    B[s][i] = 0.5*B[s-1][i-1];  // i == s+1 here
    */
  }
  rsArrayTools::copy(B[N-1], b, N+1);
  // looks plausible and seems to work in this case


  //int dummy = 0;

  /*
  // reverse to algorithm: loops run backwards, order of instrcutions reversed, increments become
  // decremets, left-hand sides and right-hand sides of assignments exchange roles, when the
  // right-hand side contains a combination (i.e. sum) we have to look at the equation, find which
  // values are already known at this point and solve for the unknown which becomes the new
  // left-hand side:
  double C[N][N+1]; // we use C to reconstruct the B matrix
  memset(C, 0, N*(N+1)*sizeof(double));
  rsCopyBuffer(b, C[N-1], N+1); // the last stage is given, we must recostruct previous stages
  s = N-1;
  for(k = 0; k <= N-2; k++)
  {
    C[s-1][s] = 2*C[s][s+1];
    for(i = s; i >= 2; i--)
      C[s-1][i-1] = 2*C[s][i] - C[s-1][i+1];
    C[s-1][0] = C[s][1] - 0.5*C[s-1][2];
    c[k] = C[s][0] - 0.5*C[s-1][1];
    s--;
  }
  c[N-1] = C[s][1];
  c[N]   = C[s][0];
  */

  // \todo: get rid of the 2D-array to store all the intermediate stages - we may re-use a 1D array
  // at each stage when using 1 (or maybe 2) temporary variable
  // todo: move this commented stuff to the "Experiments" project - we may want to have it
  // available for later reference

  /*
  rsPowersToChebychev(a, b, N);  // nope - this function is still wrong
  rsChebychevToPowers(b, c, N);

  double yp = evaluatePolynomialAt(2.0, a, N);
  double yc = rsEvaluateChebychevExpansion(2.0, b, N);

  // todo: write some more exhaustive numerical unit tests - create polynomials with random
  // coeffients, write a function to evaluate chebychev-polynomials and compare outputs of
  // power- and chebychev-expansion

  rsFillWithRandomValues(a, N+1, -9.0, 9.0, 0);
  rsPowersToChebychev(a, b, N);
  yp = evaluatePolynomialAt(2.0, a, N);
  yc = rsEvaluateChebychevExpansion(2.0, b, N);
    // !!!! error !!!!
  for(int i = -9; i <= 9; i++)
  {
    double x  = i;
    double yp = evaluatePolynomialAt(x, a, N);
    double yc = rsEvaluateChebychevExpansion(x, b, N);
    testResult &= rsIsCloseTo(yp, yc, 1.e-12);
  }

  */

  return testResult;
}


bool testPolynomialRecursion(std::string &reportString)
{
  std::string testName = "PolynomialRecursion";
  bool testResult = true;

  static const int N = 5;   // == maxOrder-1
  double a[N][N];           // polynomial coefficients

  // we need pointers to doubles instead of a 2D array for the function calls:
  double *pa[N];
  int n;
  for(n = 0; n < N; n++)
    pa[n] = &a[n][0];

  // our 3-term recurrence is generally defined as:
  // w0 * P_n(x) = (w1 + w1x * x) * P_{n-1}(x) + w2 * P_{n-2}(x)
  // where we use: w0=1, w1=2, w1x=3, w2=5, P0(x)=1, P1(x)=x, so
  // P2(x) = (2+3x)x + 5 = 5 + 2x + 3x^2, etc.

  // set recursion coefficients:
  double w0  = 1.0;
  double w1  = 2.0;
  double w1x = 3.0;
  double w2  = 5.0;

  // set coefficients of of polynomials of orders 0 and 1 to intialize the recursion:
  a[0][0] = 1;
  a[1][0] = 0;
  a[1][1] = 1;

  // compute coefficient arrays for higher orders recursively:
  for(n = 2; n < N; n++)
    rsPolynomial<double>::threeTermRecursion(pa[n], w0, n, pa[n-1], w1, w1x, pa[n-2], w2);

  // P2(x) = 5 + 2x + 3x^2:
  testResult &= a[2][0]==5 && a[2][1]==2 && a[2][2]==3;

  // P3(x) = 10 + 24x + 12x^2 + 9x^3:
  testResult &= a[3][0]==10 && a[3][1]==24 && a[3][2]==12 && a[3][3]==9;

  // P4(x) = 45 + 88x + 111x^2 + 54x^3 + 27x^4:
  testResult &= a[4][0]==45 && a[4][1]==88 && a[4][2]==111 && a[4][3]==54 && a[4][4]==27;

  // in-place application - 1st input is reused as output:
  double t1[5], t2[5];
  rsArrayTools::copy(a[2], t2, 5);
  rsArrayTools::copy(a[3], t1, 5);
  rsPolynomial<double>::threeTermRecursion(t1, w0, 4, t1, w1, w1x, t2, w2);
  testResult &= rsArrayTools::equal(a[4], t1, 5);

  // in-place application - 2nd input is reused as output:
  rsArrayTools::copy(a[3], t1, 5);
  rsPolynomial<double>::threeTermRecursion(t2, w0, 4, t1, w1, w1x, t2, w2);
  testResult &= rsArrayTools::equal(a[4], t2, 5);

  return testResult;
}

bool testJacobiPolynomials(std::string &reportString)
{
  std::string testName = "JacobiPolynomials";
  bool testResult = true;

  static const int maxOrder = 4;   // maximum order of Jacobi polynomial
  static const int N = maxOrder+1; // length of longest coefficient array
  double a = 2;                    // alpha-coefficient
  double b = 3;                    // beta-coefficient
  double c[N][N];                  // polynomial coefficients

  // we need pointers to doubles instead of a 2D array for the function calls:
  double *pc[N];
  int n;
  for(n = 0; n < N; n++)
    pc[n] = &c[n][0];

  // gnereate the coefficient arrays:
  rsPolynomial<double>::jacobiPolynomials(&pc[0], a, b, maxOrder);

  // check coefficients:
  testResult &= c[0][0]==1;
  testResult &= c[1][0]==-0.5 && c[1][1]==3.5;
  testResult &= c[2][0]==-1 && c[2][1]==-2 && c[2][2]==9;
  testResult &= c[3][0]==0.625 && c[3][1]==-5.625 && c[3][2]==-5.625 && c[3][3]==20.625;
  testResult &= c[4][0]==0.9375 && c[4][1]==3.75 && c[4][2]==-20.625 && c[4][3]==-13.75
    && c[4][4]==44.6875;

  // special case: a=b=0 leads to Legendre polynomials - test in-place generation of them:
  double L1[8], L2[8]; // Legendre polynomial for N-1 and N-2, repectively
  L1[0] = 1;
  L2[0] = 0;
  L2[1] = 1;

  // L1 and L2 now contain Legendre polynomials of orders 0 and 1, we compute Legendre polynomial
  // of successively higher orders using recursion, using the two arrays alternately for the
  // in-place computed results:
  rsPolynomial<double>::legendreRecursion(L1, 2, L2, L1);
  rsPolynomial<double>::legendreRecursion(L2, 3, L1, L2);
  rsPolynomial<double>::legendreRecursion(L1, 4, L2, L1);
  rsPolynomial<double>::legendreRecursion(L2, 5, L1, L2);
  rsPolynomial<double>::legendreRecursion(L1, 6, L2, L1);
  rsPolynomial<double>::legendreRecursion(L2, 7, L1, L2);

  // check, if the 7th order Legendre coefficients are correct:
  testResult &= L2[0]==0 && L2[1]==-2.1875 && L2[2]==0  && L2[3]==19.6875 && L2[4]==0
     && L2[5]==-43.3125 && L2[6]==0  && L2[7]==26.8125;

  return testResult;
}

bool testSpecialPolynomials()
{
  bool r = true;

  using Poly = rsPolynomial<double>;

  double y1, y2;
  double tol = 1.e-10; 
  // we need a rather high tolerance and for high-degree polynomials we need more tolerance than 
  // for lower degree - i think, it is because for the high degree, the values at the evaluation 
  // point +-2 become large, so the absolute error gets large there too - or something

  // Compare recursive vs driect evaluation of Chebychev polynomials:
  int nMax = 8;  // maximum degree
  for(int n = 0; n < nMax; n++)
  {
    y1 = Poly::chebychevRecursive(-0.5, n);
    y2 = Poly::chebychevDirect(   -0.5, n);
    r &= rsIsCloseTo(y1, y2, tol);

    y1 = Poly::chebychevRecursive(+0.5, n);
    y2 = Poly::chebychevDirect(   +0.5, n);
    r &= rsIsCloseTo(y1, y2, tol);

    y1 = Poly::chebychevRecursive(-1.0, n);
    y2 = Poly::chebychevDirect(   -1.0, n);
    r &= rsIsCloseTo(y1, y2, tol);

    y1 = Poly::chebychevRecursive(+1.0, n);
    y2 = Poly::chebychevDirect(   +1.0, n);
    r &= rsIsCloseTo(y1, y2, tol);

    y1 = Poly::chebychevRecursive(-2.0, n);
    y2 = Poly::chebychevDirect(   -2.0, n);
    r &= rsIsCloseTo(y1, y2, tol);

    y1 = Poly::chebychevRecursive(+2.0, n);
    y2 = Poly::chebychevDirect(   +2.0, n);
    r &= rsIsCloseTo(y1, y2, tol);
    //rsAssert(r);
  }

  return r;
}

//-------------------------------------------------------------------------------------------------

// rename to testPolynomialClass
bool testPolynomialOperators(std::string &reportString)
{
  std::string testName = "PolynomialOperators";
  bool testResult = true;
  typedef rsPolynomial<double> PL;



  PL p({ 7,  5,  3,  2    });   // p(x) =  7 +  5*x +  3*x^2 +  2*x^3
  PL q({23, 19, 17, 13, 11});   // q(x) = 23 + 19*x + 17*x^2 + 13*x^3 + 11*x^4


  PL r, s;
  r = p + q; testResult &= r == PL({ 30, 24, 20, 15, 11 });
  r = q - p; testResult &= r == PL({ 16, 14, 14, 11, 11 });
  p.setRoots({  1., 2., 3., 4. });    testResult &= p == PL({ 24,  -50, 35, -10, 1 });
  q.setRoots({ -1.,-2.,-3.,-4. }, 2); testResult &= q == PL({ 48,  100, 70,  20, 2 });
  q.setRoots({ -1.,-2.,-3.,-4. });    testResult &= q == PL({ 24,   50, 35,  10, 1 });

  // multiply both polynomials together:
  r = p * q;
  s = q * p;
  //testResult &= r == s && r == PL({ 576, 0, -820, 0, 273, 0, -30, 0, 1, 0 });
  testResult &= r == s && r == PL({ 576, 0, -820, 0, 273, 0, -30, 0, 1});  
  // trailing zero is automatically truncated in multiplication

  // not yet finished - these assigmenz do not work as expected:
  p = PL({ 1 });   // coeff array is wrong! [0, 0] instead of [1]
  //p = PL({ 1.0 }); // this doesn't even compile
  p.setCoeffs({ 1   }); // works
  p.setCoeffs({ 1.0 }); // works


  // test truncation of the trailing zero
  p.setCoeffs({ 1., 0., 3., 0. });
  p.truncateTrailingZeros(); testResult &= p == PL({ 1., 0., 3. });

  p.setCoeffs({ 1., 0., 3., 0., 0. });
  p.truncateTrailingZeros(); testResult &= p == PL({ 1., 0., 3. });

  p.setCoeffs({ 1., 0., 0. });
  q.setCoeffs({ 1.         });
  p.truncateTrailingZeros();
  testResult &= p == q;
  //testResult &= p == PL({ 1 });  // PL({ 1.0 }): compiler error, PL({ 1 }): test fails - wtf?

  p.setCoeffs({ 0., 0., 0. });
  q.setCoeffs({ 0.         });
  p.truncateTrailingZeros();
  testResult &= p == q;
  //testResult &= p == PL({ 0., });

  // to tes:
  // -copy constructor and assignment operator
  // -evaluation operator (for numeric result and polynomial object result)


  return testResult;
}

bool testRationalFunction(std::string& reportString)
{
  std::string testName = "RationalFunction";
  bool testResult = true;

  typedef rsPolynomial<double> PL;
  typedef rsRationalFunction<double> RF;

  std::vector<double> p({ 6,7,1 }), q({-6,-5,1}), g;
  g = polyGCD(p, q, 0.0);  // result is 1 + x

  RF r({ 1,2,3 }, { 4, 5, 6, 7 });
  RF s({ 5,6 }, { 5, 7, 11 });

  // test arithmetic operators:
  RF t;
  t = r * s; testResult &= t == RF({ 5,16,27,18 }, { 20,53,109,132,115,77 });
  t = r / s; testResult &= t == RF({ 5,17,40,43,33 }, { 20,49,60,71,42 });
  t = r + s; testResult &= t == RF({ 25,66,100,114,75 }, { 20,53,109,132,115,77 }); 
  t = r - s; testResult &= t == RF({ -15,-32,-20,-28,-9 }, { 20,53,109,132,115,77 });

  // test nesting:
  double x= 4, y1, y2;  // input and outputs
  t = r(s);             // compose/nest functions r and s (r is outer, s is inner)
  y1 = r(s(x));         // evaluate s at x, pass result to r and evaluate r at s(x)
  y2 = t(x);            // evaluate the compsed function t = r(s)
  testResult &= y1 == y2;

  // giving this to sage:
  //
  // def canonicalRational(r):
  //     return expand(r.numerator())/expand(r.denominator())
  //
  // g(x) = (2+3*x) / (5-2*x)             # inner
  // f(x) = (2-3*x+4*x^2-5*x^3)/(2+3*x)   # outer
  // h = canonicalRational(f(g(x)))
  // h, n(h(5))
  //
  // produces:
  //
  //  -(259*x^3 - 90*x^2 + 377*x - 140)/(20*x^3 - 36*x^2 - 195*x + 400)),
  //  -31.0926829268293

  double tol = 1.e-10;
  s = RF({ 2,3 }, { 5,-2 });        // inner
  r = RF({ 2,-3,4,-5 }, { 2, 3 });  // outer
  t = r(s); // 700,-2165, 1204,-1475,518;  2000,-1775,210,172,-40
  y1 = t(5);  // -31.09268...
  testResult &= rsIsCloseTo(y1, -31.0926829268293, tol);
  // our numeric result from the evaluation -31.09.. is the same as the result from sage but
  // our coefficient arrays are different. numerator and denominator still have a common
  // factor -> divide both by their gcd:

  t.reduce(tol);
  // after reducing, the coefficients are consistent with sage's - but ours are still all 
  // multiplied by a factor of two compared to sage's result - we may have to divide all 
  // coefficients by *their* gcd - but for that, we would need floating point gcd with 
  // tolerance...maybe that doesn't make much sense as we do not really assume to coeffs to be 
  // integers anyway - and the gcd of a set of real (or complex) numbers is meaningless
  // https://stackoverflow.com/questions/16628088/euclidean-algorithm-gcd-with-multiple-numbers
  // https://www.geeksforgeeks.org/gcd-two-array-numbers/




  return testResult;
}

bool testQuadraticTo3Points()
{
  bool r = true;

  using Poly = rsPolynomial<double>;

  // test prototype fitQuadratic:
  double x1 =  1, y1 = 4;
  double x2 =  2, y2 = 9;
  double x3 = -1, y3 = 6;
  double a0, a1, a2;
  fitQuadratic(x1, y1, x2, y2, x3, y3, &a0, &a1, &a2);
  double z1 = a0 + a1*x1 + a2*x1*x1;
  double z2 = a0 + a1*x2 + a2*x2*x2;
  double z3 = a0 + a1*x3 + a2*x3*x3;
  r &= z1 == y1 && z2 == y2 && z3 == y3;  // zi should be equal to yi

  // test Poly::fitQuadraticDirect:
  double a[3];
  double x[3] ={ x1,x2,x3 }, y[3] ={y1,y2,y3};
  Poly::fitQuadraticDirect(a, x, y);
  r &= a[0] == a0 && a[1] == a1 && a[2] == a2;

  // test Poly::fitQuadraticLagrange:
  Poly::fitQuadraticLagrange(a, x, y);
  r &= a[0] == a0 && a[1] == a1 && a[2] == a2;

  return r;
}

bool testBivariatePolynomial()
{
  bool r = true;

  using Poly   = rsPolynomial<double>;
  using BiPoly = rsBivariatePolynomial<double>;


  BiPoly p(2, 3, {1,2,3,4, 5,6,7,8, 9,10,11,12});
  // p(x,y) =    1     + 2*y      + 3*y^2      + 4*y^3  
  //           + 5*x   + 6*x*y    + 7*x*y^2    + 8*x*y^3
  //           + 9*x^2 + 10*x^2*y + 11*x^2*y^2 + 12*x^2*y^3

  // sage:
  // var("x y")
  // p(x,y) = 1 + 2*y + 3*y^2 + 4*y^3 + 5*x + 6*x*y + 7*x*y^2 + 8*x*y^3 + 9*x^2 + 10*x^2*y + 11*x^2*y^2 + 12*x^2*y^3
  // p(2,3), p(2,y), p(x,3)
  //
  // (2594, 68*y^3 + 61*y^2 + 54*y + 47, 462*x^2 + 302*x + 142)

  double val, val2;
  Poly   uni, uni2;
  BiPoly bi;

  val = p.evaluate(2, 3); r &= val == 2594;
  uni = p.evaluateX(2);   r &= uni == Poly({47, 54, 61, 68});
  val = uni(3);           r &= val == 2594;
  uni = p.evaluateY(3);   r &= uni == Poly({142, 302, 462});
  val = uni(2);           r &= val == 2594;
  bi = p.derivativeX(); r &= bi == BiPoly(1, 3, {5,6,7,8, 18,20,22,24});
  //    5    + 6*y    + 7*y^2    + 8*y^3 
  //  + 18*x + 20*x*y + 22*x*y^2 + 24*x*y^3

  bi = p.derivativeY(); r &= bi == BiPoly(2, 2, {2,6,12, 6,14,24, 10,22,36});
  //    2      + 6*y      + 12*y^2
  //  + 6*x    + 14*x*y   + 24*x*y^2
  //  + 10*x^2 + 22*x^2*y + 36*x^2*y^2

  // Test indefinite integration with respect to x:
  bi = p.integralX(); BiPoly(3, 3, {0,0,0,0, 1,2,3,4, 2.5,3,3.5,4, 3,10./3,11./3,4});
  //   0       + 0*y        + 0*y^2        + 0*y^3
  // + 1*x     + 2*x*y      + 3*x*y^2      + 4*x*y^3
  // + 5/2*x^2 + 3*x^2*y    + 7/2*x^2*y^2  + 4*x^2*y^3
  // + 3*x^3   + 10/3*x^3*y + 11/3*x^3*y^2 + 4*x^3*y^3

  // Test indefinite integration with respect to y:
  bi = p.integralY(); r &= bi == BiPoly(2, 4, {0,1,1,1,1, 0,5,3,7*(1./3),2, 0,9,5,11*(1./3),3 });
  //   0     + y       + y^2 + y^3 + y^4
  // + 0*x   + 5*x*y   + 3*x*y^2   + 7/3*x*y^3    + 2*x*y^4
  // + 0*x^2 + 9*x^2*y + 5*x^2*y^2 + 11/3*x^2*y^3 + 3*x^2*y^4

  // Test definite integration with respect to x:
  Poly A({1,-2,3});
  Poly B({3,-5,1,-2});
  uni = p.integralX(-2., 3.); r &= uni == Poly({122.5, 425./3, 965./6, 180});
  uni = p.integralX( A,  3.); // verify!
  uni = p.integralX(-2., B ); // verify!
  uni = p.integralX( A,  B ); // verify!
  // todo: let either of the two or both integration limits be a polynomial (in y)

  // Test definite integration with respect to y:
  uni = p.integralY(-2., 3.); r &= uni == Poly({110, 755./3, 1180./3});
  // 1180/3*x^2 + 755/3*x + 110

  // Construct bivariate polynomial as product of two univariate polynomials:
  uni  = Poly({ 2,  3,  5,  7});
  uni2 = Poly({11, 13, 17, 19, 23});
  bi   = BiPoly::multiply(uni, uni2);
  val  = bi.evaluate(2, 3);
  val2 = uni(2) * uni2(3);
  r &= val == val2;

  // Construct bivariate polynomial from a univariate polynomial in the variable a*x + b*y for 
  // given a,b
  //
  // var("x y a b")
  // p(x) = 2 + 3*x + 5*x^2 + 7*x^3
  // p(11*x + 13*y).expand()
  //
  //   2        + 39*y        + 845*y^2     + 15379*y^3
  // + 33*x     + 1430*x*y    + 39039*x*y^2 + 0
  // + 605*x^2  + 33033*x^2*y + 0           + 0
  // + 9317*x^3 + 0           + 0           + 0 
  double a = 11, b = 13;
  double x = 2,  y = 3;
  bi   = BiPoly::composeWithLinear(uni, a, b);
  val  = uni(a*x + b*y);
  val2 = bi.evaluate(x, y);
  r &= val == val2;
  bi   = BiPoly::composeWithLinearOld(uni, a, b);
  val2 = bi.evaluate(x, y);
  r &= val == val2;

  // Multiply bivariate polynomial with univariate polynomial in y:
  //
  // var("x y")
  // p(x,y) = 1 + 2*y + 3*y^2 + 4*y^3 + 5*x + 6*x*y + 7*x*y^2 + 8*x*y^3 + 9*x^2 + 10*x^2*y + 11*x^2*y^2 + 12*x^2*y^3
  // q(y) = 2 + 3*y + 5*y^2 + 7*y^3
  // (p*q).expand()
  //
  //    2     + 7*y      + 17*y^2     + 34*y^3      + 41*y^4      + 41*y^5      + 28*y^6
  // + 10*x   + 27*x*y   + 57*x*y^2   + 102*x*y^3   + 101*x*y^4   + 89*x*y^5    + 56*x*y^6
  // + 18*x^2 + 47*x^2*y + 97*x^2*y^2 + 170*x^2*y^3 + 161*x^2*y^4 + 137*x^2*y^5 + 84*x^2*y^6
  bi = p.multiplyY(uni);
  r &= bi == BiPoly(2, 6, {2,7,17,34,41,41,28, 10,27,57,102,101,89,56, 18,47,97,170,161,137,84 });

  // Test evaluation with a polynomial as argument:
  // var("x y")
  // p(x,y) = 1 + 2*y + 3*y^2 + 4*y^3 + 5*x + 6*x*y + 7*x*y^2 + 8*x*y^3 + 9*x^2 + 10*x^2*y + 11*x^2*y^2 + 12*x^2*y^3
  // q(x) = 2 + 3*x + 5*x^2 + 7*x^3
  // p(x,q).expand()
  uni2 = p.evaluateY(uni);
  r &= uni2 == Poly({49,295,1112,3091,6862,12511,18686,23377,23795,18844,11564,4116});

  // var("x y")
  // p(x,y) = 1 + 2*y + 3*y^2 + 4*y^3 + 5*x + 6*x*y + 7*x*y^2 + 8*x*y^3 + 9*x^2 + 10*x^2*y + 11*x^2*y^2 + 12*x^2*y^3
  // q(y) = 2 + 3*y + 5*y^2 + 7*y^3
  // p(q,y).expand()
  uni2 = p.evaluateX(uni);
  r &= uni2 == Poly({47,177,485,1098,1747,2375,2630,2064,1379,588});


  // Test the construction of a potential and its gradient vector field from a given divergence:
  //BiPoly D = BiPoly(2, 3, {1,2,3,4, 5,6,7,8, 9,10,11,12});   // prescribed divergence
  //BiPoly D = BiPoly(2, 2, {1,2,3, 4,5,6, 7,8,9});   // prescribed divergence
  BiPoly D = BiPoly(3, 3, {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16});   // prescribed divergence
  BiPoly P, P2;
  divergenceToPotential4(D, P);                               // the scalar potential
  BiPoly P_x  = P.derivativeX();                              // x-component of vector field
  BiPoly P_y  = P.derivativeY();                              // y-component of vector field
  BiPoly P_xx = P_x.derivativeX();
  BiPoly P_yy = P_y.derivativeY();
  BiPoly D2   = P_xx + P_yy;
  // wrong: there are some extra nonzero coeffs in D2 - seems like there are 2 extra coeffs

  //P = BiPoly(3, 3, {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16});
  //P_x  = P.derivativeX();                              // x-component of vector field
  //P_y  = P.derivativeY();                              // y-component of vector field
  //P_xx = P_x.derivativeX();
  //P_yy = P_y.derivativeY();
  //D    = P_xx + P_yy;
  //divergenceToPotential(D, P2);

  P = BiPoly(4, 5, {1,2,3,4,5,6, 7,8,9,10,11,12, 13,14,15,16,17,18, 19,20,21,22,23,24, 
                    25,26,27,28,29,30});
  P_x  = P.derivativeX();                              // x-component of vector field
  P_y  = P.derivativeY();                              // y-component of vector field
  P_xx = P_x.derivativeX();
  P_yy = P_y.derivativeY();
  D    = P_xx + P_yy;
  potentialToDivergence(P, D2); 
  r &= D2 == D;
  divergenceToPotential1(D, P2);
  potentialToDivergence(P2, D2); 
  r &= D2.isCloseTo(D);
  divergenceToPotential2(D, P2);
  potentialToDivergence(P2, D2); 
  r &= D2.isCloseTo(D);

  //r &= D2 == D;

  //BiPoly D2  = BiPoly::divergence2D(P_x, P_y);               // actual divergence of vector field
  // this does not yet work!

  return r;
}

bool testBivariatePolynomial2()
{
  // We consider the complex polynomial f(z) = z^n - 1 which has as its roots the n n-th roots of 
  // unity. We want to find its Polya vector field and a potential function for that field

  using Real    = double;  // try to use float (gives currently compiler errors)
  using Complex = std::complex<Real>;
  using PolyR   = rsPolynomial<Real>;
  using PolyC   = rsPolynomial<Complex>;
  using BiPolyR = rsBivariatePolynomial<Real>;
  using BiPolyC = rsBivariatePolynomial<Complex>;

  using C = Complex;

  // Create the complex polynomial:
  //int n = 4; PolyC p(n); p[0] = -1; p[n] =  1; // p(z) = z^n - 1
  //PolyC p({1,3,5,-2,-1,2});
  //PolyC p({2,-1,3,-3,5,2,4});
  PolyC p({C(2,1),C(-1,2),C(3,-2),C(-3,1),C(5,-2),C(2,-1),C(4,-2)});  // with complex coeffs

  // Find the Polya potential P(x,y):
  BiPolyR P = BiPolyR::getPolyaPotential(p);

  // Check, if the found potential has indeed the desired gradient by evaluating the partial 
  // derivatives at some point. (todo: implement a way to directly compare the coeff matrices):
  BiPolyR P_x = P.derivativeX();
  BiPolyR P_y = P.derivativeY();

  // Tests:
  bool ok = true;

  // Evaluate p, P_x, P_y at some point x,y and compare results:
  Real x = 2.0, y = 3.0;
  Real val1, val2;
  Complex z(x, y);
  val1 =  p(z).real(); val2 = P_x(x, y); ok &= val1 == val2;
  val1 = -p(z).imag(); val2 = P_y(x, y); ok &= val1 == val2;

  // Test the harmonic conjugate relationships between P_x, P_y
  ok &= P_x.isHarmonic();
  ok &= isHarmonic2(P_x); 
  ok &= P_y.isHarmonic();
  ok &= isHarmonic2(P_y); 
  ok &= BiPolyR::areHarmonicConjugates(P_x, -P_y); // -P_y (not P_y!) is harmonic conjugate of P_x
  ok &= BiPolyR::areHarmonicConjugates(P_y,  P_x); // the relation is antisymmetric (?)

  // Try to reconstruct the original complex univariate polynomial from P_x, P_y
  PolyC p2 = getComplexFromHarmonicUV(P_x, -P_y);
  ok &= p2 == p;

  // experimental - find an "associated polynomial" by swapping the roles of P_x, P_y in the 
  // reconstruction:
  PolyC p3 = getComplexFromHarmonicUV(P_y, P_x);
  // looks like p3 has the real and imaginary parts swapped with respect to p - but also with a 
  // sign flip: p3[i] = -Im{ p[i] } + i*Re{ p[i] } ... i think, this is a rotation by 90°, i.e. a
  // multiplication by i? let's try it:
  PolyC p4 = p; p4.scale(C(0,1)); ok &= p3 == p4;   // ...yep, indeed - nice!

  //p2 = getComplexFromHarmonicU(P_x); // should reconstruct p
  //p2 = getComplexFromHarmonicV(P_y); // seems to work except for the constant term


  // evaluate the potential at some points:
  val1 = P(+1, 0);  // -2/3 for p(z) = z^2 - 1, -4/5 for p(z) = z^4
  val2 = P(-1, 0);  // +2/3 for p(z) = z^2 - 1, +4/5 for p(z) = z^4

  // Observations:
  // -at the zeros of the original function, we seem to see saddles in the potential
  // -so far, i've not yet seen a potential that features minima or maxima - only saddles seem to
  //  occur - investigate if this is always the case and try to find a theoretical explanation
  //  ...maybe it's only because so far i have not used polynomials with complex coeffs?
  //  -could this be related to features of minimal surfaces? they also tend to from saddles
  // -saddles are easy to identify in topographic or equipotential plots - they look like crosses
  // -could the value at the zeros be given by +-(1 - 1/(n+1)) for p(z) = z^n - 1, or maybe
  //  exp(i*k) * (1 - 1/(n+1)) for k = 0...n or something?

  // todo: 
  // -plot some potentials to get a feeling for how they look like and how we can see features of
  //  complex functions in them
  // -plot them as topographic map - draw equipotentials and maybe field lines, too - obtain them 
  //  by considering the differential equation (d/dt) f(x,y) = -grad(P(x,y))...or something
  //  ...or maybe we can find the harmonic conjugate of the potential - i think, its equipotentials 
  //  are indeed the field lines...i think, from these two, the complex potential can be 
  //  constructed?
  // -maybe use a coloring according to the absolute value of the gradient? that will show extrema
  //  of the magnitude of the complex function
  // -make a unit test and test it with more complicated complex polynomials p(z)





  // test the harmonic conjugate creation function:
  BiPolyR u, v, u2, v2;

  u = BiPolyR(2,2,{0,0,-1, 0,0,0, 1,0,0});  // x^2 - y^2
  v = BiPolyR(2,2,{0,0,0,  0,2,0, 0,0,0});  // 2*x*y
  v2 = u.getHarmonicConjugate();
  u2 = v.getHarmonicConjugate();
  ok &= BiPolyR::areHarmonicConjugates(u,  v );    // make sure that the target functions are correct
  ok &= BiPolyR::areHarmonicConjugates(u,  v2);    // works

  //ok &= areHarmonicConjugates(u2, v );    // fails!
  //ok &= areHarmonicConjugates(u2, v2);
  // hmm - is it actually the case that if v is the harmonic conjugate of u then u is also the 
  // harmonic conjugate of v? or is -u the harmonic conjugate of v? this is what seems to come out
  // of the algorithm

  // u(x,y) = 2*y^3 ? 6*x^2*y + 4*x^2 + 7*x*y + 4*y^2 + 3*x + 4*y + 4
  // v(x,y) = 2*x^3 + (7/2)*x^2 + 6*x*y^2 + 8*x*y + 4*x + (7/2)*y^2 + 3*y
  u  = BiPolyR(2,3,{-4,4,-4,2, 3,-7,0,0, 4,-6,0,0}); 
  v  = BiPolyR(3,2,{0,3,-3.5, -4,8,-6, 3.5,0,0, 2,0,0}); 
  v2 = u.getHarmonicConjugate();
  ok &= v2.isCloseTo(v);

  // Example from:
  // https://math.stackexchange.com/questions/930000/calculating-a-harmonic-conjugate

  //ok &= u2.isCloseTo(u);


  //BiPolyR Q = getHarmonicConjugate(P);
  //ok &= areHarmonicConjugates(P, Q);

  Real r = 1.5;
  //plotBivariatePolynomial(P, -r, +r, 31, -r, +r, 31);

  u = BiPolyR(2,2,{0,0,1, 0,0,0, -1,0,0}); // -x^2 + y^2
  bool harm = u.isHarmonic();


  // Arithmetic operators:
  u = BiPolyR(2,2,{1,2,3, 4,5,6, 7,8,9});
  v = BiPolyR(2,2,{9,8,7, 6,5,4, 3,2,1});
  BiPolyR uv;
  uv = u+v; val1 = uv(x, y); val2 = u(x, y) + v(x, y); ok &= val1 == val2;
  uv = u-v; val1 = uv(x, y); val2 = u(x, y) - v(x, y); ok &= val1 == val2;
  uv = u*v; val1 = uv(x, y); val2 = u(x, y) * v(x, y); ok &= val1 == val2;
  // maybe use polynomials with more interesting degrees (all values different and maybe higher)

  // Composition of outer bivariate with inner univariate polynomial:
  PolyR xt({1,2,3,4});
  PolyR yt({2,3});
  PolyR ut = BiPolyR::compose(u, xt, yt);
  Real t = 5;
  val1 = u(xt(t), yt(t));
  val2 = ut(t);
  ok  &= val1 == val2;


  // Test evaluation of various types of integrals:

  // Double intgral over a rectangular region:
  // p(x,y) =   1*x^0*y^0 + 2*x^0*y^1 + 3*x^0*y^2
  //          + 4*x^1*y^0 + 5*x^1*y^1 + 6*x^1*y^2
  //          + 7*x^2*y^0 + 8*x^2*y^1 + 9*x^2*y^2
  // x0 = 1, x1 = 2, y0 = 3, y1 = 4
  //
  // var("x y")
  // p(x,y) =   1*x^0*y^0 + 2*x^0*y^1 + 3*x^0*y^2 + 4*x^1*y^0 + 5*x^1*y^1 + 6*x^1*y^2 + 7*x^2*y^0 + 8*x^2*y^1 + 9*x^2*y^2       
  // p_ix = integrate(p, x, 1, 2)
  // p_ix_iy = integrate(p_ix, y, 3, 4)
  // p_ix_iy
  //
  // 6347/12 = 528.916666666667
  val1 = u.doubleIntegralXY(1., 2., 3., 4.); ok &= val1 == 6347./12.;
  val2 = u.doubleIntegralYX(1., 2., 3., 4.); ok &= val2 == 6347./12.;

  // Path integral over a vector field:
  // var("t")
  // x(t) = 1*t^0 + 2*t^1 + 3*t^2 + 4*t^3
  // y(t) = 2*t^0 + 3*t^1
  // u = 1*x^0*y^0 + 2*x^0*y^1 + 3*x^0*y^2 + 4*x^1*y^0 + 5*x^1*y^1 + 6*x^1*y^2 + 7*x^2*y^0 + 8*x^2*y^1 + 9*x^2*y^2
  // v = 9*x^0*y^0 + 8*x^0*y^1 + 7*x^0*y^2 + 6*x^1*y^0 + 5*x^1*y^1 + 4*x^1*y^2 + 3*x^2*y^0 + 2*x^2*y^1 + 1*x^2*y^2
  // s = u*diff(x,t) + v*diff(y,t)
  // r = integrate(s, t, -1, 2)
  // r, N(r)
  //
  // 3392240031/154, 2.20275326688312e7
  Real tol = 1.e-8;
  val1 = BiPolyR::pathIntegral(u, v, xt, yt, -1., +2.);
  val2 = 3392240031./154;;
  ok  &= rsIsCloseTo(val1, val2, tol);
  // todo: maybe check the relative instead of the absolute error (absolute error is largeish)


  // Path integral over a scalar field:
  // var("t")
  // x(t) = 1*t^0 + 2*t^1
  // y(t) = 2*t^0 + 3*t^1
  // p = 1*x^0*y^0 + 2*x^0*y^1 + 3*x^0*y^2 + 4*x^1*y^0 + 5*x^1*y^1 + 6*x^1*y^2 + 7*x^2*y^0 + 8*x^2*y^1 + 9*x^2*y^2
  // dx = diff(x,t)
  // dy = diff(y,t)
  // ds = sqrt(dx*dx + dy*dy) 
  // s = p * ds 
  // r = integrate(s, t, -1, 2)
  // r, N(r)
  //
  // 102399/10*sqrt(13), 36920.4845056237
  Real err;
  tol = 1.e-15;
  xt = PolyR({1,2});
  yt = PolyR({2,3});
  val1 = BiPolyR::pathIntegral(u, xt, yt, -1., +2.);
  val2 = (102399./10.)*sqrt(13);
  err  = (val2-val1)/val2;                // relative error - almost 1.e-4 - why so large?
  ok  &= rsAbs(err) <= tol;
  //ok  &= rsIsCloseTo(val1, val2, tol);

  // https://www.youtube.com/watch?v=7mrsZzXmibg 

  // Test Green's theorem in 2D: compare the double integral of the curl over a rectangle to the 
  // path integral of the vector field itself along the boundary of the rectangle:
  tol = 1.e-13;
  Real x0 =  1, x1 = 3;
  Real y0 = -1, y1 = 2;
  BiPolyR curl = BiPolyR::curl2D(u, v);
  val1 = curl.doubleIntegralXY(x0, x1, y0, y1);
  val2 = BiPolyR::loopIntegral(u, v, x0, x1, y0, y1);
  ok  &= rsIsCloseTo(val1, val2, tol);
  // val1 seems to be more numerically precise

  // Test Gauss' theorem in 2D:
  BiPolyR divergence = BiPolyR::divergence2D(u, v);
  val1 = divergence.doubleIntegralXY(x0, x1, y0, y1);
  val2 = BiPolyR::outfluxIntegral(u, v, x0, x1, y0, y1);
  ok  &= rsIsCloseTo(val1, val2, tol);

  // Test the flux integral with more general curves:
  // see https://www.khanacademy.org/math/multivariable-calculus/integrating-multivariable-functions/line-integrals-in-vector-fields-articles/a/flux-in-two-dimensions



  // Curvature of surfaces:
  BiPolyR w = BiPolyR(2,2,{6,5,4, 9,8,7, 3,2,1});
  BiPolyR E, F, G;
  firstFundamentalForm(u, v, w, E, F, G);

  // todo:
  // -flux and circulation integrals (integrals over divergence and curl)
  // -double integral over region bounded by polynomial curves

  rsAssert(ok);
  return ok;
}

bool testTrivariatePolynomial()
{
  using Poly    = rsPolynomial<double>;
  using BiPoly  = rsBivariatePolynomial<double>;
  using TriPoly = rsTrivariatePolynomial<double>;

  TriPoly p(3, 4, 5);
  p.fillRandomly(-9.0, +9.0, 0, true); // maybe round to integers

  bool ok = true;

  double x = 5, y = 3, z = 2;
  double val1, val2, val3;
  val1 = p.evaluate(x, y, z);
  BiPoly p_yz = p.evaluateX(x);
  val2 = p_yz.evaluate(y, z);
  ok &= val1 == val2;

  // test compose:
  double u = 2, v = 3;
  BiPoly px(2,2,{1,2,3, 4,5,6, 7,8,9});           // px(u,v)
  BiPoly py(2,2,{4,5,6, 1,2,3, 7,8,9});           // py(u,v)
  BiPoly pz(2,2,{7,8,9, 4,5,6, 1,2,3});           // pz(u,v)
  p = TriPoly(2,2,2);
  p.fillRandomly(-3.0, +3.0, 0, true);
  BiPoly p_uv = TriPoly::compose(p, px, py, pz);  // p(u,v)
  x = px(u,v);
  y = py(u,v);
  z = pz(u,v);
  val1 = p.evaluate(x, y, z);
  val2 = p_uv.evaluate(u, v);
  ok &= val1 == val2;
  // ...values are getting big here!

  // Test integration and differentiation:
  TriPoly tmp;
  tmp = p.integralX().derivativeX(); ok &= p == tmp;
  tmp = p.integralY().derivativeY(); ok &= p == tmp;
  tmp = p.integralZ().derivativeZ(); ok &= p == tmp;


  // test triple integral:
  double x0 = 1, x1 = 3, y0 = 2, y1 = 4, z0 = 3, z1 = 5;
  p = TriPoly(1,1,1);
  p.coeff(1,1,1) = 3.0;  // p(x,y,z) = 3*x*y*z
  val1 = p.tripleIntegralXYZ(x0, x1, y0, y1, z0, z1);
  ok &= val1 == 576;
  // var("x y z")
  // f = 3*x*y*z
  // Fx   = integral(f,   (x, 1, 3))
  // Fxy  = integral(Fx,  (y, 2, 4))
  // Fxyz = integral(Fxy, (z, 3, 5))
  // f, Fx, Fxy, Fxyz
  //
  // 3*x*y*z, 12*y*z, 72*z, 576

  // todo: implement and test different orders of integration: XZY, YXZ, YZX, ZXY, ZYX




  tmp = p.derivativeX(); val1 = tmp.evaluate(1,2,3); ok &= val1 == 18;
  tmp = p.derivativeY(); val1 = tmp.evaluate(1,2,3); ok &= val1 == 9;
  tmp = p.derivativeZ(); val1 = tmp.evaluate(1,2,3); ok &= val1 == 6;

  // todo: test fluxIntegral
  double u0 = -1, u1 = 1, v0 = -1, v1 = 1;
  TriPoly fx(2,2,2), fy(2,2,2), fz(2,2,2);  // the vector field
  fx.fillRandomly(-3.0, +3.0, 0, true);
  fy.fillRandomly(-3.0, +3.0, 1, true);     // we need to use different seeds
  fz.fillRandomly(-3.0, +3.0, 2, true);
  //val1 = TriPoly::fluxIntegral(fx, fy, fz, px, py, pz, u0, u1, v0, v1);
  // very big value - maybe use smaller, non-integer coeffs to get more reasonable values

  // Test Gauss' theorem:
  double tol = 1.e-10;
  x0 = 2, x1 = 4, y0 = 3, y1 = 5, z0 = 4, z1 = 6;
  val1 = TriPoly::outfluxIntegral(fx, fy, fz, x0, x1, y0, y1, z0, z1);
  TriPoly divergence = TriPoly::divergence(fx, fy, fz);
  val2 = divergence.tripleIntegralXYZ(x0, x1, y0, y1, z0, z1);
  ok &= rsIsCloseTo(val1, val2, tol);
  // val2 seems to be more precise, which makes sense because the triple integral is a simpler
  // computation than the flux integral

  // Test Stokes' theorem:
  using Vec3 = rsVector3D<double>;
  //Vec3 P0(0,0,0), P1(1,0,0), P2(0,1,0);
  Vec3 P0(1,2,3), P1(2,-3,-1), P2(-1,2,1);
  std::vector<Vec3> path({P0, P1, P2, P0});
  val1 = TriPoly::pathIntegral(fx, fy, fz, path);       // circulation around triangular loop
  TriPoly cx, cy, cz;                                   // components of curl
  TriPoly::curl(fx, fy, fz, cx, cy, cz);
  val2 = TriPoly::fluxIntegral(cx, cy, cz, P0, P1, P2); // flux of curl through triangular patch
  ok &= rsIsCloseTo(val1, val2, tol);

  // ...more interesting is that the theorem holds for any surface that is bounded by the given 
  // loop - it may totally bulge or balloon away from the boundary loop (here we used just a flat 
  // (triangular) patch). We can test it by choosing a 4th point P4 and compute the 3 flux 
  // integrals through: (P0, P4, P1), (P1, P4, P2), (P2, P4, P0):

  Vec3 P4(3,1,-3);
  //Vec3 P4(0,0,0);
  //Vec3 P4 = (P0 + P1 + P2) / 3.0;
  val3  = TriPoly::fluxIntegral(cx, cy, cz, P0, P4, P1);
  val3 += TriPoly::fluxIntegral(cx, cy, cz, P1, P4, P2);
  val3 += TriPoly::fluxIntegral(cx, cy, cz, P2, P4, P0);
  // the value has a minus sign compared to val1, val2 - no matter how we choose P4 (well, tested 
  // only with a few)

  // reversing the orientations fixes this:
  val3  = TriPoly::fluxIntegral(cx, cy, cz, P1, P4, P0);
  val3 += TriPoly::fluxIntegral(cx, cy, cz, P2, P4, P1);
  val3 += TriPoly::fluxIntegral(cx, cy, cz, P0, P4, P2);
  // -> figure out why we need this

  val3  = TriPoly::fluxIntegral(cx, cy, cz, P1, P0, P4);
  val3 += TriPoly::fluxIntegral(cx, cy, cz, P2, P1, P4);
  val3 += TriPoly::fluxIntegral(cx, cy, cz, P0, P2, P4);
  // has also negative sign


  // Test Green's integral formulas (https://en.wikipedia.org/wiki/Green%27s_identities or
  // Bärwolff, pg 625):
  TriPoly f(4,4,4), g(4,4,4), gx, gy, gz;
  f.fillRandomly(-3.0, +3.0, 3, true);
  g.fillRandomly(-3.0, +3.0, 4, true);
  TriPoly::gradient(f, fx, fy, fz);
  TriPoly::gradient(g, gx, gy, gz);

  // First Green formula (Eq 8.31):
  val1 = TriPoly::outfluxIntegral(g*fx, g*fy, g*fz, x0, x1, y0, y1, z0, z1);  // lhs
  tmp = g*f.laplacian() + gx*fx + gy*fy + gz*fz;
  val2 = tmp.tripleIntegralXYZ(x0, x1, y0, y1, z0, z1);  // rhs
  double err = val2 - val1;
  ok &= err == 0.0;
  // val1 and val2 seem to be equal indeed, but the values are ridiculously large

  // Second Green formula (Eq. 8.32):
  val1 = TriPoly::outfluxIntegral(g*fx-f*gx, g*fy-f*gy, g*fz-f*gz, x0, x1, y0, y1, z0, z1);  // lhs
  tmp = g*f.laplacian() - f*g.laplacian();
  val2 = tmp.tripleIntegralXYZ(x0, x1, y0, y1, z0, z1);  // rhs
  err = val2 - val1;
  ok &= err == 0.0;

  // todo: 
  // -test 8.40 (page 629) for the 2D version of the 1st formula (but in the BiPoly unit test)
  //  is there also a 2D version of the second formula?
  // -implement and test the 1D version (also on page 625)..this has something to do with 
  //  integration by parts - maybe implement that in rsPolynomial in the 1D case, maybe also the
  //  substitution rule, if that makes sense

  // Test scalar potential:
  //tmp = getPotential(fx, fy, fz);
  tmp = TriPoly::scalarPotential(fx, fy, fz);   // maybe rename to scalarPotential
  gx = tmp.derivativeX(); ok &= fx == gx;
  gy = tmp.derivativeY(); ok &= fy == gy;
  gz = tmp.derivativeZ(); ok &= fz == gz;
  ok &= tmp == f; // works only because f's (0,0,0) coeff happens to be zero - we should really 
                  // compare only the coeffs excluding the constant term in general

  // todo: implement functions hasScalarPotential, hasVectorPotential
  // maybe implement the path-integral method to compute a (scalar) potential, too (Bärwollf pg
  // 563)

  // see
  // https://tutorial.math.lamar.edu/classes/calcIII/conservativevectorfield.aspx


  // for vector potentials, see:
  // http://galileo.math.siu.edu/Courses/251/S12/vpot.pdf
  // -integrate Fx, Fy with respect to z, call them Gx, Gy, negate Gy
  // -...

  using TP = TriPoly;

  fx.fillRandomly(-3.0, +3.0, 3, true);
  fy.fillRandomly(-3.0, +3.0, 3, true);
  fz.fillRandomly(-3.0, +3.0, 3, true);
  TP::curl(fx, fy, fz, cx, cy, cz);             // obtain curl field c using f as vector potential
  TP::vectorPotential(cx, cy, cz, gx, gy, gz);  // construct vector potential g with same curl c
  ok &= TP::isVectorPotential(gx, gy, gz, cx, cy, cz, tol);
  divergence = TP::divergence(gx, gy, gz);

  // test 2nd implementation with given divergence:
  TP d(0,0,0);  // the zero polynomial
  vectorPotential2(cx, cy, cz, d, gx, gy, gz);
  //ok &= TP::isVectorPotential(gx, gy, gz, cx, cy, cz, tol);
  divergence = TP::divergence(gx, gy, gz);
  // divergence is indeed zero, but g seems not to be a vector potential anymore - some coeffs in 
  // the computed curl match but others are totally wrong (it's not a precision issue) - i guess 
  // our integration constants a(x,y), b(x,y) do in fact depend on z? maybe we an't simply set
  // b_y = 0 as before anymore


  fx = TP(1,1,1); fx.coeff(0,1,1) = 1;          // fx(x,y,z) = y*z
  fy = TP(1,1,1); fy.coeff(1,0,1) = 1;          // fx(x,y,z) = x*z
  fz = TP(1,1,1); fz.coeff(1,1,0) = 1;          // fz(x,y,z) = x*y
  TP::vectorPotential(fx, fy, fz, gx, gy, gz);
  ok &= TP::isVectorPotential(gx, gy, gz, fx, fy, fz, tol);
  divergence = TP::divergence(gx, gy, gz);

  // -maybe test, if we can add any conservative vector field to the vector potential with 
  //  destroying its vector potential property
  // -figure out meaningful ways to add such conservative vector fields to simplify certain 
  //  calculations - this is called "fixing the gauge"





  // constant flow in z-direction with unit speed, through a unit-square plane in (x,y) at z = 0:
  double fluxX = 2, fluxY = 3, fluxZ = 4;

  double dx = x1-x0, dy = y1-y0, dz = z1-z0;
  double tgt;  // target
  fx = TriPoly(0, 0, 0, { fluxX });
  fy = TriPoly(0, 0, 0, { fluxY });
  fz = TriPoly(0, 0, 0, { fluxZ });

  // We imagine, x goes rightward, y goes forward, z goes upward

  // top and bottom plane (z constant):
  tgt  = fluxZ*dx*dy;
  px   = BiPoly(1, 1, { x0, 0, dx, 0 });
  py   = BiPoly(1, 1, { y0, dy, 0, 0 });
  pz   = BiPoly(0, 0, { z0           });
  val1 = TriPoly::fluxIntegral(fx, fy, fz, px, py, pz, 0., 1., 0., 1.);
  ok &= val1 == tgt;
  pz = BiPoly(0, 0, { z1             });
  val1 = TriPoly::fluxIntegral(fx, fy, fz, px, py, pz, 0., 1., 0., 1.);
  ok &= val1 == tgt;

  // front and back plane (y constant):
  tgt  = fluxY*dx*dz;
  px   = BiPoly(1, 1, { x0, 0,  dx, 0 });
  py   = BiPoly(0, 0, { y0            });
  pz   = BiPoly(1, 1, { z1, -dz, 0, 0 }); // z takes the role of y, direction reversed, see below
  val1 = TriPoly::fluxIntegral(fx, fy, fz, px, py, pz, 0., 1., 0., 1.);
  ok  &= val1 == tgt; 
  py = BiPoly(0, 0, { y1             });
  val1 = TriPoly::fluxIntegral(fx, fy, fz, px, py, pz, 0., 1., 0., 1.);
  ok &= val1 == tgt; 
  // We can use:
  //   px = BiPoly(1, 1, { x1, -dx, 0, 0 });
  //   py = BiPoly(0, 0, { y0            });
  //   pz = BiPoly(1, 1, { z0, 0,  dz, 0 });
  // or 
  //   px = BiPoly(1, 1, { x0, dx,  0, 0 });
  //   py = BiPoly(0, 0, { y0            });
  //   pz = BiPoly(1, 1, { z1, 0, -dz, 0 });
  // as parametrization, but not:
  //   px = BiPoly(1, 1, { x0, 0, dx, 0 });
  //   py = BiPoly(0, 0, { y0           });
  //   pz = BiPoly(1, 1, { z0, dz, 0, 0 }); 
  // because then the result will have the wrong sign. This is because when i,j,k are the unit 
  // vectors in the x,y,z directions, then we have: i x j = k, i x k = -j, j x k = i (x denoting 
  // the cross product).

  // left and right plane (x constant):
  tgt  = fluxX*dy*dz;
  px   = BiPoly(0, 0, { x0           });  // x is constant
  py   = BiPoly(1, 1, { y0, 0, dy, 0 });  // y takes the role of x
  pz   = BiPoly(1, 1, { z0, dz, 0, 0 });  // z takes the role of y
  val1 = TriPoly::fluxIntegral(fx, fy, fz, px, py, pz, 0., 1., 0., 1.);
  ok  &= val1 == tgt; 
  px   = BiPoly(0, 0, { x1           });
  val1 = TriPoly::fluxIntegral(fx, fy, fz, px, py, pz, 0., 1., 0., 1.);
  ok  &= val1 == tgt; 

  // compute outflux for simple linear fields like fx = a + b*x, fy = 0, fz = 0
  fx = TriPoly(1, 1, 1);
  fx.coeff(1, 0, 0) = fluxX;
  fy = fz = TriPoly(0,0,0);
  tgt  = fluxX*dy*dz;
  px   = BiPoly(0, 0, { x0           });
  py   = BiPoly(1, 1, { y0, 0, dy, 0 });
  pz   = BiPoly(1, 1, { z0, dz, 0, 0 });
  val1 = TriPoly::fluxIntegral(fx, fy, fz, px, py, pz, 0., 1., 0., 1.);
  val2 = tgt*x0;
  ok  &= val1 == val2; 
  px   = BiPoly(0, 0, { x1           });
  val1 = TriPoly::fluxIntegral(fx, fy, fz, px, py, pz, 0., 1., 0., 1.);
  val2 = tgt*x1;
  ok  &= val1 == val2; 
  val1 = TriPoly::outfluxIntegral(fx, fy, fz, x0, x1, y0, y1, z0, z1);
  val2 = tgt*(x1-x0);
  ok  &= val1 == val2;
  divergence = TriPoly::divergence(fx, fy, fz);
  val3 = divergence.tripleIntegralXYZ(x0, x1, y0, y1, z0, z1);
  ok  &= val1 == val3;

  fz = TriPoly(1, 1, 1);
  fz.coeff(0, 0, 1) = fluxZ;
  fx = fy = TriPoly(0,0,0);
  val1 = TriPoly::outfluxIntegral(fx, fy, fz, x0, x1, y0, y1, z0, z1);
  tgt  = fluxZ*dx*dy;
  val2 = tgt*(z1-z0);
  ok  &= val1 == val2; 
  divergence = TriPoly::divergence(fx, fy, fz);
  val3 = divergence.tripleIntegralXYZ(x0, x1, y0, y1, z0, z1);
  ok  &= val1 == val3;

  fy = TriPoly(1, 1, 1);
  fy.coeff(0, 1, 0) = fluxY;
  fx = fz = TriPoly(0,0,0);
  val1 = TriPoly::outfluxIntegral(fx, fy, fz, x0, x1, y0, y1, z0, z1);
  tgt  = fluxY*dx*dz;
  val2 = tgt*(y1-y0);
  ok  &= val1 == val2; 
  divergence = TriPoly::divergence(fx, fy, fz);
  val3 = divergence.tripleIntegralXYZ(x0, x1, y0, y1, z0, z1);
  ok  &= val1 == val3;

  // Compute path integral of vector field:
  Poly xt({1,2,3,4});
  Poly yt({2,3,4,5});
  Poly zt({3,4,5,6});
  val1 = TriPoly::pathIntegral(fx, fy, fz, xt, yt, zt, 0.0, 1.0); // 288

  rsAssert(ok);
  return ok;
}



template<class T>
bool rsIsCloseTo(const RAPT::rsPolynomial<T>& p, const RAPT::rsPolynomial<T>& q, T tol)
{
  // Compare coeffs up to the highest one of lower degree polynomial:
  int degP = p.getDegree();
  int degQ = q.getDegree();
  for(int i = 0; i <= rsMin(degP, degQ); i++)
    if( rsAbs(p.getCoeff(i)-q.getCoeff(i)) > tol )
      return false;
  if(degP == degQ)
    return true;    // if degrees match, we are done here

  // If they have different degrees, check, if the higher order coeffs of the higher degree 
  // polynomial are close to zero:
  if(degP > degQ) { 
    for(int i = degQ+1; i <= degP; i++)
      if( rsAbs(p.getCoeff(i)) > tol )
        return false; }
  else {  // degQ > degP
    for(int i = degP+1; i <= degQ; i++)
      if( rsAbs(q.getCoeff(i)) > tol )
        return false; }
  return true;
}
// move to TestTools

bool testPiecewisePolynomial1()
{
  bool r = true;

  using Poly      = RAPT::rsPolynomial<double>;
  using PiecePoly = rsPiecewisePolynomial<double>;
  Poly p({ 2,-3,5,-7 }); double pL = -1, pU = 2; // p(x) = 2 - 3*x + 5*x^2 - 7*x^3  in -1..2
  Poly q({ 3,-4,6    }); double qL =  3, qU = 4; // q(x) = 3 - 4*x + 6*x^2          in  3..4

  // Test convolving two polynomial pieces. This gives 3 segments.
  // Sage:
  // p  = piecewise([((-1,2), 2 - 3*x + 5*x^2 - 7*x^3)])
  // q  = piecewise([(( 3,4), 3 - 4*x + 6*x^2)])
  // pq = p.convolution(q)
  // pq
  // #plot([p,q,pq],xmin=-1,xmax=6)
  //
  // gives:
  // pq = 
  // -7/10*x^6 + 12/5*x^5 - 101/12*x^4 + 326*x^3 - 2271*x^2 + 92231/15*x - 28594/5 on (2, 3]
  // -441*x^3 + 5012*x^2 - 288133/15*x + 498143/20 on (3, 5] 
  // 7/10*x^6 - 12/5*x^5 + 101/12*x^4 - 767*x^3 + 14449/2*x^2 - 124586/5*x + 150942/5 on (5, 6]

  // Define target pieces:
  double tol = 1.e-11; // we need a rather high tolerance because the absolute values are quite 
                       // large -> maybe use a relative tolerance
  Poly tL({-28594./5, 92231./15, -2271, 326, -101./12, 12./5, -7./10});
  Poly tM({ 498143./20,-288133./15,5012,-441});
  Poly tR({150942./5,-124586./5,14449./2,-767.,101./12,-12./5,7./10});

  // Compute the pieces and domain boundaries and compare results to their targets:
  Poly rL, rM, rR;            // left, middle, right section of result
  double rLL, rLU, rRL, rRU;  // lower and upper limits of the sections
  PiecePoly::convolvePieces(p, pL, pU, q, qL, qU, rL, rLL, rLU, rM, rR, rRL, rRU);
  r &= rsIsCloseTo(rL, tL, tol);
  r &= rsIsCloseTo(rM, tM, tol);
  r &= rsIsCloseTo(rR, tR, tol);
  r &= rLL == 2; r &= rLU == 3;
  r &= rRL == 5; r &= rRU == 6;

  // Test role reversal of p and q:
  PiecePoly::convolvePieces(q, qL, qU, p, pL, pU, rL, rLL, rLU, rM, rR, rRL, rRU);
  r &= rsIsCloseTo(rL, tL, tol);
  r &= rsIsCloseTo(rM, tM, tol);
  r &= rsIsCloseTo(rR, tR, tol);
  r &= rLL == 2; r &= rLU == 3;
  r &= rRL == 5; r &= rRU == 6;

  // Now, let's swap the domains:
  // p  = piecewise([(( 3,4), 2 - 3*x + 5*x^2 - 7*x^3)])
  // q  = piecewise([((-1,2), 3 - 4*x + 6*x^2)])
  // pq = p.convolution(q)
  // pq
  //
  // gives:
  // -7/10*x^6 + 12/5*x^5 - 101/12*x^4 - 38*x^3 + 561*x^2 - 53089/15*x + 26206/5 on (2, 3], 
  // -3037/2*x^2 + 178022/15*x - 478007/20 on (3, 5], 
  // 7/10*x^6 - 12/5*x^5 + 101/12*x^4 - 109*x^3 - 3319/2*x^2 + 72444/5*x - 142758/5 on (5, 6]
  p  = Poly({ 2,-3,5,-7 }); pL =  3, pU = 4; // p(x) = 2 - 3*x + 5*x^2 - 7*x^3  in  3..4
  q  = Poly({ 3,-4,6    }); qL = -1, qU = 2; // q(x) = 3 - 4*x + 6*x^2          in -1..2
  tL = Poly({26206./5, -53089./15, 561, -38, -101./12, 12./5, -7./10});
  tM = Poly({-478007./20, 178022./15, -3037./2 });
  tR = Poly({-142758./5, 72444./5, -3319./2, -109, 101./12, -12./5, 7./10});

  PiecePoly::convolvePieces(p, pL, pU, q, qL, qU, rL, rLL, rLU, rM, rR, rRL, rRU);
  r &= rsIsCloseTo(rL, tL, tol);
  r &= rsIsCloseTo(rM, tM, tol);
  r &= rsIsCloseTo(rR, tR, tol);
  r &= rLL == 2; r &= rLU == 3;
  r &= rRL == 5; r &= rRU == 6;

  PiecePoly::convolvePieces(q, qL, qU, p, pL, pU, rL, rLL, rLU, rM, rR, rRL, rRU);
  r &= rsIsCloseTo(rL, tL, tol);
  r &= rsIsCloseTo(rM, tM, tol);
  r &= rsIsCloseTo(rR, tR, tol);
  r &= rLL == 2; r &= rLU == 3;
  r &= rRL == 5; r &= rRU == 6;

  // ToDo: 
  // -make tests with all possible combinations of p having higher and lower degree than q and
  //  longer, shorter, equal domain
  // -figure out, why the nominal degrees of the results are higher than the actual ones (higher
  //  order coeffs are zero) - i suppose, it's because the integral with respect to y could 
  //  potentially produce a higher degree output, but because the lower right coeffs in matrix
  //  of the bivariate polynomial are zero, the final result coeffs come out as zero too
  //  ...maybe the function should cut off trailing zeros
  // -the middle section seems to have the degree of the polynomial with the longer domain, even if
  //  that's the one with the lower degree...seems like a high-degree polynomial gets its wiggles
  //  smoothed out such that only the low degree remains in the smoothed polynomial
  // -the outer sections seem to have a degree that is the product of the input degrees? 
  //  ..nope - it's their sum plus one
  // -figure out the nominal and actual degrees of the output segments as function of the degrees
  //  (and maybe lengths of domains?) of the input segments
  // -produce B-Spline polynomials and/or Irvin-Hall distribution

  return r;
}

bool testPiecewisePolynomial2()
{
  bool r = true;

  using Poly      = RAPT::rsPolynomial<double>;
  using PiecePoly = rsPiecewisePolynomial<double>;

  Poly one({ 1.0 });  // polynomial that is constantly 1
  Poly two({ 2.0 });  // polynomial that is constantly 2

  PiecePoly p;
  p.addPiece(Poly({1,-1      }), 0, 1);  // p(x) = 1-x     in x = 0..1
  p.addPiece(Poly({3, 0,-1   }), 1, 2);  // p(x) = 3-x^2   in x = 1..2
  p.addPiece(Poly({1, 0, 0,-1}), 2, 3);  // p(x) = 1-x^3   in x = 2..3

  // Test, if the left domain boundary is included and the right one is excluded in evaluation:
  double x, y;
  x = -0.1; y = p(x); r &= y == 0;
  x =  0.0; y = p(x); r &= y == 1-x;
  x =  0.1; y = p(x); r &= y == 1-x;
  x =  0.9; y = p(x); r &= y == 1-x;
  x =  1.0; y = p(x); r &= y == 3-x*x;
  x =  1.1; y = p(x); r &= y == 3-x*x;
  x =  1.9; y = p(x); r &= y == 3-x*x;
  x =  2.0; y = p(x); r &= y == 1-x*x*x;
  x =  2.1; y = p(x); r &= y == 1-x*x*x;
  x =  2.9; y = p(x); r &= y == 1-x*x*x;
  x =  3.0; y = p(x); r &= y == 0;

  p.addPiece(Poly({0, 1, 0}), 1, 2); // middle segment should now be  p(x) = 3 + x - x^2
  x = 1.5; y = p(x); r &= y == 3 + x - x*x;

  // initialization with 6 constant segments that is repeatedly used for the following tests:
  auto init6 = [&]() 
  {  
    p.clear();
    for(int i = 0; i < 6; i++)
      p.addPiece(one, double(i), i+1.0);
  };

  // left match, right mismatch:
  init6(); p.addPiece(two, 2.0, 4.5); 
  r &= p(1.9) == 1 && p(2.1) == 3 && p(4.4) == 3 && p(4.6) == 1;
  //plot(p);

  // left match, right beyond:
  init6(); p.addPiece(two, 2.0, 6.5); 
  r &= p(1.9) == 1 && p(2.1) == 3 && p(5.9) == 3 && p(6.1) == 2 && p(6.4) == 2 && p(6.6) == 0;
  //plot(p);

  // left beyond, right match:
  init6(); p.addPiece(two, -1.0, 3.0);
  r &= p(-1.1) == 0 && p(-0.9) == 2 && p(-0.1) == 2 && p(0.1) == 3 && p(2.9) == 3 && p(3.1) == 1;
  //plot(p);

  // left mismatch, right match
  init6(); p.addPiece(two, 2.5, 5.0);
  r &= p(2.4) == 1 && p(2.6) == 3 && p(4.9) == 3 && p(5.1) == 1;
  //plot(p);

  init6(); p.addPiece(two, 3.25, 3.75);
  r &= p(3.2) == 1 && p(3.3) == 3 && p(3.7) == 3 && p(3.8) == 1;
  //plot(p);

  // left mismatch, right mismatch
  init6(); p.addPiece(two, 2.5, 4.5);
  r &= p(2.4) == 1 && p(2.6) == 3 && p(4.4) == 3 && p(4.6) == 1;
  //plot(p);

  p = p + p;
  r &= p(2.4) == 2 && p(2.6) == 6 && p(4.4) == 6 && p(4.6) == 2;
  //plot(p);

  p = p - p;
  r &= p(2.4) == 0 && p(2.6) == 0 && p(4.4) == 0 && p(4.6) == 0;

  // todo: 
  // -implement and test multiplication by a scalar (left and right)
  // -implement another unit test that adds polynomials that arise from randomly sampling noise at 
  //  random x-values and interpolating them...or maybe just use totally random polynomials

  return r;
}

bool testPiecewisePolynomial3()
{
  bool r = true;

  int N1 = 50;  // number datapoints in polynomial 1
  int N2 = 50;  // number datapoints in polynomial 2

  using Poly      = RAPT::rsPolynomial<double>;
  using PiecePoly = rsPiecewisePolynomial<double>;

  rsNoiseGenerator<double> ng;
  ng.setRange(0, 1);
  PiecePoly p1, p2;

  auto getRandomPoly = [&](int deg) -> Poly 
  {
    Poly p(deg);
    double* c = p.getCoeffPointer();
    for(int i = 0; i <= deg; i++)
      c[i] = ng.getSample();
    return p;
  };

  auto getRandomPiecePoly = [&](int numPieces) -> PiecePoly
  {
    PiecePoly pp;
    double xL = ng.getSample();
    for(int i = 0; i < numPieces; i++)
    {
      double xR = xL + ng.getSample();
      int deg = (int)ceil((8.0 * ng.getSample()));
      Poly p = getRandomPoly(deg);
      p.shiftX((xL+xR)/2);
      pp.addPiece(p, xL, xR);
      xL = xR;
    }
    return pp;
  };

  p1 = getRandomPiecePoly(N1);
  p2 = getRandomPiecePoly(N2);
  //plot(p1);
  //plot(p2);
  PiecePoly sum = p1 + p2;
  //plot(sum); 


  // evaluate at random positions
  double xMin = sum.getDomainMinimum();
  double xMax = sum.getDomainMaximum();
  double tol  = 1.e-3;
  double dMax = 0;
  for(int i = 0; i < 100; i++)
  {
    double x = xMin + (xMax-xMin) * ng.getSample();
    double y = sum.evaluate(x);
    double t = p1.evaluate(x) + p2.evaluate(x);
    dMax = rsMax(dMax, rsAbs(y-t));
  }
  r &= dMax <= tol;
  // dMax is quite large: 0.000258... that's much more than what we would expect from rounding 
  // errors - or is it? should we expect such large errors? maybe it's because the error in the 
  // coefficient gets blown up by x^N? we have polynomials of degrees up to 9 and x goes up to 
  // 25 ...oh - and we get really big coeffs at higher x - maybe it would be numerically much 
  // better to store coefficients for a polynomial in 0..1 and transform the input, i.e. work with
  // normalized x values internally - but then, we can not just simply add coefficients anymore 
  // when we add a piece, because the coeffs are incompatible - each set of coeffs requires a 
  // different affine transform of the input - but actually, this class is not really meant for
  // interpolating random datapoints anyway. it's more meant for math stuff - in particular, doing
  // computations with the Irwin-Hall distribution to shape the properties of noise


  return r;
}

bool testPiecewisePolynomial()
{
  bool r = true;
  r &= testPiecewisePolynomial1();
  r &= testPiecewisePolynomial2();
  r &= testPiecewisePolynomial3();
  //r &= testPiecewisePolynomial3();
  return r;
}

bool testPolynomial()
{
  std::string reportString = "Polynomial"; // dummy -> remove
  bool testResult = true;

  testResult &= testConvolution(                              reportString);
  testResult &= testCubicCoeffsFourPoints(                    reportString);
  testResult &= testCubicCoeffsTwoPointsAndDerivatives(       reportString);
  testResult &= testPolynomialEvaluation(                     reportString);
  testResult &= testPolynomialDivision(                       reportString);
  testResult &= testPolynomialArgumentShift(                  reportString);
  testResult &= testPolynomialDiffAndInt(                     reportString);
  testResult &= testPolynomialFiniteDifference(               reportString);
  testResult &= testPolynomialComposition(                    reportString);
  testResult &= testPolynomialWeightedSum(                    reportString);
  testResult &= testPolynomialIntegrationWithPolynomialLimits(reportString);
  testResult &= testPolynomialInterpolation(                  reportString);
  testResult &= testPolynomialRootFinder(                     reportString);
  testResult &= testPartialFractionExpansion(                 reportString);
  testResult &= testPartialFractionExpansion2(                reportString);
  testResult &= testPolynomialBaseChange(                     reportString);
  testResult &= testPolynomialRecursion(                      reportString);
  testResult &= testJacobiPolynomials(                        reportString);
  testResult &= testSpecialPolynomials();
  testResult &= testQuadraticTo3Points();

  // under construction:
  testResult &= testPowersChebychevExpansionConversion(       reportString);

  // polynomial class:
  testResult &= testPolynomialOperators(                      reportString);
    // fails!

  testResult &= testRationalFunction(reportString);


  testResult &= testBivariatePolynomial();
  testResult &= testBivariatePolynomial2();
  testResult &= testTrivariatePolynomial(); // this takes long
  testResult &= testPiecewisePolynomial();

  return testResult;
}