
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
  rsPolynomial<double>::coeffsForShiftedArgument(p, q, order, x0);

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

  rsPolynomial<double>::derivative(a, ad, 5);
  testResult &= (ad[0] == -1);
  testResult &= (ad[1] == 10);
  testResult &= (ad[2] == 21);
  testResult &= (ad[3] == -12);
  testResult &= (ad[4] == 10);

  rsPolynomial<double>::integral(a, ai, 5, 2.0);
  testResult &= (ai[0] ==  2.0);
  testResult &= (ai[1] ==  2.0/1.0);
  testResult &= (ai[2] == -1.0/2.0);
  testResult &= (ai[3] ==  5.0/3.0);
  testResult &= (ai[4] ==  7.0/4.0);
  testResult &= (ai[5] == -3.0/5.0);
  testResult &= (ai[6] ==  2.0/6.0);

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

  static const int na = 5;
  static const int nb = 4;
  static const int nc = na*nb;
  double a[na+1] = {2, -1, 5,  7, -3, 2}; // 2*x^5 - 3*x^4 + 7*x^3 + 5*x^2 - 1*x^1 + 2*x^0
  double b[nb+1] = {3,  1, 4, -5,  3};    //         3*x^4 - 5*x^3 + 4*x^2 + 1*x^1 - 3*x^0
  double c[nc+1];
  rsPolynomial<double>::compose(a, na, b, nb, c);

  // check, if the composed c-polynomial returns the same result as applying the 2nd b-polynomial
  // to the result of the 1st a-polynomial:
  double x = -3.0; // input value
  double y1, y2;
  y1 = rsPolynomial<double>::evaluate(x,  a, na);
  y1 = rsPolynomial<double>::evaluate(y1, b, nb);
  y2 = rsPolynomial<double>::evaluate(x,  c, nc);
  testResult &= (y1 == y2);

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


// these 4 should go into rsPolynomial:

/** Evaluates the Newton polynomial:
  p(x) = c[0] + c[1]*(x-r[0]) + c[2]*(x-r[0])*(x-r[1]) + ... + cN*(x-r[0])*...*(x-r[N-1])
with the Newton coeffs c[0],...,c[N] and roots r[0],...,r[N-1].  */
template<class T>
T evaluateNewton(const T& x, const T* c, const T* r, int N)
{
  T pi = T(1);   // accumulator for Newton basis polynomials
  T y  = c[0];   // accumulator for result
  for(int i = 0; i < N; i++) {
    pi *= x - r[i];
    y  += c[i+1] * pi; }
  return y;
}
// c is of length N+1, r is of length N...verify this!
// should go to Evaluation section

/** Computes coefficients for Newton basis polynomials 1,(x-x[0]),(x-x[0])*(x-x[1]),... using 
divided differences from N pairs (x[i],y[i]), i = 0,...N-1. It overwrites the y-array with the 
coefficients. @see evaluateNewton. */
template<class T>
void newtonCoeffs(const T* x, T* y, int N)
{
  for(int i = 0; i < N; i++)
    for(int j = i+1; j < N; j++)
      y[j] = (y[j] - y[i]) / (x[j] - x[i]);
}
// should go to Coefficient generation section
// https://en.wikipedia.org/wiki/Polynomial_interpolation#Non-Vandermonde_solutions
// https://en.wikipedia.org/wiki/Neville%27s_algorithm

/** Given the coefficents for the Newton basis polynomials in a[i], i = 0,..,N-1, this function 
converts the coefficients to the monomial basis, i.e. into regular polynomial coefficients. The 
result is stored in the same array a. The array x is destroyed during the process - it's re-used 
internally as temporary buffer for intermediate results to allow higher level code optimize memory
usage (otherwise, the function would need an additional buffer of length N). */
/*
template<class T>
void newtonToMonomialCoeffs(T* x, T* a, int N)
{
  T x0 = x[0]; 
  x[0] = T(1);                   // x is re-used as our convolutive accumulator
  for(int i = 1; i < N; i++) {
    T x1 = x[i];                                            // save x[i] because the next line..
    rsArrayTools::convolveWithTwoElems(x, i, -x0, T(1), x); // ..overwrites x up to x[i] but we..
    x0 = x1;                                                // ..still need it in next iteration
    for(int j = 0; j < i; j++)
      a[j] += a[i] * x[j]; }
}
// should go to Conversions section
*/

/** In place version. Overwrites x,y arrays during the process. On return, y will contain the 
polynomial coeffs and x will contain garbage (more specifically, x will contain the coefficients of
the unique monic polynomial that has roots at the given original x-values - i don't think that's 
useful for the caller, but anyway). ...verify this! i think, the last root may be missing. */
template<class T>
void interpolantNewtonInPlace(T* x, T* y, int N)
{
  newtonCoeffs(      x, y, N);   // overwrites y
  rsPolynomial<T>::newtonToMonomialCoeffs(x, y, N);   // overwrites y again and also x
}
// should got Fitting/Interpolation section


// stuff below should be replaced by convenience functions that make use of the functions above or 
// maybe removed:

/*
template<class T>
void newtonToMonomialCoeffs(const T* c, const T* x, int N, T* a, T* b)
{
  b[0] = 1;     // b is our convolutive accumulator
  a[0] = c[0];
  for(int i = 1; i < N; i++) {
    rsArrayTools::convolveWithTwoElems(b, i, -x[i-1], T(1), b);
    for(int j = 0; j < i; j++)
      a[j] += c[i] * b[j];
    a[i] = c[i]; }
}
*/
// -a is allowed to be the same array as c
// -maybe use a single array for I/O, i.e. c/a - then we can also remove the a[0] = c[0] and
//  a[i] = c[i] statements
// -is it possible to implement it in a way that allows b to be equal to x? would that be useful?
//  we could do xi = x[0] before the loop, in the loop convolve with xi instead of x[i-1] and 
//  re-assign xi = x[i] after convolution
// -this way, we could implement the whol algo using only 2 arrays of length N: x,y which would 
//  contain g,a after return (where g stands for garbage, a for the coeffs)
// -maybe get rid (or call the in-place version)

/*
template<class T>
void newtonToMonomialCoeffs(const T* c, const T* x, int N, T* a)
{
  std::vector<T> b(N);   // workspace
  newtonToMonomialCoeffs(c, x, N, a, &b[0]);
}
*/
// c: Newton basis coeffs, x: x-coordinates of evaluation/data points, a: monomial basis coeffs, 
// N: number of data points, i.e. length of x,c,a 
// todo: make order of arguments consistent: inputs first, outputs last...also make the meaning of
// N consistent with typical use in rsPolynomial
// -let the use pass a workspace - avoid allocating memory for b
// -in a higher level implementation (like interpolantNewton below), we could use the passed 
//  y-array to be used as workspace during conversion, if it's ok to destroy the y-array
// -get rid! use in place versions only!






/** Implements Newton's algorithm using diveded difference to find the polynomial interpolant 
through the N points (x[n], y[n]), n = 0,...,N-1.  */
template<class T>
void interpolantNewton(T* a, const T* x, const T* y, int N)
{
  rsArrayTools::copy(y, a, N);
  newtonPolyCoeffs(x, a, N);

  // -from this point on, the y-array is not need anymore - maybe, it can be re-used as workspace 
  //  in the following steps to avoid allocating extra space, when it's ok to destroy the y-array
  // -or we could re-use y for the c-array above, like:
  //    newtonPolyCoeffs(y, x, y, N);
  //    newtonToMonomialCoeffs(y, x, N, y);
  //  then y would contain the y-values on entry and the polynomial coeffs on exit

  newtonToMonomialCoeffs(a, x, N, a);
}
// maybe make a class rsVectorTools that just contains convenience functions for the functions from
// rsArrayTools, such that we don't need the ugly &b[0] syntax and maybe can get rid of the length
// parameters (because vectors know their lengths)...but maybe it should also include






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
  newtonCoeffs(x, c, N);
  for(n = 0; n < N; n++) {
    yc[n] = evaluateNewton(x[n], c, x, N-1);
    testResult &= rsIsCloseTo(yc[n], y[n], tol); }

  /*
  // test converting coeffs for Newton basis into monomial basis:
  rsArrayTools::fillWithNaN(a, N);
  Poly::newtonToMonomialCoeffs(c, x, N, a);
  for(n = 0; n < N; n++) {
    yc[n] = Poly::evaluate(x[n], a, N-1);
    testResult &= rsIsCloseTo(yc[n], y[n], tol); }

  // test in-place conversion:
  rsArrayTools::copy(c, a, N);
  newtonToMonomialCoeffs(a, x, N, a);
  for(n = 0; n < N; n++) {
    yc[n] = Poly::evaluate(x[n], a, N-1);
    testResult &= rsIsCloseTo(yc[n], y[n], tol); }

  // test computing the monomial coeffs via Newton basis:
  interpolantNewton(a, x, y, N);
  for(n = 0; n < N; n++) {
    yc[n] = Poly::evaluate(x[n], a, N-1);
    testResult &= rsIsCloseTo(yc[n], y[n], tol); }
    // this is the 5th time, this loop occurs -> wrap into (lambda)function
  */

  // test in-place implementation of finding polynomial coeffs via Newton expansion:
  rsArrayTools::copy(x, c, N);
  rsArrayTools::copy(y, a, N);
  interpolantNewtonInPlace(c, a, N);  // after the call, a will conatin the desired coeffs 
  for(n = 0; n < N; n++) {
    yc[n] = Poly::evaluate(x[n], a, N-1);
    testResult &= rsIsCloseTo(yc[n], y[n], tol); }
  // maybe also check, if c contains coeffs of polynomial with roots at the original x[n]


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


//bool testPolynomial(std::string &reportString)
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


  return testResult;
}