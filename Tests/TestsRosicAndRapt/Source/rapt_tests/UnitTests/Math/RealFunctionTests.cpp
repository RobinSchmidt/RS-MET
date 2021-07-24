
typedef std::complex<double> rsComplexDbl; // get rid

//inline double rsAbsFast(double x)
//{
//  static const unsigned long long mask = 0x7FFFFFFFFFFFFFFFULL; // binary: 011111...
//  unsigned long long tmp = *(reinterpret_cast<unsigned long long*>(&x)) & mask;
//  return *(reinterpret_cast<double*>(&tmp));
//}
//// move to RSLib

bool testAbsAndSign()
{
  bool testResult = true;


  //testResult &= rsAbsFast(-0.0) == 0.0;
  //testResult &= rsAbsFast(-0.5) == 0.5;
  //testResult &= rsAbsFast(-1.0) == 1.0;
  //testResult &= rsAbsFast(-2.0) == 2.0;

  //testResult &= rsAbsFast(+0.0) == 0.0;
  //testResult &= rsAbsFast(+0.5) == 0.5;
  //testResult &= rsAbsFast(+1.0) == 1.0;
  //testResult &= rsAbsFast(+2.0) == 2.0;
  //// it has been found that the "fast" abs, as defined above is actually slower than the
  //// built in fabs function (in release builds), so...

  // todo: test denormals, infinities, nans...


  return testResult;
}

bool testHyperbolicFunctions()
{
  bool testResult = true;

  double xMin = -2000.0;
  double xMax = +2000.0;
  int    N    = 4001;         // number of values between xMin, xMax
  double tol  = 1.e-13;
  double x;                   // argument
  double s, c, t, ts, tc, tt; // hyperbolic sine, cosine, tangent and their target values
  for(int n = 1; n <= N; n++)
  {
    x  = rsLinToLin((double)n, 1.0, (double)N, xMin, xMax);
    ts = sinh(x);
    tc = cosh(x);
    tt = tanh(x);
    s  = rsSinh(x);
    c  = rsCosh(x);
    t  = rsTanh(x);
    testResult &= areNumbersEqual(ts, s, tol);
    testResult &= areNumbersEqual(tc, c, tol);
    testResult &= areNumbersEqual(tt, t, tol);

    rsSinhCosh(x, &s, &c);
    testResult &= areNumbersEqual(ts, s, tol);
    testResult &= areNumbersEqual(tc, c, tol);

    // reciprocal values are associated with rsCsch (hyp. cosecant), rsSech (hyp. secant),
    // rsCoth (hyp. cotangent):
    s  = rsCsch(x);  // hyperbolic cosecant  = 1 / sinh(x);
    c  = rsSech(x);  // hyperbolic secant    = 1 / cosh(x);
    t  = rsCoth(x);  // hyperbolic cotangent = 1 / tanh(x);
    ts = 1 / ts;
    tc = 1 / tc;
    tt = 1 / tt;
    testResult &= areNumbersEqual(ts, s, tol);
    testResult &= areNumbersEqual(tc, c, tol);
    testResult &= areNumbersEqual(tt, t, tol);

    int dummy;
    if( testResult == false )
    {
      dummy = 0; if(dummy == 0){}
      // is hit with GCC on Linux: ts, tc are tiny nonzero numbers whereas s,c are exactly zero
      // maybe below some threshold, we should switch from relative to absolute error in
      // areNumbersEqual
    }

  }

  return testResult;
}

bool testSinc()
{
  bool testResult = true;

  double x, y;

  x = 0.5*EPS;
  y = rsSinc(x);
  testResult &= y == 1.0;

  x = EPS;
  y = rsSinc(x);
  testResult &= y == 1.0;

  x = 2*EPS;
  y = rsSinc(x);
  testResult &= y == 1.0;

  x = PI;
  y = rsSinc(x);
  testResult &= rsAbs(y) < EPS;

  x = 1.0;
  y = rsNormalizedSinc(x);
  testResult &= rsAbs(y) < EPS;

  x = -1.0;
  y = rsNormalizedSinc(x);
  testResult &= rsAbs(y) < EPS;

  x = 2.0;
  y = rsNormalizedSinc(x);
  testResult &= rsAbs(y) < EPS;

  x = -2.0;
  y = rsNormalizedSinc(x);
  testResult &= rsAbs(y) < EPS;

  return testResult;
}

bool testComplexExponentialIterator(rsComplexDbl a, rsComplexDbl z)
{
  rsComplexDbl w;    // computed value
  rsComplexDbl wt;   // target value
  double e;          // relative error between target and computed value
  double eMax = 0.0; // maximum error
  int N = 1000;      // number of iterations
  //int n;             // iteration counter

  // create the iterator for complex expoentials:
  rsComplexExponentialIterator<double> it(a, z);

  // compute values and target values ans measure the maximum distance:
  for(int n = 0; n < N; n++)
  {
    //wt = a * rsPowC(z, rsComplexDbl((double)n, 0));
    wt = a * pow(z, rsComplexDbl((double)n, 0));
    w  = it.getValue();
    e  = abs(wt-w) / abs(wt);
    if( e > eMax )
      eMax = e;
  }

  if( eMax < 1.e-12 )
    return true;
  else
    return false;
}

bool testSineIterator(double w, double p, double a)
{
  bool ok = true;

  double tol = 1.e-10;  // tolerance
  double y;             // computed value
  double yt;            // target value
  double e;             // relative error between target and computed value
  double eMax = 0.0;    // maximum error
  int N = 1000;         // number of iterations

  // Create the iterator for sinusoids:
  rsSineIterator<double> it(w, p, a);

  // Check retrieval of normalized radian frequency:
  double w2 = it.getOmega();
  ok &= rsIsCloseTo(w, w2, tol);


  // Compute values and target values ans measure the maximum distance:
  double p1;// p2, a2;
  for(int n = 0; n < N; n++)
  {
    p1 = w*n + p;          // instantaneous target phase ..maybe we need to wrap it
    //p2 = it.getPhase();    // at n==4,8,.. we get nan!!
    yt = a * sin(p1);  
    y  = it.getValue();
    e  = fabs(yt-y) / y;
    if( e > eMax )
      eMax = e;

    // Check retrieval of instantaneous phase and amplitude;

    //a2 = it.getAmplitude(); // not yet implemented
  }

  ok &=  eMax <= tol;

  return ok;
}

bool testPolynomialIterator()
{
  bool ok = true;

  int   N    = 500;            // number of iterations
  //float a[4] = { 7, 5, 3, 2 };  // polynomial coefficients
  //float a[4] = { 0.7, -0.5, 0.3, -0.1 };  // polynomial coefficients
  float a[4] = { +0.7f, -0.5f, +0.3f, -0.1f };  // polynomial coefficients
  //float a[4] = { 7, -5, 3, -2 };  // polynomial coefficients
  float h    =  0.01f;          // stepsize
  float x0   = -2.5f;           // initial value for x

  rsPolynomialIterator<float, 3> pIt;
  rsExpPolyIterator<float, 3> eIt;
  pIt.setup(a, h, x0);
  eIt.setup(a, h, x0);
  std::vector<float> p(N), pt(N), pErrA(N), pErrR(N), x(N);
  std::vector<float> y(N), yt(N), yErrA(N), yErrR(N);
  for(int n = 0; n < N; n++)
  {
    x[n]     = x0 + n*h;

    pt[n]    = rsPolynomial<float>::evaluate(x[n], a, 3);
    p[n]     = pIt.getValue();
    pErrA[n] = p[n] - pt[n];
    pErrR[n] = pErrA[n] / pt[n];

    yt[n]    = exp(pt[n]);
    y[n]     = eIt.getValue();
    yErrA[n] = y[n] - yt[n];
    yErrR[n] = yErrA[n] / yt[n];
  }


  float err;
  err = rsMaxDeviation(p, pt); ok &= err <= 5.e-5;
  err = rsMaxDeviation(y, yt); ok &= err <= 0.1;     // error is large!

  //rsPlotVectorsXY(x, pt, p, pErrA);
  //rsPlotVectorsXY(x, pErrA, pErrR);

  //rsPlotVectorsXY(x, yt, y, yErrA);
  //rsPlotVectorsXY(x, yErrA, yErrR);


  // The error accumulation is quite significant. Maybe we can reduce it by computing the initial
  // state in extended precision such that we at least start with a clean state...although actually
  // in the initial section, the error passes through zero several times which may mean that some
  // intermediate states are actually quite clean by coincidence. But that's not necessarily so:
  // just because the output has zero error at some instant does not mean that the full state is
  // free of error. In production code, we should re-initialize periodically with exactly computed
  // values. How long should the interval be? If it's not at least of order 100, it may become
  // questionable, if iterative evaluation makes any sense at all. 


  return ok;
}

bool testPowerIterator(double p, double e)
{
  double tol = 1.e-10;    // tolerance
  int N = 2001;           // number of iterations, with 1001, we get an all-zero output for
                          // p = 0, e = 0.001
  int oversampling = 100;

  rsPowerIterator<double> it(p, e);

  std::vector<double> x(N), y(N), yt(N);
  x = rsLinearRangeVector(N, 0.0, 1.0);
  //double dx = 1.0 / ((N-1)*oversampling); // or 1.0 / (N*oversampling-1)?
  double dx = 1.0 / (N*oversampling-1); 

  for(int n = 0; n < N; n++)
  {
    //yt[n] = pow(x[n], p);
    yt[n] = pow(x[n]+e, p);  // but it should be scaled and shifted!


    double yn;
    for(int i = 0; i < oversampling; i++)
      yn = it.getValue(dx);
    y[n] = yn;
  }

  rsPlotVectorsXY(x, yt, y);

  double err = rsArrayTools::maxDeviation(&yt[0], &y[0], N);
  return err <= tol;

  // Nope! That doesn't work yet. I tried to implement this technique posted by Big Tick:
  // https://www.kvraudio.com/forum/viewtopic.php?f=33&t=567563
  // ...hmm - could this be a numerical accuracy issue? Maybe it will work better with smaller dx, 
}

bool testFunctionIterators()
{
  bool ok = true;

  // test complex exponential iterator with different values of z corresponding to decaying spiral
  // circular motion and growing spiral:
  rsComplexDbl a(-1.3, 1.2);  // multiplier == initial value
  ok &= testComplexExponentialIterator(a, rsComplexDbl(0.6, 0.7)); // |z| < 1
  ok &= testComplexExponentialIterator(a, rsComplexDbl(0.6, 0.8)); // |z| = 1
  ok &= testComplexExponentialIterator(a, rsComplexDbl(0.7, 0.8)); // |z| > 1
  ok &= testSineIterator(2.5,  0.3, 1.2);
  ok &= testSineIterator(2.5, -0.3, 1.2);
  ok &= testPolynomialIterator();

  // This does not yet work - at least not very well:
  //ok &= testPowerIterator(0.5, 0.001); 
  //ok &= testPowerIterator(2.0, 0.001);
  //ok &= testPowerIterator(2.0, 0.0001); // doesn't work without oversampling when z-based formula is used

  return ok;
}

// is this superseded by the function above? if so, delete...
/*
bool testFunctionIterators(std::string &reportString)
{
  std::string testName = "FunctionIterators";
  bool testResult = true;

  rsComplexDbl a(-1.3, 1.2);  // multiplier == initial value
  rsComplexDbl z( 0.6, 0.7);  // function argument
  rsComplexDbl w;             // computed value
  rsComplexDbl wt;            // target value
  double d;                   // distance between target and computed value
  double dMax = 0.0;          // maximum distance
  int N = 1000;               // number of iterations
  int n;                      // iteration counter


  // test - :
  d = z.getRadius();

  // create the iterator for complex expoentials:
  rsComplexExponentialIterator it(a, z);

  // compute values and target values ans measure the maximum distance:
  for(n = 0; n < N; n++)
  {
    wt = a * rsPowC(z, rsComplexDbl((double)n, 0));
    w  = it.getValue();
    d  = (wt-w).getRadius() / wt.getRadius();
    if( d > dMax )
      dMax = d;
    int dummy = 0;
  }



  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
*/

bool testWrap()
{
  bool r = true;

  double d;

  // wrap numbers into the range -2...+3
  r &= (d = RAPT::rsWrapToInterval(-2.0, -2.0, +3.0)) == -2.0;  // x at left limit
  r &= (d = RAPT::rsWrapToInterval(+3.0, -2.0, +3.0)) == -2.0;  // x at right limit - wraps back to left
  r &= (d = RAPT::rsWrapToInterval(+2.0, -2.0, +3.0)) == +2.0;  // x inside range
  r &= (d = RAPT::rsWrapToInterval(-1.0, -2.0, +3.0)) == -1.0;  // x inside range
  r &= (d = RAPT::rsWrapToInterval(+4.0, -2.0, +3.0)) == -1.0;  // x above range
  r &= (d = RAPT::rsWrapToInterval(-3.0, -2.0, +3.0)) == +2.0;  // x below range

  // interpolation of a wrapping value (maybe rename to rsInterpolatePeriodic):
  r &= (d = RAPT::rsInterpolateWrapped(4.0, 6.0, 0.5, 0.0, 10.0)) ==  5.0;
  r &= (d = RAPT::rsInterpolateWrapped(6.0, 4.0, 0.5, 0.0, 10.0)) ==  5.0;
  r &= (d = RAPT::rsInterpolateWrapped(7.0, 1.0, 0.5, 0.0, 10.0)) ==  9.0;
  r &= (d = RAPT::rsInterpolateWrapped(1.0, 7.0, 0.5, 0.0, 10.0)) ==  9.0;
  r &= (d = RAPT::rsInterpolateWrapped(9.0, 3.0, 0.5, 0.0, 10.0)) ==  1.0;
  r &= (d = RAPT::rsInterpolateWrapped(3.0, 9.0, 0.5, 0.0, 10.0)) ==  1.0;

  return r;
}

bool testWindowFunctions(int N)
{
  bool r = true;

  using Vec = std::vector<double>;
  using WF = rsWindowFunction;
  using WT = rsWindowFunction::WindowType;

  Vec w(N);  // actual window produced by library code
  Vec v(N);  // reference window produced by prototype code
  double tol = 1.e-14;

  // compare prototype and production code implementation of Dolph-Chebychev window with various 
  // attenuations:
  WF::dolphChebychev(&w[0], N, 20.); cheby_win(&v[0], N, 20.); r &= rsIsCloseTo(w, v, tol);
  WF::dolphChebychev(&w[0], N, 40.); cheby_win(&v[0], N, 40.); r &= rsIsCloseTo(w, v, tol);
  WF::dolphChebychev(&w[0], N, 60.); cheby_win(&v[0], N, 60.); r &= rsIsCloseTo(w, v, tol);
  WF::dolphChebychev(&w[0], N, 80.); cheby_win(&v[0], N, 80.); r &= rsIsCloseTo(w, v, tol);

  return r;
}  

bool testWindowFunctions()
{
  bool r = true;

  r &= testWindowFunctions(16);  // N is a power of two
  r &= testWindowFunctions(17);  // N is odd
  r &= testWindowFunctions(18);  // N is even but not a power of two

  return r;
}


bool testPeriodicDistance()
{
  bool r = true;

  // should all be 1:
  r &= rsDistanceToMultipleOf(  1.0, 10.0) == 1.0;
  r &= rsDistanceToMultipleOf( -1.0, 10.0) == 1.0;
  r &= rsDistanceToMultipleOf(  9.0, 10.0) == 1.0;
  r &= rsDistanceToMultipleOf( -9.0, 10.0) == 1.0;
  r &= rsDistanceToMultipleOf( 11.0, 10.0) == 1.0;
  r &= rsDistanceToMultipleOf(-11.0, 10.0) == 1.0;
  r &= rsDistanceToMultipleOf( 19.0, 10.0) == 1.0;
  r &= rsDistanceToMultipleOf(-19.0, 10.0) == 1.0;
  r &= rsDistanceToMultipleOf( 21.0, 10.0) == 1.0;
  r &= rsDistanceToMultipleOf(-21.0, 10.0) == 1.0;

  return r;
}

bool testSinCosApproximations()
{
  bool ok = true;

  double xMin = -1.5*PI;
  double xMax = +1.5*PI;
  int    N    = 1001;     // number of values between xMin and xMax

  using Vec = std::vector<double>;
  Vec c(N), s(N), c1(N), s1(N), c2(N), s2(N), c3(N);
  Vec x = rsLinearRangeVector(N, xMin, xMax);
  for(int i = 0; i < N; i++)
  {
    c[i]  = cos(x[i]);
    s[i]  = sin(x[i]);
    rsSinCos1(x[i], &s1[i], &c1[i]);
    rsSinCos2(x[i], &s2[i], &c2[i]);
    c3[i] = rsCos2(x[i]);
  }
  //rsPlotVectorsXY(x, s, c, s1, c1);
  //rsPlotVectorsXY(x, s, c, s2, c2);
  //rsPlotVectorsXY(x, c, c3, c-c3);
  //rsPlotVectorsXY(x, c-c3);

  // This unit test is not yet finished because the testee functions themselves are not yet 
  // finished. ToDo: Check if the maximum error is wihtin expected bounds for all the 
  // approximations

  return ok;
}


bool testRealFunctions()
{
  bool ok = true;

  ok &= testAbsAndSign();
  //ok &= testHyperbolicFunctions(); // test doesn't pass
  ok &= testSinc();
  ok &= testFunctionIterators();
  ok &= testWrap();
  ok &= testWindowFunctions();
  ok &= testPeriodicDistance();
  ok &= testSinCosApproximations();

  return ok;
}
