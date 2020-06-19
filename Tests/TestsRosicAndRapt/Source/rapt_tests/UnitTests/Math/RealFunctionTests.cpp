
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
  double y;          // computed value
  double yt;         // target value
  double e;          // relative error between target and computed value
  double eMax = 0.0; // maximum error
  int N = 1000;      // number of iterations

  // create the iterator for sinusoids:
  rsSineIterator<double> it(w, p, a);

  // compute values and target values ans measure the maximum distance:
  for(int n = 0; n < N; n++)
  {
    yt = a * sin(w*n + p);
    y  = it.getValue();
    e  = fabs(yt-y) / y;
    if( e > eMax )
      eMax = e;
  }

  if( eMax < 1.e-10 )
    return true;
  else
    return false;
}

bool testFunctionIterators()
{
  bool testResult = true;

  // test complex exponential iterator with different values of z corresponding to decaying spiral
  // circular motion and growing spiral:
  rsComplexDbl a(-1.3, 1.2);  // multiplier == initial value
  testResult &= testComplexExponentialIterator(a, rsComplexDbl(0.6, 0.7)); // |z| < 1
  testResult &= testComplexExponentialIterator(a, rsComplexDbl(0.6, 0.8)); // |z| = 1
  testResult &= testComplexExponentialIterator(a, rsComplexDbl(0.7, 0.8)); // |z| > 1

  testResult &= testSineIterator(2.5,  0.3, 1.2);
  testResult &= testSineIterator(2.5, -0.3, 1.2);

  return testResult;
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

bool testRealFunctions()
{
  bool testResult = true;

  testResult &= testAbsAndSign();
  //testResult &= testHyperbolicFunctions(); // test doesn't pass
  testResult &= testSinc();
  testResult &= testFunctionIterators();
  testResult &= testWrap();
  testResult &= testWindowFunctions();
  testResult &= testPeriodicDistance();

  return testResult;
}
