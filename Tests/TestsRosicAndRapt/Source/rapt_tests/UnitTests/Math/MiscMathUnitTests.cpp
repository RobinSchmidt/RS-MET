
bool testExponentialCurveFitting(std::string &reportString)
{
  std::string testName = "ExponentialCurveFitting";
  bool testResult = true;

  // example from Hamming's "Numerical Methods ..." page 620
  static const int N = 4;
  static const int k = N/2;
  //double x[N] = { 0.0,  1.0,  2.0,  3.0};
  double y[N] = {32.0, 20.0, 14.0, 11.0};

  // find exponents and weights:
  double a[k];   // exponents
  double A[k];   // weights
  bool success = rsCurveFitter::fitExponentialSum(y, N, A, a, k);
  rsAssert(success);

  // now we should have: y[n] = A[0]*exp(a[0]*n) + A[1]*exp(a[1]*n) - verify this:
  double yc[N];
  int n;
  for(n = 0; n < N; n++)
    yc[n] = A[0]*exp(a[0]*n) + A[1]*exp(a[1]*n);
  for(n = 0; n < N; n++)
    testResult &= (yc[n] == y[n]);

  return testResult;
}

bool testRootFinding(std::string &reportString)
{
  std::string testName = "RootFinding";
  bool testResult = true;

  double r;
  double tol = std::numeric_limits<double>::epsilon();

  UnivariateScalarFunctionViaPointer<double> sine(&sin, &cos);

  // Newton iteration:
  r = sine.findRootViaNewtonNonRobust(3.0);
  testResult &= rsIsCloseTo(r, PI, tol);
  r = sine.findRootViaNewtonNonRobust(6.0);
  testResult &= rsIsCloseTo(r, 2*PI, tol);

  // Chebychev method (generalization of Newton iteration using also 2nd derivative):
  // ...stopping criterion still not working - so it's commented:
  //r = sine.findRootViaChebychevNonRobust(3.0);
  //testResult &= rsIsCloseTo(r, PI, tol);
  //r = sine.findRootViaChebychevNonRobust(6.0);
  //testResult &= rsIsCloseTo(r, 2*PI, tol);

  // Ridders' method does not yet converge here, so it's commented
  //r = sine.findRootViaRidders(6.2, 6.3);
  //r = sine.findRootViaRidders(3.1, 3.2);
  //r = sine.findRootViaRidders(2.5, 3.5);

  return testResult;
}

bool testGradientBasedOptimization(std::string &reportString)
{
  std::string testName = "GradientBasedOptimization";
  bool testResult = true;

  // set up the minimizer:
  double tol = 0.00001; // tolerance - later, pass this to the minimizer
  GradientBasedMinimizer<double> minimizer;
  minimizer.setBetaFormula(GradientBasedMinimizer<double>::POLAK_RIBIERE);
  //minimizer.setBetaFormula(GradientBasedMinimizer::FLETCHER_REEVES);
  //minimizer.setBetaFormula(GradientBasedMinimizer::HESTENES_STIEFEL);
  //minimizer.setPrintInfo(true);

  // create the error function object:
  QuadraticTestErrorFunction<double> error;
  double aMin[2]  = { 2, -2};  // x = [ 2, -2] is the desired minimum
  double aInit[2] = {-2, -2};  // x = [-2, -2] is the intial guess
  rsVectorDbl xMin( 2, aMin);
  rsVectorDbl xInit(2, aInit);
  rsVectorDbl xFinal;

  // do the minimization using different algorithms:
  minimizer.setAlgorithm(GradientBasedMinimizer<double>::GRADIENT_DESCENT);
  xFinal = minimizer.minimizeFunction(&error, xInit);
  testResult &= (xFinal-xMin).getEuclideanNorm() < tol;

  minimizer.setAlgorithm(GradientBasedMinimizer<double>::BOLD_DRIVER_WITH_MOMENTUM);
  xFinal = minimizer.minimizeFunction(&error, xInit);
  testResult &= (xFinal-xMin).getEuclideanNorm() < tol;

  minimizer.setAlgorithm(GradientBasedMinimizer<double>::CONJUGATE_GRADIENT);
  xFinal = minimizer.minimizeFunction(&error, xInit);
  testResult &= (xFinal-xMin).getEuclideanNorm() < tol;

  minimizer.setAlgorithm(GradientBasedMinimizer<double>::SCALED_CONJUGATE_GRADIENT);
  xFinal = minimizer.minimizeFunction(&error, xInit);
  testResult &= (xFinal-xMin).getEuclideanNorm() < tol;

  return testResult;
}

bool testMinSqrDifFixSum(std::string &reportString)
{
  // code moved to experiments - todo: implement actual unit-tests
  std::string testName = "MinSqrDifFixSum";
  bool testResult = true;


  return testResult;
}

bool testPhaseUnwrapStuff(std::string &reportString)  // rename to testUnwrapping
{
  bool r = true;  // test result

  // we consider the range -3..+7 with wrap-around - all numbers ineach column from an equivalence
  // class and numbers in the leftmost column are identified with numbers in the rightmost column 
  // (for example -PI and PI would represent the same angle):
  // -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -13   k = -2
  // -13 -12 -11 -10 -09 -08 -07 -06 -05 -04 -03   k = -1
  // -03 -02 -01 +00 +01 +02 +03 +04 +05 +06 +07   k =  0 (base range)
  // +07 +08 +09 +10 +11 +12 +13 +14 +15 +16 +17   k = +1
  // +17 +18 +19 +20 +21 +22 +23 +24 +25 +26 +27   k = +2
   

  double rangeMin = -3.0;
  double rangeMax = +7.0;
  double val;

  val = rsConsistentUnwrappedValue(-21.0, 5.0, rangeMin, rangeMax); r &= val == -25.0;

  val = rsConsistentUnwrappedValue(-20.0, 5.0, rangeMin, rangeMax); r &= val == -15.0;
  val = rsConsistentUnwrappedValue(-19.0, 5.0, rangeMin, rangeMax); r &= val == -15.0;
  val = rsConsistentUnwrappedValue(-18.0, 5.0, rangeMin, rangeMax); r &= val == -15.0;
  val = rsConsistentUnwrappedValue(-17.0, 5.0, rangeMin, rangeMax); r &= val == -15.0;
  val = rsConsistentUnwrappedValue(-16.0, 5.0, rangeMin, rangeMax); r &= val == -15.0;
  val = rsConsistentUnwrappedValue(-15.0, 5.0, rangeMin, rangeMax); r &= val == -15.0;
  val = rsConsistentUnwrappedValue(-14.0, 5.0, rangeMin, rangeMax); r &= val == -15.0;
  val = rsConsistentUnwrappedValue(-13.0, 5.0, rangeMin, rangeMax); r &= val == -15.0;
  val = rsConsistentUnwrappedValue(-12.0, 5.0, rangeMin, rangeMax); r &= val == -15.0;
  val = rsConsistentUnwrappedValue(-11.0, 5.0, rangeMin, rangeMax); r &= val == -15.0;

  val = rsConsistentUnwrappedValue(-10.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsConsistentUnwrappedValue(- 9.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsConsistentUnwrappedValue(- 8.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsConsistentUnwrappedValue(- 7.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsConsistentUnwrappedValue(- 6.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsConsistentUnwrappedValue(- 5.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsConsistentUnwrappedValue(- 4.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsConsistentUnwrappedValue(- 3.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsConsistentUnwrappedValue(- 2.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsConsistentUnwrappedValue(- 1.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;

  // interestingly, this block has 11 entries, the others have 10:
  val = rsConsistentUnwrappedValue(  0.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsConsistentUnwrappedValue(  1.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsConsistentUnwrappedValue(  2.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsConsistentUnwrappedValue(  3.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsConsistentUnwrappedValue(  4.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsConsistentUnwrappedValue(  5.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsConsistentUnwrappedValue(  6.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsConsistentUnwrappedValue(  7.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsConsistentUnwrappedValue(  8.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsConsistentUnwrappedValue(  9.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsConsistentUnwrappedValue( 10.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;

  val = rsConsistentUnwrappedValue(11.0, 5.0, rangeMin, rangeMax); r &= val == 15.0;
  val = rsConsistentUnwrappedValue(12.0, 5.0, rangeMin, rangeMax); r &= val == 15.0;
  val = rsConsistentUnwrappedValue(13.0, 5.0, rangeMin, rangeMax); r &= val == 15.0;
  val = rsConsistentUnwrappedValue(14.0, 5.0, rangeMin, rangeMax); r &= val == 15.0;
  val = rsConsistentUnwrappedValue(15.0, 5.0, rangeMin, rangeMax); r &= val == 15.0;
  val = rsConsistentUnwrappedValue(16.0, 5.0, rangeMin, rangeMax); r &= val == 15.0;
  val = rsConsistentUnwrappedValue(17.0, 5.0, rangeMin, rangeMax); r &= val == 15.0;
  val = rsConsistentUnwrappedValue(18.0, 5.0, rangeMin, rangeMax); r &= val == 15.0;
  val = rsConsistentUnwrappedValue(19.0, 5.0, rangeMin, rangeMax); r &= val == 15.0;
  val = rsConsistentUnwrappedValue(20.0, 5.0, rangeMin, rangeMax); r &= val == 15.0;

  val = rsConsistentUnwrappedValue(21.0, 5.0, rangeMin, rangeMax); r &= val == 25.0;
  // ...and so on...


  // now, we consider the range 0..5 and use the function that supposes a zero lower limit

  // 00 01 02 03 04 05    k = 0 (base range)
  // 05 06 07 08 09 10    k = 1
  // 10 11 12 13 14 15    k = 2
  // 15 16 17 18 19 20    k = 3

  rangeMax = 5.0;
  val = rsConsistentUnwrappedValue0( 0.0, 2.0, rangeMax); r &= val ==  2.0;
  val = rsConsistentUnwrappedValue0( 1.0, 2.0, rangeMax); r &= val ==  2.0;
  val = rsConsistentUnwrappedValue0( 2.0, 2.0, rangeMax); r &= val ==  2.0;
  val = rsConsistentUnwrappedValue0( 3.0, 2.0, rangeMax); r &= val ==  2.0;
  val = rsConsistentUnwrappedValue0( 4.0, 2.0, rangeMax); r &= val ==  2.0;

  val = rsConsistentUnwrappedValue0( 5.0, 2.0, rangeMax); r &= val ==  7.0;
  val = rsConsistentUnwrappedValue0( 6.0, 2.0, rangeMax); r &= val ==  7.0;
  val = rsConsistentUnwrappedValue0( 7.0, 2.0, rangeMax); r &= val ==  7.0;
  val = rsConsistentUnwrappedValue0( 8.0, 2.0, rangeMax); r &= val ==  7.0;
  val = rsConsistentUnwrappedValue0( 9.0, 2.0, rangeMax); r &= val ==  7.0;

  val = rsConsistentUnwrappedValue0(10.0, 2.0, rangeMax); r &= val == 12.0;
  val = rsConsistentUnwrappedValue0(11.0, 2.0, rangeMax); r &= val == 12.0;
  val = rsConsistentUnwrappedValue0(12.0, 2.0, rangeMax); r &= val == 12.0;
  val = rsConsistentUnwrappedValue0(13.0, 2.0, rangeMax); r &= val == 12.0;
  val = rsConsistentUnwrappedValue0(14.0, 2.0, rangeMax); r &= val == 12.0;

  // maybe try with a range of not exactly representable numbers such as -pi..pi

  // tests various functions that have to do with phase-unwrapping

  // todo: test wrapped interpolation, rsWarpTointerval, unwrap, ...


  return r;
}

// s: stencil offsets, d: derivative order, t: target coeffs, tol: tolerance
bool testStencil(const std::vector<double>& s, int d, const std::vector<double>& t, 
  double tol = 1.e-13)
{
  int N = (int) s.size();
  std::vector<double> c(N);
  RAPT::rsNumericDifferentiator<double>::stencilCoeffs(&s[0], N, d, &c[0]); 
  return rsEquals(c, t, tol);
}
bool testNumDiffStencils()
{
  // Tests the computation of coefficients for arbitrary finite difference stencils according to:
  // http://web.media.mit.edu/~crtaylor/calculator.html. We compare the results returned by
  // RAPT::getNumDiffStencilCoeffs to the results produced by that website. There's a copy of that 
  // html file in my private repo, just in case, the page disappears - the html has the javascript 
  // code embedded

  bool r = true;
  typedef std::vector<double> Vec;
  Vec s;

  // symmetric, equidistant 3-point stencil -1,0,1:
  s = {-1, 0, 1};
  r &= testStencil(s, 1, Vec({-1.,   0.,  1}) / 2.);
  r &= testStencil(s, 2, Vec({ 1.,  -2.,  1}) / 1.);

  // symmetric, equidistant 5-point stencil -2,-1,0,+1,+2:
  s = {-2, -1, 0, 1, 2};
  r &= testStencil(s, 1, Vec({ 1., -8.,   0.,  8., -1.}) / 12.);
  r &= testStencil(s, 2, Vec({-1., 16., -30., 16., -1.}) / 12.);
  r &= testStencil(s, 3, Vec({-1.,  2.,   0., -2.,  1.}) /  2.);
  r &= testStencil(s, 4, Vec({ 1., -4.,   6., -4.,  1.}) /  1.);

  // symmetric, equidistant 7-point stencil -3,-2,-1,0,+1,+2,+3:
  s = {-3, -2, -1, 0, 1, 2, 3};
  r &= testStencil(s, 1, Vec({-1.,   9., -45.,    0.,   45.,  -9.,  1.}) /  60.);
  r &= testStencil(s, 2, Vec({ 2., -27., 270., -490.,  270., -27.,  2.}) / 180.);
  r &= testStencil(s, 3, Vec({ 1.,  -8.,  13.,    0.,  -13.,   8., -1.}) /   8.);
  r &= testStencil(s, 4, Vec({-1.,  12., -39.,    56., -39.,  12., -1.}) /   6.);
  r &= testStencil(s, 5, Vec({-1.,   4.,  -5.,     0.,   5.,  -4.,  1.}) /   2.);
  r &= testStencil(s, 6, Vec({ 1.,  -6.,  15.,   -20.,  15.,  -6.,  1.}) /   1.);

  return r;

  // todo: try some weird stencils (asymmetric and/or non-equidistant, even  number of points,...)
}

// Compute Hessian using NumDiff::hessianTimesVector - for each unit vector, we compute one row
// of the Hessian. ...maybe move to ScratchPad or Prototypes
template<class T, class F>
static void hessian2(const F& f, T* x, int N, T* pH, const T* h, T k)
{
  using Vec = std::vector<double>;
  using NumDiff = rsNumericDifferentiator<double>;
  Vec wrk(2*N);
  Vec e(N);  // should we init this with zeros or can we rely on being initialized?
  for(int n = 0; n < N; n++)  {
    e[n] = 1;
    NumDiff::hessianTimesVector(f, x, &e[0], N, &pH[N*n], h, k, &wrk[0]);
    e[n] = 0; }
  // # evals: N * 4*N = 4*N^2
}

bool testNumericGradientAndHessian()
{
  // Tests numeric gradient and Hessian matrix computation by comparing the results to analytically
  // computed ones.

  bool r = true;

  // Our example is the trivariate scalar function:
  //   f(x,y,z) = x^2 * y^3 * z^4
  // where v = (x y z) is the 3D input vector.
  std::function<double(double*)> f = [=](double* v)->double
  { 
    double x = v[0], y = v[1], z = v[2];
    return x*x * y*y*y * z*z*z*z; 
  };
  // We need to wrap f into a std::function object and cannot use a raw lambda function like:
  //   auto f = [=](double* v)->double  ...
  // because the function NumDiff::hessian requires a std::function. This is because it's defined
  // in the cpp file and therefore not inlined, so we need an explicit instantiation - and i don't 
  // know how to make explicit instantiations for lambda functions (or if that is even possible) 
  // and that would not be desirable anyway because it would presumbly need a separate 
  // instantiation for each lambda function that is defined somewhere. I've found somewhere the 
  // statement, that lambdas can decay into function pointers, iff they have no capture variables, 
  // so maybe in this case, and instantiation for a C-style function-pointer would work. However, 
  // std::function seems to be the most convenient and flexible way to do it, so that's what i 
  // opted for for instantiation.

  // Computes the gradient of f:
  //   g(f) = (f_x  f_y  f_z)
  // at the given position v analytically and writes the result into g:
  auto gf = [=](double* v, double* g)
  {
    double x = v[0], y = v[1], z = v[2];
    g[0] = 2*x * y*y*y * z*z*z*z;  // f_x
    g[1] = x*x * 3*y*y * z*z*z*z;  // f_y
    g[2] = x*x * y*y*y * 4*z*z*z;  // f_z
  };

  // Computes the Hessian matrix of f:
  //          f_xx  f_xy  f_xz
  //   H(f) = f_yx  f_yy  f_yz
  //          f_zx  f_zy  f_zz
  // at the given position v analytically.
  auto Hf = [=](double* v)->rsMatrix<double>
  {
    double x = v[0], y = v[1], z = v[2];
    rsMatrix<double> H(3,3);
    H(0,0)          = 2*1 * y*y*y * z*z*z*z;  // f_xx
    H(0,1) = H(1,0) = 2*x * 3*y*y * z*z*z*z;  // f_xy = f_yx
    H(0,2) = H(2,0) = 2*x * y*y*y * 4*z*z*z;  // f_xz = f_zx
    H(1,1)          = x*x * 3*2*y * z*z*z*z;  // f_yy
    H(1,2) = H(2,1) = x*x * 3*y*y * 4*z*z*z;  // f_yz = f_zy
    H(2,2)          = x*x * y*y*y * 4*3*z*z;  // f_zz
    return H;
  };


  using Vec = std::vector<double>;
  using Mat = rsMatrix<double>;
  using NumDiff = rsNumericDifferentiator<double>;
  //using ND  = rsNumericDifferentiator<double, double>;

  double hx = pow(2, -18);    // from 2^-18, maxErr in the gradient becomes 0
  double hy = hx/2;
  double hz = hx/4;
  double h[3] = {hx, hy, hz};  // h as array

  Vec v({5,3,2});         // point at which we evaluate gradient and Hessian - maybe rename to x
  double vf = f(&v[0]);   // compute function value at v

  // compute gradient analytically and numerically and compare results:
  Vec ga(3); gf(&v[0], &ga[0]);
  Vec gn = NumDiff::gradient(f, v, h);
  Vec err = ga - gn;
  double maxErr = rsMaxAbs(err);
  r &= maxErr == 0.0;

  // compute Hessian matrix analytically and numerically and compare results:
  Mat Ha = Hf(&v[0]);
  Mat Hn = NumDiff::hessian(f, v, h);
  Mat He = Ha - Hn;  // error matrix
  maxErr = He.getAbsoluteMaximum();
  r &= maxErr == 0.0;

  // compute 1st and 2nd partial derivatives:
  double f0,f1,f2;
  NumDiff::partialDerivativesUpTo2(f, &v[0], 3, 0, h[0], &f0, &f1, &f2); 
  r &= f0 == vf && f1 == ga[0] && f2 == Ha(0, 0);
  NumDiff::partialDerivativesUpTo2(f, &v[0], 3, 1, h[1], &f0, &f1, &f2);
  r &= f0 == vf && f1 == ga[1] && f2 == Ha(1, 1);
  NumDiff::partialDerivativesUpTo2(f, &v[0], 3, 2, h[2], &f0, &f1, &f2);
  r &= f0 == vf && f1 == ga[2] && f2 == Ha(2, 2);

  // compute Hessian times a vector w:
  Vec w({-2,1,3});    // fixed vector
  Vec Hwa = Ha * w;   // analytic matrix vector product H*w
  Hwa = w * Ha;
  double k = pow(2, -13);
  Vec Hwn(3);         // numeric matrix vector product H*w
  NumDiff::hessianTimesVector(f, &v[0], &w[0], 3, &Hwn[0], h, k);
  maxErr = rsMaxAbs(Hwa-Hwn);
  r &= maxErr < 0.004;  // it doesn't get any better by tweaking k - but maybe we could tweak
                        // the h-values as well...


  // not a unit-test test - just trying to figure out what we get, when computing the Hessian times 
  // the coordinate unit vectors:
  w = Vec({1,0,0}); NumDiff::hessianTimesVector(f, &v[0], &w[0], 3, &Hwn[0], h, k);
  w = Vec({0,1,0}); NumDiff::hessianTimesVector(f, &v[0], &w[0], 3, &Hwn[0], h, k);
  w = Vec({0,0,1}); NumDiff::hessianTimesVector(f, &v[0], &w[0], 3, &Hwn[0], h, k);
  // ...when w is the i-th unit vector, we get the i-th row (or column) of the Hessian, so
  // maybe another way to numerically approximate the Hessian is to use hessianTimesVector with
  // the N coordinate basis vectors? maybe the H(i,j) element of the Hessian says how much the i-th
  // component of the gradient changes when wiggling the j-th coordinate? or vice versa? it 
  // shouldn't make a difference because the Hessian is symmetric - but why would the j-th gradient
  // element g[i] change by the same same amount when wiggling x[i] as the g[j] when wiggling x[j]
  // ...intutively, these seem to be unrelated quantities -> figure out!
  // ....and implement another algorithm to estimate the Hessian baseed on the above observation 
  // and compare it to the algo we already have in terms off accuracy and number of function 
  // evaluations

  hessian2(f, &v[0], 3, Hn.getDataPointer(), h, k);
  // The 6th value changes a bit (i.e. gets less accurate), so this algo seems to be less accurate 
  // than NumDiff::hessian (but maybe that can be improved by using a different k). It also uses 
  // 4*N^2 evaluations of f compared to 2*N^2 + 1 of NumDiff::hessian, and it needs a workspace 
  // which NumDiff::hessian doesn't. It also needs the user to specify this additional step-size k. 
  // So, this way of computing the Hessian is clearly inferior in all sorts of ways, so there's no 
  // point in adding it to the library. Maybe move this code elsewhere - it doesn't really belong 
  // in a unit test - but here we already have all the variables for the call available...

  return r;
}


template<class T, class F>
static void jacobian(const F& f, T* x, int N, T* pJ, const T* h)
{
  using NumDif = rsNumericDifferentiator<T>;
  int M = (int) f.size();  // number of outputs, maybe use rsSize
  for(int m = 0; m < M; m++)
    NumDif::gradient(f[m], x, N, &pJ[m*N], h);
}
// is this the actual jacobian or the transposed jacobian?

bool testNumericJacobian()
{
  bool r = true;

  // create a function from R^2 to R^3, i.e. a parametric surface

  using Func = std::function<double(double*)>;
  using Mat   = rsMatrix<double>;

  int N = 2;  // number of input dimensions
  int M = 3;

  // input vector:
  double x = 7, y = 5;
  std::vector<double> v({x,y});

  // approximation stepsizes:
  double hx = pow(2, -18);    // from 2^-18, maxErr in the gradient becomes 0
  double hy = hx/2;
  double hz = hx/4;
  double h[3] = {hx, hy, hz};  // h as array

  // Our example functions are:
  //   f1(x,y) =   x^2 *   y^3
  //   f2(x,y) = 2*x^3 * 3*y^2
  //   f3(x,y) = x^2 + y^2 - 2*x*y
  // and the Jacobian is:
  //   df1/dx = 2*x   *   y^3      df1/dy = x^2   * 3*y^2
  //   df2/dx = 6*x^2 * 3*y^2      df2/dy = 2*x^3 * 6*y
  //   df3/dx = 2*x * y^2 - 2*y    df3/dy = x^2 * 2*y - 2*x
  Func f1 = [=](double* v)->double { double x = v[0], y = v[1]; return x*x     * y*y*y;   };
  Func f2 = [=](double* v)->double { double x = v[0], y = v[1]; return 2*x*x*x * 3*y*y;   };
  Func f3 = [=](double* v)->double { double x = v[0], y = v[1]; return x*x * y*y - 2*x*y; };
  std::vector<Func> f({f1,f2,f3});  // array/vector of functors
  Mat Ja(3, 2);                     // analytical Jacobian
  Ja(0, 0) = 2*x * y*y*y; 
  Ja(0, 1) = x*x * 3*y*y;
  Ja(1, 0) = 6*x*x * 3*y*y;
  Ja(1, 1) = 2*x*x*x * 6*y;
  Ja(2, 0) = 2*x * y*y - 2*y;
  Ja(2, 1) = x*x * 2*y - 2*x;
  // what if we use auto instead of Func for f1,f2,f3?

  // Compute numerical Jacobian and compare with analytic result:
  Mat Jn(3, 2);  
  double *pJ = Jn.getDataPointer();
  jacobian(f, &v[0], N, pJ, h);
  Mat Je = Ja - Jn;  // error matrix
  double maxErr = Je.getAbsoluteMaximum();
  r &= maxErr == 0.0;

  return r;
}

bool testMultiLayerPerceptronOld(std::string &reportString)
{
  std::string testName = "MultiLayerPerceptron";
  bool testResult = true;

  // this should actually be in the "Experiments" suite

  //int hiddenNeuronsInLayers[2] = {2, 4};
  //MultiLayerPerceptron mlp(3, 2, 2, hiddenNeuronsInLayers);

  //int hiddenNeuronsInLayers[3] = {2, 4, 3};
  //MultiLayerPerceptron mlp(1, 1, 3, hiddenNeuronsInLayers);

  //int hiddenNeuronsInLayers[3] = {4, 3, 5};
  //MultiLayerPerceptron mlp(1, 1, 3, hiddenNeuronsInLayers);

  int hiddenNeuronsInLayers[1] = {10};
  MultiLayerPerceptron<double> mlp(1, 1, 1, hiddenNeuronsInLayers);

  //int hiddenNeuronsInLayers[2] = {3, 5};
  //MultiLayerPerceptron mlp(1, 1, 2, hiddenNeuronsInLayers);

  //int hiddenNeuronsInLayers[2] = {2,4};
  //MultiLayerPerceptron mlp(1, 1, 2, hiddenNeuronsInLayers);

  // create the error-function object and pass some training data to it:
  MultiLayerPerceptronErrorFunction<double> mlpError(&mlp);
  static const int N = 200;
  rsVectorDbl x[N];
  rsVectorDbl y[N];
  for(int n=0; n<N; n++)
  {
    x[n].setDimensionality(1);
    y[n].setDimensionality(1);
  }
  x[0].v[0] = rsRandomUniform(-3.0, 3.0, 7);
  for(int n=0; n<N; n++)
  {
    x[n].v[0] = rsRandomUniform(-3.0, 3.0);
    //x[n].v[0] = (double) n / (double) (N-1);
    //x[n].v[0] = linToLin(x[n].v[0], 0.0, 1.0, -3.0, 3.0);
    y[n].v[0] = x[n].v[0] * x[n].v[0];            // y = x^2
    //y[n].v[0] = fabs(x[n].v[0]);                    // y = |x|
    //y[n].v[0] = sign(x[n].v[0]);                    // y = sign(x)
    //y[n].v[0] = sin(x[n].v[0] * x[n].v[0]);            // y = sin(x^2)
  }
  mlpError.setTrainingData(x, y, N);

  // create the minimizer and minimize the training error:
  GradientBasedMinimizer<double> mlpTrainer;
  mlpTrainer.setPrintInfo(true);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::GRADIENT_DESCENT);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::BOLD_DRIVER_WITH_MOMENTUM);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::CONJUGATE_GRADIENT);
  mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::SCALED_CONJUGATE_GRADIENT);
  //mlpTrainer.setBetaFormula(GradientBasedMinimizer<double>::POLAK_RIBIERE);
  //mlpTrainer.setBetaFormula(GradientBasedMinimizer<double>::FLETCHER_REEVES);
  //mlpTrainer.setBetaFormula(GradientBasedMinimizer<double>::HESTENES_STIEFEL);
  rsVectorDbl w = mlpTrainer.minimizeFunction(&mlpError, mlp.getWeightsAsVector());
  mlp.setWeightVector(w);

  int dummy = 0;


  //double *test = new double[0];
  //test[0] = 1.0;

  //double x[3]  = {3,2,1};
  //double yt[2] = {1,-1};
  //double y[2];
  //mlp.computeNetworkOutput(x, y);
  //mlp.printWeights();
  //mlp.printActivations();

  //mlpError.computePatternGradientByWeightPerturbation(Vector(2, yt));
  //mlpError.printPatternGradient();

  //mlpError.computePatternGradient(Vector(2, yt));
  //mlpError.printPatternGradient();



  //Vector wv = mlp.getWeightsAsVector();
  //wv.print();
  //mlp.setWeightVector(wv);
  //mlp.printWeights();
  //Vector gv = mlpTrainer.getGradientVector();
  //gv.print();


  double xTest[N];
  double yTest[N];
  for(int n=0; n<N; n++)
  {
    xTest[n] = (double) n / (double) (N-1);
    xTest[n] = rsLinToLin(xTest[n], 0.0, 1.0, -6.0, 6.0);
    mlp.computeNetworkOutput(&xTest[n], &yTest[n]);
    rsVectorDbl yVec = mlp.getOutput();
    int dummy = 0;
  }
  //Plotter::plotData(N, xTest, yTest);
    
  printf("%s", "\n Press Key \n");
  getchar();

  return testResult;
}


rsVectorDbl mlpTestFunction(rsVectorDbl x, double noiseAmplitude)
{
  rsAssert(x.dim == 3);
  rsVectorDbl y(2);

  // compute output values according to deterministic function:
  y[0]  =  2*x[0] + 3*x[1] - 1*x[2] - 2*x[0]*x[1] + 1*x[0]*x[2] + 3*x[1]*x[2];
  y[1]  = -3*x[0] + 1*x[1] + 2*x[2] + 3*x[0]*x[1] - 2*x[0]*x[2] + 1*x[1]*x[2];

  // add some noise:
  y[0] += noiseAmplitude * rsRandomUniform(-1.0, 1.0);
  y[1] += noiseAmplitude * rsRandomUniform(-1.0, 1.0);

  return y;
}
rsVectorDbl randomVector(int numDimensions, double min, double max)
{
  rsVectorDbl x(numDimensions);
  for(int i = 0; i < numDimensions; i++)
    x[i] = rsRandomUniform(min, max);
  return x;
}
bool testMultiLayerPerceptron(std::string &reportString)
{
  std::string testName = "MultiLayerPerceptron";
  bool testResult = true;

  // we try to approximate the function f: R^3 -> R^2, defined by:
  // y1 =  2*x1 + 3*x2 - 1*x3 - 2*x1*x2 + 1*x1*x3 + 3*x2*x3
  // y2 = -3*x1 + 1*x2 + 2*x3 + 3*x1*x2 - 2*x1*x3 + 1*x2*x3
  // by using noisy datapoints

  // ...that particular function seems to be quite difficult to approximate by an MLP
  // maybe use a function that is easier to represent - possibly by using outputs of an MLP
  // itself (with randomized weights) - see if the weights can be retrieved by training (or an 
  // equivalent weight-vector - there are weight-space symmetries)

  // create test data:
  int n;
  static const int numPatterns    = 200;
  double noiseAmplitude = 0.2;
  rsVectorDbl x[numPatterns];
  rsVectorDbl y[numPatterns];
  for(n = 0; n < numPatterns; n++)
  {
    x[n].setDimensionality(3);
    y[n].setDimensionality(2);
  }
  rsRandomUniform(-1.0, 1.0, 7); // init PRNG
  for(n = 0; n < numPatterns; n++)
  {
    x[n] = randomVector(3, -2.0, 2.0);
    y[n] = mlpTestFunction(x[n], noiseAmplitude);
  }

  // create an MLP with 3 hidden layers with 4, 3 and 5 neurons respectively:
  //int hiddenNeuronsInLayers[3] = {4, 3, 5};
  //MultiLayerPerceptron mlp(3, 2, 3, hiddenNeuronsInLayers);

  int hiddenNeuronsInLayers[1] = {10};
  MultiLayerPerceptron<double> mlp(3, 2, 1, hiddenNeuronsInLayers);

  // create the error-function object and pass the training data to it:
  MultiLayerPerceptronErrorFunction<double> mlpError(&mlp);
  mlpError.setTrainingData(x, y, numPatterns);

  // create the minimizer for the error-function find the weights that minimize the training error:
  GradientBasedMinimizer<double> mlpTrainer;
  mlpTrainer.setPrintInfo(true);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::GRADIENT_DESCENT);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::BOLD_DRIVER_WITH_MOMENTUM);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::CONJUGATE_GRADIENT);
  mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::SCALED_CONJUGATE_GRADIENT);
  rsVectorDbl w = mlpTrainer.minimizeFunction(&mlpError, mlp.getWeightsAsVector());
  mlp.setWeightVector(w);  // set up the network with the optimal weight-vector

  return testResult;
}


bool testMiscMath()
{
  std::string dummy;    // get rid

  bool testResult = true;

  testResult &= testExponentialCurveFitting(  dummy);
  testResult &= testRootFinding(              dummy);
  testResult &= testGradientBasedOptimization(dummy);
  testResult &= testMinSqrDifFixSum(          dummy);
  testResult &= testPhaseUnwrapStuff(         dummy);
  testResult &= testNumDiffStencils();
  testResult &= testNumericGradientAndHessian();
  testResult &= testNumericJacobian();
  //testResult &= testNumericMinimization();  // has been moved to experiments

  //testResult &= testMultiLayerPerceptronOld(  dummy); // produces verbose output
  //testResult &= testMultiLayerPerceptron(     dummy); // maybe move to experiments

  return testResult;
}