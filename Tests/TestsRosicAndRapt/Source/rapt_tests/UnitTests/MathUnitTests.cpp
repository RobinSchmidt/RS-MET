
// these files are not all compiled all by themselves separately in order to reduce the number of
// compilation unit to improve build time:
#include "Math/LinearAlgebraUnitTests.cpp" // produces linker errors
#include "Math/PolynomialUnitTests.cpp"
#include "Math/VectorUnitTests.cpp"
#include "Math/MatrixUnitTests.cpp"
#include "Math/MiscMathUnitTests.cpp"

// under construction:
#include "Math/DifferentialEquationTests.cpp"
#include "Math/IntegerFunctionTests.cpp"
#include "Math/MoebiusTransformTests.cpp"
#include "Math/NumberTheoryTests.cpp"
#include "Math/RealFunctionTests.cpp"
#include "Math/TransformsTests.cpp"
//#include "Math/StatisticsTests.cpp"  // there is no such file - why?
#include "Math/Point2DTests.cpp"
#include "Math/Polygon2DTests.cpp"
#include "Math/Triangle2DTests.cpp"

#include "Math/MultiArrayTests.cpp"

bool coordinateMapperUnitTest()
{
  bool r = true;      // test result

  rsCoordinateMapper2D<float> mapper;

  mapper.setInputRange( 0.125f, 8.f, -50.f, +50.f);
  mapper.setOutputRange(0.f,  400.f, 200.f,   0.f);
  mapper.mapperX.setLogScaled(true);

  float x, y;
  x = mapper.mapX(1.f);      r &= x == 200.f;
  x = mapper.unmapX(200.f);  r &= x == 1.f;
  y = mapper.mapY(0.f);      r &= y == 100.f;
  y = mapper.unmapY(100.f);  r &= y == 0.f;

  return r;
}

bool correlationUnitTest()
{
  bool ok = true;              // test result

  static const int N = 11;    // signal length (todo: make this a variable in a loop to test
                              // with different lengths)
  static const int M = 2*N-1; // length of correlation sequence

  using AT = rsArrayTools;

  double x[N], y[N], xr[N], yr[N];    // inputs and reversed versions
  double c1[M], c2[M], c3[M], c4[M];  // correlation sequences (results)
  //AT::fillWithRandomValues(x, N, -1, +1, 0);
  //AT::fillWithRandomValues(y, N, -1, +1, 1);
  AT::fillWithRandomIntegers(x, N, -9, +9, 0);
  AT::fillWithRandomIntegers(y, N, -9, +9, 1);
  AT::reverse(x, xr, N);
  AT::reverse(y, yr, N);

  // Obtain cross-correlation sequences via various algorithms:
  rsCrossCorrelationDirect(x, y, N, c1);
  rsCrossCorrelationDirect(y, x, N, c2);  // should it be commutative? it's not!
  rsCrossCorrelationFFT(   x, y, N, c3);
  AT::convolve(x, N, yr, N, c4);

  // Results should all be the same up to roundoff:
  double err;
  double tol = 1.e-13;
  //err = AT::maxDeviation(c1, c2, N); ok &= err <= tol; // nope!
  err = AT::maxDeviation(c1, c3, N); ok &= err <= tol;
  //err = AT::maxDeviation(c1, c4, N); ok &= err <= tol;
  // c4 is completely different - wtf? ...oh - it seems, the rsCrossCorrelation functions only
  // return the 2nd half of the array - well, yeah, the results are only of length N whereas
  // rsArrayTools::convolve produces a result of length 2*N-1. c4 is also completely different
  // except for the first coeff

  // Test autocorrelation:
  rsCrossCorrelationDirect(x, x, N, c1);
  rsAutoCorrelationFFT(    x,    N, c2);
  err = AT::maxDeviation(c1, c2, N); ok &= err <= tol;




  // try de-biasing a convolution result:
  AT::fillWithValue(x, N, 1.0);
  AT::fillWithValue(y, N, 1.0);
  AT::convolve(x, N, y, N, c3);
  for(int n = 0; n < N; n++)
  {
    double scale = double(N)/double(N-n);
    c3[N-1+n] *= scale;
    c3[N-1-n] *= scale;
    int dummy = 0;
  }
  // ok - looks good


  return ok;
}

// move to RAPT - CurveFitting ...hmm...but actually it's not a curve-fit in the least-squares
// sense but an interpolation...maybe it belongs into a class rsInterpolator
void fitOddRatFunc4(double *x, double* y, double *coeffs)
{
  // establish matrix and rhs vector:
  double** A;
  double b[4];
  RAPT::rsArrayTools::allocateSquareArray2D(A, 4);
  for(int i = 0; i < 4; i++) {
    double xi2 = x[i]*x[i];     //  xi^2
    A[i][0] = x[i];             //  xi
    A[i][1] = xi2*x[i];         //  xi^3
    A[i][2] = -xi2*y[i];        // -xi^2*yi
    A[i][3] = -xi2*xi2*y[i];    // -xi^4*yi
    b[i]    = y[i];
  }

  // solve linear system and clean up:
  RAPT::rsLinearAlgebra::rsSolveLinearSystemInPlace(A, coeffs, b, 4);
  RAPT::rsArrayTools::deAllocateSquareArray2D(A, 4);
}
double oddRatFunc4(double x, double *c) // evaluates (c0*x + c1*x^3) / (1 + c2*x^2 + c3*x^4)
{
  double x2 = x*x;  // x^2
  return (c[0]*x + c[1]*x*x2) / (1 + c[2]*x2 + c[3]*x2*x2);
}
bool fitRationalUnitTest()
{
  bool r = true;      // test result

  // create data to match and allocate coeffs:
  double x[4] = { 0.25, 0.5,  0.75,  1.0 };  // x-values
  double y[4] = { 0.5,  0.85, 0.975, 1.0 };  // y-values
  double c[4]; fitOddRatFunc4(x, y, c);      // coeffs

  // evaluate function at the data point and verify, if it matches them indeed:
  double Y[4];
  double tol = 1.e-13;
  for(int i = 0; i < 4; i++) {
    Y[i] = oddRatFunc4(x[i], c);
    r &= fabs(y[i]-Y[i]) < tol;
  }

  // the function values match indeed - but the function may be useless as an approximation to the
  // data as it features a pole within the range of data points, see the plot:
  // https://www.desmos.com/calculator/5le2p3qsvb
  // -can this be fixed by using x^5 instead of x^3 in the numerator? or replace x by x^5?
  // -can the fitting routine be suitably generalized to allow the user to specify which
  //  exponents he wants in the numerator and denominator?
  // -can we detect problems by checking if there are poles with in the x-range of datapoints?
  // -can the algorithm be extended to choose itself a set of exponents that works? ..just
  //  systematically try a bunch of combinations of exponents and return the first one that works?
  //  ...or maybe try a few more and return one that works best - ...but how do we define best?

  return r;
}
/*
todo: use te same approach to interpolate with arbitrary basis functions. we have N+M basis
functions, N in the numerator and M in the denominator and use the relation
yi = f(xi) = sum(i=0..N-1, ci*fi(xi)) / ( 1 + sum(i=N..L, ci*fi(xi)) ), where L = N+M-1
this gives the matrix equation:

|f0(x0) f1(x0) ... fN-1(x0) -y0*fN(x0) ... -y0*fL(x0)|   |c0|   |y0|
|                                                    | * |  | = |  |
|f0(xL) f1(xL) ... fN-1(xL) -yL*fN(xL) ... -yL*fL(xL)|   |cL|   |yL|

which can solved just the same way. When M=0 and fi = x^i, we interpolate with an N-1 th order
polynomial - test with the results are the same as with regular polynomial interpolation. The
fitting routine should get as inputs the arrays of x and y values and an array of function-pointers
for the basis functions (maybe std::function) and give back the coeffs as usual */

bool interpolatingFunctionUnitTest()
{
  bool r = true;      // test result

  rsNodeBasedFunction<float> ipf; // rename to nbf

  std::vector<float> xa, ya; // x,y value arrays
  float /*x,*/ y;                // x,y values
  size_t i;                  // in dex

  y = ipf.getValue(10.f);
  r &= y == 0.f;

  // test adding datapoints:

  // i: 0
  // x: 0
  // y: 2
  i = ipf.addNode(0.f, 2.f); r &= i == 0;
  y = ipf.getValue(-1.f);    r &= y == 2.f;
  y = ipf.getValue( 0.f);    r &= y == 2.f;
  y = ipf.getValue(+1.f);    r &= y == 2.f;

  // i: 0  1
  // x: 0  1
  // y: 2  3
  i = ipf.addNode(1.f, 3.f); r &= i == 1;
  y = ipf.getValue( 0.25f);  r &= y == 2.25f;

  // i: 0  1  2
  // x: 0  1  5
  // y: 2  3  1
  i = ipf.addNode(5.f, 1.f); r &= i == 2;
  y = ipf.getValue( 3.f);    r &= y == 2.f;

  // i: 0  1  2  3
  // x: 0  1  4  5
  // y: 2  3  0  1
  i = ipf.addNode(4.f, 0.f); r &= i == 2;
  y = ipf.getValue( 2.f);    r &= y == 2.f;
  y = ipf.getValue( 3.f);    r &= y == 1.f;

  // i: 0  1  2  3  4
  // x: 0  1  3  4  5
  // y: 2  3  2  0  1
  i = ipf.addNode(3.f, 2.f); r &= i == 2;
  y = ipf.getValue( 2.f);    r &= y == 2.5f;

  // i: 0  1  2  3  4  5
  // x: 0  1  2  3  4  5
  // y: 2  3  4  2  0  1
  i = ipf.addNode(2.f, 4.f); r &= i == 2;
  y = ipf.getValue( 1.25f);  r &= y == 3.25f;

  // i: 0  1  2  3  4  5  6
  // x: 0  1  2  2  3  4  5
  // y: 2  3  4  6  2  0  1
  i = ipf.addNode(2.f, 6.f); r &= i == 3;
  y = ipf.getValue( 2.5f);   r &= y == 4.f;

  // i:  0  1  2  3  4  5  6  7
  // x: -1  0  1  2  2  3  4  5
  // y:  0  2  3  4  6  2  0  1
  i = ipf.addNode(-1.f, 0.f); r &= i == 0;
  y = ipf.getValue( -.5f);    r &= y == 1.f;

  // test moving datapoints:

  // i:  0  1  2  3  4  5  6  7
  // x: -1  0  1  2  3  3  4  5
  // y:  0  2  3  4  2  6  0  1
  i = ipf.moveNode(4, 3.f, 6.f); r &= i == 5;


  // test removing datapoints:

  // i:  0  1  2  3  4  5  6
  // x: -1  0  1  2  3  4  5
  // y:  0  2  3  4  6  0  1
  ipf.removeNode(4);
  y = ipf.getValue( 2.5f); r &= y == 5.f;

  return r;
}

bool resampleNonUniform()
{
  bool r = true;

  static const int inLength = 5;
  double xIn[] = { 2, 4, 5, 7, 8 };
  double yIn[] = { 1, 3, 1, 2, 3 };

  static const int outLength = 9;
  double xOut[] = { 1.5, 2.5, 3, 3.5, 4.5, 6, 7.5, 8.5, 9.5 };
  double yOut[outLength];

  resampleNonUniformLinear(xIn, yIn, inLength, xOut, yOut, outLength);

  r &= yOut[0] == 0.5;
  r &= yOut[1] == 1.5;
  r &= yOut[2] == 2.0;
  r &= yOut[3] == 2.5;
  r &= yOut[4] == 2.0;
  r &= yOut[5] == 1.5;
  r &= yOut[6] == 2.5;
  r &= yOut[7] == 3.5;
  r &= yOut[8] == 4.5;

  return r;
}

template<class T>
bool splineSlopesUnitTest(T tol)
{
  // Example taken from Meister - Numerik, page 62. It's the Runge function 1/(1+x^2). Maybe 
  // evaluate it at more points

  bool ok = true;

  using Vec = std::vector<T>;
  int N = 5;

  /*
  Vec x({ -1,0,1,2,3 }), y(N); 
  for(int i = 0; i < N; i++)
    y[i] = T(1) / (1 + x[i]*x[i]);// ({ 0.5,1,0.5,0.2,0.1 });

  Vec t({39, -3, -27, -9, 3}); t = t / T(50);   // target values for slope
  Vec s = splineSlopes(x, y, true, T(0), T(0)); // natural spline
  Vec r = s - t;                                // error
  ok &= rsIsCloseTo(s, t, tol);

  t = Vec({175, 6, -199, -50, -21}) / T(350);
  s = splineSlopes(x, y, false, T(1)/2, -T(6)/100);
  r = s - t;
  ok &= rsIsCloseTo(s, t, tol);
  */

  Vec x, y, t, s, r;

  // We create some unequally spaced data using a cubic polynomial. The cubic spline should exactly
  // match the original polynomial at all x, iff we also prescribe the correct values for the (1st 
  // or 2nd) derivative at the boundary points. In this case, the computed slopes should of course 
  // also match the actual slopes of the polynomial at the sampling points:
  x = Vec({ -5,-3,-2,-1,2,3,5,6,8 });
  N = (int) x.size();;
  y.resize(N);
  t.resize(N);
  rsPolynomial<T> p({ 7,5,3,2 });              // our polynomial is y(x) = 7 + 5*x + 3*x^2 + 2*x^3
  for(int i = 0; i < N; i++) {
    y[i] = p.evaluate(x[i]);                   // evaluate polynomial at the x[i]
    t[i] = p.derivativeAt(x[i]); }             // target values for the computed slopes
  s = splineSlopes(x, y, false, t[0], t[N-1]); // prescribe correct slopes at boundary points
  ok &= rsIsCloseTo(s, t, tol);                // computed slopes should match target slopes
  // works only when T=double, apparently, we get integer overflow with fractions

  // Now, we prescribe 2nd derivatives at the endpoints instead of 1st derivatives:
  T k0 = p.derivativeAt(x[0],   2);
  T k1 = p.derivativeAt(x[N-1], 2);
  s = splineSlopes(x, y, true, k0, k1);
  r = s - t;
  ok &= rsIsCloseTo(s, t, tol);

  //rsPlotVector(r);
  // error starts small at left and grows big at right end...maybe the equation for the last 
  // datapoint is wrong

  // ToDo: 
  // -Test it using data obtained from a cubic. In this case, the computed spline derivatives 
  //  should exactly match those of the original polynomial, regardless where the datapoints are.
  //  ...but only if we completely specify the spline, i.e. prescribe either 1st or 2nd derivative
  //  values at the end which should match those of the original polynomial

  return ok;
}

bool interpolationUnitTest()
{
  bool ok = true;

  ok &= fitRationalUnitTest();  // fails on linux ("illegal instruction") - encounters singular matrix
  ok &= interpolatingFunctionUnitTest();
  ok &= resampleNonUniform();
  //ok &= splineSlopesUnitTest<rsFraction<int>>(0);
  ok &= splineSlopesUnitTest<double>(1.e-13);

  return ok;
};



// For testing the root-finder, we use a 3rd oder polynomial as example function with roots at
// -1, +1, +2. We define a function and a functor that implements that function in order to pass it
// to the root finder:
float testFunction(float x)
{
  return (x+1)*(x-1)*(x-2);
}
class TestFunctor
{
public:
  inline float operator()(const float x)
  {
    return (x+1)*(x-1)*(x-2); // the roots could be adjustable member data
  }
};
bool testRootFinding(std::function<float(float)>& func, float xL, float xR, float targetRoot,
  float targetY = 0.f)
{
  bool result = true;
  float root;
  float tol = std::numeric_limits<float>::epsilon();

  root = rsRootFinder<float>::bisection(func, xL, xR);
  result &= abs(root-targetRoot) <= tol*targetRoot;   // should it be < or <=?

  root = rsRootFinder<float>::falsePosition(func, xL, xR);
  result &= abs(root-targetRoot) <= tol*targetRoot;

  return result;
}
// todo: templatize this function and test with double and float, maybe include also a tolerance
// parameter - maybe the root-finder itself also needs a tolerance parameter...or maybe two - one
// for dx and one for dy...but maybe that should be optional for convenience

bool rootFinderUnitTest()
{
  bool r = true;                 // test result
  float x, y;                    // function in/out values
  std::function<float(float)> f; // the function to find the roots of
  //float root;



  // moved temporarily up for debugging:
  // that's the evil function from here https://github.com/RobinSchmidt/RS-MET/issues/249
  float b = 5;
  float a = 1.f/tanh(b);
  f = [=] (float x)->float { return ((atanh(x/a)/b)+1.f)*0.5f; };
  //r &= testRootFinding(f, 0.f, 1.f, 1.f, 1.f);

  f = [] (float x)->float { return 1.f/(x-1); };        // has a pole at x=1
  //r &= testRootFinding(f, -1.f, +2.f, 1.f, 0.f);  // biscetion works, false position fails

  f = [] (float x)->float { return 1.f/x; };        // has a pole at x=0
  //r &= testRootFinding(f, -1.f, +1.f, 0.f, 0.f);  // this fails

  // make sure to test special cases where the root is exactly at one of the bracket values
  // for example y = x^2 - 1 has roots at -1 and +1
  // make sure it works also when the bracketing interval does not actually bracket the root





  // create example (lambda) function with roots at -1,+1,+2 and verify positions of roots:
  f = [] (float x)->float { return (x+1)*(x-1)*(x-2); };
  y = f(-1.f); r &= y == 0.f;
  y = f( 1.f); r &= y == 0.f;
  y = f( 2.f); r &= y == 0.f;

  // find the roots via bisection:
  x = rsRootFinder<float>::bisection(f, -1.3f, -0.8f); r &= x == -1.f;
  x = rsRootFinder<float>::bisection(f,  0.8f,  1.3f); r &= x ==  1.f;
  x = rsRootFinder<float>::bisection(f,  1.7f,  2.2f); r &= x ==  2.f;

  // find the roots via false position:
  x = rsRootFinder<float>::falsePosition(f, -1.3f, -0.8f); r &= x == -1.f;
  x = rsRootFinder<float>::falsePosition(f,  0.8f,  1.3f); r &= x ==  1.f;
  x = rsRootFinder<float>::falsePosition(f,  1.7f,  2.2f); r &= x ==  2.f;

  // use a function pointer:
  f = &testFunction;
  x = rsRootFinder<float>::falsePosition(f, -1.3f, -0.8f); r &= x == -1.f;
  x = rsRootFinder<float>::falsePosition(f,  0.8f,  1.3f); r &= x ==  1.f;
  x = rsRootFinder<float>::falsePosition(f,  1.7f,  2.2f); r &= x ==  2.f;

  // use a function object (functor):
  //TestFunctor functor; f = functor;
  f = TestFunctor();  // shorter syntax
  x = rsRootFinder<float>::falsePosition(f, -1.3f, -0.8f); r &= x == -1.f;
  x = rsRootFinder<float>::falsePosition(f,  0.8f,  1.3f); r &= x ==  1.f;
  x = rsRootFinder<float>::falsePosition(f,  1.7f,  2.2f); r &= x ==  2.f;

  // OK - this was to test that it works calling it with different types of paremeters for the
  // function (funtion-pointer, lambda function, functor). Now we want to test the robustness if
  // the actual algorithm by throwing some evil functions at it (we use lambda functions for this)

  f = [] (float x)->float { return sin(x); };       // nice function - has a root at x=0
  r &= testRootFinding(f, -1.f, +1.f, 0.f, 0.f);


  // it should probably return 0 as result - the location of jump from -inf to inf

  //f = [] (float x)->float { return (x*x)/x; };     // is undefined at x=0...but has a root there

  // can we somehow also check the numbers of iterations taken?

  return r;
}

//// cube root
//float cbrt(float x)
//{
//  return pow(x, float(1.0/3.0));
//}

// finds one solution to the equation z^3 = p*z + q
float cubicRootPQ(float p, float q)
{
  float q2 = q/2;
  float p3 = p/3;
  float r  = sqrt(q2*q2 - p3*p3*p3);
  return cbrt(q2+r) + cbrt(q2-r);
}

// finds one solution to the equation a * z^3 + b * z^2 + c*z + d = 0
float cubicRoot(float d, float c, float b, float a)
{
  float ai = 1/a; b *= ai; c *= ai; d *= ai;  // make monic (normalize a to 1)
  float b2 = b*b;
  return cubicRootPQ((b2-3*c)*(1.f/3.f), -b2*b*(2.f/27.f) + 0.5f*b*c - d);

  // Formulas for p,q  were found with sage via:
  // var("a b c d z")                # declare symbolic variables
  // f(x) = a*x^3 + b*x^2 + c*x + d  # define a symbolic function
  // g = f.subs(x = z - b/(3*a))     # substitue z-b/(3*a) for x, removes x^2 term
  // h = g/a                         # divide resulting expression by a
  // h = h.collect(z)                # collect terms with respect to z
  // h.coefficients(z)               # gives coefficients for powers of z
}

void cubicRoots(float d, float c, float b, float a, std::complex<float>* r1,
  std::complex<float>* r2, std::complex<float>* r3)
{
  float rr = cubicRoot(d, c, b, a); // find real root
  *r2 = rr;
  float cof[4] = {d, c, b, a };     // collect coeffs into array
  //float rem;                        // dummy for remainder
  //RAPT::rsPolynomial<float>::dividePolynomialByMonomialInPlace(cof, 4, rr, &rem); // causes crash
  //RAPT::rsPolynomial<float>::rootsQuadraticComplex(cof[0], cof[1], cof[2], r1, r3);
  int dummy = 0;
}

bool polynomialRootsUnitTest()
{
  bool r = true;

  typedef RAPT::rsPolynomial<float> P;
  typedef std::complex<float> C;

  float rr1, rr2;                     // real roots
  std::complex<float> cr1, cr2, cr3;  // complex roots
  float d; // discriminant

  // 6 - 9*x + 3*x^2 , roots: 1, 2
  P::rootsQuadraticReal(6.f, -9.f, 3.f, &rr1, &rr2);
  r &= rr1 == 1.f;
  r &= rr2 == 2.f;

  // 15 - 12*x + 3*x^2 , roots: 2-i, 2+i
  P::rootsQuadraticComplex(C(15.f), C(-12.f), C(3.f), &cr1, &cr2);
  r &= cr1 == C(2, -1);
  r &= cr2 == C(2, +1);
  P::rootsQuadraticReal(15.f, -12.f, 3.f, &rr1, &rr2); // should give real parts
  r &= rr1 == 2;
  r &= rr2 == 2;

  // x^3 - x, roots: -1, 0, +1, d = 4
  d = P::cubicDiscriminant(0.f,  -1.f, 0.f, 1.f);
  P::rootsCubicComplex(C(0), C(-1), C(0), C(1), &cr1, &cr2, &cr3);
  r &= d   == 4;
  r &= cr1 == C(-1, 0);
  r &= cr2 == C( 0, 0);
  r &= cr3 == C(+1, 0);
  // has roundoff error, otherwise ok

  // x^3 + x, roots: 0, -i, +i, d = -4
  d = P::cubicDiscriminant( 0.f, 1.f, 0.f, 1.f);
  P::rootsCubicComplex(C(0), C(1), C(0), C(1), &cr1, &cr2, &cr3);
  r &= d   == -4;
  r &= cr1 == C(0,  0);
  r &= cr2 == C(0, -1);
  r &= cr3 == C(0, +1);
  // totally wrong complex roots

  float test;
  test = cubicRootPQ(-1, 0);    // should find root at 0
  test = cubicRoot(0, 1, 0, 1); // this too

  test = cubicRootPQ(2, 4);  // should be 2
  cubicRoots(-8, -4, 0, 2, &cr1, &cr2, &cr3); // roots: 2, -1-i, -1+i


  // -18 + 33*x - 18*x^2 + 3*x^3, roots: 1, 2, 3, d = 324
  d = P::cubicDiscriminant( -18.f, 33.f, -18.f, 3.f);
  P::rootsCubicComplex(C(-18), C(33), C(-18), C(3), &cr1, &cr2, &cr3);
  r &= d   == 324;
  r &= cr1 == C(1, 0);
  r &= cr2 == C(2, 0);
  r &= cr3 == C(3, 0);

  // -6 + 15*x - 12*x^2 + 3*x^3, roots: 1, 1, 2, d = 0
  d = P::cubicDiscriminant( -6.f, 15.f, -12.f, 3.f);
  P::rootsCubicComplex(C(-6), C(15), C(-12), C(3), &cr1, &cr2, &cr3);
  r &= d   == 0;
  r &= cr1 == C(1, 0);
  r &= cr2 == C(1, 0);
  r &= cr3 == C(2, 0);

  // -15 + 27*x - 15*x^2 + 3*x^3, roots: 1, 2-i, 2+i, d = -1296
  d = P::cubicDiscriminant( -15.f, 27.f, -15.f, 3.f);
  P::rootsCubicComplex(C(-15), C(27), C(-15), C(3), &cr1, &cr2, &cr3);
  r &= d   == -1296;
  r &= cr1 == C(1,  0);
  r &= cr2 == C(2, -1);
  r &= cr3 == C(2, +1);
  // gives the wrong result - it's not even a real root and a complex conjugate pair



  // todo: fix warnings, implement cubic qnd quartic formulas

  // quartic p(x) = x^4 - 7x^3 + 21*x^2 - 23*x - 52,  roots: 2+3i, 2-3i, -1, 4

  return r;
}


double minusSquare(double x)
{
  return -x*x;
}

bool minimizerUnitTest()
{
  bool  ok   = true;

  using Real = double;
  using MM   = RAPT::rsMinimizer1D<Real>;

  Real x, y;                      // function in/out values
  std::function<Real(Real)> f;    // the function to find the minimum of


  f = [](Real x){ return x*x; };  // f(x) = x^2 with minimum at x = 0


  Real xMin;

  xMin = goldenRatioMethodMax(&minusSquare, -1.0, +1.0);
  xMin = goldenRatioMethodMax(&minusSquare, -1.0, +2.0);


  //xMin = MM::goldenSection(f, -1.0, +1.0);
  //ok &= xMin == 0.0;  // maybe we need a tolerance

  xMin = MM::goldenSectionMin(f, -2.0, +3.0);
  ok &= xMin == 0.0;  








  return ok;
}

