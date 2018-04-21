#include "MathUnitTests.h"

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
//#include "Math/NumberTheoryTests.cpp"
//#include "Math/RealFunctionFunctionTests.cpp"
//#include "Math/TransformsTests.cpp"
//#include "Math/StatisticsTests.cpp"  // there is no such file - why?

bool coordinateMapperUnitTest()
{
  bool r = true;      // test result

  rsCoordinateMapper2DF mapper;

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

bool interpolatingFunctionUnitTest()
{	
  bool r = true;      // test result

  rsNodeBasedFunctionF ipf; // rename to nbf

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
bool rootFinderUnitTest()
{
  bool r = true;                 // test result
  float x, y;                    // function in/out values
  std::function<float(float)> f; // the function to find the roots of

  // create example (lambda) function with roots at -1,+1,+2 and verify positions of roots:
  f = [] (float x)->float { return (x+1)*(x-1)*(x-2); };
  y = f(-1.f); r &= y == 0.f;
  y = f( 1.f); r &= y == 0.f;
  y = f( 2.f); r &= y == 0.f;

  // find the roots via bisection:
  x = rsRootFinderF::bisection(f, -1.3f, -0.8f); r &= x == -1.f;
  x = rsRootFinderF::bisection(f,  0.8f,  1.3f); r &= x ==  1.f;
  x = rsRootFinderF::bisection(f,  1.7f,  2.2f); r &= x ==  2.f;

  // find the roots via false position:
  x = rsRootFinderF::falsePosition(f, -1.3f, -0.8f); r &= x == -1.f;
  x = rsRootFinderF::falsePosition(f,  0.8f,  1.3f); r &= x ==  1.f;
  x = rsRootFinderF::falsePosition(f,  1.7f,  2.2f); r &= x ==  2.f;

  // use a function pointer:
  f = &testFunction;
  x = rsRootFinderF::falsePosition(f, -1.3f, -0.8f); r &= x == -1.f;
  x = rsRootFinderF::falsePosition(f,  0.8f,  1.3f); r &= x ==  1.f;
  x = rsRootFinderF::falsePosition(f,  1.7f,  2.2f); r &= x ==  2.f;

  // use a function object (functor):
  //TestFunctor functor; f = functor;
  f = TestFunctor();  // shorter syntax
  x = rsRootFinderF::falsePosition(f, -1.3f, -0.8f); r &= x == -1.f;
  x = rsRootFinderF::falsePosition(f,  0.8f,  1.3f); r &= x ==  1.f;
  x = rsRootFinderF::falsePosition(f,  1.7f,  2.2f); r &= x ==  2.f;

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
  P::rootsQuadraticReal(6, -9, 3, &rr1, &rr2);
  r &= rr1 == 1.f;
  r &= rr2 == 2.f;

  // 15 - 12*x + 3*x^2 , roots: 2-i, 2+i
  P::rootsQuadraticComplex(15, -12, 3, &cr1, &cr2);
  r &= cr1 == C(2, -1);
  r &= cr2 == C(2, +1);
  P::rootsQuadraticReal(15, -12, 3, &rr1, &rr2); // should give real parts
  r &= rr1 == 2;
  r &= rr2 == 2;

  // x^3 - x, roots: -1, 0, +1, d = 4
  d = P::discriminant(0,  -1, 0, 1);
  P::rootsCubicComplex(0, -1, 0, 1, &cr1, &cr2, &cr3);
  r &= d   == 4;
  r &= cr1 == C(-1, 0);
  r &= cr2 == C( 0, 0);
  r &= cr3 == C(+1, 0);
  // has roundoff error, otherwise ok

  // x^3 + x, roots: 0, -i, +i, d = -4
  d = P::discriminant( 0, 1, 0, 1);
  P::rootsCubicComplex(0, 1, 0, 1, &cr1, &cr2, &cr3);
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
  d = P::discriminant( -18, 33, -18, 3);
  P::rootsCubicComplex(-18, 33, -18, 3, &cr1, &cr2, &cr3);
  r &= d   == 324;
  r &= cr1 == C(1, 0);
  r &= cr2 == C(2, 0);
  r &= cr3 == C(3, 0);

  // -6 + 15*x - 12*x^2 + 3*x^3, roots: 1, 1, 2, d = 0
  d = P::discriminant( -6, 15, -12, 3);
  P::rootsCubicComplex(-6, 15, -12, 3, &cr1, &cr2, &cr3);
  r &= d   == 0;
  r &= cr1 == C(1, 0);
  r &= cr2 == C(1, 0);
  r &= cr3 == C(2, 0);

  // -15 + 27*x - 15*x^2 + 3*x^3, roots: 1, 2-i, 2+i, d = -1296
  d = P::discriminant( -15, 27, -15, 3);
  P::rootsCubicComplex(-15, 27, -15, 3, &cr1, &cr2, &cr3);
  r &= d   == -1296;
  r &= cr1 == C(1,  0);
  r &= cr2 == C(2, -1);
  r &= cr3 == C(2, +1);
  // gives the wrong result - it's not even a real root and a complex conjugate pair



  // todo: fix warnings, implement cubic qnd quartic formulas
  
  // quartic p(x) = x^4 - 7x^3 + 21*x^2 - 23*x - 52,  roots: 2+3i, 2-3i, -1, 4

  return r;
}
