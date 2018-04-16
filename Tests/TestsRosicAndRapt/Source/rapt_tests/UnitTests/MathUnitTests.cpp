#include "MathUnitTests.h"

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

bool polynomialRootsUnitTest()
{
  bool r = true;

  typedef RAPT::rsPolynomial<float> P;
  typedef std::complex<float> C;

  float rr1, rr2;                     // real roots
  std::complex<float> cr1, cr2, cr3;  // complex roots

  // 3*x^2 - 9*x + 6, roots: 1, 2
  P::rootsQuadraticReal(3, -9, 6, &rr1, &rr2);
  r &= rr1 == 1.f;
  r &= rr2 == 2.f;

  // 3*x^2 - 12*x + 15, roots: 2-i, 2+i
  P::rootsQuadraticComplex(3, -12, 15, &cr1, &cr2);
  r &= cr1 == C(2, -1);
  r &= cr2 == C(2, +1);
  P::rootsQuadraticReal(3, -12, 15, &rr1, &rr2); // should give real parts
  r &= rr1 == 2;
  r &= rr2 == 2;

  // -18 + 33*x - 18*x^2 + 3*x^3, roots: 1, 2, 3
  P::rootsCubicComplex(-18, 33, -18, 3, &cr1, &cr2, &cr3);
  r &= cr1 == C(1, 0);
  r &= cr2 == C(2, 0);
  r &= cr3 == C(3, 0);

  // -15 + 27*x - 15*x^2 + 3*x^3, roots: 1, 2-i, 2+i
  P::rootsCubicComplex(-15, 27, -15, 3, &cr1, &cr2, &cr3);
  r &= cr1 == C(1, 0);
  r &= cr2 == C(2, -1);
  r &= cr3 == C(2, +1);



  // todo: fix warnings, implement cubic qnd quartic formulas
  
  // quartic p(x) = x^4 - 7x^3 + 21*x^2 - 23*x - 52,  roots: 2+3i, 2-3i, -1, 4

  return r;
}
