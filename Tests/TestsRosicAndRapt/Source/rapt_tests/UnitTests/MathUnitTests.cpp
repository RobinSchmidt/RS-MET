#include "MathUnitTests.h"

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

bool rootFinderUnitTest()
{
  bool r = true;  // test result
  float x, y;     // function in/out values

  // create example function with roots at -1,+1,+2 and verify positions of roots:
  std::function<float(float)> f = [] (float x)->float { return (x+1)*(x-1)*(x-2); };
  y = f(-1.f); r &= y == 0.f;
  y = f( 1.f); r &= y == 0.f;
  y = f( 2.f); r &= y == 0.f;

  // find the roots via bisection:
  x = rsRootFinderFF::bisection(f, -1.3f, -0.8f); r &= x == -1.f;
  x = rsRootFinderFF::bisection(f,  0.8f,  1.3f); r &= x ==  1.f;
  x = rsRootFinderFF::bisection(f,  1.7f,  2.2f); r &= x ==  2.f;

  return r;
}


