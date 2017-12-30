#include "MathUnitTests.h"

bool interpolatingFunctionUnitTest()
{	
  bool r = true;      // test result

  rsInterpolatingFunctionF ipf;

  float y;

  y = ipf.getValue(10.f);
  r &= y == 0.f;

  ipf.addDataPoint(0.f, 2.f);
  y = ipf.getValue(-1.f);  r &= y == 2.f;
  y = ipf.getValue( 0.f);  r &= y == 2.f;
  y = ipf.getValue(+1.f);  r &= y == 2.f;



  return r;
}


