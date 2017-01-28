#include "UnitTests.h"

bool runUnitTests()
{
  bool r = true;  // test result

  r &= imagePainterUnitTest();
  //...
  //...more to come...

  return r;
}
