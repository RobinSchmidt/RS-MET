#include "UnitTests.h"

bool runUnitTest(bool (*test)(), const string& name)
{
  bool r = test();

  // todo: print result

  return r;
}

bool runAllUnitTests()
{
  bool r = true;  // test result

  // Filters:
  r &= runUnitTest(&prototypeDesignUnitTest, "rsPrototypeDesigner");

  // Visualization:
  r &= runUnitTest(&imagePainterUnitTest, "rsImagePainter");

  /*r &= imagePainterUnitTest();*/
  //...
  //...more to come...

  return r;
}
