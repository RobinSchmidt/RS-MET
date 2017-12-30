#include "UnitTests.h"

bool runUnitTest(bool (*test)(), const string& name)
{
  //cout << "Testing: " + name + ": ";
  cout << name + ": ";
  bool passed = test();
  rsAssert(passed); // break, if test fails
  if(passed)
    cout << "Passed\n";
  else
    cout << "Failed\n";
  return passed;
}

bool runAllUnitTests()
{
  bool passed = true;  // test result

  // Filters:
  passed &= runUnitTest(&interpolatingFunctionUnitTest, "rsInterpolatingFunction");

  // Filters:
  passed &= runUnitTest(&prototypeDesignUnitTest, "rsPrototypeDesigner");

  // Visualization:
  passed &= runUnitTest(&imagePainterUnitTest, "rsImagePainter");

  //...
  //...more to come...

  return passed;
}
