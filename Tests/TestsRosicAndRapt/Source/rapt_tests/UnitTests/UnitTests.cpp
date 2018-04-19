#include "UnitTests.h"

using namespace RAPT;

/*
template<class T>
template<class S>
void Foo<T>::bar(S x)
{

}
*/

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

  // Data:
  passed &= runUnitTest(&float64x2UnitTest,        "rsFloat64x2");
  passed &= runUnitTest(&complexFloat64x2UnitTest, "std::complex<rsFloat64x2>");

  // Math:
  passed &= runUnitTest(&coordinateMapperUnitTest,      "rsCoordinateMapper2D");
  passed &= runUnitTest(&interpolatingFunctionUnitTest, "rsInterpolatingFunction");
  passed &= runUnitTest(&rootFinderUnitTest,            "rsRootFinder");

  passed &= runUnitTest(&testLinearAlgebra,             "rsLinearAlgebra");
  passed &= runUnitTest(&testPolynomial,                "rsPolynomial");

  passed &= runUnitTest(&polynomialRootsUnitTest,       "rsPolynomial: root finding"); // absorb in rsPolynomial



  // Filters:
  //passed &= runUnitTest(&prototypeDesignUnitTest, "rsPrototypeDesigner");

  // Visualization:
  passed &= runUnitTest(&imagePainterUnitTest, "rsImagePainter");

  //...
  //...more to come...

  return passed;
}
