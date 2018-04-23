#include "UnitTests.h"

using namespace RAPT;

bool runUnitTest(bool (*test)(), const string& name)
{
  //cout << "Testing: " + name + ": ";
  cout << name + ": ";
  bool passed = test();
  rsAssert(passed); // break, if test fails
  if(passed)
    cout << "Passed\n";
  else
    cout << "!!!!----> F A I L E D <----!!!!\n";
  return passed;
}

bool runAllUnitTests()
{
  bool passed = true;  // test result

  // Data:
  passed &= runUnitTest(&float64x2UnitTest,        "rsFloat64x2");
  passed &= runUnitTest(&complexFloat64x2UnitTest, "std::complex<rsFloat64x2>");

  // Math:
  passed &= runUnitTest(&coordinateMapperUnitTest,       "rsCoordinateMapper2D");
  passed &= runUnitTest(&interpolatingFunctionUnitTest,  "rsInterpolatingFunction");
  passed &= runUnitTest(&rootFinderUnitTest,             "rsRootFinder");

  passed &= runUnitTest(&testVector,                     "rsVector");
  passed &= runUnitTest(&testMatrix,                     "rsMatrix");
  passed &= runUnitTest(&testMiscMath,                   "misc math");
  passed &= runUnitTest(&testLinearAlgebra,              "rsLinearAlgebra");
  passed &= runUnitTest(&testPolynomial,                 "rsPolynomial");
  //passed &= runUnitTest(&polynomialRootsUnitTest,       "rsPolynomial: root finding"); // absorb in rsPolynomial  

  passed &= runUnitTest(&testDifferentialEquationSystem, "rsDifferentialEquationSystem");
  passed &= runUnitTest(&testIntegerFunctions,           "integer functions");
  passed &= runUnitTest(&testMoebiusTransform,           "rsMoebiusTransform");
  passed &= runUnitTest(&testNumberTheory,               "number theory"); 
  passed &= runUnitTest(&testRealFunctions,              "real functions"); 
  passed &= runUnitTest(&testTransforms,                 "transforms");






  // Filters:
  //passed &= runUnitTest(&prototypeDesignUnitTest, "rsPrototypeDesigner");

  // Visualization:
  passed &= runUnitTest(&imagePainterUnitTest, "rsImagePainter");

  //...
  //...more to come...

  return passed;
}
