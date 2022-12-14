#include "UnitTestsRapt.h"

#include <rs_testing/rs_testing.h>


#include "DataUnitTests.cpp"
#include "FilterUnitTests.cpp"
#include "GeneratorUnitTests.cpp"
#include "ImageUnitTests.cpp"
#include "MathUnitTests.cpp"
#include "DrawingUnitTests.cpp"
#include "MiscUnitTests.cpp"
#include "SortAndSearchTests.cpp"
#include "BufferFunctionTests.cpp"


bool runUnitTestsRapt()
{
  bool ok = true;  // test result

  std::cout << "Running unit tests for RAPT\n";

  // Test for the currently developed class - it's also again run down there below, but I want
  // the test for the code i'm currently working on to go first for faster edit/build/test cycles, 
  // because some of the test take longer to perfom. So this line is volatile:
  //ok &= runUnitTest(&colorUnitTest,  "rsColor");
  //ok &= runUnitTest(&ladderUnitTest, "rsLadder");

  //// these tests should go into UnitTestsRosic.cpp:
  //ok &= runUnitTest(&stateVariableFilterUnitTest,"rsStateVariableFilter"); 
  //ok &= runUnitTest(&analysisUnitTest,      "Analysis");
  ok &= runUnitTest(&samplerEngineUnitTest,       "rsSamplerEngine");
  return ok;


  // ToDo:
  // ok &= runUnitTest(&testUtilsTest, "Test utilities");  
  // not yet implemented - this is supposed to test the test utility functions and classes such
  // as rsLoggingVector

  // Data:
  ok &= runUnitTest(&arrayUnitTest,            "rsArrayTools and std::vector stuff");
  ok &= runUnitTest(&testBufferFunctions,      "BufferFunctions");  // merge with rsArrayTools tests
  ok &= runUnitTest(&testSortAndSearch,        "SortAndSearch");
  ok &= runUnitTest(&binaryHeapUnitTest,       "rsBinaryHeap");
  ok &= runUnitTest(&ringBufferUnitTest,       "rsDelayBuffer");
  ok &= runUnitTest(&doubleEndedQueueUnitTest, "rsDoubleEndedQueue");
  ok &= runUnitTest(&simdUnitTest,             "SIMD Types");

  // Math:
  ok &= runUnitTest(&coordinateMapperUnitTest,       "rsCoordinateMapper2D");
  ok &= runUnitTest(&interpolationUnitTest,          "Interpolation and curve fitting");
  ok &= runUnitTest(&rootFinderUnitTest,             "rsRootFinder");
  ok &= runUnitTest(&correlationUnitTest,            "correlation");
  ok &= runUnitTest(&testVector,                     "rsVector");
  ok &= runUnitTest(&testMatrix,                     "rsMatrix");
  //ok &= runUnitTest(&testRationalNumber,             "rsRationalNumber"); // is in misc math tests
  ok &= runUnitTest(&testMiscMath,                   "misc math");  // fails on linux ("illegal instruction") - encounters a singular matrix
  ok &= runUnitTest(&testLinearAlgebra,              "rsLinearAlgebra");  // fails on linux ("illegal instruction")
  ok &= runUnitTest(&testDifferentialEquationSystem, "rsDifferentialEquationSystem");
  ok &= runUnitTest(&testIntegerFunctions,           "integer functions");
  ok &= runUnitTest(&testMoebiusTransform,           "rsMoebiusTransform");
  ok &= runUnitTest(&testNumberTheory,               "number theory");
  ok &= runUnitTest(&testRealFunctions,              "real functions");
  ok &= runUnitTest(&testTransforms,                 "transforms");
  ok &= runUnitTest(&testTriangle2D,                 "rsTriangle2D");
  ok &= runUnitTest(&testPoint2D,                    "rsPoint2D");
  ok &= runUnitTest(&testPolygon2D,                  "rsPolygon2D");
  ok &= runUnitTest(&testMultiArray,                 "rsMultiArray");
  ok &= runUnitTest(&testPolynomial,                 "rsPolynomial");

  //ok &= runUnitTest(&polynomialRootsUnitTest,       "rsPolynomial: root finding"); 
  // FAILS! Fix it and the absorb the test in testPolynomial


  // Filters:
  //ok &= runUnitTest(&prototypeDesignUnitTest, "rsPrototypeDesigner"); // why commented?
  ok &= runUnitTest(&filterSpecUnitTest,         "rsFilterSpecification (BA/ZPK)");
  ok &= runUnitTest(&movingMaximumUnitTest,      "moving maximum filter");
  ok &= runUnitTest(&movingQuantileUnitTest,     "moving quantile filter"); // under construction
  ok &= runUnitTest(&ladderUnitTest,             "rsLadder"); 
  ok &= runUnitTest(&stateVariableFilterUnitTest,"rsStateVariableFilter"); 

  // Visualization:
  ok &= runUnitTest(&imagePainterUnitTest,   "rsImagePainter");
  ok &= runUnitTest(&triangleRasterization,  "Triangle Rasterization");
  //ok &= runUnitTest(&triangleRasterization2, "Triangle Rasterization 2"); // merge
  ok &= runUnitTest(&colorUnitTest,          "rsColor");
  // move down later

  // Generators:
  ok &= runUnitTest(&samplerEngineUnitTest,       "rsSamplerEngine");

  // Misc:
  ok &= runUnitTest(&blepUnitTest,  "Blit/Blep/Blamp");  // move to generator unit tests
  ok &= runUnitTest(&syncUnitTest,  "Osc Sync");         // dito


  ok &= runUnitTest(&spectrogramUnitTest,   "rsSpectrogramProcessor");
  ok &= runUnitTest(&sineModelingUnitTest,  "SineModeling");
  ok &= runUnitTest(&analysisUnitTest,      "Analysis");
  // todo: test rsCycleMarkFinder: give it a sine-sweep as input and check, if mark-deltas are 
  // within sane limits ...but maybe that's better for an experiment
  //...
  //...more to come...

  if(ok) std::cout << "RAPT: OK\n";
  else   std::cout << "RAPT: !!!!----> F A I L E D <----!!!!\n";
  std::cout << "\n";
  return ok;
}

// Currently, a unit test is just a function that returns a bool: true if the test has passed and 
// false if it has failed. Maybe at some stage, switch to some sort of unit testing framework. This
// looks interesting (but rquires C++20):
// https://www.youtube.com/watch?v=-qAXShy1xiE
// https://github.com/boost-ext/ut/blob/master/include/boost/ut.hpp
