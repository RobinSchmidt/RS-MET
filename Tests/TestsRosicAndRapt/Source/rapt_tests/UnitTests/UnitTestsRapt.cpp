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

  // these tests should go into UnitTestsRosic.cpp:
  ok &= runUnitTest(&samplerDataUnitTest,         "rsSamplerData");
  ok &= runUnitTest(&samplerEngineUnitTest,       "rsSamplerEngine");
  ok &= runUnitTest(&samplerEngineUnitTestFileIO, "rsSamplerEngine - File I/O");
  return ok;




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
  ok &= runUnitTest(&fitRationalUnitTest,            "fit rational");   // fails on linux ("illegal instruction") - encounters singular matrix
  ok &= runUnitTest(&interpolatingFunctionUnitTest,  "rsInterpolatingFunction");
  ok &= runUnitTest(&resampleNonUniform,             "resampleNonUniform");
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
  //ok &= runUnitTest(&polynomialRootsUnitTest,       "rsPolynomial: root finding"); // fails! ...absorb in rsPolynomial


  // Filters:
  //ok &= runUnitTest(&prototypeDesignUnitTest, "rsPrototypeDesigner"); // why commented?
  ok &= runUnitTest(&filterSpecUnitTest,     "rsFilterSpecification (BA/ZPK)");
  ok &= runUnitTest(&movingMaximumUnitTest,  "moving maximum filter");
  ok &= runUnitTest(&movingQuantileUnitTest, "moving quantile filter"); // under construction
  ok &= runUnitTest(&ladderUnitTest,         "rsLadder"); 


  // Visualization:
  ok &= runUnitTest(&imagePainterUnitTest,   "rsImagePainter");
  ok &= runUnitTest(&triangleRasterization,  "Triangle Rasterization");
  //ok &= runUnitTest(&triangleRasterization2, "Triangle Rasterization 2"); // merge
  ok &= runUnitTest(&colorUnitTest,          "rsColor");
  // move down later

  // Generators:
  ok &= runUnitTest(&samplerDataUnitTest,         "rsSamplerData");
  ok &= runUnitTest(&samplerEngineUnitTest,       "rsSamplerEngine");
  ok &= runUnitTest(&samplerEngineUnitTestFileIO, "rsSamplerEngine - File I/O");

  // Misc:
  ok &= runUnitTest(&blepUnitTest,  "Blit/Blep/Blamp");  // move to generator unit tests
  ok &= runUnitTest(&syncUnitTest,  "Osc Sync");         // dito


  ok &= runUnitTest(&spectrogramUnitTest,   "rsSpectrogramProcessor");
  ok &= runUnitTest(&sineModelingUnitTest,  "SineModeling");
  // todo: test rsCycleMarkFinder: give it a sine-sweep as input and check, if mark-deltas are 
  // within sane limits ...but maybe that's better for an experiment
  //...
  //...more to come...

  if(ok) std::cout << "RAPT: OK\n";
  else   std::cout << "RAPT: !!!!----> F A I L E D <----!!!!\n";
  std::cout << "\n";
  return ok;
}
