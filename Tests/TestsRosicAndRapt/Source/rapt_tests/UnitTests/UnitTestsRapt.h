#ifndef UNITTESTS_INCLUDED
#define UNITTESTS_INCLUDED

// The function that runs all the unit tests:
bool runUnitTestsRapt(); 

// Data:
bool arrayUnitTest();
bool testBufferFunctions();      // rename or absorb in arrayUnitTests
bool testSortAndSearch();        // rename to sortAndSearchUnitTest
bool binaryHeapUnitTest();
bool doubleEndedQueueUnitTest();
bool float64x2UnitTest();
bool float32x4UnitTest();
bool complexFloat64x2UnitTest();

// Drawing:
bool triangleRasterization();
bool triangleRasterization2();  // absorb in function above
bool imagePainterUnitTest();    // implementation is in extra file - maybe mae a GraphicsUnitTests.cpp file for all
bool colorUnitTest();

// Filter:
bool prototypeDesignUnitTest();
bool filterSpecUnitTest();
bool movingMaximumUnitTest();

// Math:
bool coordinateMapperUnitTest();
bool correlationUnitTest();
bool fitRationalUnitTest();
bool interpolatingFunctionUnitTest();
bool resampleNonUniform();
bool rootFinderUnitTest();
bool polynomialRootsUnitTest(); // the new explicit formulas - move to PolynomialUnitTests

// Misc:
bool blepUnitTest();  // move to GeneratorTests
bool syncUnitTest();
bool spectrogramUnitTest();
bool sineModelingUnitTest();

#endif