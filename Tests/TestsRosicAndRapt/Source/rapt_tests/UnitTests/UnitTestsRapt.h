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
bool simdUnitTest();

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
bool interpolationUnitTest();
bool rootFinderUnitTest();
bool polynomialRootsUnitTest(); // the new explicit formulas - move to PolynomialUnitTests
bool minimizerUnitTest();       // 1D function minimization


// Misc:
bool blepUnitTest();            // move to GeneratorTests
bool syncUnitTest();
bool spectrogramUnitTest();
bool sineModelingUnitTest();
bool analysisUnitTest();        // maybe move to AnalysisTests

#endif