#pragma once


// make similar files PerformanceTestsRosic, ExperimentsRosic, ExamplesRosic


//#include "PortedFromRSLib/UnitTestsRSLib.h" 
// Clean up and copy contents to here - hmm -it's actually not used anymore - figure out, if the 
// tests from there are already integrated somewhere here - if not, do it


namespace rotes  // maybe get rid of this namespace
{

// todo: many of the tests here seem to not be unit tests but rather experiments (they don't return
// a bool) - maybe sort the tests according to their nature ...hmm - some work with raising 
// assertions - that was the old way - they should be adapted to return a bool instead of breaking
// in case of failure

// string
bool testRosicString(); // all tests for rosic::String
//void testStringComparisons();   // why commented?
bool testCharacterComparisons();
bool testStringBufferCopying();
bool testStringIntConversions(int numIterations = 10000);
bool testStringDoubleConversions();
bool testStringDoubleConversionsRandom(int numIterations = 10000); 
bool testStringDoubleConversionsSpecialValues(); 
bool testStringDoubleConversionsDenormals(); 
bool testStringDoubleConversionsLarge(); 
bool testStringDoubleConversionsGeometricProgression(double start, double factor);
//void testStringConcatenation();
//void testStringComparison();...
rosic::rsString createStringWithAllCharacters();
rosic::rsString createStringWithAllPrintableCharacters();

// file:
bool testFileTextReadWrite();  // tests, if we can write a string into a file and retrieve it 
                               // again, the string must not contain non-printable characters

// effects:
bool testFastGeneralizedHadamardTransform();
bool testFeedbackDelayNetwork();
bool testMultiComp();

// filters:
void testLadderFilter();                   // is esperiment
void testModalFilter();                    // is esperiment
void testModalFilterWithAttack();          // is esperiment
void testBiquadPhasePlot();                // is esperiment
void testFiniteImpulseResponseDesigner();  // is esperiment
bool testConvolverPartitioned();
bool testFiniteImpulseResponseFilter();
void testFilterAnalyzer();
void testBiquadCascade();
void testCrossover4Way();
void testCrossover4Way2();
void testCrossoverNewVsOld();
void testSlopeFilter();
void testPrototypeDesigner();
void testLowpassToLowshelf();
void testBesselPrototypeDesign();
void testPapoulisPrototypeDesign();
void testEngineersFilter();
void testPoleZeroMapping();
void highOrderFilterPolesAndZeros();

// generators:
void testOscillatorStereo();
void testLorentzSystem();
bool testSnowflake();
bool testResetter();
void testTurtleReverse();
void testTurtleSource();

// modulators:
void testConsecutiveExponentialDecay();

// analysis:
void testOscilloscopeBuffer();

// basics:
bool testBinomialCoefficients();    // obsolete thx to rapt, now?
bool testMathFunctions();
void testWindowFunctions();
void testInterpolation();
void testHermiteTwoPoint1();
void testHermiteTwoPoint2();
void testHermiteTwoPoint3();
void testHermiteTwoPointM();
rosic::Matrix createHermiteInterpolatorImpulseResponses(int inLength, int oversampling, const int M[5], double shape);
void plotOneSidedInterpolatorContinuousResponses(int M[5], double shape);
void plotOneSidedInterpolatorPolyphaseResponses(int M, double shape, double d[5]);
void testAsymmetricPolynomialInterpolatorsOld();


// math:
bool testComplexSqrt();
bool testCubicCoeffsTwoPointsAndDerivatives();
bool testPolynomialDiffAndInt();
bool testPolynomialComposition();
bool testPolynomialWeightedSum();
bool testPolynomialIntegrationWithPolynomialLimits();
void testPolynomialRootFinder();
void testLinLogEquationSolver();
void testLinLogEquationSolverOld();
void testLinearSystemSolver();

// numerical (maybe merge with math - or get rid):
void testUnivariateScalarFunction();
void testUnivariateRootFinder();

// non-realtime:
bool testMinimumPhaseReconstruction();

// others:
void testSlewRateLimiterLinear();

// unit test drivers:
bool testAllRosicClasses();
void testRosicAnalysis();
bool testRosicBasics(bool runExperiments = false);
bool testRosicFile();
void testRosicEffects();
bool testRosicFilter(bool runExperiments = false);
void testRosicGenerators(bool runExperiments = false);
void testRosicModulators(bool runExperiments = false);
void testRosicMath();
void testRosicNumerical();
bool testRosicNonRealTime();
void testRosicOthers();
}



bool runUnitTestsRosic();

bool testFilterPolynomials();
bool testHighOrderFilter();

bool testModalFilter2();
bool testModalSynth();

bool testBandwidthConversions();
bool testSincInterpolation();
bool testSineParameters();
bool testZeroCrossingFinder();

bool testNumberManipulations();
bool testDoubleIntConversions();  // get rid - is subtest of  testNumberManipulations
bool testExponentExtraction();    // dito


bool testAutoCorrelationPitchDetector();
bool testTypeSizes();