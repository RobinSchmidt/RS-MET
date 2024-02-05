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

// file:
bool testFileText();
bool testFileWave();  

// effects:
void allpassDisperser();           // is experiment
void allpassDelay();               // is experiment
void allpassDelayChain();          // is experiment
bool testFeedbackDelayNetwork();   // is experiment
bool testFastGeneralizedHadamardTransform();
bool testAllpassDelay();
bool testFreqShifter();
bool testMultiComp();
void spectralFilter();
void formantShifter();
void spectralShifter();

// filters:
void testLadderFilter();                   // is experiment
void testModalFilter();                    // is experiment
void testModalFilterWithAttack();          // is experiment
void testBiquadPhasePlot();                // is experiment
void testFiniteImpulseResponseDesigner();  // is experiment
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
void testCombustionEngine();
bool testSnowflake();
bool testResetter();
void testTurtleReverse();
void testTurtleSource();
void testSamplerEngine();

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
bool testPolynomialRootFinder();
void testLinLogEquationSolver();
void testLinLogEquationSolverOld();
bool testLinearSystemSolver();

// numerical (maybe merge with math - or get rid):
bool testUnivariateScalarFunction();
//void testUnivariateRootFinder();    // empty

// non-realtime:
bool testMinimumPhaseReconstruction();

// others:
void testSlewRateLimiterLinear();

// unit test drivers:
//void testRosicAnalysis();
bool testRosicBasics();
bool testRosicFile();
bool testRosicEffects();
bool testRosicFilter();
//void testRosicGenerators();
//void testRosicModulators();
bool testRosicMath();
bool testRosicNumerical();
bool testRosicNonRealTime();
//void testRosicOthers();
// these runExperiments parameters should go away when the disentanglement of the unit tests from 
// the experiments is complete
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