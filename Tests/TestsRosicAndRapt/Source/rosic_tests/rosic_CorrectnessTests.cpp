//#include "rosic_CorrectnessTests.h"
//using namespace rotes;

// try to get rid of this file

bool rotes::testAllRosicClasses()
{
  // ToDo: split into unit tests and experiments, unit tests should return a bool
  bool ok = true;

  //printf("Warning: in testAllRosicClasses, not all tests are updated yet.\n");
  rsWarning("In testAllRosicClasses(), not all tests are updated yet.");

  ok &= testRosicBasics();
  ok &= testRosicFile();
  ok &= testRosicFilter();
  testRosicGenerators();
  testRosicModulators();
  testRosicNonRealTime();
  testRosicOthers();
  testRosicAnalysis();
  testRosicEffects();
  testRosicNumerical();
  testRosicMath();

  ok &= testRosicString();  // fails! reason: double/string roundtrip..see comments in the tests


  return ok;
}

void rotes::testRosicAnalysis()
{
  testOscilloscopeBuffer();  // creates plot
}

bool rotes::testRosicBasics()
{
  bool ok = true;
  ok &= testBinomialCoefficients();  // obsolete thx to rapt, now?
  ok &= testMathFunctions();
  //testWindowFunctions();           // is experiment, not unit test
  //testInterpolation();             // ditto  
  return ok;
}

bool rotes::testRosicFile()
{
  bool ok = true;
  ok &= testFileTextReadWrite();
  // todo: test wavefile read/write
  return ok;
}

void rotes::testRosicEffects()
{
  testFastGeneralizedHadamardTransform(); // returns bool 
  testFeedbackDelayNetwork();             // writes wave file
}

bool rotes::testRosicFilter()
{
  // To disentangle unit tests from experiments - eventually, the experiments should go elsewhere,
  // probably they should be alled directly in Main.cpp and only the unit tests shall remain here:
  auto runExperiments = []()
  {
    testLadderFilter();
    testModalFilter();
    testModalFilterWithAttack();
    testBiquadPhasePlot();
    testFiniteImpulseResponseDesigner();
    testFilterAnalyzer();
    testBiquadCascade();
    testSlopeFilter();
    testPrototypeDesigner();
    testLowpassToLowshelf();
    testBesselPrototypeDesign();
    testPapoulisPrototypeDesign();
    testEngineersFilter();
    testPoleZeroMapping();
    highOrderFilterPolesAndZeros();  // reference output production for RSLib (obsolete?)
    testCrossover4Way();
    testCrossover4Way2();
  };
  auto runUnitTests = []()
  {
    bool ok = true;
    ok &= testConvolverPartitioned();
    ok &= testFiniteImpulseResponseFilter();  // fails! convolver imp-resp update is conditional
    return ok;
  };

  //runExperiments();
  bool ok = runUnitTests();
  return ok;
}

void rotes::testRosicGenerators()
{
  testOscillatorStereo();
  testLorentzSystem();  // creates a plot
}

void rotes::testRosicModulators()
{
  //testConsecutiveExponentialDecay(); // creates a plot
}

void rotes::testRosicMath()
{
  testComplexSqrt();
  testCubicCoeffsTwoPointsAndDerivatives();
  testPolynomialDiffAndInt();
  testPolynomialComposition();
  testPolynomialWeightedSum();
  testPolynomialIntegrationWithPolynomialLimits();
  testPolynomialRootFinder();
  testLinLogEquationSolver();  // creates plot
  testLinearSystemSolver();
}

void rotes::testRosicNumerical()
{
  testUnivariateScalarFunction();
  testUnivariateRootFinder();
}

void rotes::testRosicNonRealTime()
{
  testMinimumPhaseReconstruction();
}

void rotes::testRosicOthers()
{
  testSlewRateLimiterLinear();  // creates plot
}
