//#include "rosic_CorrectnessTests.h"
//using namespace rotes;

// try to get rid of this file

bool rotes::testAllRosicClasses()
{
  // ToDo: split into unit tests and experiments, unit tests should return a bool
  bool ok = true;

  //printf("Warning: in testAllRosicClasses, not all tests are updated yet.\n");
  rsWarning("In testAllRosicClasses(), not all tests are updated yet.");

  ok &= testRosicString();
  ok &= testRosicBasics();
  testRosicFile();
  testRosicFilter();
  testRosicGenerators();
  testRosicModulators();
  testRosicNonRealTime();
  testRosicOthers();
  testRosicAnalysis();
  testRosicEffects();
  testRosicNumerical();
  testRosicMath();

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
  testWindowFunctions();
  testInterpolation();
  return ok;
}

void rotes::testRosicFile()
{
  testFileTextReadWrite();
}

void rotes::testRosicEffects()
{
  testFastGeneralizedHadamardTransform(); // returns bool 
  testFeedbackDelayNetwork();             // writes wave file
}

void rotes::testRosicFilter()
{
  testLadderFilter();
  testModalFilter();
  testModalFilterWithAttack();
  testBiquadPhasePlot();
  testFiniteImpulseResponseDesigner();
  testConvolverPartitioned();
  testFiniteImpulseResponseFilter();
  testFilterAnalyzer();
  testBiquadCascade();
  //testCrossover4Way(); // linker error
  testCrossover4Way2();
  testSlopeFilter();
  testPrototypeDesigner();
  testLowpassToLowshelf();
  testBesselPrototypeDesign();
  testPapoulisPrototypeDesign();
  testEngineersFilter();
  testPoleZeroMapping();

  // reference output production for RSLib:
  highOrderFilterPolesAndZeros();

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
