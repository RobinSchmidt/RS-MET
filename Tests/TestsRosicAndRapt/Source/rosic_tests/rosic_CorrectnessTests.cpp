#include "rosic_CorrectnessTests.h"
using namespace rotes;

void rotes::testAllRosicClasses()
{
  testRosicString();
  testRosicBasics();
  testRosicFile();
  testRosicFilter();
  testRosicGenerators();
  testRosicModulators();
  testRosicNonRealTime();
  testRosicOthers();
}

void rotes::testRosicAnalysis()
{
  testOscilloscopeBuffer();
}

void rotes::testRosicBasics()
{
  //testBinomialCoefficients();
  //testMathFunctions();
  //testWindowFunctions();
  //testInterpolation();
}

void rotes::testRosicFile()
{
  //testFileTextReadWrite();
}

void rotes::testRosicEffects()
{
  testFastGeneralizedHadamardTransform();
  testFeedbackDelayNetwork();
}

void rotes::testRosicFilter()
{
  //testLadderFilter();
  //testModalFilter();
  //testModalFilterWithAttack();
  //testBiquadPhasePlot();
  //testFiniteImpulseResponseDesigner();
  //testConvolverPartitioned();
  //testFiniteImpulseResponseFilter();
  //testFilterAnalyzer();
  //testBiquadCascade();
  //testCrossover4Way();
  //testCrossover4Way2();
  //testSlopeFilter();
  //testPrototypeDesigner();
  //testLowpassToLowshelf();
  //testBesselPrototypeDesign();
  //testPapoulisPrototypeDesign();
  //testEngineersFilter();
  testPoleZeroMapping();

  // reference output production for RSLib:
  highOrderFilterPolesAndZeros();

}

void rotes::testRosicGenerators()
{
  //testOscillatorStereo();
  testLorentzSystem();
}

void rotes::testRosicModulators()
{
  testConsecutiveExponentialDecay();
}

void rotes::testRosicMath()
{
  testComplexSqrt();
  testCubicCoeffsTwoPointsAndDerivatives();
  testPolynomialDiffAndInt();
  testPolynomialComposition();
  testPolynomialWeightedSum();
  testPolynomialIntegrationWithPolynomialLimits();
  //testPolynomialRootFinder();
  testLinLogEquationSolver();
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
  testSlewRateLimiterLinear();
}
