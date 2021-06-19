// ToDo: Try to get rid of this file by merging the unit tests and experiments here with those for
// rapt. At the moment, all these testRosic*() functions call a mixture of unit tests (returning a 
// bool) and experiments (returning nothing, but popping up plots and/or writing wavefiles, etc). 
// This should be disentangled. Currently, it is controlled by a switch, whether the experiments
// should run....

bool rotes::testAllRosicClasses()
{
  bool runExperiments = false;
  // This is the switch that decides whether or not the experiments should run (see comment above).
  // When the experiments are disentangled from the unit tests, this can go away.

  rsWarning("In testAllRosicClasses(), not all tests are updated yet.");


  bool ok = true;
  ok &= testRosicBasics(runExperiments);
  ok &= testRosicFile();
  ok &= testRosicFilter(runExperiments);
  testRosicGenerators(runExperiments);    // has no unit tests yet
  testRosicModulators(runExperiments);    // ditto
  ok &= testRosicNonRealTime();
  testRosicOthers();
  testRosicAnalysis();
  testRosicEffects();
  testRosicNumerical();
  testRosicMath();

  ok &= testRosicString();  // Fails! reason: double/string roundtrip..see comments in the tests

  return ok;
}

void rotes::testRosicAnalysis()
{
  testOscilloscopeBuffer();  // creates plot
}

bool rotes::testRosicBasics(bool runExperiments)
{
  if(runExperiments)
  {
    testWindowFunctions();
    testInterpolation();
  }

  bool ok = true;
  ok &= testBinomialCoefficients();  // obsolete thx to rapt, now?
  ok &= testMathFunctions();
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

bool rotes::testRosicFilter(bool runExperiments)
{
  if(runExperiments)
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
  }
  bool ok = true;
  ok &= testConvolverPartitioned();
  ok &= testFiniteImpulseResponseFilter();  // fails! convolver imp-resp update is conditional
  return ok;
}

void rotes::testRosicGenerators(bool runExperiments)
{
  if(runExperiments)
  {
    testOscillatorStereo();
    testLorentzSystem();
  }
  // has not unit tests yet.
}

void rotes::testRosicModulators(bool runExperiments)
{
  if(runExperiments)
  {
    testConsecutiveExponentialDecay();
  }
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

bool rotes::testRosicNonRealTime()
{
  bool ok = true;
  ok &= testMinimumPhaseReconstruction(); // maybe move to testRosicMath
  return ok;
}

void rotes::testRosicOthers()
{
  testSlewRateLimiterLinear();  // creates plot
}
