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
  ok &= testRosicNonRealTime();
  ok &= testRosicEffects();
  ok &= testRosicMath();
  ok &= testRosicNumerical(); // should go to rapt
  ok &= testRosicString();    // Fails! reason: double/string roundtrip..see comments in the tests
  return ok;
}

/*
void rotes::testRosicAnalysis()
{

}
*/

bool rotes::testRosicBasics(bool runExperiments)
{
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

bool rotes::testRosicEffects(bool runExperiments)
{
  bool ok = true;
  ok &= testFastGeneralizedHadamardTransform();
  ok &= testMultiComp();
  return ok;
}

bool rotes::testRosicFilter(bool runExperiments)
{
  bool ok = true;
  ok &= testConvolverPartitioned();
  ok &= testFiniteImpulseResponseFilter();  // fails! convolver imp-resp update is conditional
  return ok;
}

/*
void rotes::testRosicGenerators(bool runExperiments)
{

}

void rotes::testRosicModulators(bool runExperiments)
{

}
*/

bool rotes::testRosicMath(bool runExperiments)
{
  bool ok = true;
  ok &= testComplexSqrt();
  ok &= testCubicCoeffsTwoPointsAndDerivatives();
  ok &= testPolynomialDiffAndInt();
  ok &= testPolynomialComposition();
  ok &= testPolynomialWeightedSum();
  ok &= testPolynomialIntegrationWithPolynomialLimits();
  ok &= testPolynomialRootFinder();
  ok &= testLinearSystemSolver();
  return ok;
}

bool rotes::testRosicNumerical()
{
  bool ok = true;
  ok &= testUnivariateScalarFunction();
  return ok;
  // ToDo: Numerical algos should all be templatized and go to rapt. 
}

bool rotes::testRosicNonRealTime()
{
  bool ok = true;
  ok &= testMinimumPhaseReconstruction(); 
  return ok;
  // maybe move to testRosicMath or better: move the code into rapt
}

/*
void rotes::testRosicOthers(bool runExperiments)
{

}
*/
