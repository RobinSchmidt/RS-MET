// ToDo: Try to get rid of this file by merging the unit tests and experiments here with those for
// rapt. 

/*
void rotes::testRosicAnalysis()
{

}
*/

bool rotes::testRosicBasics()
{
  bool ok = true;
  ok &= testBinomialCoefficients();  // obsolete thx to rapt, now?
  ok &= testMathFunctions();
  return ok;
}

bool rotes::testRosicFile()
{
  bool ok = true;
  ok &= testFileText();
  ok &= testFileWave();
  return ok;
}





bool rotes::testRosicEffects()
{
  bool ok = true;
  ok &= testFastGeneralizedHadamardTransform();
  //ok &= testFeedbackDelayNetwork();  // This is not an actual unit test! ToDo: move/rename/etc.
  ok &= testAllpassDelay();
  ok &= testFreqShifter();
  ok &= testMultiComp();
  return ok;
}

bool rotes::testRosicFilter()
{
  bool ok = true;
  ok &= testConvolverPartitioned();
  ok &= testFiniteImpulseResponseFilter();  // fails! convolver imp-resp update is conditional
  return ok;
}

/*
void rotes::testRosicGenerators()
{

}

void rotes::testRosicModulators()
{

}
*/

bool rotes::testRosicMath()
{
  bool ok = true;
  ok &= testComplexSqrt();
  ok &= testCubicCoeffsTwoPointsAndDerivatives();
  ok &= testPolynomialDiffAndInt();
  ok &= testPolynomialComposition();
  ok &= testPolynomialWeightedSum();
  ok &= testPolynomialIntegrationWithPolynomialLimits();

  ok &= testPolynomialRootFinder();
  // I had to set it to a quite high tolerance to make it pass. -> Figure out why and try to 
  // reduce the tolerance. It seems that it did pass with lower tolerances at some point (there
  // are code remnants of running the test with lower tolerance)

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
void rotes::testRosicOthers()
{

}
*/
