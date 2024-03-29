#include "UnitTestsRosic.h"

#include <rs_testing/rs_testing.h>  // get rid - or maybe not?

#include "rosic_AnalysisTests.cpp"
#include "rosic_BasicsTests.cpp"
#include "rosic_CorrectnessTests.cpp"
#include "rosic_EffectsTests.cpp"
#include "rosic_FileTests.cpp"
#include "rosic_FilterTests.cpp"
#include "rosic_GeneratorsTests.cpp"
#include "rosic_MathTests.cpp"
#include "rosic_ModulatorsTests.cpp"
#include "rosic_NonRealTimeTests.cpp"
#include "rosic_NumericalTests.cpp"
#include "rosic_OthersTests.cpp"
#include "rosic_StringTests.cpp"
#include "PortedFromRSLib/UnitTests/FilterTests.cpp"
#include "PortedFromRSLib/UnitTests/MiscAudioTests.cpp"
#include "PortedFromRSLib/UnitTests/ModalTests.cpp"
#include "PortedFromRSLib/UnitTests/NumberManipulationsTests.cpp"
#include "PortedFromRSLib/UnitTests/PitchDetectorTests.cpp"
#include "PortedFromRSLib/UnitTests/TypeSizeTests.cpp"


bool runUnitTestsRosic()
{
  bool ok = true;  // test result

  std::cout << "Running unit tests for rosic\n";

  ok &= runUnitTest(&testRosicBasics,         "Basics");
  ok &= runUnitTest(&testRosicFile,           "File");
  ok &= runUnitTest(&testRosicFilter,         "Filter");
  ok &= runUnitTest(&testRosicNonRealTime,    "NonRealTime");
  ok &= runUnitTest(&testRosicEffects,        "Effects");
  ok &= runUnitTest(&testRosicMath,           "Math");
  ok &= runUnitTest(&testRosicNumerical,      "Numerical");
  ok &= runUnitTest(&testRosicString,         "String");            // Fails! reason: double/string roundtrip..see comments in the tests
  ok &= runUnitTest(&testTypeSizes,           "TypeSizes");
  ok &= runUnitTest(&testNumberManipulations, "NumberManipulations");
  ok &= runUnitTest(&testFilterPolynomials,   "FilterPolynomials");
  ok &= runUnitTest(&testHighOrderFilter,     "HighOrderFilter");   // takes long
  ok &= runUnitTest(&testModalFilter2,        "ModalFilter2");
  ok &= runUnitTest(&testModalSynth,          "ModalSynth");     // doesn't do anything useful
  ok &= runUnitTest(&testAutoCorrelationPitchDetector, "AutoCorrPitchDetect");

  if(ok) std::cout << "rosic: OK\n";
  else   std::cout << "rosic: !!!!----> F A I L E D <----!!!!\n";
  std::cout << "\n";
  return ok;
}