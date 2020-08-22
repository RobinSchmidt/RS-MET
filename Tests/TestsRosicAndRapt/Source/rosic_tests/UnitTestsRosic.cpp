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
  bool passed = true;  // test result

  std::cout << "Running unit tests for rosic\n";

  passed &= runUnitTest(&testTypeSizes,          "TypeSizes");
  passed &= runUnitTest(&testExponentExtraction, "ExponentExtraction"); // is part of numberManipulations
  passed &= runUnitTest(&testNumberManipulations,"NumberManipulations"); // fails due to rounding -> figure out


  passed &= runUnitTest(&testFilterPolynomials,  "FilterPolynomials");
  passed &= runUnitTest(&testHighOrderFilter,    "HighOrderFilter");
  passed &= runUnitTest(&testModalFilter2,        "ModalFilter2");
  passed &= runUnitTest(&testModalSynth,          "ModalSynth");  // doesn't do anything useful
  passed &= runUnitTest(&testAutoCorrelationPitchDetector, "AutoCorrPitchDetect");


  // these need to be adapted
  //testAllRosicClasses();
  //testRosicAnalysis();
  //testRosicBasics();
  //testRosicFile();
  //testRosicEffects();
  //testRosicGenerators();
  //testRosicFilter();
  //testRosicNumerical();
  //testRosicMath();
  //testRosicNonRealTime();
  //testRosicOthers();

  std::cout << "\n";
  return passed;
}