#include "UnitTestsRosic.h"
#include <rs_testing/rs_testing.h>  // get rid

bool runUnitTestsRosic()
{
  bool passed = true;  // test result


  passed &= runUnitTest(&testTypeSizes,          "TypeSizes");
  passed &= runUnitTest(&testExponentExtraction, "ExponentExtraction"); // is oart of numberManipulations
  passed &= runUnitTest(&testFilterPolynomials,  "FilterPolynomials");
  passed &= runUnitTest(&testHighOrderFilter,    "HighOrderFilter");
  passed &= runUnitTest(&testModalFilter2,        "ModalFilter2");
  passed &= runUnitTest(&testModalSynth,          "ModalSynth");  // doesn't do anything useful
  passed &= runUnitTest(&testAutoCorrelationPitchDetector, "AutoCorrPitchDetect");
  passed &= runUnitTest(&testNumberManipulations,          "NumberManipulations"); // fails due to rounding -> figure out



  return passed;
}