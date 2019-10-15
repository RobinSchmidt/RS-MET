#pragma once

// get rid, move content of those files here:
#include <string> // get rid - we don't need the report strings anymore

#include <rs_testing/rs_testing.h>


#include "PortedFromRSLib/UnitTestsRSLib.h" 
#include "rosic_CorrectnessTests.h"


bool runUnitTestsRosic();


bool testFilterPolynomials();
bool testHighOrderFilter();

bool testModalFilter2();
bool testModalSynth();


bool testNumberManipulations();
bool testDoubleIntConversions();  // get rid - is subtest of  testNumberManipulations
bool testExponentExtraction();    // dito


bool testAutoCorrelationPitchDetector();

bool testTypeSizes();