#pragma once

// get rid, move content of those files here:
#include <string> // get rid - we don't need the report strings anymore

//#include <rs_testing/rs_testing.h>



#include "PortedFromRSLib/UnitTestsRSLib.h" // clean up and copy contents to here



#include "rosic_CorrectnessTests.h" // copy contents her, get rid of file

namespace rotes  // maybe get rid of this namespace
{




// unit test drivers:
void testAllRosicClasses();
void testRosicAnalysis();
void testRosicBasics();
void testRosicFile();
void testRosicEffects();
void testRosicFilter();
void testRosicGenerators();
void testRosicModulators();
void testRosicMath();
void testRosicNumerical();
void testRosicNonRealTime();
void testRosicOthers();
}






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