#pragma once


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