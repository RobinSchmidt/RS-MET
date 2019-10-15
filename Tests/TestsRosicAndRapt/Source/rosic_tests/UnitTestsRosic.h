#pragma once

// get rid, move content of those files here:
#include <string> // get rid - we don't need the report strings anymore

//#include <rs_testing/rs_testing.h>



#include "PortedFromRSLib/UnitTestsRSLib.h" // clean up and copy contents to here



#include "rosic_CorrectnessTests.h" // copy contents her, get rid of file

namespace rotes  // maybe get rid of this namespace
{

// string
void testRosicString(); // all tests for rosic::String
//void testStringComparisons();
void testCharacterComparisons();
void testStringBufferCopying();
void testStringIntConversions(int numIterations = 10000); 
void testStringDoubleConversions(); 
void testStringDoubleConversionsRandom(int numIterations = 10000); 
void testStringDoubleConversionsSpecialValues(); 
void testStringDoubleConversionsDenormals(); 
void testStringDoubleConversionsLarge(); 
void testStringDoubleConversionsGeometricProgression(double start, double factor);
//void testStringConcatenation();
//void testStringComparison();...
rosic::rsString createStringWithAllCharacters();
rosic::rsString createStringWithAllPrintableCharacters();

// file:
void testFileTextReadWrite();  // tests, if we can write a string into a file and retrieve it 
                               // again, the string must not contain non-printable characters

// effects:
bool testFastGeneralizedHadamardTransform();
bool testFeedbackDelayNetwork();
bool testMultiComp();

// filters:
void testLadderFilter();
void testModalFilter();
void testModalFilterWithAttack();
void testBiquadPhasePlot();
void testFiniteImpulseResponseDesigner();
void testConvolverPartitioned();
void testFiniteImpulseResponseFilter();
void testFilterAnalyzer();
void testBiquadCascade();
void testCrossover4Way();
void testCrossover4Way2();
void testCrossoverNewVsOld();
void testSlopeFilter();
void testPrototypeDesigner();
void testLowpassToLowshelf();
void testBesselPrototypeDesign();
void testPapoulisPrototypeDesign();
void testEngineersFilter();
void testPoleZeroMapping();
void highOrderFilterPolesAndZeros();


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