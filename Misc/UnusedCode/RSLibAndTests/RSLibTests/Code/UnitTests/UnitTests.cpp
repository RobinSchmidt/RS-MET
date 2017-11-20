#include "RSCore/RSCoreTests.h"
#include "RSMath/RSMathTests.h"
#include "RSAudio/RSAudioTests.h"
#include "RSGraphics/RSGraphicsTests.h"

bool runTestsCore(std::string &reportString)
{  
  bool testResult = true;

  testResult &= testTypeSizes(reportString);

  reportString += "\nContainers:\n";
  testResult &= testArray(reportString);
  testResult &= testKeyValueMap(reportString);
  testResult &= testList(reportString);
  testResult &= testTree(reportString);

  reportString += "\nUtilities:\n";
  testResult &= testBufferFunctions(reportString);
  testResult &= testNumberManipulations(reportString);
  testResult &= testSortAndSearch(reportString);

  reportString += "\nText:\n";
  testResult &= testString(reportString);
  testResult &= testKeyValueStringTree(reportString);

  reportString += "\nFiles:\n";
  testResult &= testFile(reportString);

  reportString += "\nMisc:\n";
  testResult &= testCallbacks(reportString);

  return testResult;
}

bool runTestsMath(std::string &reportString)
{  
  bool testResult = true;

  reportString += "\nFunctions:\n";
  testResult &= testIntegerFunctions(reportString);
  testResult &= testMoebiusTransform(reportString);
  testResult &= testRealFunctions(reportString);

  reportString += "\nGeometry:\n";
  testResult &= testPoint2D(reportString);
  testResult &= testPolygon2D(reportString);
  testResult &= testTriangle2D(reportString);

  reportString += "\nMath-Misc:\n";
  testResult &= testDifferentialEquationSystem(reportString);
  testResult &= testTransforms(reportString);
  testResult &= testLinearAlgebra(reportString);
  testResult &= testPolynomial(reportString);
  testResult &= testNumberTheory(reportString);

  reportString += "\nMath-Types:\n";
  testResult &= testArbitraryPrecision(reportString);
  testResult &= testComplex(reportString);
  testResult &= testMultiArray(reportString);
  testResult &= testVector(reportString);
  testResult &= testMatrix(reportString);

  testResult &= testMiscMath(reportString);
  // some stuff in Math/Misc depends on rsMatrix - todo: we want to run the tests in the order of
  // dependencies such that the more independent classes get tested first but we also want to have 
  // the tests ordered as in the directory structure - these goals are in conflict - resolve by 
  // changing directory structure. maybe an Algorithms folder could be created to hold the core
  // algorithms (which do not yet depend on any of our math-types)

  return testResult;
}

bool runTestsAudio(std::string &reportString)
{  
  bool testResult = true;

  reportString += "\nAudio-Analysis:\n";
  testResult &= testAutoCorrelationPitchDetector(reportString);

  reportString += "\nAudio-Filters:\n";
  testResult &= testFilterPolynomials(reportString);
  testResult &= testHighOrderFilter(  reportString);

  reportString += "\nMisc-Audio:\n";
  testResult &= testBandwidthConversions(reportString);
  testResult &= testSincInterpolation(   reportString);
  testResult &= testSineParameters(      reportString);
  testResult &= testZeroCrossingFinder(  reportString);

  reportString += "\nAudio-Synthesis:\n";
  testResult &= testModalFilter(reportString);
  testResult &= testModalSynth(reportString);

  return testResult;
}

bool runTestsGraphics(std::string &reportString)
{  
  bool testResult = true;

  reportString += "\nGraphics:\n";
  testResult &= testImage(reportString);
  testResult &= testFonts(reportString);
  testResult &= testRendering(reportString);

  return testResult;
}

void runTests()
{
  std::string reportString;
  bool testResult = true;
  
  testResult &= runTestsCore(reportString);
  testResult &= runTestsMath(reportString);
  testResult &= runTestsAudio(reportString);
  testResult &= runTestsGraphics(reportString);

  if( testResult == true )
    reportString += "\nAll Tests passed";
  else
    reportString += "\nAt least one test !!! FAILED !!!";

  std::cout << reportString;
}

int main(int argc, char** argv)
{
  runTests();

  //rsString *leakString = new rsString("LeakString"); 
    // uncomment this when you want to trigger the leak-detection deliberately (for testing the 
    // leak-detection itself)

  if( detectMemoryLeaks() )
    std::cout << "\n\n!!! Memory leaks detected !!! \n";

  getchar();
  return(EXIT_SUCCESS);
}
