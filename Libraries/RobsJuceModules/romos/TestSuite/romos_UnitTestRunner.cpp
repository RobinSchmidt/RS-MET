#include "romos_UnitTestRunner.h"
//using namespace rsTestRomos;

namespace rsTestRomos
{

UnitTestRunner::UnitTestRunner()
{

}

bool UnitTestRunner::runAllTestsAndPrintResultsToConsole()
{
  bool result = true;

  result &= runMiscTests();
  result &= runGlobalFrameworkTests();
  result &= runProcessingTests();  // causes memleak, fails with gcc
  result &= runContainerManipulationTests();
  result &= runSystemTests();

  return result;
}

bool UnitTestRunner::runGlobalFrameworkTests()
{
  const char* testName = "GlobalFrameworkTests";
  printf("%s %s", testName, ":\n");

  UnitTest* test;
  bool testsPassed = true;

  test = new TopLevelModuleTest();  testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new VoiceAllocatorTest();  testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  //test = new TriggerAndKillTest();  testsPassed &= test->runTestAndPrintResultToConsole(); delete test;

  printTestResultToConsole(testsPassed, testName);
  return testsPassed;
}


bool UnitTestRunner::runProcessingTests()
{
  const char* testName = "ProcessingTests";
  printf("%s %s", testName, ":\n");

  UnitTest* test;
  bool testsPassed = true;

  test = new Formula_N_1Test();                testsPassed &= test->runTestAndPrintResultToConsole(); delete test;

  test = new IdentityTest();                   testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new AdderTest();                      testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new Adder3Test();                     testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new Adder4Test();                     testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new Adder5Test();                     testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new SubtractorTest();                 testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new UnitDelayTest();                  testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new NoiseGeneratorTest();             testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new WrappedAdderTest();               testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new SumDiffProdTest();                testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new WrappedSumDiffProdTest();         testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new WrappedAdderNTest();              testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new SummedDiffsTest();                testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new MovingAverageTest();              testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new DelayedConnectionTest();          testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new LeakyIntegratorTest();            testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new LeakyIntegratorDoubleDelayTest(); testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new TestFilter1Test();                testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new BiquadMacroTest();                testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new BiquadAtomicTest();               testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new BiquadFormulaTest();              testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new BlipTest();                       testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new MonoToPolyTest();                 testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new VoiceCombinerTest();              testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new Formula1In1OutTest();             testsPassed &= test->runTestAndPrintResultToConsole(); delete test;

  //test = new PolyBlipStereoTest();             testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  //test = new GateAndKillTest();  testsPassed &= test->runTestAndPrintResultToConsole(); delete test;

  printTestResultToConsole(testsPassed, testName);
  return testsPassed;
}

bool UnitTestRunner::runContainerManipulationTests()
{
  const char* testName = "ContainerManipulationTests";
  printf("%s %s", testName, ":\n");

  UnitTest* test;
  bool testsPassed = true;


  test = new Containerize01();                     testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new Containerize02();                     testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new OutputModuleDeletion();               testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new ContainerizationAddedConstantsTest(); testsPassed &= test->runTestAndPrintResultToConsole(); delete test; // fails
  test = new PinSortingTest();                     testsPassed &= test->runTestAndPrintResultToConsole(); delete test;


  printTestResultToConsole(testsPassed, testName);
  return testsPassed;
}


bool UnitTestRunner::runSystemTests()
{
  const char* testName = "SystemTests";
  printf("%s %s", testName, ":\n");

  UnitTest* test;
  bool testsPassed = true;

  test = new BypassTest();          testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new BypassWithChildTest(); testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new PolyBlipStereoTest();  testsPassed &= test->runTestAndPrintResultToConsole(); delete test;
  test = new NoiseFluteTest();      testsPassed &= test->runTestAndPrintResultToConsole(); delete test;

  printTestResultToConsole(testsPassed, testName);
  return testsPassed;
}

bool UnitTestRunner::runMiscTests()
{
  bool testsPassed = true;
  const char* testName = "MiscTests";
  printf("%s %s", testName, ":\n");

  testsPassed &= testFormulaModules();


  printTestResultToConsole(testsPassed, testName);
  return testsPassed;
}


void UnitTestRunner::printTestResultToConsole(bool testPassed, const char* testName)
{
  if(testPassed)
    printf("%s %s", testName, " passed.\n\n");
  else
    printf("%s %s %s", "!!! ", testName, " FAILED !!!\n\n");
}


}