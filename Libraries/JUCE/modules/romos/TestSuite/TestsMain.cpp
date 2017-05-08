//#include "AutomaticTests.h"
//#include "InteractiveTests.h"
//#include "PerformanceTests.h"

#include "../romos.h"

bool runUnitTests()
{
  romos::UnitTestRunner testRunner;
  bool testsPassed = testRunner.runAllTestsAndPrintResultsToConsole();


  if( testsPassed == false )
    printf("%s", "At least one unit test FAILED !!! FAILED !!! FAILED !!! FAILED !!! FAILED !!!\n");
  else
    printf("%s", "All unit tests passed\n");
  printf("%s", "\n");


  return testsPassed;
}

void runPerformanceTests(bool createLogFile)
{
  romos::PerformanceTestRunner testRunner;
  testRunner.runAllTestsAndPrintResultsToConsole(createLogFile);
}

/*
void runAutomatedCorrectnessTests()
{
  bool verboseOutput    = false;
  int  numVoicesToCheck = 3;

  runUnitTests();
  //runPerformanceTests(false);

  //double *leak = new double[100];
  //String *leakString = new String("LeakString"); // uncomment this when you want to trigger the leak-detection deliberately
}
*/



void testCodeGenerator()
{
  // here we choose what kind of module to create (and genreate code for):
  //romos::Module *testModule = createSumDiffModule(NULL);
  //romos::Module *testModule = createWrappedSumDiffModule(NULL);
  romos::Module *testModule = romos::TestModuleBuilder::createTestFilter1("TestFilter1", 0, 0, false);

  rosic::String codeForModule = romos::ModuleBuildCodeGenerator::getCodeForModule(testModule);
  codeForModule.printToStandardOutput();
  romos::ModuleFactory::deleteModule(testModule);
}

void runInteractiveTests()
{
  romos::InteractiveTestRunner testRunner;
  testRunner.runTests();


  //testCodeGenerator();  // move this into the testrunner, create a class for the codegenerator test etc.
}

/*
void runPerformanceTests()
{
  printModuleByteSizes();
  //printf("%s %.3f %s", "Performance hit for biquad: ", measureFrameworkOverheadWithBiquad(), "\n" );
}
*/

void runTests()
{
  /*
  // move this into a test-class of it own - maybe we need some reliable, invertible and pretty 
  // double->string->double roundtrip implementations in romos
  rosic::String s01 = rosic::String(6467325487632587632.0234);
  rosic::String s02 = rosic::String(0.2);
  rosic::String s03 = rosic::String(0.7);
  rosic::String s04 = rosic::String(0.9);
  rosic::String s05 = rosic::String(0.0002);
  rosic::String s06 = rosic::String(0.9999);
  rosic::String s07 = rosic::String(1.2);
  rosic::String s08 = rosic::String(1.9);
  rosic::String s09 = rosic::String(10.2);
  rosic::String s10 = rosic::String(10.9);
  rosic::String s11 = rosic::String(10.79);
  rosic::String s12 = rosic::String(10.99);
  rosic::String s13 = rosic::String(11.99);
  rosic::String s14 = rosic::String(117.99);
  rosic::String s15 = rosic::String(9.9999999999999999);
  */

  double someDouble = -0.4;
  double absOfSomeDouble = fabs(someDouble);
  double abs2 = rAbs(someDouble);
  //unsigned long long intAbsValue = *((unsigned long long*) &someDouble) & 0x7FFFFFFFFFFFFFFFULL;
  //abs2 = *((double*) &intAbsValue);


  //runUnitTests();
  runInteractiveTests();
  //runPerformanceTests(false);

  romos::ModuleTypeRegistry::deleteSoleInstance(); // deletes the singleton object
  romos::BlitIntegratorInitialStates::deleteStateValueTables();
}





int main(int argc, char** argv)
{


  runTests();

  if( detectMemoryLeaks() )
    printf("%s", "Memory leaks detected\n");

  printf("%s", "Tests done.");
  //getchar();
  return 0;
}
