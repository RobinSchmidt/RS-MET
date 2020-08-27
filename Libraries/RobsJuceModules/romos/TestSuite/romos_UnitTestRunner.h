#ifndef romos_UnitTestRunner_h
#define romos_UnitTestRunner_h

//#include "romos_GlobalFrameworkTests.h"
//#include "romos_ConcreteProcessingTests.h"
//#include "romos_ContainerManipulationTests.h"
//#include "romos_ConcreteModularSystemTests.h"

namespace rsTestRomos
{

/** This class serves as the high-level interface for running all the unit-tests for the romos 
library. To run all the tests, you simply create an insteance of this class and call 
runAllTestsAndPrintResultsToConsole() on it. You can also selectively run only certain tests by 
using the other runXxxxTests methods. */

class UnitTestRunner
{
public:

  /** Constructor.  */
  UnitTestRunner();

  /** This function runs all the test and informs whether or not all tests have passed. */
  virtual bool runAllTestsAndPrintResultsToConsole();


  /** Runs tests that check global things like voice-allocation, etc. */
  virtual bool runGlobalFrameworkTests();

  /** Runs tests that involve polyphonic (and possibly also monophonic) processing functions of modules. */
  virtual bool runProcessingTests();

  /** Runs the tests that involve manipulations of container modules such as moving around child-modules, containerizing/uncontainerizing,
  etc. */
  virtual bool runContainerManipulationTests();

  /** Runs the tests that involve the whole ModularSystem. */
  virtual bool runSystemTests();


  virtual bool runMiscTests();


  // separate into runTestsWithoutEvents, rundTestsWithEvents, runProcessingStatusTests, etc.

protected:

  /** Prints the result of the test to the console. You should pass the name of the test and the result (whether or not the test has
  passed. The function will them print a message like "testName passed.\n" or "!!! testName FAILED !!!\n". */
  virtual void printTestResultToConsole(bool testPassed, const char *testName);


};

} // end namespace romos

#endif 
