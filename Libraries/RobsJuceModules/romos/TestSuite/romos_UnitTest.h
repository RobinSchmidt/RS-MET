#ifndef romos_UnitTest_h
#define romos_UnitTest_h

//#include "../framework/romos_ContainerModule.h"

namespace rsTestRomos
{

  /**

  This class serves as baseclass for all unit-tests.

  */

  class UnitTest
  {
  public:

    /** Constructor. You should pass the name of the test case to the constructor - this name will be used in report strings that will be
    retrieved from the object. */
    UnitTest(const char *testName);

    /** Destructor. */
    virtual ~UnitTest();

    /** This function runs the test and informs whether or not the test has passed. */
    virtual bool runTestAndPrintResultToConsole();

    //---------------------------------------------------------------------------------------------
    // information output:

    /** Prints the structure of the moduleToTest member to the console. */
    virtual void printModuleStructure(romos::Module *module, int indent = 0);

  protected:

    /** This function must be overriden by concrete UnitTest subclasses and return true when the 
    test has passed and false when it failed. */
    virtual bool runTest() = 0;

    /** A hook function that is called after a test has been performed with the result of the 
    test. Subclasses may override this, for example to output debug information in case of a failed 
    test. */
    virtual void handleTestResult(bool didTestPass)
    {
      if(!didTestPass)
        RS_DEBUG_BREAK;
    }

    char *name;

  };

} // end namespace romos

#endif 
