#ifndef romos_PerformanceTest_h
#define romos_PerformanceTest_h

//#include "romos_ProcessingTest.h"

namespace rsTestRomos
{

  /**

  Baseclass for all performance measurement tests.

  */

  class PerformanceTest : public UnitTest
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    PerformanceTest(const char *testName);

    /** Destructor. */
    virtual ~PerformanceTest();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Sets up all the members for either polyphonic or monophonic processing. */
    virtual void setTestPolyphonic(bool shouldBePolyphonic);



    //-------------------------------------------------------------------------------------------------------------------------------------
    // information output

    /** Prints the result of the test to the console. */
    //virtual void printResultToConsole();



    virtual rosic::rsString runTestsAndGetReport();


    // data for formatting the report string:
    static const int  nameLength   = 32;
    static const int  numberLength = 12;
    static const char padCharacter = ' ';


    //=====================================================================================================================================

  protected:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // running the test:

    /** This function runs the tests and returns true when the tests. the return value can be ignored as we are not concerned about 
    pass/fail here.  ...nahh not needed - but we must implement it anyway because it's purely virtual in the baseclass */
    virtual bool runTest();

    virtual rosic::rsString runFrameWiseTestAndGetReport();



    virtual rosic::rsString runBlockWiseTestAndGetReport();




    //-------------------------------------------------------------------------------------------------------------------------------------
    // data:

    romos::Module *moduleToTest;

    double randomValues[10000];

    int numFramesPerRun, numRuns;
    
    rosic::ProcessorCycleCounter counter;



  };

} // end namespace romos

#endif 
