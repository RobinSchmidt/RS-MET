#ifndef romos_PerformanceTestRunner_h
#define romos_PerformanceTestRunner_h

//#include "romos_ConcretePerformanceTests.h"

namespace rsTestRomos
{

  /**

  This class serves as the high-level interface for running all the performance-tests for the romos library.

  */

  class PerformanceTestRunner
  {
  public:

    /** Constructor.  */
    PerformanceTestRunner();

    virtual void runAllTestsAndPrintResultsToConsole(bool createLogFile);




    /** Runs tests that check global the efficiency of the gloabl framework using interconnected modules inside containers where the 
    modules themselves do not do much proceesing - i.e. we use simple modules like adders, identities, etc. These tests can be used to 
    optimize global framwork things like the way that connections are computed, buffering, memory access etc. */
    virtual rosic::rsString runFrameworkPerformanceTests();

    /** Runs tests that check the efficiency of individual atomic modules. */
    virtual rosic::rsString runAtomicModulePerformanceTests();


    /** runs performance tests for simple standard functions like sin, cos tan, exp and self-defined functions liek rMax, sign, etc. */
    virtual rosic::rsString runInternalFunctionPerformanceTests();

  protected:





  };

  double dummyFunction(double x); 

} // end namespace romos

#endif 
