#ifndef romos_InteractiveTestRunner_h
#define romos_InteractiveTestRunner_h

//#include "romos_InteractivePlottingTests.h"

namespace rsTestRomos
{

  /**

  This class serves as the high-level interface for running all the interactive tests for the romos library.

  */

  class InteractiveTestRunner
  {
  public:

    /** Constructor.  */
    InteractiveTestRunner();



    /** Runs tests and outputs the result. */
    virtual void runTests();

  protected:





  };

} // end namespace romos

#endif 
