#ifndef rosic_ProcessorCycleCounter_h
#define rosic_ProcessorCycleCounter_h

#ifdef _MSC_VER
#include <intrin.h> // Or #include <ia32intrin.h> etc.
#endif

// rosic-indcludes:
#include "../basics/GlobalDefinitions.h"

namespace rosic
{

  /**

  This class implements a CPU cycle counter which is useful for performance tests. It does not 
  literally count but instead measure the time stamps before and after a sequence of commands. 
  Works only with MSVC at the moment.

  \todo rename into something more appropriate

  */

  class ProcessorCycleCounter  
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ProcessorCycleCounter();

    /** Destructor */
    ~ProcessorCycleCounter();

    //---------------------------------------------------------------------------------------------
    // initialization and readout:

    /** Resets the counter to zero. */
    INLINE void init();

    /** Returns the number of CPU-cycles since the last call to init(). */
    INLINE INT64 getNumCyclesSinceInit();

    //=============================================================================================

  protected:

    INT64 initTime;    // remembers the time at which the call to init() occured 

    /** Function for reading the time-stamp counter - this is taken from Agner Fog's optimization 
    tutorials. */
    INT64 ReadTSC()  // Returns time stamp counter
    {
#ifdef _MSC_VER
      int   dummy[4];      // For unused returns
      INT64 clock;         // Time 
      __cpuid(dummy, 0);   // Serialize
      clock = __rdtsc();   // Read time
      __cpuid(dummy, 0);   // Serialize again
      return clock;
#else
			return 0;
#endif
    }

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void ProcessorCycleCounter::init()
  {
    initTime = ReadTSC();
  }

  INLINE INT64 ProcessorCycleCounter::getNumCyclesSinceInit()
  {
    return ReadTSC() - initTime;
  }

} 

#endif 
