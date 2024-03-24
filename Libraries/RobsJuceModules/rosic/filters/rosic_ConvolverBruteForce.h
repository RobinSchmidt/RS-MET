#ifndef rosic_ConvolverBruteForce_h
#define rosic_ConvolverBruteForce_h

//// standard includes
//#include <string.h> // for memmove
//
//// rosic-indcludes:
//#include "../basics/GlobalFunctions.h"
//#include "../infrastructure/rosic_MutexLock.h"

namespace rosic
{

  /**

  This class implements a convolution via directly implementing the convolution sum:

  \f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] + ... + b[M-1] x[M-1] \f]

  where M is the length of the impulse response.

  */

  class ConvolverBruteForce
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ConvolverBruteForce();

    /** Destructor. */
    ~ConvolverBruteForce();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Sets up the impulse response to be used. */
    void setImpulseResponse(double *newImpulseResponse, int newLength);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a single convolved output-sample. */
    INLINE double getSample(double in);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    /** Sets the impulse response to an impulse (1 0 0 0....). */
    void clearImpulseResponse();

    /** Sets the buffer for the previous input samples to zero. */
    void clearInputBuffer();

    //=====================================================================================================================================

  protected:

    doubleA *x;       // buffer for past input samples
    doubleA *h;       // array containing the impulse response
    int M;            // length of the impulse response
    MutexLock mutex;  // mutex-lock for accessing the buffers in a thread safe way - 
                      // \todo get rid of that, thread-safety should be handled elsewhere

  };

  //---------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double ConvolverBruteForce::getSample(double in)
  {
    if( x == NULL || h == NULL )
    {
      DEBUG_BREAK;
      return 0.0;
    }

    mutex.lock();

    // shift all the values in the vector which contains the past inputs down by one position and
    // write the new input sample into the slot with index 0:
    memmove(&x[1], &x[0], (M-1)*sizeof(double)); // same as: for(i=(M-1); i>0; i--)  x[i] = x[i-1];
    x[0] = in;                                   // x[k] represents the input delayed by k samples

    // do the sum:
    doubleA tmp = 0.0;
    for(int k=0; k<M; k++)
      tmp += h[k] * x[k];

    mutex.unlock();

    return tmp;
  }
  // ToDo: implement this using a circular buffer - get rid of the memmove - that's wasteful. But 
  // measure performance! Maybe keep both implementations

} 

#endif 
