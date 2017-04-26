#ifndef rosic_IntegerDelayLine_h
#define rosic_IntegerDelayLine_h

// rosic-indcludes:
#include "rosic_BasicIntegerDelayLine.h"

namespace rosic
{

  /**

  In addition to its baseclass BasicIntegerDelayLine, this class keeps information about the 
  sample-rate and the delay in seconds.

  \todo: facilitate temp-sync by maintaining a bpm-value and a sync-flag 
  ...maybe do this in a subclass...

  */

  class IntegerDelayLine : public BasicIntegerDelayLine
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a delay-line with a given maximum number of samples delay. This
    has to be a power of two minus 1 - otherwise the next power of two minus 1 will be used. */
    IntegerDelayLine(int maximumDelayInSamples = 65536);

    /** Destructor */
    ~IntegerDelayLine();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the delay-time in samples. */
    int setDelayInSamples(int newDelayInSamples);

    /** Sets the delay-time in seconds. */
    void setDelayInSeconds(double newDelayInSeconds);

    /** Sets the delay-time in milliseconds. */
    void setDelayInMilliseconds(double newDelayInMilliseconds);

    //---------------------------------------------------------------------------------------------
    // inquiry (get-, is-, etc. functions):

    /** Returns the delay-time in seconds. */
    double getDelayInSeconds() const { return delayInSeconds; }

    /** Returns the delay-time in milliseconds. */
    double getDelayInMilliseconds() const { return 1000.0 * delayInSeconds; }

    //=============================================================================================

  protected:

    double delayInSeconds;
    double sampleRate;

  };

} // end namespace rosic

#endif // #ifndef rosic_IntegerDelayLine_h
