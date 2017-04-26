#ifndef rosic_Limiter_h
#define rosic_Limiter_h

// rosic-indcludes:
#include "rosic_DynamicsProcessorBase.h"

namespace rosic
{

  /**

  This class implements a brickwall-limiter.

  */

  class Limiter : public DynamicsProcessorBase
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a dynamics processor with a given maximum number of samples 
    lookahead. */
    Limiter(int newLookAheadBufferSize = DEFAULT_LOOKAHEAD_SIZE);

    /** Destructor */
    ~Limiter();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the limit which shall not be exceeded (in decibels). */
    void setLimit(double newLimit) { limit = dB2amp(newLimit); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //=============================================================================================

  protected:

    double limit;  // limit as raw amplitude value

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void Limiter::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    // establish the mono sum of the (non-delayed) input signal and obtain its envelope:
    double tmpL = inputGainFactor * (*inOutL);
    double tmpR = inputGainFactor * (*inOutR);  
    double inM  = rmax(fabs(tmpL), fabs(tmpR));
    double env  = levelDetector.attackReleaseFollower.getSample(inM);

    double limiterGain;
    if( env > limit )
      limiterGain = limit / env;
    else
      limiterGain = 1.0;
    limiterGain *= outputGainFactor;

    // apply delay, gain an mix dry/wet (this can be almalgameted into a single operation):
    applyLookAheadDelay(&tmpL, &tmpR, &tmpL, &tmpR);
    tmpL    *= (dry + wet*limiterGain);
    tmpR    *= (dry + wet*limiterGain);
    *inOutL  = tmpL; 
    *inOutR  = tmpR; 
  }

} // end namespace rosic

#endif // #ifndef rosic_Limiter_h
