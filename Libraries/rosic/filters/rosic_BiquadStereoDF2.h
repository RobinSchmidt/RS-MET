#ifndef rosic_BiquadStereoDF2_h
#define rosic_BiquadStereoDF2_h

// rosic-indcludes:
#include "rosic_BiquadBase.h"

namespace rosic
{

  /**

  This class implements a biquad filter which realize the pair of difference equations:

  \f[ w[n] = x[n] + a_1 y[n-1] + a_2 y[n-2]     \f]
  \f[ y[n] = b_0 w[n] + b_1 w[n-1] + b_2 w[n-2] \f]

  note the positive sign of the feedback part, as opposed to negative sign seen in many DSP textbooks - this has been chosen to allow for 
  better potential for parallelization.

  ...performance tests indicate that DF1 seems to be more efficient

  */

  class BiquadStereoDF2 : public BiquadBase
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    BiquadStereoDF2();

    /** Destructor. */
    ~BiquadStereoDF2();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output sample at a time. */
    INLINE double getSample(double in);

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    /** Sets the buffers for the previous input and output samples to zero. */
    void reset();

    //=====================================================================================================================================

  protected:

    // internal state buffers:
    double lw1, lw2, rw1, rw2;

  };

  //---------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double BiquadStereoDF2::getSample(double in)
  {
    double w, y;

    w   = (in+TINY) + a1*lw1 + a2*lw2;
    y   = b0*w      + b1*lw1 + b2*lw2;
    lw2 = lw1;
    lw1 = w;

    return y;
  }

  INLINE void BiquadStereoDF2::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double w;

    // left:
    w       = (*inOutL+TINY) + a1*lw1 + a2*lw2;
    *inOutL = b0*w           + b1*lw1 + b2*lw2;
    lw2     = lw1;
    lw1     = w;

    // right:
    w       = (*inOutR+TINY) + a1*rw1 + a2*rw2;
    *inOutR = b0*w           + b1*rw1 + b2*rw2;
    rw2     = rw1;
    rw1     = w;
  }

} 

#endif 
