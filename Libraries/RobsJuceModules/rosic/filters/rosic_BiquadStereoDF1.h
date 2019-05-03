#ifndef rosic_BiquadStereoDF1_h
#define rosic_BiquadStereoDF1_h

//// rosic-indcludes:
//#include "rosic_BiquadBase.h"

namespace rosic
{

  /**

  This class implements a biquad filter which realize the difference equation:

  \f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] + a_1 y[n-1] + a_2 y[n-2] \f]

  note the positive sign of the feedback part, as opposed to negative sign seen in many DSP
  textbooks - this has been chosen to allow for better potential for parallelization.

  */

  class BiquadStereoDF1 : public BiquadBase
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    BiquadStereoDF1();

    /** Destructor. */
    ~BiquadStereoDF1();

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a single filtered output-sample in place. */
    INLINE void getSampleInPlace(double &inOut);

    /** Calculates a single filtered output-sample. */
    INLINE double getSample(double in);

    /** Calculates a single filtered output-sample. */
    INLINE float getSample(float in);

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(float *inOutL, float *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Sets the buffers for the previous input and output samples to zero. */
    void reset();

    //=============================================================================================

  protected:

    // internal state buffers:
    double lx1, lx2, ly1, ly2, rx1, rx2, ry1, ry2;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void BiquadStereoDF1::getSampleInPlace(double &inOut)
  {
    doubleA tmp = inOut;

    tmp = (b0*tmp+TINY) + (b1*lx1 + b2*lx2) + (a1*ly1 + a2*ly2); // parentheses facilitate out-of-order execution

    lx2 = lx1;
    lx1 = inOut;
    ly2 = ly1;
    ly1 = tmp;

    inOut = tmp;
  }

  INLINE double BiquadStereoDF1::getSample(double in)
  {
    doubleA tmp;

    tmp = (b0*in+TINY) + (b1*lx1 + b2*lx2) + (a1*ly1 + a2*ly2);

    lx2 = lx1;
    lx1 = in;
    ly2 = ly1;
    ly1 = tmp;

    return tmp;
  }

  INLINE float BiquadStereoDF1::getSample(float in)
  {
    return (float) getSample((double) in);
  }

  INLINE void BiquadStereoDF1::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    doubleA tmp;

    // left:
    tmp     = (b0*(*inOutL)+TINY) + (b1*lx1 + b2*lx2) + (a1*ly1 + a2*ly2);
    lx2     = lx1;
    lx1     = *inOutL;
    ly2     = ly1;
    ly1     = tmp;
    *inOutL = tmp;

    // right:
    tmp     = (b0*(*inOutR)+TINY) + (b1*rx1 + b2*rx2) + (a1*ry1 + a2*ry2);
    rx2     = rx1;
    rx1     = *inOutR;
    ry2     = ry1;
    ry1     = tmp;
    *inOutR = tmp;
  }

  INLINE void BiquadStereoDF1::getSampleFrameStereo(float *inOutL, float *inOutR)
  {
    doubleA tmp;

    // left:
    tmp     = (b0*(*inOutL)+TINY) + (b1*lx1 + b2*lx2) + (a1*ly1 + a2*ly2);
    lx2     = lx1;
    lx1     = *inOutL;
    ly2     = ly1;
    ly1     = tmp;
    *inOutL = (float) tmp;

    // right:
    tmp     = (b0*(*inOutR)+TINY) + (b1*rx1 + b2*rx2) + (a1*ry1 + a2*ry2);
    rx2     = rx1;
    rx1     = *inOutR;
    ry2     = ry1;
    ry1     = tmp;
    *inOutR = (float) tmp;
  }

} // end namespace rosic

#endif // rosic_BiquadStereoDF1_h
