#ifndef rosic_BiquadMonoDF1_h
#define rosic_BiquadMonoDF1_h

// rosic-indcludes:
#include "rosic_BiquadBase.h"

namespace rosic
{

  /**

  This class implements a biquad filter which realize the difference equation:

  \f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] + a_1 y[n-1] + a_2 y[n-2] \f]

  note the positive sign of the feedback part, as opposed to negative sign seen in many DSP
  textbooks - this has been chosen to allow for better potential for parallelization.

  */

  class BiquadMonoDF1 : public BiquadBase
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    BiquadMonoDF1();

    /** Destructor. */
    ~BiquadMonoDF1();

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a single filtered output-sample in place. */
    INLINE void getSampleInPlace(double &inOut);

    /** Calculates a single filtered output-sample. */
    INLINE double getSample(double in);

    /** Calculates a single filtered output-sample. */
    INLINE float getSample(float in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Sets the buffers for the previous input and output samples to zero. */
    void reset();

    //=============================================================================================

  protected:

    // internal state buffers:
    double x1, x2, y1, y2;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void BiquadMonoDF1::getSampleInPlace(double &inOut)
  {
    doubleA tmp = inOut;

    tmp = (b0*tmp+TINY) + (b1*x1 + b2*x2) + (a1*y1 + a2*y2);
      // parentheses facilitate out-of-order execution

    x2 = x1;
    x1 = inOut;
    y2 = y1;
    y1 = tmp;

    inOut = tmp;
  }

  INLINE double BiquadMonoDF1::getSample(double in)
  {
    doubleA tmp;

    tmp = (b0*in+TINY) + (b1*x1 + b2*x2) + (a1*y1 + a2*y2);
      // parentheses facilitate out-of-order execution

    x2 = x1;
    x1 = in;
    y2 = y1;
    y1 = tmp;

    return tmp;
  }

  INLINE float BiquadMonoDF1::getSample(float in)
  {
    return (float) getSample((double) in);
  }

} // end namespace rosic

#endif // rosic_BiquadMonoDF1_h
