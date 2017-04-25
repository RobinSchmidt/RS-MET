#ifndef rosic_ChannelMatrix2x2_h
#define rosic_ChannelMatrix2x2_h

// rosic-indcludes:
#include "GlobalFunctions.h"

namespace rosic
{

  /**

  This class implements a matrix mutliplication of channels from two signals to yield another pair
  of two signals according to yL = gLL*xL + gRL*xR, yR = gLR*xL + gRR*xR.

  */

  class ChannelMatrix2x2
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ChannelMatrix2x2();   

    /** Destructor. */
    ~ChannelMatrix2x2();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the amount by which the left input goes to the right output. */
    void setLeftToLeftGain(double newGain) { gLL = newGain; };

    /** Sets the amount by which the right input goes to the left output. */
    void setRightToLeftGain(double newGain) { gRL = newGain; }

    /** Sets the amount by which the left input goes to the right output. */
    void setLeftToRightGain(double newGain) { gLR = newGain; }

    /** Sets the amount by which the right input goes to the right output. */
    void setRightToRightGain(double newGain) { gRR = newGain; }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the amount by which the left input goes to the right output. */
    double getLeftToLeftGain() { return gLL; }

    /** Returns the amount by which the right input goes to the left output. */
    double getRightToLeftGain() { return gRL; }

    /** Returns the amount by which the left input goes to the right output. */
    double getLeftToRightGain() { return gLR; }

    /** Returns the amount by which the right input goes to the right output. */
    double getRightToRightGain() { return gRR; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a stereo-ouput frame. */
    INLINE void getSampleFrameStereo(double* inOutL, double* inOutR);

    //=============================================================================================

  protected:

    double gLL, gRL, gLR, gRR;

  };

  //-----------------------------------------------------------------------------------------------
  // inlined member-functions:

  INLINE void ChannelMatrix2x2::getSampleFrameStereo(double* inOutL,  double* inOutR)
  {
    double tmpL  = *inOutL;
    double tmpR  = *inOutR;
    *inOutL      = gLL*tmpL + gRL*tmpR;
    *inOutR      = gLR*tmpL + gRR*tmpR;
  }

} // end namespace rosic

#endif // rosic_ChannelMatrix2x2_h
