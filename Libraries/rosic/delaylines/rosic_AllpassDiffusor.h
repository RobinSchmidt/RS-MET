#ifndef rosic_AllpassDiffusor_h
#define rosic_AllpassDiffusor_h

// rosic-indcludes:
#include "rosic_IntegerDelayLine.h"

namespace rosic
{

  /**

  This class implements an allpass diffusor based on a delay-line with an integer number of 
  samples delay (it is derived from class IntegerDelayLine). It introduces a new parameter - namely
  the allpass gain g which is related to the amount of diffusion (between -100...+100 percent) via
  g = amount/100 where g is used as in the DAFX book, page 177. If the amount (and therefore g) is 
  zero, the allpass reduces to a plain delay (as in IntegerDelayLine). 

  */

  class AllpassDiffusor : public IntegerDelayLine
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a delay-line with a given maximum number of samples delay. */
    AllpassDiffusor(int maximumDelayInSamples = 65536);

    /** Destructor */
    ~AllpassDiffusor();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the amount of diffusion in percent. It can take on values between -100.0...+100.0 and
    will map to the allpass gain like g = amount/100 where g is used as in the DAFX book, page 177. 
    If the amount (and therefore g) is zero, the allpass reduces to a plain delay. */
    void setDiffusionAmount(double newAmount);

    //---------------------------------------------------------------------------------------------
    // inquiry (get-, is-, etc. functions):

    /** Returns the amount of diffusion in percent. */
    double getDiffusionAmount(); 

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Overrides the inherited method from the IntegerDelayLine base class. */
    INLINE double getSample(double in);

    /** Overrides the inherited method from the IntegerDelayLine base class. */
    INLINE double getSampleSuppressTapIncrements(double in);

    //=============================================================================================

  protected:

    double g;
      // The allpass gain coefficient (using the notation from DAFX) - if it is zero, the allpass
      // reduces to a plain delay. */

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double AllpassDiffusor::getSample(double in)
  {
    incrementTapPointers();
    return getSampleSuppressTapIncrements(in);
  }

  INLINE double AllpassDiffusor::getSampleSuppressTapIncrements(double in)
  {
    // get the output-sample: 
    double out = -g*in + delayLine[tapOut];

    // write the incoming sample plus the feedback-part into the delay-line:
    delayLine[tapIn] = in + g*out;

    return out;
  }

} // end namespace rosic

#endif // #ifndef rosic_AllpassDiffusor_h
