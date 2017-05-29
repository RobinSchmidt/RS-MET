#ifndef rosic_DelayPhaser_h
#define rosic_DelayPhaser_h

//// rosic-indcludes:
//#include "rosic_Phaser.h"
//#include "../delaylines/rosic_PingPongEcho.h"

namespace rosic
{

  /**

  This class combines two Phasers and a PingPongEcho into one effect unit. Moreover, it allows
  for feedback around phaser1->delay, delay->phaser2 and phaser1->delay->phaser2.

  */

  class DelayPhaser
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    DelayPhaser();

    /** Destructor */
    ~DelayPhaser();

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets up the tempo in  beats per minute. */
    void setTempoInBPM(double newTempoInBPM);

    /** Sets the amount of feedback (as raw factor) from the output of the delay to the input of 
    the first phaser. */
    void setFeedback1(double newFeedback) { fb1 = newFeedback; }

    /** Sets the amount of feedback (as raw factor) from the output of the second phaser 
    to the input of the delay. */
    void setFeedback2(double newFeedback) { fb2 = newFeedback; }

    /** Sets the amount of feedback (as raw factor) from the output of the second phaser 
    to the input of the first phaser. */
    void setFeedback3(double newFeedback) { fb3 = newFeedback; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one stereo sample frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state. */
    void reset();

    //---------------------------------------------------------------------------------------------
    // embedded objects:

    Phaser phaser1, phaser2;
    PingPongEcho delay;

    //=============================================================================================

  protected:

    // feedback factors:
    double fb1, fb2, fb3;

    // left and right output signals from delay and phaser2 (for feedback):
    double dL, dR, p2L, p2R;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void DelayPhaser::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double tmpL = *inOutL + fb1*dL + fb3*p2L;
    double tmpR = *inOutR + fb1*dR + fb3*p2R;

    phaser1.getSampleFrameStereo(&tmpL, &tmpR);
    tmpL += fb2*p2L;
    tmpR += fb2*p2R;

    delay.getSampleFrameStereo(&tmpL, &tmpR);
    dL = tmpL;
    dR = tmpR;

    phaser2.getSampleFrameStereo(&tmpL, &tmpR);
    p2L = tmpL;
    p2R = tmpR;

    *inOutL = tmpL;
    *inOutR = tmpR;
  }

} // end namespace rosic

#endif // #ifndef rosic_DelayPhaser_h
