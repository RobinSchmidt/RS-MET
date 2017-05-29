#ifndef rosic_AlgoVerb_h
#define rosic_AlgoVerb_h

//// rosic-indcludes:
//#include "rosic_FeedbackDelayNetwork16.h"

namespace rosic
{

  /**

  This class implements an algorithmic reverb using a tapped delayline for creating early 
  reflections and a feedback delay network for creating the late reverb.

  ToDo: add input diffusors which smear the input signal into the FDN over some amount of time
  (should be within the integration time of human hearing, that is within 50 ms or so)
  -a diffusion parameter may simultanously control the length and the feedback of the allpass so 
   as to keep the smearing time constant
  -replace the delaylines in the FDN with notchpass filters
  -add modulation to the feedback matrix by applying rotation-matrices to any pair of adjacent 
   elements in the signal vector (maybe 2 passes, 1st pass starting at index 0, 2nd at 1)
   modulating signal is a sinusoid wich serves as angle-input to the sin/cos function 
   - angle signal can be scaled ('modulation-amount') and synced ('cycle-length' in beats)


  */

  class AlgoVerb
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    AlgoVerb();  

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the sample-rate. */
    void setSampleRate(float newSampleRate);

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(float newDryWet) 
    { equalPowerGainFactors(newDryWet, &dryVol, &wetVol, 0.0, 1.0); }

    /** Sets the level of the dry signal (in dB). */
    //void setDryLevel(double newLevel){ dryVol = (float) dB2amp(newLevel); }

    /** Sets the level of the early reflections (in dB). */
    void setEarlyReflectionLevel(double newLevel) { earlyVol = (float) dB2amp(newLevel); }

    /** Sets the level of the late reverb (in dB). */
    void setLateReverbLevel(double newLevel){ lateVol = (float) dB2amp(newLevel); }


    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Feeds in an impulse to the  reverb unit for auditioning the impulse response. */
    void feedInImpulse(float amplitude = 1.f) { impulseAmplitude = amplitude; }

    /** Generates one ouput sample at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state (content of the delaylines, etc.). */
    void reset();

    //---------------------------------------------------------------------------------------------
    // embedded modules

    FeedbackDelayNetwork16 fdn;

    //=============================================================================================

  protected:

    float impulseAmplitude; // if nonzero we feed in an impulse in the next call to getSample

    double dryVol, wetVol, earlyVol, lateVol;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void AlgoVerb::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    // possibly add the impulse (if nonzero) for auditioning purposes:
    double lateL = *inOutL + impulseAmplitude;
    double lateR = *inOutR + impulseAmplitude;
    impulseAmplitude = 0.f;

    fdn.getSampleFrameStereo(&lateL, &lateR);

    double earlyL = 0.0;
    double earlyR = 0.0;
     // todo: insert the early reflection generator here





    *inOutL = dryVol*(*inOutL) + wetVol * (lateVol*lateL + earlyVol*earlyL);
    *inOutR = dryVol*(*inOutR) + wetVol * (lateVol*lateR + earlyVol*earlyR);
  }

} // end namespace rosic

#endif // #ifndef rosic_AlgoVerb_h
