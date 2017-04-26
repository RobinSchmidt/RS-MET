#ifndef rosic_Noisifier_h
#define rosic_Noisifier_h

// rosic-indcludes:
#include "rosic_ModulationEffect.h"
#include "../generators/rosic_NoiseGenerator.h"
#include "../analysis/rosic_EnvelopeFollower.h"

namespace rosic
{

  /**

  This class implements an effect that adds colored noise to an incoming signal, thereby possibly 
  applying the signals amplitude envelope to the generated noise (partly or fully, controlled via
  setEnvelopeAmount().

  */

  class Noisifier
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Noisifier();   

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the amount of the input signal's envelope to be applied to the generated noise. */
    void setEnvelopeAmount(double newAmount) { envelopeAmount = newAmount; }

    /** Sets the attack time constant for the envelope follower in milliseconds. */
    void setAttackTime(double newAttack) { envFollower.setAttackTime(newAttack); }

    /** Sets the release time constant for the envelope follower in milliseconds. */
    void setReleaseTime(double newRelease) { envFollower.setReleaseTime(newRelease); }

    /** Sets the level with which the input signal is passed through. */
    void setPassThroughLevel(double newLevel) { passGain = dB2amp(newLevel); }

    /** Sets the level with which the generated noise is mixed in. */
    void setNoiseLevel(double newLevel) { noiseGain = dB2amp(newLevel); }

    /** sets the spectral slope for the noise generator in dB/oct. */
    void setNoiseSpectralSlope(double newSlope) { noiseGenerator.setSpectralSlope(newSlope); }
        
    /** Sets the lowest frequency in the noise to be generated (in Hz) */
    void setLowestFrequency(double newFreq) { noiseGenerator.setLowestFrequency(newFreq); }

    /** Sets the highest frequency in the noise to be generated (in Hz) */
    void setHighestFrequency(double newFreq) { noiseGenerator.setHighestFrequency(newFreq); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the effect. */
    void reset();

    //=============================================================================================

  protected:

    // embedded modules:
    NoiseGenerator   noiseGenerator;
    EnvelopeFollower envFollower; 

    double envelopeAmount;
    double passGain, noiseGain;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void Noisifier::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double env         = envFollower.getSample(*inOutL + *inOutR);
    double noiseFactor = noiseGain * (1.0 - envelopeAmount + envelopeAmount*env);
    double noise       = noiseGenerator.getSampleThreadSafe();   

    *inOutL = passGain*(*inOutL) + noiseFactor*noise;
    *inOutR = passGain*(*inOutR) + noiseFactor*noise;
  }

} // end namespace rosic

#endif // #ifndef rosic_Noisifier_h
