#ifndef rosic_AudioToMidi_h
#define rosic_AudioToMidi_h

// rosic-indcludes:
#include "../analysis/rosic_PitchDetector.h"

namespace rosic
{

  /**

  This is a simple grain-based pitch shifter. It uses a delayline with two output tap pointers 
  which are spaced apart by half of the delayline length. The individual tap-outputs are multiplied 
  by a cos^2 shaped window-function (or grain-envelope) which has its maximum when the pointer 
  passes through the point of half the delayline length and its minima at the two 
  (forward/backward) wraparound points. These two output tap signals are then summed. For upward 
  pitch-shifting, the input can optionally be lowpass-filtered before entering the delayline in 
  order to avoid aliasing.

  */

  class AudioToMidi
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    AudioToMidi();
    ///< Constructor.

    ~AudioToMidi();
    ///< Destructor.

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    void setSampleRate(double newSampleRate);
    /**< Sets the sample-rate. */

    void setMinFrequency(double newMinFrequency);
    /**< Sets the minimum expected fundamental frequency in Hz. */

    void setMaxFrequency(double newMaxFrequency);
    /**< Sets the maximum expected fundamental frequency in Hz. */





    //---------------------------------------------------------------------------------------------
    // inquiry:

    double getMinFrequency();
    /**< Returns the minimum expected fundamental frequency in Hz. */

    double getMaxFrequency();
    /**< Returns the maximum expected fundamental frequency in Hz. */




    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double* inL, double* inR, double* outL, double* outR);
    /**< Calculates a stereo-ouput frame. */

    //---------------------------------------------------------------------------------------------
    // others:

    void reset();
    /**< Resets the content of the delaylines to all zeros. */

    //=============================================================================================

    PitchDetector pitchDetector;

  protected:

    double sampleRate;
    double minFreq, maxFreq;
    double pitchSmooth;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE void AudioToMidi::getSampleFrameStereo(double* inL,  double* inR,
                                                 double* outL, double* outR)
  { 

    // test:
    //*outL = *outR = pitchDetector.predictor.getSample(*inL + *inR);
    //*outL = *outR = pitchDetector.estimatePeriod(*inL + *inR);

    
    pitchDetector.estimatePeriod(*inL + *inR);

    // mix dry wet and store the result:
    *outL = *inL;
    *outR = *inR;
    
  }

} // end namespace rosic

#endif // rosic_AudioToMidi_h
