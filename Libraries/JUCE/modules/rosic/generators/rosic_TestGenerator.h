#ifndef rosic_TestGenerator_h
#define rosic_TestGenerator_h

// rosic-indcludes:
#include "rosic_SineOscillatorStereo.h"

namespace rosic
{

  /**

  This is a digital TestGenerator/distortion unit with sample rate decimation and 
  re-quantization (the latter is better known as bitcrushing).

  */

  class TestGenerator
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    TestGenerator();   ///< Constructor.
    ~TestGenerator();  ///< Destructor.

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    void setSampleRate(double newSampleRate);
    /**< Sets the sample-rate. */

    void setFrequency(double newFrequency);
    /**< Sets the frequency for the oscillator. */

    void setLevel(double newLevel);
    /**< Sets the amplitude level (in dB) for the oscillator. */

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double* inL, double* inR, double* outL, double* outR);
    /**< Calculates a stereo-ouput frame. */

    //=============================================================================================

  protected:

    //double sampleRate;
    //double frequency;
    //double level;

    SineOscillatorStereo sineOscillator;

  };

  //---------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE void TestGenerator::getSampleFrameStereo(double* inL,  double* inR, 
                                                  double* outL, double* outR)
  {
    sineOscillator.getSampleFrameStereo(outL, outR);
  }

} // end namespace rosic

#endif // rosic_TestGenerator_h
