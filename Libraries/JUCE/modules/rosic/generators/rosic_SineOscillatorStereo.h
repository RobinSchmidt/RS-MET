#ifndef rosic_SineOscillatorStereo_h
#define rosic_SineOscillatorStereo_h

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This class implements a pair of a sine and a cosine oscillator based on a self oscillating (free 
  running) second order recursive filter. 

  */

  class SineOscillatorStereo  
  {

    friend class FunctionNodeSineOscillator;

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    SineOscillatorStereo();   

    /** Destructor */
    ~SineOscillatorStereo();  

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the frequency of the sine/cosine to be generated. */
    void setFrequency(double newFrequency);

    /** Sets the amplitude (as raw multiplier) of the sine/cosine to be generated. */
    void setAmplitude(double newAmplitude);

    //---------------------------------------------------------------------------------------------
    // event handling:

    /** Resets the oscillator to its start phase. */
    void trigger();

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output-sample frame at a time. */
    INLINE void getSampleFrameStereo(double *outL, double *outR);

    //=============================================================================================

  protected:

    double a1;         /**< The filter coefficient for the recursion. */
    double s1, s2;     /**< Past sine outputs (delayed by one or two samples respectively). */
    double c1, c2;     /**< Past cosine outputs (delayed by one or two samples respectively). */

    double frequency;  /**< The frequency of the oscillator in Hz. */
    double amplitude;  /**< The amplitude of the oscillator as raw multiplier. */
    double startPhase; /**< The phase to which the oscillator is reset on trigger(). */
    double sampleRate; /**< The sample-rate. */
    double omega;      /**< Normalized radian frequency of the oscillator. */

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void SineOscillatorStereo::getSampleFrameStereo(double *outL, double *outR)
  {
    *outL  = a1*s1 - s2;
    s2     = s1;
    s1     = *outL;

    *outR  = a1*c1 - c2;
    c2     = c1;
    c1     = *outR;

    *outL *= amplitude;
    *outR *= amplitude;
  }

} // end namespace rosic

#endif // #ifndef rosic_SineOscillatorStereo_h
