#ifndef rosic_OscillatorBank_h
#define rosic_OscillatorBank_h

// rosic-indcludes:
#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This class implements bank of oscillators of the sort as defined in SineOscillatorStereo.
  
  \todo implement linear frequency ramping for each voice - this
  can use the same recursion for the evaluation of cosines of successive radian frequencies (which
  are needed for the filter coefficient).

  */

  class OscillatorBank  
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    OscillatorBank(int newMaxNumSines = 512);   
    /**< Constructor. The input parameter chooses, how many voices this OscillatorBank will be able
    to generate at maximum. */

    ~OscillatorBank();  
    /**< Destructor */

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    void setSampleRate(double newSampleRate);
    /**< Sets the sample-rate. */

    void setNumSines(int newNumSines);
    /**< Sets the number of voices (up to maxNumSines which can be set only from the 
    constructor). */

    void setFundamentalFrequency(double newFundamentalFrequency);
    /**< Sets the fundamental frequency. */

    void setMasterAmplitude(double newMasterAmplitude);
    /**< Sets the master amplitude. */

    void setRelativeFrequencies(double* newRelativeFrequencies);
    /**< Sets the relative frequencies of the sinusoids with respect to the fundamental 
    frequency as factors. The array passed must have (at least) maxNumSines entries. */

    void setRelativeAmplitudes(double* newRelativeAmplitudes);
    /**< Sets the relative amplitudes (as raw multiplier) of the sinusoids with respect to the 
    master amplitude. The array passed must have (at least) maxNumSines entries.  */

    void setPhases(double* newPhases);
    /**< Sets the phases (in radians) of the sinusoids to be generated.  The array passed must 
    have (at least) maxNumSines entries.  */

    //---------------------------------------------------------------------------------------------
    // event handling:

    void triggerAll();
    /**< Resets all the oscillator-voices up to maxnumVoices to their start phase. */

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double *outL, double *outR);
    /**< Calculates one output-sample frame at a time. */

    //---------------------------------------------------------------------------------------------
    // others:

    void initToSawWave();
    /**< Initializes the relative frequencies and amplitudes and the start phases such that a 
    sawtooth wave will be produced (on the left channel and its hilbert transfom on the right 
    channel). */

    //=============================================================================================

  protected:

    void calculateAllCoefficients();
    /**< Triggers a calculation of all recursion coefficients up to maxNumSines. */

    double *a1;           /**< The filter coefficients for the recursions. */
    double *s1, *s2;      /**< Past sine outputs (delayed by one or two samples respectively). */
    double *c1, *c2;      /**< Past cosine outputs (delayed by one or two samples respectively). */


    double fundamentalFrequency; /**< The fundamental frequency in Hz. */
    double masterAmplitude;      /**< The  master amplitude as raw multiplier. */
    double sampleRate;           /**< The sample-rate in Hz. */
    int    numSines;             /**< Number of sines to generate. */
    int    maxNumSines;          /**< Maximum number of sines to generate. */

    double *relativeFrequencies; /**< The relative frequencies of the partials (as multiplier). */
    double *relativeAmplitudes;  /**< The relative amplitudes of the partials (as multiplier).   */
    double *phases;              /**< The start-phases of the oscillators in radians. */
  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void OscillatorBank::getSampleFrameStereo(double *outL, double *outR)
  {
    *outL = 0.0;
    *outR = 0.0;
    double tmp1, tmp2;
    for(int i=0; i<numSines; i++)
    {
      tmp1   = a1[i]*s1[i] - s2[i];
      s2[i]  = s1[i];
      s1[i]  = tmp1;
      *outL += relativeAmplitudes[i] * tmp1;

      tmp2   = a1[i]*c1[i] - c2[i];
      c2[i]  = c1[i];
      c1[i]  = tmp2;
      *outR += relativeAmplitudes[i] * tmp2;
    }
    *outL *= masterAmplitude;
    *outR *= masterAmplitude;
  }


} // end namespace rosic

#endif // #ifndef rosic_OscillatorBank_h
