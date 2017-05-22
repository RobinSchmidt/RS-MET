#ifndef rosic_CombResonator_h
#define rosic_CombResonator_h

//// rosic-indcludes:
//#include "rosic_CombFilter.h"
//#include "rosic_DampingFilter.h"

namespace rosic
{

  /**

  This class extends the CombFilter class by incorporating a biquadratic low-/high- shelving filter
  into the feedback loop which allows the user to set up decay times for mid, low and high 
  fequencies seperately. These decay-times provide a more musical parametrization of the feedback
  gain. Moreover, it will be possible to set the comb into a mode in which it only generates odd
  harmonics - this supersedes the negative feedback option of the baseclass (which also went an 
  octave lower on negative feedback, which we avoid here).

  */

  class CombResonator : public CombFilter
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - will allocate a delay-buffer with a given maximum number of samples delay. */
    CombResonator(int bufferLengthToAllocate = 65536);

    /** Destructor */
    ~CombResonator();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the fundamental frequency of the comb filter (in case of positive feedback). */
    void setFrequency(double newFrequency);

    /** Sets the time for the comb to decay to -60 dB (in seconds). */
    void setDecayTime(double newDecayTime);

    /** Scales the decay-time for the high-frequency band. */
    void setHighDecayScale(double newScale);

    /** Scales the decay-time for the low-frequency band. */
    void setLowDecayScale(double newScale);

    /** Sets the crossover-frequency between mid and high frequencies. */
    void setHighCrossoverFreq(double newFreq);

    /** Sets the crossover-frequency between low and mid frequencies. */
    void setLowCrossoverFreq(double newFreq);

    /** Switches the resonator into a mode where it produces only odd harmonics. */
    void setOddOnlyMode(bool shouldCreateOnlyOddHarmonics);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE double getSample(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the content of the delayline to all zeros and resets the feedback filter. */
    void clearBuffer();

    //=============================================================================================

  protected:

    /** Sets up the delay-time in samples according to the chosen frequency and sample-rate user 
    parameters. */
    void setupDelayInSamples();

    /** Sets up the feedback- and correction filter such as to obtain the desired frequency 
    dependent decay characteristic. */
    void setupFilters();

    double decayTime, highDecayScale, lowDecayScale, highCrossOver, lowCrossOver;
    bool   oddOnlyMode;

    DampingFilter feedbackFilter, correctionFilter;

  private:

    // make assignment operator and copy constructor unavailable because this class contains 
    // pointer members:
    CombResonator& operator=(const CombResonator &other) { return *this = other; }
    CombResonator(const CombResonator& other)            {        *this = other; }

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  double CombResonator::getSample(double in)
  {
    // write the incoming sample into the delay-line:
    buffer[tapIn] = in + feedbackFilter.getSample(yOld);

    // if the tap-pointer is smaller than the interpolator-margin, we have to write the sample
    // behind the end of the used length of the delay-line as well:
    if( tapIn < interpolatorMargin )
      buffer[length+tapIn] = buffer[tapIn];

    // calculate the output-sample by invoking the embedded Interpolator-object:
    yOld = interpolator.getSample(frac, &(buffer[tapOut]));

    // increment tap-pointers:
    tapIn  = wrapAround(tapIn+1,  length);
    tapOut = wrapAround(tapOut+1, length);

    return (1.0-dryWetRatio)*in + dryWetRatio*wetPolarity*correctionFilter.getSample(yOld);
  }

} // end namespace rosic

#endif // #ifndef rosic_CombResonator_h
