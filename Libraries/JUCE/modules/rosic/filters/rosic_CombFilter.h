#ifndef rosic_CombFilter_h
#define rosic_CombFilter_h

//// rosic-indcludes:
//#include "../basics/rosic_Interpolator.h"
//#include <new>

namespace rosic
{

  /**

  This class implements a comb filter which may be used in flanger effects and the likes.

  */

  class CombFilter
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - will allocate a delay-buffer with a given maximum number of samples delay. */
    CombFilter(int bufferLengthToAllocate = 65536);

    /** Destructor */
    ~CombFilter();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the fundamental frequency of the comb filter (in case of positive feedback). */
    void setFrequency(double newFrequency);

    /** Sets the feedback around the allpass chain as raw factor. */
    void setFeedbackFactor(double newFeedbackFactor) 
    { feedbackFactor = newFeedbackFactor; } 

    /** Sets the polarity of the wet signal negative (or positive, if false). */
    void setNegativePolarity(bool shouldBeNegative) 
    { wetPolarity = - (double) (2*(int)shouldBeNegative-1); }

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) { dryWetRatio = newDryWet; }

    /** Sets the interpolation method. */
    void setInterpolationMethod(int newMethod);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE double getSample(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the content of the delayline to all zeros. */
    void clearBuffer();

    //=============================================================================================

  protected:

    /** Sets up the delay-time in samples according to the chosen frequency and sample-rate user 
    parameters. */
    void setupDelayInSamples();

    static const int interpolatorMargin = 1;
    //static const int interpolatorMargin = 8;
      // The allocated memory will be a bit larger than the required delayline-length in order to
      // make life easier for the interpolator (such that the interpolator is not concerned with
      // buffer-wraparounds). This is the number of samples which the buffer is longer.
      // a large margin imposes long minimum delay time (minimum = margin-1), but allows for higher
      // oder interpolation

    int    tapIn, tapOut;
    double *buffer;
    int    length;           // nominal length (excluding the interpolator margin)
    double yOld;             // old output sample to be used for feedback
    double frac;             // actual readout-position is this (fractional) number of samples ahead 
                             // the tapOut-pointer position. 
    double delayInSamples;
    double frequency;
    double feedbackFactor;
    double dryWetRatio;      // dry/wet ratio 0....1
    double sampleRate;
    double wetPolarity;

    Interpolator interpolator;

  private:

    // make assignment operator and copy constructor unavailable because this class contains 
    // pointer members:
    CombFilter& operator=(const CombFilter &other) { return *this = other;  }
    CombFilter(const CombFilter& other)            { *this = other; }

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  double CombFilter::getSample(double in)
  {
    // write the incoming sample into the delay-line:
    buffer[tapIn] = in+feedbackFactor*yOld;

    // if the tap-pointer is smaller than the interpolator-margin, we have to write the sample
    // behind the end of the used length of the delay-line as well:
    if( tapIn < interpolatorMargin )
      buffer[length+tapIn] = buffer[tapIn];

    // calculate the output-sample by invoking the embedded Interpolator-object:
    yOld = interpolator.getSample(frac, &(buffer[tapOut]));

    // increment tap-pointers:
    tapIn  = RAPT::rsWrapAround(tapIn+1,  length);
    tapOut = RAPT::rsWrapAround(tapOut+1, length);

    return (1.0-dryWetRatio)*in + dryWetRatio*wetPolarity*yOld;
  }

} // end namespace rosic

#endif // #ifndef rosic_CombFilter_h
