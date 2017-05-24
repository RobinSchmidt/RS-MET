#ifndef rosic_FractionalDelayLineStereo_h
#define rosic_FractionalDelayLineStereo_h

//// rosic-indcludes:
//#include "../basics/rosic_Interpolator.h"

namespace rosic
{

  /**

  This class implements a basic delay-line with various interpolation methods. 

  \todo: make a mono version

  */

  class FractionalDelayLineStereo  
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a delay-line with a given maximum number of samples delay. */
    FractionalDelayLineStereo(int maximumDelayInSamples = 65536);

    /** Destructor */
    ~FractionalDelayLineStereo();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the interpolation method. */
    void setInterpolationMethod(int newMethod);

    /** Sets the delay-time in seconds. */
    void setDelayTime(double newDelayTime); 

    /** Sets the delay-time in seconds. */
    void setDelayTimeInMilliseconds(double newDelayTimeInMilliseconds) { setDelayTime(0.001*newDelayTimeInMilliseconds); }

    //---------------------------------------------------------------------------------------------
    // inquiry (get-, is-, etc. functions):

    /** Returns the delay-time in seconds. */
    double getDelayTime() const { return delayInSeconds; }

    /** Returns the interpolation method. */
    int getInterpolationMethod() const { return interpolatorL.getInterpolationMethod(); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output-sample at a time. */
    INLINE double getSample(double in);

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the content of the delayline to all zeros. */
    void clearDelayBuffers();

    //---------------------------------------------------------------------------------------------
    // embedded public audio-modules:

    Interpolator interpolatorL, interpolatorR;

  //===============================================================================================

  protected:

    static const int interpolatorMargin = 8;
      // The allocated memory will be a bit larger than the required delayline-length in order to 
      // make life easier for the interpolator (such that the interpolator is not concerned with 
      // buffer-wraparounds). This is the number of samples which the buffer is longer. */

    int    tapIn, tapOut;
    double *delayBufferL, *delayBufferR;

    double frac; 
      // The actual readout-position is this (fractional) number of samples ahead the 
      // tapOut-pointer position. It is given by  1.0 - delayInSampleFractionalPart. */

    double delayInSamplesFractionalPart;
    double delayInSeconds;
    double delayInSamples;
    int    delayInSamplesIntegerPart;
    int    maxDelayInSamples;
    double sampleRate;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  double FractionalDelayLineStereo::getSample(double in)
  {
    double out;

    // write the incoming sample into the delay-line:
    delayBufferL[tapIn] = in;

    // if the tap-pointer is smaller than the interpolator-margin, we have to write the sample 
    // behind the end of the used length of the delay-line as well:
    if( tapIn < interpolatorMargin )
      delayBufferL[maxDelayInSamples+tapIn] = in;

    // calculate the output-sample by invoking the embedded Interpolator-object:
    out = interpolatorL.getSample(frac, &(delayBufferL[tapOut]));

    // increment tap-pointers:
    tapIn  = wrapAround(tapIn+1,  maxDelayInSamples);
    tapOut = wrapAround(tapOut+1, maxDelayInSamples);

    return out;
  }

  void FractionalDelayLineStereo::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    // write the incoming samples into the delay-lines:
    delayBufferL[tapIn] = *inOutL;
    delayBufferR[tapIn] = *inOutR;

    // if the tap-pointer is smaller than the interpolator-margin, we have to write the sample 
    // behind the end of the used length of the delay-line as well:
    if( tapIn < interpolatorMargin )
    {
      delayBufferL[maxDelayInSamples+tapIn] = *inOutL;
      delayBufferR[maxDelayInSamples+tapIn] = *inOutR;
    }

    // calculate the output-samples by invoking the embedded Interpolator-object:
    *inOutL = interpolatorL.getSample(frac, &(delayBufferL[tapOut]));
    *inOutR = interpolatorR.getSample(frac, &(delayBufferR[tapOut]));

    // increment tap-pointers:
    tapIn  = wrapAround(tapIn+1,  maxDelayInSamples);
    tapOut = wrapAround(tapOut+1, maxDelayInSamples);
  }

} // end namespace rosic

#endif // #ifndef rosic_FractionalDelayLineStereo_h
