#ifndef rosic_FractionalDelayLine_h
#define rosic_FractionalDelayLine_h

//// rosic-indcludes:
//#include "../basics/rosic_Interpolator.h"

namespace rosic
{

  /**

  This class implements a basic delay-line with various interpolation methods.

  \todo: define copy constructor to create a deep copy of the delayBuffer

  */

  class FractionalDelayLine
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a delay-line with a given maximum number of samples delay. */
    FractionalDelayLine(int maximumDelayInSamples = 65536);

    /** Destructor */
    ~FractionalDelayLine();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the delay-time in seconds or beats (depending on whether sync is active). */
    void setDelayTime(double newDelayTime);

    /** Switches the tempo-sync on or off. */
    void setSyncMode(bool shouldTempoSync);

    /** Sets up the tempo in  beats per minute. */
    void setTempoInBPM(double newTempoInBPM);

    /** Sets the interpolation method. */
    void setInterpolationMethod(int newMethod);

    //---------------------------------------------------------------------------------------------
    // inquiry (get-, is-, etc. functions):

    /** Returns the delay-time in seconds or beats (depending on whether sync is active). */
    double getDelayTime() const { return delayTime; }

    /** Returns true when tempo-sync is active, false otherwise. */
    int isInSyncMode() const { return tempoSync; }

    /** Returns the interpolation method. */
    int getInterpolationMethod() const { return interpolator.getInterpolationMethod(); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output-sample at a time. */
    INLINE double getSample(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the content of the delayline to all zeros. */
    void clearDelayBuffer();

    //=============================================================================================

  protected:

    /** Wraps an integer (read/write) position into the permitted range (0...length-1). */
    INLINE int wrapAround(int position);

    /** Sets up the delay-time in samples according to the chosen delayTime, sync-mode and
    sample-rate user parameters. */
    void setupDelayInSamples();

    static const int interpolatorMargin = 1;
    //static const int interpolatorMargin = 8;
      // The allocated memory will be a bit larger than the required delayline-length in order to
      // make life easier for the interpolator (such that the interpolator is not concerned with
      // buffer-wraparounds). This is the number of samples which the buffer is longer.
      // a large margin imposes long minimum delay time (minimum = margin-1), but allows for higher
      // oder interpolation

    int    tapIn, tapOut;

    double *delayBuffer;

    int    length;
      // nominal length (excluding the interpolator margin, maximum delay will be length-1

    double frac;
      // The actual readout-position is this (fractional) number of samples ahead the
      // tapOut-pointer position. It is given by  1.0 - delayInSampleFractionalPart. */


    double delayInSamples;
    double delayTime;      // in seconds or beats
    double sampleRate;
    double bpm;
    bool   tempoSync;

    Interpolator interpolator;

  private:

    // make assignment operator and copy constructor unavailable because this class contains
    // pointer members:
    FractionalDelayLine& operator=(const FractionalDelayLine& /*other*/) { return *this; }
    FractionalDelayLine(const FractionalDelayLine& /*other*/) { }

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE int FractionalDelayLine::wrapAround(int position)
  {
    while( position >= length )
      position =- length;
    while( position < 0 )
      position += length;
    return position;
  }

  double FractionalDelayLine::getSample(double in)
  {
    doubleA out;

    // write the incoming sample into the delay-line:
    delayBuffer[tapIn] = in;

    // if the tap-pointer is smaller than the interpolator-margin, we have to write the sample
    // behind the end of the used length of the delay-line as well:
    if( tapIn < interpolatorMargin )
      delayBuffer[length+tapIn] = in;

    // calculate the output-sample by invoking the embedded Interpolator-object:
    out = interpolator.getSample(frac, &(delayBuffer[tapOut]));

    // increment tap-pointers:
    tapIn  = wrapAround(tapIn+1);
    tapOut = wrapAround(tapOut+1);

    return out;
  }

} // end namespace rosic

#endif // #ifndef rosic_FractionalDelayLine_h
