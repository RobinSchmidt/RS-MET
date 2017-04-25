#ifndef rosic_DelayLineStereo_h
#define rosic_DelayLineStereo_h

// rosic-indcludes:
#include "../math/rosic_ElementaryFunctionsReal.h"
#include "../basics/rosic_Interpolator.h"
#include <new>

namespace rosic
{

  /**

  This class implements serves as baseclass for effects that need two delaylines (for left and 
  right channel) such as true-stereo versions of vibrato, chorus, etc.

  */

  class DelayLineStereo
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - will allocate a delay-buffer with a given maximum number of samples 
    delay. */
    DelayLineStereo(int bufferLengthToAllocate);

    /** Destructor */
    ~DelayLineStereo();

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Feeds a stereo input sample-frame into the delayline (taking duplication of samples for the 
    interpolator into account) - should be wrapped into acquireLock()/releaseLock(). */
    INLINE void feedIn(double inL, double inR);

    /** Increments the write pointer by one sample (taking wraparound into account) - should be 
    wrapped into acquireLock()/releaseLock(). */
    INLINE void incrementWritePointer();

    /** Returns the output of the left delayline at some given amount of delay, using cubic hermite
    interpolation - should be wrapped into acquireLock()/releaseLock(). */
    INLINE double getLeftOutputHermiteAt(double delayInSamples);

    /** Returns the output of the right delayline at some given amount of delay, using cubic hermite
    interpolation - should be wrapped into acquireLock()/releaseLock(). */
    INLINE double getRightOutputHermiteAt(double delayInSamples);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Sets the content of the delayline to zeros. */
    void clearBuffers();

    //=============================================================================================

  protected:

    /** Converts a desired fractional delay (in samples) into integer and fractional part of the 
    delayline readout position. */
    INLINE void delayToIntAndFrac(double delayInSamples, int &iPart, double &fPart);

    static const int interpolatorMargin = 3;

    int     tapIn;               // write pointer
    int     length;              // delayline length
    double  *bufferL, *bufferR;  // buffers for left and right channel 

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void DelayLineStereo::feedIn(double inL, double inR)
  {
    // feed the signals into the delay-buffers:
    bufferL[tapIn] = inL;
    bufferR[tapIn] = inR;

    // if the tap-pointer is smaller than the interpolator-margin, we have to write the sample
    // behind the end of the used length of the delay-line as well:
    if( tapIn < interpolatorMargin )
    {
      bufferL[length+tapIn] = inL;
      bufferR[length+tapIn] = inR;
    }
  }

  INLINE void DelayLineStereo::incrementWritePointer()
  {
    tapIn = wrapAround(tapIn+1, length);
  }
 
  INLINE double DelayLineStereo::getLeftOutputHermiteAt(double delayInSamples)
  {
    int    iPart;
    double fPart;
    delayToIntAndFrac(delayInSamples, iPart, fPart);
    return Interpolator::getSampleHermite4p3o(fPart, &(bufferL[iPart]));
  }

  INLINE double DelayLineStereo::getRightOutputHermiteAt(double delayInSamples)
  {
    int    iPart;
    double fPart;
    delayToIntAndFrac(delayInSamples, iPart, fPart);
    return Interpolator::getSampleHermite4p3o(fPart, &(bufferR[iPart]));
  }

  INLINE void DelayLineStereo::delayToIntAndFrac(double delayInSamples, int &iPart, double &fPart)
  {
    iPart = floorInt(delayInSamples);  // integer part of the delay
    fPart = delayInSamples - iPart;    // fractional part of the delay
    fPart = 1.0 - fPart;               // fractional part of the read-position
    if( fPart >= 1.0 )                 // wrap to open interval [0.0...1.0[
    {
      fPart  = 0.0;
      iPart -= 1;
    }
    iPart = wrapAround(tapIn-iPart-1, length);   // integer part of the read-position
  }

} // end namespace rosic

#endif // #ifndef rosic_DelayLineStereo_h
