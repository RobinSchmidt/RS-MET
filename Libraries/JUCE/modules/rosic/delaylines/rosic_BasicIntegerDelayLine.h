#ifndef rosic_BasicIntegerDelayLine_h
#define rosic_BasicIntegerDelayLine_h

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This class implements a basic delay-line which allows only for integer delays. No interpolation
  algorithm is involved - that makes it especially efficient.

  \todo: use only one pointer for tapIn and tapOut (see Julius Smith's pasp-book)

  */

  class BasicIntegerDelayLine
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a delay-line with a given maximum number of samples delay. This
    has to be a power of two otherwise the next power of two will be used. */
    BasicIntegerDelayLine(int maximumDelayInSamples = 65536);

    /** Destructor */
    ~BasicIntegerDelayLine();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the delay-time in samples and returns the value again. If the passed value exceeds the 
    maximum possible delay, the delayline will set itself to this maximum value and return this 
    maximum value. */
    int setDelayInSamples(int newDelayInSamples);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output-sample at a time and handles all the tap-pointer increments. */
    INLINE double getSample(double in);

    /** Calculates one output-sample at a time but suppresses the incrementation of the tapIn
    tapOut pointers and should be used in conjunction with incrementTapPointers(). This is useful
    when something must be added to the input which is not yet available at the time of/before
    calling getSample() -  for example, when the delayline is part of a feedback loop. The adding
    itself can the be done via addToInput(). */
    INLINE double getSampleSuppressTapIncrements(double in);

    /** Adds some signal value to the current tapIn-position in the delayLine - useful for
    feedback and crossfeedback stuff. */
    INLINE void addToInput(double signalToAdd);

    /** Does the increment for the tap pointers and wraps them around if necesarray - should be
    used in conjunction with getSampleSuppressTapIncrements(). */
    INLINE void incrementTapPointers();

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the content of the delayline to all zeros. */
    void clearDelayBuffer();

    //=============================================================================================

  protected:

    double *delayLine; 
    int    tapIn, tapOut, bitMask;  // the maximum possible delay will be bitMask+1
    
  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double BasicIntegerDelayLine::getSample(double in)
  {
    incrementTapPointers();
    return getSampleSuppressTapIncrements(in);
  }

  INLINE double BasicIntegerDelayLine::getSampleSuppressTapIncrements(double in)
  {
    delayLine[tapIn] = in;
    return delayLine[tapOut];
  }

  INLINE void BasicIntegerDelayLine::addToInput(double signalToAdd)
  {
    delayLine[tapIn] += signalToAdd;
  }

  INLINE void BasicIntegerDelayLine::incrementTapPointers()
  {
    //tapIn  = (tapIn++  & bitMask);
    //tapOut = (tapOut++ & bitMask); 

    // ...strange: this leads to an access violation when plugging out the StereoDelay plugin. 
    // the code below works...

    tapIn++;
    tapOut++;

    // wraparound:
    if( tapIn > bitMask )
      tapIn = 0;
    if( tapOut > bitMask )
      tapOut = 0;
  }

} // end namespace rosic

#endif // #ifndef rosic_BasicIntegerDelayLine_h
