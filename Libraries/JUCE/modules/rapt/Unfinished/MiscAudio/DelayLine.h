#ifndef RAPT_DELAYLINE_H
#define RAPT_DELAYLINE_H

  // this file contains a couple of delayline classes with increasingly complex functionality

/** This class implements a basic delay-line which allows only for integer delays. No interpolation
algorithm is involved - that makes it especially efficient.

\todo: use only one pointer for tapIn and tapOut (see Julius Smith's pasp-book) */

class rsBasicDelayLine
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsBasicDelayLine();

  /** Destructor */
  ~rsBasicDelayLine();


  /** \name Setup */

  /** Sets the maximum delay in samples that can be uesed - i.e. the length of the internal
  delayline. For effciiency reasion, this should be a power-of-two-minus-one (such that the
  length of the delayline can be the respective power-of-two itself) - if it isn't, the next
  power-of-two-minus-one will be used. */
  void setMaximumDelayInSamples(int newMaxDelay);

  /** Sets the delay-time in samples. If the passed value exceeds the length of the delayline,
  new memory will be allocated which is large enough to support the desired delay. You probably
  want to avoid this (this could introduce artifacts) by calling setMaximumDelayInSamples in some
  safe place in your client code with a value that is larger than the largest delay, you
  expect. */
  void setDelayInSamples(int newDelay);


  /** \name Inquiry */

  int getDelayInSamples() const
  {
    int delay = tapIn - tapOut;
    if(delay < 0)
      delay += maxDelay+1;
    return delay;
  };


  /** \name Audio Processing */

  /** Calculates one output-sample at a time and handles all the tap-pointer increments. */
  RS_INLINE double getSample(double in);

  /** Calculates one output-sample at a time but suppresses the incrementation of the tapIn
  tapOut pointers and should be used in conjunction with incrementTapPointers(). This is useful
  when something must be added to the input which is not yet available at the time of/before
  calling getSample() - for example, when the delayline is part of a feedback loop. The adding
  itself can the be done via addToInput(). */
  RS_INLINE double getSampleSuppressTapIncrements(double in);

  /** Adds some signal value to the current tapIn-position in the delayLine - useful for
  feedback and crossfeedback stuff. */
  RS_INLINE void addToInput(double signalToAdd);

  /** Does the increment for the tap pointers and wraps them around if necesarray - should be
  used in conjunction with getSampleSuppressTapIncrements(). */
  RS_INLINE void incrementTapPointers();


  /** \name Misc */

  /** Resets the content of the delayline contents to all zeros. */
  void reset();

protected:

  /** \name Data */

  double *delayLine;
  int    tapIn, tapOut, maxDelay;

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

RS_INLINE double rsBasicDelayLine::getSample(double in)
{
  incrementTapPointers();
  return getSampleSuppressTapIncrements(in);
}

RS_INLINE double rsBasicDelayLine::getSampleSuppressTapIncrements(double in)
{
  delayLine[tapIn] = in;
  return delayLine[tapOut];
}

RS_INLINE void rsBasicDelayLine::addToInput(double signalToAdd)
{
  delayLine[tapIn] += signalToAdd;
}

RS_INLINE void rsBasicDelayLine::incrementTapPointers()
{
  tapIn  = (tapIn+1)  & maxDelay;
  tapOut = (tapOut+1) & maxDelay;

  //tapIn  = (++tapIn)  & maxDelay;
  //tapOut = (++tapOut) & maxDelay;
    // it's crucial to use pre-increment "++tapIn" rather than post-increment "tapIn++" because
    // with post-increment, the bitmask will be applied before incrementing which results in
    // reading/writing one sample behind the allocated memory
}

//===============================================================================================

/**

Extends BasicDelayLine by keeping information about the sample rate and the delay in
seconds.

\todo: facilitate temp-sync by maintaining a bpm-value and a sync-flag
...maybe do this in a subclass...

*/

class rsDelayLine : public rsBasicDelayLine
{

public:

  /** \name Construction/Destruction */

  /** Constructor - constructs a delay-line with a given maximum number of samples delay. This
  has to be a power of two minus 1 - otherwise the next power of two minus 1 will be used. */
  rsDelayLine();

  /** Destructor */
  ~rsDelayLine();


  /** \name Setup */

  /** Sets the sample-rate. */
  void setSampleRate(double newSampleRate);

  /** Sets the delay-time in samples. */
  void setDelayInSamples(int newDelayInSamples);

  /** Sets the delay-time in seconds. */
  void setDelayInSeconds(double newDelayInSeconds);

  /** Sets the delay-time in milliseconds. */
  void setDelayInMilliseconds(double newDelayInMilliseconds);


  /** \name Inquiry */

  /** Returns the delay-time in seconds. */
  RS_INLINE double getDelayInSeconds() const { return delayInSeconds; }

  /** Returns the delay-time in milliseconds. */
  RS_INLINE double getDelayInMilliseconds() const { return 1000.0 * delayInSeconds; }

protected:

  /** \name Data */

  double delayInSeconds;
  double sampleRate;

};

//===============================================================================================

/**

This class implements a basic delay-line with various interpolation methods.

\todo: define copy constructor to create a deep copy of the delayBuffer
\todo: maybe get rid of the tempo-sync stuff - this is actually something for a higher level
       it would probably be best, if the delay time is set up in samples such that the class can
       be agnostic of the samplerate

*/

class rsFractionalDelayLine
{

public:

  /** \name Construction/Destruction */

  /** Constructor - constructs a delay-line with a given maximum number of samples delay. */
  rsFractionalDelayLine(int maximumDelayInSamples = 65536);

  /** Destructor */
  ~rsFractionalDelayLine();


  /** \name Setup */

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


  /** \name Inquiry */

  /** Returns the delay-time in seconds or beats (depending on whether sync is active). */
  double getDelayTime() const { return delayTime; }

  /** Returns true when tempo-sync is active, false otherwise. */
  int isInSyncMode() const { return tempoSync; }

  /** Returns the interpolation method. */
  int getInterpolationMethod() const { return interpolator.getInterpolationMethod(); }


  /** \name Audio Processing */

  /** Calculates one output-sample at a time. */
  RS_INLINE double getSample(double in);


  /** \name Misc */

  /** Resets the content of the delayline to all zeros. */
  void clearDelayBuffer();


protected:

  /** Wraps an integer (read/write) position into the permitted range (0...length-1). */
  RS_INLINE int wrapAround(int position);

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

  rsInterpolator interpolator;

private:

  // make assignment operator and copy constructor unavailable because this class contains 
  // pointer members:
  rsFractionalDelayLine& operator=(const rsFractionalDelayLine &other) { return *this; }
  rsFractionalDelayLine(const rsFractionalDelayLine& other) { }

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

RS_INLINE int rsFractionalDelayLine::wrapAround(int position)
{
  while(position >= length)
    position =-length;
  while(position < 0)
    position += length;
  return position;
  // optimize using a bitmask
}

double rsFractionalDelayLine::getSample(double in)
{
  doubleA out;

  // write the incoming sample into the delay-line:
  delayBuffer[tapIn] = in;

  // if the tap-pointer is smaller than the interpolator-margin, we have to write the sample
  // behind the end of the used length of the delay-line as well:
  if(tapIn < interpolatorMargin)
    delayBuffer[length+tapIn] = in;

  // calculate the output-sample by invoking the embedded Interpolator-object:
  out = interpolator.getSample(frac, &(delayBuffer[tapOut]));

  // increment tap-pointers:
  tapIn  = wrapAround(tapIn+1);
  tapOut = wrapAround(tapOut+1);

  return out;
}

#endif
