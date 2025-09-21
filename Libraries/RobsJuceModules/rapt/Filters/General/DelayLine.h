#ifndef RAPT_DELAYLINE_H
#define RAPT_DELAYLINE_H

// This file contains a couple of delayline classes with increasingly complex functionality

/** This class implements a basic delay-line which allows only for integer delays. ...TBC...  */

template<class T>
class rsDelay   // Maybe rename to rsDelayLine, rsDelayLineInteger ..or just to rsDelay, or rsDelayInteger
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Constructor. */
  rsDelay();

  /** Destructor */
  ~rsDelay();


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the maximum delay in samples that can be uesed - i.e. the length of the internal
  delayline. For effciiency reasion, this should be a power-of-two-minus-one (such that the
  length of the delayline can be the respective power-of-two itself) - if it isn't, the next
  power-of-two-minus-one will be used. */
  void setMaxDelayInSamples(int newMaxDelay);

  /** Sets the delay-time in samples. If the passed value exceeds the length of the delayline,
  new memory will be allocated which is large enough to support the desired delay. You probably
  want to avoid this (this could introduce artifacts) by calling setMaxDelayInSamples in some
  safe place in your client code with a value that is larger than the largest delay, you
  expect. */
  void setDelayInSamples(int newDelay);


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the delay in samples that this delayline produces. */
  int getDelayInSamples() const
  {
    int delay = tapIn - tapOut;
    if(delay < 0)
      delay += maxDelay+1;
    return delay;

    // ToDo: Document why we need to add maxDelay+1. One might expect that we should add maxDelay 
    // and suspect a bug here. But apparently, it's actually correct that way. There's a unit test
    // for this and it passes just fine. I think it's related to the maxDelay always being a power 
    // of two minus one.
  };

  /** Returns the maximum delay that this delayline can produce */
  int getMaxDelayInSamples() const { return maxDelay; }
  // Maybe rename to getMaxDelay. Same for the other delay setters/getters. Maybe adopt the 
  // convention that in higher level classes that accept a delay time in a physical unit (i.e. 
  // seconds), the functions should be named setDelayTime() etc. But then it seems a bit tedious to
  // apply that convention consistently. What about envelope followers and compressors, etc.? Should
  // we then also call their setters setAttackTime() rather than setAttack() etc.? For consistency,
  // we really should. We'll see....

  /** Returns the value of the transfer function H(z) at the given value of z. If M is the delay in
  samples, then H(z) = z^-M. The function has its own template parameter TArg because the argument 
  type will typically be different from the type T with which the class is instantiated. A higher 
  level class that uses template parameters TSig, TPar for parameters and signals may instantiate 
  the rsDelay class with T = TSig and call getTransferFunctionAt() with an argument of type 
  rsComplex<TPar>, for example. */
  template<class TArg>
  TArg getTransferFunctionAt(const TArg& z) const
  {
    int M = getDelayInSamples();     // M is our delay
    return rsPow(z, TArg(-M));       // H(z) = z^-M
  }
  // Needs unit tests!
  

  template<class TTol>
  void getTransferFunction(rsSparseTransferFunction<T, TTol>* tf) const
  {
    int M = getDelayInSamples();

    //// Old:
    //tf->num._setNumTerms(1); tf->num._setTerm(0, T(1), M);
    //tf->den._setNumTerms(1); tf->den._setTerm(0, T(1), 0);

    // New:
    tf->_setNumTerms(1, 1);
    tf->_setNumeratorTerm(  0, T(1), M);
    tf->_setDenominatorTerm(0, T(1), 0);
  }
  // Needs unit tests! 
  // 
  // Maybe drag this function out of the class as a free function like:
  // 
  //   rsGetTransferFunction(const rsDelay<T>& delay, rsSparseTransferFunction<T, TTol>* tf)
  // 
  // I don't really like to introduce a coupling to class rsSparseTransferFunction here. I think,
  // that generally thorughout the library, when we need functions that glue together two classes
  // (like here rsDelay and rsSparseTransferFunction), they should preferably implemented as
  // free functions (maybe) unless one of the two classes *very* basic - so basic that the 
  // introduction of the coupling doesn't matter because we can safely assume the basic class to be
  // avaibable (like #included or #imported or something) anyway. We actually can here because
  // rsSparseTransferFunction is part of the Math submodule which is more basic than the 
  // AudioBasics submodule (to which rsDelayLine belongs) - but still... Maybe someday I want to
  // reorganize things. But: As a free function, it will be harder to discover. But maybe its
  // ok when the free function resides in this file. It's less discovereable than when it's a 
  // member function but perhaps still divcoverable enough.




  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Calculates one output-sample at a time and handles all the tap-pointer increments. */
  RS_INLINE T getSample(T in);

  /** Calculates one output-sample at a time but suppresses the incrementation of the tapIn
  tapOut pointers and should be used in conjunction with incrementTapPointers(). This is useful
  when something must be added to the input which is not yet available at the time of/before
  calling getSample() - for example, when the delayline is part of a feedback loop. The adding
  itself can the be done via addToInput(). */
  RS_INLINE T getSampleSuppressTapIncrements(T in);
  // rename to getSampleNoTapIncrement or getSampleNoUpdate


  RS_INLINE void writeInput(T in);

  /** Adds some signal value to the current tapIn-position in the delayLine - useful for
  feedback and crossfeedback stuff. */
  RS_INLINE void addToInput(T signalToAdd);

  RS_INLINE void addToInputAt(T signalToAdd, int delay);




  /** Does the increment for the tap pointers and wraps them around if necesarray - should be
  used in conjunction with getSampleSuppressTapIncrements(). */
  RS_INLINE void incrementTapPointers();
  // rename to updateTaps




  // New functions - more convenient in certain situations:

  /** Reads the content of the delayline at the current tapOut pointer which is always M samples
  behind the tapIn pointer where M is the delay in samples. This readout triggers no update action.
  It's sometimes convenient to control the read/write/update steps from outside. That's why this
  function exists. */
  inline T readOutput() const { return delayLine[tapOut]; }

  /** Reads the content of the delayline at an arbitrary delay, i.e. at a position that is 
  independent from our current tapOut pointer. This can be used to implement multitap delaylines.
  We can just read the delayline whereever we want. */
  inline T readOutputAt(int delay) const 
  { 
    int readPos = tapIn - delay;   // Compute nominal read index
    readPos = readPos & maxDelay;  // Apply bit masking for wrap-around behavior
    return delayLine[readPos];     // Read out the delayline
  }

  /** Reads the content of the delayline at a position of the current tapOut pointer minus some
  additional delay. That means, when the delayline is set up to give a delay of M samples and the
  normal readOutput function would return x[n-M], this function here returns x[n-(M+k)] = x[n-M-k]
  where k is the desired additional delay. This function facilitates the implementation of 
  fractional delaylines on top of the integer one. For example, for a delay of 10.3 samples, one
  could set up the delayline to 10 samples, retrieve x[n-10] via readOutput and x[n-11] via
  readOutputWithAdditionalDelay(1) and then compute the output va linear interpolation: 
  y = 0.7*x[n-10] + 0.3*x[n-11].  */
  inline T readOutputWithAdditionalDelay(int additionalDelay) const
  {
    int readPos = tapOut - additionalDelay;

    readPos = readPos & maxDelay;
    // Is this correct? If so, document why. I think this function is already successfully in use 
    // in the CombAllpass delays (Verify!), so it is supposedly indeed correct. Conceptually, I 
    // think, it should behave like  readPos += maxDelay  whenever  readPos < 0. Does the 2s 
    // complement representation take care of this behavior?

    return delayLine[readPos];
  }
  // Needs tests


  inline void writeInputNoUpdate(T in) { delayLine[tapIn] = in; }

  inline void writeInputAndUpdate(T in)
  {
    writeInputNoUpdate(in);
    incrementTapPointers();
  }
  // rename to writeInputAndUpdate

  // Verify this and add it to the documentation:
  //
  // The caller should either use:
  //   y = dl.getSample(x);
  // or
  //   y = dl.readOutput();
  //   dl.writeInput(x);
  // or 
  //   y = dl.getSampleSuppressTapIncrements(0); 
  //   dl.addToInput(x); 
  //   dl.incrementTapPointers();
  // or
  //   y = dl.readOutput();
  //   dl.writeInputNoIncrement(x);
  //   dl.incrementTapPointers();



  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Resets the content of the delayline contents to all zeros. */
  void reset();

protected:


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  T* delayLine = nullptr;
  int tapIn = 0, tapOut = 0, maxDelay = 0;
  // ToDo: use std::vector for the delayLine. We may then get rid of maxDelay because it's stored
  // in the vector's size. ...or maybe capacity - depends on how we implement it.

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

template<class T>
RS_INLINE T rsDelay<T>::getSample(T in)
{
  incrementTapPointers();
  return getSampleSuppressTapIncrements(in);
}

template<class T>
RS_INLINE T rsDelay<T>::getSampleSuppressTapIncrements(T in)
{
  delayLine[tapIn] = in;    // Maybe use writeInput(in) instead
  return delayLine[tapOut];
}

template<class T>
RS_INLINE void rsDelay<T>::writeInput(T in)
{
  delayLine[tapIn] = in;
}

template<class T>
RS_INLINE void rsDelay<T>::addToInput(T signalToAdd)
{
  delayLine[tapIn] += signalToAdd;
}

template<class T>
RS_INLINE void rsDelay<T>::addToInputAt(T signalToAdd, int delay)
{
  //int p = tapIn - delay;
  //if(p < 0)
  //  p += maxDelay;  
    // Is that correct or should we use p += getDelayInSamples? That would be more costly, though
    // so I actually hope, it is correct this way. Compare to implementation of 
    // readOutputWithAdditionalDelay() and readOutputAt(). There, we do 
    // readPos = readPos & maxDelay;

  int writePos = tapIn - delay;
  writePos = writePos & maxDelay;
  delayLine[writePos] += signalToAdd;

  // Maybe optimize this to just:
  // 
  //   delayLine[(tapIn-delay) & maxDelay] += signalToAdd;
  //
  // and benchmark if it makes any difference.
}
// Needs unit tests!


template<class T>
RS_INLINE void rsDelay<T>::incrementTapPointers()
{
  tapIn  = (tapIn+1)  & maxDelay;
  tapOut = (tapOut+1) & maxDelay;

  //tapIn  = (++tapIn)  & maxDelay;
  //tapOut = (++tapOut) & maxDelay;
    // it's crucial to use pre-increment "++tapIn" rather than post-increment "tapIn++" because
    // with post-increment, the bitmask will be applied before incrementing which results in
    // reading/writing one sample behind the allocated memory. OK - it seems the parentheses solve
    // it also. But maybe switch back to pre-increment for potential performance gains (unlikely
    // but still)
}
// Maybe rename to incrementTaps() or updateTaps(). It's shorter, more descriptive and also more
// abstract in the desireable sense of hiding more implementation details that are suppsoed to be 
// irrelevant to the user. The user isn't interested in the questions whether or not we use 
// pointers (and strictly speaking, we don't - we use integers as pointer-offsets) and also whether
// we increment or decrement something or doing something even more weird.






//=================================================================================================

/**

This class implements a basic delay-line with various interpolation methods.

\todo: define copy constructor to create a deep copy of the delayBuffer
\todo: maybe get rid of the tempo-sync stuff - this is actually something for a higher level
       it would probably be best, if the delay time is set up in samples such that the class can
       be agnostic of the samplerate


ToDo:

- Try to get rid! We want to reclaim the name for a different implementation. But we can't because
  it's actually used in rsFakeResonanceFilter. But maybe we can turn this implementation into the 
  one, we want

- Or rename this one into rsTimeBasedDelayLine and factor out a class rsDelayTempoSynced that
  only has the stuff that is needed to set it up in terms of a fractional delay in samples. None
  of that sampleRate, bpm, tempoSync stuff   ...or rsDelayTempoSynced

*/

template<class TSig, class TPar>
class rsDelayTempoSynced
{

public:

  /** \name Construction/Destruction */

  /** Constructor - constructs a delay-line with a given maximum number of samples delay. */
  rsDelayTempoSynced(int maximumDelayInSamples = 65536);

  /** Destructor */
  ~rsDelayTempoSynced();


  /** \name Setup */

  /** Sets the sample-rate. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the delay-time in seconds or beats (depending on whether sync is active). */
  void setDelayTime(TPar newDelayTime);

  /** Switches the tempo-sync on or off. */
  void setSyncMode(bool shouldTempoSync);

  /** Sets up the tempo in  beats per minute. */
  void setTempoInBPM(TPar newTempoInBPM);

  /** Sets the interpolation method. */
  void setInterpolationMethod(int newMethod);


  /** \name Inquiry */

  /** Returns the delay-time in seconds or beats (depending on whether sync is active). */
  TPar getDelayTime() const { return delayTime; }

  /** Returns true when tempo-sync is active, false otherwise. */
  int isInSyncMode() const { return tempoSync; }

  /** Returns the interpolation method. */
  int getInterpolationMethod() const { return interpolator.getInterpolationMethod(); }


  /** \name Audio Processing */

  /** Calculates one output-sample at a time. */
  RS_INLINE TSig getSample(TSig in);


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

  TSig *delayBuffer;

  int    length;
  // nominal length (excluding the interpolator margin, maximum delay will be length-1

  TPar frac;
  // The actual readout-position is this (fractional) number of samples ahead the
  // tapOut-pointer position. It is given by  1.0 - delayInSampleFractionalPart. */


  TPar delayInSamples;
  TPar delayTime;      // in seconds or beats
  TPar sampleRate;
  TPar bpm;
  bool tempoSync;

  rsInterpolator<TSig> interpolator;

private:

  // Make assignment operator and copy constructor unavailable because this class contains 
  // pointer members:
  rsDelayTempoSynced& operator=(const rsDelayTempoSynced& /*other*/) { return *this; }
  rsDelayTempoSynced(const rsDelayTempoSynced& /*other*/) { }
  // ToDo: use a macro for this

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
RS_INLINE int rsDelayTempoSynced<TSig, TPar>::wrapAround(int position)
{
  while(position >= length)
    position =-length;
  while(position < 0)
    position += length;
  return position;
  // optimize using a bitmask
}

template<class TSig, class TPar>
TSig rsDelayTempoSynced<TSig, TPar>::getSample(TSig in)
{
  TSig out;

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
