#pragma once

//=================================================================================================

/** Implements a realtime moving maximum filter, that is: a filter that looks at audio samples
x[n], x[n-1],..., x[n-k] and extracts the maximum value in this range of samples. The algorithm is
based on a double ended queue and has an amortized complexity of O(1) per sample (a naive
implementation would have complexity O(k) per sample). For an explanation of the algorithm,
see: https://www.nayuki.io/page/sliding-window-minimum-maximum-algorithm  

see also here:
https://www.kvraudio.com/forum/viewtopic.php?f=33&t=465829 

todo: 
-maybe implement detection of intersample-peaks by fitting a parabola to triplets of samples
 and sloving for its peaks in cases, where a peak is detected - that could actually be a useful
 feature for any kind of envelope detector

*/

template<class T>
class rsMovingMaximumFilter  // rename to rsExtremumFilter, maybe move to Scientific
{

public:

  rsMovingMaximumFilter(size_t maxLength);
  // todo: give it a default value (maybe 0 or 1)

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the length of the filter, i.e. the number of samples within which a maximum is
  searched.
  Currently, this will also cause the delayline to be flushed - todo: allow for dynamic
  length changes. */
  void setLength(size_t newLength)
  {
    delayLine.setLength(newLength);
    reset();
    // preliminary - we need to put the filter into its initial state on length changes, so the
    // length cannot yet be modulated
    // todo: we may have to update the content of the maxDeque - when this is called during running
    // the filter and new length is shorter than the old, some values that are currently in the 
    // deque will never be removed - so we should remove them in the moment, when the length 
    // changes
  }

  // todo: add function setMaxLength, setLength with non-integer parameter (should crossfade between
  // floor(length) and ceil(length)) .oh - i think, it's not so easily possible to implement 
  // non-integer lengths



  /** Sets the "greater-than" comparison function. Note that you can actually also pass a function
  that implements a less-than comparison in which case the whole filter turns into a moving-minimum
  filter. */
  void setGreaterThanFunction(bool (*greaterThan)(const T&, const T&))
  { greater = greaterThan; }

  // maybe use std::function instead of plain function pointer:
  //void setComparisonFunction(const std::function<bool(const T&, const T&)>& greaterThan)
  //{ greater = greaterThan; }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */


  /** Returns the length of the filter, i.e. the number of samples within which a maximum is
  searched. */
  size_t getLength() const { return delayLine.getLength(); }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes and returns an output sample.  */
  inline T getSample(T in)
  {
    // debug - so we can see the length when we hit an assert in pushBack:
    size_t L = maxDeque.getLength();

    // accept new incoming sample - this corresponds to
    // Nayuki's Step 2 - "increment the array range’s right endpoint"
    while(!maxDeque.isEmpty() && greater(in, maxDeque.readTail()) )
      maxDeque.popBack();
    maxDeque.pushBack(in);

    // update delayline (forget oldest sample and remember current sample for later), this
    // corresponds to Nayuki's Step 3 - "increment the array range’s left endpoint":
    T oldest = delayLine.getSample(in);
    if(!maxDeque.isEmpty()) {            // happens when length is set to zero
      T maxVal = maxDeque.readHead();
      if(maxVal == oldest)
        maxDeque.popFront();
      return maxVal;
    }
    else
      return in;
  }

  /** Computes an output sample using a naive algorithm that scans the whole ringbuffer for its
  maximum value. The algorithm has a per sample complexity O(k) where k is the length of the
  filter. The implementation is mainly for testing purposes and should probably not be used in
  production code. */
  inline T getSampleNaive(T in)
  {
    delayLine.getSample(in);  // output of getSample not needed here
    T maxVal = delayLine.getMaximum(); // maybe pass greater function to the maximum finder
    return maxVal;
  }

  /** Resets the filter to its initial state. */
  void reset()
  {
    delayLine.reset();
    maxDeque.clear();
  }

protected:

  rsDelayBuffer<T> delayLine;
  rsDoubleEndedQueue<T> maxDeque;

  bool (*greater)(const T&, const T&) = &rsGreater;
  //std::function<bool(const T&, const T&)> greater; // = &rsGreater;

};

// maybe call it rsMovingSelector, movingExtremum ...an even more general concept would be a
// moving aggregator which could also be something like a moving median (or r-th quantile)
//
// https://www.nayuki.io/page/sliding-window-minimum-maximum-algorithm

//=================================================================================================

/** A variant of the rsMovingMaximumFilter that also extracts the minimum. It's slightly more
economic to extract min and max with a single filter object rather than using separate filters for
min and max (specifically, the delayline can be shared by both filters). */

template<class T>
class rsMovingMinMaxFilter : public rsMovingMaximumFilter<T>  // rename to rsMinMaxFilter
{

public:

  rsMovingMinMaxFilter(size_t maxLength) : rsMovingMaximumFilter<T>(maxLength), minDeque(maxLength) {}

  void setLessThanFunction(bool (*lessThan)(const T&, const T&)) { less = lessThan; }

  inline void getMinMax(T in, T* minVal, T* maxVal)
  {
    // update deque tails:
    //while(!this->maxDeque.isEmpty() && greater(in, this->maxDeque.readTail()) )  // old, fails on gcc
    //  this->maxDeque.popBack();
    while(!this->maxDeque.isEmpty() && less(this->maxDeque.readTail(), in) )    // new, seems ok
      this->maxDeque.popBack();
    this->maxDeque.pushBack(in);
    while(!this->minDeque.isEmpty() && less(in, this->minDeque.readTail()) )
      this->minDeque.popBack();
    this->minDeque.pushBack(in);

    // update delayline and init outputs:
    T oldest = this->delayLine.getSample(in);
    *maxVal = *minVal = in;

    // update deque heads:
    if(!this->maxDeque.isEmpty()) {
      *maxVal = this->maxDeque.readHead();
      if(*maxVal == oldest)
        this->maxDeque.popFront(); }
    if(!minDeque.isEmpty()) {
      *minVal = minDeque.readHead();
      if(*minVal == oldest)
        minDeque.popFront(); }
  }
  // hmm...maybe that shouldn't be inlined...we'll see

  void reset()
  {
    rsMovingMaximumFilter<T>::reset();
    minDeque.clear();
  }

protected:

  rsDoubleEndedQueue<T> minDeque;
  bool (*less)(const T&, const T&) = &rsLess;

};

//=================================================================================================

template<class T>
class rsSlewRateLimiterLinear2 // we already have another class of that name - clean up!
{                              // maybe use this here as baseclass for the other

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the maximum upward change from one sample to the next. */
  void setUpwardLimit(const T& newLimit)
  {
    upwardLimit = newLimit;
  }

  /** Sets the maximum downward change from one sample to the next. */
  void setDownwardLimit(const T& newLimit)
  {
    downwardLimit = newLimit;
  }

  /** Convenience function to set both, upward and downward limits, to the same value. */
  void setLimits(const T& newLimit)
  {
    setUpwardLimit(newLimit);
    setDownwardLimit(newLimit);
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  T getSample(T in)
  {
    y1 += rsClip(in-y1, -downwardLimit, upwardLimit);
    return y1;
  }

  void reset() { y1 = T(0); }

protected:

  T upwardLimit = RS_INF(T), downwardLimit = RS_INF(T);
  T y1 = T(0);  // previous output sample

};

//=================================================================================================

/** A class for applying a special kind of nonlinear smoothing algorithm. It is based on min/max
filtering over a certain number of samples and taking an adjustable linear combination of both and
then applying a linear slew rate limiter to the result where the limit on the slew rate is
adaptively updated according to the current values of min and max such that that it could go from
min to max in the same given number of samples (or some fraction or multiple thereof).

It is good for post-processing the raw output of an envelope follower to make it smoother. The raw
envelope follower output will typically show parasitic oscillations with the signal's frequency.
Setting up a smoother with a length equal to the cycle-length of the (enveloped) input wave will
give optimal smoothing for the signal - in an ideal situation, the oscillations in the detected
envelope will be completely flattened out. */

template<class T>
class rsMinMaxSmoother
{

public:

  rsMinMaxSmoother(size_t maxLength) : minMaxFilter(maxLength) {}


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the smoothing length in samples. */
  void setLength(size_t newLength)
  {
    minMaxFilter.setLength(newLength);
    updateSlewRateLimitingFactor();
  }

  /** Adjusts the scaling factor for the slew-rate. With a factor of one, the slew rate is set up
  such that a transition between the currently measured min and max can occur in between L samples
  where L is the length of the min-max filter. An amount of zero turns slew rate limitig off and
  amounts higher than one make the whole thing sluggish and are probably not useful. */
  void setSlewRateLimiting(T newAmount)
  {
    slewLimitingAmount = newAmount;
    updateSlewRateLimitingFactor();
  }

  /** Sets the mix between min and max which is used as output signal (before it goes into the
  adaptive, linear slew rate limiter) where 0 corresponds to min only, 1 corresponds to max only
  and 0.5 to the arithmetic mean between min and max (todo: maybe a generalized mean could be
  useful?). */
  void setMinMaxMix(T newMix)
  {
    mix = newMix;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  size_t getLength() const { return minMaxFilter.getLength(); }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Processes one sample at a time. */
  T getSample(T in)
  {
    T minVal, maxVal;
    minMaxFilter.getMinMax(in, &minVal, &maxVal);
    if(limitingFactor == RS_INF(T))
      slewLimiter.setLimits(limitingFactor); // avoids 0 * inf = NaN in case maxVal == minVal
    else
      slewLimiter.setLimits((maxVal-minVal) * limitingFactor);
    return slewLimiter.getSample( (1-mix)*minVal + mix*maxVal );
  }

  /** Resets the filter to its initial state. */
  void reset()
  {
    minMaxFilter.reset();
    slewLimiter.reset();
  }

protected:

  void updateSlewRateLimitingFactor()
  {
    limitingFactor = 1.0 / (slewLimitingAmount * minMaxFilter.getLength());
  }

  T limitingFactor = 0.0;
  T slewLimitingAmount = 1.0;
  T mix = 0.5; // mixing between minimum an maximum output, 0: min only, 1: max only
  rsMovingMinMaxFilter<T> minMaxFilter;
  rsSlewRateLimiterLinear2<T> slewLimiter;

};
