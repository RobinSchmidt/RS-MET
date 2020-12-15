#ifndef RAPT_QUANTILEFILTER_H
#define RAPT_QUANTILEFILTER_H

/** Class to exctract moving quantiles (such as the median) from a signal in realtime. If the
quantile runs over N samples, the filter takes O(log(N)) operations per sample. It achieves this
by using an rsDoubleHeap together with a circular buffer of indices into that double-heap. The
process that takes place in getSample is to replace the oldest sample in the double-heap with the
new incoming sample. The potential re-ordering of the heaps due to such an replacement is kept
track of by the circular buffer, such that it always points to the oldest sample in the
double-heap.  */

template<class T>
class rsQuantileFilterCore
{

public:


  rsQuantileFilterCore(int maxLength = 2)
  {
    rsAssert(maxLength >= 2, "A maxLength of at least 2 is needed");
    auto swapNodes = [&](Node& a, Node& b) { this->swapNodes(a, b); };
    dblHp.setSwapFunction(swapNodes);
    setMaxLength(maxLength);
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the maximum length of the filter. This may re-allocate memory. */
  void setMaxLength(int newMaxLength)
  {
    // newMaxLength = rsNextPowerOfTwo(newMaxLength); // use later for optimizing wraparounds
    // mask = newMaxLength-1;

    small.resize(newMaxLength);
    large.resize(newMaxLength);
    keyBuf.resize(newMaxLength);
    setLengthAndReadPosition(L, p, true);
  }
  // maybe try to optimize the memory usage - we actually just need one nodes array of length
  // newMaxLength of which one part is used for the min-heap and the other for the max-heap - but
  // maybe that would make modulating the length more difficult - if we use only one buffer, we may
  // have to move around more data, when the length changes (i guess) - we'll see - hmm maybe not

  /** Sets the length of the filter, i.e. the number of past samples that we look at. */
  void setLength(int newLength, bool hard = false)
  { setLengthAndReadPosition(newLength, p, hard); }

  /** Sets the read position in the sorted array of stored past values. This array does not exist
  literally but only conceptually (in the naive implementation in the prototypes section, this
  actually exists literally). In practice it's the largest of the small values, i.e. the front
  element of the min-heap of large alues in our double-heap. So what this function actually does is
  to update the sizes of the max-heap (of small values) and the min-heap (of large values) in the
  double-heap. But conceptually, think about it simply as the readout index in an array of sorted
  past values. */
  void setReadPosition(int newPosition, bool hard = false)
  { setLengthAndReadPosition(L, newPosition, hard); }

  /** Sets length and read-position at once. */
  void setLengthAndReadPosition(int newLength, int newPosition, bool hard = false);

  /** When a higher level class wants to implement non-integer readout positions, we need to do a
  linear interpolation between the values at sorted-array positions to the left and to the right of
  actual requested non-integer position. This sets the weight w for the value to the right, the
  weight for the value to the left is then 1-w. By default, w = 1, so we give full weight to the
  sample to the right which is the smallest of the large values. */
  void setRightWeight(T newWeight) { w = newWeight; }

  void setLengthAndQuantile(int newLength, T newQuantile, bool hard = false)
  {
    lengthAndQuantileToPositionAndWeight(newLength, newQuantile, &p, &w);
    setLengthAndReadPosition(newLength, p, hard);
  }

  /** Sets a pointer to a signal buffer that is driven by client code to facilitate artifact-free
  modulation of the length. If you leave this unassigned, you may see artifacts when the length
  of the filter is switched from a lower to a higher value, due to the fact that we don't know the
  older samples here and just have to make up some values heuristically to fill up the heaps. When
  the true old samples are available via such a buffer, these artifacts can be avoided. As said,
  client code is supposed to feed/drive this buffer because this will also facilitate sharing of
  such buffers between a bunch of parallel quantile filters or to create a highpass version. That
  means before calling setLength (and then getSample) on this object, client code should have
  called getSample on the buffer object. */
  void setDelayBuffer(rsDelayBuffer<T>* newBuffer) { sigBuf = newBuffer; }
  // maybe rename to setSignalBuffer, setDelayBuffer or setInputBuffer - it's used for other 
  // things, too


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the length of the filter in samples. */
  int getLength() const { return L; }

  /** Returns the readout position in the (conceptual, not actually existent) array of sorted past
  input samples. */
  int getReadPosition() const { return p; }

  /** Returns the maximum length (in samples) that the filter currently supports, i.e. the capacity
  of the allocated memory for the buffers and heaps. */
  int getMaxLength() const { return (int) small.size(); } // large.size() == buf.size();

  /** Returns the quantile that this filter produces. */
  T getQuantile() const { return (p-1+w) / (L-1); }
  // formula needs to be verified
  // maybe have also a setLengthAndQuantile function that sets up L,p,w from L,q - more convenient
  // for the user - i'm currently doing the required computations in the embedding class 
  // rsQuantileFilter, but maybe they should be moved into this class


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes an output sample from a given input sample x. */
  T getSample(T x)
  {
    storeInput(x);       // updates double-heap and key-buffer
    return readOutput(); // sample readout
  }
  // maybe move to cpp file (also the two functions that are called) - it's perhaps better to not
  // inline this to keep binaries small (it does enough stuff to make the calling-overhead probably
  // preferable over the inlining-bloat) - maybe do benchmarks

  /** Resets the filter into its initial state. */
  void reset();

  /** Produces a sample that would have been produced, if the length of the filter would be
  longer by one sample, i.e. L+1 instead of L. This can be used to implement non-integer length 
  filters at a higher level by crossfading between the outputs of two filters whose lengths are L 
  and L+1 without literally running a second filter of length L+1. The output of the L+1 filter is 
  simulated by doing some trickery. From the values returned by the regular getSample call and the 
  call to this function afterwards, a non-integer length filter sample can be computed by 
  crossfading. To use this feature, the input buffer (delayline) must be assigned and properly 
  driven by client code, because the x[n-L] sample is not in the heaps, so we must retrieve it from
  the delayline. */
  T getElongatedOutput()
  {
    rsAssert(this->sigBuf != nullptr, "To use this feature, the input buffer must be assigned.");
    T xL = (*this->sigBuf)[this->L];   // should be x[n-L], client code must assure this
    return getElongatedOutput(xL);
  }
  // maybe also implement a function getShortenedOutput that produces the sample that would have 
  // resulted when the filter would be 1 sample shorter. this has the advantage that it doesn't 
  // need to retrieve the older sample from the delayline...could we even go down to a filter of 
  // length 1 by crossfading between a length 2 filter and the original input? That would allow
  // a neutral setting of the filter...and it may also be simpler to implement...maybe

  /** After calling getSample, this function may be called to produce an output that getSample 
  would have produced when the length would have been one sample longer, i.e. L+1 instead of L and 
  at some time within this larger time interval, the value x would have been fed into the filter. 
  Used internally by getElongatedOutput() in which case the input sample from L samples ago is 
  passed as x (which is the first (most recent, newest) value which we don't have in the heaps 
  anymore). */
  T getElongatedOutput(T x)
  {
    int p1;                                                 // read position
    T w1, xS, xL;                                           // weight, xLarge, xSmall
    T q = getQuantile();
    lengthAndQuantileToPositionAndWeight(L+1, q, &p1, &w1);
    T S0 = getS0(), L0 = getL0();
    if(p1 == p) {                                           // additional slot is in the large heap
      T S1 = getS1(x);
      if(     x > L0) { xS = S0; xL = L0; }
      else if(x > S0) { xS = S0; xL = x;  }
      else if(x > S1) { xS = x;  xL = S0; }
      else            { xS = S1; xL = S0; } }
    else {                                                  // additional slot is in the small heap
      rsAssert(p1 == p+1);                                  // sanity check
      T L1 = getL1(x);
      if(     x < S0) { xS = S0; xL = L0; }
      else if(x < L0) { xS = x;  xL = L0; }
      else if(x < L1) { xS = L0; xL = x;  }
      else            { xS = L0; xL = L1; } }
    return (T(1)-w1)*xS + w1*xL;
  }
  // I get "unresolved external symbol" linker errors in visual studio (msc) when trying to move 
  // this into the cpp file. I think, it has to do with member functions only being instantiated 
  // when they are actually called - which this function obviously is (otherwise the linker would 
  // have nothing to complain about) - but maybe it's called in a compilation unit that is compiled 
  // too late...or something
  // -> move code to .cpp and fix this - the function is already there but commented out
  // maybe when we call it in rsQuantileFilter, the problem will disappear - we'll see.

  /** After calling getSample, this function may be called to produce an output that getSample 
  would have produced when the length would have been one sample shorter, i.e. L-1 instead of L. 
  This is complementary to getElongatedOutput and can also be used to implement filters with 
  non-integer length by configuring the filter to have a length floor(length)+1 and blend the 
  regular output the filter obtained from getSample with an output obtained from this function. 
  The advantage is that it doesn't need the delay buffer to be assigned. */
  T getShortenedOutput()
  {
    int p1;                                     // read position
    T w1, xS, xL;                               // weight, xLarge, xSmall
    T q = getQuantile();
    lengthAndQuantileToPositionAndWeight(L-1, q, &p1, &w1);

    int k = keyBuf[bufIdx];                     // heap-key of oldest sample xOld
    bool kInUpper = dblHp.isKeyInLargeHeap(k);  // indicates, if oldest sample is in upper heap
    if(kInUpper)
      k = dblHp.toLargeHeapIndex(k);

    if(p1 == p) {
      if(kInUpper) {
        if(k == 0)  { xS = getS0();  xL = getL1(0); }   // xOld == L0
        else        { xS = getS0();  xL = getL0();  }}  // xOld >  L0
      else         {  xS = getL0();  xL = getL1(0); } } // xOld <= S0
    else {
      rsAssert(p1 == p-1);                              // sanity check
      if(kInUpper) {  xS = getS1(0); xL = getS0();  }   // xOld >= L0
      else         {
        if(k == 0)  { xS = getS1(0); xL = getL0(); }    // xOld == S0
        else        { xS = getS0();  xL = getL0(); }}}  // xOld <  S0

    return (T(1)-w1)*xS + w1*xL;
  }
  // todo: benchmark getElongatedOutput vs getShortenedOutput and use the more efficient version to
  // implement non-integer lengths


  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  static void lengthAndQuantileToPositionAndWeight(int L, T q, int* p, T* w)
  {
    T P = q * (L-1);            // non-integer read position
    *w  = P - floor(P);         // weight for value right to P
    *p  = (int)floor(P) + 1;    // +1 because p is actually the value to the right of P
    if(*p > L-1) {              // quantile q == 1 (maximum) needs special care
      *p = L-1; *w = T(1);  }
  }
  // same linker error as for getElongatedOutput
  // optimize: avoid calling floor twice

protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Internals */

  /** Accepts a new input sample and updates our internal buffers accordingly. This will have
  (conceptually) the effect that the oldest input sample will be remvoved from the double-heap and
  the new input will be inserted (what actually happens is a replacement, but that's an
  implementation detail). */
  void storeInput(T x)
  {
    //rsAssert(isStateConsistent(), "inconsistent state");
    int k = keyBuf[bufIdx];        // heap-key of oldest sample
    int w = wrap(bufIdx + L);      // (write) index of new node in keyBuf
    keyBuf[w] = k;                 // store preliminary key (== old node's key) in kexBuf at w
    dblHp.replace(k, Node(x, w));  // replace the old node, reshuffles heaps and keyBuf
    bufIdx = wrap(bufIdx + 1);     // update position in circular buffer of keys
    //rsAssert(isStateConsistent(), "inconsistent state");
  }

  /** Produces the current ouput sample by reading out the largest of the small values and the
  smallest of the large values and forming a linear combination. Has been factored out from
  getSample because it's also needed for modulating the length when no delay-buffer of old samples
  is available (i.e. when sigBuffer == nullptr because the user has not assigned it). */
  T readOutput()
  {
    T yS = dblHp.getLargestSmallValue().value;   // smaller value
    T yL = dblHp.getSmallestLargeValue().value;  // larger value
    T y  = (T(1)-w)*yS + w*yL;                   // linear interpolation
    return y;
  }

  /** This is called from setLengthAndReadPosition when its "hard" parameter is false. Instead of
  setting up a new length and position immediately and doing a hard reset, this function adapts the
  sizes of the two heaps by removing or inserting data into the heaps and/or moving data between
  the small and large heap. If N is the number of nodes that need to be removed/moved/inserted,
  the cost of this function is O(N*log(N)). */
  void modulateLengthAndReadPosition(int newLength, int newPosition);

  /** Moves the first (smallest) value of the large heap over to the small heap, such that it
  becomes the new largest element of the small heap. This decreases the size of the large heap by
  one and increases the size of the small heap by one. The return value informs whether this was
  successful - it will fail when the large heap has only one element in it because both heaps must
  always contain at least one element. */
  bool moveFirstLargeToSmall();

  /** Analog to moveFirstLargeToSmall. */
  bool moveFirstSmallToLarge();

  /** Discards the oldest sample. This will shrink one of the two heaps, wherever the oldest sample
  happens to be. This will also shrink the circular buffer by one. */
  void discardOldestSample();

  /** Adds a sample to one of the heaps as new oldest sample. The sample to be added is either a
  true older sample (if the user has assigned our past-signal-values buffer via
  setModulationBuffer and correctly keeps it up to date) or else, we will just heuristically invent
  a sample value to insert (choosing the most recent output value for that). This will also grow
  the circular buffer by one. */
  void addOlderSample();

  /** Wraparound for the bufIdx (handles only forward wraparound, i.e. input should be
  non-negative). */
  int wrap(int i) { rsAssert(i >= 0); return i % (int)keyBuf.size(); }
  // todo: optimize using a bitmask -> restrict buffer-sizes to powers of two

  /** Returns either the 2nd largest value of the small values (in the typical case) or - in the 
  edge case where the length of the small heap is 1 such that it doesn't have such a value - the 
  passed input value x. It's used in getElongatedOutput to handle these edge cases. I'm 
  not totally sure, why it works to return x in these cases, but i think, it's because of the 
  conditionals that are done in this function. (ToDo: figure out and explain) */
  T getS1(T x) const
  {
    if(this->dblHp.isSmallSizeOne()) return x;
    return this->dblHp.get2ndLargestSmallValue().value;
  }

  /** Analogous to get2ndLargestSmallOrX */
  T getL1(T x) const
  {
    if(this->dblHp.isLargeSizeOne()) return x;
    return this->dblHp.get2ndSmallestLargeValue().value;
  }

  T getS0() const { return dblHp.getLargestSmallValue().value;  }
  T getL0() const { return dblHp.getSmallestLargeValue().value; }

  
  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  /** A node stores an incoming signal value together with its index in the circular buffer. The
  circular buffer, in turn, stores the double-heap key of the node. Via these links, we can refer
  back and forth from buffer-indices and node-keys. The "less-than" comparison is based on the
  signal value. */
  struct Node
  {
    T value = T(0);
    int bufIdx = 0;   // maybe use size_t for better compatibility with rsDelayBuffer

    Node(T v = T(0), int i = 0) { value = v; bufIdx = i; }

    bool operator<(const Node& b) const
    { return this->value < b.value; }

    bool operator==(const Node& b) const
    { return this->value == b.value && this->bufIdx == b.bufIdx; }
  };

  /** The swapping for nodes must do the actual swap of a and b as usual but also let the circular
  buffer keep track of what gets swapped. */
  void swapNodes(Node& a, Node& b)
  {
    rsSwap(a, b);
    rsSwap(keyBuf[a.bufIdx], keyBuf[b.bufIdx]);
  }

  std::vector<Node>  small, large;     // storage arrays of the nodes
  rsDoubleHeap<Node> dblHp;            // maintains large/small as double-heap
  std::vector<int>   keyBuf;           // circular buffer of heap keys
  rsDelayBuffer<T>*  sigBuf = nullptr; // (possibly shared) buffer of delayed input samples

  int bufIdx = 0;    // index into keyBuf, mayb rename to keyIdx
  int L      = 2;    // total length of filter
  int p      = 1;    // readout position in imagined sorted array, 1 <= p <= L-1
  T   w      = T(1); // weight for smallest large value in the linear interpolation
  // Maybe use size_t instead of int. That would be consistent with rsMovingMaximumFilter and
  // better compatible with rsDelayBuffer - but: size_t is 64 bit and int only 32, so int has lower
  // memory consumption...so maybe not..but wait: we dont use arrays of this type, so that may be
  // irrelevant. For the heap-keys, we'll use int anyway.


  //-----------------------------------------------------------------------------------------------
  /** \name Scaffolding. For self-tests for debugging. May be removed when class is finished.
  ...but maybe let's keep them because they are nice for documenting the algorithm as well. */

  bool isStateConsistent();               // O(N)
  bool isNodeConsistent(const Node& n);   // O(N)
  bool isBufferSlotConsistent(int i);     // O(1)

};

//=================================================================================================

/** Subclass of rsQuantileFilterCore that facilitates smoother sweeps of the filter length by 
supporting non-integer filter lengths. These are implemented by means of crossfading between a 
filter of length L and L+1, where L is the floor of the desired length and the crossfade is done
via its fractional part. The filter of length L+1 is not actually a second filter - instead we 
simulate such a second filter with the heaps we already have and some additional trickery. It 
supports also filter lengths less than 2. This is implemented by crossfading between a length 2
quantile filter and the input. */

template<class T>
class rsQuantileFilterCore2 : protected rsQuantileFilterCore<T>
{

public:

  using Base = rsQuantileFilterCore<T>;

  /** Sets the new length and quantile. The length can be a non-integer number >= 1 (smaller 
  values are clipped from below at 1). The quantile should be a number between 0 and 1 (both 
  ends inclusive). */
  void setLengthAndQuantile(T L, T q)
  {
    rsAssert(q >= T(0) && q <= T(1), "Quantiles must be values between 0 and 1.");
    if(L < T(2)) {
      L = rsMax(L, T(1));     // 1 is really the lower limit for the length
      blend = L - floor(L);
      frac  = q;
      Base::setLengthAndQuantile(2, q); }
    else {
      blend = T(1);
      frac  = L - floor(L);
      Base::setLengthAndQuantile((int)L, q); }
  }

  /** Prodcues one sample at a time. */
  T getSample(T x)
  {
    T y;
    if(this->sigBuf) {
      T y0 = Base::getSample(x);          // output of filter of length L
      T y1 = Base::getElongatedOutput();  // output of filter of length L+1
      y = (T(1)-frac)*y0 + frac*y1; }     // crossfade between length L and L+1
    else
      y = Base::getSample(x);
    return (T(1)-blend)*x + blend*y;      // crossfade betwenn input and output for L < 2 case
  }

  // some delegations to the basclass:
  void setMaxLength(int newMaxLength) { Base::setMaxLength(newMaxLength); }
  void setDelayBuffer(rsDelayBuffer<T>* newBuffer) { Base::setDelayBuffer(newBuffer); }
  void reset() { Base::reset(); }


protected:

  // additional member variables to avoid recomputation of them in readOutputWithOneMoreInput

  // algo parameters:
  T frac  = T(0);  // fractional part of length
  T blend = T(1);  // blend between filter output and input (1.0: only output)
  //int p1;      // readout position used in getElongatedOutput
  //T w1;        // weight used in getElongatedOutput
  // cached w1, p1 values can later be used for optimization purposes (they need to be recomputed
  // only when the settings change, not necessarily per sample as is currently done)

};

//=================================================================================================

/** This is a more user-friendly, audio-oriented, plugin-ready wrapper around the core algorithm
of the moving quantile filter. It uses a more intuitive, audio-friendly parametrization and also
adds features such as highpass mode, .... */

template<class T>
class rsQuantileFilter
{

public:

  rsQuantileFilter() : dirty(true)
  {
    allocateResources();
    core.setDelayBuffer(&delayLine);
    //dirty = true;
    //updateInternals(); // hangs - why?
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the sample rate at which this filter should operate. May re-allocate memory. */
  void setSampleRate(T newSampleRate)
  {
    sampleRate = newSampleRate;
    allocateResources();
    dirty = true;
  }

  /** Sets the maximum length (in seconds) for this filter. May re-allocate memory. */
  void setMaxLength(T newMaxLength)
  {
    maxLength = newMaxLength;
    allocateResources();
    dirty = true;
  }

  /** Sets sample rate and maximum length at the same time. May re-allocate memory. This may avoid
  some re-allocations compared to using setSampleRate followed by setMaxLength (or vice versa), 
  depending on the situation - so, if possible, it's recommended to set both at the same time. */
  void setSampleRateAndMaxLength(T newSampleRate, T newMaxLength)
  {
    sampleRate = newSampleRate;
    maxLength  = newMaxLength;
    allocateResources();
    dirty = true;
  }

  /** Sets the length of the filter in seconds, i.e. the time interval of the samples that we
  look at. */
  void setLength(T newLength) { length = newLength; dirty = true; }

  /** Sets some sort of pseudo cutoff frequency which we just define as the reciprocal of the
  length. The filter is actually very nonlinear, so the notion of a cutoff frequency or even just
  of a frequency response does not really apply. Instead, if you consider a periodic waveform and
  pass it through a quantile filter whose length matches the wave's period, the output will be a
  constant. This can be understood from the fact, that the array of sorted input values will be
  the same, no matter, at which phase you pick the starting point of the cycle. (todo: verify this
  empricially) */
  void setFrequency(T newFrequency) { setLength(T(1) / newFrequency); }

  /** Sets the quantile that should be extracted. 0.0 gives the minimum, 1.0 gives the maximum,
  0.5 gives the median, 0.25 the lower quartile, 0.75 the upper quartile, etc. */
  void setQuantile(T newQuantile1) { quantile = newQuantile1; dirty = true; }

  /** Sets the gain (as raw scale factor) for the actual filter output (which is lowpass'ish in
  character, which is why we call this parameter lowpass gain). */
  void setLowpassGain(T newGain) { loGain = newGain; }

  /** Sets the gain for a highpass signal that is obtained by subtracting the (lowpass) filter
  output from an appropriately delayed input signal. Highpass and lowpass gain can be used to
  obtain different response types. The default is lo = 1, hi = 0 giving a lowpass filter,
  lo = 0, hi = 1 gives a highpass filter, lo = 1, hi = 1 gives a pure delay, lo = 1, hi = 2 gives
  a high frequency boost, etc. */
  void setHighpassGain(T newGain) { hiGain = newGain; }

  // void setDelayScaler(T newScaler) { delayScl = newScaler; }
  // typical range should be 0..2....i'm not yet sure, if that's useful as a user parameter...

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the delay that this filter introduces (in samples). */
  T getDelayInSamples() const { return delay; }
  // hmm...actually, the returned value is not yet up to date when this is called before getSample
  // because getSample may recalculate the delay before actually producing a sample - so if client
  // code calls this before getSample, it will always work with a value that is lagging behind
  // by one sample...hmmm...if client code wants to avoid this, it could call updateInternals
  // itself before calling getSample...but that's very inconvenient, API wise
  // maybe we should add to the documentation, that the returned value applies to the sample that
  // *was* produced most recently, not to the sample that *will be* produced next


  T getFrequency() const { return T(1) / length; }

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Produces one output sample from a given input sample. */
  T getSample(T x)
  {
    delayLine.getSample(x);
    // Must be called *before* updateInternals because in updateInternals, we may modulate the
    // length and the artifact-avoidance strategy for the modulation assumes that we have already
    // called getSample on the delay buffer. At this point, we are not (yet) interested in the
    // output of the delayline - we will retrieve its ouput later via random access using the []
    // operator. Actually, it would not make sense to try to retrieve an output sample before
    // calling updateInternals anyway, because that call also updates the delay which we need to
    // read out the delayline at the correct position.

    if(dirty) updateInternals();
    T yL = core.getSample(x);               // lowpass part
    T yH = delayLine[delayScl*delay] - yL;  // highpass part
    return loGain * yL + hiGain * yH;
  }
  // maybe factor out a function to produce lowpass and highpass getSampleLoHi or something at the
  // same time - client code may find that useful - or maybe getOutputs to be consistent with
  // rsStateVariableFilter

  /** Resets the filter into its initial state. */
  void reset()
  {
    core.reset();
    delayLine.reset();
    //y = T(0);
  }

  /** Updates the internal algorithm parameters and embedded objects according to the user
  parameters. This is called in getSample, if the state is dirty but sometimes it may be
  convenient to call it from client code, too. */
  virtual void updateInternals()
  {
    // compute internal and set up core parameters:
    double L = length*sampleRate;  // length in samples
    core.setLengthAndQuantile(L, quantile);
    delay = T(0.5)*(L-1);
    dirty = false;
  }


protected:

  /** Computes filter algorithm parameters length L, readout point p (both in samples), weight
  w for the linear interpolation and the delay d from user parameters length (in seconds), quantile
  (normalized 0..1) and sampleRate. */
  //static void convertParameters(T length, T quantile, T sampleRate, int* L, int* p, T* w, T* d);
  // obsolete - delete soon

  /** Returns the maximum length in samples that may be needed for the delayline and core 
  buffers. */
  int getMaxRequiredLengthInSamples()
  {
    return (int) ceil(maxLength * sampleRate); // maxLength is a user parameter in seconds
  }

  /** Allocates the memory used for the delay-buffers, heaps, etc. */
  virtual void allocateResources()
  {
    int mL = getMaxRequiredLengthInSamples();
    core.setMaxLength(mL);
    delayLine.setCapacity(mL);
  }

  // user parameters:
  //T feedback   = 0.0; // ...might be interesting to experiment with later
  T sampleRate = 44100; // in Hz
  T maxLength  = 0.1;   // in seconds
  T length     = 0.01;  // in seconds
  T quantile   = 0.5;   // in 0..1, 0.5 is median
  T loGain     = 1.0;   // gain for the filter output (which kinda lowpass)
  T hiGain     = 0.0;   // gain for delayed input minus filter output (kinda highpass)
  T delayScl   = 1.0;   // scaler for the delay for the input signal

  // filter state, algo parameters, infrastructure:
  //T y     = 0.0;   // previous output sample (for feedback)
  T delay = 0.0;   // required delay to implement highpass by subtraction
  std::atomic<bool> dirty = true;  // flag to indicate that algo params must be recalculated

  // embedded objects:
  rsQuantileFilterCore2<T> core;
  rsDelayBuffer<T> delayLine;

};

// Notes:
// I'm trying to make the class thread-safe by itself in a typical plugin scenario by using an
// atomic flag to indicate that internal algorithm parameters need to be re-calculated, which is
// then done in getSample (supposed to be called on the audio thread). However, this does not apply
// to the loGain and hiGain variables - they are set directly. However, if the type T is a type
// for which load and store operations actually are atomic, modifying them (from a different
// thread, like the gui thread) should be safe, too. I'm not yet sure what to do - if i use
// std::atomic<T>, it may lock mutexes which is also undesirable...so, even though i attempt to
// give thread safety, at the time being, it's better to assume, it's not and client code should
// implement its own synchronization scheme, if it deems it necessarry. See:
// https://stackoverflow.com/questions/36624881/why-is-integer-assignment-on-a-naturally-aligned-variable-atomic-on-x86/36685056#36685056




#endif
