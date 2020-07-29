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
  literally but only conceptually (in the naive implementation, this actually exists literally). In 
  practice it's the largest of the small values, i.e. the front element of the min-heap of large 
  values in our double-heap. So what this function actually does is to update the sizes of the 
  max-heap (of small values) and the min-heap (of large values) in the double-heap. But 
  conceptually, think about it simply as the readout index in an array of stored past values. */
  void setReadPosition(int newPosition, bool hard = false) 
  { setLengthAndReadPosition(L, newPosition, hard); }

  /** Sets length and read-position at once. */
  void setLengthAndReadPosition(int newLength, int newPosition, bool hard = false);

  /** When a higher level class wants to implement non-integer readout positions, we need to do a 
  linear interpolation between the values at sorted-array positions to the left and to the right of 
  actual requested non-integer position. This sets the weight w for the value to the right, the 
  weight for the value to the left is the 1-w. By default, w = 1, so we give full weight to the 
  sample to the right which is the smallest of the large values. */
  void setRightWeight(T newWeight) { w = newWeight; }

  /** Sets a pointer to a signal buffer that is driven by client code to facilitate artifact-free
  modulation of the length. If you leave this unassigned, you may see artifacts when the length
  of the filter is switched from a lower to a higher value, due to the fact that we don't know the 
  older samples here and just have to make up some values to fill up the heaps. When the true old 
  samples are available via such a buffer, these artifacts can be avoided. As said, client code 
  is supposed to feed/drive this buffer because this will also facilitate sharing of such buffers 
  between a bunch of parallel quantile filters or to create a highpass version. That means before
  calling setLength (and then getSample) on this object, client code should have called getSample
  on the buffer object. */
  void setModulationBuffer(rsDelayBuffer<T>* newBuffer) { sigBuf = newBuffer; }


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


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes an output sample from a given input sample x. */
  T getSample(T x)
  {
    storeInput(x);       // updates double-heap and key-buffer
    return readOutput(); // sample readout
  }

  /** Resets the filter into its initial state. */
  void reset()
  {
    for(int n = 0; n < L; n++) {
      dblHp.atIndex(n).value  = T(0);
      int k = dblHp.indexToKey(n);
      dblHp.atIndex(n).bufIdx = n;
      keyBuf[n] = k; }
    bufIdx = 0;
  }


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

  /** Adds a (dummy) sample to one of the heaps as new oldest sample. This will also grow the 
  circular buffer by one. */
  void addOlderSample();

  /** Conversion from a coneceptual buffer index in the range 0..L-1 to the actual index in the 
  range 0..buf.size()-1 */
  //int convertBufIndex(int indexFromOldest)
  //{ return (indexFromOldest + bufIdx) % (int)buf.size(); }
  // we don't really need this anymore - this was meant for scaffold code

  /** Wraparound for the bufIdx (handles only forward wraparound, i.e. input should be 
  non-negative). */
  int wrap(int i) { rsAssert(i >= 0); return i % (int)keyBuf.size(); }
  // todo: optimize using a bitmask -> restrict buffer-sizes to powers of two
 

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
  std::vector<int>   keyBuf;           // circular buffer of heap keys - rename to keyBuf
  rsDelayBuffer<T>*  sigBuf = nullptr; // (possibly shared) buffer of delayed input samples

  int bufIdx = 0;    // index into keyBuf, mayb rename to keyIdx
  int L      = 2;    // total length of filter
  int p      = 1;    // readout position, 1 <= p <= L-1
  T   w      = T(1); // weight for smallest large value in the linear interpolation
  // maybe use size_t instead of int


  //-----------------------------------------------------------------------------------------------
  /** \name Scaffolding (may be removed when class is finished) */

  // some scaffolding code for development:
  //int getNodeKey(Node node);  // O(N)
  //void fixInconsistentBufferKeys(Node startNode);
  //void makeBufferConsistent();

  // self-tests for debugging:
  bool isStateConsistent(); 
  bool isNodeConsistent(const Node& n);
  bool isBufferSlotConsistent(int i);
  // maybe don't delete them, these are nice for documentation purposes

};


#endif