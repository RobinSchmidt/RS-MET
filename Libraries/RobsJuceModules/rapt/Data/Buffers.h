#pragma once

/** Implements (the basis for) a circular buffer */

template<class T>
class rsRingBuffer
{

public:

  /** Creates a buffer with a given initial capacity. The actual capacity will be the next power of
  two of given value. */
  rsRingBuffer(size_t capacity);

  /** Sets a new capacity for this buffer. The actual capacity will be the next power of two of
  given value. */
  void setCapacity(size_t newCapacity);

  /** Initializes all the buffer elements with given value (default is zero). */
  void initBufferValues(T value = T(0));

  /** Returns the capacity of this buffer. */
  size_t getCapacity() const { return data.size(); }
  // is this correct? or do we need to subtract 1? no - i don't think so. -> do unit test!


  std::vector<T> data;
  // temporarily moved to public - but maybe it should really be public

protected:

  /** Wraparound */
  inline size_t wrap(size_t i) const { return i & mask; }

  size_t mask;

};
// todo:
// -maybe implement [] operator as: return data[i & mask]

//=================================================================================================

/** Subclass of rsRingBuffer that adds some features to make it more convenient to use as
delay-buffer or delay-line in signal processing applications. */

template<class T>
class rsDelayBuffer : public rsRingBuffer<T>
{

public:


  rsDelayBuffer(size_t capacity = 0) : rsRingBuffer<T>(capacity) {}
  // todo: make a default constructor creating a buffer with capacity 0 and provide a setCapacity
  // function



  /** Sets up a new buffer length and adjusts the left index accordingly, leaving the right index
  where it is. The reason to do it that way and not the other way around is that in a delayline,
  the right index represents the write-head and the left index represents the read-head and on a
  change of the delay time, we want to move the read-head. */
  void setLength(size_t newLength)
  {
    rsAssert(newLength <= this->getCapacity(), "Desired length exceeds capacity");
    // maybe automatically increase capacity in such a case - maybe do this optionally controlled
    // by a 2nd boolean parameter "increaseCapacityIfNeeded" or "reallocateIfNeeded" - it should
    // clearly communicate that memory allocation may take place

    //length = newLength; // old
    length = rsMin(newLength, this->getCapacity()); // new
    adjustReadIndex();
  }
  // maybe rename to setDelay

  inline size_t getLength() const { return length; }



  inline T getSample(T in)
  {
    // We increment the indices before the actual read/write operations to allow further operations
    // after getSample (where the indices must be still valid):
    incrementIndices(1);
    // BUT: if client code actually relies on this, this points to a design flaw. We should
    // provide a method to modify the most recently stored value - we already have such a method:
    // setNewest. Client code should be forced to used this, if it wants to modify the stored
    // value after calling getSample. Then we are free here to use pre- or post-increment however
    // we want without affecting client code - we just would have to adapt getNewest if we switch
    // to post-increment

    // The actual read/write operations
    this->data[writeIndex] = in;           // right index is write index
    T out = this->data[readIndex];         // left index is read index
    return out;
  }
  // rename to pushRightPopLeft maybe make a subclass delayline that defines an alias function
  // name getSample...any maybe merge class with rsDoubleEndedQueue

  /** Advances our read- and write pointers by the given number of positions. */
  inline void incrementIndices(size_t amount = 1)
  {
    writeIndex = this->wrap(writeIndex + amount);
    readIndex  = this->wrap(readIndex  + amount);
  }
  // maybe rename to incrementIndices, make protected

  /*
  inline void retractPointers(size_t amount = 1)
  {
    if(amount > rightIndex) rightIndex += getCapacity();
    if(amount > leftIndex)  leftIndex  += getCapacity();
    rightIndex = this->wrap(rightIndex - amount);
    leftIndex  = this->wrap(leftIndex  - amount);
  }
  // needs test
  */


  // Return samples from the buffer without updating anything

  T getOldest() const { return this->data[this->wrap(readIndex+1)]; }  // why +1?
  T getNewest() const { return this->data[writeIndex];        }


  /** Sets/replaces the newest sample stored in the buffer. If you need to modify the most recently
  stored sample after you have already called getSample(), this function is your friend. */
  void setNewest(T x) { this->data[writeIndex] = x; }

  //void setOldest(T x) { this->data[wrap(leftIndex+1)] = x } // seems useless

  /** Sets the sample that should be returned in the next call to getSample or in the delay-th call
  after that. */
  void setNextReturnSample(T x, size_t delay = 0) { this->data[this->wrap(readIndex+delay+1)] = x; }
  // the +1 is due to the use of pre-increment in getSample

  /** Sets the object up in such a way that the in the call to getSample, data[0] will be
  returned. */
  void setNextReadIndex(size_t i)
  {
    if(i > 0) readIndex = this->wrap(i-1);
    else      readIndex = this->data.size()-1; // == mask?
    adjustWriteIndex();
  }



  size_t getReadIndex()  const { return readIndex; }
  size_t getWriteIndex() const { return writeIndex; }


  size_t getIndexFromOldest(size_t i) const
  {
    // is this correct? -> make unit tests:
    //if(i > length)    i -= length;  // should never happen
    //if(i < leftIndex) i += getCapacity();

    return this->wrap(readIndex + i + 1); // why do we need +1?
  }
  // needs test - maybe we have to do something similar as in getIndexFromNewest ..maybe
  // if(i < leftIndex) i += getCapacity()  ?
  // and/or maybe we have to take the length into account? it seems to work only when
  // length == capacity
  // but actually, we do the i += stuff below only to avoid negative numbers when we do
  // rightIndex-i - but negative number cannot happen here because we do not subtract

  /** Returns the index j into our internal data array such that data[j] contains the value from
  i samples ago. */
  size_t getIndexFromNewest(size_t i) const
  {
    if(i > writeIndex)
      i += this->getCapacity();
    // do this branchless: i += (i > rightIndex) * getCapacity();
    // is this even needed or will underflow in the subtraction below produce the correct result as
    // well? -> figure out and test. we do it to get negative numbers in the rightIndex-i operation

    return this->wrap(writeIndex - i);
  }
  // maybe make protected



  T& fromNewest(size_t i) { return this->data[getIndexFromNewest(i)]; }
  T& fromOldest(size_t i) { return this->data[getIndexFromOldest(i)]; }

  // readout at given delay d
  T operator[](T d)
  {
    size_t i = (size_t) d;
    T f  = d-i;
    T y0 = fromNewest(i);
    T y1 = fromNewest(i+1);
    return (1-f)*y0 + f*y1;  // lerp
  }



  /** Writes the content of this circular buffer into the given linear buffer in either
  chronological or reverse chronological order. */
  void copyTo(T* buffer, bool newestFirst)  // maybe provide false as default
  {
    if(newestFirst)
      for(size_t i = 0; i < getLength(); i++)
        buffer[i] = fromNewest(i);
    else
      for(size_t i = 0; i < getLength(); i++)
        buffer[i] = fromOldest(i);
  }


  /** Returns the maximum value in the range between the two pointers. Mostly for testing
  purposes (the moving-max filter) - not efficient, takes time O(L) where L is the length. */
  T getMaximum()
  {
    size_t i = writeIndex;
    T ex = this->data[i];
    while(i != readIndex) {
      i = this->wrap(i-1);
      if(this->data[i] > ex)  // or >= or maybe use a comparison function and rename function to
        ex = this->data[i];   // getOptimum/Extremum
    }
    return ex;
  }

  void reset();

protected:

  /** Updates the current left index according right index and length. */
  void adjustReadIndex()  { readIndex = this->wrap(writeIndex - length); }
  // rename to updateReadIndex or adjustReadIndex - does the formula always work, even when
  // writeIndex - length is nagtive? we are using unsigned integers here

  void adjustWriteIndex()  { writeIndex = this->wrap(readIndex + length); }

  /** Updates the current right index according left index and length. */
  //void updateRightIndex() { rightIndex = wrap(leftIndex + length); }
  // maybe uncomment if needed - it's not typically used

  size_t writeIndex = 0, readIndex = 0;
  size_t length = 0; // rename to capacity ..or use data.size() - nope it's the current length
  // maybe name them generally rightEnd, leftEnd or something ... rgt, lft. L,R - then we may
  // also make rsDoubleEndedQueue a subclass and inherit the data and wrap function

  // maybe use w,r,L for writeIndex, readIndex, length


  //template<class U> friend class rsMovingQuantileFilter; // try to get rid

};

//=================================================================================================

/** Implements the double-ended queue data structure. This implementation provides a static, finite
capacity that has to be passed at construction time. Due to implementation details, the actual
capacity is always a power of two minus one, so when you pass an arbitrary number to the
constructor that is not of the form 2^k - 1 the smallest possible k will be chosen such that
2^k - 1 is greater than or equal to your requested capacity.

The implementation is based on a circular buffer. Wikipedia says, the C++ std::deque class uses a
multiple-array implementation which is probably not so good for realtime purposes:
https://en.wikipedia.org/wiki/Double-ended_queue */

template<class T>
class rsDoubleEndedQueue : public rsRingBuffer<T>
{

public:

  /** Constructor. Allocates enough memory for an internal buffer to hold the "capacity" number of
  queued values. The memory allocated for the buffer will be the smallest power of two that is
  greater or equal to the desired capacity plus 1 (an additional index value is required to make
  the index arithmetic work right). */
  rsDoubleEndedQueue(size_t capacity) : rsRingBuffer<T>(capacity+1) {}

  //-----------------------------------------------------------------------------------------------
  /** \name Data Access */

  /** Appends a value to right/head side of the queue. */
  inline void pushFront(T value)
  {
    RAPT::rsAssert(!isFull(), "Trying to push onto full deque");
    this->data[head] = value;
    head = this->wrap(head+1);
  }

  /** Appends a value to left/tail side of the queue. */
  inline void pushBack(T value)
  {
    RAPT::rsAssert(!isFull(), "Trying to push onto full deque");
    this->data[tail] = value;
    tail = this->wrap(tail-1);
  }

  /** Removes a value from the right/head side of the queue and returns it. */
  inline T popFront()
  {
    RAPT::rsAssert(!isEmpty(), "Trying to pop from empty deque");
    head = this->wrap(head-1);
    return this->data[head];
  }

  /** Removes a value from the left/tail side of the queue and returns it. */
  inline T popBack()
  {
    RAPT::rsAssert(!isEmpty(), "Trying to pop from empty deque");
    tail = this->wrap(tail+1);
    return this->data[tail];
  }

  /** Returns the value at the head of the queue, i.e. the rightmost value. */
  inline T readHead() const
  {
    RAPT::rsAssert(!isEmpty(), "Trying to read from empty deque");
    return this->data[this->wrap(head-1)];
  }

  /** ...UNDER CONSTRUCTION...
  Returns a linearly interpolated value between the rightmost value and its immediate left 
  neighbor using the given interpolation coefficient in the range [0,1]. 0 means rightmost only, 
  1 means left neighbor only. This is used to facilitate non-integer lengths in min/max filters. */
  /*
  inline T readHeadInterpolated(T coeff = T(0)) const
  {
    RAPT::rsAssert(getLength() >= 2, "Deque too short for interpolation");
    return (1-coeff) * this->data[this->wrap(head-1)] + coeff * this->data[this->wrap(head-2)];
  }
  */
  // needs test, todo: implement a similar function readTailInterpolated

  /** Returns the value at the tail of the queue, i.e. the leftmost value. */
  inline T readTail() const
  {
    RAPT::rsAssert(!isEmpty(), "Trying to read from empty deque");
    return this->data[this->wrap(tail+1)];
  }
  // or maybe readFirst/Last Front/Back or peekFront/peekBack, maybe also provide
  // writeFront/writeBack

  /** Clears the queue and optionally resets the data in the underlying buffer to all zeros. The
  clearing of the queue itself just resets some pointers/indices, turning the data into
  inaccessible garbage - but it should have no effect for the functionality of the queue whether
  this data buffer is uninitialized/garbage or all zeros. */
  void clear(bool clearGarbage = false)
  {
    head = 1;
    tail = 0;
    if(clearGarbage)
      this->initBufferValues(0);
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the number of values that are currently in the queue. */
  inline size_t getLength() const { return this->wrap(head - tail - 1); }

  /** Returns the maximum allowed length for the queue. Due to the way, the head- and tail pointers
  operate, the usable length is 2 less than the size of the underlying data buffer. */
  inline size_t getMaxLength() const { return this->data.size()-1; } //

  /** Returns true if the queue is empty. */
  inline bool isEmpty() const { return getLength() == 0; }

  /** Returns true, if the queue is full. */
  inline bool isFull() const { return getLength() >= getMaxLength(); }


protected:

  size_t head = 1, tail = 0; // a.k.a front/back, first/last
  // the actual data is stored at data[wrap(head-1)] and data[wrap(tail+1)]

};

/*
ToDo:
-test it more thoroughly when it's completely full - there are some hickups with the moving-max
 filter when it operates at full capacity. at capacity-1, it still works well. but the problem
 is more likely the moving-max filter itself
-maybe provide push-functions that take a vector argument and pop-functions that pop a range of
 values at once
-maybe provide some sort of random access functions:
 T fromFront(size_t i) where i = 0 returns the front element, 1 the one after the front, etc.
 or fromHead - analogously for fromBack/fromTail/fromRight
-maybe binary index search functions (at which index, counted from front (or back) does the value
 in the que get greater than some value
 size_t searchFromFront...this may become relevant for optimizing the moving maximum filter
-maybe factor out an rsDoubleEndedQueueView that just stores a pointer to an array somewhere
 ..but may that's not so useful after all
*/

// maybe make convenience subclasses rsQueue and rsStack - both datastructures have different
// subsets of functionality of the double ended queue actually, but the convenience functions could
// have more stacky or queuey names like push/pop/top, enqueue/dequeue/front/back - they should
// use protected inheritance and use delegation, like Stack::push delegates to Base::pushFront,
// etc.




