#pragma once

template<class T>
class rsBuffer
{

public:

  /** Creates a buffer with a given initial capacity. The actual capacity will be the next power of
  two of given value. */
  rsBuffer(size_t capacity);

  /** Sets a new capacity for this buffer. The actual capacity will be the next power of two of 
  given value. */
  void setCapacity(size_t newCapacity);

  /** Initializes all the buffer elements with given value (default is zero). */
  void initBufferValues(T value = T(0));

  /** Returns the capacity of this buffer. */
  size_t getCapacity() const { return data.size(); }
  // is this correct? or do we need to subtract 1? -> do unit test!


  std::vector<T> data;
  // temporarily moved to public - but maybe it should really be public

protected:

  /** Wraparound */
  inline size_t wrap(size_t i) const { return i & mask; }

  size_t mask;

};
// maybe this class should be named rsRingBuffer or rsCircularBuffer or rsRingBufferBase and the
// class below should be named rsDelayLine or rsDelayBuffer

//=================================================================================================

template<class T>
class rsRingBuffer : public rsBuffer<T>
{

public:


  rsRingBuffer(size_t capacity = 0) : rsBuffer<T>(capacity) {}
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

  T getOldest() const { return this->data[wrap(readIndex+1)]; }  // why +1?
  T getNewest() const { return this->data[writeIndex];        }


  /** Sets/replaces the newest sample stored in the buffer. If you need to modify the most recently
  stored sample after you have already called getSample(), this function is your friend. */
  void setNewest(T x) { this->data[writeIndex] = x; }

  //void setOldest(T x) { this->data[wrap(leftIndex+1)] = x } // seems useless

  /** Sets the sample that should be returned in the next call to getSample or in the delay-th call
  after that. */
  void setNextReturnSample(T x, size_t delay = 0) { this->data[wrap(readIndex+delay+1)] = x; }
  // the +1 is due to the use of pre-increment in getSample





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
      i += getCapacity();  
    // do this branchless: i += (i > rightIndex) * getCapacity();
    // is this even needed or will underflow in the subtraction below produce the correct result as
    // well? -> figure out and test. we do it to get negative numbers in the rightIndex-i operation

    return this->wrap(writeIndex - i);
  }
  // maybe make protected



  T& fromNewest(size_t i) { return data[getIndexFromNewest(i)]; }
  T& fromOldest(size_t i) { return data[getIndexFromOldest(i)]; }





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

template<class T>
bool rsGreater(const T& a, const T& b);

template<class T>
bool rsLess(const T& a, const T& b);  // merge with function in RAPT...SortAndSearch

/** Implements the double-ended queue data structure. This implementation provides a static, finite
capacity that has to be passed at construction time. Due to implementation details, the actual
capacity is always a power of two minus two, so when you pass an arbitrary number to the
constructor that is not of the form 2^k - 2 the smallest possible k will be chosen such that
2^k - 2 is greater than or equal to your requested capacity.

The implementation is based on a circular buffer. Wikipedia says, the C++ std::deque class uses a
multiple-array implementation which is probably not so good for realtime purposes:
https://en.wikipedia.org/wiki/Double-ended_queue */

template<class T>
class rsDoubleEndedQueue : public rsBuffer<T>
{

public:

  /** Constructor. Allocates enough memory for an internal buffer to hold the "capacity" number of
  queued values. The memory allocated for the buffer will be the smallest power of two that is
  greater or equal to the desired capacity plus 2 (two additional and unusable index values are
  required to make the index arithmetic work right). */
  rsDoubleEndedQueue(size_t capacity) : rsBuffer<T>(capacity+2) {}

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

  /** Returns the value at the tail of the queue, i.e. the leftmost value. */
  inline T readTail() const
  {
    RAPT::rsAssert(!isEmpty(), "Trying to read from empty deque");
    return this->data[this->wrap(tail+1)];
  }
  // or maybe readFirst/Last Front/Back or peekFront/peekBack

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
  inline size_t getLength() const 
  { 
    size_t L = wrap(head - tail - 1);   // for debug
    rsAssert(L <= getMaxLength());

    return this->wrap(head - tail - 1); 
  }

  /** Returns the maximum allowed length for the queue. Due to the way, the head- and tail pointers
  operate, the usable length is 2 less than the size of the underlying data buffer. */
  inline size_t getMaxLength() const { return this->data.size()-1; } //
  //inline size_t getMaxLength() const { return this->data.size()-2; } //
  // hmm..i think, we can actually use size()-1 memory locations

  /** Returns true if the queue is empty. */
  inline bool isEmpty() const { return getLength() == 0; }

  /** Returns true, if the queue is full. */
  //inline bool isFull() const { return getLength() >= getMaxLength(); }
  inline bool isFull() const { return getLength() > getMaxLength(); }
  // shouldn't that be >= ? ...somehow it seems, it can actually take one more entry ...but only
  // under certain conditions? more tests needed....
  // it seems >= is the right one, but if we use it, we trigger an assert in movingMaximumUnitTest
  // although the test actually passes


protected:

  size_t head = 1, tail = 0; // check out correct terminology in the literature
  // The head-pointer is always one position to the right of rightmost value and the tail pointer
  // is always one position to the left of the leftmost value (both with wraparound). This is
  // required in order to compute the length of the queue from the head/tail indices. Placing them
  // directly ON the leftmost/rightmost values would not allow to distiguish between an empty
  // and one-element queue (in both cases, we would have head == tail). So the actual queued values
  // are at indices i with: tail < i < head, not: tail <= i <= head and when head = tail + 1, the
  // queue is empty because no index i fits in between tail and head
};

/*
ToDo:
-allow to dynamically resize the capacity at runtime - but such resize operations should be an
 absolute exception in realtime code (only a last resort - they actually should be considered user
 errors)
-maybe provide push-functions that take a vector argument and pop-functions that pop a range of
 values at once
-maybe provide some sort of random access functions:
 T fromFront(size_t i) where i = 0 returns the front element, 1 the one after the front, etc.
 or fromHead - analogously for fromBack/fromTail/fromRight
-maybe binary index search functions (at which index, counted from front (or back) does the value
 in the que get greater than some value
 size_t searchFromFront...this may become relevant for optimizing the moving maximum filter
*/

// maybe make convenience subclasses rsQueue and rsStack - both datastructures have different
// subsets of functionality of the double ended queue actually, but the convenience functions could
// have more stacky or queuey names like push/pop/top, enqueue/dequeue/front/back - they should
// use protected inheritance and use delegation, like Stack::push delegates to Base::pushFront, 
// etc.




