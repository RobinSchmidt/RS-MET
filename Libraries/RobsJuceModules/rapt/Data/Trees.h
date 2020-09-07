#ifndef RAPT_TREES_H
#define RAPT_TREES_H

/** This file contains tree based data structures. */


/** Baseclass for almost complete binary tree data-structures that are stored in an external array,
i.e. this class acts like "view" or "wrapper" around an already existing data array. An example is
the subclass rsBinaryHeap.   */

template<class T>
class rsBinaryTree  // maybe rename to rsImplicitBinaryTree
{

public:

  rsBinaryTree(T* newData = nullptr, int newSize = 0, int newCapacity = 0)
  {
    setData(newData, newSize, newCapacity);
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the data that should be treated as tree. The object does not take ownership of the data.
  Instead, it acts pretty much like a "view" (as in rsMatrixView) - it just operates on an existing
  data array whose lifetime is managed elsewhere. */
  void setData(T* newData, int newSize, int newCapacity)
  {
    data = newData;
    size = newSize;
    capacity = newCapacity;
  }

  void setData(std::vector<T>& newData)
  { setData(&newData[0], (int)newData.size(), (int)newData.size()); }

  /** Sets the comparison function to be used. If it implements "less-than", you'll get a max-heap
  and if it implements "greater-than", you'll get a min-heap. By default, it's assigned to a
  less-than function based on the < operator. */
  void setCompareFunction(const std::function<bool(const T&, const T&)>& newFunc)
  { less = newFunc; }
  // maybe that should be in the subclass rsBinaryHeap - there may be trees with other properties
  // that do not rely on a less-than comparison
  // question: what happens, if we use a less-or-equal or greater-or-equal function for this?

  /** Sets the function used to swap two elements. By default, it just uses rsSwap but you may want
  to do something extra whenever a swap takes place in certain circumstances, for example, to keep
  track of when items get moved around (see rsMovingQuantileFilter for an example). */
  void setSwapFunction(const std::function<void(T&, T&)>& newFunc)
  { this->swap = newFunc; }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the number of nodes in the tree. */
  int getSize() const { return size; }
  // maybe make a version that takes a node-index as parameter

  /** Returns the capacity, i.e. the maxumum number of nodes that can be stored in the underlying
  data array. */
  int getCapacity() const { return capacity; }

  /** Read/write access to elements. Warning: overwriting elements may destroy the defining
  property (heap, binary-search, etc.) of the tree, so it should be used only if you know what you
  are doing. */
  T& operator[](int i)
  {
    rsAssert(i >= 0 && i < size, "Index out of range");
    return data[i];
  }

  T& operator[](size_t i)
  {
    rsAssert((int) i < size, "Index out of range");
    return data[i];
  }

  const T& at(int i) const
  {
    rsAssert((int) i < size, "Index out of range");
    return data[i];
  }
  // get rid of the duplication

  /** Returns true, if a is less than b as defined by our less function as et by
  setCompareFunction. */
  bool isLess(const T& a, const T& b) const { return less(a, b); }
  // maybe rename to comesBefore to reflect that the comparision function not necessarily a
  // less-than function


protected:

  /** Index of parent of node i. */
  static inline int parent(int i) { return (i-1) / 2; }

  /** Index of left child of node i. */
  static inline int left(int i)   { return 2*i + 1; }

  /** Index of right child of node i. */
  static inline int right(int i)  { return 2*i + 2; }

  /** Returns true, iff index i is the index of a left child node. */
  static inline bool isLeft(int i) { return i == left(parent(i)); }
  // needs test - i think, we can do it simpler: the odd indices are left children and the even
  // indices are right children...verify that

  /** Returns true, iff index i is the index of a right child node. */
  static inline bool isRight(int i) { return i == right(parent(i)); }
  // needs test

  // The actual data array:
  T* data = nullptr;
  int size = 0;
  int capacity = 0;

  // Comparison and swapping functions:
  std::function<bool(const T&, const T&)>
    less = [](const T& a, const T& b)->bool { return a < b;  };
  std::function<void(T&, T&)>
    swap = [](T& a, T& b) { rsSwap(a, b); };
  // maybe these should be references? hmm..that would make the class inconveninient to use

};

//=================================================================================================

/** Class for representing an array A of data as binary heap with functions for establishing and
maintaining the heap-property. To understand what that property means, we must first interpret the
flat array A as a binary tree in the following way:

  left(i)   = 2*i + 1          left child index of index i
  right(i)  = 2*i + 2          right child index of index i
  parent(i) = (i-1) / 2        parent index of index i

Given that, the heap property says that for every node with index i (except the root), it holds
that:

  A[parent(i)] >= A[i]

Instead of >=, we could have also used <=. In the former case, we are dealing with a max-heap in
the latter with a min-heap.

References:
  (1) Introduction to Algorithms, 2nd Ed. (Cormen, Leiserson, Rivest, Stein)  */

template<class T>
class rsBinaryHeap : public rsBinaryTree<T>
{

public:

  using rsBinaryTree<T>::rsBinaryTree;  // inherit constructors

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns true, iff the underyling data array currently satisfies the heap property. The
  function is meant mostly for testing and debugging purposes and is currently implemented
  recursively (i.e. inefficiently). */
  bool isHeap(int i = 0) const;




  //-----------------------------------------------------------------------------------------------
  /** \name Data Manipulation */

  /** Reorders the underlying data array so as to satisfy the heap property. Runs in O(N) time.
  From the code, it would appear as having O(N*log(N)) complexity because we call an O(log(N))
  function inside the loop. However, the runtime of floatDown depends on the argument i in such
  a way to give an overall O(N) behavior (see reference (1)). The memory complexity is O(1). */
  void buildHeap()
  {
    for(int i = this->size/2-1; i >= 0; i--)  // or should we use (size-1)/2 ?
      floatDown(i);
  }

  /** Replaces the element at index i with the new given element x and rebalances the heap to
  maintain the heap-property which amounts to floating the new element x up or down. The return
  value is the array index, where the new element actually ended up. */
  int replace(int i, const T& x) { this->data[i] = x; return floatIntoPlace(i); }

  /** Inserts the element x into the heap at the correct position to maintain the heap-property and
  returns the array-index where it was inserted. */
  int insert(const T& x);

  /** Removes the element at given index i from the heap and re-orders the remaining elements to
  maintain the heap-property. */
  void remove(int i);

  /** Removes the first element from the heap */
  void removeFirst() { remove(0); }
  // maybe we should have a special implementation to remove the front element? i think, in this
  // case, we may simply call floatDown instead of floatIntoPlace

  /** Removes the first element and returns it */
  T extractFirst() { T first = this->data[0]; removeFirst(); return first; }


protected:

  /** Functions to establish or maintain the heap-property of the underlying data array. */

  /** Lets the node at index i float up or down into its appropriate place and returns the
  resulting new array index. */
  int floatIntoPlace(int i) { return floatUp(floatDown(i)); }

  /** Lets the node at index i float up the tree if its parent violates the heap property. Returns
  the new array index. */
  int floatUp(int i);

  /** Assuming that the subtrees rooted at left(i) and right(i) satisfy the heap property, this
  function makes sure that node/index i also satifies the heap property. If it doesn't already
  satisfy it, the function lets the value at i float down the subtree rooted at i. It has a time
  complexity of O(log(N)) and memory complexity of O(1). The return value is the new array index of
  the value that was previously located at index i. */
  int floatDown(int i);


  // todo: implement functions:  T peek(int i), extractMax, increaseKey, heapMax, see (1)
  // todo: implement int find(const T& x); ...should be possible in log(N) time complexity, or 
  // maybe not? if both children are larger than our searched value, we would have to traverse
  // two subtrees..but maybe at least O(log(N)) on average?

  template<class U> friend class rsDoubleHeap;

};

//=================================================================================================

/** Data structure that combines a max-heap with a min-heap for the purpose of splitting a bunch
of values into two groups: the small values and the large values. The two heaps satisfy the
property that all values in the small heap are less-or-equal to all values in the large heap. The
small values are held in a max-heap such that the largest of them can be extracted in O(1) time and
the large values are held in a min-heap such that the smallest of them can also be extracted in
O(1) time. This facilitates to extract the median or other quantiles.

Heap items can referred to by an integer that serves as key. We use the convention that when the
sign bit is zero, the number is treated directly as index into the small heap and if it's one,
it's an index into the large heap given by whatever number results when the sign bit is masked
away. But client code should really not care about this implementation detail and just treat the
keys as opaque identifiers. But note that any manipulations such as insertion, removal or
replacement of heap elements will invalidate the keys. If you want to refer to an item later, you
need to keep track of how its key changes due to any heap operations that may have occured in
the meantime. You may do this by implementing the swap function in a way that allow you to track
where an item floats to in the double-heap (and therefore, how its key changes). Basically, what
you need to do when swapping items is to just swap your stored keys, too. For this, an item needs
a way to know, where its key is stored. See rsQuantileFilterCore for an example how this may work.
It's a bit tricky (and maybe you don't need that feature anyway). */

#undef small // some silly header #defines this as macro - who would do such a thing? pure evil!

template<class T>
class rsDoubleHeap
{

public:

  rsDoubleHeap()
  {
    large.less = [](const T& a, const T& b)->bool { return b < a;  };
    // large.less is actually a greater-than function which we obtain from the less-than operator
    // by swapping the arguments - this turns large into a min-heap rather than a max-heap. todo:
    // use the more generic name comp for "compare" in rsBinaryTree for the comparison function
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the data arrays that actually hold the underlying data. The lifetime of these is the
  responsibility of client code. */
  void setData(T* newSmallData, int newSmallSize, int newSmallCapacity,
    T* newLargeData, int newLargeSize, int newLargeCapacity)
  {
    small.setData(newSmallData, newSmallSize, newSmallCapacity);
    large.setData(newLargeData, newLargeSize, newLargeCapacity);
  }

  /** Convenience function. */
  void setData(std::vector<T>& newSmall, std::vector<T>& newLarge)
  {
    setData(&newSmall[0], (int) newSmall.size(), (int) newSmall.size(),
      &newLarge[0], (int) newLarge.size(), (int) newLarge.size());
  }

  /** Sets the function to be used to swap two heap elements. */
  void setSwapFunction(const std::function<void(T&, T&)>& newFunc)
  {
    small.setSwapFunction(newFunc);
    large.setSwapFunction(newFunc);
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  bool isKeyValid(int k) const
  {
    if(isKeyInLargeHeap(k))
    {
      // for debuging a weird bug in gcc:
      int lhi = toLargeHeapIndex(k);
      int lhs = large.getSize();
      int shs = small.getSize();  // just for info
      rsAssert(lhi < lhs);
      // we need a test for rsDoubleHeap that fills up both heaps up to their capacity. somehow,
      // with gcc, it tries to add an element to the large heap when it's already full. this
      // happens also when the size of the heap is well below the allocated capacity. It happens
      // only with the random int and linearly descending vector, not with the linearly ascending
      // one
      // could it be that rsQuantileFilterCore::wrap behaves differently in different compilers, i.e.
      // the modulo operator behaves differently for negative inputs?

      return toLargeHeapIndex(k) < large.getSize();
    }
    else
      return k < small.getSize();
  }

  bool isIndexValid(int i) { return i >= 0 && i < small.getSize() + large.getSize(); }


  int getNumSmallValues() const { return (int) small.getSize(); }
  int getNumLargeValues() const { return (int) large.getSize(); }
  int getNumValues()      const { return (int) (small.getSize() + large.getSize()); }
  bool isSmallSizeOne()   const { return small.getSize() == 1; }
  bool isLargeSizeOne()   const { return large.getSize() == 1; }


  /** Returns the largest (i.e. first) element from our max-heap of small values. */
  T getLargestSmallValue()  const { return small.at(0); }

  /** Returns the smallest (i.e. first) element from our min-heap of large values  */
  T getSmallestLargeValue() const { return large.at(0); }

  /** Returns the second largest element from our max-heap of small values. */
  T get2ndLargestSmallValue(/*T retValWhenHeapLengthIs1*/) const
  {
    if(this->small.getSize() >= 3) return rsMaxViaLess(this->small.at(1), this->small.at(2));
    if(this->small.getSize() >= 2) return this->small.at(1);
    return small.at(0);  // does this make sense? ...nope!
    //return retValWhenHeapLengthIs1;
  }

  /** Returns the second smallest element from our min-heap of large values. */
  T get2ndSmallestLargeValue(/*T retValWhenHeapLengthIs1*/) const
  {
    if(this->large.getSize() >= 3) return rsMin(this->large.at(1), this->large.at(2));
    if(this->large.getSize() >= 2) return this->large.at(1);
    return large.at(0);
    //return retValWhenHeapLengthIs1;
  }
  // ..returning large.at(0) when the length is 1 does not seem to be the right thing to do in
  // rsQuantileFilterCore::readOutputWithOneMoreInput - instead, maybe the caller should specify
  // what should be returned in this case
  // maybe return const references instead of variables by value

  /** Returns the last element from our max-heap of small values (this is not  necessarily the
  smalles value, because in a heap, the order of the children of a node is unspecified). */
  T getLastSmallValue()     const { return small.at(small.getSize()-1); }


  T getLastLargeValue()     const { return large.at(large.getSize()-1); }
  // rename to getFirstSmallValue getLast.. the last is not necessarily the largest or smallest
  // the children are both larger/smaller but their order is not defined...hmm - can we define
  // it to give a heap even more structure? could that be usefu?


  /** Returns a preliminary key which is a key from the end of either the small or large heap.
  This is needed in rsQuantileFilterCore for appending nodes. */
  int getPreliminaryKey(const T& newValue)
  {
    if(small.less(newValue, large[0]))
      return small.getSize();
    else
      return large.getSize() | firstBitOnly;
  }
  // maybe rename to something more meaningful - but what?

  /** Returns true, iff this object satisfies the double-heap property. Meant mostly for testing
  and debugging (is costly to call!). */
  bool isDoubleHeap()
  {
    return small.isHeap() && large.isHeap()
      && !(small.less(large[0], small[0]));  // last condition means small[0] <= large[0]
  }


  static bool isKeyInLargeHeap(int k) { return k & firstBitOnly;  }

  /** Returns the raw large-heap index, which is the key k with the first bit shaved off. */
  static inline int rawLargeHeapIndex(int k) // rename to largeHeapKeyToIndex
  { return k & allBitsButFirst; }

  int toLargeHeapIndex(int k)  const { return rawLargeHeapIndex(k); }
  // get rid

  //-----------------------------------------------------------------------------------------------
  /** \name Data access */

  /** Replaces the element at given key with the given new value and returns the the key where the
  new element actually ended up after doing all the floating up/down and potential swapping
  business. */
  int replace(int key, const T& newValue);

  /** Inserts a new item and returns the key that indicates, where it ended up. */
  int insert(const T& newValue)
  {
    if(small.less(newValue, large[0]))
      return small.insert(newValue);
    else
      return large.insert(newValue) | firstBitOnly;
  }

  /** Removes the item at the given key. */
  void remove(int key)
  {
    rsAssert(isKeyValid(key), "Key out of range");
    if(isKeyInLargeHeap(key))
      large.remove(toLargeHeapIndex(key));
    else
      small.remove(key);
  }

  /** Element access via an integer key.  */
  T& atKey(int k)
  {
    rsAssert(isKeyValid(k), "Key out of range");
    if(isKeyInLargeHeap(k))
      return large[toLargeHeapIndex(k)];
    else
      return small[k];
  }

  /** Element access via an integer index. If nS is the number of values in the small heap, indices
  i < nS refer directly to samples in the small heap with that same small-heap-index i, whereas
  indices >= nS are interpreted as large-heap-index i-nS into the large heap. */
  T& atIndex(int i)
  {
    rsAssert(isIndexValid(i), "Index out of range");
    int nS = small.getSize();
    if(i < nS)
      return small[i];
    else
      return large[i-nS];
  }

  /** Transfers the smallest of the large values to the heap of small values and returns the key
  where it ended up. This will shrink the large heap and grow the large heap by one element. */
  int moveFirstLargeToSmall() { T e = large.extractFirst(); return small.insert(e); }
  // hmm - return value should be always 0 anyway

  int moveFirstSmallToLarge()
  { T e = small.extractFirst(); return large.insert(e) + small.getSize(); }
  // return value should always be small.getSize(), size measured after the operation


  /** Conversion from index to key. */
  int indexToKey(int i)
  {
    int nS = small.getSize();
    if(i < nS) return i;
    else       return (i-nS) | firstBitOnly;
  }
  // takes the double-heap index from 0....nS+nL-1
  // needs test

  /** Conversion from key to index. */
  int keyToIndex(int k)
  {
    if(isKeyInLargeHeap(k)) {
      int i = toLargeHeapIndex(k);
      return i + small.getSize(); }
    else
      return k;
  }
  // returns the double-heap index from 0....nS+nL-1
  // needs test


  // when rsQuantileFilterCore has been move to rapt, move this into protected and uncomment
  // friend:
  rsBinaryHeap<T> small, large;  // the two heaps for the small and large numbers
  //template<class U> friend class ::rsQuantileFilterCore; // temporary - try to get rid

protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Internal */






  //template<class U> friend class rsQuantileFilterCore2; // preliminary - try to get rid

};
// clean this up - order the functions according to the usual scheme and add documentation


#endif
