#ifndef RAPT_TREES_H
#define RAPT_TREES_H

/** This file contains tree based data structures. */


/** Baseclass for amost complete binary tree data-structures such as rsBinaryHeap (and later maybe
rsBinarySearchTree, which is still under construction in the prototypes section).  */

template<class T>
class rsBinaryTree
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
  {
    setData(&newData[0], (int)newData.size(), (int)newData.size());
  }

  /** Sets the comparison function to be used. If it implements "less-than", you'll get a max-heap
  and if it implements "greater-than", you'll get a min-heap. By default, it's assigned to a 
  less-than function based on the < operator. */
  void setCompareFunction(const std::function<bool(const T&, const T&)>& newFunc) 
  { less = newFunc; }
  // question: what happens, if we use a less-or-equal or greater-or-equal function for this?

  /** Sets the function used to swap two elements. By default, it just uses rsSwap but you may want
  to do something extra whenever a swap takes place in certain circumstances, for example, to keep 
  track of when items get moved around (see rsMovingQuantileFilter for an example). */
  void setSwapFunction(const std::function<void(T&, T&)>& newFunc)
  { swap = newFunc; }


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
  // maybe these should be references?

  // Some related classes need direct acces to our members:
  //template<class U> friend class rsDoubleHeap;

};



#endif