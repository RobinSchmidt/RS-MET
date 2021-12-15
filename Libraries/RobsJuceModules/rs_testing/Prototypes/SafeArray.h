// todo: rename file to rsArrays - we have different variants of arrays here..



//=================================================================================================

/** Under Construction

A class that can be used like std::vector but has the feature that it doesn't re-allocate memory 
when it needs to grow. Instead, the already allocated memory is retained and additional memory is 
allocated elsewhere. So, the data is not necessarily in a contiguous block of memory. The lack
of reallcoation is useful for arrays of objects to which a pointer or reference exists elsewhere
in the program. Re-allocation would require to invalidate all these pointers but with this class, 
they remain valid ...tbc.... 

ToDo:
-Maybe rename to rsRetainedArray. It's not true that it never ever re-allocates. Its just that it 
 doesn't implicitly reallocate in push_back when it needs to grow. Later we can add funtions that 
 explicitly may trigger reallocation such as "defragmentMemory".
-Maybe then it could make sense to also have a class rsContiguous array which would basically be a
 re-implementation of std::vector. Maybe not for production use, more as "excercise" and as 
 reference for implementing more complex STL compliant containers later.
*/

template<class T>
class rsNonReAllocatingArray
{

public:


  rsNonReAllocatingArray() {}



  size_t size() const { return totalSize; }

  size_t capacity() const { return totalCapacity; }


  bool empty() const { return size() == 0; }

  void reserve(size_t numElems);

  void push_back(const T& elem);


  void clear() { totalSize = 0; totalCapacity = 0; chunks.clear(); }



  /** Random access to element at given index i. The complexity is linear in the number of chunks.
  Try to avoid using this when iterating over all the elements in the array. If possible, use an 
  iterator which gives constant complexity for sequential access. ...tbc... */
  T& operator[](const size_t i);

  // todo: have an init(size_t initialCapacity) function


  /** Returns reference to element k in chunk j. */
  T& at(size_t j, size_t k) { return chunks[j][k]; }
  // this belongs into the private section! it's only public during development for convenience!


  // inc/dec should be done in constant time, +=, -= etc. may be linear in thenumber of chunks
  // ...but maybe it can be made constant time, too? let's at least try it!
  class iterator
  {
  public:  // try to get rid and make as much as possible private

    iterator(rsNonReAllocatingArray& iteratee, size_t chunkIndex, size_t elemIndex) 
      : a(iteratee), j(chunkIndex), k(elemIndex) 
    {
      if(a.capacity() == 0) {
        rsError("not yet implemented"); }
      else {
        c0 = a.chunks[0].size(); }
    }

    iterator& operator++() 
    { 
      rsError("not yet implemented"); 
      return *this;
    } // pre-inc

    iterator& operator--()    
    {
      rsError("not yet implemented"); 
      return *this;
    } // pre-dec

    iterator  operator++(int) { iterator tmp = *this; ++*this; return tmp; } // post-inc
    iterator  operator--(int) { iterator tmp = *this; --*this; return tmp; } // post-dec


    // Dereferencing:
    T& operator*()   { return   a.at(j, k);  }
    T* operator->()  { return &(a.at(j, k)); }


  private:

    rsNonReAllocatingArray& a; // the array to iterate over
    size_t j, k;               // chunk- and element index within chunk
    size_t c0;                 // initial capacity

  };

  iterator begin() { return iterator(*this, 0, 0); }

  iterator end() 
  { 
    size_t j, k;
    flatToChunkAndElemIndex(totalSize);
    return iterator(*this, j, k); 
  }



protected:



  /** Converts a flat index as seen by the user to the chunk-index and element index within the
  chunk */
  void flatToChunkAndElemIndex(size_t flatIndex, size_t& chunkIndex, size_t& elemIndex) const;

  std::vector<std::vector<T>> chunks;
  /**< chunks[0] contains the 1st chunk of the whole array, chunks[1] the 2nd and so on. 
  Initially, there's only 1 chunk. */

  size_t totalSize = 0;
  /**< Currently used total size. */

  size_t totalCapacity = 0;
  /**< Sum of the sizes of the chunks, cached for quick access. The internal vectors are always
  initially resized (and not just reserved) to their desired sizes/capacities, i.e. for the chunks,
  there's no distinction between their internal size and capacity .*/

};

template<class T>
void rsNonReAllocatingArray<T>::flatToChunkAndElemIndex(size_t i, size_t& j, size_t& k) const
{
  j = 0;   // chunk index
  k = i;   // element index within chunk
  size_t s = chunks[0].size();            // s is the initial capacity
  if(   k >= s) { k -= s; ++j; }          // first 2 chunks have size s
  while(k >= s) { k -= s; ++j; s *= 2; }  // after that, they grow by factor 2

  // If the initial capacity s is 4, the chunk sizes are: 4,4,8,16,32,64,... The next chunk size
  // is always equal to the sum of the sizes that came before it. That's why the first needs to
  // violate the 2^k pattern. It has to be *some* power of two that the user can pick via the 
  // initial call to reserve.
}

template<class T>
T& rsNonReAllocatingArray<T>::operator[](const size_t i) 
{ 
  RAPT::rsAssert(i < totalSize);
  size_t j, k;
  flatToChunkAndElemIndex(i, j, k);
  return chunks[j][k];
}


template<class T>
void rsNonReAllocatingArray<T>::reserve(size_t numElems)
{
  rsAssert(rsIsPowerOfTwo(numElems), "Must be a power of two!");
  // todo: relax that - just round up to the next power of two


  if(numElems <= totalCapacity)
    return;  // nothing to do
  if(chunks.size() == 0)
  {
    // This branch is for the initial reservation shortly after creation or before usage.
    chunks.resize(1);
    chunks[0].resize(numElems);
    totalCapacity = numElems; 
  }
  else
  {
    // This branch is for later reservation of more memory than was initially reserved.
    size_t needed = totalCapacity - numElems;

    // ToDo:
    // -figure out, how many more chunks are needed - we must follow the pattern: 
    //  C,C,2C,4C,8C,16C etc. where C is the initial capacity, i.e. chunks[0].size()

    rsError("not yet implemented");
  }
}

template<class T>
void rsNonReAllocatingArray<T>::push_back(const T& elem)
{
  if(totalCapacity == 0)
    reserve(1);
  size_t j, k;
  flatToChunkAndElemIndex(totalSize, j, k);
  if(j >= chunks.size() )
  {
    // A new chunk needs to be added. The size of the new chunk is equal to the current total size.
    // This ensures that the capacity is always a power of two given that the initial capacity was
    // a power of two (which we ensure in reserve).
    rsAssert(j == chunks.size() && k == 0);  // if j >= chunks.size() happens, then like that
    chunks.push_back(std::vector<T>(totalSize));
  }
  chunks[j][k] = elem;
  ++totalSize;
}


//=================================================================================================

/** just a STUB at the moment, experimental

An attempt to wrap a raw C-array to give it (partially) the interface of std::vector. The goal
is to make std::algorithms like std::sort, std::partial_sort, std::transform, std::accumulate, etc.
available to raw arrays. If you have a C-array in your code like so:

  int N = 10;
  float *a = malloc(N*sizeof(float));
  // ...

you could wrap an rsArrayView around it and apply algorithms from the standard library to it:

  rsArrayView<float> v(a, N);
  std::sort(v.begin(), v.end(), [](float a, float b){ return a < b; });

...that's the goal, at least. The rationale is that much of the time, we have data as raw arrays or
wrapped in some data structure other than std::vector or similar and we want to harness the 
std:algorithms in these cases, too. 

Maybe later, I can add my own algorithms, like:

  rsConvolve(const rsArrayView<T>& x, const rsArrayView<T>& h, rsArrayView<T>& y);

*/

template<class T>
class rsArrayView
{

public:


  using difference_type = std::ptrdiff_t;




  rsArrayView(T* data, size_t numElements) : _data(data), _numElems(numElements) {}


  // functions for compatibility with std::vector


  T* data() { return _data; }

  size_t size() const { return _numElems; }





  // todo: iterators, dereferencing, etc

  class iterator
  {
  public:  // try to get rid and make as much as possible private
    iterator(T* data, size_t index) : _data(data), _index(index) {}

    iterator& operator++()    { _index++; return *this; }                    // pre-inc
    iterator  operator++(int) { iterator tmp = *this; ++*this; return tmp; } // post-inc
    iterator& operator--()    { _index--; return *this; }                    // pre-dec
    iterator  operator--(int) { iterator tmp = *this; --*this; return tmp; } // post-dec


    //iterator operator+(const iterator& r) const { return iterator(_data, _index + r._index); }
    //iterator operator-(const iterator& r) const { return iterator(_data, _index - r._index); }
    // wait..no...that makes no sense: r should be a size_t, not another irterator
    // std::sort requires a binary - operator that takes two iterators as inputs and returns
    // some sort of "iterator-difference" type
    // https://stackoverflow.com/questions/46695349/how-to-handle-iteratordifference-type-when-you-have-no-way-of-measuring-the-di
    // https://en.cppreference.com/w/cpp/iterator/iterator_traits
    // https://www.py4u.net/discuss/66628
    // https://www.cplusplus.com/reference/iterator/iterator_traits/
    // https://www.cplusplus.com/reference/iterator/iterator/

    // https://www.fluentcpp.com/2018/05/08/std-iterator-deprecated/

    // hmmm - none of this make std::sort work with the iterators:

    //typedef std::ptrdiff_t difference_type;
    //difference_type operator-(const iterator& r) const { return _index - r._index; }

    //using difference_type = size_t;
    //typedef size_t difference_type;
    //difference_type operator-(const iterator& r) const { return _index - r._index; }


    using difference_type = std::ptrdiff_t;
    difference_type operator-(const iterator& r) const { return _data - r._data; }
    //difference_type operator-(const iterator& r) const { return _index - r._index; }


    //std::ptrdiff_t operator-(const iterator& r) const { return _index - r._index; }

    //size_t operator-(const iterator& r) const { return _index - r._index; }

    //_Iter_diff_t<iterator> operator-(const iterator& r) const { return _index - r._index; }

    //template<class iterator>
    //std::iterator_traits<iterator>::difference_type 
    //  operator-(const iterator& r) const { return _index - r._index; }

    //template<class iter>
    //std::iterator_traits<iter>::difference_type 
    //  operator-(const iterator& r) const { return _index - r._index; }


    // todo: 
    // -binary +,- and +=, -=
    // -dereferencing -> (...i think)


    T* _data = nullptr;
    size_t _index = 0;
  };


  iterator begin() { return iterator(_data, 0);         }
  iterator end()   { return iterator(_data, _numElems); }  // or -1?

  // resize can perhaps only be done for shrinking the size...unless the caller passes on 
  // construction a pointer that points to more allocated memory than N indicates....hmm...


  T& operator[](const int i) { return _data[i]; }
  const T& operator[](const int i) const { return _data[i]; }
  // todo: overload for size_t and, importantly, iterator


  T& operator[](const iterator i) { return _data[i._index]; }



protected:

  T* _data = nullptr;
  size_t _numElems = 0;

};

// see:
// https://en.cppreference.com/w/cpp/container/vector
// https://en.cppreference.com/w/cpp/algorithm/sort

// https://stackoverflow.com/questions/3182843/writing-stl-compatible-iterators
// https://www.fluentcpp.com/2018/04/24/following-conventions-stl/

// https://docs.microsoft.com/en-us/cpp/cpp/increment-and-decrement-operator-overloading-cpp?view=msvc-170


//=================================================================================================

template<class T>
class rsSafeArray1
{

public:


  // todo: push_back, pop_back, insert, remove, etc...

  /** Read and write access to i-th element */
  //T& operator[](int i);



  /** Read-only access to i-th element */
  const T& operator[](int i) const
  {
    rsAssert(isValidIndex(i), "Invalid array index");
    return data[i];
  }


  bool isValidIndex(int i) const { return i >= 0 && i < numElems; }


private:

  std::vector<T> vec1, vec2;   // we always use either vec1 or vec2
  T* data = nullptr;           // points to &vec1[0] or &vec2[0] or is nullptr
  size_t numElems = 0;         // equals vec1.size() or vec2.size()
  std::recursive_mutex mutex;  // mutex for modifying the array

};

/*
template<class T>
T& rsSafeArray1<T>::operator[](int i)
{
  std::lock_guard<std::recursive_mutex> lock(mutex);
}
*/