// todo: rename file to rsArrays - we have different variants of arrays here..


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


    // hmmm - none of this compiles:
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