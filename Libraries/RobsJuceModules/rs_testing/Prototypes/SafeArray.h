


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