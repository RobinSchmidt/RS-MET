#ifndef RAPT_STANDARDCONTAINER_H_INCLUDED
#define RAPT_STANDARDCONTAINER_H_INCLUDED

/** A collection of convenience functions for the container classes of the C++ standard template 
library (STL), such as std::vector, std::map, etc. */

//=================================================================================================
// functions for std::vector

template<class T>
inline std::vector<T> rsConstantVector(size_t size, T value)
{
  std::vector<T> v(size);
  for(size_t i = 0; i < size; i++)
    v[i] = value;
  return v;
}

template<class T>
inline size_t rsSize(const std::vector<T>& v)
{
  return v.size();
}

template<class T>
inline void rsAppend(std::vector<T>& v, T newElement)
{
  v.push_back(newElement);
}

/** Appends vector w to vector v. */
template<class T>
inline void rsAppend(std::vector<T>& v, const std::vector<T>& w)
{
  size_t nv = v.size();  // old size of v
  size_t nw = w.size();  // 
  v.resize(v.size() + w.size()); // if v and w are the same, this will also change the size of w, 
  rsArray::copyBuffer(&w[0], &v[nv], (int)nw);  // ...that's why we needed to figure out nw before

  // another implementation - looks safer but the above works, too
  //if(&v[0] == &w[0]) {       // appending a vector to itself, we need a temporary...
  //  std::vector<T> tmp(w);
  //  rsAppend(v, tmp);
  //}
  //else {
  //  size_t nv = v.size();    // old size of v
  //  v.resize(v.size() + w.size());
  //  rsArray::copyBuffer(&w[0], &v[nv], (int)w.size());
  //}
}


template<class T>
inline std::vector<T> rsConcatenate(std::vector<T>& v, const std::vector<T>& w)
{
  std::vector<T> r(v.size() + w.size());
  for(size_t i = 0; i < v.size(); i++) r[i]            = v[i];
  for(size_t i = 0; i < w.size(); i++) r[i + v.size()] = w[i];
  return r;
}

template<class T>
inline bool rsAreVectorsEqual(const std::vector<T>& v, const std::vector<T>& w, double tolerance)
{
  if(v.size() != w.size())
    return false;
  for(size_t i = 0; i < v.size(); i++) {
    if(abs(v[i]-w[i]) > tolerance)
      return false;
  }
  return true;
}

template<class T>
inline void rsInsert(std::vector<T>& v, const T& newElement, size_t index)
{
  v.insert(v.begin() + index, newElement);
}

template<class T>
inline void rsInsertValue(std::vector<T>& v, T newElement, size_t index)
{
  v.insert(v.begin() + index, newElement);
}

template<class T>
inline void rsPrepend(std::vector<T>& v, const T& newElement)
{
  v.insert(v.begin(), newElement);
}

template<class T>
inline void rsRemove(std::vector<T>& v, size_t index)
{
  v.erase(v.begin() + index);
}

template<class T>
inline size_t rsFind(std::vector<T>& v, T elementToFind)
{
  for(size_t i = 0; i < v.size(); i++)
    if(v[i] == elementToFind)
      return i;
  return v.size(); // as convention, return vector-length, if element is not found
}

template<class T>
inline T rsLast(const std::vector<T>& v)
{
  return v[v.size()-1];
} // make also an rsFirst

template<class T>
inline bool rsRemoveFirstOccurrence(std::vector<T>& v, T elementToRemove)
{
  for(size_t i = 0; i < v.size(); i++)
    if(v[i] == elementToRemove){
      rsRemove(v, i);
      return true; }
  return false;
}

template<class T>
inline void rsReverse(std::vector<T>& v)
{
  rsArray::reverse(&v[0], (int) v.size());
}

template<class T>
inline bool rsContains(std::vector<T>& v, T elementToCheckFor)
{
  for(size_t i = 0; i < v.size(); i++)
    if(v[i] == elementToCheckFor)
      return true;
  return false;
}

template<class T>
inline void rsAppendIfNotAlreadyThere(std::vector<T>& v, T newElement)
{
  if(!rsContains(v, newElement))
    rsAppend(v, newElement);
}

template<class T>
inline T rsGetAndRemoveLast(std::vector<T>& v)
{
  T result = v[v.size()-1];
  v.pop_back();
  return result;
}

/** Converts C-array to std::vector. */
template<class T>
inline std::vector<T> toVector(T* theArray, size_t size) // rename to rsToVector
{
  std::vector<T> v(size);
  for(size_t i = 0; i < size; i++)
    v[i] = theArray[i];
  return v;
}

/** Multiplies a scalar and a vector. */
template<class T>
inline std::vector<T> operator*(const T& x, const std::vector<T>& v)
{
  std::vector<T> result(v.size());
  for(int i = 0; i < v.size(); i++)
    result[i] = x * v[i];
  return result;
}

/** Adds a scalar to a vector. */
template<class T>
inline std::vector<T> operator+(const std::vector<T>& v, const T& x)
{
  std::vector<T> result(v.size());
  for(int i = 0; i < v.size(); i++)
    result[i] = v[i] + x;
  return result;
}


/** Adds two vectors element wise. */
template<class T>
inline std::vector<T> operator+(const std::vector<T>& x, const std::vector<T>& y)
{
  size_t Nmax = std::max(x.size(), y.size());
  size_t Nmin = std::min(x.size(), y.size());
  std::vector<T> result(Nmax);
  for(int i = 0; i < Nmin; i++)
    result[i] = x[i] + y[i];
  return result;
}

/** Subtracts two vectors element wise. */
template<class T>
inline std::vector<T> operator-(const std::vector<T>& x, const std::vector<T>& y)
{
  size_t Nmax = std::max(x.size(), y.size());
  size_t Nmin = std::min(x.size(), y.size());
  std::vector<T> result(Nmax);
  for(int i = 0; i < Nmin; i++)
    result[i] = x[i] - y[i];
  return result;
}

//=================================================================================================
// functions for std::map

/** Checks, if a map contains a given key. */
template<class Key, class Value>
inline bool rsContains(const std::map<Key, Value>& map, const Key& key)
{
  auto iterator = map.find(key);
  if( iterator == map.end() )
    return false;
  return true;
}



#endif