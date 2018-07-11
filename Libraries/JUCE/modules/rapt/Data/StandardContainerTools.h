#ifndef RAPT_STANDARDCONTAINER_H_INCLUDED
#define RAPT_STANDARDCONTAINER_H_INCLUDED

/** A collection of (convenience) functions for the container classes of the C++ standard library 
such as std::vector, etc. */

//=================================================================================================
// functions for std::vector

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
inline void rsRemove(std::vector<T>& v, size_t index)
{
  v.erase(v.begin() + index);
}

template<class T>
inline size_t rsFind(std::vector<T>& v, T elementToFind)
{
  for(size_t i = 0; i < size(v); i++)
    if(v[i] == elementToFind)
      return i;
  return size(v); // as convention, return vector-length, if element is not found
}

template<class T>
inline T rsLast(const std::vector<T>& v)
{
  return v[v.size()-1];
}

template<class T>
inline bool rsRemoveFirstOccurrence(std::vector<T>& v, T elementToRemove)
{
  for(size_t i = 0; i < size(v); i++)
    if(v[i] == elementToRemove){
      remove(v, i);
      return true; }
  return false;
}

template<class T>
inline bool rsContains(std::vector<T>& v, T elementToCheckFor)
{
  for(size_t i = 0; i < size(v); i++)
    if(v[i] == elementToCheckFor)
      return true;
  return false;
}

template<class T>
inline void rsAppendIfNotAlreadyThere(std::vector<T>& v, T newElement)
{
  if(!contains(v, newElement))
    append(v, newElement);
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



#endif