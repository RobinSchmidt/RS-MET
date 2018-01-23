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

//=================================================================================================



#endif