#ifndef RAPT_STANDARDCONTAINER_H_INCLUDED
#define RAPT_STANDARDCONTAINER_H_INCLUDED

/** A collection of convenience functions for the container classes of the C++ standard template
library (STL), such as std::vector, std::map, etc. */

//=================================================================================================
// Conveniennce functions for std::vector
// maybe wrap into a class rsStdVectorTools

template<class T>
inline std::vector<T> rsAbs(const std::vector<T>& v)
{
  std::vector<T> va(v.size());
  for(size_t i = 0; i < v.size(); i++)
    va[i] = rsAbs(v[i]);
  return va;
}

/** Wraps iterator syntax to simplify calls to std::all_of. */
template<class T, class UnaryPredicate >
bool rsAllOf(const std::vector<T>& v, UnaryPredicate p)
{
  return std::all_of(v.cbegin(), v.cend(), p);
}

/** Wraps iterator syntax to simplify calls to std::any_of. */
template<class T, class UnaryPredicate >
bool rsAnyOf(const std::vector<T>& v, UnaryPredicate p)
{
  return std::any_of(v.cbegin(), v.cend(), p);
}

/** Returns a new vector which conatins elements where the function f is applied to each 
corresponding element in the input vector v. The function f should be a function pointer. 
Usage example:

  std::vector<double> y = rsApplyFunction(x, &log);

to get the logarithms of all values in a vector x of double precision numbers. */
template<class T>
std::vector<T> rsApplyFunction(const std::vector<T>& v, T (*f) (T))
{
  std::vector<T> r(v.size());
  for(size_t i = 0; i < r.size(); i++)
    r[i] = f(v[i]);
  return r;
}
// todo: can this be generalized to any kind of callable f by using a 2nd template parameter F for
// the function type?

template<class T>
inline std::vector<T> rsConstantVector(size_t size, T value)
{
  std::vector<T> v(size);
  for(size_t i = 0; i < size; i++)
    v[i] = value;
  return v;
}

/** Copies data from scr to dst. Both must have the same size. */
template<class T>
void rsCopy(const std::vector<T>& src, std::vector<T>& dst)
{
  rsAssert(src.size() == dst.size());
  rsArrayTools::copy(&src[0], &dst[0], (int)src.size());
}

/** For a vector v of pointers to objects, this functions deletes all the objects. */
template<class T>
inline void rsDeleteObjects(const std::vector<T*>& v)
{
  for(size_t i = 0; i < v.size(); i++)
    delete v[i];
}

/** Calls rsDeleteObjects and clears the vector. */
template<class T>
inline void rsDeleteObjectsAndClear(const std::vector<T*>& v)
{
  rsDeleteObjects(v);
  v.clear();
}

/** Like rsDeleteObjects but additionally sets the pointers to nullptr. */
template<class T>
inline void rsDeleteObjectsAndNull(const std::vector<T*>& v)
{
  for(size_t i = 0; i < v.size(); i++) {
    delete v[i];
    v[i] = nullptr; }
}


template<class T>
inline size_t rsSize(const std::vector<T>& v)
{
  return v.size();
}
// todo: return an int

template<class T>
inline void rsAppend(std::vector<T>& v, T newElement)
{
  v.push_back(newElement);
}

/** Appends vector w to vector v. */
template<class T>
inline void rsAppend(std::vector<T>& v, const std::vector<T>& w)
{
  if(w.size() == 0)
    return;
  size_t nv = v.size();  // old size of v
  size_t nw = w.size();  //
  v.resize(v.size() + w.size()); // if v and w are the same, this will also change the size of w,
  rsArrayTools::copy(&w[0], &v[nv], (int)nw);  // ...that's why we needed to figure out nw before

  // another implementation - looks safer but the above works, too
  //if(&v[0] == &w[0]) {       // appending a vector to itself, we need a temporary...
  //  std::vector<T> tmp(w);
  //  rsAppend(v, tmp);
  //}
  //else {
  //  size_t nv = v.size();    // old size of v
  //  v.resize(v.size() + w.size());
  //  rsArrayTools::copy(&w[0], &v[nv], (int)w.size());
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
// remove or derprecate - use rsEquals instead

template<class T>
inline bool rsEquals(const std::vector<T>& x, const std::vector<T>& y, T tol = T(0))
{
  if( x.size() != y.size() )
    return false;
  for(size_t i = 0; i < x.size(); i++) {
    if( rsAbs(x[i]-y[i]) > tol )
      return false; }
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
inline void rsInsert(std::vector<T>& v, const std::vector<T>& w, size_t index)
{
  v.insert(v.begin() + index, w.begin(), w.end());
}

/** Inserts value x into sorted vector v at its appropriate position. */
template<typename T>
typename std::vector<T>::iterator rsInsertSorted(std::vector<T>& v, T const& x)
{
  return v.insert(std::upper_bound(v.begin(), v.end(), x), x);
}
// https://stackoverflow.com/questions/15843525/how-do-you-insert-the-value-in-a-sorted-vector

template<class T>
inline void rsNegate(std::vector<T>& v)
{
  rsArrayTools::negate(&v[0], &v[0], (int) v.size());
}

/** Wraps iterator syntax to simplify calls to std::none_of. */
template<class T, class UnaryPredicate >
bool rsNoneOf(const std::vector<T>& v, UnaryPredicate p)
{
  return std::none_of(v.cbegin(), v.cend(), p);
}
// todo: make similar functions for any_of, all_of

template<class T>
inline void rsPrepend(std::vector<T>& v, const T& newElement)
{
  v.insert(v.begin(), newElement);
}

template<class T>
std::vector<T> rsRangeLinear(T min, T max, int N)
{
  std::vector<T> r(N);
  rsArrayTools::fillWithRangeLinear(&r[0], N, min, max);
  return r;
}
// maybe rename to rsLinSpace - consistent with numpy

template<class T>
inline void rsRemove(std::vector<T>& v, size_t index)
{
  v.erase(v.begin() + index);
}

/** Removes the range indexed from first to last, both ends inclusive! Attention: in
std::vector::erase, the "last" would be excluded - but that's really counter-intuitive and provokes
off-by-one bugs, that's why i use a both-ends-inclusive convention. */
template<class T>
inline void rsRemoveRange(std::vector<T>& v, size_t first, size_t last)
{
  v.erase(v.begin() + first, v.begin() + last + 1);
}

template<class T>
inline void rsPopFront(std::vector<T>& v)
{
  rsRemoveRange(v, 0, 0);
}



template<class TVal, class TIdx>
inline std::vector<TVal> rsSelect(const TVal* x, const std::vector<TIdx>& indices)
{
  std::vector<TVal> r(indices.size());
  for(TIdx i = 0; i < (TIdx) indices.size(); i++)
    r[i] = x[indices[i]];
  return r;
}

/** Returns a new vector that contains only those elements from the original vector v at the given
indices. TVal is the type of the values in v and TIdx is the type of the index-values - for example 
size_t or int. Note that it doesn't check, if the entries in the indicies array are actually valid
indices for v-array - making that sure is the responsibility of the caller. */
template<class TVal, class TIdx>
inline std::vector<TVal> rsSelect(const std::vector<TVal>& v, const std::vector<TIdx>& indices)
{
  return rsSelect(&v[0], indices);
  //std::vector<TVal> r(indices.size());
  //for(TIdx i = 0; i < (TIdx) indices.size(); i++)
  //  r[i] = v[indices[i]];
  //return r;
}

/*
template<class T>
inline std::vector<T> rsSelect(std::vector<T>& v, std::vector<size_t> indices)
{
  std::vector<T> r(indices.size());
  for(size_t i = 0; i < indices.size(); i++)
    r[i] = v[indices[i]];
  return r;
}
*/




template<class T>
inline void rsFill(std::vector<T>& v, T value)
{
  for(size_t i = 0; i < v.size(); i++)
    v[i] = value;
}

template<class T>
inline size_t rsFind(const std::vector<T>& v, T elementToFind)
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
  rsArrayTools::reverse(&v[0], (int) v.size());
}

template<class T>
inline void rsSetAllValues(std::vector<T>& v, T value)
{
  for(size_t i = 0; i < v.size(); i++)
    v[i] = value;
}

template<class T>
inline void rsSetZero(std::vector<T>& v)
{
  rsSetAllValues(v, T(0));
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

/** Resizes the vector to the given new size and if this new size is greater than the old size,
intializes all the new entries at the end with "value". */
template<class T>
inline void rsResizeWithInit(std::vector<double>& v, size_t newSize, T value)
{
  size_t oldSize = v.size();
  v.resize(newSize);
  for(size_t i = oldSize; i < newSize; i++)
    v[i] = value;
}

/** Prepends "amount" copies of "value" to the vector, i.e. increases the size of the vector by
"amount", shifts all values "amount" places to the right and inserts "amount" number of "value"s
at the front. */
template<class T>
inline void rsPadLeft(std::vector<double>& v, size_t amount, T value)
{
  if(amount == 0)
    return;
  size_t oldSize = v.size();
  v.resize(oldSize+amount);
  for(size_t i = v.size()-1; i >= amount; i--)
    v[i] = v[i-amount];
  for(size_t i = 0; i < amount; i++)
    v[i] = value;
  // maybe, it can be done with memmove and memset more efficiently?
}

template<class T>
std::vector<T> rsDifference(const std::vector<T> x)
{
  if(x.size() < 2)
    return std::vector<T>();  // result is empty
  std::vector<T> d(x.size()-1);
  for(size_t i = 0; i < d.size(); i++)
    d[i] = x[i+1] - x[i];
  return d;
}

/** Converts a vector of amplitudes to decibel values with a given floor decibel value to avoid
negative infinity for zero amplitudes. */
template<class T>
std::vector<T> rsAmpToDb(const std::vector<T>& a, T floorDb = -std::numeric_limits<T>::infinity())
{
  std::vector<T> db(a.size());
  for(size_t i = 0; i < a.size(); i++)
    db[i] = std::max(floorDb, 20*log10(a[i]));
  return db;
}

template<class T>
std::vector<T> rsDecimate(const std::vector<T>& x, int factor)
{
  int Ny = (int) x.size() / factor;
  std::vector<T> y(Ny);
  RAPT::rsArrayTools::decimate(&x[0], (int)x.size(), &y[0], factor);
  return y;
}

template<class T>
std::vector<T> rsDecimateViaMean(const std::vector<T>& x, int factor)
{
  int Ny = (int) x.size() / factor;
  std::vector<T> y(Ny);
  RAPT::rsArrayTools::decimateViaMean(&x[0], (int)x.size(), &y[0], factor);
  return y;
}

// these implementations cause compiler errors with gcc
///** Iterator to minimum element. */
//template<class T>
////auto rsMinIter(const std::vector<T>& x) { return std::min_element(x.cbegin(), x.cend()); } // gcc complains about auto
//std::vector<T>::const_iterator rsMinIter(const std::vector<T>& x)
//{ return std::min_element(x.cbegin(), x.cend()); }
//
///** Minimum element. */
//template<class T>
//T rsMinValue(const std::vector<T>& x) { return *rsMinIter(x); }
//
//template<class T>
//auto rsMaxIter(const std::vector<T>& x) { return std::max_element(x.cbegin(), x.cend()); }
//
//template<class T>
//T rsMaxValue(const std::vector<T>& x) { return *rsMaxIter(x); }


template<class T>
T rsMinValue(const std::vector<T>& x) { return rsArrayTools::minValue(&x[0], (int) x.size()); }

template<class T>
T rsMaxValue(const std::vector<T>& x) { return rsArrayTools::maxValue(&x[0], (int) x.size()); }

template<class T>
T rsMaxAbs(const std::vector<T>& x) { return rsArrayTools::maxAbs(&x[0], (int) x.size()); }

template<class T>
T rsMaxDeviation(const std::vector<T>& x, const std::vector<T>& y)
{
  rsAssert(x.size() == y.size());
  return rsArrayTools::maxDeviation(&x[0], &y[0], (int)x.size());
}

template<class T>
bool rsIsCloseTo(const std::vector<T>& x, const std::vector<T>& y, T tol)
{
  if(x.size() != y.size()) 
    return false;
  return rsArrayTools::almostEqual(&x[0], &y[0], (int) x.size(), tol);
}

//template<class T>
//T rsMaxValue(const std::vector<T>& x)
//{
//  return *std::max_element(x.cbegin(), x.cend());
//}

  //T max = std::numeric_limits<T>::min(); // we should instead use -inf for double/float? -> make explicit specilizations
  //for(size_t i = 0; i < x.size(); i++) {
  //  if(x[i] > max)
  //    max = x[i];
  //}
  //return max;

/** Returns the first index for which vector v has a nonzero value or -1 if all values are zero. */
template<class T>
int rsIndexOfFirstNonZero(const std::vector<T>& v)
{
  return rsArrayTools::firstIndexWithNonZeroValue(&v[0], (int)v.size()); 
  // todo: rename that to rsIndexOfFirstNonZero, too -> shorter and consistent
}

/** Returns the number of nonzero values in the vector v. */
template<class T>
int rsNumNonZeros(const std::vector<T>& v)
{
  return rsArrayTools::numNonZeros(&v[0], (int)v.size()); 
}

template<class T>
T rsSum(const std::vector<T>& x)
{
  T sum = T(0);
  for(size_t i = 0; i < x.size(); i++)
    sum += x[i];
  return sum;
}

template<class T>
T rsMean(const std::vector<T>& x)
{
  return rsSum(x) / T(x.size());
}

template<class T>
void rsScale(std::vector<T>& x, T scaler)
{
  for(size_t i = 0; i < x.size(); i++)
    x[i] *= scaler;
}


/** Converts C-array to std::vector. */
template<class T>
inline std::vector<T> toVector(const T* theArray, const size_t size) // rename to rsToVector
{
  std::vector<T> v(size);
  for(size_t i = 0; i < size; i++)
    v[i] = theArray[i];
  return v;
}

template<class T>
inline std::vector<T> rsToVector(const T* theArray, const int size) 
{
  std::vector<T> v(size);
  for(int i = 0; i < size; i++)
    v[i] = theArray[i];
  return v;
}


/** Copies data from existing C-array into an existing std::vector */
template<class T>
inline void rsCopyToVector(const T* a, int N, std::vector<T>& v)
{
  v.resize(N);
  for(int i = 0; i < N; i++)
    v[i] = a[i];
}

// ...hmm - defining operators for std::vector could clash with other libraries that do the same
// ...maybe i should use a wrapper rsArray or rsVector that wraps std::vector and has the same 
// interface (plus soem extra stuff) - maybe we should wrap it into an #ifdef 
// RS_USE_STD_VECTOR_OPERATORS


/** Unary minus for std::vector. Negates all elements. */
template<class T>
inline std::vector<T> operator-(const std::vector<T>& v)
{
  std::vector<T> result(v.size());
  for(size_t i = 0; i < v.size(); i++)
    result[i] = -v[i];
  return result;
}

/** Multiplies a scalar and a vector. */
template<class T>
inline std::vector<T> operator*(const T& x, const std::vector<T>& v)
{
  std::vector<T> result(v.size());
  for(size_t i = 0; i < v.size(); i++)
    result[i] = x * v[i];
  return result;
}

/** Multiplies a vector and a scalar. */
template<class T>
inline std::vector<T> operator*(const std::vector<T>& v, const T& x)
{
  std::vector<T> result(v.size());
  for(size_t i = 0; i < v.size(); i++)
    result[i] = x * v[i];
  return result;
}


/** Divides a vector by a scalar. */
template<class T>
inline std::vector<T> operator/(const std::vector<T>& v, const T& x)
{
  std::vector<T> result(v.size());
  for(size_t i = 0; i < v.size(); i++)
    result[i] = v[i] / x;
  return result;
}


/** Divides a scalar by a vector. */
template<class T>
inline std::vector<T> operator/(const T& x, const std::vector<T>& v)
{
  std::vector<T> result(v.size());
  for(int i = 0; i < v.size(); i++)
    result[i] = x / v[i];
  return result;
}


/** Adds a scalar to a vector. */
template<class T>
inline std::vector<T> operator+(const std::vector<T>& v, const T& x)
{
  std::vector<T> result(v.size());
  for(size_t i = 0; i < v.size(); i++)
    result[i] = v[i] + x;
  return result;
}

/** Subtracts a scalar from a vector. */
template<class T>
inline std::vector<T> operator-(const std::vector<T>& v, const T& x)
{
  std::vector<T> result(v.size());
  for(size_t i = 0; i < v.size(); i++)
    result[i] = v[i] - x;
  return result;
}

/** Adds two vectors element wise. */
template<class T>
inline std::vector<T> operator+(const std::vector<T>& x, const std::vector<T>& y)
{
  size_t Nmax = std::max(x.size(), y.size());
  size_t Nmin = std::min(x.size(), y.size());
  std::vector<T> result(Nmax);
  for(size_t i = 0; i < Nmin; i++)
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
  for(size_t i = 0; i < Nmin; i++)
    result[i] = x[i] - y[i];
  return result;
}

/** Multiplies two vectors element wise. */
template<class T>
inline std::vector<T> operator*(const std::vector<T>& x, const std::vector<T>& y)
{
  size_t Nmax = std::max(x.size(), y.size());
  size_t Nmin = std::min(x.size(), y.size());
  std::vector<T> result(Nmax);
  for(int i = 0; i < Nmin; i++)
    result[i] = x[i] * y[i];
  return result;
}

/*
template<class T>
inline bool operator==(const std::vector<T>& x, const std::vector<T>& y)
{
  if(x.size() != y.size())
    return false;
  for(size_t i = 0; i < x.size(); i++)
    if(x[i] != y[i])
      return false;
  return true;
}
*/


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
