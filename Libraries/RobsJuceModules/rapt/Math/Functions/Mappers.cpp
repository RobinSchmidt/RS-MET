template<class T>
void rsCoordinateMapper<T>::setInputRange(T newMin, T newMax)
{
  inMin = newMin;
  inMax = newMax;
}

template<class T>
void rsCoordinateMapper<T>::setOutputRange(T newMin, T newMax)
{
  outMin = newMin;
  outMax = newMax;
}

template<class T>
void rsCoordinateMapper<T>::setLogScaled(bool shoulBeLogScaled)
{
  logScaled = shoulBeLogScaled;
}

template<class T>
T rsCoordinateMapper<T>::map(T x) const
{
  if(logScaled)
    return (T)rsExpToLin(x, inMin, inMax, outMin, outMax);
  return (T)rsLinToLin(x, inMin, inMax, outMin, outMax);
}

template<class T>
T rsCoordinateMapper<T>::unmap(T x) const
{
  if(logScaled)
    return (T)rsLinToExp(x, outMin, outMax, inMin, inMax);
  return (T)rsLinToLin(x, outMin, outMax, inMin, inMax);
}

// todo: get rid of code duplication and optimize (precompute coeffs)

//=================================================================================================

template<class T>
void rsCoordinateMapper2D<T>::setInputRange(T minX, T maxX, T minY, T maxY)
{
  mapperX.setInputRange(minX, maxX);
  mapperY.setInputRange(minY, maxY);
}

template<class T>
void rsCoordinateMapper2D<T>::setOutputRange(T minX, T maxX, T minY, T maxY)
{
  mapperX.setOutputRange(minX, maxX);
  mapperY.setOutputRange(minY, maxY);
}

template<class T>
void rsCoordinateMapper2D<T>::map(T *x, T *y) const
{
  *x = mapX(*x);
  *y = mapY(*y);
}

template<class T>
void rsCoordinateMapper2D<T>::unmap(T *x, T *y) const
{
  *x = unmapX(*x);
  *y = unmapY(*y);
}
