
template<class T>
void rsCoordinateMapper2D<T>::map(T *x, T *y)
{
  *x = mapX(*x);
  *y = mapY(*y);
}

template<class T>
void rsCoordinateMapper2D<T>::unmap(T *x, T *y)
{
  *x = unmapX(*x);
  *y = unmapY(*y);
}

template<class T>
T rsCoordinateMapper2D<T>::mapX(T x)
{
  return x;
}

template<class T>
T rsCoordinateMapper2D<T>::mapY(T y)
{
  return y;
}

template<class T>
T rsCoordinateMapper2D<T>::unmapX(T x)
{
  return x;
}

template<class T>
T rsCoordinateMapper2D<T>::unmapY(T y)
{
  return y;
}