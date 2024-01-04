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



/*=================================================================================================

Ideas:
-Make similar classes for 3D
-Make a class ParameterMapper. This is supposed to implement a monotonic function that maps the
 normalized parameter range 0..1 to some other range min..max in some possibly nonlinear 
 fashion. I already have such a class un jura_framework. Maybe the implementation should be 
 based on that and then the class in jura can use the new class in RAPT - or maybe it then 
 becomes obsolete entirely.
-A nice function to map 0..1 to 0..inf that feels natural should probably go through the points:
   x:  0.0  0.25  0.5  0.75  1.0
   y:  0.0  0.5   1.0  2.0   inf
 In a chat on the discord server "The Audio Programmer" a suggested function that works as 
 approximation (with some deviations from the desired (0.25,0.5) and (0.75,2) but otherwise exact)
 was  f(x) = (x / (1-x))^a  where  a = log_3(2) ~= 0.63. I don't know, how the value log_3(2) 
 arises and personally think a value close to a = 0.65 could be better in terms of error at the two
 non-exact points. This exponent can perhaps be numerically optimized to minimize
 (f(0.25) - 0.5)^2 + (f(0.75) - 2.0)^2

*/