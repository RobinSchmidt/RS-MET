template<class T>
rsComplexExponentialIterator<T>::rsComplexExponentialIterator(std::complex<T> a, std::complex<T> z)
{
  this->w = a;
  this->z = z;
}

template<class T>
void rsComplexExponentialIterator<T>::resetValue(std::complex<T> initialValue)
{ 
  w = initialValue;
}

template<class T>
void rsComplexExponentialIterator<T>::setZ(std::complex<T> newZ)
{
  z = newZ;
}

//=================================================================================================

/*
template<class T>
rsSineIterator<T>::rsSineIterator()
{
  a1 =  1.0806046117362795;
  s1 = -0.84147098480789650;
  s2 = -0.90929742682568171;
    // calling setup(1, 0, 1) would compute these values, but that would be more costly.
}

template<class T>
rsSineIterator<T>::rsSineIterator(T w, T p, T a)
{
  setup(w, p, a);
}
*/

template<class T>
void rsSineIterator<T>::setup(T w, T p, T a)
{
  a1 = 2.0*cos(w);
  s1 = a*sin(p-    w);
  s2 = a*sin(p-2.0*w);
  // Try to optimize computations of s1,s2 using addition theorems - i think, we may save one call
  // to sin
}
