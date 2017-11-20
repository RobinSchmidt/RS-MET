#ifndef RS_POLYNOMIAL_INL
#define RS_POLYNOMIAL_INL

using namespace RSLib;

// construction/destruction:

template<class T>
rsPolynomial<T>::rsPolynomial(T *coeffs, int order)
{
  N = order;
  a = new T[N+1];
  for(int n = 0; n <= N; n++)
    a[n] = coeffs[n];
}

template<class T>
rsPolynomial<T>::rsPolynomial(const rsPolynomial<T>& other)
{
  N = other.N;
  a = new T[N+1];
  for(int n = 0; n <= N; n++)
    a[n] = other.a[n];
}

template<class T>
rsPolynomial<T>::~rsPolynomial()
{
  delete[] a;
}



/*
// operators:

template<class T>
rsPolynomial<T> rsPolynomial<T>::operator-() const
{
  if( value == 0 )
    return *this;
  else
    return rsPolynomial<T>(modulus-value, modulus);
}

template<class T>
bool rsPolynomial<T>::operator==(const rsPolynomial<T>& other) const
{
  return value == other.value && modulus == other.modulus;
}

template<class T>
bool rsPolynomial<T>::operator!=(const rsPolynomial<T>& other) const
{
  return !(*this == other);
}

template<class T>
rsPolynomial<T> rsPolynomial<T>::operator+(const rsPolynomial<T> &other)
{
  rsAssert( modulus == other.modulus );
  T r = this->value + other.value;
  if( r >= modulus )
    r -= modulus;
  return rsPolynomial<T>(r, modulus);
}

template<class T>
rsPolynomial<T> rsPolynomial<T>::operator-(const rsPolynomial<T> &other)
{
  rsAssert( modulus == other.modulus );
  T r;
  if( other.value > this->value )
    r = modulus + this->value - other.value;
  else
    r = this->value - other.value;
  return rsPolynomial<T>(r, modulus);
}

template<class T>
rsPolynomial<T> rsPolynomial<T>::operator*(const rsPolynomial<T> &other)
{
  rsAssert( modulus == other.modulus );
  T r = (this->value * other.value) % modulus;
  return rsPolynomial<T>(r, modulus);
}

template<class T>
rsPolynomial<T> rsPolynomial<T>::operator/(const rsPolynomial<T> &other)
{
  rsAssert( modulus == other.modulus );
  return *this * rsModularInverse(other.value, modulus);
}

template<class T>
rsPolynomial<T>& rsPolynomial<T>::operator+=(const rsPolynomial<T> &other)
{
  *this = *this + other;
  return *this;
}
template<class T>
rsPolynomial<T>& rsPolynomial<T>::operator-=(const rsPolynomial<T> &other)
{
  *this = *this - other;
  return *this;
}
template<class T>
rsPolynomial<T>& rsPolynomial<T>::operator*=(const rsPolynomial<T> &other)
{
  *this = *this * other;
  return *this;
}
template<class T>
rsPolynomial<T>& rsPolynomial<T>::operator/=(const rsPolynomial<T> &other)
{
  *this = *this / other;
  return *this;
}
*/

#endif
