#ifndef RS_MODULARINTEGER_INL
#define RS_MODULARINTEGER_INL

#include "ModularInteger.h" // why is this include needed? maybe it's not
using namespace RSLib;

// construction/destruction:

template<class T>
rsModularInteger<T>::rsModularInteger(rsUint64 initialValue, rsUint64 modulusToUse)
{
  modulus = modulusToUse;
  value   = initialValue;
  rsAssert( value >= T(0) && value < modulus );
}

template<class T>
rsModularInteger<T>::rsModularInteger(const rsModularInteger<T>& other)
{
  modulus = other.modulus;
  value   = other.value;
}

// operators:

template<class T>
rsModularInteger<T> rsModularInteger<T>::operator-() const
{
  if( value == 0 )
    return *this;
  else
    return rsModularInteger<T>(modulus-value, modulus);
}

template<class T>
bool rsModularInteger<T>::operator==(const rsModularInteger<T>& other) const
{
  return value == other.value && modulus == other.modulus;
}

template<class T>
bool rsModularInteger<T>::operator!=(const rsModularInteger<T>& other) const
{
  return !(*this == other);
}

template<class T>
rsModularInteger<T> rsModularInteger<T>::operator+(const rsModularInteger<T> &other)
{
  rsAssert( modulus == other.modulus );
  T r = this->value + other.value;
  if( r >= modulus )
    r -= modulus;
  return rsModularInteger<T>(r, modulus);
}

template<class T>
rsModularInteger<T> rsModularInteger<T>::operator-(const rsModularInteger<T> &other)
{
  rsAssert( modulus == other.modulus );
  T r;
  if( other.value > this->value )
    r = modulus + this->value - other.value;
  else
    r = this->value - other.value;
  return rsModularInteger<T>(r, modulus);
}

template<class T>
rsModularInteger<T> rsModularInteger<T>::operator*(const rsModularInteger<T> &other)
{
  rsAssert( modulus == other.modulus );
  T r = (this->value * other.value) % modulus;
  return rsModularInteger<T>(r, modulus);
}

template<class T>
rsModularInteger<T> rsModularInteger<T>::operator/(const rsModularInteger<T> &other)
{
  rsAssert( modulus == other.modulus );
  return *this * rsModularInverse(other.value, modulus);
}

template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator+=(const rsModularInteger<T> &other)
{
  *this = *this + other;
  return *this;
}
template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator-=(const rsModularInteger<T> &other)
{
  *this = *this - other;
  return *this;
}
template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator*=(const rsModularInteger<T> &other)
{
  *this = *this * other;
  return *this;
}
template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator/=(const rsModularInteger<T> &other)
{
  *this = *this / other;
  return *this;
}
template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator++()
{
  *this = *this + rsModularInteger<T>(1, modulus);
  return *this;
}
template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator--()
{
  *this = *this - rsModularInteger<T>(1, modulus);
  return *this;
}

#endif
