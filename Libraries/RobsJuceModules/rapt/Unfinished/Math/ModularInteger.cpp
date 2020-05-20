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
// what, if there is no modular inverse (i think, it exists only if the value is coprime with the 
// modulus - verify)

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

/*
Ideas:
-templatize on the integer type to use - allow for arbitrary size integers
-make it work also for negative values (-> verify, if the % operator works correctly, when the 
 left operand is negative - if not, use a custom function instead)
-how about negative moduli?
-currently, the arithmetic operations make sense only when the two operands have the same modulus
 -generalize this to a sort of "multi-modular" or "mixed-modular" arithmetic
 -the modulus of the result should be the lowest common multiple of the moduli of the operands 
 -would that make sense?
-does the notion of a modular rational number make any sense? i.e. numerator and/or denominator are 
 modular integers?

*/