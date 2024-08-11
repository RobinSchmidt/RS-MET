#ifndef RAPT_MODULARINTEGER_H
#define RAPT_MODULARINTEGER_H

/** This is a class for representing integers modulo some modulus m. The value will always be in
the range 0...m-1. If you pass a value outside this range to a constructor or to set(), then the 
class will automatically wrap it into the allowed range such that even in this case, the invariant
is maintained.

Note that, if you weant to use the division operator, you should instatiate the template-class
with a signed integer type. This is because computation of the modular inverse relies on the
extended Euclid algorithm, which operates with signed values.


ToDo:

- Try to make it work also with unsigned integer types for T

- Add unit tests for division when the modulus is not a prime. In this case, there are zero 
  divisors, i.e. number a != 0 with the property that a * b = 0 for some b != 0. That makes 
  division by such numbers impossible - it basically behaves like division by zero. The resulting
  mathematical structure is then only a ring, not a field. Document this.

- Add more unit tests in general.    

- Compare implementation to rsFraction and maybe use the same style/conventions here with regard to
  what goes into the .h and what goes into the .cpp file, naming operator parameters, etc.

*/

template<class T>
class rsModularInteger
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Default constructor. On default construction, the value will be initialized to 0 and the 
  modulus to 2. That's the smallest modulus that makes sense. */
  rsModularInteger() {}

  /** Constructor. Creates a modular integer with given value and modulus. */
  rsModularInteger(const T& initialValue, const T& modulusToUse);

  /** Copy constructor. */
  rsModularInteger(const rsModularInteger& other);

  // ToDo: Copy/Move-assigment, etc.


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets up this modular integer to a new value and modulus. The modulus must be > 1. The value
  will be stored in its canonical representation. For example, when the modulus is 5, it doesn't 
  matter if you pass 3 or 8 or 13,... or -2 or -7 or -12,... for the value. It will always be 
  stored as 3. */
  void set(T newValue, T newModulus);


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry


  T getValue() const { return value; }

  T getModulus() const { return modulus; }

  /** Returns true, iff this modular integer has a multiplicative inverse. */
  bool hasInverse() const;
  // ToDo: Document the mathematical conditions for when this is the case. I think, the value must
  // be coprime with the modulus or something...


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  rsModularInteger& operator=(const rsModularInteger& b)
  { value = b.value; modulus = b.modulus; return *this; }

  rsModularInteger operator-() const;
  // Maybe add unary + as well. It's just the identity, i.e. trivial...but still. Some code 
  // explicitly writes out the plusses

  bool operator==(const rsModularInteger& other) const;
  bool operator!=(const rsModularInteger& other) const;

  // The size comparison operators have questionable semantics but are needed for some things:
  bool operator> (const rsModularInteger& other) const;
  bool operator< (const rsModularInteger& other) const;
  bool operator>=(const rsModularInteger& other) const;
  bool operator<=(const rsModularInteger& other) const;

  rsModularInteger operator+(const rsModularInteger& other) const;
  rsModularInteger operator-(const rsModularInteger& other) const;
  rsModularInteger operator*(const rsModularInteger& other) const;
  rsModularInteger operator/(const rsModularInteger& other) const;

  rsModularInteger& operator+=(const rsModularInteger& other);
  rsModularInteger& operator-=(const rsModularInteger& other);
  rsModularInteger& operator*=(const rsModularInteger& other);
  rsModularInteger& operator/=(const rsModularInteger& other);

  rsModularInteger& operator++();
  rsModularInteger& operator--();


protected:


  /** Makes sure that our value is in the range 0...modulus-1. */
  void canonicalize() { value = rsModulo(value, modulus); }

  /** Checks, if the representation is canonical. It's protected because we use it only internally 
  for assertions and unit tests. Client code is supposed to assume that the value is always 
  represented canonically or it shouldn't even care how it is represented, so there shouldn't be a
  situation where client code wants to check that condition. If it isn't, we have a bug. */
  bool isCanonical() const { return value >= 0 && value < modulus; }


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  T value   = T(0);
  T modulus = T(2);

};


// Explicit template instantiations, to be used by the rsPow template-function and also in 
// rsPolynomial etc.:

template<class T>
rsModularInteger<T> rsZeroValue(rsModularInteger<T> value)
{ 
  return rsModularInteger<T>(T(0), value.getModulus()); 
}

template<class T>
rsModularInteger<T> rsUnityValue(rsModularInteger<T> value)
{ 
  return rsModularInteger<T>(T(1), value.getModulus()); 
}

template<class T> 
rsModularInteger<T> rsIntValue(int value, rsModularInteger<T> targetTemplate) 
{ 
  return rsModularInteger<T>(T(value), targetTemplate.getModulus());
}

/*
template <class T>
bool rsGreaterAbs(const rsModularInteger<T>& lhs, const rsModularInteger<T>& rhs)
{
  rsAssert(lhs.getModulus() == rhs.getModulus());
  return lhs.getValue() > rhs.getValue();
}
*/

// ToDo: Implement the default, templatized rsZeroValue and rsUnityValue functions in terms of 
// rsIntValue such that rsModularInteger and all other similar classes need to provide only an 
// implementation of rsConstantValue

#endif
