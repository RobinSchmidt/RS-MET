#ifndef RAPT_MODULARINTEGER_H
#define RAPT_MODULARINTEGER_H

/** This is a class for representing integers modulo some modulus m. The arithmetic operators 
always return values between 0...m-1 (inclusive) and it is assumed (and not checked), that the user
always initializes the value inside that range as well (for example, in constructors and
assignments).

Note that, if you weant to use the division operator, you should instatiate the template-class
with a signed integer type (this is because computation of the modular inverse relies on the
extended Euclid algorithm, which operates with signed values) - maybe we should get rid of the
division operator anyway...

class is not yet tested  */

template<class T>
class rsModularInteger
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Default constructor. Leaves value and modulus uninitialized. */
  rsModularInteger() {}

  /** Constructor. You may initialize the number by passing some unsigned 64-bit integer. */
  rsModularInteger(rsUint64 initialValue, rsUint64 modulusToUse);

  /** Copy constructor. */
  rsModularInteger(const rsModularInteger& other);


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets up this modular integer to a new value and modulus. The modulus must be > 1. The value
  will be stored in its canonical representation. For example, when the modulus is 5, it doesn't 
  matter if you pass 3 or 8 or 13,... or -2 or -7 or -12,... for the value. It will always be 
  stored as 3. */
  void set(T newValue, T newModulus);








  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  rsModularInteger& operator=(const rsModularInteger& b)
  { value = b.value; modulus = b.modulus; return *this; }

  rsModularInteger operator-() const;

  bool operator==(const rsModularInteger& other) const;
  bool operator!=(const rsModularInteger& other) const;

  rsModularInteger operator+(const rsModularInteger& other);
  rsModularInteger operator-(const rsModularInteger& other);
  rsModularInteger operator*(const rsModularInteger& other);
  rsModularInteger operator/(const rsModularInteger& other);

  rsModularInteger& operator+=(const rsModularInteger& other);
  rsModularInteger& operator-=(const rsModularInteger& other);
  rsModularInteger& operator*=(const rsModularInteger& other);
  rsModularInteger& operator/=(const rsModularInteger& other);

  rsModularInteger& operator++();
  rsModularInteger& operator--();


  //-----------------------------------------------------------------------------------------------
  /** \name Misc */


  /** Implements the modulo operation in a way that works also for when the value is negative.
  The result of the C++ operator % is equal to modulo-result only when the left operand is >= 0.
  ToDo: Maybe allow m to be negative, too.  */
  static T modulo (T x, T m)
  {
    rsAssert(m > 1);   // modulus must be positive integer >= 2
    T r = x % m;       // division remainder
    if(r < 0)          // r < 0 happens for x < 0
      return r + m;
    return r;
    // See https://stackoverflow.com/questions/11720656/modulo-operation-with-negative-numbers
    // There's also code for when m is negative, but we don't need that here. What would that even
    // mean? Maybe move the function into the library as rsModulo.
  };
  // maybe move to .cpp



  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  T value, modulus;  
  // ToDo: make protected, provide accessors. Reason: Setters should canonicalize the 
  // representation, similar to rsFraction

protected:

  void canonicalize();


};

/** Explicit template instantiation, to be used by the rsPow template-function and also in 
rsPolynomial etc.. */
template<class T>
rsModularInteger<T> rsUnityValue(rsModularInteger<T> value)
{ 
  return rsModularInteger<T>(T(1), value.modulus); 
}

template<class T> 
rsModularInteger<T> rsConstantValue(T value, rsModularInteger<T> targetTemplate) 
{ 
  return rsModularInteger<T>(value, targetTemplate.modulus);
}

// todo: add also rsZeroValue

#endif
