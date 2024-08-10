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
  //rsModularInteger(rsUint64 initialValue, rsUint64 modulusToUse);
  // This is weird! Why does the constructor not just take a pair of type T? Try to get rid!



  rsModularInteger(const T& initialValue, const T& modulusToUse);

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
  // \name Inquiry


  T getValue() const { return value; }

  T getModulus() const { return modulus; }

  bool hasInverse() const;


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  rsModularInteger& operator=(const rsModularInteger& b)
  { value = b.value; modulus = b.modulus; return *this; }

  rsModularInteger operator-() const;

  bool operator==(const rsModularInteger& other) const;
  bool operator!=(const rsModularInteger& other) const;

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


  //-----------------------------------------------------------------------------------------------
  /** \name Misc */


  /** Implements the modulo operation in a way that works also for when the left operand (or first
  function parameter) is negative. The result of the C++ operator % is equal to modulo-result only 
  when the left operand is >= 0. ToDo: Maybe allow m to be negative, too - but what would a 
  negative modulus even mean mathematically? Does it even make sense? */
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
  // But what if T is unsigned?
  // maybe move to .cpp



  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  T value, modulus;  
  // ToDo: make protected, provide accessors. Reason: Setters should canonicalize the 
  // representation, similar to rsFraction

protected:

  void canonicalize() { value = modulo(value, modulus); }

  /** Checks, if the representation is canonical. It's protected because we use it only internally 
  for assertions and unit tests. Client code is supposed to assume that the value is always 
  represented canonically or it shouldn't even care how it is represented, so there shouldn't be a
  situation where client code wants to check that condition. If it isn't, we have a bug. */
  bool isCanonical() const { return value >= 0 && value < modulus; }
  // Maybe move out of the class - it hinders instantiation of the class for T = rsPolynomial. But
  // that alone will not solve the problem because we actually use it in an assertion in the 
  // copy-constructor. Maybe the solution is to implement the comparison operators in rsPolynomial
  // in some meaningful way. But how? Or get rid of the function and the assertion. Or maybe 
  // compile the assertion only conditionally if the type T implements <, <=. See:
  // https://en.cppreference.com/w/cpp/header/type_traits
  // An instantiation with rsPolynomial may be useful for implementing Galois fields (I think).
  // OK...done...needs tests...

};

/** Explicit template instantiations, to be used by the rsPow template-function and also in 
rsPolynomial etc.. */

template<class T>
rsModularInteger<T> rsZeroValue(rsModularInteger<T> value)
{ 
  return rsModularInteger<T>(T(0), value.modulus); 
}

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

// ToDo: Implement the default, templatized rsZeroValue and rsUnityValue functions in terms of 
// rsConstantValue such that rsModularInteger and all other similar classes need to provide only an 
// implementation of rsConstantValue

#endif
