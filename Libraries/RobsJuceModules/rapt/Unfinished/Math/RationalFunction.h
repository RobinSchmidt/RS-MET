#pragma once

//namespace RAPT
//{

/** A class for representing rational functions R(x) = P(x) / Q(x) where P and Q are both 
polynomials in x. Rational functions are important in signal processing because the transfer 
functions of linear time invariant systems (a.k.a. filters) are of that type. The class provides 
facilities for performing arithmetic with rational functions (including composition), evaluation,
partial fraction expansion, etc. 

Note that instances of this class are not suitable for realtime code - there's a lot of memory 
allocation going on in the operators (std::vectors are created, resized, assigned, etc.) and the
algorithms are not optimized either.

to verify:
I think, the rational functions over a field (like the real or complex numbers) form themselves a 
field - as opposed to the polynomials, which form only a ring. That's why we can always do a 
division of two rational functions without having to worry about a remainder (which may occur in 
polynomial division). */

template<class T>
class rsRationalFunction
{

public:


  rsRationalFunction() {}

  rsRationalFunction(
    const std::vector<T>& numeratorCoeffs, const std::vector<T>& denominatorCoeffs)
    : num(numeratorCoeffs), den(denominatorCoeffs)
  {

  }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  // setNumerator, setDenominator // for polynomial and std::vector and maybe plain arrays

  bool reduce(T tol);



  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Adds two rational functions. */
  rsRationalFunction<T> operator+(const rsRationalFunction<T>& q) const
  {
    rsRationalFunction<T> r;
    ratAdd(num.coeffs, den.coeffs, q.num.coeffs, q.den.coeffs, r.num.coeffs, r.den.coeffs);
    return r;
  }

  /** Subtracts two rational functions. */
  rsRationalFunction<T> operator-(const rsRationalFunction<T>& q) const
  {
    rsRationalFunction<T> r;
    ratAdd(num.coeffs, den.coeffs, q.num.coeffs, q.den.coeffs, r.num.coeffs, r.den.coeffs, 
      T(0), T(1), T(-1));
    return r;
  }

  /** Multiplies two rational functions. */
  rsRationalFunction<T> operator*(const rsRationalFunction<T>& q) const
  {
    rsRationalFunction<T> r;
    ratMul(num.coeffs, den.coeffs, q.num.coeffs, q.den.coeffs, r.num.coeffs, r.den.coeffs);
    return r;
  }

  /** Divides two rational functions. */
  rsRationalFunction<T> operator/(const rsRationalFunction<T>& q) const
  {
    rsRationalFunction<T> r;
    ratDiv(num.coeffs, den.coeffs, q.num.coeffs, q.den.coeffs, r.num.coeffs, r.den.coeffs);
    return r;
    // implement conversion operator in rsPolynomial to get rid of accessing coeffs
  }

  /** Returns the rational function that results from nesting/composing the given inner rational 
  function with this function as outer function. You can use it like h = f(g) where h,f,g are all 
  rsRationalFunction objects. */
  rsRationalFunction<T> operator()(rsRationalFunction<T> inner) const
  {
    rsRationalFunction<T> r;
    ratNest(inner.num.coeffs, inner.den.coeffs, num.coeffs, den.coeffs, 
      r.num.coeffs, r.den.coeffs);
    return r;
  }

  /** Evaluates the function at the given input x. */
  T operator()(T x) const
  {
    return num(x) / den(x);
  }

  /** Compares this rational function to another one for equality. */
  bool operator==(const rsRationalFunction<T>& q) const 
  { 
    return num == q.num && den == q.den;
  }




  //===============================================================================================
  /** \name Computations on raw coefficient arrays */

  static void partialFractionExpansion(
    std::complex<T>* numerator, int numeratorDegree,
    std::complex<T>* denominator, int denominatorDegree,
    std::complex<T>* poles, int* multiplicities, int numDistinctPoles,
    std::complex<T>* pfeCoeffs, std::complex<T>* polyCoeffs = nullptr);
  // allocates heap memory
  // ToDo:
  // -comment the format of the output - how does it deal with repeated poles? see unit tests
  // -maybe have another complex array *polyCoeffs which will contain the polynomial part when
  //  numeratorDegree >= denominatorDegree
  // -have a higher-level version of the function that doesn't require the poles to be passed (the
  //  function should find them itself via a root finder)



protected:

  rsPolynomial<T> num, den;  // numerator and denominator polynomials

};


//}

/*
todo:
-implement a conversion constructor that can take a polynomial and make a rationla function from it

*/