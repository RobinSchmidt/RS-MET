#pragma once

namespace RAPT
{

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
      0, 1, -1);
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

  /** Evaluates the function at the given input x. */
  T operator()(T x) const
  {
    return num(x) / den(x);
  }

  /** Returns the rational function that results from nesting/composing the given inner rational 
  function this function as outer function. You can use it like h = f(g) where h,f,g are all 
  rsRationalFunction objects. */
  rsRationalFunction<T> operator()(rsRationalFunction<T> inner) const
  {
    rsRationalFunction<T> r;
    ratNest(inner.num.coeffs, inner.den.coeffs, num.coeffs, den.coeffs, 
      r.num.coeffs, r.den.coeffs);
    return r;
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
    std::complex<T>* pfeCoeffs);
  // allocate heap memory


protected:

  //T tol = 0.0;
  rsPolynomial<T> num, den;  // numerator and denominator polynomials

};


}

/*
todo:
-implement a conversion constructor that can take a polynomial and make a rationla function from it

*/