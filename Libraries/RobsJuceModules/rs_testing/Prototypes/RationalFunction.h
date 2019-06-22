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