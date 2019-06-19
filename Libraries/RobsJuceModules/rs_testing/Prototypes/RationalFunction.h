#pragma once

template<class T>
class rsRationalFunction
{

public:


  /** Evaluates the function at the given input x. */
  T operator()(T x) const
  {
    return num(x) / den(x);
  }





  //===============================================================================================
  /** \name Computations on raw coefficient arrays */

  static void partialFractionExpansion(
    std::complex<T> *numerator, int numeratorDegree,
    std::complex<T> *denominator, int denominatorDegree,
    std::complex<T> *poles, int *multiplicities, int numDistinctPoles,
    std::complex<T> *pfeCoeffs);
  // allocate heap memory


protected:

  rsPolynomial<T> num, den;  // numerator and denominator polynomials

};

/*
todo:
-implement a conversion constructor that can take a polynomial and make a rationla function from it

*/