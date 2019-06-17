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

protected:

  rsPolynomial<T> num, den;  // numerator and denominator polynomials

};