#ifndef RAPT_SPARSETRANSFERFUNCTION_H
#define RAPT_SPARSETRANSFERFUNCTION_H

//=================================================================================================

/** A subclass of rsSparseRationalFunction that is meant to deal specifically with transfer 
functions of digital filters. A general transfer function for a digital filter looks like:

          b0 + b1 * z^-1 + b2 * z^-2 + ... + bN * z^-N
  H(z) = ----------------------------------------------
          1  + a1 * z^-1 + a2 * z^-2 + ... + aM * z^-M

Note that this is actually a rational function not in z itself but in z^-1. We override the 
function evaluation operator () to take care of this reciprocation of z. We also implement some 
additional functionality on top of the baseclass that is specific to such transfer functions. For 
example, digital filter transfer functions are usually normalized to a0 = 1, as seen above. We 
implement a check for that condition (and a few others) in isCanonical(). We also provide functions
to invert the transfer function (basically, swapping numerator and denominator but maintaining the 
a0 = 1 condition by appropriate pre- and post scaling), reflecting the zeros about the unit circle 
(turning minimum phase filters into maximum phase ones), etc. ...TBC...   */

template<class T, class TTol = rsEmptyType> // Maybe rename T to TCoef
class rsSparseTransferFunction : public rsSparseRationalFunction<T, TTol>
{

public:

  using SparseTransFunc = rsSparseTransferFunction<T, TTol>;
  using Base            = rsSparseRationalFunction<T, TTol>;
  using Base::Base;                                             // Inherit constructors


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Adds an overall predelay to the whole filter by shifting all exponents of z^-1 by the given
  amount which must be a nonnegative integer. */
  void addPreDelay(int amountInSamples);

  /** Removes the predelay from this filter, if any is present. This makes sure that the lowest
  exponent of z^-1 in the numerator is zero. see getPreDelay(), addPreDelay()  */
  void removePreDelay() { num.shiftPowers(-getPreDelay()); }

  /** Turns the filter into its inverse. This basically amounts to swapping numerator and
  denominator and possibly applying some scaling of the coefficients if b0 != 1. */
  void invert();

  /** Reflects the zeros of the filter about the unit circle. This will turn a minimum phase
  filter into a maximum phase one and vice versa. For mixed phase filters, it inverts the mix.
  It doesn't affect stability or filter order. */
  void reflectZeros();


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the predelay introduced by this filter. A predelay is characterized by the fact that
  the lowest exponent of z^-1 in the numerator is not zero. That means the numerator does not look
  like b0 + b1*z-^1 + b2*z^-2 + ... but rather something like b7*z^-7 + b8*z^-8 + b9*z^-9 + ...
  for a predelay of 7 samples, for example. */
  int getPreDelay() const { return num.getPower(0); }

  /** Returns the order of the filter. This is the maximum exponent of z^-1 that occurs in the
  transfer function. */
  int getFilterOrder() const { return rsMax(num.getDegree(), den.getDegree()); }

  /** Checks if this transfer function is in canonical representation. A transfer function for a 
  digital filter in canonical representation has canonical numerator and denominator and the 
  constant term of the denominator is 1. Note how these criteria are different from a canoncial
  representation of general (sparse) rational functions which require the denominator to be monic,
  i.e. the leading coeff rather than the constant coeff must be 1. Also, we do not require the 
  numerator and denominator to be coprime here. That means we allow pole-zero cancellations in the
  transfer function and still consider it canonical. The rationale is that such pole-zero 
  cancellations may indeed occur in practice for certain perfectly valid parameter settings of user
  adjustable filters. Not having to check for coprimality has the nice side effect that we can 
  perform the test for canonicalness without the need to allocate temporary heap memory (which the
  coprimality test needs for computing the gcd). The function here overrides a non-virtual(!) 
  baseclass method, i.e. provides compile time polymorphism for the _isCanonical() member 
  function. There is no runtime polymorphism, though - so take care! */
  bool _isCanonical() const;

  /** Computes the density of the numerator defined as the number of actual nonzero coeffs divided
  by the number of potentially nonzero coeffs given the degree of the numerator. */
  double getNumeratorDensity() const
  { return double(num.getNumTerms()) / double(num.getDegree()+1); }

  /** Computes the density of the denominator defined as the number of actual nonzero coeffs 
  divided by the number of potentially nonzero coeffs given the degree of the denominator. */
  double getDenominatorDensity() const
  { return double(den.getNumTerms()-1) / double(den.getDegree()); }

  /** Returns the "combined density" defined as the number of actual nonzero coeffs of the filter 
  divided by the number of potential nonzero coeffs for the given filter order. The a0 coeff 
  doesn't count because it's always 1. */
  double getCombinedDensity() const
  {
    int numPossibleCoeffs = 2*getFilterOrder() + 1;
    int numActualCoeffs   = num.getNumTerms() + (den.getNumTerms()-1);
    return double(numActualCoeffs) / double(numPossibleCoeffs);
  }

  /** Returns the "separated density" defined as the as number of actual nonzero coeffs in 
  numerator and denominator divided by the number of potential nozero coeffs for the given orders
  of numerator and denominator. */
  double getSeparatedDensity() const
  {
    int numPossibleCoeffs = num.getDegree()+1 + den.getDegree();
    int numActualCoeffs   = num.getNumTerms() + (den.getNumTerms()-1);
    return double(numActualCoeffs) / double(numPossibleCoeffs);
  }
  // ToDo: Explain this better! We consider numerator and denominator as separate filters that can
  // have their own orders and therefore the computation of the number of possibly nonzero coeffs
  // is different. For example, in an Nth order allpole filter, we have a 0th order numerator. In 
  // getCombinedDensity, we would assume that it could potentially have N+1 coeffs. Here, we assume
  // that it can have only one because the order of the numerator is zero.

  // Maybe this could be moved into the baseclass. But I'm not sure, if the +1 fo the num and
  // no +1 for the den also applies there...well...I think, it does when we assume a canonical
  // representation with monic denomionator. Here we normalize the denominator to a0=1 - but
  // whatever way we normalize the function, the denominator has always one degree of freedom
  // less than the numerator (for the same degree). A degree 1 polynomial has 2 coeffs and a degree
  // N polynomial has N+1 coeffs. That number applies to the numerator as is. But the denominator
  // is normalized so we lose one degree of freedom and subtract 1 again.

  // ToDo: Maybe return the densities as rsFraction<int>. They are rational numbers so maybe we 
  // should treat them as such.

  // ToDo: Document use cases for these getDensity() functions.


  /** Overrides the function evaluation operator in order to reciprocate the input z before 
  applying the rational function to it. This is needed because we store the coeffs of H(z^-1) 
  rather than of H(z). That means the numerator and denominator polynomials have coeffs that 
  multiply powers of z^-1 such as z^-1, z^-2, z^-3, etc. and not z^1, z^2, z^3, etc.. Note that 
  the override works only at compile time. We don't support runtime polymorphism here. */
  template<class TArg>
  TArg operator()(TArg z) const { TArg zr = TArg(1) / z; return num(zr) / den(zr); }


  // This boilerplate is needed to have the desired arithmetic operators available also for the 
  // derived class. They can not be inherited from the baseclass because their parameter and
  // return types are different. It may work with pointer-types but not with value-types (I guess):


  SparseTransFunc operator-() const { return SparseTransFunc(-num, den); }

  SparseTransFunc operator+(const SparseTransFunc& q) const 
  { SparseTransFunc r; weightedSum(*this, T(1), q, T(1), &r, T(0)); return r; }

  SparseTransFunc operator-(const SparseTransFunc& q) const 
  { SparseTransFunc r; weightedSum(*this, T(1), q, T(-1), &r, T(0)); return r; }

  SparseTransFunc operator*(const SparseTransFunc& q) const 
  { return SparseTransFunc(num * q.num, den * q.den); }

  SparseTransFunc operator/(const SparseTransFunc& q) const 
  { return SparseTransFunc(num * q.den, den * q.num); }

  // ToDo: Check if we need to reduce the results to lowest terms or maybe to "canonicalize" ..but
  // we should not use the baseclass implementation of canonicalize here. I think, we should leave
  // out the reduce step and use a modified normaliation step. Also: implement +=, -=, *=, /=. The
  // implementations should avoid allocations of temporary objects.


};

/** Multiplies a coefficient and a sparse digital transfer function. */
template<class T, class TTol>
inline rsSparseTransferFunction<T, TTol> operator*(
  const T& s, const rsSparseTransferFunction<T, TTol>& p)
{
  rsSparseTransferFunction<T, TTol> r(p);
  r._scale(s);
  return r;
}


// Under construction:

// ToDo: Maybe define these for the baseclass instead - then verify, that the right functions are
// called:

template<class T, class TTol>
inline bool rsIsBetterPivot(
  const rsSparseTransferFunction<T, TTol>& x,
  const rsSparseTransferFunction<T, TTol>& y)
{ 
  // A zero x is never a better pivot than any y:
  if(x.isZero())
    return false;

  // Any nonzero x is always better than a zero y:
  if(y.isZero())
    return true;

  // An x with simpler denominator is better than a y with more complex denominator:
  if(x.getNumDenominatorTerms() < y.getNumDenominatorTerms())
    return true;

  // An x with more complex denominator is worse that a y with simpler denominator:
  if(x.getNumDenominatorTerms() > y.getNumDenominatorTerms())
    return false;

  // When x and y have the same denominator, we compare the numerators. A simpler numerator is
  // better (except when it's zero - but this has already been ruled out):
  return x.getNumNumeratorTerms() < y.getNumNumeratorTerms();
}

template<class T, class TTol>
inline bool rsIsInvalidDivisor(
  const rsSparseTransferFunction<T, TTol>& x,
  const rsSparseTransferFunction<T, TTol>& /*tol*/)
{
  return x.isZero(); 
  // Maybe we should pass in a tolerance into x.isZero()? But then this tolerance should 
  // probably be just passed on from our second parameter. But for this, we need to make the 2nd
  // parameter of type T (or more likely TTol) rather than 
  // rsSparseTransferFunction<T, TTol>. I think, we generally need to change the API to 
  // admit different types for x and tol - perhaps even the API auf the Gaussian eliminatio algo.
  // Maybe it also needs to take a tolerance parameter of type TTol.
}

template<class T, class TTol>
inline rsSparseTransferFunction<T, TTol> rsGetPivotingTolerance(
  const rsMatrixView<rsSparseTransferFunction<T, TTol>>& /*A*/)
{
  return rsSparseTransferFunction<T, TTol>();
  // I think, this should return a TTol. Maybe the return type should be set to "auto"
}





#endif