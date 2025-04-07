#ifndef RAPT_SPARSERATIONALFUNCTION_H
#define RAPT_SPARSERATIONALFUNCTION_H

//=================================================================================================

// Under construction

/** Implements a sparse rational function. ...TBC... */

//template<class T>
template<class T, class TTol = rsEmptyType>
class rsSparseRationalFunction
{

public:


  // For convenience:
  using SparsePoly    = rsSparsePolynomial<T, TTol>; 
  using SparseRatFunc = rsSparseRationalFunction<T, TTol>;


  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Default constructor. Creates the rational function that is constantly zero. The zero function
  is canonically represented as the rational function  f(x) = 0 * x^0 / 1 * x^0  where the zero in
  the numerator is canonically represented as empty coefficient array. */
  rsSparseRationalFunction() { den._appendTerm(T(1), 0); }

  /** Constructor that converts a number c to the constnat function that just produces c for any 
  input. */
  rsSparseRationalFunction(const T& c)
  {
    if(!rsIsZero(c))           // The zero sparse polynomial is canonically represented as empty,
      num._appendTerm(c, 0);   // ..so we append the c*x^0 term only if c is nonzero.
    den._appendTerm(T(1), 0);  // The denominator is always one.
  }

  rsSparseRationalFunction(
    const SparsePoly& numerator, const SparsePoly& denominator) 
    : num(numerator), den(denominator)  {}

  // Maybe implement it for const rsSparsePolynomial<T>&&, too. Or maybe that one is then enough, 
  // i.e. can also accept lvalue references? I think, an rvalue reference parameter can accept 
  // both kinds of arguments: rvalue- and lvalue references but lvalue reference parameters can
  // only accept lvalue reference arguments. Verify!


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the numerical tolerance that is used to determine if a coefficient should be 
  considered zero, i.e. with in the numerical roundoff noise. */
  void setRoundoffTolerance(TTol newTolerance) 
  { num.setRoundoffTolerance(newTolerance); den.setRoundoffTolerance(newTolerance); }

  /** Initializes this rational function to the zero function: f(x) = 0. We represent this as
  f(x) = 0*x^0 / 1*x^0. The array for the numerator will be empty and the array for the 
  denominator will have a single coefficient of 1 for the monomial x^0. */
  void setToZero()
  {
    _clear();                  // f(x) = 0/0. That's indeterminate!
    den._appendTerm(T(1), 0);  // f(x) = 0/1. That's much better.
  }

  /** Initializes this rational function that is constantly one: f(x) = 1. We represent this as
  f(x) = 1*x^0 / 1*x^0. */
  void setToOne()
  {
    setToZero();               // f(x) = 0/1
    num._appendTerm(T(1), 0);  // f(x) = 1/1
  }

  void setToIdentity()
  {
    setToZero();               // f(x) = 0/1
    num._appendTerm(T(1), 1);  // f(x) = x/1
  }





  /*
  void multiplyBy(rsMonomial<T> factor) 
  { 
    _multiplyBy(factor);
    _canonicalize();
    // Calling this always may be overkill! I think, we need it only when the denominator has no
    // constant term. Also, we may not need all steps of the canonicalization
    //
    // This now triggers an assert!
  }
  */












  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  bool isCloseTo(const SparseRatFunc& q, T tol) const
  { return q.num.isCloseTo(num, tol) && q.den.isCloseTo(den, tol); }


  bool isZero() const { return num.isZero(); }

  rsSparsePolynomial<T, TTol>& getNumerator() { return num; }

  rsSparsePolynomial<T, TTol>& getDenominator() { return den; }

  const rsSparsePolynomial<T, TTol>& getNumeratorConst() const { return num; }

  const rsSparsePolynomial<T, TTol>& getDenominatorConst() const { return den; }

  int getNumNumeratorTerms()   const { return num.getNumTerms(); }

  int getNumDenominatorTerms() const { return den.getNumTerms(); }



  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Evaluates the function at the given input x. */
  T operator()(T x) const { return num(x) / den(x); }

  /** Evaluates the function at the given input z whose type may be different from the 
  coefficient type, for example, for evaluating functions with real coeffs at complex arguments.
  WARNING: the same considerations as for @see rsPolynomial::operator(TArg) apply. */
  template<class TArg>
  TArg operator()(TArg z) const { return num(z) / den(z); }

  /** Negates this rational function. */
  SparseRatFunc operator-() const { return SparseRatFunc(-num, den); }

  /** Adds two rational functions. */
  SparseRatFunc operator+(const SparseRatFunc& q) const 
  { SparseRatFunc r; weightedSum(*this, T(1), q, T(1), &r, T(0)); return r; }

  /** Subtracts two rational functions. */
  SparseRatFunc operator-(const SparseRatFunc& q) const 
  { SparseRatFunc r; weightedSum(*this, T(1), q, T(-1), &r, T(0)); return r; }

  /** Multiplies two rational functions. */
  SparseRatFunc operator*(const SparseRatFunc& q) const 
  { return SparseRatFunc(num * q.num, den * q.den); }

  /** Divides two rational functions. */
  SparseRatFunc operator/(const SparseRatFunc& q) const 
  { return SparseRatFunc(num * q.den, den * q.num); }


  //-----------------------------------------------------------------------------------------------
  /** \name Boilerplate */

  SparseRatFunc& operator+=(const SparseRatFunc& b) { return *this = (*this) + b; }





  //-----------------------------------------------------------------------------------------------
  /** \name Low level API.  
  ToDo: Maybe for some of these functions, provide versions without underscore. They should call
  the underscore version and then take appropriate action to ensure a canonical representation.
  In the simplest case, this may mean to just call canonicalize(). But this is expensive and in 
  certain cases, it may be possible to get away with a cheaper method so we need to figure out 
  what is strictly necessarry in each case and then do only that. */

  static void weightedSum(const SparseRatFunc& p, T wp, const SparseRatFunc& q, T wq, 
    SparseRatFunc* r, T tol);
  // Get rid of tol.

  static void weightedSumDestructive(SparseRatFunc* p, T wp, SparseRatFunc* q, T wq,
    SparseRatFunc* r, T tol);
  // Get rid of tol.
  // Document that The first parameter p may alias to the result r.

  /** Clears numerator and denominator. Note that this puts the object into an invalid state. It 
  would formally represent the indeterminate expression 0/0. So, this function should be used with
  great care and perhaps only internally in low level code. That's why it's marked with an 
  underscore. In higher level code, consider using setToZero() instead which sets the function to
  the zero function f(x) = 0 which is quite probably what you actually want to achieve anyway. */
  void _clear() { num.clear(); den.clear(); }

  void _copyDataFrom(const SparseRatFunc& q) { num = q.num; den = q.den; }
  // Maybe get rid of it and use (default) assignment operator instead.

  void _setupFromDenseCoeffs(
    const std::vector<T>& newNumeratorCoeffs,
    const std::vector<T>& newDenominatorCoeffs)
  {
    num.setupFromDenseCoeffs(newNumeratorCoeffs);
    den.setupFromDenseCoeffs(newDenominatorCoeffs);
  }

  void _setupFromDenseCoeffs(
    const T* newNumeratorCoeffs,   int newNumNumeratorTerms, 
    const T* newDenominatorCoeffs, int newNumDenominatorTerms)
  {
    num.setupFromDenseCoeffs(newNumeratorCoeffs,   newNumNumeratorTerms);
    den.setupFromDenseCoeffs(newDenominatorCoeffs, newNumDenominatorTerms);
  }

  void _setNumTerms(int newNumNumeratorTerms, int newNumDenominatorTerms)
  { num._setNumTerms(newNumNumeratorTerms); den._setNumTerms(newNumDenominatorTerms); }

  void _setNumeratorTerm(int index, const T& newCoeff, int power)
  { num._setTerm(index, newCoeff, power); }

  void _setDenominatorTerm(int index, const T& newCoeff, int power)
  { den._setTerm(index, newCoeff, power); }


  /** Applies a scaling factor to this rational function. This basically means to scale all 
  numerator coeffs by that factor. */
  void _scale(T scaler) { num._scaleCoeffs(scaler); }
  // May decanonicalize when scaler is zero

  /** Adds the given constant c to the rational function. This has the effect of adding a scaled
  version of the denominator to the numerator because N/D + c = N/D + c*D/D = (N + c*D)/D. */
  void _addConstant(T c) { num.addScaled(den, c); }
  // ToDo: Document, if this works in place (I think so)
  // Maybe it needs an underscore? Could it destroy the "reduced" property? I think so - but figure
  // out hwo and document this!



  void _multiplyBy(rsMonomial<T> factor) { num._multiplyBy(factor); }
  // Give it an underscore - why? could it possibly destroy canonicalness? I guess, it could 
  // destroy the no-common-factors (aka irreducibility) property. Maybe we should call a reduce()
  // function
  // Yes - it can destroy the reduced property. Consider R(x) = ((x+1)*(x+2)) / ((x+3)*x). It's
  // canonical. Num and den have no common factors. But if we multiply the numerator by the 
  // monomial x, they will have the common factor x. I think, this occurs whenever the denominator
  // has a monomial as factor, i.e. a factor of x^p, i.e. a root (possibly with multiplicity) at 
  // x = 0. This is equivalent to den not having a constant term. 


  void _multiplyBy(const SparseRatFunc& factor) 
  { 
    num.multiplyBy(factor.num);
    den.multiplyBy(factor.den);

    // I think, num and den may now have a common factor, so we potentially need to divide that
    // out:
    //canonicalize();
    // But maybe we should not automatically reduce the result by default because doing so may 
    // require memory allocations (because the GCD algo needs temporaries) and we really need this
    // function to be realtime safe (it's used in rsDampedCombAllpass::getCombTransferFunction(), 
    // for example - and that is used in rsDampedMultiCombAllpass::updateFilters() which could 
    // potentially be called on an audio thread). Maybe we should have a boolean parameter 
    // "reduce"? that deafults to true but that we set to false in a realtime context?)

    // Consider 14/15 * 3/4 = 42/60 = 7/10. Although both factors are in lowest terms, their 
    // product is not. The same thing could happen with rational functions.
  }
  // Get rid of tol, give it an underscore
  // I think, it will destroy the reduced feature if den and factor have a common factor, i.e.
  // theri gcd isn't 1. ...Verify this!

  void _multiplyByDenseCoeffs(const T* numeratorCoeffs,   int numNumeratorTerms,
    const T* denominatorCoeffs, int numDenominatorTerms)
  {
    num.multiplyByDenseCoeffs(numeratorCoeffs,   numNumeratorTerms);
    den.multiplyByDenseCoeffs(denominatorCoeffs, numDenominatorTerms);
  }
  // Needs test, get rid of tol param, give it an underscore


  // Maybe also make a divideByDenseCoeffs function by just calling multiplyByDenseCoeffs with
  // swapped arguments...or maybe not - client code can do that itself - no need to increase the
  // API surface area






  /** Reduces this rational function to lowest terms. That means, it divides out the greatest 
  common divisor of numerator and denominator from both. This doesn't change the represented 
  rational function mathematically.*/
  void _reduce();

  /** Canonicalizes numerator and denominator. This doesn't change the represented rational 
  function mathematically. */
  void _canonicalizeNumAndDen() { num._canonicalize(); den._canonicalize(); }

  /** Makes our denominator monic by dividing out the leading coefficient of the denominator from
  both, numerator and denominator. This doesn't change the represented rational function 
  mathematically. */
  void _makeDenominatorMonic()
  { T s = T(1) / den.getLeadingCoeff(); num._scaleCoeffs(s); den._scaleCoeffs(s); }

  /** Puts this rational function into its canonical representation. That means it will be reduced
  to lowest terms, numerator and denominator will be in canonical representation and the 
  denominator will be monic. */
  void _canonicalize() { _reduce(); _canonicalizeNumAndDen(); _makeDenominatorMonic(); }

  /** Returns true iff this rational function is in canonical representation. A canonical 
  representation has canonical numerator and denominator with no common factors (i.e. they are 
  coprime) and the denominator is monic (i.e. has leading coeff 1). */
  bool _isCanonical();
  // Needs tests.



protected:


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  rsSparsePolynomial<T, TTol> num, den;  // Numerator and denominator polynomials.

};


/** Multiplies a coefficient and a sparse rational function. */
template<class T, class TTol>
inline rsSparseRationalFunction<T, TTol> operator*(
  const T& s, const rsSparseRationalFunction<T, TTol>& p)
{
  rsSparseRationalFunction<T, TTol> r(p);
  r._scale(s);
  return r;
}


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
implement a check for that condition (and a few others) isCanonical(). We also provide functions to
invert the transfer function (basically, swapping numerator and denominator but maintaining the 
a0 = 1 condition by appropriate pre- and post scaling), reflecting the zeros about the unit circle 
(turning minimum phase filters into maximum phase ones), etc. ...TBC... 

ToDo: Move this class into its own dedicated pair of .h/.cpp files in Filters/General

*/


template<class T, class TTol>  // ToDo: Let TTol default to rsEmptyType
class rsSparseDigitalTransferFunction : public rsSparseRationalFunction<T, TTol>
{

public:

  using Base = rsSparseRationalFunction<T, TTol>;  // For convenience
  using Base::Base;                                // Inherit constructors



  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Adds an overall predelay to the whole filter by shifting all exponents of z^-1 by the given
  amount. */
  void addPreDelay(int amountInSamples)
  {
    if(amountInSamples < 0)
    {
      rsError("Negative predelay is not allowed");
      return;
      // We could allow it if we already have some predelay and the given amount would just 
      // reduce it. Maybe we can relax the restriction to amountInSamples >= -num.getMinPower() or
      // something. If the minimum power is 3, we could allow a predelay amount of -3. This would 
      // then just reduce the predelay to zero.
    }

    num.shiftPowers(amountInSamples);


    // Maybe do num.shiftPowers(rsMax(amountInSamples, -getPreDelay()) ); and write unit tests for
    // this
  }

  /** Removes the predelay from this filter, if any is present. This makes sure that the lowest
  exponent of z^-1 in the numerator is zero. see getPreDelay(), addPreDelay()  */
  void removePreDelay() { num.shiftPowers(-getPreDelay()); }

  /** Turns the filter into its inverse. This basically amounts to swapping numerator and
  denominator and possibly applying some scaling of the coefficients if b0 != 1. */
  void invert()
  {
    rsAssert(isCanonical());

    // A filter with predelay cannot be inverted (at least not if it operates in realtime). The 
    // best thing we can do in this case is to invert the filter up to the predelay. Removing the
    // predelay ensures that the 0-th term in the numerator has power of 0, i.e. it's a  b0 * z^-0  
    // term:
    removePreDelay();
    rsAssert(num.getPower(0) == 0); // Numerator is already asserted to be non-empty in isCanoncial
    rsAssert(num.getCoeff(0) != 0); // ..so we can access the 0-th element without risk here

    // Swap numerator and denominator while maintaining the a0 = 0 normalization condition by 
    // appropriately scaling the numerator before and after the swap:
    T s = T(1) / num.getCoeff(0);  // Desired scaler s is 1/b0
    _scale(s);                     // Scale old numerator to achieve b0 = 1 before swap
    std::swap(num, den);           // Swap numerator and denominator. b0 is now 1 because a0 was.
    _scale(s);                     // Scale new numerator to achieve desired overall gain

    // I'm pretty sure it doesnt' allocate. The swap of the underyling std::vectors should use move 
    // semantics. Verify and document this. We need to be able to call this function on a realtime 
    // thread in some damped comb allpass filters, so it is important that this function is 
    // non-allocating.
  }

  /** Reflects the zeros of the filter about the unit circle. This will turn a minimum phase
  filter into a maximum phase one and vice versa. For mixed phase filters, it inverts the mix.
  It doesn't affect stability or filter order. */
  void reflectZeros()
  {
    int deg = num.getDegree();                   // ToDo: use canonical getDegree
    for(int i = 0; i < num.getNumTerms(); i++)
      num._setPower(i, deg - num.getPower(i));
    num._reverse();                               // Order array by ascending powers again

    // How about a reflectPoles() function? But that would turn stable filters into unstable ones,
    // so it's usefulness is questionable. For the time being, we can do without. Maybe we should 
    // also apply complex conjugation of the coeffs in case of complex coeffs? If so, maybe
    // do it in a function num.conjugateCoeffs() which just calls rsConj() on each coeff (which is
    // an empty function for real types)
  }


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

  /** Performs some sanity checks. Is meant for debug assertions. */
  bool isCanonical() const
  {
    bool ok = true;

    // Numerator and denominator polynomials should not be empty:
    ok &= num.getNumTerms() > 0 && den.getNumTerms() > 0;
    // But is that right? Maybe an empty numerator is admissible to represent the zero function? 
    // It's not very useful as a transfer function, but still...

    // We assume the filter polynomials to be in canonical shape:
    ok &= num._isCanonical();
    ok &= den._isCanonical();

    // Filter should satisfy the a0 == 1 normalization property:
    ok &= den.getPower(0) == 0 && den.getCoeff(0) == T(1);

    return ok;
  }


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





  T operator()(T x) const 
  { 
    T xr = T(1) / x;
    return num(xr) / den(xr); 
  }
  // Do we really need this when we already have the variant with TArg below? Maybe we need it only
  // because the baseclass also has it and we need to override that? If so, maybe delete it here 
  // *and* in the baseclass and keep only the variant with TArg in both


  template<class TArg>
  TArg operator()(TArg z) const 
  { 
    TArg zr = TArg(1) / z;
    return num(zr) / den(zr);
  }
  // Reciprocation of z needed because we store the coeffs of H(z^-1)

  // We also need to override the evaluateAt functions...or maybe get rid of them in the baseclass.
  // Oh - looks like, we don't have such functions there. OK - good.

  // This boilerplate is needed to have the desired arithmetic operators available also for the 
  // derived class. They can not be inherited from the baseclass because their parameter and
  // return types are different. It may work with pointer-types but not with value-types (I guess):


  rsSparseDigitalTransferFunction<T, TTol> operator-() const
  {
    return rsSparseDigitalTransferFunction<T, TTol>(-num, den);
  }


  rsSparseDigitalTransferFunction<T, TTol> operator+(const rsSparseDigitalTransferFunction<T, TTol>& q) const 
  { rsSparseDigitalTransferFunction<T, TTol> r; weightedSum(*this, T(1), q, T(1), &r, T(0)); return r; }

  rsSparseDigitalTransferFunction<T, TTol> operator-(const rsSparseDigitalTransferFunction<T, TTol>& q) const 
  { rsSparseDigitalTransferFunction<T, TTol> r; weightedSum(*this, T(1), q, T(-1), &r, T(0)); return r; }

  rsSparseDigitalTransferFunction<T, TTol> operator*(const rsSparseDigitalTransferFunction<T, TTol>& q) const 
  { return rsSparseDigitalTransferFunction(num * q.num, den * q.den); }

  rsSparseDigitalTransferFunction<T, TTol> operator/(const rsSparseDigitalTransferFunction<T, TTol>& q) const 
  { return rsSparseDigitalTransferFunction(num * q.den, den * q.num); }

  // ToDo: Check if we need to reduce the results to lowest terms or maybe to "canonicalize".


};

/** Multiplies a coefficient and a sparse digital transfer function. */
template<class T, class TTol>
inline rsSparseDigitalTransferFunction<T, TTol> operator*(
  const T& s, const rsSparseDigitalTransferFunction<T, TTol>& p)
{
  rsSparseDigitalTransferFunction<T, TTol> r(p);
  r._scale(s);
  return r;
}




// Under construction:

// ToDo: Maybe define these for the baseclass instead - then verify, that the right functions are
// called:

template<class T, class TTol>
inline bool rsIsBetterPivot(
  const rsSparseDigitalTransferFunction<T, TTol>& x,
  const rsSparseDigitalTransferFunction<T, TTol>& y)
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
  const rsSparseDigitalTransferFunction<T, TTol>& x,
  const rsSparseDigitalTransferFunction<T, TTol>& tol)
{
  return x.isZero(); 
  // Maybe we should pass in a tolerance into this function? But then this tolerance should 
  // probably be just passed on from our second parameter. But for this, we need to make the 2nd
  // parameter of type T rather than rsSparseDigitalTransferFunction. I think, we generally need to
  // change the API to admit different types for x and tol
}

template<class T, class TTol>
inline rsSparseDigitalTransferFunction<T, TTol> rsGetPivotingTolerance(
  const rsMatrixView<rsSparseDigitalTransferFunction<T, TTol>>& A)
{
  return rsSparseDigitalTransferFunction<T, TTol>();
}



#endif