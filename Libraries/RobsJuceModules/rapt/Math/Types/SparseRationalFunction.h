#ifndef RAPT_SPARSERATIONALFUNCTION_H
#define RAPT_SPARSERATIONALFUNCTION_H

//=================================================================================================

/** Implements a sparse rational function, i.e. a quotient of two sparse polynomials. To represent
these, this class has two members of type rsSparsePolynomial. Just like rsSparsePolynomial, the 
class rsSparserationalFunction maintains the representation in a canonical form as long as you use
only the high level interface. The low level interface is again signaled by method names that start
with an underscore. In a canonical representation, we require that the two constituting polynomials
are themselves in canonical representation and moreover, we require that numerator and denominator 
have no common factors (i.e. they are coprime) and the denominator is monic (i.e. has leading 
coefficient of 1).  ...TBC...   */

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
    : num(numerator), den(denominator)  
  {
  
  }
  // Maybe it should take an (optional?) tolerance parameter and call canonicalize()


  // Maybe implement it for const rsSparsePolynomial<T>&&, too. Or maybe that one is then enough, 
  // i.e. can also accept lvalue references? I think, an rvalue reference parameter can accept 
  // both kinds of arguments: rvalue- and lvalue references but lvalue reference parameters can
  // only accept lvalue reference arguments. Verify!

  rsSparseRationalFunction(
    const std::vector<T>& numeratorCoeffs,
    const std::vector<T>& denominatorCoeffs,
    const TTol& tolerance)
    : num(numeratorCoeffs, tolerance), den(denominatorCoeffs, tolerance)
  {

  }


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the numerical tolerance that is used to determine if a coefficient should be 
  considered zero, i.e. with in the numerical roundoff noise. */
  void setRoundoffTolerance(TTol newTolerance) 
  { num.setRoundoffTolerance(newTolerance); den.setRoundoffTolerance(newTolerance); }

  /** Sets this rational function to the zero function: f(x) = 0. We represent this as
  f(x) = 0*x^0 / 1*x^0. The array for the numerator will be empty and the array for the 
  denominator will have a single coefficient of 1 for the monomial x^0. */
  void setToZero() { _clear(); den._appendTerm(T(1), 0); }

  /** Sets this rational function to the function that is constantly one: f(x) = 1. We represent 
  this as f(x) = 1*x^0 / 1*x^0. */
  void setToOne() { setToZero(); num._appendTerm(T(1), 0); }

  /** Sets this rational function the identity function: f(x) = x. We represent this as
  f(x) = 1*x^1 / 1*x^0. */
  void setToIdentity() { setToZero(); num._appendTerm(T(1), 1); }
  // Needs test

  /** Sets this rational function the power function: f(x) = x^p. We represent this as
  f(x) = 1*x^p / 1*x^0. */
  void setToPower(int p) 
  { 
    //rsAssert(p >= 0);
    setToZero(); 
    if(p < 0)
      den._appendTerm(T(1), -p);
    else
      num._appendTerm(T(1),  p); 
  }
  // Needs test. 




  /*
  void multiplyBy(rsMonomial<T> factor) 
  { 
    _multiplyBy(factor);
    _canonicalize();
    // Calling this always may be overkill! I think, we need it only when the denominator has no
    // constant term. Only then we may factor out an x from the denominator (which is also a factor
    // of "factor"). Also, even in this case, we may not need all steps of the canonicalization. We
    // May actually use a custom optimized algorithm that just determines the lowest power of x in
    // num and den and reduces all powers by that number, I think.
    //
    // This now triggers an assert!
  }
  */


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  TTol getRoundoffTolerance() const
  { return rsMax(num.getRoundoffTolerance(), den.getRoundoffTolerance()); }

  bool isCloseTo(const SparseRatFunc& q, TTol tol) const
  { return q.num.isCloseTo(num, tol) && q.den.isCloseTo(den, tol); }

  bool isCloseTo(const SparseRatFunc& q) const
  { return isCloseTo(q, rsMax(getRoundoffTolerance(), q.getRoundoffTolerance())); }
  // Needs tests

  bool isZero() const { return num.isZero(); }




  const rsSparsePolynomial<T, TTol>& getNumeratorConst() const { return num; }

  const rsSparsePolynomial<T, TTol>& getDenominatorConst() const { return den; }

  int getNumNumeratorTerms()   const { return num.getNumTerms(); }

  int getNumDenominatorTerms() const { return den.getNumTerms(); }


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Evaluates the function at the given input x. */
  //T operator()(T x) const { return num(x) / den(x); }

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

  void _setupFromDenseCoeffs(const std::vector<T>& newNumeratorCoeffs,
                             const std::vector<T>& newDenominatorCoeffs)
  {
    num.setupFromDenseCoeffs(newNumeratorCoeffs);
    den.setupFromDenseCoeffs(newDenominatorCoeffs);
  }

  void _setupFromDenseCoeffs(const T* newNumeratorCoeffs,   int newNumNumeratorTerms, 
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



  // Prefix by underscore:
  rsSparsePolynomial<T, TTol>& _getNumerator() { return num; }

  rsSparsePolynomial<T, TTol>& _getDenominator() { return den; }





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
  // Needs test


  // Maybe also make a divideByDenseCoeffs function by just calling multiplyByDenseCoeffs with
  // swapped arguments...or maybe not - client code can do that itself - no need to increase the
  // API surface area


  /** Reduces this rational function to lowest terms. That means, it divides out the greatest 
  common divisor of numerator and denominator from both. This doesn't change the represented 
  rational function mathematically.*/
  void _reduce();
  // Needs tests

  /** Canonicalizes numerator and denominator. This doesn't change the represented rational 
  function mathematically. */
  void _canonicalizeNumAndDen() { num._canonicalize(); den._canonicalize(); }
  // Needs tests

  /** Makes our denominator monic by dividing out the leading coefficient of the denominator from
  both, numerator and denominator. This doesn't change the represented rational function 
  mathematically. */
  void _makeDenominatorMonic()
  { T s = T(1) / den.getLeadingCoeff(); num._scaleCoeffs(s); den._scaleCoeffs(s); }
  // Needs tests

  /** Puts this rational function into its canonical representation. That means it will be reduced
  to lowest terms, numerator and denominator will be in canonical representation and the 
  denominator will be monic. It will also make sure that the roundoff tolerances of numerator and
  denominator match (if they don't match already, it will pick the maximum of both). */
  void _canonicalize();
  // Needs tests

  /** Returns true iff this rational function is in canonical representation. A canonical 
  representation has canonical numerator and denominator with no common factors (i.e. they are 
  coprime) and the denominator is monic (i.e. has leading coeff 1). */
  bool _isCanonical() const;
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


#endif