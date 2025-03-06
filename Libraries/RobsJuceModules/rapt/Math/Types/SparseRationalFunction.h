#ifndef RAPT_SPARSERATIONALFUNCTION_H
#define RAPT_SPARSERATIONALFUNCTION_H

//=================================================================================================

// Under construction

/** Implements a sparse rational function. ...TBC... */

template<class T>
class rsSparseRationalFunction
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  /** Numerator and denominator polynomials. These data members are public because it's really more
  convenient that way. We could do some sort of facade pattern and delegation but it would 
  literally just be boilerplate - here and in client code - and a lot of it. We would have to 
  implement functions like:

    void setNumeratorCoeff(int index, T newCoeff) { num.setCoeff(index, newCoeff); }

  and then client code would call things like:

    r.setNumeratorCoeff(...);

  instead of:

    r.num.setCoeff(...);

  We would have to do this basically for all setters and getters of rsSparsePolynomial - twice. 
  And there are a lot of setters and getters. Nope. Just nope! Let's make num and den public 
  instead. Yes, I'm fully aware that it goes against OOP encapsulation practices. I know the rules 
  and break them deliberately here. We do not really have to maintain any class invariants or 
  anything like that so it's ok to let client code directly access and manipulate the numerator and 
  denominator. The only donwside may be that the variable names "num" and "den" now become part of
  the public API of the class and can't be changed later. I can live with that. */
  rsSparsePolynomial<T> num, den;
  // ...well...wait: There actually is a class invariant that (maybe) should be maintained: The 
  // denominator should be nonzero...hmmm...well...or maybe we just take the position that the onus 
  // is on the client to avoid divisions by zero. That's actually also how it works for int and 
  // float. Such variables (and also rsFraction) also do not nanny the programmer that way. So why
  // should we? Or maybe I'm just too lazy to write the boilerplate and trying to rationalize it? 
  // But it's not just about writing the boilerplate. It's also about readability and bloat - not 
  // on the binary code side (the delegations would be inlined) but on the source code side. 
  //
  // Hmm...but maybe the assumption that we really want expose all the setters and getters for the
  // two polynomials is wrong? Maybe we actually want to deal with a higher level interface here?
  // If really access to the full functionality of rsSparsPolynomial is needed, we could provide
  // getters like getNumerator/DenominatorReference() for that.
  //
  // We'll see.....



  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */


  rsSparseRationalFunction() 
  {
    den._appendTerm(T(1), 0);
  }


  rsSparseRationalFunction(
    const rsSparsePolynomial<T>& numerator, const rsSparsePolynomial<T>& denominator) 
    : num(numerator), den(denominator)  {}

  // Maybe implement it for const rsSparsePolynomial<T>&&, too. Or maybe that one is enough?


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Clears numerator and denominator. Note that this puts the object into an invalid state. It 
  would formally represent the indeterminate expression 0/0. So, this function should be used with
  great care and perhaps only internally in low level code. That's why it's marked with an 
  underscore. In higher level code, consider using initToZero() instead which sets the function to
  the zero function f(x) = 0 which is quite probably what you actually want to achieve anyway. */
  void _clear()
  {
    num.clear();
    den.clear();
  }
  // Maybe it should have an underscore? It puts the function itno an invalid state representing
  // the function 0/0 (I think)

  /** Initializes this rational function to the zero function: f(x) = 0. We represent this as
  f(x) = 0*x^0 / 1*x^0. The array for the numerator will be empty and the array for the 
  denominator will have a single coefficient of 1 for the monomial x^0. */
  void initToZero()
  {
    _clear();                  // f(x) = 0/0. That's indeterminate!
    den._appendTerm(T(1), 0);  // f(x) = 0/1. That's much better.
  }

  /** Initializes this rational function that is constantly one: f(x) = 1. We represent this as
  f(x) = 1*x^0 / 1*x^0. */
  void initToOne()
  {
    initToZero();              // f(x) = 0/1
    num._appendTerm(T(1), 0);  // f(x) = 1/1
  }

  /*
  void initToIdentity()
  {
    initToZero();              // f(x) = 0/1
    num._appendTerm(T(1), 1);  // f(x) = x/1
  }
  */


  void _setNumTerms(int newNumNumeratorTerms, int newNumDenominatorTerms)
  {
    num._setNumTerms(newNumNumeratorTerms);
    den._setNumTerms(newNumDenominatorTerms);
  }
  // Maybe it should have an underscore! It's a low level method. It may put the object into an
  // undefined state

  void setupFromDenseCoeffs(
    const std::vector<T>& newNumeratorCoeffs,
    const std::vector<T>& newDenominatorCoeffs,
    T tol)
  {
    num.setupFromDenseCoeffs(newNumeratorCoeffs,   tol);
    den.setupFromDenseCoeffs(newDenominatorCoeffs, tol);
  }

  void setupFromDenseCoeffs(
    const T* newNumeratorCoeffs,   int newNumNumeratorTerms, 
    const T* newDenominatorCoeffs, int newNumDenominatorTerms,
    T tol)
  {
    num.setupFromDenseCoeffs(newNumeratorCoeffs,   newNumNumeratorTerms,   tol);
    den.setupFromDenseCoeffs(newDenominatorCoeffs, newNumDenominatorTerms, tol);
  }


  void copyDataFrom(const rsSparseRationalFunction<T>& q)
  {
    num.copyDataFrom(q.num);
    den.copyDataFrom(q.den);
  }



  /** Applies a scaling factor to this rational function. This basically means to scale all 
  numerator coeffs by that factor. */
  void scale(T scaler) { num.scale(scaler); }


  void multiplyBy(rsMonomial<T> factor) { num.multiplyBy(factor); }


  void multiplyBy(const rsSparseRationalFunction<T>& factor, T tol) 
  { 
    num.multiplyBy(factor.num, tol);
    den.multiplyBy(factor.den, tol);

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


  void multiplyByDenseCoeffs(const T* numeratorCoeffs,   int numNumeratorTerms,
                             const T* denominatorCoeffs, int numDenominatorTerms, T tol)
  {
    num.multiplyByDenseCoeffs(numeratorCoeffs,   numNumeratorTerms,   tol);
    den.multiplyByDenseCoeffs(denominatorCoeffs, numDenominatorTerms, tol);
  }
  // Needs test


  // Maybe also make a divideByDenseCoeffs function by just calling multiplyByDenseCoeffs with
  // swapped arguments...or maybe not - client code can do that itself - no need to increase the
  // API surface area


  





  void addConstant(T constant, T tol)
  {
    num.addScaledPolynomial(den, constant, tol);
  }
  // ToDo: Document why this formula is right:
  //
  // N(x) / D(x) + k = N(x) / D(x) + k * D(x) / D(x) = (N(x) + k*D(x)) / D(x)
  //
  // Document, if this works in place (I think so)


  //void canonicalize();
  // Should: (1) Divide out the GCD of num and den. (2) Canonicalize num and den. 
  // (3) Divide num and den by the leading coeff of den (i.e. make den monic)
  //
  // But maybe it could also make sense to define a canonical representation as one with monic
  // numerator? But no! This can't represent the zero function.


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  bool isCloseTo(const rsSparseRationalFunction<T>& q, T tol) const
  {
    return q.num.isCloseTo(num, tol) && q.den.isCloseTo(den, tol);
  }


  // isCanonical()
  // A canonical representation has canonical numerator and denominator with no common factors
  // and the denominator is monic



  rsSparsePolynomial<T>& getNumerator() { return num; }


  rsSparsePolynomial<T>& getDenominator() { return den; }




  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Evaluates the function at the given input x. */
  T operator()(T x) const { return num(x) / den(x); }

  /** Evaluates the function at the given input z whose type may be different from the 
  coefficient type, for example, for evaluating functions with real coeffs at complex arguments.
  WARNING: the same considerations as for @see rsPolynomial::operator(TArg) apply. */
  template<class TArg>
  TArg operator()(TArg z) const { return num(z) / den(z); }

  /** Adds two rational functions. */
  rsSparseRationalFunction<T> operator+(const rsSparseRationalFunction<T>& q) const 
  { rsSparseRationalFunction<T> r; weightedSum(*this, T(1), q, T(1), &r, T(0)); return r; }

  /** Subtracts two rational functions. */
  rsSparseRationalFunction<T> operator-(const rsSparseRationalFunction<T>& q) const 
  { rsSparseRationalFunction<T> r; weightedSum(*this, T(1), q, T(-1), &r, T(0)); return r; }

  /** Multiplies two rational functions. */
  rsSparseRationalFunction<T> operator*(const rsSparseRationalFunction<T>& q) const 
  { return rsSparseRationalFunction(num * q.num, den * q.den); }

  /** Divides two rational functions. */
  rsSparseRationalFunction<T> operator/(const rsSparseRationalFunction<T>& q) const 
  { return rsSparseRationalFunction(num * q.den, den * q.num); }


  //-----------------------------------------------------------------------------------------------
  /** \name Low level API.  */


  static void weightedSum(
    const rsSparseRationalFunction<T>& p, T wp,
    const rsSparseRationalFunction<T>& q, T wq,
    rsSparseRationalFunction<T>* r, T tol);

  static void weightedSumDestructive(
    rsSparseRationalFunction<T>* p, T wp,
    rsSparseRationalFunction<T>* q, T wq,
    rsSparseRationalFunction<T>* r, T tol);
  // The first parameter p may alias to the result r. 


};

/** Multiplies a coefficient and a sparse rational function. */
template<class T>
inline rsSparseRationalFunction<T> operator*(const T& s, const rsSparseRationalFunction<T>& p)
{
  rsSparseRationalFunction<T> r(p);
  r.scale(s);
  return r;
}

template<class T>
void rsSparseRationalFunction<T>::weightedSum(
  const rsSparseRationalFunction<T>& p, T wp,
  const rsSparseRationalFunction<T>& q, T wq,
  rsSparseRationalFunction<T>* r, T tol)
{
  r->den = p.den * q.den;
  r->num = wp * p.num * q.den  +  wq * q.num * p.den;

  // This can probably be optimized with respect to avoid unnecessary temporary objects and heap
  // allocations. We may also use the gcd instead of just cross-mutiplying the denominators.
}

template<class T>
void rsSparseRationalFunction<T>::weightedSumDestructive(
  rsSparseRationalFunction<T>* p, T wp,
  rsSparseRationalFunction<T>* q, T wq,
  rsSparseRationalFunction<T>* r, T tol)

{
  rsAssert(rsAreAddressesDistinct(*p, *q));
  //rsAssert(rsAreAddressesDistinct(*r, *p));  // We may actually allow this!
  rsAssert(rsAreAddressesDistinct(*r, *q));
  // Maybe we can relax this? It would be really nice if r could be equal to at least one of p or
  // q. Requiring p and q to be distinct is not such a big problem. I think, we can allow this, if
  // SP::weightedSum can work in place. But at the moment. I think, it can't. But maybe it can be
  // made so. Looking at the code, it seems like it could work when the result aliases to the 1st 
  // argument. Test and document this! A test indicates that this may indeed work out. Investigate
  // this further and document! We could perhaps make it work to also allow r == q by swapping p 
  // and q (and wp and wq) in this case. But what if r == p == q? ...well...in that case, we could 
  // leave the denominator of r (and p and q) alone and just multiply the numerator by the scaler 
  // (wp+wq), I think. I think, the p == q != r case could possibly also be handled - just copy p
  // or q into r and then scale by (wp+wq)

  using SP = rsSparsePolynomial<T>;

  //              arg1        arg2        result
  SP::multiply(   p->num,     q->den,     &p->num, tol);  // Replace p->num by p->num * q->den
  SP::multiply(   q->num,     p->den,     &q->num, tol);  // Replace q->num by q->num * p->den
  SP::weightedSum(p->num, wp, q->num, wq, &r->num, tol);  // Establish r->num
  SP::multiply(   p->den,     q->den,     &r->den, tol);  // Establish r->den

  // ToDo:
  //
  // - Document exactly, how it can be used with respect to which pointers must be distinct and 
  //   which one may alias (and to what)
}
// Needs tests











#endif