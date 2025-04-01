#ifndef RAPT_SPARSEPOLYNOMIAL_H
#define RAPT_SPARSEPOLYNOMIAL_H


//=================================================================================================

/** A class for representing (univariate) monomials, i.e. expressions of the form  c * x^p  for
some coefficient c and integer power (or exponent) p. Strictly speaking, we should require p to be
nonnegative, but we don't really enforce this here.

See: https://en.wikipedia.org/wiki/Monomial   */

template<class T> 
class rsMonomial
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  explicit rsMonomial(T newCoeff = T(0), int newPower = 0) : coeff(newCoeff), power(newPower) { }
  // Marked as explicit because we want to avoid hidden automatic conversions from type T


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setup(T newCoeff, int newPower) { coeff = newCoeff; power = newPower; }

  void setCoeff(T newCoeff)   { coeff = newCoeff; }

  void setPower(int newPower) { power = newPower; }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the coefficient c in the expression c * x^p. */
  T getCoeff() const { return coeff; }

  /** Returns the power (aka exponent) p in the expression c * x^p. */
  int getPower() const { return power; }


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Evaluates the expression c * x^p at the given x. The data type of the argument x may be 
  different from the type T with which the class is instantiated. Can be used, for example, to 
  evaluate monomials with real coefficients at complex arguments. */
  template<class TArg>
  TArg operator()(TArg x) const { return TArg(coeff) * rsPow(x, TArg(power)); }
  // Preliminary. We may want to use rsPowInt for integer exponents. That may be more efficient.
  // Here, we explicitly first convert the exponent to type TArg and then call 
  // rsPow(TArg x, TArg y). But currently rsPowInt is not suitably defined. It expects two 
  // (unsigned?) ints. rsPowInt is currently only defined for x and power being both integers. We 
  // really need a function where the base is an arbitrary type and the epxonent is an integer.

  /** Returns the negative of this monomial. */
  rsMonomial<T> operator-() const { return rsMonomial<T>(-getCoeff(), getPower()); }

  /** Multiplies two monomials. */
  rsMonomial<T> operator*(const rsMonomial<T>& q) const
  { return rsMonomial<T>(getCoeff() * q.getCoeff(), getPower() + q.getPower()); }
  // Needs tests

  /** Divides two monomials. */
  rsMonomial<T> operator/(const rsMonomial<T>& q) const
  { return rsMonomial<T>(getCoeff() / q.getCoeff(), getPower() - q.getPower()); }
  // Needs tests. 
  // If q.power > this->power, this will lead to a negative power in the result. Should we do 
  // something about this like triggering an rsAssert? And what if q.getCoeff() returns zero?


protected:

  T   coeff = T(0);
  int power = T(0);

};


//=================================================================================================

/** A class for representing sparse polynomials, i.e. polynomials that have many zero coefficients.
Think of something like p(x) = 2*x^37 - 5*x^129 + 3*x^435. Such polynomials occur in filter 
transfer functions that use delaylines instead of unit delays. We represent such sparse polynomials
basically as a std::vector of monomials which we call terms in this context. We say that a sparse 
polynomial is in canonical representation if the powers of the terms are strictly increasing as 
function of array index (implying that no power appears more than once) and there are no terms with
a coefficient of zero (up to some tolerance - i.e. the absolute values of all coeffs should be 
greater than the tolerance). If the array of terms is empty, we treat that as the canonical 
representaion of the zero polynomial. This is different from how we represent the zero polynomial 
in the class rsPolynomial for dense polynomials. There, the coeff array of the zero polynomial has 
a single element whose value is zero. But because we want to avoid (near) zero coeffs here, this 
representation would not be a good fit because then we would have to make a special rule for the 
zeroth coeff which may not even exist in certain nonzero sparse polynomials like p(x) = 3*x^100.

For many purposes, it is convenient to assume a canonical representation and many setters will 
maintain such a representation - but not all of them. Sometimes, one needs to - at least 
temporarily - violate such a canonical representation for performance reasons. Therefore, the API 
has two levels. A higher level that assumes and maintains canonical representations and a lower 
level that makes no such assumption and gives no such maintenance guarantee. The lower level member
functions are prefixed with an underscore _ to indicate at the call site that now some low-level 
stuff is going on and special care should be taken. ...TBC...


ToDo:

- Remove the tol parameter from the non-static member functions. These should now use the tol 
  member for that purpose.

- Document clearly under which circumstances the user can assume the polynomial to be in a 
  canonical representation (and what that even means). I'm still not quite sure myself, whether or 
  not the API should always enforce a canonical representation as class invariant. Maintaining that
  at all times - in particluar when adding or modifying terms - is costly. On the other hand, 
  certain other operations (like extracting the leading term) are cheaper when we can assume a 
  canonical representation. At the moment a canonical representation is not enforced. ...TBC...

- Maybe prefix the low-level functions that may destroy a canonical representation by an 
  underscore. This signals at the call site that now the low-level API is being used and special 
  care is required, if maintaining a canonical representation is desired. ...done...

- Or: maybe add suffixes _c and _n to functions that work with canonical and non-canonical 
  representations. For setters, _c should mean that the function assumes the polynomial in 
  canonical representation as precondition *and* ensures that this still holds when the function 
  returns, i.e. as postcondition. For getters, only the precondition is relevant because they don't
  change the object. Functions with suffix _n do not assume such a precondition and even in the 
  case that the condition is met, they do not assure to maintain it. ...hmm...or maybe only mark
  the decanonicalizing methods somehow

- Sort the high-level and low-level access functions, i.e. let them have their own category like
  Setup (high level, maintaining canonical representation), Setup (low level, may destroy canonical
  representation), Inquiry (assuming canonical representation), Inquiry (not assuming canonical
  representation)

*/

template<class T, class TTol = rsEmptyType>
class rsSparsePolynomial
{

public:

  // ToDo: Maybe use this abbreviation for convenience in the member function declarations:
  using SparsePoly = rsSparsePolynomial<T, TTol>;

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Default constructor. Constructs an empty sparse polynomial. This represents, by definition,
  the zero polynomial. */
  rsSparsePolynomial() {}

  /** Creates a polynomial from an initializer list for the terms. */
  rsSparsePolynomial(std::initializer_list<rsMonomial<T>> initList, TTol tolerance) 
    : terms(initList), tol(tolerance) {}

  rsSparsePolynomial(const std::vector<T>& coefficients, TTol tolerance) 
  { 
    setupFromDenseCoeffs(coefficients, tolerance);

    // Maybe make the tolerance parameter optional. I'm not sure about that, though. It may invite
    // forgetting to set it when it's really needed. But on the other hand, some types T don't need
    // any tolerance at all.
  }


  // What about copy- and move constructors and copy- and move assignment operators? Do we need to
  // define them or can we rely on the auto-generated ones? It's important that swapping two
  // sparse polynomials can be done allocation free. This is needed for inverting sparse filters by
  // swapping numerator and denominator of their transfer functions (plus some extra stuff to 
  // maintain the a0 = 1 normalization). This is an an operation that we need to do in a realtime 
  // safe manner. Verify and document this!


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Reserves memory for the given number of terms. Can be called before calling functions like 
  addTerm() to pre-allocate the desired amount of memory beforehand when multiple terms are being
  added in a sequence. */
  void reserve(size_t numTerms) { terms.reserve(numTerms); }

  /** Clears the array of terms. */
  void clear() { terms.clear(); }

  /** Sets up the numerical tolerance that is used to determine if a coefficient should be 
  considered zero, i.e. with in the numerical roundoff noise. */
  void setRoundoffTolerance(TTol newTolerance) { tol = newTolerance; }

  /** Sets up the polynomial from a dense arrays of polynomial coeffs. When a coefficient in the 
  dense representation is zero, we not create a term for that. */
  void setupFromDenseCoeffs(const std::vector<T>& newCoeffs, TTol newTol)
  { setupFromDenseCoeffs(&newCoeffs[0], (int) newCoeffs.size(), newTol); }

  /** Like setupFromDenseCoeffs(const std::vector<T>&, ...) but for raw C-arrays. */
  void setupFromDenseCoeffs(const T* newCoeffs, int newNumTerms, TTol newTol);
  // ToDo: Provide methods that don't require a tol parameter. They should do the same thing 
  // except setting our tol member. 
  
  /** Adds the term c * x^p with coeff c and power p to the polynomial. If a term with the same 
  power already exists, this will just shift its coefficient. If the cofficient happens to be zero 
  after shift (up to the given tolerance), the term will be removed. */
  void addTerm(T coeff, int power);
  // This function assumes that the polynomial is in canonical representation! Document this and 
  // maybe reflect it in the function name.

  /** Adds the given monomial to the polynomial. */
  void addTerm(const rsMonomial<T>& newTerm)
  { addTerm(newTerm.getCoeff(), newTerm.getPower()); }

  /** Subtracts the given monomial from the polynomial. */
  void subtractTerm(const rsMonomial<T>& newTerm)
  { addTerm(-newTerm.getCoeff(), newTerm.getPower()); }

  /** Adds a scaled version of the given polynomial p to this polynomial. */
  void addScaledPolynomial(const SparsePoly& p, T scaler)
  {
    tol = rsMax(tol, p.tol);
    for(int i = 0; i < p.getNumTerms(); i++)
      addTerm(scaler * p.getCoeff(i), p.getPower(i));

    // Maybe it would be better to just append a scaled version and then canonicalize? This may 
    // result in less data movement - but it may blow up the required memory temporarily. So: no -
    // let's not do that in general. It may even lead to allocations when we really don't want 
    // them.
  }

  /** Scales all coeffs by the given scaler. */
  void scaleCoeffs(T scaler) { rsAssert(scaler != T(0)); _scaleCoeffs(scaler); }

  /** Alias for scaleCoeffs() for compatibility with API of rsPolynomial. */
  void scale(T scaler) { scaleCoeffs(scaler); }

  /** Makes the polynomial monic by dividing all coeffs by the leading coeff. A monic polynomial is
  a polynomial in which the leading coefficient is unity (aka one).*/
  void makeMonic() { scale(T(1) / _getLeadingCoeff()); }
  // Maybe use canonical getLeadingCoeff()

  /** Shifts all powers by the given amount. If the amount is p, this corresponds to multiplying 
  the polynomial by a monomial factor with unit coefficient, i.e. by x^p. */
  void shiftPowers(int amount)
  {
    for(int i = 0; i < getNumTerms(); i++)
      _shiftPower(i, amount);
  }
  // Shifting all powers by the same amount should be unproblematic with regard to 
  // decanonicalization. That's why this function doesn't need an underscore

  /** Multiplies this polynomial by the given monomial factor. This results in all coeffs being 
  multiplied by the coeff of the monomial and all powers being increased by the pwer of the 
  monomial. */
  void multiplyBy(const rsMonomial<T>& factor)
  { scaleCoeffs(factor.getCoeff()); shiftPowers(factor.getPower()); }

  /** Multiplies this polynomial by the given other polynomial factor. Works in place and 
  re-allocates only when the capacity is too low (VERIFY!). */
  void multiplyBy(const SparsePoly& factor, TTol tol) { multiply(*this, factor, this); }
  // I think, this may also decanonicalize! We may get multiple terms with same exponent. But we
  // may actually repair this inside the function. But no! It calls canonicalize at the end, so 
  // even if it temporarily decanonicalizes, it cleans everything up at the end.

  /** Multiplies this polynomial by a desne polynomial represented by the given array of 
  coefficients. Works in place and re-allocates only when the capacity is too low. */
  void multiplyByDenseCoeffs(const T* coeffs, int numTerms);
  // I think, this may also decanonicalize! See comment above. It's the same here

  void divideBy(const rsMonomial<T>& divisor)
  {
    scaleCoeffs(T(1) / divisor.getCoeff());
    shiftPowers(     - divisor.getPower());
  }
  // Needs tests. 
  // Maybe assert that this->getMinPower() >= divisor.getPower() to avoid producing negative 
  // powers.

  void addScaled(const SparsePoly& summand, const rsMonomial<T>& scaler);
  // ToDo: implement add(summand), i.e. the same thing but without the scaler.
  // ...and maybe one with the scaler being a simple coeff

  /** Turns the representation of the polynomial into a canonical one. A canonical representation 
  has the following properties: (1) The powers are strictly increasing as function of index. 
  (2) No power appears more than once. (3) No zero coefficients appear. We achieve this by 
  first sorting the terms, then consolidating multiple terms with equal exponents into single
  terms and finally deleting all terms that have a coefficient zero (up to the given tolerance). */
  void canonicalize();
  // Maybe it should also go into the low level interface section and get and underscore. Or maybe
  // it should even go into the protected section. ...but maybe client code sometimes needs it.


  void _copyDataFrom(const SparsePoly& other)
  {
    tol = other.tol;
    _setNumTerms(other.getNumTerms());
    for(int i = 0; i < getNumTerms(); i++)
      _setTerm(i, other.getCoeff(i), other.getPower(i));
  }
  // Maybe this should have an _ at the start. It will decanonicalize this polynomial, iff the 
  // other polynomial is in non-canonical representation. Hmmm...this is a gray area. if the client
  // code uses only other non-underscored function, this here may get away without underscore, too.
  // But maybe client code should use the assignment operator anyway (which we need to define - for
  // copy and move assignment)


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns true, iff this polynomial is empty, i.e. has no terms. */
  bool isEmpty() const { return terms.empty(); }
  // Maybe get rid. Clients should use isZero()

  /** Returns the numerical tolerance that is used to determine if a coefficient should be 
  considered zero, i.e. with in the numerical roundoff noise. */
  TTol getRoundoffTolerance() const { return tol; }

  /** Returns true iff this polynomial is the zero polynomial. */
  bool isZero() const { rsAssert(isCanonical()); return isEmpty(); }

  /** Returns true, iff the rhs polynomial equals this polynomial up to the given tolerance. */
  bool isCloseTo(const SparsePoly& rhs, TTol tol) const;

  /** Like isCloseTo above but without the tol parameter. It uses the maximum of our tol member 
  and the rhs's tol member. */
  bool isCloseTo(const SparsePoly& rhs) const { return isCloseTo(rhs, rsMax(tol, rhs.tol)); }

  /** Return true, iff the given index is valid, i.e. the object has a term with given index. */
  bool isValidIndex(int i) const { return i >= 0 && i < getNumTerms(); }

  /** Returns the number of terms in this polynomial. The i-th term is a monomial of the form 
  ci * x^pi with a coefficient ci and a power/exponent pi. */
  int getNumTerms() const { return (int) terms.size(); }

  /** Returns the term (i.e. the monomial) at the given index. */
  rsMonomial<T> getTerm(int index) const { rsAssert(isValidIndex(index)); return terms[index]; }

  /** Returns the coefficient of the term with given index. */
  T getCoeff(int index) const {  rsAssert(isValidIndex(index)); return terms[index].getCoeff(); }
  // ToDo: Return the coeff as const reference! We intend this class to be potentially used with
  // large coeff types like matrices or arbitrary precision floats, so this optimization may make 
  // sense

  /** Returns the power of the term with given index. */
  int getPower(int index) const { rsAssert(isValidIndex(index)); return terms[index].getPower(); }

  /** Checks if this sparse polynomial is in canonical representation. A representation is 
  canonical if it has no zero coefficients (up to a given tolerance) and if the powers are strictly
  increasing (as function of term-index). The empty polynomial is also accepted as a canonical 
  representation. It represents the zero polynomial. */
  bool isCanonical() const;


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */


  /** Copy assignment operator. Copies data from rhs into this object. */
  SparsePoly& operator=(const SparsePoly& rhs) 
  { 
    _copyDataFrom(rhs); 
    return *this; 
  }
  // ToDo: Implement move assignment


  /** Evaluates the function at the given input z whose type may be different from the 
  coefficient type T. This may be used, for example, for evaluating polynomials with real coeffs at
  complex arguments. 
  WARNING: the same considerations as for @see rsPolynomial::operator(TArg) apply. */
  template<class TArg>
  TArg operator()(TArg z) const
  { 
    TArg y(0);
    for(auto& t_i : terms)
      y += t_i(z);             // t_i(z) = coeffs[i] * z^powers[i]
    return y;
  }

  /** Implements the unary minus operator. */
  SparsePoly operator-() const 
  { 
    SparsePoly r;
    r._copyDataFrom(*this);
    r.scale(T(-1));  // Maybe use a special negate() function. It may be more efficient.
    return r;
  }
  // ToDo: Implement unary plus, too. It's trivial but sometimes, we may want to use it for 
  // clarity. But maybe it should return a (const?) reference rather than a value? Is that even 
  // possible?

  /** Adds two polynomials. */
  SparsePoly operator+(const SparsePoly& q) const 
  { SparsePoly r; add(*this, q, &r); return r; }

  /** Subtracts two polynomials. */
  SparsePoly operator-(const SparsePoly& q) const 
  { SparsePoly r; subtract(*this, q, &r); return r; }

  /** Multiplies two polynomials. */
  SparsePoly operator*(const SparsePoly& q) const 
  { SparsePoly r; multiply(*this, q, &r); return r; }

  /** Divides two polynomials. */
  SparsePoly operator/(const SparsePoly& q) const 
  { SparsePoly quot, rem; divide(*this, q, &quot, &rem); return quot; }

  /** Computes remainder of polynomial division, i.e. implements the modulo operation. */
  SparsePoly operator%(const SparsePoly& q) const 
  { SparsePoly quot, rem; divide(*this, q, &quot, &rem); return rem; }


  //-----------------------------------------------------------------------------------------------
  /** \name Static member functions */

  /** Computes the greatest common divisor of the polynomials p and q. The monic parameter defines 
  if the returned gcd should be normalized to be monic. The gcd of polynomials is unique only up to
  a constant scale factor, so it may make sense to make it well defined by requiring it to be
  monic. If monic is false, the returned gcd may be scaled by some arbitrary scale factor which
  depends on the details of the algorithm but has no mathematical significance (I think). But maybe
  it has? Figure out! */
  static SparsePoly greatestCommonDivisor(
    const SparsePoly& p, const SparsePoly& q, bool monic = true)
  {
    SparsePoly a = p, b = q, tmp1, tmp2;
    SparsePoly::greatestCommonDivisorInPlace(&a, &b, &tmp1, &tmp2, monic);
    return a;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Low level API. These functions operate on pre-allocated output parameters passed by 
  pointer (to make it obvious at the call site that the parameter may be modified). Using these 
  functions with pre-allocated sparse polynomials may potentially avoid heap allocations which the 
  more convenient functions that return sparse polynomials do. This includes the +,-,*,/,% 
  operators. So, for real time code, these operators are actually forbidden and one has to resort 
  to the low level API. */

  // These functions should not take a tol parameter. Instead the should assign the tolerance of
  // the result to the max of the tolerances of the operands:

  // Maybe move the functions with underscore into this area, too

  static void add(const SparsePoly& p, const SparsePoly& q, SparsePoly* r);

  static void subtract(const SparsePoly& p, const SparsePoly& q, SparsePoly* r);

  static void weightedSum(const SparsePoly& p, T wp, const SparsePoly& q, T wq, SparsePoly* r);

  /** Multiplies polynomials p and q and stores the result in r. It may be used in place, i.e. the
  result polynomial r can point to the memory location of the arguments p and/or q. */
  static void multiply(const SparsePoly& p, const SparsePoly& q, SparsePoly* r);

  /** Implements polynomial division with remainder. ...TBC... */
  static void divide(const SparsePoly& numerator, const SparsePoly& denominator,
    SparsePoly* quotient, SparsePoly* remainder);

  /** Computes the greatest common divisor of two polynomials. It works in place meaning that it
  allocates no temporary sparse polynomials internally. The first parameter is an input/output 
  parameter. On input it should contain the first argument. On return, it will contain the result,
  i.e. the GCD. The second parameter is the second argument and will also be used internally for 
  temporary data such that on return, it will be destroyed (it will be zeroed out by the 
  algorithm). The algorithm also needs two additional temporaries that you need to pass. Their 
  content on output is undefined. On input, they may contain anything - it doesn't matter. For 
  example usage, see the greatestCommonDivisor() function which basically serves as convenience 
  function for the in-place version. */
  static void greatestCommonDivisorInPlace(SparsePoly* FirstArgAndResult, SparsePoly* SecondArg,
    SparsePoly* temp1, SparsePoly* temp2, bool makeResultMonic);
  // I think, if all passed polynomials have large enough capacity, then the function should not
  // (re)allocate any heap memory. Verify and document this! How large is "large enough"?



  /** Appends a term with given coeff and power to the end of our terms array. Beware that this 
  may decanonicalize the representation. */
  void _appendTerm(T coeff, int power) { terms.emplace_back(rsMonomial<T>(coeff, power)); } 
  // May decanonicalize

  /** Sets the number of terms. If the new number is less than the current number, it will just 
  cut off terms from the end. If the new number is greater than the current number, it will just
  extend our vector of terms and the added terms at the end are uninitialized, i.e. may contain 
  garbage. This function should only be used if you intend to set up the new terms via e.g. 
  _setTerm() after calling _setNumTerms(). So, it's a function that needs a lot of care to be used
  properly. */
  void _setNumTerms(int newNumTerms) { terms.resize(newNumTerms); }

  /** Directly sets the coefficient and power of the term with given index with no regard for 
  maintaining a canonical representation. This is intended to be used in a sequence of calls with a 
  subsequent manual call to canonicalize when performance matters. */
  void _setTerm(int index, T coeff, int power) 
  { rsAssert(isValidIndex(index));  terms[index].setup(coeff, power); }
  // May decanonicalize in various ways.

  /** Directly sets the power of the term with given index with no regard for maintaining a 
  canonical representation. This is intended to be used in a sequence of calls with a subsequent 
  manual call to canonicalize when performance matters. */
  void _setPower(int index, int newPower)
  { rsAssert(isValidIndex(index)); terms[index].setPower(newPower); }
  // May decanonicalize by destroying the "all powers appear only once" property.

  /** Directly sets the coefficient of the term with given index with no regard for maintaining a 
  canonical representation. This is intended to be used in a sequence of calls with a subsequent 
  manual call to canonicalize when performance matters. */
  void _setCoeff(int index, T newCoeff)
  { rsAssert(isValidIndex(index)); terms[index].setCoeff(newCoeff); }

  /** Scales the coefficient with the given index by the given scaler. */
  void _scaleCoeff(int index, T scaler) { _setCoeff(index, scaler * getCoeff(index)); }
  // May decanonicalize if the scaler is zero.

  /** Scales all coefficients by the given scaler. */
  void _scaleCoeffs(T scaler)
  {
    for(int i = 0; i < getNumTerms(); i++)
      _scaleCoeff(i, scaler);
  }
  // Decanonicalizes when scaler == 0

  /** Shifts the coefficient with the given index by the given amount, i.e. adds the given amount 
  to the coeff. It may decanonicalize the representation by leading to a zero coeff. */
  void _shiftCoeff(int index, T amount) { _setCoeff(index, amount + getCoeff(index)); }

  /** Shifts the power at the given index by the given amount. It may decanonicalize the 
  representation by introducing two terms with equal power. */
  void _shiftPower(int index, int amount) { _setPower(index, amount + getPower(index)); }

  /** Reverses the array of terms. It may appear to be a weird thing to do on polynomials but this 
  operation is needed when transforming minimum phase filters into maximum phase ones (or vice 
  versa) and when producing allpass filters from allpole filters. */
  void _reverse() { rsReverse(terms); }

  /** Returns the minimum power that occurs in this polynomial. */
  int _getMinPower() const;
  // Implement a getMinPower() for canonical representations that just returns 0 or the power of
  // the 0-th term

  /** Returns the maximum power that occurs in this polynomial. In mathematical jargon, the 
  highest power in a polynomial is also known as the degree or order of the polynomial. */
  int _getMaxPower() const;
  // dito

  /** Returns the index of the maximum power or -1 in the case of an empty array of terms. */
  int _getMaxPowerIndex() const;
  // dito

  /** Alias for getMaxPower() for compatibility with API of rsPolynomial. Returns the degree of
  the polynomial. This is mathematical term for the term with the highest power/exponent that has 
  a nonzero coefficient. */
  int _getDegree() const { return _getMaxPower(); }
  // This is basically an alias name for getMaxPower(). I'm not sure, if it's a good idea to have 
  // two functions that do the exact same thing. Maybe get rid of it. But on the other hand, it's 
  // nice to have to be consistent with the API of class rsPolynomial. 

  /** Returns the leading coefficient, i.e. the coefficient that multiplies the highest power of
  the input variable x. */
  T _getLeadingCoeff() const;

  /** Returns the leading term in this polynomial, i.e. the monomial  cn x^n  that has the highest
  exponent n. */
  rsMonomial<T> _getLeadingTerm() const;




  // ToDo: 
  //
  // - Implement compose (see free function rsComposeNaive() in Prototypes.h file), 
  //   lowestCommonMultiple
  //
  // - Maybe try to get rid of some of the _ functions like _getDegree() etc. I don't think, we 
  //   will ever need them.


protected:

  std::vector<rsMonomial<T>> terms;
  TTol tol = TTol(0);

};


// ToDo: Move these implementations below into the class

/** Multiplies a coefficient and a sparse polynomial. */
template<class T, class TTol>
inline rsSparsePolynomial<T, TTol> operator*(const T& s, const rsSparsePolynomial<T, TTol>& p)
{
  rsSparsePolynomial<T, TTol> r;
  r._copyDataFrom(p);
  r.scale(s);
  return r;

  // ToDo: Write an operator that takes a monomial as left operand. It should scale r by the 
  // monomial's coeff as above and shift the powers of r by the monomial's power. Maybe the 
  // "copyDataFrom" function should already include the possible scaling and shifting. But then
  // we should call it copyScaledDataFrom and/or copyScaledAndShiftedDataFrom.
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::setupFromDenseCoeffs(
  const T* newCoeffs, int newNumTerms, TTol newTol)
{
  tol = newTol;
  terms.clear();
  terms.reserve(newNumTerms);
  for(int i = 0; i < newNumTerms; i++)
    if( !rsIsNegligible(newCoeffs[i], tol) ) 
      terms.emplace_back(rsMonomial<T>(newCoeffs[i], i));

  //canonicalize();  // Superfluous!
  // The result is actually ensured to be canonical already anyway. The dense coeffs are always in
  // the right order and we take care of not appending negligible coeffs.
}


/** Specializes rsMaxNorm() for rsSparsePolynomial. The max norm of a sparse polynomial is defined
as the maximum norm of all the coefficients. */
template<class T, class TTol>
auto rsMaxNorm(const rsSparsePolynomial<T, TTol>& p)
{
  auto max = rsMaxNorm(T(0));
  for(int i = 0; i < p.getNumTerms(); i++)
    max = rsMax(max, rsMaxNorm(p.getCoeff(i)));
  return max;
}


#endif