#ifndef RAPT_SPARSEPOLYNOMIAL_H
#define RAPT_SPARSEPOLYNOMIAL_H


//=================================================================================================

/** A class for representing (univariate) monomials, i.e. expressions of the form  c * x^p  for
some coefficient c and integer power (or exponent) p. Strictly speaking, we should require p to be
nonnegative, but we don't really enforce this here. So, it can later be used to represent more 
general terms of this form with possibly negative integer powers.

See: https://en.wikipedia.org/wiki/Monomial   */

template<class T>  // Maybe rename to TCoef 
class rsMonomial
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  explicit rsMonomial(T newCoeff = T(0), int newPower = 0) : coeff(newCoeff), power(newPower) { }
  // Marked as explicit because we want to avoid hidden automatic conversions from type T
  // Maybe pass newCoeff by const ref


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setup(T newCoeff, int newPower) { coeff = newCoeff; power = newPower; }

  void setCoeff(T newCoeff)   { coeff = newCoeff; }

  // ToDo: Pass newCoeff by const reference. Rationale: At some point we may want to use types T 
  // that are expensive to copy such as arbitrary precision floating point numbers.

  void setPower(int newPower) { power = newPower; }

  void negate() { coeff = -coeff; }

  void shiftPower(int amount) { power += amount; }

  void scaleCoeff(const T& scaler) { coeff *= scaler; }

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the coefficient c in the expression c * x^p. */
  T getCoeff() const { return coeff; }
  // Maybe return by const ref? ...not sure if this is a good idea, though.

  /** Returns the power (aka exponent) p in the expression c * x^p. */
  int getPower() const { return power; }


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Evaluates the expression c * x^p at the given x. The data type of the argument x may be 
  different from the coefficient type T with which the class is instantiated. Can be used, for 
  example, to evaluate monomials with real coefficients at complex arguments. */
  template<class TArg>
  TArg operator()(TArg x) const { return TArg(coeff) * rsPow(x, TArg(power)); }
  // Preliminary. We may want to use rsPowInt for integer exponents. That may be more efficient.
  // Here, we explicitly first convert the exponent to type TArg and then call 
  // rsPow(TArg x, TArg y). But currently rsPowInt is not suitably defined. It expects two 
  // (unsigned?) ints. We really need a function where the base is an arbitrary type and the 
  // exponent is an integer. I think, the fast exponentiation algorithm should work in these 
  // cases, too. It is important to not require x^n to be defined for n being of type T or TArg 
  // because the Type T or TArg may be a matrix type and we don't know how to raise a matrix to 
  // the power of another matrix - but we do know how to raise a matrix to an integer power. So,
  // to make the class as flexibly instantiatable as possible, we should require only integer
  // powers to be defined. Actually, it would be enough to require non-negative integer powers,
  // i.e. not require inversion/reciprocation to be defined. Or we could use an unsigned int type
  // for the power.

  /** Returns the negative of this monomial. */
  rsMonomial<T> operator-() const { return rsMonomial<T>(-getCoeff(), getPower()); }

  /** Multiplies two monomials. */
  rsMonomial<T> operator*(const rsMonomial<T>& q) const
  { return rsMonomial<T>(getCoeff() * q.getCoeff(), getPower() + q.getPower()); }
  // Needs tests. Maybe use coeff and power directly instead of using the getters. That makes the
  // code a bit shorter.

  /** Divides two monomials. */
  rsMonomial<T> operator/(const rsMonomial<T>& q) const
  { return rsMonomial<T>(getCoeff() / q.getCoeff(), getPower() - q.getPower()); }
  // Needs tests. 
  // If q.power > this->power, this will lead to a negative power in the result. Should we do 
  // something about this like triggering an rsAssert? And what if q.getCoeff() returns zero?
  // Maybe warn about this in the documentation and add assertions.

  // Addition and subtraction cannot be generally defined. These operations would only make sense
  // when both operands have the same power which is more an exceptional case rather than the rule
  // so it's better to not define these operators at all. At least, for the moment. But maybe we 
  // can define multiplication and division by constants for convenience. Maybe we should also 
  // define *= and /= in terms of setup(). It might avoid allocations for complicated types T such
  // as matrices.

protected:

  T   coeff = T(0);
  int power = T(0);

};


//=================================================================================================

/** A class for representing sparse polynomials, i.e. polynomials that have many zero coefficients.
Think of something like p(x) = 2*x^37 - 5*x^129 + 3*x^435. With a dense representation (as in e.g.
rsPolynomial), we would have to store 436 coefficients 433 of which would be zero. This is, of
course, inacceptably wasteful. And such polynomials actually do occur a lot in filter transfer 
functions that use delaylines instead of unit delays. Think of Schroeder allpasses, feedback delay 
networks (FDNs), etc. i.e. the building blocks of reverb algorithms. 

We represent such sparse polynomials basically as a std::vector of monomials which we call terms in
this context. We say that a sparse polynomial is in canonical representation if the powers of the 
terms are strictly increasing as function of array index (implying that no power appears more than 
once) and there are no terms with a coefficient of zero (up to some tolerance - i.e. the absolute 
values of all coeffs should be greater than the tolerance). If the array of terms is empty, we 
treat that as the canonical representation of the zero polynomial. This is different from how we 
represent the zero polynomial in the class rsPolynomial for dense polynomials. There, the coeff 
array of the zero polynomial has a single element whose value is zero. But because we want to avoid
(near) zero coeffs here, this representation would not be a good fit because then we would have to 
make a special rule for the zeroth coeff which may not even exist in certain nonzero sparse 
polynomials like p(x) = 3*x^100.

For many purposes, it is convenient to assume a canonical representation and many setters will 
maintain such a representation - but not all of them. Sometimes, one needs to - at least 
temporarily - violate such a canonical representation for performance reasons. Therefore, the API 
has two levels. A higher level that assumes and maintains canonical representations and a lower 
level that makes no such assumption and gives no such maintenance guarantee. The lower level member
functions are prefixed with an underscore _ to indicate at the call site that now some low-level 
stuff is going on and special care should be taken. When using the potentially decanonicalizing 
setup methods (prefixed by an underscore), there are 2 options: (1) You know exactly what you are 
doing and that this is in fact ok, i.e. doesn't actually decanonicalize. (2) You re-canonicalize 
after you have finished with your operations by calling e.g. canonicalize(). Option (2) is quite 
expensive, though. If you opt for option (1), it might be good practice to put something like
rsAssert(p._isCanonical()) after your operations on a polynomial p such that, if you inadvertently 
mess up the canonical representation of p, you will fail fast. As long as you don't use the 
underscore methods, you don't need to worry about this. But maybe in such cases you should worry 
about performance if you exclusively stick to the high-level API - especially if you call setup 
functions in loops over terms. */

template<class T, class TTol = rsEmptyType>
class rsSparsePolynomial
{

public:

  using SparsePoly = rsSparsePolynomial<T, TTol>;  // For convenience


  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Default constructor. Constructs an empty sparse polynomial. This represents, by definition,
  the zero polynomial. */
  rsSparsePolynomial() {}

  /** Creates a polynomial from an initializer list for the terms. */
  rsSparsePolynomial(std::initializer_list<rsMonomial<T>> initList, TTol tolerance) 
    : terms(initList), tol(tolerance) 
  { _canonicalize(); }

  /** Creates a polynomial from a std::vector that is interpreted as containing the coefficients 
  for a dense polynomial. */
  rsSparsePolynomial(const std::vector<T>& coefficients, TTol tolerance) 
  { tol = tolerance; setupFromDenseCoeffs(coefficients);  }
  // Maybe replace this by a static factory function fromDenseCoeffs()


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Reserves memory for the given number of terms. Can be called before calling functions like 
  addTerm() to pre-allocate the desired amount of memory beforehand when multiple terms are being
  added in a sequence. */
  void reserve(size_t numTerms) { terms.reserve(numTerms); }

  /** Clears the array of terms. */
  void clear() { terms.clear(); }

  /** Sets up the numerical tolerance that is used to determine if a coefficient should be 
  considered zero, i.e. within the numerical roundoff error. The new setting will immediately take
  effect. That is: If the polynomial currently contains any terms that fall below the new 
  threshold, they will be removed. */
  void setRoundoffTolerance(TTol newTolerance) { tol = newTolerance; _removeTermsWithZeroCoeff(); }

  /** Sets up the polynomial from a dense arrays of polynomial coeffs. When a coefficient in the 
  dense representation is zero, we will not create a term for that. The comparison to zero is to be
  understood as an inexact comparison with some tolerance that can be set up via 
  setRoundoffTolerance() to accomodate for floating point roundoff errors. */
  void setupFromDenseCoeffs(const std::vector<T>& newCoeffs)
  { setupFromDenseCoeffs(&newCoeffs[0], (int) newCoeffs.size()); }

  /** Like setupFromDenseCoeffs(const std::vector<T>&, ...) but for raw C-arrays. */
  void setupFromDenseCoeffs(const T* newCoeffs, int newNumTerms);
  
  /** Adds the term c * x^p with coeff c and power p to the polynomial. If a term with the same 
  power already exists, this will just shift its coefficient. If the cofficient happens to be zero 
  after shift (up to the roundoff tolerance), the term will be removed. */
  void addTerm(T coeff, int power);

  /** Adds the given monomial to the polynomial. */
  void addTerm(const rsMonomial<T>& newTerm)
  { addTerm(newTerm.getCoeff(), newTerm.getPower()); }

  /** Subtracts the given monomial from the polynomial. */
  void subtractTerm(const rsMonomial<T>& newTerm)
  { addTerm(-newTerm.getCoeff(), newTerm.getPower()); }

  /** Negates this polynomial, i.e. multiplies all coeffs by -1. */
  void negate() { for(auto& t : terms) t.negate(); }

  /** Makes the polynomial monic by dividing all coeffs by the leading coeff. A monic polynomial is
  a polynomial in which the leading coefficient is unity (aka one).*/
  void makeMonic() { _scaleCoeffs(T(1) / getLeadingCoeff()); }
  // ToDo: Maybe we should now again apply the tolerance threshold to the coeffs? But it would 
  // really be weird to remove coeffs in a makeMonic operation! Maybe we should treat the threshold
  // as relative (to the largest coeff) anyway? Then, the relative size of the coeffs wouldn't 
  // change due to scaling. A coeff that was above the relative threshold before scaling would 
  // still be above the relative threshold after scaling. ...at least if we ignore possible 
  // roundoff error changes.

  /** Shifts all powers by the given amount. If the amount is p, this corresponds to multiplying 
  the polynomial by a monomial factor with unit coefficient, i.e. by x^p. */
  void shiftPowers(int amount) { for(auto& t : terms) t.shiftPower(amount); }
  // ToDo: Maybe assert that the amount is <= the power of our smallest terms because otherwise,
  // we'll produce negative powers. It may at some point make sense to allow negative powers, 
  // though (for example, to represent truncated Laurent series), but at the moment, we assume to
  // deal with just normal polynomials.

  /** Multiplies this polynomial by the given other polynomial factor. Works in place and 
  re-allocates only when the capacity is too low (VERIFY!). */
  void multiplyBy(const SparsePoly& factor) { multiply(*this, factor, this); }

  /** Multiplies this polynomial by a dense polynomial represented by the given array of 
  coefficients. Works in place and re-allocates only when the capacity is too low. */
  void multiplyByDenseCoeffs(const T* coeffs, int numTerms);

  /** Adds the polynomial p scaled by the given scaler monomial to this polynomial. */
  void addScaled(const SparsePoly& p, const rsMonomial<T>& scaler);

  /** Adds the polynomial p scaled by the given scaling factor to this polynomial. */
  void addScaled(const SparsePoly& p, T scaler);


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the numerical tolerance that is used to determine if a coefficient should be 
  considered zero, i.e. within the numerical roundoff error. */
  TTol getRoundoffTolerance() const { return tol; }

  /** Returns the leading term in this polynomial, i.e. the monomial  a_n * x^p[n]  that has the 
  highest exponent p[n]. */
  rsMonomial<T> getLeadingTerm() const
  {
    rsAssert(_isCanonical());
    if(terms.empty())
      return rsMonomial<T>(T(0), 0);
    return terms[terms.size()-1];
  }
  // Maybe move to .cpp file

  /** Returns the leading coefficient of this polynomial, i.e. the coefficient in front of the 
  highest power of x. */
  T getLeadingCoeff() const { return getLeadingTerm().getCoeff(); }
  // Needs unit test. ..I think, it should be covered already because it's called by makeMonic 
  // which is called in the gcd algo - and gcd has a unit test. -> Verify this!

  /** Returns the degree of the polynomial, i.e. the exponent of the highest power of x that 
  occurs. */
  int getDegree() const { return getLeadingTerm().getPower(); }
  // Verify if this has unit tests that cover the edge cases (zero polynomial, constant polynomial,
  // polynomial with single term, ...)

  /** Returns true iff this polynomial is the zero polynomial. */
  bool isZero() const { rsAssert(_isCanonical()); return terms.empty(); }

  /** Returns true iff this polynomial is a constant polynomial, i.e. has a degree of zero. */
  bool isConstant() const { return getDegree() == 0; }
  // Needs test

  /** Returns true, iff the rhs polynomial equals this polynomial up to the given tolerance. */
  bool isCloseTo(const SparsePoly& rhs, TTol tol) const;

  /** Like isCloseTo above but without the tol parameter. It uses the maximum of our tol member 
  and the rhs's tol member. */
  bool isCloseTo(const SparsePoly& rhs) const { return isCloseTo(rhs, rsMax(tol, rhs.tol)); }

  /** Returns true, iff this polynomial is monic, i.e. if its leading coefficient is 1 (up to the
  numeric roundoff tolerance). */
  bool isMonic() const { return rsIsCloseTo(getLeadingCoeff(), T(1), tol); }

  /** Returns true, iff the greatest common divisor of this polynomial and the other one is just
  a constant. If this is the case, the two polynomials are said to be "mutually prime" or 
  "coprime". This is analoguous to the notion of coprimality of integer numbers. */
  bool isCoprimeTo(const SparsePoly& other) const
  { SparsePoly gcd = greatestCommonDivisor(*this, other); return gcd.isConstant(); }
  // Allocates! A non-allocating version would need temporary workspace parameters. We can't 
  // compute the gcd without temporary memory.
  // Needs tests. ToDo: Maybe the criterion is stricter and requires the constant to be one? But
  // no - I don't think so. Verify!

  /** Return true, iff the given index is valid, i.e. the object has a term with given index. */
  bool isValidIndex(int i) const { return i >= 0 && i < getNumTerms(); }

  /** Returns the number of terms in this polynomial. The i-th term is a monomial of the form 
  ci * x^pi with a coefficient ci and a power/exponent pi. */
  int getNumTerms() const { return (int) terms.size(); }

  /** Returns the term (i.e. the monomial) at the given index. */
  rsMonomial<T> getTerm(int index) const { rsAssert(isValidIndex(index)); return terms[index]; }

  /** Returns the coefficient of the term with given index. */
  T getCoeff(int index) const { rsAssert(isValidIndex(index)); return terms[index].getCoeff(); }
  // ToDo: Return the coeff as const reference! We intend this class to be potentially used with
  // large coeff types like matrices or arbitrary precision floats, so this optimization may make 
  // sense. Hmm...I tried it and it breaks a unit test. Figure out why that happens and document 
  // it!

  /** Returns the power of the term with given index. */
  int getPower(int index) const { rsAssert(isValidIndex(index)); return terms[index].getPower(); }


  //-----------------------------------------------------------------------------------------------
  /** \name Operators. The arithmetic operators all allocate. */

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
  // Maybe try to use std::accumulate() if possible. It might not be because of mismatch of types
  // of vector elements and output. Also, it may not use the potentially more efficient += 
  // operator. ...soo - maybe it's better to leave the code as is.

  /** Implements the unary minus operator. */
  SparsePoly operator-() const 
  { SparsePoly r(*this); r.negate(); return r; }

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
  /** \name Low level API. The static functions operate on pre-allocated output parameters passed 
  by pointer (to make it obvious at the call site that the parameter may be modified). Using these 
  functions with pre-allocated sparse polynomials may potentially avoid heap allocations which the 
  more convenient functions that return sparse polynomials do. This includes the +,-,*,/,% 
  operators. So, for real time code, these operators are actually forbidden and one has to resort 
  to the low level API. 
  
  The low level non-static member functions starting with an underscore are meant for low level 
  manipulations that may temporarily destroy the canonical representation. If you use these (for 
  example, for performance reasons), it's your own responsibility to maintain a canonical 
  representation or when you destroy it in some process, to restore it when your are finished, for
  example by calling _canonicalize(). */

  /** Computes the sum r = p + q of the polynomials p and q and stores the result in r. */
  static void add(const SparsePoly& p, const SparsePoly& q, SparsePoly* r);
  // ToDo: document whether or not it can be used in place.

  /** Computes the difference r = p - q of the polynomials p and q and stores the result in r. */
  static void subtract(const SparsePoly& p, const SparsePoly& q, SparsePoly* r);
  // ToDo: document whether or not it can be used in place.

  /** Computes the weighted sum r = wp * p + wq * q of the polynomials p and q and stores the 
  result in r. */
  static void weightedSum(const SparsePoly& p, T wp, const SparsePoly& q, T wq, SparsePoly* r);
  // ToDo: Document whether or not it can be used in place.

  /** Multiplies polynomials p and q and stores the result in r. It may be used in place, i.e. the
  result polynomial r can point to the memory location of the arguments p and/or q. */
  static void multiply(const SparsePoly& p, const SparsePoly& q, SparsePoly* r);
  // ToDo: document whether or not it can be used in place.
  // ...Done... ToDo: Verify, if we have unit tests for in place usage. Maybe document that, too.
  // But maybe not as part of the doxygen documentation

  /** Implements polynomial division with remainder. ...TBC... */
  static void divide(const SparsePoly& numerator, const SparsePoly& denominator,
    SparsePoly* quotient, SparsePoly* remainder);

  /** Computes the greatest common divisor of the polynomials p and q. The "monic" parameter 
  selects if the returned gcd should be normalized to be monic (that is: have leading coeff 1).
  The gcd of polynomials is unique only up to a constant scale factor, so it may make sense to 
  make it uniquely defined by requiring it to be monic. If "monic" is false, the returned gcd may 
  be scaled by some arbitrary scale factor which I'm not sure about whether or not it could be 
  useful to retain. That's why I make the normalization to a monic polynomial optional. After all, 
  this normalization would throw away a piece of information. See comments in the .cpp file for an 
  idea for what the leading coeff could mean. */
  static SparsePoly greatestCommonDivisor(
    const SparsePoly& p, const SparsePoly& q, bool monic = true);

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


  /** Turns the representation of the polynomial into a canonical one. A canonical representation 
  has the following properties: (1) The powers are strictly increasing as function of array index.
  (2) No power appears more than once. (3) No zero coefficients appear. We achieve this by 
  first sorting the terms, then combining multiple terms with equal exponents into single terms 
  and finally deleting all terms that have a coefficient zero (up to the roundoff tolerance). */
  void _canonicalize();

  /** Removes all the terms that have a coefficient of zero (up to the roundofff tolerance). */
  void _removeTermsWithZeroCoeff();

  /** Appends a term with given coeff and power to the end of our terms array. This may 
  decanonicalize the representation by appending a term of a power lower than the current degree
  and/or by duplicating one of the existing exponents and/or by having a zero coefficient. So, if 
  you use this function, you need to either ensure that none of these things happen or else call 
  canonicalize() at some point after your manipulations. */
  void _appendTerm(T coeff, int power) { terms.emplace_back(rsMonomial<T>(coeff, power)); } 

  /** Sets the number of terms. If the new number is less than the current number, it will just 
  cut off terms from the end. If the new number is greater than the current number, it will just
  extend our vector of terms and the added terms at the end are uninitialized, i.e. may contain 
  garbage. This function should only be used if you intend to set up the new terms via e.g. 
  _setTerm() after calling _setNumTerms(). So, it's a function that needs a lot of care to be used
  properly. */
  void _setNumTerms(int newNumTerms) { terms.resize(newNumTerms); }

  /** Directly sets the coefficient and power of the term with given index. This may 
  decanonicalize the representation in all sorts of ways: by destroying the strict order, by 
  destroying the "powers appear at most once" property, by setting a coeff to zero. */
  void _setTerm(int index, T coeff, int power) 
  { rsAssert(isValidIndex(index));  terms[index].setup(coeff, power); }

  /** Directly sets the power of the term with given index. It may decanonicalize by destroying 
  the "all powers appear only once" property and/or the "powers are in strict ascending order"
  property. */
  void _setPower(int index, int newPower)
  { rsAssert(isValidIndex(index)); terms[index].setPower(newPower); }

  /** Directly sets the coefficient of the term with given index. It may destroy the canonical
  representation by setting a coeff to zero. */
  void _setCoeff(int index, T newCoeff)
  { rsAssert(isValidIndex(index)); terms[index].setCoeff(newCoeff); }

  /** Scales the coefficient with the given index by the given scaler. It may destroy the 
  canonical representation by setting a coeff to zero. */
  void _scaleCoeff(int index, T scaler) { _setCoeff(index, scaler * getCoeff(index)); }

  /** Scales all coefficients by the given scaler. It may destroy the canonical representation 
  by setting the coeffs to zero. */
  void _scaleCoeffs(T scaler) { for(auto& t : terms) t.scaleCoeff(scaler); }

  /** Shifts the coefficient with the given index by the given amount, i.e. adds the given amount 
  to the coeff. It may decanonicalize the representation by leading to a zero coeff. */
  void _shiftCoeff(int index, T amount) { _setCoeff(index, amount + getCoeff(index)); }

  /** Shifts the power at the given index by the given amount. It may decanonicalize the 
  representation by introducing two terms with equal power and/or destroying the increasing order
  of terms. */
  void _shiftPower(int index, int amount) { _setPower(index, amount + getPower(index)); }

  /** Multiplies this polynomial by the given monomial factor. This results in all coeffs being 
  multiplied by the coeff of the monomial and all powers being increased by the power of the 
  monomial. */
  void _multiplyBy(const rsMonomial<T>& factor)
  { _scaleCoeffs(factor.getCoeff()); shiftPowers(factor.getPower()); }
  // How could it possibly decanonicalize the representation? I think, it can't. Maybe it's a safe
  // method and should be moved up into the high level API and lose its underscore. Ah - wait!
  // With the current implementation, it can indeed decanonicalize. However - it would be easy to 
  // implement it in a safe way. Just check, if factor.getCoeff() is zero. If it is, set "this"
  // canoncially to zero. If it isn't, the function will not decanonicalize. ...I think....

  /** Divides this polynomial by the given monomial factor. This results in all coeffs being 
  multiplied by reciprocal of the coeff of the monomial and all powers being decreased by the 
  power of the monomial. If you are not careful, this may result in polynomials that contain 
  negative powers. Currently, this would be seen as a bug (and may trigger an rsAssert at some 
  point later) - but if needed, we can easily lift this restriction at some point. */
  void _divideBy(const rsMonomial<T>& divisor)
  { _scaleCoeffs(T(1) / divisor.getCoeff()); shiftPowers(-divisor.getPower()); }
  // Needs tests. 
  // Maybe assert that this->getMinPower() >= divisor.getPower() to avoid producing negative 
  // powers. Also: what if divisor.coeff is zero (division by zero error) or infinite (makes all
  // our coeffs zero and thereby destroys canonical representation)

  /** Reverses the array of terms. It may appear to be a weird thing to do on polynomials but this 
  operation is needed when transforming minimum phase filters into maximum phase ones (or vice 
  versa) and when producing allpass filters from allpole filters. This destroys the canonical 
  order from lower to higher powers of x. */
  void _reverse() { rsReverse(terms); }

  /** Checks if this sparse polynomial is in canonical representation. A representation is 
  canonical if it has no zero coefficients (up to the roundoff tolerance) and if the powers are 
  strictly increasing (as function of array index) and if no power occurs more than once. The empty
  polynomial is also accepted as a canonical representation. It represents the zero polynomial. */
  bool _isCanonical() const;

  /** Returns true iff the powers of our terms are strictly increasing as function of array index.
  This strict monotonicity also entails uniqueness of the powers. That means that this function 
  serves two purposes at the same time: Making sure that the terms are sorted by power and that no
  power appears more than once. These are two of the requirements for a canonical 
  representation. */
  bool _areTermsStrictlySorted() const;
  // Needs test

  /** Returns true iff any of our terms array has a coefficient of zero (up to the roundoff 
  tolerance). In a canonical representation, this is forbidden. */
  bool _hasZeroCoeffs() const;
  // Needs test

  /** Returns true iff any of our terms has a negative power. We currently consider this as not
  allowed, i.e. a bug - but this restriction can be lifted later, if needed. */
  bool _hasNegativePowers() const;
  // Needs test


protected:

  std::vector<rsMonomial<T>> terms;  // Terms of the form c_i * x^p[i]
  TTol tol = TTol(0);                // Roundoff error tolerance (relevant for e.g. T = float)

};


/** Multiplies a coefficient and a sparse polynomial. */
template<class T, class TTol>
inline rsSparsePolynomial<T, TTol> operator*(
  const T& s, const rsSparsePolynomial<T, TTol>& p)
{
  rsSparsePolynomial<T, TTol> r(p);
  r._scaleCoeffs(s);
  return r;
}

/** Multiplies a monomial and a sparse polynomial. */
template<class T, class TTol>
inline rsSparsePolynomial<T, TTol> operator*(
  const rsMonomial<T>& mon, const rsSparsePolynomial<T, TTol>& p)
{
  rsSparsePolynomial<T, TTol> r(p);
  r._scaleCoeffs(mon.getCoeff());
  r.shiftPowers(mon.getPower());
  return r;
}
// Needs tests!


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