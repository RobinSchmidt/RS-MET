#ifndef RAPT_POLYNOMIAL_H
#define RAPT_POLYNOMIAL_H

/** A class for representing polynomials and doing computations with them. */

// todo: 

// -implement greatest common divisor algorithm (maybe the one for the integers can be
//  used as is?)
// -implement missing operators:
//  -operator() that takes a polynomial and returns another polynomial as result (implements 
//   nesting)
//  -arithmetic operators that take a number as second (left or right) argument
//  -maybe we could meaningfully define <,<=, ...? look at how python or other scientific libraries
//   handle that - in my own python polynoamial class, i'm taking the asymptotic behavior


template<class T>
class rsPolynomial
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction */

  /** Creates a polynomial of given degree. Allocates memory for the coefficients and optionally
  intializes them with zeros. */
  rsPolynomial(int degree = 0, bool initWithZeros = true);

  /** Creates a polynomial from a std::vector of coefficients. */
  //rsPolynomial(const std::vector<T>& coefficients) : coeffs(coefficients) {}
  rsPolynomial(const std::vector<T>& coefficients) 
  {
    setCoeffs(coefficients);
  }

  /** Promotes a number to a 0th degree polynomial. */
  rsPolynomial(const T& number)
  { coeffs.resize(1); coeffs[0] = number; }


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setCoeffs(const T* newCoeffs, int newDegree);

  void setCoeffs(const std::vector<T>& newCoeffs)
  { setCoeffs(newCoeffs.data(), (int)newCoeffs.size()-1); }

  void setRoots(const T* newRoots, int numRoots, T scaler = T(1));

  void setRoots(const std::vector<T>& newRoots, T scaler = T(1))
  { setRoots(newRoots.data(), (int)newRoots.size(), scaler); }

  void truncateTrailingZeros(const T& threshold = T(0));

  void negate()
  { rsArrayTools::negate(&coeffs[0], &coeffs[0], (int) coeffs.size()); }
  // the conversion to int may be avoided



  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the maximum order that this poloynomial may have which is the length of the
  coefficient array minus one. When there are trailing zero coefficients, the actual degree of
  the polynomial is lower. */
  //int getMaxOrder() const { return (int)coeffs.size()-1; }
  // deprecate this 

  /** Returns the degree of the polynomial, defined as... */
  int getDegree() const { return (int)coeffs.size()-1; }
  // should take into account traling zeros ..or maybe have a boolean flag
  // "takeZeroCoeffsIntoAccount" which defaults to false...or maybe it shouldn't have any default
  // value - client code must be explicit...or maybe have functions getAllocatedDegree, 
  // getActualDegree(tolerance)...or getDegree has an optional parameter for the tolerance 
  // defaulting to 0

  /** Returns the leading coefficient, i.e. the coefficient that multiplies the highest power of 
  x. */
  T getLeadingCoeff() const { return rsLast(coeffs); }
  // what if we have trailing zeros in the coeff array?

  /** Returns a pointer to our coefficient array - breaks encapsulation - use with care! */
  T* getCoeffPointer() { return &coeffs[0]; }
  // nope! when we really need low-level access to the coeff-array, we declare the 
  // functions/classes that need it as friends. comment this out later

  const T* getCoeffPointerConst() const { return &coeffs[0]; }

  /** Returns true, iff this polynomial is monic, i.e. the coefficient for the highest power (the
  leading coefficient) is unity. Monic polynomials are important because they arise when 
  multiplying out the product form. */
  bool isMonic() const { return getLeadingCoeff() == T(1); }
  // what if we have trailing zeros in the coeff array? should we have a tolerance?

  /** Evaluates the first derivative of this polynomial at the given x. */
  T derivativeAt(const T& x) 
  { return evaluateDerivative(x, &coeffs[0], getDegree()); }

  /** Evaluates the order-th derivative of this polynomial at the given x. Works also for the 0th 
  derivative, which is the function value itself. ...but the order must be non-negative. */
  T derivativeAt(const T& x, int order) 
  { return evaluateDerivative(x, &coeffs[0], getDegree(), order); }
  // todo: maybe make it also work for negative orders (in which case the antiderivative of 
  // given order will be evaluated (setting integration constants to zero))


  //T definiteIntegral(const T& lowerLimit, const T& upperLimit);



  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Adds two polynomials. */
  rsPolynomial<T> operator+(const rsPolynomial<T>& q) const {
    rsPolynomial<T> r(rsMax(getDegree(), q.getDegree()), false);
    weightedSum(coeffs.data(), getDegree(), T(1),
      q.coeffs.data(), q.getDegree(), T(1),
      r.coeffs.data());
    return r;
  }

  /** Subtracts two polynomials. */
  rsPolynomial<T> operator-(const rsPolynomial<T>& q) const {
    rsPolynomial<T> r(rsMax(getDegree(), q.getDegree()), false);
    weightedSum(coeffs.data(), getDegree(), T(+1),
      q.coeffs.data(), q.getDegree(), T(-1),
      r.coeffs.data());
    return r;
  }

  /** Multiplies two polynomials. */
  rsPolynomial<T> operator*(const rsPolynomial<T>& q) const {
    //rsPolynomial<T> r(getDegree() + q.getDegree() + 1, false);
    rsPolynomial<T> r(getDegree() + q.getDegree(), false);
    multiply(coeffs.data(), getDegree(), q.coeffs.data(), q.getDegree(), r.coeffs.data());
    return r;
  }

  /** Divides this polynomial by the given divisor polynomial d and returns the quotient. */
  rsPolynomial<T> operator/(const rsPolynomial<T>& d) const {
    rsPolynomial<T> q(getDegree(), false);  // quotient
    rsPolynomial<T> r(getDegree(), false);  // remainder
    divide(coeffs.data(), getDegree(), 
      d.coeffs.data(), d.getDegree(), q.coeffs.data(), r.coeffs.data());
    return q;
  }

  /** Divides this polynomial by the given divisor polynomial d and returns the remainder. */
  rsPolynomial<T> operator%(const rsPolynomial<T>& d) const {
    rsPolynomial<T> q(getDegree(), false);  // quotient
    rsPolynomial<T> r(getDegree(), false);  // remainder
    divide(coeffs.data(), getDegree(),
      d.coeffs.data(), d.getDegree(), q.coeffs.data(), r.coeffs.data());
    return r;
  }

  // maybe move implementations into cpp file, implement +=, *=, etc.
  // +,-,*,/ for scalar second arguments (left and right), unary -

  bool operator==(const rsPolynomial<T>& p) const { return coeffs == p.coeffs; }

  bool operator!=(const rsPolynomial<T>& p) const { return coeffs != p.coeffs; }

  rsPolynomial<T>& operator=(const rsPolynomial<T>& p)
  {
    setCoeffs(p.coeffs);
    return *this;
  }


  /** Multiplies a polynomial and a number. */
  rsPolynomial<T> operator*(const T& k) const 
  {
    return rsPolynomial<T>(k*coeffs);
  }

  /** Raises a polynomial to a (non-negative) integer power. ...Hmm...it's perhaps no good idea to 
  implement the operator that way because it has lower precedence than *, / in C++ - this can be 
  confusing - maybe write a "pow" function instead. */
  rsPolynomial<T> operator^(int k) const 
  {
    // optimize: allocate once and repeatedly convolve the coeff-array of r with the coeff array 
    // of this (in place)
    rsPolynomial<T> r(std::vector<T>({T(1)}));
    for(int i = 1; i <= k; i++)
      r = r * (*this);
    return r;
  }

  /** Allow the rsPolynomial object to be converted to a (const) std::vector of coefficients. */
  operator const std::vector<T>&() const { return coeffs; }
  // todo: maybe have a non-const version?


  // how to deal with the trailing zeros in quotient and/or remainder? should we cut them
  // off...if so, what should be the numerical threshold? maybe there should be a member function
  // removeTrailingZeros that client code must explicitly call
  // maybe we should also cutoff in add/subtract, if the trailing coeffs happen to add/subtract to
  // zero? what about multiplication (i.e. convolution) - can trailing zeros happen there, too - or
  // is that impossible? ...maybe try to get a zero with two very low (1st) order polynomials
  // no - it can't happen - the highest power coeff is always the product of the two highest power
  // coeffs of the input polynomials - and if they are assumed to be nonzero, so is their product


  /** Evaluates the polynomial at the given input x. */
  T operator()(T x) const { return evaluate(x, &coeffs[0], getDegree()); }
  // todo: have an overloaded operator() that takes a polynomial as input and returns another 
  // polynomial -> implement nesting/composition

  /** Read and write access to i-th coefficient (breaks encapsulation - use with care). */
  T& operator[](int i) { return coeffs[i]; }

  //===============================================================================================
  /** \name Computations on raw coefficient arrays */

  //-----------------------------------------------------------------------------------------------
  /** \name Evaluation */

  /** Evaluates the polynomial defined by the array of coefficients "a" at argument "x".  The array
  of coefficients must be of length degree+1 and is interpreted as follows: a[0] is taken to be the
  constant term, a[1] is the multiplier for x^1, a[2] the multiplier for x^2 and so on until
  a[degree] which is the multiplier for a^degree. */
  static T evaluate(const T& x, const T *a, int degree);

  /** Evaluates the polynomial defined by the array of roots "r" at argument "x". If infinite roots
  are encountered, they are skipped - this is consistent with what we need when evaluating filter
  transfer functions that have zeros at infinity. */
  static std::complex<T> evaluateFromRoots(const std::complex<T>& x,
    const std::complex<T>* roots, int numRoots);

  /** Like evaluateFromRoots but leaves one root out in the evaluation. This is equivalent to 
  evaluating the polynomial and divide the result by the linear factor corresponding to the left
  out root (at the given x), i.e. evaluating g(x) = f(x)/(x-r_i) where r_i is the i-th root that
  has been left out. */
  static std::complex<T> evaluateFromRootsOneLeftOut(const std::complex<T>& x,
    const std::complex<T>* roots, int numRoots, int leaveOutIndex);

  /** Evaluates the first derivative of the polynomial a at the given x.  */
  static T evaluateDerivative(const T& x, const T *a, int degree);
  //{ T y, yd; evaluateWithDerivative(x, a, degree, &y, &yd); return yd; }
  // todo: implement it similar to evaluateDerivative that allows the order to be passed - the only
  // difference is that the product is replaced by a single term

  /** Evaluates the order-th derivative at the given x. */
  static T evaluateDerivative(const T& x, const T *a, int degree, int order);
  //static inline T evaluateDerivative(const T& x, const T *a, int degree, int order)
  //{ T y[32]; evaluateWithDerivatives(x, a, degree, y, order); return y[order]; }
  // todo: use algorithm in polyDerivative in MathExperiments.cpp - maybe compare performance to
  // the current one

  /** Evaluates the polynomial defined by the array of coefficients 'a' and its first derivative at
  argument 'x'. The value of the polynomial will be stored in y and the value of the derivative
  will be stored in yd. @see: evaluate() */
  static void evaluateWithDerivative(const T& x, const T *a, int degree, T *y, T *yd);

  /** Evaluates the polynomial defined by the array of coefficients 'a' and a given number of
  derivatives at argument 'x'. The value of the polynomial will be stored in results[0], the 1st
  derivative in results[1] and so on. numDerivatives should be <= 31.
  @see: evaluatePolynomialAt() */
  static void evaluateWithDerivatives(const T& x, const T *a, int degree, T *results, 
    int numDerivatives);

  /** Evaluates a complex polynomial with coeffs "a" and its first and second derivative at the 
  input "z", stores the results in P[0],P[1],P[2] and returns an error estimate for the evaluated
  P[0] (verify this). This is used in the Laguerre root-finding algorithm. */
  static T evaluateWithTwoDerivativesAndError(const std::complex<T>* a, int degree,
    std::complex<T> z, std::complex<T>* P);
  // rename "P" to "y" ...or "w2 as is common in complex functions

  /** Evaluates the cubic polynomial a[0] + a[1]*x + a[2]*x^2 + a[3]*x^3 at the given x. */
  static inline T evaluateCubic(const T& x, const T* a)
  {
    return a[0] + (a[1] + (a[2] + a[3]*x)*x)*x;
    //T x2 = x*x; return a[0] + a[1]*x + a[2]*x2 + a[3]*x*x2;  // alternative implementation
  }

  /** Evaluates the cubic polynomial a + b*x + c*x^2 + d*x^3 at the given x. */
  static inline T evaluateCubic(const T x, const T& a, const T& b, const T& c, const T& d)
  {
    return a + (b + (c + d*x)*x)*x;
  }
  // make consistent with rootCubic

  /** Evaluates the Hermite polynomial of given degree. This is the physicist's version of the
  Hermite polynomials with leading coefficient 2^n. The probabilist would use those with leading 
  coeff 1. */
  static T evaluateHermite(const T& x, int degree);

  // todo: evaluateDerivative, evaluateIntegral (or AntiDerivative)

  /** Given a coefficient array p of length maxDegree+1, this function returns the actual degree
  of the polynomial - meaning that trailing zeros in the array don't count. So it's the degree 
  that takes only into account the non-zero coefficients. */
  //static int actualDegree(T* p, int maxDegree, T tol = T(0));
  // maybe move to an "Inquiry" section ...or "Misc"

  //-----------------------------------------------------------------------------------------------
  /** \name Arithmetic */

  /** Forms a weighted sum of two polynomials p(x) and q(x) with weights wp and wq respectively and
  stores the coeffficients of the resulting polynomial r(x) in the r-array. The polynomials p(x)
  and q(x) do not need to be of the same degree and the resulting polynomial will have an degree of
  max(pN, qN). However, in this weighted sum, it may happen that the higher order terms sum to zero
  in which case the *actual* degree of the resulting polynomial may be less than that. In an
  extreme case, you could even end up with the polynomial that is identically zero. If it's
  important to know the *actual* degree rather than the length of the coeff-array, the caller needs
  to check for (non)-zero-ness of the coeff values in the result. */
  static void weightedSum(const T* p, int pN, const T& wp, const T* q, int qN, const T& wq, T* r);

  /** Adds polynomial q(x) from polynomial p(x) and stores the coeffficients of the resulting
  polynomial in r(x) which is of degree max(pN, qN). */
  static void add(const T* p, int pN, const T* q, int qN, T* r)
  { weightedSum(p, pN, T(1), q, qN, T(1), r); }

  /** Subtracts polynomial q(x) from polynomial p(x) and stores the coeffficients of the resulting
  polynomial in r(x) which is of degree max(pN, qN). */
  static void subtract(const T* p, int pN, const T* q, int qN, T* r) 
  { weightedSum(p, pN, T(1), q, qN, T(-1), r); }

  /** Multiplies the polynomials represented by the coefficient vectors 'a' and 'b' and stores the
  resulting coefficients in 'result'. The resulting polynomial will be or order aDegree+bDegree and
  the coefficient vector should have allocated space for
  (aDegree+1)+(bDegree+1)-1 = aDegree+bDegree+1 = aLength+bLength-1 elements.
  Can work in place, i.e. result may point to the same array as a and/or b.   */
  static void multiply(const T *a, int aDegree, const T *b, int bDegree, T *result)
  {
    rsArrayTools::convolve(a, aDegree+1, b, bDegree+1, result);
  }

  /** Divides the polynomials represented by the coefficient arrays 'dividend' and 'divisor' and
  stores the resulting coefficients for the quotient and remainder in the respective output arrays.
  ?? The resulting quotient polynomial will be of degree dividendDegree-divisorDegree and the
  remainder polynom will be at most of degree... ??
  However, the output arrays must have the same length as the dividend, where
  remainder[divisorDegree...dividendDegree] and
  quotient[dividendDegree-divisorDegree+1...dividendDegree] will be filled with zeros.
  ...\todo check this comment */
  static void divide(const T* dividend, int dividendDegree, const T* divisor, int divisorDegree,
    T* quotient, T* remainder);
  // todo: complete this comment ...figure out, what the degrees of the quotient and remainder are
  // ...supposedly, they can be inferred by client code by checking from which index in the array
  // there are only zeros?

  /** Divides the dividend by the monomial factor (x-x0) and stores the result in the same array
  again. The remainder (which is just a numerical value) will be stored in 'remainder'. This
  function is useful for splitting off a linear factor, if x0 is a root of the dividend, in which
  case the remainder will come out as zero (check this). */
  template<class S>
  static void divideByMonomialInPlace(S* dividendAndResult, int dividendDegree, S x0, S* remainder);
  // maybe make a version that can store the result in a different array - make it so that the
  // result array may or may not point to the same array as the dividend
  // shouldn't x0 be const S& ? ...comment, what S is and why we don't use T ..it has to do with
  // real vs complex


  static void greatestCommonDivisor(const T* p, int pDegree, const T* q, int qDegree, 
    T* gcd, int* gcdDegree, T tolerance = T(0));
  // under construction
  // maybe implement extended GCD


  /** Creates an array of arrays with polynomial cofficients that represent the polynomial with
  coefficients a[] raised to successive powers up to and including "highestPower". aPowers[0]
  will have a single entry equal to unity (representing a[]^0), aPowers[1] will contain a copy
  of a[] itself (representing a[]^1), aPowers[2] will contain a[] convolved with itself
  (representing a[]^2), aPowers[3] will contain a[]^2 convolved with a[] and so on. Thus, we
  repeatedly convolve the result of the previous iteration with the a[] array. The function fills
  up the arrays only up to the point where the array of polynomial coefficients actually ends, it
  doesn't fill additional zeros. That means, you may pass an array of arrays in aPowers, where
  each sub-array has exactly the length that is strictly required. If you allocate more memory (for
  convenience or whatever), make sure to initialize the sub-arrays with zeros. */
  static void powers(const T* a, int N, T** aPowers, int highestPower);

  /** Let A(x) and B(x) be polynomials represented by their coefficient arrays a[] and b[]
  respectively. This function creates the coefficients of a polynomial C(x), represented by the
  coefficient array c[], that results from composing the polynomials A(x) and B(x), that is: first
  the polynomial A(x) is applied to the value x, and then the polynomial B(x) is applied to the
  result of the first polynomial, such that C(x) = B(A(x)). This nesting or composition of two
  polynomials can itself be seen as a polynomial in its own right. This resulting polynomial has
  a degree of cN = aN*bN, where aN and bN are the degrees of the a[] and b[] polynomials,
  respectively, so the caller has to make sure that the c[] array has at least a length of
  aN*bN+1. */
  static void compose(const T* a, int aN, const T* b, int bN, T* c);
  // allocates heap memory

  /** Given an array of polynomial coefficients "a" such that
  p(x) = a[0]*x^0 + a[1]*x^1 + ... + a[N]*x^N, this function returns (in "am") the coefficients for
  a polynomial q(x) such that q(x) = p(-x). This amounts to sign-inverting all coefficients which
  multiply odd powers of x. */
  static void coeffsForNegativeArgument(const T *a, T *am, int N);
  // rename to negateArgument

  /** Given an array of polynomial coefficients "a" such that
  p(x) = a[0]*x^0 + a[1]*x^1 + ... + a[N]*x^N, this function returns (in "aShifted") the coefficients
  for a polynomial q(x) such that q(x) = p(x-x0). */
  static void coeffsForShiftedArgument(const T *a, T *aShifted, int N, T x0);
  // allocates heap memory
  // rename to shiftArgument, maybe move into "Conversions" section


  //-----------------------------------------------------------------------------------------------
  /** \name Calculus */

  /** Finds the coefficients of the derivative of the N-th degree polynomial with coeffs in "a" and
  stores them in "ad". The degree of the polynomial represented by the coeffs in "ad" will be
  N-1. The "ad" array may point to the same array as the "a" array, so you can use the same array
  for input and output. */
  static void derivative(const T *a, T *ad, int N);

  /** Finds the coefficients of the indefinite integral of the N-th degree polynomial with coeffs
  in "a" and stores them in "ai". The degree of the polynomial represented by the coeffs in "ai"
  will be N+1. The constant term in the ai[] polynomial is the arbitrary integration constant
  which may be passed in as parameter "c" - this parameter is optional, it defaults to zero.  */
  static void integral(const T *a, T *ai, int N, T c = T(0));
  // maybe rename to antiderivative?

  /** Computes the definite integral of the polynomial "p" where the lower integration limit is
  given by the polynomial "a" and the upper limit is given by the polynomial "b". "p", "a", "b"
  are assumed to be of degrees "pN", "aN" and "bN" respectively and the result will be stored in
  as polynomial "q" which will be of degree pN*max(aN, bN). */
  static void integrateWithPolynomialLimits(const T *p, int pN, const T *a, int aN, const T *b, 
    int bN, T *q);
  // allocates heap memory

  /** Given a polynomial p(x) via its coefficient array a, this function computes the coefficients
  of a polynomial q(x) = p(x+h) - p(x) if direction = 1, or q(x) = p(x) - p(x-h) if direction = -1.
  This represents the first forward or backward difference polynomial in finite difference
  calculus and is analog to the derivative in infinitesimal calculus. As in infinitesimal calculus,
  the resulting polynomial is one degree less than p(x) itself. Scaled by 1/h, this can be seen as 
  an approximation to the derivative using stepsize h.
  \todo provide a finite central difference q(x) = p(x+h/2) - p(x-h/2) when direction = 0
   ...could this be just the average between forward and backward difference? ...research! */
  static void finiteDifference(const T *a, T *ad, int N, int direction = 1, T h = 1);
  // allocates heap memory


  //-----------------------------------------------------------------------------------------------
  /** \name Roots */

  /** Finds all complex roots of a polynomial by Laguerre's method and returns them in "roots". */
  static void roots(const std::complex<T>* a, int degree, std::complex<T>* roots);
  // allocates heap memory

  static void roots(const T *a, int degree, std::complex<T> *roots);
  // allocates heap memory

  /** Converges to a complex root of a polynomial by means of Laguerre's method using the
  "initialGuess" as first estimate. */
  static std::complex<T> convergeToRootViaLaguerre(const std::complex<T> *a, int degree,
    std::complex<T> initialGuess = std::complex<T>(0.0, 0.0));
  // allocates heap memory

    /** Computes the root of the linear equation: \f[ a x + b = 0 \f] which is simply given by
  \f[ x_0 = -\frac{b}{a} \f] */
  static T rootLinear(const T& a, const T& b);
  // rename inputs to a0,a1 (change their order)

  /** Computes the two roots of the quadratic equation: \f[ a x^2 + b x + c = 0 \f] which are
  given by: \f[ x_{1,2} = \frac{-b \pm \sqrt{b^2-4ac}}{2a} \f] and stores the result in two-element
  array which is returned. When the qudratic is degenerate (i.e, a == 0), it will fall back to the
  rootsLinear() function, and return a one-element array.  */
  static std::vector<std::complex<T>> rootsQuadratic(const T& a, const T& b, const T& c);
  // rename inputs to a0,a1,a2 (change their order) ..or maybe deprecate this function

    /** Computes the two roots of the quadratic equation: \f[ a_0 + a_1 x + a_2 x^2 = 0 \f] and
  stores them in r1, r2. When the equation has two distinct real roots, they will be returned in
  ascending order, i.e. r1 < r2. In case of a (real) double root, we'll have r1 == r2 and when the
  roots of the equation are actually complex, the outputs will also be equal and contain the real
  part of the complex conjugate pair. */
  template<class R>
  static void rootsQuadraticReal(const R& a0, const R& a1, const R& a2, R* r1, R* r2);

  static void rootsQuadraticComplex(
    const std::complex<T>& a0, const std::complex<T>& a1, const std::complex<T>& a2,
    std::complex<T>* x1, std::complex<T>* x2);
    // todo: make optimized version for real coefficients (but complex outputs)

  /** Computes the three roots of the cubic equation: \f[ a x^3 + b x^2 + c x + d = 0 \f] and
  stores the result in the three-element array which is returned. When the cubic is degenerate
  (i.e, a == 0), it will fall back to the getRootsOfQuadraticEquation() function, and return a
  two-element array (or a one-element array, when b is also zero). */
  static std::vector<std::complex<T>> rootsCubic(const T& a, const T& b, const T& c, const T& d);
  // todo: make the order of the arguments consistent with evaluateCubic - but careful - this will
  // break client code! ...rename parameters to a0,a1,a2,a3 to make it clear, how it's meant

  /** Discriminant of cubic polynomial \f[ a_0 + a_1 x + a_2 x^2 + a_3 x^3 = 0 \f].
  D > 0: 3 distinct real roots, D == 0: 3 real roots, 2 or 3 of which may coincide,
  D < 0: 1 real root and 2 complex conjugate roots */
  static T cubicDiscriminant(const T& a0, const T& a1, const T& a2, const T& a3);
  // rename to discriminantCubic

  // todo: write function quadraticDiscriminant

  /** under construction - does not yet work */
  static void rootsCubicComplex(
    std::complex<T> a0, std::complex<T> a1, 
    std::complex<T> a2, std::complex<T> a3, 
    std::complex<T>* r1, std::complex<T>* r2, std::complex<T>* r3);


  // implement rootsQuadraticReal, rootsQuadraticComplex, rootsCubicReal, rootsCubicComplex
  // rootsQuarticReal, rootsQuarticComplex
  //

  /** Iteratively improves an initial estimate for the root of the cubic equation:
  \f[ a x^3 + b x^2 + c x + d = 0               \f]
  by means of the Newton-Raphson iteration:
  \f[ x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}    \f]
  where f and f' are calcutated as:
  \f[ f(x) = ax^3+bx^2+cx+d =  ((ax+b)x+c)x+d   \f]
  \f[ f(x) = 3ax^2+2bx+c    =  (3ax+2b)x+c      \f]
  the arguments min and max give upper and lower bounds for the root (which will be returned in
  cases where the iteration diverges, which the caller should avoid in the first place) and
  maxIterations gives the maximum number of iteration steps. */
  static T cubicRootNear(T x, const T& a, const T& b, const T& c, const T& d, const T& min, 
    const T& max, int maxIterations = 10);
  // todo: rename to rootCubicNear, change order of variables, maybe use bisection, if 
  // newton-iteration diverges

  /** Iteratively improves an initial estimate for the root of the polynomial equation:
  \f[ a[order] x^order + ... + a[1] x + a[0] = 0   \f]
  by means of the Newton-Raphson iteration:
  \f[ x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}       \f]
  the arguments min and max give upper and lower bounds for the root (which will be returned in
  cases where the iteration diverges, which you should avoid in the first place) and maxIterations
  gives the maximum number of iteration steps. */
  static T rootNear(T x, const T* a, int order, const T& min, const T& max, 
    int maxIterations = 32);

  /** Same as above but accepts real coefficients. */
  //void findPolynomialRootsInternal(T *a, int order, Complex *roots, bool polish = true);

  //void findPolynomialRootsNew(Complex *a, int order, Complex *roots);


  //-----------------------------------------------------------------------------------------------
  /** \name Conversions */
  // changes for the representation of the polynomial

  /** Given expansion coefficients a[k] of an arbitrary polynomial P(x) with given degree in terms
  of a set of N+1 basis polynomials Q0(x), Q1(x), ..., QN(x) such that:
  P(x) = a[0]*Q0[x] + a[1]*Q1[x] + ... + a[N]*QN[x], this function returns the expansion
  coefficients b[k] with respect to another set of basis polynomials R, such that:
  P(x) = b[0]*R0[x] + b[1]*R1[x] + ... + b[N]*RN[x]
  The basis polynomials Q and R are passed as 2-dimensional arrays where the k-th row represents
  the coefficients of the k-th basis polynomial. If R is not a basis, the function will not succeed
  and return false, otherwise true. */
  static bool baseChange(T **Q, T *a, T **R, T *b, int degree)
  {
    return rsLinearAlgebra::rsChangeOfBasisRowWise(Q, R, a, b, degree+1);
  }
  // make const-correct - first make functions in rsLinearAlgebra const-correct

  /** Computes polynomial coefficients from the roots. \todo: get rid of that - replace by function
  below */
  static std::vector<std::complex<T>> rootsToCoeffs(const std::vector<std::complex<T>>& roots);
  // allocates heap memory

  /** Computes polynomial coefficients from the roots. The roots should be passed in the array "r"
  of length "N", the coefficients will be returned in the array "a" of length "N" + 1. The
  coefficient for the highest power a[N] will be normalized to unity. */
  static void rootsToCoeffs(const std::complex<T> *r, std::complex<T> *a, int N);
  // allocates heap memory
  // rename to finiteRootsToCoeffs

  /** Similar to rootsToCoeffs(Complex *r, Complex *a, int N), but assumes that the roots are
  either real or occur in complex conjugate pairs. This means that the polynomial has purely real
  coefficients, so the type of the coefficient-array is T instead of Complex. You should use
  this function only if you know in advance that the coefficients will indeed come out as purely
  real */
  static void rootsToCoeffs(const std::complex<T> *r, T *a, int N);
  // allocates heap memory

  static void rootsToCoeffs(const T* r, T* a, int N, T scaler = T(1));

  // drag the ..shiftArgument function in this group


  //-----------------------------------------------------------------------------------------------
  /** \name Fitting/Interpolation */

  /** Computes coefficients a[0], a[1], a[2], a[3] for the cubic polynomial that goes through the
  points (x[0], y[0]) and (x[1], y[1]) and has first derivatives of dy[0] and dy[1] at these points
  respectively. */
  static void cubicCoeffsTwoPointsAndDerivatives(T *a, const T *x, const T *y, const T *dy);
  // rename to fitCubic or cubicHermiteCoeffs

  /** Simplified version of cubicCoeffsTwoPointsAndDerivatives(T *a, T *x, T *y, T *dy)
  which assumes that x[0] = 0, x[1] = 1 - so we don't need to actually pass any x-array. The 
  formulas are far simpler. */
  static void cubicCoeffsTwoPointsAndDerivatives(T *a, const T *y, const T *dy);
  // rename to fitCubic_0_1 or cubicHermiteCoeffs_0_1

  // make a function fitQuartic with the same constraints as the cubic Hermite with integral 
  // normalized to that of the linear interpolant

  // \todo void cubicCoeffsFourPoints(T *a, T *x, T *y);

  /** Computes coefficients a[0], a[1], a[2], a[3] for the cubic polynomial that goes through the
  points (-1, y[-1]), (0, y[0]), (1, y[1]), (2, y[2]). NOTE, that the y-array is accessed at values
  y[-1]...y[2] - the caller should make sure, these values exist. */
  static void cubicCoeffsFourPoints(T *a, const T *y);
  // rename to fitCubic_m1_0_1_2 or cubicLagrange_m1_0_1_2 or cubicThrough

  /** Allocates and fills an NxN matrix A wher A[i][j] are given by x[i]^j. The caller is
  responsible for deallocation. So it's used like:
  T **A = rsVandermondeMatrix(x, N);
  // ...do stuff with matrix A
  rsDeAllocateSquareArray2D(A, N);  */
  static T** vandermondeMatrix(const T *x, int N);
    // move to rsMatrixOld

  /** Computes coefficients a[0],..., a[N-1] for a polynomial of degree N-1 that goes through the N
  data points (x[0], y[0]),...,(x[N-1], y[N-1]). */
  static void interpolant(T *a, const T *x, const T *y, int N);
    // maybe move to Interplation

  /** Like rsInterpolatingPolynomial(T *a, T *x, T *y, int N), but instead of
  passing an x-array, you should pass a start value x0 and an increment dx and it will use x-values
  given by x0, x0+dx, x0+2*dx,...,x0+(N-1)*dx, so this function assumes equidisant abscissa
  values. */
  static void interpolant(T *a, const T& x0, const T& dx, const T *y, int N);
  // allocates heap memory

  // \todo void quinticCoeffsTwoPointsAndDerivatives(T *a, T *x, T *y, T *dy,
  //                                                 T *d2y);

  /** Fits the quadratic parabola defined by y(x) = a[2]*x^2 + a[1]*x + a[0] to the
  3 points (x[0],y[0]), (x[1],y[1]), (x[2],y[2]). */
  static void fitQuadratic(T *a, const T *x, const T *y);

  /** Fits the quadratic parabola defined by y(x) = a[2]*x^2 + a[1]*x + a[0] to the
  3 points (0,y[0]), (1,y[1]), (2,y[2]). */
  static void fitQuadratic_0_1_2(T *a, const T *y);

  /** Fits the quadratic parabola defined by y(x) = a[2]*x^2 + a[1]*x + a[0] to the
  3 points (-1,y[0]), (0,y[1]), (1,y[2]). */
  static void fitQuadratic_m1_0_1(T *a, const T *y);

  /** Returns the position (i.e. the x-coordinate) of the extremum (minimum or maximum) of the 
  quadratic parabola y(x) = a[2]*x^2 + a[1]*x + a[0]. Note that a[2] must be nonzero, otherwise the
  parabola degenerates to a line and there will be no extremum - which will lead to a division by
  zero in the formula. */
  static T quadraticExtremumPosition(const T *a);
  // maybe inline this, move to Evaluation group

  /** Fits the quartic defined by y(x) = a[4]*x^4 + a3*x^3 + a[2]*x^2 + a[1]*x + a[0] to the
  3 points (0,y[0]), (1,y[1]), (2,y[2]) and also matches the derivatives (slopes)
  y'(0) = s0, y'(2) = s2 */
  static void fitQuarticWithDerivatives(T *a, const T *y, const T& s0, const T& s2);

  /** Given coefficients of a polynomial a2*x^2 + a1*x + a0, this function determines whether its
  roots are on or inside the unit circle.
   \todo: write and run a unit-test for this function. */
  static bool areRootsOnOrInsideUnitCircle(const T& a0, const T& a1, const T& a2);
  // move to Evaluation

  // \todo fitPolynomial(T *a, int order, T *x, T *y, int numValues);
  // degree+1 == numValues: exact fit
  // degree+1 >  numValues: exact fit, some higher coeffs unused -> maybe via recursive call
  // degree+1 <  numValues: least-squares fit


  //-----------------------------------------------------------------------------------------------
  /** \name Coefficient generation */
  // functions to generate coefficient arrays for certain special polynomials

  /** Computes polynomial coefficients of a polynomial that is defined recursively by
  w0 * P_n(x) = (w1 + w1x * x) * P_{n-1}(x) + w2 * P_{n-2}(x)
  where n is the degree of the polynomial represented by the "a" array, the degree of a1 is n-1, 
  the degree of a2 is n-2. The lengths of the corresponding arrays equals their respective degree 
  plus 1. w0, w1, w2 are the weighting coeffients of the linear 3-term recurrence relation. The 
  pointer for result "a" may point to the same memory location as either of the input argument 
  arrays "a1", "a2", so the function may be used in place. */
  static void threeTermRecursion(T *a, const T& w0, int degree, const T *a1, const T& w1, 
    const T& w1x, const T*a2, const T& w2);

  /** Fills the array with coefficients for a Bessel-polynomial of given degree. */
  static void besselPolynomial(T *a, int degree);
  // allocates heap memory
  // todo: maybe use rsPolynomialRecursion inside

  /** Fills the array with coefficients for a Legendre-polynomial (of the 1st kind) of given
  degree. */
  static void legendrePolynomial(T *a, int degree);
  // allocates heap memory
    // todo: maybe use rsPolynomialRecursion - or maybe get rid of the function
    // (move to prototypes)

  /** Computes the recursion coefficients (as used in rsPolynomialRecursion) for the Jacobi
  polynomial of degree n (n >= 2) with parameters a and b. */
  static void jacobiRecursionCoeffs(int n, T a, T b, T *w0, T *w1, T *w1x, T *w2);

  /** Given the coefficients of the Jacobi polynomials of degrees n-1 and n-2 in c1 and c2, this
  function computes the coefficients of the next Jacobi polynomial of degree n by using the 3 term
  recurrence relation with parameters a and b. */
  static void jacobiRecursion(T *c, int n, T *c1, T *c2, T a, T b);

  /** Computes the coefficients of the Jacobi polynomials of degrees up to maxDegree with 
  parameters a and b (usually alpha and beta in formulas). The 2D array "c" will contain the 
  coefficients on return. The first index in the c-array runs over indices 0...maxDegree inclusive,
  so the outer dimension should be maxDegree+1. The "c" array may actually be triangular, with the 
  1st inner array of length 1 , the 2nd of length 2, etc. where the last one should have length
  maxDegree+1 (same as the outer). You may also use a square matrix for convenience - then unused
  elements will not be touched in this case. */
  static void jacobiPolynomials(T **c, T a, T b, int maxDegree);

  /** Analog to jacobiRecursion */
  static void legendreRecursion(T *a, int n, T *a1, T *a2);

  // todo: void chebychevPolynomial(T *a, int degree);

  /** Constructs a polynomial p(x) of degree 2*N+1 with the following properties:
  p(0) = 0, p(1) = 1, p'(x) >= 0 for all x (monotonically increasing), p'(1) = maximum possible
  when monotonicity is assumed. \todo: check if these properties are actually true. Such
  polynomials are used in Papoulis filters. */
  static void maxSlopeMonotonic(T *a, int N);
  // allocates heap memory

  // \todo for Halpern filters:
  //void jacobiPolynomial(T *a, int degree); // the U-polynomials
  //void maximallyDivergingMonotonicPolynomial(T *a, int degree); // the T-polynomial



protected:

  std::vector<T> coeffs;   // array of coefficients - index correpsonds to power of x

  // Some functions and classes need low-level access to the coefficient array. These are 
  // declared as friends here. The friends should all themselves be part of the library. We don't
  // want external friends (except maybe during development):
  template<class U> friend class rsRationalFunction;
  
  //friend rsPolynomial<T> fitPolynomial(int numDataPoints, T* x, T* y, int degree);
  //template<class U> friend rsPolynomial<U> ::fitPolynomial(int numDataPoints, U* x, U* y, int degree);

};

// todo: implement a function that determines the number of real roots of a polynomial in an
// interval by means of Sturmian sequences (see Einführung in die computerorientierte Mathematik
// mit Sage, p.163ff) - may be useful for rational interpolants to figure out, if there's a pole
// inside the interpolation interval

// make all the functions const correct

// todo: implement polynomial GCD (i.e. adapt the integer gcd for polynomials)

// \todo implement funcitons to construct coefficient arrays for certain recursively defined
// polynomials, such as Chebychev, etc. provide functions to evaluate linear combinations of
// them efficiently (maybe look up Clenshaw algorithm)

// wrap all (or many) of these functions into a class rsPolynom as static member functions
// let class polynom have members T *a; int order; where a is the coefficient array an T is the
// type for the coefficients. implement operators +,-,*,/,% where the two latter should call
// a general rsDivMod(T n, T d, T *q, T *r); function n: numerator, d: denominator, q: quotient
// r: remainder
// make a class rsRational<T>, where T can be an integer type or class rsPolynomial. In the
// division of rsRational, the rsDivMod<T> function may be used
// make an subclass rsRationalFunction : public rsRational<rsPolynom> that implements
// integretaion and differentiation
// maybe both rsPolynomial and rsRationalFunction should also be subclasses of
// rsUnivariateScalarFunction

#endif
