#ifndef RAPT_POLYNOMIAL_H
#define RAPT_POLYNOMIAL_H

/** A class for representing polynomials and doing computations with them. Much of the code is 
implemented as static member functions that operate directly on arrays of the type T, so you don't
have to create an instance of class rsPolynomial to use its functionality. However, for 
convenience, you may also instantiate polynomial objects and then you can do arithemetic operations
with these objects directly, for example writing code like:

  rsPolynomial<double> r = p*q;

where p and q are both polynomials and so is r - and the * operator implements polynomial 
multiplication. But doing it this way is recommended mostly for prototyping only because creating
polynomials involves dynamic memory allocations. For production (especially real-time) code, it's 
better to operate on pre-allocated arrays with the static functions - this makes the code harder to 
read but more efficient.

There are some static functions that use their own template parameter types independent from the 
"T" that is used to instantiate the class. That's the case, for example, for functions that expect 
real inputs and produce complex outputs (like root-finders) which then use "R" for the real type
and std::complex<R> for the complex type. This is done because this class template should be able
to be instantiated for real and complex types "T", so using the same template parameter could lead
to confusion like the compiler using a nested complex type which makes no sense.
....under construction....tbc...  */

// todo: 
// -in the low-level function interfaces, use consistently "degree" or "aDeg", "bDeg" etc. 
//  instead of the generic N, aN, bN etc. to make it clear to the caller that the degree must be 
//  passed and *NOT* the length of the coefficient array (which is one more than the degree). 
//  Conflating the two is a constant source of confusion and off-by-one bugs and even heap 
//  corruptions.
// -use consistently input arrays as first and output arrays as last parameters for example in
//  interpolant, fitQuadraticDirect -> this will silenty break client code, so be extra careful to
//  adapt the code at *every* call site - it should be a consistent pattern through the library:
//  inputs first, then outputs...hmm - but that doesn't really work well when we want to have some
//  inputs optional - optional parameters must come last...hmmmm....maybe use:
//  required inputs, outputs, optional inputs - what about optional outputs? (for example, when we
//  fill an output-array only when a non-nullptr is passed?)...see what rsArrayTools does...
//  ...or maybe such an enforced consistency might not be a good idea, after all?
// -implement greatest common divisor algorithm (maybe the one for the integers can be
//  used as is?)
// -implement missing operators:
//  -operator() that takes a polynomial and returns another polynomial as result (implements 
//   nesting)
//  -arithmetic operators that take a number as second (left or right) argument
//  -maybe we could meaningfully define <,<=, ...? look at how python or other scientific libraries
//   handle that - in my own python polynomial class, i'm taking the asymptotic behavior


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
  //rsPolynomial(const std::vector<T>& coefficients) : coeffs(coefficients) {}  // why not?
  rsPolynomial(const std::vector<T>& coefficients) { setCoeffs(coefficients); }

  /** Creates a polynomial from an initializer list for the coefficients. */
  rsPolynomial(std::initializer_list<T> l) : coeffs(l) {}
  //  needs tests
  // maybe take a reference? hmm...no, the examples here also use a value:
  // https://en.cppreference.com/w/cpp/utility/initializer_list
  // https://www.learncpp.com/cpp-tutorial/10-7-stdinitializer_list/

  /** Promotes a number to a 0th degree polynomial. */
  rsPolynomial(const T& number) { coeffs.resize(1); coeffs[0] = number; }

  // make a constructor that initializes from a raw array


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

  /** Stretches (or compresses) the polynomial along the x axis. */
  void stretch(T factor)
  { scaleArgument(&coeffs[0], &coeffs[0], getDegree(), T(1)/factor); }

  /** Scales the output of the whole polynomial by the given factor. */
  void scale(T factor)
  { rsArrayTools::scale(&coeffs[0], (int) coeffs.size(), factor); }


  void shiftX(T dx) 
  {
    std::vector<T> tmp = coeffs;
    shiftArgument(&tmp[0], &coeffs[0], getDegree(), dx);
  }
  // needs tests

  /** Shifts the polynomial up and down in the y direction by the given dy. */
  void shiftY(T dy) { coeffs[0] += dy; }

  /** Turns this polynomial into the indefinite integral of itself with integration constant c 
  (this c becomes the coeff fo x^0 = 1). */
  void integrate(T c = T(0))
  { 
    coeffs.resize(coeffs.size()+1);
    integral(&coeffs[0], &coeffs[0], getDegree()-1, c); // -1 bcs resize has increased degree
  }

  /** Sets the allcoated degree of this polynomial. The "allocated" qualifier means, that we are 
  talking about the length of the coeff array (minus 1) regardless whether the highest coeff in 
  this array is zero or not. If the new degree is less than the old one, we'll just cut off the 
  coefficient array at the given new endpoint. If the new degree is greater than the old one, we
  take over the lower coefficients and fill higher coefficients with zero. */
  void setAllocatedDegree(int newDegree) { coeffs.resize(newDegree+1); }

  /** Adds the givne polynomial q multiplied by a scalar weight into this one. If q has higher 
  degree thatn this one, the degree of this will be increased to accomodate for the higher 
  coefficients in q - so the function may reallocate memory. */
  void addWithWeight(const rsPolynomial<T>& q, T w)
  {
    int qDeg = q.getDegree();
    if(qDeg > getDegree())
      setAllocatedDegree(qDeg);
    for(int i = 0; i <= qDeg; i++)
      coeffs[i] += w * q.coeffs[i];
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the maximum order that this poloynomial may have which is the length of the
  coefficient array minus one. When there are trailing zero coefficients, the actual degree of
  the polynomial is lower. */
  //int getMaxOrder() const { return (int)coeffs.size()-1; }
  // deprecate this 

  /** Returns the degree of the polynomial. Mathematically, this is defined as the exponent of the 
  highest power of x which has a nonzero coefficient....but the function currently just returns the 
  degree as it is determined by the length of the coefficient array, not checking if that last 
  value is zero because doing so would need some tolerance when using floating point numbers and 
  i'm not yet sure, how to best handle that...typically, client code wants to know, how long the
  coefficient array is anyway....tbc... */
  int getDegree() const { return (int)coeffs.size()-1; }
  // should take into account trailing zeros ..or maybe have a boolean flag
  // "takeZeroCoeffsIntoAccount" which defaults to false...or maybe it shouldn't have any default
  // value - client code must be explicit...or maybe have functions getAllocatedDegree, 
  // getActualDegree(tolerance)...or getDegree has an optional parameter for the tolerance 
  // defaulting to 0...but no - the calls may still be ambiguous from client code, when nothing is 
  // passed...it's also confusing to pass a number into such a getter

  /** Returns the number of coefficients in this polynomial. */
  int getNumCoeffs() const { return (int)coeffs.size(); }
  // maybe client code should preferably use this, when it wants to know the length of the coeff 
  // array and not getDegree because of the ambiguity

  /** Returns the leading coefficient, i.e. the coefficient that multiplies the highest power of 
  x. */
  T getLeadingCoeff() const { return rsLast(coeffs); }
  // what if we have trailing zeros in the coeff array?

  /** Returns the i-th coefficient, i.e. the coefficient for x^i. */
  T getCoeff(int i) const { rsAssert(i >= 0 && i <= getDegree()); return coeffs[i]; }

  /** Returns the i-th coefficient, i.e. the coefficient for x^i or zero if i is greater than the 
  degree of this polynomial. */
  T getCoeffPadded(int i) const 
  { 
    rsAssert(i >= 0); 
    if(i > getDegree())
      return T(0);
    return coeffs[i]; 
  }

  /** Returns a pointer to our coefficient array - breaks encapsulation - use with care! */
  T* getCoeffPointer() { return &coeffs[0]; }
  // Try to get rid of this - when we really need low-level access to the coeff-array, we declare 
  // the functions/classes that need it as friends. Comment this out later

  const T* getCoeffPointerConst() const { return &coeffs[0]; }

  /** Returns true, iff this polynomial is monic, i.e. the coefficient for the highest power (the
  leading coefficient) is unity. Monic polynomials are important because they arise when 
  multiplying out the product form. */
  bool isMonic() const { return getLeadingCoeff() == T(1); }
  // what if we have trailing zeros in the coeff array? should we have a tolerance?






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

  /** Adds the polynomial b to "this" polynomial. */
  rsPolynomial<T>& operator+=(const rsPolynomial<T>& b) 
  { 
    return *this = (*this) + b;
    // this implementation is preliminary - todo: optimize: if deg(b) > deg(this) -> resize, then
    // loop through the coeffs up to min(deg(b), deg(this)) using +=
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
  //T operator()(T x) const { return evaluate(x, &coeffs[0], getDegree()); }
  T operator()(T x) const { return evaluate(x); }

  /** Overloaded evaluation operator () that takes a polynomial as input and returns another 
  polynomial. This implements nesting/composition. The given x is the inner polynomial and "this" 
  is the outer polynomial */
  rsPolynomial<T> operator()(const rsPolynomial<T>& p) const
  {
    rsPolynomial<T> r(getDegree() * p.getDegree());
    compose(&p.coeffs[0], p.getDegree(), &coeffs[0], getDegree(), &r.coeffs[0]);
    return r;
  }

  /** Read and write access to i-th coefficient (breaks encapsulation - use with care). */
  T& operator[](int i) { return coeffs[i]; }





  //-----------------------------------------------------------------------------------------------
  /** \name Evaluation (High Level) */

  /** Evaluates the polynomial at the given input x. */
  T evaluate(T x) const { return evaluate(x, &coeffs[0], getDegree()); }

  /** Evaluates the first derivative of this polynomial at the given x. */
  T derivativeAt(const T& x) 
  { return evaluateDerivative(x, &coeffs[0], getDegree()); }

  /** Evaluates the order-th derivative of this polynomial at the given x. Works also for the 0th 
  derivative, which is the function value itself. ...but the order must be non-negative. */
  T derivativeAt(const T& x, int order) 
  { return evaluateDerivative(x, &coeffs[0], getDegree(), order); }
  // todo: maybe make it also work for negative orders (in which case the antiderivative of 
  // given order will be evaluated (setting integration constants to zero))

  T integralAt(const T& x, const T c = T(0))
  { return evaluateIntegral(x, &coeffs[0], getDegree(), c); }

  T definiteIntegral(const T& lowerLimit, const T& upperLimit)
  { return integralAt(upperLimit) - integralAt(lowerLimit); }



  //-----------------------------------------------------------------------------------------------
  /** \name Calculus (High Level) */

  rsPolynomial<T> derivative() const
  { 
    if(getDegree() == 0) return rsPolynomial<T>(); // a constant polynomial has zero derivative
    rsPolynomial<T> d(getDegree()-1); derivative(&coeffs[0], &d.coeffs[0], getDegree()); 
    return d;
  }



  //===============================================================================================
  /** \name Computations on raw coefficient arrays (low level functions) */

  //-----------------------------------------------------------------------------------------------
  /** \name Evaluation (Low Level) */

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
  template<class R>
  static R evaluateWithTwoDerivativesAndError(const std::complex<R>* a, int degree,
    std::complex<R> z, std::complex<R>* P);
  // rename "P" to "y" ...or "w2 as is common in complex functions

  /** Evaluates the indefinite integral of the polynomial at given x with given integration 
  constant c. */
  static T evaluateIntegral(const T& x, const T *a, int degree, T c = T(0));


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
  // maybe also implement evaluateHermiteMonic which should be the probabilist's version

  // todo: evaluateDerivative, evaluateIntegral (or AntiDerivative)

  /** Evaluates the Newton polynomial with the Newton expansion coeffs c[0],...,c[N] and roots 
  r[0],...,r[N-1] given by:
    p(x) = c[0] + c[1]*(x-r[0]) + c[2]*(x-r[0])*(x-r[1]) + ... + c[N]*(x-r[0])*...*(x-r[N-1])
  at the given input x. Array c is of length N+1, r is of length N. */
  static T evaluateNewton(const T& x, const T* c, const T* r, int N);


  /** Given a coefficient array p of length maxDegree+1, this function returns the actual degree
  of the polynomial - meaning that trailing zeros in the array don't count. So it's the degree 
  that takes only into account the non-zero coefficients. */
  //static int actualDegree(T* p, int maxDegree, T tol = T(0));
  // maybe move to an "Inquiry" section ...or "Misc"

  // maybe abbreviate evaluate with eval, add evaluation functions for evaluating polynomials in
  // different bases, for example, using chebychev polynomials as basis (instead of the regular
  // monomial basis x^0,x^1,x^2,x^3,... use T0(x),T1(x),T2(x),T3(x),...) - also for Newton 
  // polynomials (they are additionally parametrized by a set of roots)

  //-----------------------------------------------------------------------------------------------
  /** \name Arithmetic (Low Level) */

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

  /** Like powers(const T* a, int N, T** aPowers, int highestPower) but using a flat array for the
  output ...tbc... */
  static void powers(const T* a, int N, T* aPowers, int highestPower, int stride);

  /** Let A(x) and B(x) be polynomials represented by their coefficient arrays a[] and b[]
  respectively. This function creates the coefficients of a polynomial C(x), represented by the
  coefficient array c[], that results from composing the polynomials A(x) and B(x), that is: first
  the polynomial A(x) is applied to the value x, and then the polynomial B(x) is applied to the
  result of the first polynomial, such that C(x) = B(A(x)). This nesting or composition of two
  polynomials can itself be seen as a polynomial in its own right. This resulting polynomial has
  a degree of cN = aN*bN, where aN and bN are the degrees of the a[] and b[] polynomials,
  respectively, so the caller has to make sure that the c[] array has at least a length of
  aN*bN+1. The workspace must also be of length aN*bN+1. */
  static void compose(const T* a, int aN, const T* b, int bN, T* c, T* workspace);
  // i think, the complexity is O(aN^2 * bN)...verify!

  /** Convenience function that allocates a workspace internally. */
  static void compose(const T* a, int aN, const T* b, int bN, T* c);

  /** Composes (nests) the outer polynomial a(x) = a0 + a1*x + a2*x^2 + a3*x^3 with the inner
  polynomial b(x) = b0 + b1*x and writes the resulting cofficients into c (which may point to the
  same array as a for in-place use). */
  static void composeLinearWithCubic(T* a, T* c, T b0, T b1);

  /** Given an array of polynomial coefficients "a" such that
  p(x) = a[0]*x^0 + a[1]*x^1 + ... + a[N]*x^N, this function returns (in "am") the coefficients for
  a polynomial q(x) such that q(x) = p(-x). This amounts to sign-inverting all coefficients which
  multiply odd powers of x. */
  static void negateArgument(const T *a, T *am, int N);

  /** Given an array of polynomial coefficients "a" such that
  p(x) = a[0]*x^0 + a[1]*x^1 + ... + a[N]*x^N, this function returns (in "as") the coefficients for
  a polynomial q(x) such that q(x) = p(scaler*x). This amounts to scaling all coefficients with the
  scaler raised to the same power as the respective x. */
  static void scaleArgument(const T *a, T *as, int N, T scaler);
  // maybe rename to stretch or scaleX

  /** Given an array of polynomial coefficients "a" such that
  p(x) = a[0]*x^0 + a[1]*x^1 + ... + a[N]*x^N, this function returns (in "aShifted") the coefficients
  for a polynomial q(x) such that q(x) = p(x-x0). */
  static void shiftArgument(const T *a, T *aShifted, int N, T x0);
  // allocates heap memory and algo is O(N^2) - can this be done better? ..i think so ...yes - we 
  // should use compose and a workspace of length N (or N+1) should suffice, see comments in
  // shiftPolynomial in MathExperiments.cpp
  // maybe move into "Conversions" section
  // maybe rename to shiftX


  //-----------------------------------------------------------------------------------------------
  /** \name Calculus (Low Level) */

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
  /** \name Roots (Low Level) */

  /** Finds all complex roots of a polynomial by Laguerre's method and returns them in "roots". */
  template<class R>
  static void roots(const std::complex<R>* a, int degree, std::complex<R>* roots);
  // allocates heap memory

  /** Same, but for real coefficients */
  template<class R>
  static void roots(const R *a, int degree, std::complex<R> *roots);
  // allocates heap memory

  /** Converges to a complex root of a polynomial by means of Laguerre's method using the
  "initialGuess" as first estimate. */
  template<class R>
  static std::complex<R> convergeToRootViaLaguerre(const std::complex<R> *a, int degree,
    std::complex<R> initialGuess = std::complex<R>(0.0, 0.0));
  // allocates heap memory

    /** Computes the root of the linear equation: \f[ a x + b = 0 \f] which is simply given by
  \f[ x_0 = -\frac{b}{a} \f] */
  static T rootLinear(const T& a, const T& b);
  // rename inputs to a0,a1 (change their order)

  /** Computes the two roots of the quadratic equation: \f[ a x^2 + b x + c = 0 \f] which are
  given by: \f[ x_{1,2} = \frac{-b \pm \sqrt{b^2-4ac}}{2a} \f] and stores the result in two-element
  array which is returned. When the qudratic is degenerate (i.e, a == 0), it will fall back to the
  rootsLinear() function, and return a one-element array.  */
  template<class R>
  static std::vector<std::complex<R>> rootsQuadratic(const R& a, const R& b, const R& c);
  // rename inputs to a0,a1,a2 (change their order) ..or maybe deprecate this function

    /** Computes the two roots of the quadratic equation: \f[ a_0 + a_1 x + a_2 x^2 = 0 \f] and
  stores them in r1, r2. When the equation has two distinct real roots, they will be returned in
  ascending order, i.e. r1 < r2. In case of a (real) double root, we'll have r1 == r2 and when the
  roots of the equation are actually complex, the outputs will also be equal and contain the real
  part of the complex conjugate pair. */
  template<class R>
  static void rootsQuadraticReal(const R& a0, const R& a1, const R& a2, R* r1, R* r2);

  template<class R>
  static void rootsQuadraticComplex(
    const std::complex<R>& a0, const std::complex<R>& a1, const std::complex<R>& a2,
    std::complex<R>* x1, std::complex<R>* x2);
    // todo: make optimized version for real coefficients (but complex outputs)

  /** Computes the three roots of the cubic equation: \f[ a x^3 + b x^2 + c x + d = 0 \f] with real
  coefficients and stores the result in the three-element array which is returned. When the cubic 
  is degenerate (i.e, a == 0), it will fall back to the getRootsOfQuadraticEquation() function, and 
  return a two-element array (or a one-element array, when b is also zero). */
  template<class R>
  static std::vector<std::complex<R>> rootsCubic(const R& a, const R& b, const R& c, const R& d);
  // todo: make the order of the arguments consistent with evaluateCubic - but careful - this will
  // break client code! ...rename parameters to a0,a1,a2,a3 to make it clear, how it's meant

  /** Discriminant of cubic polynomial \f[ a_0 + a_1 x + a_2 x^2 + a_3 x^3 = 0 \f].
  D > 0: 3 distinct real roots, D == 0: 3 real roots, 2 or 3 of which may coincide,
  D < 0: 1 real root and 2 complex conjugate roots */
  template<class R>
  static R cubicDiscriminant(const R& a0, const R& a1, const R& a2, const R& a3);
  // rename to discriminantCubic

  // todo: write function quadraticDiscriminant

  /** under construction - does not yet work */
  template<class R>
  static void rootsCubicComplex(
    std::complex<R> a0, std::complex<R> a1, 
    std::complex<R> a2, std::complex<R> a3, 
    std::complex<R>* r1, std::complex<R>* r2, std::complex<R>* r3);


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
  template<class R>
  static R cubicRootNear(R x, const R& a, const R& b, const R& c, const R& d, const R& min, 
    const R& max, int maxIterations = 10);
  // todo: rename to rootCubicNear, change order of variables, maybe use bisection, if 
  // newton-iteration diverges

  /** Iteratively improves an initial estimate for the root of the polynomial equation:
  \f[ a[order] x^order + ... + a[1] x + a[0] = 0   \f]
  by means of the Newton-Raphson iteration:
  \f[ x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}       \f]
  the arguments min and max give upper and lower bounds for the root (which will be returned in
  cases where the iteration diverges, which you should avoid in the first place) and maxIterations
  gives the maximum number of iteration steps. */
  template<class R>
  static R rootNear(R x, const R* a, int order, const R& min, const R& max, 
    int maxIterations = 32);

  /** Same as above but accepts real coefficients. */
  //void findPolynomialRootsInternal(T *a, int order, Complex *roots, bool polish = true);

  //void findPolynomialRootsNew(Complex *a, int order, Complex *roots);


  //-----------------------------------------------------------------------------------------------
  /** \name Conversions (Low Level) */
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

  /** Computes polynomial coefficients from the roots. 
  \todo: get rid of that - replace by function below */
  static std::vector<std::complex<T>> rootsToCoeffs(const std::vector<std::complex<T>>& roots);
  // allocates heap memory - todo: avoid that

  /** Computes polynomial coefficients from the roots. The roots should be passed in the array "r"
  of length "N", the coefficients will be returned in the array "a" of length "N" + 1. The
  coefficient for the highest power a[N] will be normalized to unity. */
  static void rootsToCoeffs(const std::complex<T> *r, std::complex<T> *a, int N);
  // allocates heap memory - todo: avoid that
  // rename to finiteRootsToCoeffs

  /** Similar to rootsToCoeffs(Complex *r, Complex *a, int N), but assumes that the roots are
  either real or occur in complex conjugate pairs. This means that the polynomial has purely real
  coefficients, so the type of the coefficient-array is T instead of Complex. You should use
  this function only if you know in advance that the coefficients will indeed come out as purely
  real */
  template<class R>
  static void complexRootsToRealCoeffs(const std::complex<R> *r, R *a, int N);
  // allocates heap memory

  static void rootsToCoeffs(const T* r, T* a, int N, T scaler = T(1));

  /** Given the coefficents for the Newton basis polynomials 1,(x-x[0]),(x-x[0])*(x-x[1]),... in 
  a[i], i = 0,..,N-1, this function converts the coefficients to the monomial basis, i.e. into 
  regular polynomial coefficients. The result is stored in the same array a. The array x is 
  destroyed during the process - it's re-used internally as temporary buffer for intermediate 
  results to allow higher level code optimize memory usage (otherwise, the function would need an 
  additional buffer of length N - if desired, the caller can create this additional buffer itself 
  and copy the x-values into it and then use this function). */
  static void newtonToMonomialCoeffs(T* x, T* a, int N);


  // drag the ..shiftArgument function in this group

  //-----------------------------------------------------------------------------------------------
  /** \name Fitting/Interpolation (Low Level) */

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
  // -rename to fitCubic_m1_0_1_2 or cubicLagrange_m1_0_1_2 or cubicThrough 
  //  -document why do we adopt the convention that we start x at -1 and access the y[-1] element - i 
  //   think, this is more convenient in the content of interpolating in realtime as in the 
  //   pitch-detector?

  /** Allocates and fills an NxN matrix A where A[i][j] are given by x[i]^j. The caller is
  responsible for deallocation. So it's used like:
  T **A = rsVandermondeMatrix(x, N);
  // ...do stuff with matrix A
  rsDeAllocateSquareArray2D(A, N);  */
  static T** vandermondeMatrix(const T *x, int N);
  // move to rsMatrixOld, deprecate! ..implement the algo with the new rsMatrix/rsLinearAlgebra
  // implementation - the implementation is actually just for reference anyway - for production 
  // code, the algorithm using the Newton polynomials should be used

  /** Computes coefficients a[0],..., a[N-1] for a polynomial of degree N-1 that goes through the N
  data points (x[0], y[0]),...,(x[N-1], y[N-1]). 
  Note that N here is the number of datapoints, not the degree (which is one less)...that's a bit 
  quirky....  */
  static void interpolant(T *a, const T *x, const T *y, int numDataPoints);
  // maybe move to Interplation
  // the meaning of N here is inconsistent with the rest of the class - it's the number of 
  // datapoints - maybe move interpolant into some interpolator class where its conventional to 
  // pass the number of datapoints.

  /** Version of interpolant that avoids memory allocation by letting use pass a workspace - this
  must be at least of length N+1. */
  static void interpolant(T* a, const T* x, const T* y, int N, T* workspace);

  /** Like rsInterpolatingPolynomial(T *a, T *x, T *y, int N), but instead of
  passing an x-array, you should pass a start value x0 and an increment dx and it will use x-values
  given by x0, x0+dx, x0+2*dx,...,x0+(N-1)*dx, so this function assumes equidisant abscissa
  values. */
  static void interpolant(T *a, const T& x0, const T& dx, const T *y, int numDataPoints);
  // allocates heap memory

  /** Implements Newton's algorithm using diveded difference to find the polynomial interpolant 
  through the N points (x[n], y[n]), n = 0,...,N-1. It stores the polynomial coefficients in "a" 
  and it needs a workspace array of length N. All arrays should be of length N. */
  static void interpolantViaNewton(T* a, const T* x, const T* y, int N, T* workspace);
  // make similar functions interpolantViaLagrange, interpolantViaVandermode - interpolant may call 
  // any of these ...maybe the Lagrange version should have 2 variants - one using the "master"
  // polynomial and dividing out one root at a time - this leads to an O(N^2) algo as well - we may
  // then compare it to the Newton algo

  /** In place version of interpolantViaNewton. Overwrites the x,y arrays during the process. On 
  return, y will contain the polynomial coeffs and x will contain garbage. (More specifically, x 
  will contain the coefficients of the unique monic polynomial of degree N-1 that has roots at the 
  given original x-values, except the last one. I don't think that's any useful for the caller, 
  but just for info). The time complexity is O(N^2) and space complexity is O(1). */
  static void interpolantViaNewtonInPlace(T* x, T* y, int N);


  // \todo void quinticCoeffsTwoPointsAndDerivatives(T *a, T *x, T *y, T *dy,
  //                                                 T *d2y);

  /** Fits the quadratic parabola defined by y(x) = a[2]*x^2 + a[1]*x + a[0] to the
  3 points (x[0],y[0]), (x[1],y[1]), (x[2],y[2]). */
  static void fitQuadratic(T* a, const T* x, const T* y)
  { fitQuadraticDirect(a, x, y); }
  // todo: figure which formula is better, the direct or the Lagrange-formula - choose the better
  // one

  /** Uses a formula that resulted from setting up the 3x3 linear system of equations 
        y[i] = a0 + a1*x[i] + a2*x[i]^2   i = 0,1,2 
  and solving it directly. */
  static void fitQuadraticDirect(T *a, const T *x, const T *y);

  /** Computes the same thing as fitQuadraticDirect but uses formulas derived from setting up 3 
  polynomials in product form, where each has zeros at all but one of the datapoints, say xi, and
  to have value yi at xi and then adding them up (idea due to Lagrange):
     p1(x) = k1*(x-x2)*(x-x3)       p1 has zeros at at x2,x3
     p2(x) = k2*(x-x1)*(x-x3)       p2 has zeros at at x1,x3
     p3(x) = k3*(x-x1)*(x-x2)       p3 has zeros at at x1,x2
  Require:
     p1(x1) = y1, p2(x2) = y2, p3(x3) = y3
  Solve these for the ki, i.e. k1,k2,k3. For example, k1 = y1 / ((x1-x2)*(x1-x3)). Plug, for 
  example, k1 back into the p1 equation and multiply it out to obtain its coeffs - do the same 
  for p2 and p3 and then obtain the final polynomial coeffs by adding the corresponding  coeffs 
  of each of the partial polynomials. This organizes the computations in a more orderly manner but
  i'm not yet sure if it's more numerically accurate and/or efficient than the direct way - tests 
  are needed. */
  static void fitQuadraticLagrange(T *a, const T *x, const T *y);

  /** Fits the quadratic parabola defined by y(x) = a[2]*x^2 + a[1]*x + a[0] to the
  3 points (0,y[0]), (1,y[1]), (2,y[2]). */
  static void fitQuadratic_0_1_2(T *a, const T *y);
    // can this be used in place, i.e. a==y? if so, document that! ...i think so

  /** Fits the quadratic parabola defined by y(x) = a[2]*x^2 + a[1]*x + a[0] to the
  3 points (-1,y[0]), (0,y[1]), (1,y[2]). */
  static void fitQuadratic_m1_0_1(T *a, const T *y);
    // can this be used in place, i.e. a==y? if so, document that! ...i don't think so but maybe it
    // can be modified for use in place

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

  /** Given coefficients of a polynomial a2*x^2 + a1*x + a0 with real coefficients, this function 
  determines whether its roots are on or inside the unit circle.
   \todo: write and run a unit-test for this function. */
  template<class R>
  static bool areRootsOnOrInsideUnitCircle(const R& a0, const R& a1, const R& a2);
  // move to Evaluation

  // \todo fitPolynomial(T *a, int order, T *x, T *y, int numValues);
  // degree+1 == numValues: exact fit
  // degree+1 >  numValues: exact fit, some higher coeffs unused -> maybe via recursive call
  // degree+1 <  numValues: least-squares fit


  //-----------------------------------------------------------------------------------------------
  /** \name Coefficient generation (Low Level) */
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
  // rename to bessel (we are already inside class rsPolynomial, so the "Polynomial" part is 
  // redundant), or maybe besselCoeffs or coeffsBessel (having coeffs first orders the functions
  // more meaningfully alphabetically)

  /** Fills the array with coefficients for a Legendre-polynomial (of the 1st kind) of given
  degree. */
  template<class R>
  static void legendrePolynomial(R *a, int degree);
  // allocates heap memory
    // todo: maybe use rsPolynomialRecursion - or maybe get rid of the function
    // (move to prototypes)
  // rename to legendre ("Polynomial" is redundant)

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
  template<class R>
  static void maxSlopeMonotonic(R *a, int N);
  // rename to coeffsMaxSlopeMonotonic
  // allocates heap memory

  // \todo for Halpern filters:
  //void jacobiPolynomial(T *a, int degree); // the U-polynomials
  //void maximallyDivergingMonotonicPolynomial(T *a, int degree); // the T-polynomial

  /** Computes coefficients for Newton basis polynomials 1,(x-x[0]),(x-x[0])*(x-x[1]),... using 
  divided differences from N pairs (x[i],y[i]), i = 0,...N-1. It overwrites the y-array with the 
  coefficients. @see evaluateNewton. */
  static void coeffsNewton(const T* x, T* y, int N);
  // todo: make a version that doesn't overwrite the y-array and instead accepts another array
  // for the coeffs to write into
  // maybe rename to coeffsNewtonInPlace

  //-----------------------------------------------------------------------------------------------
  // Evaluation of special polynomials (Low Level)
  // maybe move into the Evaluation section

  // move the evaluateHermite function here - use consistent naming - either they should all start
  // with "evaluate" or none should

  /** Evaluates the N-th degree Chebychev polynomial T_N(x) at x by recursion. */
  static T chebychevRecursive(T x, int N)
  {
    rsAssert(N >= 0, "polynomial degree must be non-negative");
    T t0 = T(1); T t1 = x; T tn = T(1);
    for(int i = 0; i < N; i++) {
      tn = T(2)*x*t1 - t0; t0 = t1; t1 = tn; }
    return t0;
  }
  // maybe rename to evalChebyRecursive, have also a function coeffsCheby

  /** Evaluates the N-th degree Chebychev polynomial T_N(x) at x by means of acos and cos or 
  acosh and cosh. */
  template<class U>
  static U chebychevDirect(U x, int N)
  {
    rsAssert(N >= 0, "polynomial degree must be non-negative");
    if(rsAbs(x) <= U(1)) return  cos( U(N)*acos ( x));
    if(      x  >  U(1)) return  cosh(U(N)*acosh( x));
    if(rsIsEven(N))      return  cosh(U(N)*acosh(-x));
    else                 return -cosh(U(N)*acosh(-x));
  }
  // we can't use template parameter T here because of compiler errors when instantiating 
  // rsPolynomial for std::complex

  // todo: figure out for which N which of the two functions is faster and/or more accurate - maybe
  // provide a dispatcher function - it seems, at least for lower degrees, the recursion is more 
  // accurate - well, at least for inputs that are exactly representable (like not-too-large 
  // integers), which is expected because it just does basic arithmetic - as long as every 
  // intermediate result is exactly representable, the recursion will give exact results

  // We have: cos(n*a) = T_n(cos(a)), see: https://www.youtube.com/watch?v=VOM3giwqMJw&t=14m5s
  // ...how about sin(n*a) = S_n(sin(a)) for some other sort of polynomial? is this possible?
  // maybe S_n(x) = sin(n*asin(x)) in analogy with the code above (for -1 <= x <= +1)?
  // see https://en.wikipedia.org/wiki/Chebyshev_polynomials#Trigonometric_definition
  // ...so the U_n polynomials do not work that way (but similar, so maybe also useful)
  // https://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html
  // https://mathworld.wolfram.com/ChebyshevPolynomialoftheSecondKind.html



protected:

  std::vector<T> coeffs;   // array of coefficients - index correpsonds to power of x

  // Some functions and classes need low-level access to the coefficient array. These are 
  // declared as friends here. The friends should all themselves be part of the library. We don't
  // want external friends (except maybe during development):
  template<class U> friend class rsRationalFunction;
  
  //friend rsPolynomial<T> fitPolynomial(int numDataPoints, T* x, T* y, int degree);
  //template<class U> friend rsPolynomial<U> ::fitPolynomial(int numDataPoints, U* x, U* y, int degree);

};

/** Multiplies a number and a polynomial. */
template<class T>
inline rsPolynomial<T> operator*(const T& s, const rsPolynomial<T>& p)
{
  rsPolynomial<T> q = p;
  q.scale(s);
  return q;
}


// todo: implement a function that determines the number of real roots of a polynomial in an
// interval by means of Sturmian sequences (see Einführung in die computerorientierte Mathematik
// mit Sage, p.163ff) - may be useful for rational interpolants to figure out, if there's a pole
// inside the interpolation interval

// make all the functions const correct

// todo: implement polynomial GCD (i.e. adapt the integer gcd for polynomials)

// \todo implement functions to construct coefficient arrays for certain recursively defined
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
