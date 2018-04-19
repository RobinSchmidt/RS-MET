#ifndef RAPT_POLYNOMIAL_H
#define RAPT_POLYNOMIAL_H

/** A class for representing polynomials and doing computations with them. */

template<class T>
class rsPolynomial
{

public:

  // todo:  have a std::vector<T> member for the coeffs and define operators


  //===============================================================================================
  /** \name Static functions for computing with coefficient arrays */

  /** Evaluates the polynomial defined by the array of roots "r" at argument "x". */
  static std::complex<T> evaluatePolynomialWithRoots(std::complex<T> x, std::complex<T> *r,
    int numRoots);
  // rename to evaluateFromRoots

  /** Evaluates the polynomial defined by the array of coefficients "a" at argument "x".  The array
  of coefficients must be of length order+1 and is interpreted as follows: a[0] is taken to be the
  constant term, a[1] is the multiplier for x^1, a[2] the multiplier for x^2 and so on until
  a[order] which is the multiplier for a^order. */
  static T evaluatePolynomialAt(T x, T *a, int order);
  // rename to evaluateAt

  /** Evaluates the polynomial defined by the array of coefficients 'a' and its first derivative at
  argument 'x'. The value of the polynomial will be stored in y and the value of the derivative
  will be stored in yd. @see: evaluatePolynomialAt() */
  static void evaluatePolynomialAndDerivativeAt(T x, T *a, int order, T *y, T *yd);
  // rename to evaluateWithDerivative

  /** Evaluates the polynomial defined by the array of coefficients 'a' and a given number of
  derivatives at argument 'x'. The value of the polynomial will be stored in results[0], the 1st
  derivative in results[1] and so on. numDerivatives should be <= 31.
  @see: evaluatePolynomialAt() */
  static void evaluatePolynomialAndDerivativesAt(T x, T *a, int order, T *results, 
    int numDerivatives);
  // rename to evaluateWithDerivatives

  /** Multiplies the polynomials represented by the coefficient vectors 'a' and 'b' and stores the
  resulting coefficients in 'result'. The resulting polynom will be or order aOrder+bOrder and the
  coefficient vector should have allocated space for
  (aOrder+1)+(bOrder+1)-1 = aOrder+bOrder+1 = aLength+bLength-1 elements. */
  static void multiplyPolynomials(T *a, int aOrder, T *b, int bOrder, T *result);
  // rename to multiply

  /** Divides the polynomials represented by the coefficient arrays 'dividend' and 'divisor' and
  stores the resulting coefficients for the quotient and remainder in the respective output arrays.
  ?? The resulting quotient polynom will be of order dividendOrder-divisorOrder and the
  remainder polynom will be at most of order ??
  However, the output arrays must have the same length as the dividend, where
  remainder[divisorOrder...dividendOrder] and 
  quotient[dividendOrder-divisorOrder+1...dividendOrder] will be filled with zeros.  
  ...\todo check this comment */
  static void dividePolynomials(T *dividend, int dividendOrder, T *divisor, int divisorOrder, 
    T *quotient, T *remainder);
  // rename to divide

  /** Divides the dividend by the monomial factor (x-x0) and stores the result in the same array
  again. The remainder (which is just a numerical value) will be stored in 'remainder'. */
  template<class S>
  static void dividePolynomialByMonomialInPlace(S *dividendAndResult, int dividendOrder, S x0, 
    S *remainder);

  //static void dividePolynomialByMonomialInPlace(T *dividendAndResult, int dividendOrder, T x0, 
  //  T *remainder);


  /** Given an array of polynomial coefficients "a" such that
  p(x) = a[0]*x^0 + a[1]*x^1 + ... + a[N]*x^N, this function returns (in "am") the coefficients for
  a polynomial q(x) such that q(x) = p(-x). This amounts to sign-inverting all coefficients which
  multiply odd powers of x. */
  static void polyCoeffsForNegativeArgument(T *a, T *am, int N);

  /** Given an array of polynomial coefficients "a" such that
  p(x) = a[0]*x^0 + a[1]*x^1 + ... + a[N]*x^N, this function returns (in "aShifted") the coefficients
  for a polynomial q(x) such that q(x) = p(x-x0). */
  static void polyCoeffsForShiftedArgument(T *a, T *aShifted, int N, T x0);

  /** Finds the coefficients of the derivative of the N-th order polynomial with coeffs in "a" and
  stores them in "ad". The order of the polynomial represented by the coeffs in "ad" will be
  N-1. The "ad" array may point to the same array as the "a" array, so you can use the same array
  for input and output. */
  static void polyDerivative(T *a, T *ad, int N);

  /** Given a polynomial p(x) via its coefficient array a, this function computes the coefficients
  of a polynomial q(x) = p(x+h) - p(x) if direction = 1, or q(x) = p(x) - p(x-h) if direction = -1.
  This represents the first forward or backward difference polynomial in finite difference
  calculus and is analog to the derivative in infinitesimal calculus. As in infinitesimal calculus,
  the resulting polynomial is one order less than p(x) itself. Scaled by 1/h, this can be seen as an
  approximation to the derivative using stepsize h.
  \todo provide a finite central difference q(x) = p(x+h/2) - p(x-h/2) when direction = 0
   ...could this be just the average between forward and backward difference? ...research! */
  static void polyFiniteDifference(T *a, T *ad, int N, int direction = 1, T h = 1);

  /** Finds the coefficients of the indefinite integral of the N-th order polynomial with coeffs
  in "a" and stores them in "ai". The order of the polynomial represented by the coeffs in "ai"
  will be N+1. The constant term in the ai[] polynomial is the arbitrary integration constant
  which may be passed in as parameter "c" - this parameter is optional, it defaults to zero.  */
  static void polyIntegral(T *a, T *ai, int N, T c = T(0));

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
  static void createPolynomialPowers(T *a, int N, T **aPowers, int highestPower);

  /** Let A(x) and B(x) be polynomials represented by their coefficient arrays a[] and b[]
  respectively. This function creates the coefficients of a polynomial C(x), represented by the
  coefficient array c[], that results from composing the polynomials A(x) and B(x), that is: first
  the polynomial A(x) is applied to the value x, and then the polynomial B(x) is applied to the
  result of the first polynomial, such that C(x) = B(A(x)). This nesting or composition of two
  polynomials can itself be seen as a polynomial in its own right. This resulting polynomial has
  an order of cN = aN*bN, where aN and bN are the orders of the a[] and b[] polynomials,
  respectively, so the caller has to make sure that the c[] array has at least a length of
  aN*bN+1. */
  static void composePolynomials(T *a, int aN, T *b, int bN, T *c);

  /** Forms a weighted sum of two polynomials p(x) and q(x) with weights wp and wq respectively and
  stores the coeffficients of the resulting polynomial r(x) in the r-array. The polynomials p(x)
  and q(x) do not need to be of the same order and the resulting polynomial will have an order of
  max(pN, qN). */
  static void weightedSumOfPolynomials(T *p, int pN, T wp, T *q, int qN, T wq, T *r);

  /** Subtracts polynomial q(x) from polynomial p(x) and stores the coeffficients of the resulting
  polynomial in r(x) which is of order max(pN, qN). */
  static void subtractPolynomials(T *p, int pN, T *q, int qN, T *r);

  /** Computes the definite integral of the polynomial "p" where the lower integration limit is
  given by the polynomial "a" and the upper limit is given by the polynomial "b". "p", "a", "b"
  are assumed to be of orders "pN", "aN" and "bN" respectively and the result wil be stored in
  as polynomial "q" which will be of order pN*max(aN, bN). */
  static void integratePolynomialWithPolynomialLimits(T *p, int pN, T *a, int aN, T *b, 
    int bN, T *q);

  /** Given expansion coefficients a[k] of an arbitrary polynomial P(x) with given order in terms
  of a set of N+1 basis polynomials Q0(x), Q1(x), ..., QN(x) such that:
  P(x) = a[0]*Q0[x] + a[1]*Q1[x] + ... + a[N]*QN[x], this function returns the expansion
  coefficients b[k] with respect to another set of basis polynomials R, such that:
  P(x) = b[0]*R0[x] + b[1]*R1[x] + ... + b[N]*RN[x]
  The basis polynomials Q and R are passed as 2-dimensional arrays where the k-th row represents
  the coefficients of the k-th basis polynomial. If R is not a basis, the function will not succeed
  and return false, otherwise true. */
  static bool rsPolynomialBaseChange(T **Q, T *a, T **R, T *b, int order);

  /** Converges to a complex root of a polynomial by means of Laguerre's method using the
  "initialGuess" as first estimate. */
  static std::complex<T> convergeToRootViaLaguerre(std::complex<T> *a, int order,
    std::complex<T> initialGuess = std::complex<T>(0.0, 0.0));

  /** Finds all complex roots of a polynomial by Laguerre's method and returns them in "roots". */
  static void findPolynomialRoots(std::complex<T> *a, int order, std::complex<T> *roots);

  static void findPolynomialRoots(T *a, int order, std::complex<T> *roots);

  /** Same as above but accepts real coefficients. */
  //void findPolynomialRootsInternal(T *a, int order, Complex *roots, bool polish = true);

  //void findPolynomialRootsNew(Complex *a, int order, Complex *roots);

  /** Computes polynomial coefficients from the roots. \todo: get rid of that - replace by function
  below */
  static std::vector<std::complex<T>> getPolynomialCoefficientsFromRoots(
    std::vector<std::complex<T>> roots);

  /** Computes polynomial coefficients from the roots. The roots should be passed in the array "r"
  of length "N", the coefficients will be returned in the array "a" of length "N" + 1. The
  coefficient for the highest power a[N] will be normalized to unity. */
  static void rootsToCoeffs(std::complex<T> *r, std::complex<T> *a, int N);

  /** Similar to rootsToCoeffs(Complex *r, Complex *a, int N), but assumes that the roots are
  either real or occur in complex conjugate pairs. This means that the polynomial has purely real
  coefficients, so the type of the coefficient-array is T instead of Complex. You should use
  this function only if you know in advance that the coefficients will indeed come out as purely
  real */
  static void rootsToCoeffs(std::complex<T> *r, T *a, int N);

  /** Computes the root of the linear equation: \f[ a x + b = 0 \f] which is simply given by
  \f[ x_0 = -\frac{b}{a} \f] */
  static T getRootOfLinearEquation(T a, T b);
    // rename to rootLinear

  /** Computes the two roots of the quadratic equation: \f[ a x^2 + b x + c = 0 \f] which are
  given by: \f[ x_{1,2} = \frac{-b \pm \sqrt{b^2-4ac}}{2a} \f] and stores the result in two-element
  array which is returned. When the qudratic is degenerate (i.e, a == 0), it will fall back to the
  getRootsOfLinearEquation() function, and return a one-element array.  */
  static std::vector<std::complex<T>> getRootsOfQuadraticEquation(T a, T b, T c);

  /** Computes the three roots of the cubic equation: \f[ a x^3 + b x^2 + c x + d = 0 \f] and
  stores the result in the three-element array which is returned. When the cubic is degenerate
  (i.e, a == 0), it will fall back to the getRootsOfQuadraticEquation() function, and return a
  two-element array (or a one-element array, when b is also zero). */
  static std::vector<std::complex<T>> getRootsOfCubicEquation(T a, T b, T c, 
    T d);


  /** Computes the two roots of the quadratic equation: \f[ a_0 + a_1 x + a_2 x^2 = 0 \f] and 
  stores them in r1, r2. When the equation has two distinct real roots, they will be returned in 
  ascending order, i.e. r1 < r2. In case of a (real) double root, we'll have r1 == r2 and when the 
  roots of the equation are actually complex, the outputs will also be equal and contain the real 
  part of the complex conjugate pair. */
  static void rootsQuadraticReal(T a0, T a1, T a2, T* r1, T* r2);

  static void rootsQuadraticComplex(std::complex<T> a0, std::complex<T> a1, std::complex<T> a2, 
    std::complex<T>* x1, std::complex<T>* x2);
  // todo: make optimized version for real coefficients

  /** Discriminant of cubic polynomial \f[ a_0 + a_1 x + a_2 x^2 + a_3 x^3 = 0 \f]. 
  D > 0: 3 distinct real roots, D == 0: 3 real roots, 2 or 3 of which may coincide, 
  D < 0: 1 real root and 2 complex conjugate roots */
  static T discriminant(T a0, T a1, T a2, T a3);

  static void rootsCubicComplex(std::complex<T> a0, std::complex<T> a1, std::complex<T> a2, 
    std::complex<T> a3, std::complex<T>* r1, std::complex<T>* r2, std::complex<T>* r3);
  // under construction - does not yet work

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
  static T getCubicRootNear(T x, T a, T b, T c, T d, T min,
    T max, int maxIterations = 10);

  /** Iteratively improves an initial estimate for the root of the polynomial equation:
  \f[ a[order] x^order + ... + a[1] x + a[0] = 0   \f]
  by means of the Newton-Raphson iteration:
  \f[ x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}       \f]
  the arguments min and max give upper and lower bounds for the root (which will be returned in
  cases where the iteration diverges, which you should avoid in the first place) and maxIterations
  gives the maximum number of iteration steps. */
  static T getRootNear(T x, T *a, int order, T min, T max,
    int maxIterations = 32);

  /** Computes coefficients a[0], a[1], a[2], a[3] for the cubic polynomial that goes through the
  points (x[0], y[0]) and (x[1], y[1]) and has first derivatives of dy[0] and dy[1] at these points
  respectively. */
  static void cubicCoeffsTwoPointsAndDerivatives(T *a, T *x, T *y, T *dy);

  /** Simplified version of
  cubicCoeffsTwoPointsAndDerivatives(T *a, T *x, T *y, T *dy)
  which assumes that x[0] = 0, x[1] = 1 - so we don't need to actually pass the x-array. */
  static void rsCubicCoeffsTwoPointsAndDerivatives(T *a, T *y, T *dy);

  // \todo void cubicCoeffsFourPoints(T *a, T *x, T *y);

  /** Computes coefficients a[0], a[1], a[2], a[3] for the cubic polynomial that goes through the
  points (-1, y[-1]), (0, y[0]), (1, y[1]), (2, y[2]). NOTE, that the y-array is accessed at values
  y[-1]...y[2] - the caller should make sure, these values exist. */
  static void rsCubicCoeffsFourPoints(T *a, T *y);

  /** Allocates and fills an NxN matrix A wher A[i][j] are given by x[i]^j. The caller is
  responsible for deallocation. So it's used like:
  T **A = rsVandermondeMatrix(x, N);
  // ...do stuff with matrix A
  rsDeAllocateSquareArray2D(A, N);  */
  static T** rsVandermondeMatrix(T *x, int N);

  /** Computes coefficients a[0],..., a[N-1] for a polynomial of order N-1 that goes through the N
  data points (x[0], y[0]),...,(x[N-1], y[N-1]). */
  static void rsInterpolatingPolynomial(T *a, T *x, T *y, int N);
    // maybe move to Interplation

  /** Like rsInterpolatingPolynomial(T *a, T *x, T *y, int N), but instead of
  passing an x-array, you should pass a start value x0 and an increment dx and it will use x-values
  given by x0, x0+dx, x0+2*dx,...,x0+(N-1)*dx, so this function assumes equidisant abscissa
  values. */
  static void rsInterpolatingPolynomial(T *a, T x0, T dx, T *y, int N);

  // \todo void quinticCoeffsTwoPointsAndDerivatives(T *a, T *x, T *y, T *dy,
  //                                                 T *d2y);

  /** Fits the quadratic parabola defined by y(x) = a[2]*x^2 + a[1]*x + a[0] to the
  3 points (x[0],y[0]), (x[1],y[1]), (x[2],y[2]). */
  static void fitQuadratic(T *a, T *x, T *y);

  /** Fits the quadratic parabola defined by y(x) = a[2]*x^2 + a[1]*x + a[0] to the
  3 points (0,y[0]), (1,y[1]), (2,y[2]).
  \todo - maybe it's simpler to compute and more convenient to use if we fit a quadratic to
  (-1,y[-1]), (0,y[0]), (1,y[1]). ...but changing this behaviour will affect client code in a way,
  that might not be immediately apparent - we should keep this function and give the other function
  a different name [....] - do this as soon a all code from the legacy codebase is copied over.
  ..OK done - renamed to fitQuadratic_0_1_2 - the other one can be named just fitQuadratic
  */
  static void fitQuadratic_0_1_2(T *a, T *y);

  /** Fits the quartic defined by y(x) = a[4]*x^4 + a3*x^3 + a[2]*x^2 + a[1]*x + a[0] to the
  3 points (0,y[0]), (1,y[1]), (2,y[2]) and also matches the derivatives (slopes)
  y'(0) = s0, y'(2) = s2 */
  static void fitQuarticWithDerivatives(T *a, T *y, T s0, T s2);

  /** Given coefficients of a polynomial a2*x^2 + a1*x + a0, this function determines whether its
  roots are on or inside the unit circle.
   \todo: write and run a unit-test for this function.
  */
  static bool areRootsOnOrInsideUnitCircle(T a0, T a1, T a2);

  // \todo fitPolynomial(T *a, int order, T *x, T *y, int numValues);
  // order+1 == numValues: exact fit
  // order+1 >  numValues: exact fit, some higher coeffs unused -> maybe via recursive call
  // order+1 <  numValues: least-squares fit

  /** Computes polynomial coefficients of a polynomial that is defined recursively by
  w0 * P_n(x) = (w1 + w1x * x) * P_{n-1}(x) + w2 * P_{n-2}(x)
  where n is the order of the polynomial represented by the "a" array, the order of a1 is n-1, the
  order of a2 is n-2. The lengths of the corresponding arrays equals their respective order plus 1.
  w0, w1, w2 are the weighting coeffients of the linear 3-term recurrence relation. The pointer for
  result "a" may point to the same memory location as either of the input argument arrays "a1",
  "a2", so the function may be used in place. */
  static void rsPolynomialRecursion(T *a, T w0, int order, T *a1, T w1, T w1x, T *a2, T w2);

  /** Fills the array with coefficients for a Bessel-polynomial of given order. */
  static void besselPolynomial(T *a, int order);
    // todo: maybe use rsPolynomialRecursion inside

  /** Fills the array with coefficients for a Legendre-polynomial (of the 1st kind) of given
  order. */
  static void legendrePolynomial(T *a, int order);
    // todo: maybe use rsPolynomialRecursion - or maybe get rid of the function 
    // (move to prototypes)

  /** Computes the recursion coefficients (as used in rsPolynomialRecursion) for the Jacobi
  polynomial of order n (n >= 2) with parameters a and b. */
  static void rsJacobiPolynomialRecursionCoeffs(int n, T a, T b, T *w0,
    T *w1, T *w1x, T *w2);

  /** Given the coefficients of the Jacobi polynomials of orders n-1 and n-2 in c1 and c2, this
  function computes the coefficients of the next Jacobi polynomial of order n by using the 3 term
  recurrence relation with parameters a and b. */
  static void rsJacobiPolynomialRecursion(T *c, int n, T *c1, T *c2, T a,
    T b);

  /** Computes the coefficients of the Jacobi polynomials of orders up to maxOrder with parameters
  a and b (usually alpha and beta in formulas). The 2D array "c" will contain the coefficients on
  return. The first index in the c-array runs over indices 0...maxOrder inclusive, so the outer
  dimension should be maxOrder+1. The "c" array may actually be triangular, with the 1st inner
  array of length 1 , the 2nd of length 2, etc. where the last one should have length
  maxOrder+1 (same as the outer). You may also use a square matrix for convenience - then unused
  elements will not be touched in this case. */
  static void rsJacobiPolynomials(T **c, T a, T b, int maxOrder);
    // rename to jacobi

  /** Analog to rsJacobiPolynomialRecursion */
  static void rsLegendrePolynomialRecursion(T *a, int n, T *a1, T *a2);
   // rename to legendreRecursion

  // todo: void chebychevPolynomial(T *a, int order);

  /** Constructs a polynomial p(x) of order 2*N+1 with the following properties: 
  p(0) = 0, p(1) = 1, p'(x) >= 0 for all x (monotonically increasing), p'(1) = maximum possible 
  when monotonicity is assumed. \todo: check if these properties are actually true. Such 
  polynomials are used in Papoulis filters. */
  static void maximumSlopeMonotonicPolynomial(T *a, int N);
    // rename to maxSlopeMonotonic


  // \todo for Halpern filters:
  //void jacobiPolynomial(T *a, int order); // the U-polynomials
  //void maximallyDivergingMonotonicPolynomial(T *a, int order); // the T-polynomial

  // comment this function, maybe use a more efficent algorithm if all
  // poles are simple, (see also Experiments - there's something said about that)
  static void rsPartialFractionExpansion(
    std::complex<T> *numerator, int numeratorOrder,
    std::complex<T> *denominator, int denominatorOrder,
    std::complex<T> *poles, int *multiplicities, int numDistinctPoles,
    std::complex<T> *pfeCoeffs);

  //===============================================================================================
  /** \name Non-static member functions and operators */

  rsPolynomial(int order = 0, bool initWithZeros = true);


  /** Returns the maximum order that this poloynomial may have which is the length of the 
  coefficient array minus one. When there a re trailing zero coefficients, the actual degree of
  the polynomial is lower. */
  int getMaxOrder() const { return (int)coeffs.size()-1; }
  // int getOrder()
  // should take into account traling zeros ..or maybe have a boolean flag 
  // "takeZeroCoeffsIntoAccount" which defaults to false...or maybe it shouldn't have any default
  // value - client code must be explicit

  /** Adds two polynomials */
  rsPolynomial<T> operator+(const rsPolynomial<T>& q)
  { 
    rsPolynomial<T> r(rsMax(getMaxOrder(), q.getMaxOrder()), true);

    // something to do

    return r;
  }



    // maybe we whould take into account the possibility of trailing zero coeffs?
    // maybe have two functions: degreeMax, 



protected:

  std::vector<T> coeffs;

};

// todo: implement a function that determines the number of real roots of a polynomial in an 
// interval by means of Sturmian sequences (see Einführung in die computerorientierte Mathematik 
// mit Sage, p.163ff)

// \todo: implement differentiation and integration, make a class, implement some algos for
// rational functions (partial fraction expansion, differentiation, integration)

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
