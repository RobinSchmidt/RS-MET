#ifndef RS_POLYNOMIALALGORITHMS_H
#define RS_POLYNOMIALALGORITHMS_H

namespace RSLib
{

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

  /** Evaluates the polynomial defined by the array of roots "r" at argument "x". */
  RSLib_API rsComplexDbl evaluatePolynomialWithRoots(rsComplexDbl x, rsComplexDbl *r,
    int numRoots);

  /** Evaluates the polynomial defined by the array of coefficients "a" at argument "x".  The array
  of coefficients must be of length order+1 and is interpreted as follows: a[0] is taken to be the
  constant term, a[1] is the multiplier for x^1, a[2] the multiplier for x^2 and so on until
  a[order] which is the multiplier for a^order. */
  template <class T>
  T evaluatePolynomialAt(T x, T *a, int order);

  /** Evaluates the polynomial defined by the array of coefficients 'a' and its first derivative at
  argument 'x'. The value of the polynomial will be stored in y and the value of the derivative
  will be stored in yd. @see: evaluatePolynomialAt() */
  template <class T>
  void evaluatePolynomialAndDerivativeAt(T x, T *a, int order, T *y, T *yd);

  /** Evaluates the polynomial defined by the array of coefficients 'a' and a given number of
  derivatives at argument 'x'. The value of the polynomial will be stored in results[0], the 1st
  derivative in results[1] and so on. numDerivatives should be <= 31. 
  @see: evaluatePolynomialAt() */
  template <class T>
  void evaluatePolynomialAndDerivativesAt(T x, T *a, int order, T *results, int numDerivatives);

  /** Multiplies the polynomials represented by the coefficient vectors 'a' and 'b' and stores the
  resulting coefficients in 'result'. The resulting polynom will be or order aOrder+bOrder and the
  coefficient vector should have allocated space for
  (aOrder+1)+(bOrder+1)-1 = aOrder+bOrder+1 = aLength+bLength-1 elements. */
  template <class T>
  void multiplyPolynomials(T *a, int aOrder, T *b, int bOrder, T *result);

  /** Divides the polynomials represented by the coefficient arrays 'dividend' and 'divisor' and
  stores the resulting coefficients for the quotient and remainder in the respective output arrays.
  ?? The resulting quotient polynom will be of order dividendOrder-divisorOrder and the
  remainder polynom will be at most of order ??
  However, the output arrays must have the same length as the dividend, where
  remainder[divisorOrder...dividendOrder] and quotient[dividendOrder-divisorOrder+1...dividendOrder]
  will be filled with zeros.  ...\todo check this comment */
  template <class T>
  void dividePolynomials(T *dividend, int dividendOrder, T *divisor, int divisorOrder, T *quotient,
    T *remainder);

  /** Divides the dividend by the monomial factor (x-x0) and stores the result in the same array
  again. The remainder (which is just a numerical value) will be stored in 'remainder'. */
  template <class T>
  void dividePolynomialByMonomialInPlace(T *dividendAndResult, int dividendOrder, T x0, T *remainder);

  /** Given an array of polynomial coefficients "a" such that
  p(x) = a[0]*x^0 + a[1]*x^1 + ... + a[N]*x^N, this function returns (in "am") the coefficients for
  a polynomial q(x) such that q(x) = p(-x). This amounts to sign-inverting all coefficients which
  multiply odd powers of x. */
  template <class T>
  void polyCoeffsForNegativeArgument(T *a, T *am, int N);

  /** Given an array of polynomial coefficients "a" such that
  p(x) = a[0]*x^0 + a[1]*x^1 + ... + a[N]*x^N, this function returns (in "aShifted") the coefficients 
  for a polynomial q(x) such that q(x) = p(x-x0). */
  template <class T>
  void polyCoeffsForShiftedArgument(T *a, T *aShifted, int N, T x0);

  /** Finds the coefficients of the derivative of the N-th order polynomial with coeffs in "a" and
  stores them in "ad". The order of the polynomial represented by the coeffs in "ad" will be
  N-1. The "ad" array may point to the same array as the "a" array, so you can use the same array
  for input and output. */
  template <class T>
  void polyDerivative(T *a, T *ad, int N);

  /** Given a polynomial p(x) via its coefficient array a, this function computes the coefficients 
  of a polynomial q(x) = p(x+h) - p(x) if direction = 1, or q(x) = p(x) - p(x-h) if direction = -1.
  This represents the first forward or backward difference polynomial in finite difference 
  calculus and is analog to the derivative in infinitesimal calculus. As in infinitesimal calculus, 
  the resulting polynomial is one order less than p(x) itself. Scaled by 1/h, this can be seen as an
  approximation to the derivative using stepsize h. 
  \todo provide a finite central difference q(x) = p(x+h/2) - p(x-h/2) when direction = 0
   ...could this be just the average between forward and backward difference? ...research! */
  template <class T>
  void polyFiniteDifference(T *a, T *ad, int N, int direction = 1, T h = 1);

  /** Finds the coefficients of the indefinite integral of the N-th order polynomial with coeffs
  in "a" and stores them in "ai". The order of the polynomial represented by the coeffs in "ai"
  will be N+1. The constant term in the ai[] polynomial is the arbitrary integration constant
  which may be passed in as parameter "c" - this parameter is optional, it defaults to zero.  */
  template <class T>
  void polyIntegral(T *a, T *ai, int N, T c = T(0));

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
  template<class T>
  void createPolynomialPowers(T *a, int N, T **aPowers, int highestPower);

  /** Let A(x) and B(x) be polynomials represented by their coefficient arrays a[] and b[]
  respectively. This function creates the coefficients of a polynomial C(x), represented by the
  coefficient array c[], that results from composing the polynomials A(x) and B(x), that is: first
  the polynomial A(x) is applied to the value x, and then the polynomial B(x) is applied to the
  result of the first polynomial, such that C(x) = B(A(x)). This nesting or composition of two
  polynomials can itself be seen as a polynomial in its own right. This resulting polynomial has
  an order of cN = aN*bN, where aN and bN are the orders of the a[] and b[] polynomials,
  respectively, so the caller has to make sure that the c[] array has at least a length of
  aN*bN+1. */
  template<class T>
  void composePolynomials(T *a, int aN, T *b, int bN, T *c);

  /** Forms a weighted sum of two polynomials p(x) and q(x) with weights wp and wq respectively and
  stores the coeffficients of the resulting polynomial r(x) in the r-array. The polynomials p(x)
  and q(x) do not need to be of the same order and the resulting polynomial will have an order of
  max(pN, qN). */
  template<class T>
  void weightedSumOfPolynomials(T *p, int pN, T wp, T *q, int qN, T wq, T *r);

  /** Subtracts polynomial q(x) from polynomial p(x) and stores the coeffficients of the resulting
  polynomial in r(x) which is of order max(pN, qN). */
  template<class T>
  void subtractPolynomials(T *p, int pN, T *q, int qN, T *r);

  /** Computes the definite integral of the polynomial "p" where the lower integration limit is
  given by the polynomial "a" and the upper limit is given by the polynomial "b". "p", "a", "b"
  are assumed to be of orders "pN", "aN" and "bN" respectively and the result wil be stored in
  as polynomial "q" which will be of order pN*max(aN, bN). */
  template<class T>
  void integratePolynomialWithPolynomialLimits(T *p, int pN, T *a, int aN, T *b, int bN, T *q);

  /** Given expansion coefficients a[k] of an arbitrary polynomial P(x) with given order in terms 
  of a set of N+1 basis polynomials Q0(x), Q1(x), ..., QN(x) such that:
  P(x) = a[0]*Q0[x] + a[1]*Q1[x] + ... + a[N]*QN[x], this function returns the expansion 
  coefficients b[k] with respect to another set of basis polynomials R, such that:
  P(x) = b[0]*R0[x] + b[1]*R1[x] + ... + b[N]*RN[x]
  The basis polynomials Q and R are passed as 2-dimensional arrays where the k-th row represents 
  the coefficients of the k-th basis polynomial. If R is not a basis, the function will not succeed
  and return false, otherwise true. */
  template<class T>
  bool rsPolynomialBaseChange(T **Q, T *a, T **R, T *b, int order);

  /** Converges to a complex root of a polynomial by means of Laguerre's method using the
  "initialGuess" as first estimate. */
  RSLib_API rsComplexDbl convergeToRootViaLaguerre(rsComplexDbl *a, int order,
    rsComplexDbl initialGuess = rsComplexDbl(0.0, 0.0));

  /** Finds all complex roots of a polynomial by Laguerre's method and returns them in "roots". */
  RSLib_API void findPolynomialRoots(rsComplexDbl *a, int order, rsComplexDbl *roots);

  RSLib_API void findPolynomialRoots(double  *a, int order, rsComplexDbl *roots);


  /** Same as above but accepts real coefficients. */
  //void findPolynomialRootsInternal(double *a, int order, Complex *roots, bool polish = true);

  //void findPolynomialRootsNew(Complex *a, int order, Complex *roots);


  /** Computes polynomial coefficients from the roots. \todo: get rid of that - replace by function
  below */
  RSLib_API rsArray<rsComplexDbl> getPolynomialCoefficientsFromRoots(rsArray<rsComplexDbl> roots);


  /** Computes polynomial coefficients from the roots. The roots should be passed in the array "r"
  of length "N", the coefficients will be returned in the array "a" of length "N" + 1. The
  coefficient for the highest power a[N] will be normalized to unity. */
  RSLib_API void rootsToCoeffs(rsComplexDbl *r, rsComplexDbl *a, int N);

  /** Similar to rootsToCoeffs(Complex *r, Complex *a, int N), but assumes that the roots are
  either real or occur in complex conjugate pairs. This means that the polynomial has purely real
  coefficients, so the type of the coefficient-array is double instead of Complex. You should use
  this function only if you know in advance that the coefficients will indeed come out as purely
  real */
  RSLib_API void rootsToCoeffs(rsComplexDbl *r, double *a, int N);

  /** Computes the root of the linear equation: \f[ a x + b = 0 \f] which is simply given by
  \f[ x_0 = -\frac{b}{a} \f] */
  RSLib_API double getRootOfLinearEquation(double a, double b);

  /** Computes the two roots of the quadratic equation: \f[ a x^2 + b x + c = 0 \f] which are
  given by: \f[ x_{1,2} = \frac{-b \pm \sqrt{b^2-4ac}}{2a} \f] and stores the result in two-element
  array which is returned. When the qudratic is degenerate (i.e, a == 0), it will fall back to the
  getRootsOfLinearEquation() function, and return a one-element array.  */
  RSLib_API rsArray<rsComplexDbl> getRootsOfQuadraticEquation(double a, double b, double c);

  /** Computes the three roots of the cubic equation: \f[ a x^3 + b x^2 + c x + d = 0 \f] and
  stores the result in the three-element array which is returned. When the cubic is degenerate
  (i.e, a == 0), it will fall back to the getRootsOfQuadraticEquation() function, and return a
  two-element array (or a one-element array, when b is also zero). */
  rsArray<rsComplexDbl> getRootsOfCubicEquation(double a, double b, double c, double d);

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
  RSLib_API double getCubicRootNear(double x, double a, double b, double c, double d, double min,
    double max, int maxIterations = 10);

  /** Iteratively improves an initial estimate for the root of the polynomial equation:
  \f[ a[order] x^order + ... + a[1] x + a[0] = 0   \f]
  by means of the Newton-Raphson iteration:
  \f[ x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}       \f]
  the arguments min and max give upper and lower bounds for the root (which will be returned in
  cases where the iteration diverges, which you should avoid in the first place) and maxIterations
  gives the maximum number of iteration steps. */
  RSLib_API double getRootNear(double x, double *a, int order, double min, double max,
    int maxIterations = 32);

  /** Computes coefficients a[0], a[1], a[2], a[3] for the cubic polynomial that goes through the
  points (x[0], y[0]) and (x[1], y[1]) and has first derivatives of dy[0] and dy[1] at these points
  respectively. */
  RSLib_API void cubicCoeffsTwoPointsAndDerivatives(double *a, double *x, double *y, double *dy);

  /** Simplified version of 
  cubicCoeffsTwoPointsAndDerivatives(double *a, double *x, double *y, double *dy)
  which assumes that x[0] = 0, x[1] = 1 - so we don't need to actually pass the x-array. */
  RSLib_API void rsCubicCoeffsTwoPointsAndDerivatives(double *a, double *y, double *dy);

  // \todo void cubicCoeffsFourPoints(double *a, double *x, double *y);

  /** Computes coefficients a[0], a[1], a[2], a[3] for the cubic polynomial that goes through the
  points (-1, y[-1]), (0, y[0]), (1, y[1]), (2, y[2]). NOTE, that the y-array is accessed at values
  y[-1]...y[2] - the caller should make sure, these values exist. */
  RSLib_API void rsCubicCoeffsFourPoints(double *a, double *y);

  /** Allocates and fills an NxN matrix A wher A[i][j] are given by x[i]^j. The caller is 
  responsible for deallocation. So it's used like:
  double **A = rsVandermondeMatrix(x, N);
  // ...do stuff with matrix A
  rsDeAllocateSquareArray2D(A, N);  */
  RSLib_API double** rsVandermondeMatrix(double *x, int N);

  /** Computes coefficients a[0],..., a[N-1] for a polynomial of order N-1 that goes through the N
  data points (x[0], y[0]),...,(x[N-1], y[N-1]). */
  RSLib_API void rsInterpolatingPolynomial(double *a, double *x, double *y, int N);
    // maybe move to Interplation

  /** Like rsInterpolatingPolynomial(double *a, double *x, double *y, int N), but instead of 
  passing an x-array, you should pass a start value x0 and an increment dx and it will use x-values
  given by x0, x0+dx, x0+2*dx,...,x0+(N-1)*dx, so this function assumes equidisant abscissa 
  values. */
  RSLib_API void rsInterpolatingPolynomial(double *a, double x0, double dx, double *y, int N);



  // \todo void quinticCoeffsTwoPointsAndDerivatives(double *a, double *x, double *y, double *dy,
  //                                                 double *d2y);

  /** Fits the quadratic parabola defined by y(x) = a[2]*x^2 + a[1]*x + a[0] to the
  3 points (x[0],y[0]), (x[1],y[1]), (x[2],y[2]). */
  RSLib_API void fitQuadratic(double *a, double *x, double *y);

  /** Fits the quadratic parabola defined by y(x) = a[2]*x^2 + a[1]*x + a[0] to the
  3 points (0,y[0]), (1,y[1]), (2,y[2]). 
  \todo - maybe it's simpler to compute and more convenient to use if we fit a quadratic to
  (-1,y[-1]), (0,y[0]), (1,y[1]). ...but changing this behaviour will affect client code in a way, 
  that might not be immediately apparent - we should keep this function and give the other function
  a different name [....] - do this as soon a all code from the legacy codebase is copied over. 
  ..OK done - renamed to fitQuadratic_0_1_2 - the other one can be named just fitQuadratic
  */
  RSLib_API void fitQuadratic_0_1_2(double *a, double *y);

  /** Fits the quartic defined by y(x) = a[4]*x^4 + a3*x^3 + a[2]*x^2 + a[1]*x + a[0] to the
  3 points (0,y[0]), (1,y[1]), (2,y[2]) and also matches the derivatives (slopes)
  y'(0) = s0, y'(2) = s2 */
  RSLib_API void fitQuarticWithDerivatives(double *a, double *y, double s0, double s2);




  /** Given coefficients of a polynomial a2*x^2 + a1*x + a0, this function determines whether its
  roots are on or inside the unit circle.
   \todo: write and run a unit-test for this function.
  */
  RSLib_API bool areRootsOnOrInsideUnitCircle(double a0, double a1, double a2);

  // \todo fitPolynomial(double *a, int order, double *x, double *y, int numValues);
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
  template <class T>
  void rsPolynomialRecursion(T *a, T w0, int order, T *a1, T w1, T w1x, T *a2, T w2);

  /** Fills the array with coefficients for a Bessel-polynomial of given order. */
  RSLib_API void besselPolynomial(double *a, int order);
    // todo: maybe use rsPolynomialRecursion inside

  /** Fills the array with coefficients for a Legendre-polynomial (of the 1st kind) of given
  order. */
  RSLib_API void legendrePolynomial(double *a, int order);
    // todo: maybe use rsPolynomialRecursion - or maybe get rid of the function 
    // (move to prototypes)

  /** Computes the recursion coefficients (as used in rsPolynomialRecursion) for the Jacobi 
  polynomial of order n (n >= 2) with parameters a and b. */
  RSLib_API void rsJacobiPolynomialRecursionCoeffs(int n, double a, double b, double *w0, 
    double *w1, double *w1x, double *w2);

  /** Given the coefficients of the Jacobi polynomials of orders n-1 and n-2 in c1 and c2, this
  function computes the coefficients of the next Jacobi polynomial of order n by using the 3 term 
  recurrence relation with parameters a and b. */
  RSLib_API void rsJacobiPolynomialRecursion(double *c, int n, double *c1, double *c2, double a, 
    double b);

  /** Computes the coefficients of the Jacobi polynomials of orders up to maxOrder with parameters 
  a and b (usually alpha and beta in formulas). The 2D array "c" will contain the coefficients on 
  return. The first index in the c-array runs over indices 0...maxOrder inclusive, so the outer 
  dimension should be maxOrder+1. The "c" array may actually be triangular, with the 1st inner 
  array of length 1 , the 2nd of length 2, etc. where the last one should have length 
  maxOrder+1 (same as the outer). You may also use a square matrix for convenience - then unused 
  elements will not be touched in this case. */
  RSLib_API void rsJacobiPolynomials(double **c, double a, double b, int maxOrder);

  /** Analog to rsJacobiPolynomialRecursion */
  RSLib_API void rsLegendrePolynomialRecursion(double *a, int n, double *a1, double *a2);


  // todo: void chebychevPolynomial(double *a, int order);

  // \todo for Halpern filters:
  //void jacobiPolynomial(double *a, int order); // the U-polynomials
  //void maximallyDivergingMonotonicPolynomial(double *a, int order); // the T-polynomial

  // comment this function, maybe use a more efficent algorithm if all
  // poles are simple, (see also Experiments - there's something said about that)
  RSLib_API void rsPartialFractionExpansion(
    rsComplexDbl *numerator,   int numeratorOrder,
    rsComplexDbl *denominator, int denominatorOrder, 
    rsComplexDbl *poles, int *multiplicities, int numDistinctPoles,
    rsComplexDbl *pfeCoeffs);

  //-----------------------------------------------------------------------------------------------
  // template function definitions:

  template <class T>
  T evaluatePolynomialAt(T x, T *a, int order)
  {
    if( order < 0 )
      return T(0);
    T y = a[order];
    for(int i = order-1; i >= 0; i--)
      y = y*x + a[i];
    return y;
  }

  template <class T>
  void evaluatePolynomialAndDerivativeAt(T x, T *a, int order, T *y, T *yd)
  {
    *y  = a[order];
    *yd = 0.0;
    for(int i = order-1; i >= 0; i--)
    {
      *yd = *yd * x + *y;
      *y  = *y  * x + a[i];
    }
  }

  template <class T>
  void evaluatePolynomialAndDerivativesAt(T x, T *a, int order, T *results, int numDerivatives)
  {
    results[0] = a[order];
    rsFillWithZeros(&results[1], numDerivatives);
    for(int i = order-1; i >= 0; i--)
    {
      int n = rsMin(numDerivatives, order-1);
      for(int j = n; j >= 1; j--)
        results[j] = results[j]*x + results[j-1];
      results[0] = results[0]*x + a[i];
    }
    rsMultiply(&results[2], &rsFactorials[2], &results[2], numDerivatives-1);
  }

  template <class T>
  void multiplyPolynomials(T *a, int aOrder, T *b, int bOrder, T *result)
  {
    rsConvolve(a, aOrder+1, b, bOrder+1, result);
  }

  template <class T>
  void dividePolynomials(T *p, int pOrder, T *d, int dOrder, T *q, T *r)
  {
    rsCopyBuffer(p, r, pOrder+1); // init remainder with p
    rsFillWithZeros(q, pOrder+1); // init quotient with zeros
    for(int k = pOrder-dOrder; k >= 0; k--)
    {
      q[k] = r[dOrder+k] / d[dOrder];
      for(int j = dOrder+k-1; j >= k; j--)
        r[j] -= q[k] * d[j-k];
    }
    rsFillWithZeros(&r[dOrder], pOrder-dOrder+1);
  }

  template <class T>
  void dividePolynomialByMonomialInPlace(T *dividendAndResult, int dividendOrder, T x0,
    T *remainder)
  {
    *remainder = dividendAndResult[dividendOrder];
    dividendAndResult[dividendOrder] = T(0);
    for(int i=dividendOrder-1; i>=0; i--)
    {
      T swap               = dividendAndResult[i];
      dividendAndResult[i] = *remainder;
      *remainder           = swap + *remainder*x0;
    }
  }



  template <class T>
  void polyCoeffsForNegativeArgument(T *a, T *am, int N)
  {
    double s = 1.0;
    for(int n = 0; n <= N; n++)
    {
      am[n]  = s*a[n];
      s     *= -1.0;
    }
  }

  // todo: polyCoeffsForScaledArgument: aScaled[n] = a[n] * scaler^n - when the scaler equals -1, 
  // it reduces to polyCoeffsForNegativeArgument - this function is superfluous then

  template <class T>
  void polyCoeffsForShiftedArgument(T *a, T *as, int N, T x0)
  {
    int numLines = N+1;
    int length   = (numLines*(numLines+1))/2;
    rsUint32 *pt = new rsUint32[length];
    rsCreatePascalTriangle(pt, numLines);
    T *x0n = new T[N+1];  // +- x0^n
    x0n[0] = 1.0;
    for(int n = 1; n <= N; n++)
      x0n[n] = -x0 * x0n[n-1];
    for(int n = 0; n <= N; n++)
    {
      as[n] = 0.0;
      for(int k = n; k <= N; k++)
        as[n] += rsPascalTriangle(pt, k, k-n) * x0n[k-n] * a[k];
    }
    delete[] pt;
    delete[] x0n;
  }

  template <class T>
  void polyDerivative(T *a, T *ad, int N)
  {
    for(int n = 1; n <= N; n++)
      ad[n-1] = n * a[n];
  }

  template <class T>
  void polyFiniteDifference(T *a, T *ad, int N, int direction, T h)
  {
    // (possibly alternating) powers of the stepsize h:
    T *hk = new T[N+1];
    T hs  = direction*h;
    hk[0] = T(1);
    for(int k = 1; k <= N; k++)
      hk[k] = hk[k-1] * hs;

    // binomial coefficients:
    int numCoeffs    = N+1;
    int triangleSize = (numCoeffs*(numCoeffs+1))/2;
    rsUint32 *binomCoeffs = new rsUint32[triangleSize];
    rsCreatePascalTriangle(binomCoeffs, numCoeffs);

    // actual coefficient computation for ad:
    rsFillWithZeros(ad, N);
    for(int n = 0; n <= N; n++)
    {
      for(int k = 1; k <= n; k++)
        ad[n-k] += a[n] * rsPascalTriangle(binomCoeffs, n, k) * hk[k];
    }
    if( direction == -1 )
      rsScale(ad, N, -1);

    delete[] hk;
    delete[] binomCoeffs;
  }

  template <class T>
  void polyIntegral(T *a, T *ai, int N, T c)
  {
    for(int n = N+1; n >= 1; n--)
      ai[n] = a[n-1] / n;
    ai[0] = c;
  }

  template <class T>
  void createPolynomialPowers(T *a, int N, T **aPowers, int highestPower)
  {
    aPowers[0][0] = 1;
    if( highestPower < 1 )
      return;
    rsCopyBuffer(a, aPowers[1], N+1);
    for(int k = 2; k <= highestPower; k++)
      rsConvolve(aPowers[k-1], (k-1)*N+1, a, N+1, aPowers[k]);
  }

  template <class T>
  void composePolynomials(T *a, int aN, T *b, int bN, T *c)
  {
    int cN = aN*bN;
    T *an  = new T[cN+1];  // array for the successive powers of a[]
    an[0]  = T(1);         // initialize to a[]^0

    // accumulation:
    rsFillWithZeros(c, cN+1);
    c[0] = b[0];
    int K = 1;
    for(int n = 1; n <= bN; n++)
    {
      rsConvolveInPlace(an, K, a, aN+1);
      K += aN;
      for(int k = 0; k < K; k++)
        c[k] += b[n] * an[k];
    }

    delete[] an;
  }

  template <class T>
  void rsPolynomialRecursion(T *a, T w0, int order, T *a1, T w1, T w1x, T *a2, T w2)
  {
    rsAssert( order >= 2 );
    int n = order;
    a[n] = (w1x*a1[n-1]) / w0;
    n--;
    a[n] = (w1*a1[n] + w1x*a1[n-1]) / w0;
    for(n = n-1; n > 0; n--)
      a[n] = (w1*a1[n] + w1x*a1[n-1] + w2*a2[n]) / w0;
    a[0] = (w1*a1[0] + w2*a2[0]) / w0;
  }

  template<class T>
  void weightedSumOfPolynomials(T *p, int pN, T wp, T *q, int qN, T wq, T *r)
  {
    int i;
    if( pN >= qN )
    {
      for(i = 0; i <= qN; i++)
        r[i] = wp*p[i] + wq*q[i];
      for(i = qN+1; i <= pN; i++)
        r[i] = wp*p[i];
    }
    else
    {
      for(i = 0; i <= pN; i++)
        r[i] = wp*p[i] + wq*q[i];
      for(i = pN+1; i <= qN; i++)
        r[i] = wq*q[i];
    }
  }

  template<class T>
  void subtractPolynomials(T *p, int pN, T *q, int qN, T *r)
  {
    weightedSumOfPolynomials(p, pN, T(1), q, qN, T(-1), r);
  }

  template<class T>
  void integratePolynomialWithPolynomialLimits(T *p, int pN, T *a, int aN, T *b, int bN, T *q)
  {
    int PN = pN+1;
    int AN = aN*PN;
    int BN = bN*PN;

    T *P = new T[PN+1];
    T *A = new T[AN+1];
    T *B = new T[BN+1];

    polyIntegral(p, P, pN);               // P(x) is the antiderivative of p(x)
    composePolynomials(a, aN, P, PN, A);  // A(x) = P(a(x))
    composePolynomials(b, bN, P, PN, B);  // B(x) = P(b(x)) 
    subtractPolynomials(B, BN, A, AN, q); // q(x) = B(x) - A(x)

    delete[] P;
    delete[] A;
    delete[] B;
  }

  template<class T>
  bool rsPolynomialBaseChange(T **Q, T *a, T **R, T *b, int order)
  {
    return rsChangeOfBasisRowWise(Q, R, a, b, order+1);
  }

}

#endif
