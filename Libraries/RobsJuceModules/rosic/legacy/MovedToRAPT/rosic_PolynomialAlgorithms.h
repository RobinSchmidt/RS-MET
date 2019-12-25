#ifndef rosic_PolynomialAlgorithms_h
#define rosic_PolynomialAlgorithms_h

// commented code is obsolete and superseded by templatized versions in RAPT
// todo: replace also the functions involving complex ins/outs - maybe replace higher level code
// that uses them before - for example, in rosic_FilterAnalyzer.h/cpp

namespace rosic
{

  // \todo: implement differentiation and integration, make a class, implement some algos for 
  // rational functions (partial fraction expansion, differentiation, integration)

#ifdef NOT_DEFINED // faux commenting out

  /** Evaluates the polynomial defined by the array of roots "r" at argument "x". */
  Complex evaluatePolynomialWithRoots(Complex x, Complex *r, int numRoots);

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
  derivative in results[1] and so on. @see: evaluatePolynomialAt() */
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
  a polynomial q(x) such that q(x) = p(-x). */
  template <class T>
  void polyCoeffsForNegativeArgument(T *a, T *am, int N);

  /** Finds the coefficients of the derivative of the N-th order polynomial with coeffs in "a" and
  stores them in "ad". The order of the polynomial represented by the coeffs in "ad" will be 
  N-1. */
  template <class T>
  void polyDerivative(T *a, T *ad, int N);

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


  /** Converges to a complex root of a polynomial by means of Laguerre's method using the 
  "initialGuess" as first estimate. */
  Complex convergeToRootViaLaguerre(Complex *a, int order, 
    Complex initialGuess = Complex(0.0, 0.0));

  /** Finds all complex roots of a polynomial by Laguerre's method and returns them in "roots". */
  void findPolynomialRoots(Complex *a, int order, Complex *roots);

  void findPolynomialRoots(double  *a, int order, Complex *roots);

  /** Same as above but accepts real coefficients. */
  //void findPolynomialRootsInternal(double *a, int order, Complex *roots, bool polish = true);

  //void findPolynomialRootsNew(Complex *a, int order, Complex *roots);


  /** Computes polynomial coefficients from the roots. \todo: get rid of that - replace by function 
  below */
  rsArrayTools<Complex> getPolynomialCoefficientsFromRoots(rsArrayTools<Complex> roots);


  /** Computes polynomial coefficients from the roots. The roots should be passed in the array "r" 
  of length "N", the coefficients will be returned in the array "a" of length "N" + 1. The 
  coefficient for the highest power a[N] will be normalized to unity. */
  void rootsToCoeffs(Complex *r, Complex *a, int N);

  /** Similar to rootsToCoeffs(Complex *r, Complex *a, int N), but assumes that the roots are 
  either real or occur in complex conjugate pairs. This means that the polynomial has purely real 
  coefficients, so the type of the coefficient-array is double instead of Complex. You should use 
  this function only if you know in advance that the coefficients will indeed come out as purely 
  real */
  void rootsToCoeffs(Complex *r, double *a, int N);

  /** Computes the root of the linear equation: \f[ a x + b = 0 \f] which is simply given by 
  \f[ x_0 = -\frac{b}{a} \f] */
  double getRootOfLinearEquation(double a, double b);

  /** Computes the two roots of the quadratic equation: \f[ a x^2 + b x + c = 0 \f] which are 
  given by: \f[ x_{1,2} = \frac{-b \pm \sqrt{b^2-4ac}}{2a} \f] and stores the result in two-element 
  array which is returned. When the qudratic is degenerate (i.e, a == 0), it will fall back to the 
  getRootsOfLinearEquation() function, and return a one-element array.  */
  rsArrayTools<Complex> getRootsOfQuadraticEquation(double a, double b, double c);

  /** Computes the three roots of the cubic equation: \f[ a x^3 + b x^2 + c x + d = 0 \f] and 
  stores the result in the three-element array which is returned. When the cubic is degenerate 
  (i.e, a == 0), it will fall back to the getRootsOfQuadraticEquation() function, and return a 
  two-element array (or a one-element array, when b is also zero). */
  rsArrayTools<Complex> getRootsOfCubicEquation(double a, double b, double c, double d);

#endif

#ifdef NOT_DEFINED // faux commenting out

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
  double getCubicRootNear(double x, double a, double b, double c, double d, double min, double max, 
    int maxIterations = 10);

  /** Iteratively improves an initial estimate for the root of the polynomial equation:
  \f[ a[order] x^order + ... + a[1] x + a[0] = 0   \f]
  by means of the Newton-Raphson iteration:
  \f[ x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}       \f]
  the arguments min and max give upper and lower bounds for the root (which will be returned in 
  cases where the iteration diverges, which you should avoid in the first place) and maxIterations 
  gives the maximum number of iteration steps. */
  double getRootNear(double x, double *a, int order, double min, double max, 
    int maxIterations = 32);

  /** Computes coefficients a[0], a[1], a[2], a[3] for the cubic polynomial that goes through the 
  points (x[0], y[0]) and (x[1], y[1]) and has first derivatives of dy[0] and dy[1] at these points 
  respectively. */
  void cubicCoeffsTwoPointsAndDerivatives(double *a, double *x, double *y, double *dy);

  // \todo void cubicCoeffsFourPoints(double *a, double *x, double *y);
  // \todo void quinticCoeffsTwoPointsAndDerivatives(double *a, double *x, double *y, double *dy, double *d2y);

  // \todo fitPolynomial(double *a, int order, double *x, double *y, int numValues); 
  // order+1 == numValues: exact fit 
  // order+1 >  numValues: exact fit, some higher coeffs unused -> maybe via recursive call
  // order+1 <  numValues: least-squares fit


  /** Fills the array with coefficients for a Bessel-polynomial of given order. */
  void besselPolynomial(double *a, int order);

  /** Fills the array with coefficients for a Legendre-polynomial (of the 1st kind) of given 
  order. */
  void legendrePolynomial(double *a, int order);

  /** Constructs a polynomial p(x) of order 2*N+1 with the following properties: 
  p(0) = 0, p(1) = 1, p'(x) >= 0 for all x (monotonically increasing), p'(1) = maximum possible 
  when monotonicity is assumed. \todo: check if these properties are actually true. Such 
  polynomials are used in Papoulis filters. */
  void maximumSlopeMonotonicPolynomial(double *a, int N);

  // \todo for Halpern filters (see Paarmann, page 255-259):
  //void jacobiPolynomial(double *a, int order); // the U-polynomials
  //void maximallyDivergingMonotonicPolynomial(double *a, int order); // the T-polynomial

#endif


  //-----------------------------------------------------------------------------------------------
  // template function definitions:
  /*
  template <class T>
  T evaluatePolynomialAt(T x, T *a, int order)
  {
    if( order < 0 )
      return T(0);
    T y = a[order];
    for(int i=order-1; i>=0; i--)
      y = y*x + a[i];
    return y;
  }

  template <class T>
  void evaluatePolynomialAndDerivativeAt(T x, T *a, int order, T *y, T *yd)
  {
    *y  = a[order];
    *yd = 0.0;
    for(int i=order-1; i>=0; i--)
    {
      *yd = *yd * x + *y;
      *y  = *y  * x + a[i];
    }
  }

  template <class T>
  void evaluatePolynomialAndDerivativesAt(T x, T *a, int order, T *results, int numDerivatives)
  {
    int nnd, j, i;  // rename nnd variable, change declaration order
    T   constant(1);
    results[0] = a[order];
    for(j=1; j<=numDerivatives; j++) 
      results[j] = T(0);
    for(i=order-1; i>=0; i--) 
    {
      nnd = (numDerivatives < (order-i) ? numDerivatives : order-i); // use if-statement

      for(j=nnd; j>=1; j--)
        results[j] = results[j]*x + results[j-1];
      results[0] = results[0]*x + a[i];
    }
    for(i=2; i<=numDerivatives; i++) 
    {
      constant   *= i;
      results[i] *= constant;
    }
  }

  template <class T>
  void multiplyPolynomials(T *a, int aOrder, T *b, int bOrder, T *result)
  {
    convolve(a, aOrder+1, b, bOrder+1, result);
  }

  template <class T>
  void dividePolynomials(T *dividend, int dividendOrder, T *divisor, int divisorOrder, T *quotient, 
    T *remainder)
  {
    int k, j;
    for(j=0; j<=dividendOrder; j++) 
    {
      remainder[j] = dividend[j];
      quotient[j]  = T(0);
    }
    for(k=dividendOrder-divisorOrder; k>=0; k--) 
    {
      quotient[k] = remainder[divisorOrder+k] / divisor[divisorOrder];
      for(j=divisorOrder+k-1; j>=k; j--) 
        remainder[j] -= quotient[k] * divisor[j-k];
    }
    for(j=divisorOrder; j<=dividendOrder; j++) 
      remainder[j] = T(0);
  }

  template <class T>
  void dividePolynomialByMonomialInPlace(T *dividendAndResult, int dividendOrder, T x0, 
    T *remainder)
  {
    *remainder                       = dividendAndResult[dividendOrder];
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
    // we invert the sign for all coefficients that multiply odd powers of x:
    double s = 1.0; // the sign
    for(int n = 0; n <= N; n++)
    {
      am[n]  = s*a[n];
      s     *= -1.0;
    }
  }

  template <class T>
  void polyDerivative(T *a, T *ad, int N)
  {
    for(int n = 1; n <= N; n++)
      ad[n-1] = n * a[n];
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
    copy(a, aPowers[1], N+1);
    for(int k = 2; k <= highestPower; k++)
      convolve(aPowers[k-1], (k-1)*N+1, a, N+1, aPowers[k]);
  }

  template <class T>
  void composePolynomials(T *a, int aN, T *b, int bN, T *c)
  {
    int cN = aN*bN;
    T *an  = new T[cN+1];  // array for the successive powers of a[]
    an[0]  = T(1);         // initialize to a[]^0

    // accumulation:
    fillWithZeros(c, cN+1);
    c[0] = b[0];
    int K = 1;
    for(int n = 1; n <= bN; n++)
    {
      convolveInPlace(an, K, a, aN+1);
      K += aN;
      for(int k = 0; k < K; k++)
        c[k] += b[n] * an[k];
    }

    delete[] an;
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

    RAPT::rsPolynomial<double>::polyIntegral(p, P, pN);
    composePolynomials(a, aN, P, PN, A);
    composePolynomials(b, bN, P, PN, B);
    subtractPolynomials(B, BN, A, AN, q);

    delete[] P;
    delete[] A;
    delete[] B;
  }
  */
}

#endif
