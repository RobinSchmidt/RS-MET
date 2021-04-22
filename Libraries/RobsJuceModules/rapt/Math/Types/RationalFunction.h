#pragma once

//namespace RAPT
//{

/** A class for representing rational functions R(x) = P(x) / Q(x) where P and Q are both 
polynomials in x. Rational functions are important in signal processing because the transfer 
functions of linear time invariant systems (a.k.a. filters) are of that type. The class provides 
facilities for performing arithmetic with rational functions (including composition), evaluation,
partial fraction expansion, etc. 

Note that instances of this class are not suitable for realtime code - there's a lot of memory 
allocation going on in the operators (std::vectors are created, resized, assigned, etc.) and the
algorithms are not optimized either. ToDo: factor out low-level implementations that operate on
pre-allocated raw arrays.

to verify:
I think, the rational functions over a field (like the real or complex numbers) form themselves a 
field - as opposed to the polynomials, which form only a ring. That's why we can always do a 
division of two rational functions without having to worry about a remainder (which may occur in 
polynomial division). */

template<class T>
class rsRationalFunction
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  rsRationalFunction() {}

  rsRationalFunction(
    const std::vector<T>& numeratorCoeffs, 
    const std::vector<T>& denominatorCoeffs,
    const T& reductionTolerance = T(0) )
    : num(numeratorCoeffs), den(denominatorCoeffs), tol(reductionTolerance)
  {

  }

  rsRationalFunction(const T& number) : num(number), den(T(1))
  { 

  }


  /*
  rsRationalFunction(const T& number, const T& reductionTolerance = T(0)) 
    : num(number), den(T(1)), tol(reductionTolerance)
  { 

  }
  */
  // get rid of this constructor - when calling it like
  // r = rsRationalFunction<double>({1}, {1});  we don't get the constant "1" function but one with
  // a tolerance of 1 - the pattern matching does not work as one might want or expect
  // ...but damn - we need a constructor that takes an int - it gets called, for example, in 
  // rsMatrix::getDiagonalProduct - but this construtor cannot take an optional tolerance - that
  // messes up the pattern matching that decides which constructor is called (and which implicit 
  // conversions are made)

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  // setNumerator, setDenominator // for polynomial and std::vector and maybe plain arrays

  bool reduce(T tol);


  void scaleNumerator(const T& scaler) { num.scale(scaler); }

  void scaleDenominator(const T& scaler) { den.scale(scaler); }

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  rsPolynomial<T> getNumerator()   const { return num; }
  rsPolynomial<T> getDenominator() const { return den; }

  int getNumeratorDegree()   const { return num.getDegree(); }
  int getDenominatorDegree() const { return den.getDegree(); }


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Returns a new rational function that is the negative of this one. */
  rsRationalFunction<T> operator-() const 
  {
    rsRationalFunction<T> r = *this;
    r.num.negate();
    return r;
  }

  /** Adds two rational functions. */
  rsRationalFunction<T> operator+(const rsRationalFunction<T>& q) const
  {
    rsRationalFunction<T> r;
    ratAdd(num.coeffs, den.coeffs, q.num.coeffs, q.den.coeffs, r.num.coeffs, r.den.coeffs, 
      rsBiggest(tol, q.tol)); 
    //ratAdd(num.coeffs, den.coeffs, q.num.coeffs, q.den.coeffs, r.num.coeffs, r.den.coeffs);
    return r;
  }

  /** Subtracts two rational functions. */
  rsRationalFunction<T> operator-(const rsRationalFunction<T>& q) const
  {
    rsRationalFunction<T> r;
    ratAdd(num.coeffs, den.coeffs, q.num.coeffs, q.den.coeffs, r.num.coeffs, r.den.coeffs, 
      rsBiggest(tol, q.tol), T(1), T(-1));
    return r;
  }

  /** Multiplies two rational functions. */
  rsRationalFunction<T> operator*(const rsRationalFunction<T>& q) const
  {
    rsRationalFunction<T> r;
    ratMul(num.coeffs, den.coeffs, q.num.coeffs, q.den.coeffs, r.num.coeffs, r.den.coeffs, 
      rsBiggest(tol, q.tol));
    return r;
  }

  /** Divides two rational functions. */
  rsRationalFunction<T> operator/(const rsRationalFunction<T>& q) const
  {
    rsRationalFunction<T> r;
    ratDiv(num.coeffs, den.coeffs, q.num.coeffs, q.den.coeffs, r.num.coeffs, r.den.coeffs, 
      rsBiggest(tol, q.tol));
    return r;
    // implement conversion operator in rsPolynomial to get rid of accessing coeffs
  }

  /** Adds another rational function to this one. */
  rsRationalFunction<T>& operator+=(const rsRationalFunction<T>& q)
  { 
    ratAdd(num.coeffs, den.coeffs, q.num.coeffs, q.den.coeffs, num.coeffs, den.coeffs, 
      rsBiggest(tol, q.tol));
    return *this;
  }

  /** Multiplies this rational function by another one. */
  rsRationalFunction<T>& operator*=(const rsRationalFunction<T>& q)
  { 
    ratMul(num.coeffs, den.coeffs, q.num.coeffs, q.den.coeffs, num.coeffs, den.coeffs, 
      rsBiggest(tol, q.tol));
    return *this;
  }


  rsRationalFunction<T> operator*(const T& a) const
  { return rsRationalFunction<T>(num*a, den); } // write a*num




  /** Returns the rational function that results from nesting/composing the given inner rational 
  function with this function as outer function. You can use it like h = f(g) where h,f,g are all 
  rsRationalFunction objects. */
  rsRationalFunction<T> operator()(rsRationalFunction<T> inner) const
  {
    rsRationalFunction<T> r;
    ratNest(inner.num.coeffs, inner.den.coeffs, num.coeffs, den.coeffs, 
      r.num.coeffs, r.den.coeffs, rsBiggest(tol, inner.tol));
    return r;
  }

  /** Evaluates the function at the given input x. */
  T operator()(T x) const
  {
    return num(x) / den(x);
  }

  /** Evaluates the function at the given input z whose type may be different from the 
  coefficient type, for example, for evaluating functions with real coeffs at complex arguments. */
  template<class TArg>
  TArg operator()(TArg z) const 
  { 
    return num(z) / den(z);
  }

  /** Compares this rational function to another one for equality. */
  bool operator==(const rsRationalFunction<T>& q) const 
  { 
    return num == q.num && den == q.den;
  }








  //===============================================================================================
  /** \name Computations on std::vector */

  /** Evaluates the polynomial p a the given x using Horner's algorithm */
  static T polyEval(std::vector<T>& p, T x);

  /** Truncates trailing zeros of the list p */
  static void polyTrunc(std::vector<T>& p, T tol = 0.0);

  /** Makes the polynomial a monic, i.e. divides all coefficients by the  leading coefficient to 
  make the leading coefficient 1. Will result in division by zero error, if p is the zero 
  polynomial. It works in place and will return the leading coefficient (which may or may not be 
  of interest to the caller) */
  static T makeMonic(std::vector<T>& p);

  /** Forms a weighted sum of the two coefficient lists p and q with weights wp and wq 
  respectively. If the resulting list will have trailing zeros, these will be truncated. */
  static std::vector<T> polyAdd(const std::vector<T>& p, const std::vector<T>& q,  
    T tol = 0.0, T wp = 1, T wq = 1);

  /** Subtracts the coefficient list q from the coefficient list p. If the result has trailing 
  zeros, these will be truncated. */
  static std::vector<T> polySub(const std::vector<T>& p, const std::vector<T>& q, T tol = 0.0);

  /** Multiplies two lists of polynomial coefficients by convolution. */
  static std::vector<T> polyMul(const std::vector<T>& p, const std::vector<T>& q, T tol = 0.0);

  /** Divides polynomial p (product) by polynomial d (divisor) and returns the quotient in q and 
  remainder in r */
  static void polyDivMod(std::vector<T> p, std::vector<T> d, 
    std::vector<T>& q, std::vector<T>& r, T tol = 0.0);

  /** Quotient of polynomial division - this corresponds to the integer part of the division of 
  natural numbers. */
  static std::vector<T> polyDiv(std::vector<T> p, std::vector<T> d, T tol);

  /** Remainder of polynomial division */
  static std::vector<T> polyMod(std::vector<T> p, std::vector<T> d, T tol);

  /** Checks, if vector v contains only zeros. */
  static bool isAllZeros(const std::vector<T>& v, T tol); 
  // move to StdContainerTools - maybe it's already there? -> check and delete, if so

  /** Computes the greatest common divisor of polynomials p and q which is defined as the 
  polynomial of highest degree that divides both p and q. Such a polynomial is unique only up to 
  multiplication by a constant, so it is often additionally required to be a monic polynomial to 
  make it unique. This normalization can be controlled by the monic parameter. */
  static std::vector<T> polyGCD(
    const std::vector<T>& p, const std::vector<T>& q, T tol, bool monic = true);

  /** Given the coefficient lists of two polynomials a(x) and b(x), this function computes the 
  coefficient list of the polynomial c(x) that results from nesting a(x) and b(x) where a(x) is the
  inner and b(x) the outer polynomial such that: c(x) = b(a(x)) */
  static std::vector<T> polyNest(const std::vector<T>& a, const std::vector<T>& b);

  /** Reduces rational function p/q to the lowest possible denominator. */
  static void ratReduce(const std::vector<T>& pIn, const std::vector<T>& qIn, 
    std::vector<T>& pOut, std::vector<T>& qOut, T tol);

  /** Multiplies two rational functions represented as lists of coefficients for	numerator and 
  denominator. Computes u/v = (p/q) * (r/s). By default, it will reduce the result to the lowest 
  possible denominator but you can turn that off via the reduced parameter. */
  static void ratMul(const std::vector<T>& p, const std::vector<T>& q, const std::vector<T>& r, 
    const std::vector<T>& s, std::vector<T>& u, std::vector<T>& v, T tol = 0.0, 
    bool reduced = true);

  /** Divides two rational functions */
  static void ratDiv(const std::vector<T>& p, const std::vector<T>& q, const std::vector<T>& r, 
    const std::vector<T>& s, std::vector<T>& u, std::vector<T>& v, T tol = 0.0, 
    bool reduced = true);

  /** Adds two rational functions represented as lists of coefficients for numerator and 
  denominator with optional weighting. It computes numerator and denominator of 
  nr/dr = w1*n1/d1 + w2*n2/d2. */
  static void ratAdd(const std::vector<T>& n1, const std::vector<T>& d1, 
    const std::vector<T>& n2, const std::vector<T>& d2, std::vector<T>& nr, std::vector<T>& dr, 
    T tol = T(0), T w1 = T(1), T w2 = T(1));

  /** Nesting of an inner rational function ni/di with an outer polynomial po. */
  static void ratPolyNest(const std::vector<T>& ni, const std::vector<T>& di, 
    const std::vector<T>& po, std::vector<T>& nr, std::vector<T>& dr, T tol = 0.0);

  /** Nesting of two rational functions. The inner functions numerator and denominator are given by
  nI, dI, likewise for the outer and nO, dO. The result is returned in nR, dR. */
  static void ratNest(const std::vector<T>& nI, const std::vector<T>& dI, const std::vector<T>& nO,
    const std::vector<T>& dO, std::vector<T>& nR, std::vector<T>& dR, T tol = 0.0);


  //===============================================================================================
  /** \name Computations on raw coefficient arrays */

  /** A general routine for performing a partial fraction expansion of a rational function with 
  known poles. You must pass the coefficient arrays of numerator and denominator, the array of 
  poles along with an array of their respective multiplicities (these two arrays are both of length
  "numDistinctPoles". On return, the pfeCoeffs will contain the residues corresponding the the 
  given poles. The array is of length equal to the total number of poles, each counted with its 
  multiplicity, which is the same number as the degree of the denominator. The ordering of the 
  pfeCoeffs corresponds to the ordering of the poles and if a pole p_i has a  multiplicity k, then 
  you will first get the residue for r_i/(x-p_i), then r_i/(x-p_i)^2, etc. If the rational function
  is not strictly proper (i.e. numeratorDegree >= denominatorDegree), the expansion will also 
  feature a polynomial part which is returned in polyCoeffs - which must be of length 
  numeratorDegree+1 which is the *maximum potential* number of coefficients in the polynomial part.
  It's *not* enough if it is just long enough to hold the *actual* number of coeffs. If the actual 
  number of polynomial coeffs is less than the maximum possible number, you'll notice this by 
  getting (numerically close to) zero-valued coeffs for the higher order terms. Depending on 
  whether or not all poles are distinct, it dipsatches between the two functions 
  partialFractionExpansionDistinctPoles and partialFractionExpansionMultiplePoles for the actual
  work. It's the high-level interface function meant to be called from client code.

  Note that function may manipulate the incoming numerator and denominator arrays in place (in 
  order to make the denominator monic and/or divide out the polynomial part from the numerator). 
  So, if you still need them in their original form after the function returns, you should create a
  local copy to pass into the function. */
  template<class R>
  static void partialFractionExpansion(
    std::complex<R>* numerator, int numeratorDegree,
    std::complex<R>* denominator, int denominatorDegree,
    const std::complex<R>* poles, const int* multiplicities, int numDistinctPoles,
    std::complex<R>* pfeCoeffs, std::complex<R>* polyCoeffs = nullptr);
  // -may allocate heap memory (in case of multiple poles)
  // ToDo:
  // -have a higher-level version of the function that doesn't require the poles to be known and
  //  passed in by the caller (the function should find them itself via a root finder)
  // maybe rename to partialFractions
  // maybe needs a separate template parameter R for real numbers


  /** A routine to perform a partial fraction expansion of a strictly proper rational function when
  all poles are distinct. In this common special case, a much more efficient and numerically more 
  precise (supposedly - verify that) algorithm can be used than in the general case where poles may 
  have multiplicities. This function implements the cover-up method, see:
  https://en.wikipedia.org/wiki/Heaviside_cover-up_method  
  called interbally by partialFractionExpansion  */
  template<class R>
  static void partialFractionExpansionDistinctPoles(
    std::complex<R> *num, int numDeg, std::complex<R> *den, int denDeg,
    const std::complex<R> *poles, std::complex<R> *pfeCoeffs);

  /** A routine to perform a partial fraction expansion of a strictly proper rational function when
  some poles may have mutliplicities. The algorithm implemented here solves the linear system
  of equations that results from equating the original rational function to a partial fraction
  expansion with undetermined coefficients, multiplying both sides by the denominator and equating
  coefficients of the polynomials on both sides. called interbally by partialFractionExpansion */
  template<class R>
  static void partialFractionExpansionMultiplePoles(
    const std::complex<R>* num, int numDeg, const std::complex<R>* den, int denDeg,
    const std::complex<R>* poles, const int* multiplicities, int numDistinctPoles,
    std::complex<R>* pfeCoeffs);
  // allocates heap memory


  // convenience functions to work on std::vector (not recommended for production code due to 
  // extra memory allocations - mostly for making tests and experiments more convenient):

  template<class R>
  static std::vector<std::complex<R>> partialFractions(
    const std::vector<std::complex<R>>& numerator,
    const std::vector<std::complex<R>>& denominator,
    const std::vector<std::complex<R>>& poles);

  template<class R>
  static std::vector<std::complex<R>> partialFractions(
    const std::vector<std::complex<R>>& numerator,
    const std::vector<std::complex<R>>& denominator,
    const std::vector<std::complex<R>>& poles,
    const std::vector<int>& multiplicities);
  // maybe have another parameter later to switch between algorithms




protected:

  rsPolynomial<T> num, den;  // numerator and denominator polynomials
  T tol = T(0);              // tolerance for reduction in the arithmetic operators

};

/** Multiplies a number and a rational function. */
template<class T>
inline rsRationalFunction<T> operator*(const T& s, const rsRationalFunction<T>& p)
{
  rsRationalFunction<T> q = p;
  q.scaleNumerator(s);
  return q;
}



//}

/*
todo:
-implement a conversion constructor that can take a polynomial and make a rationla function from it

*/