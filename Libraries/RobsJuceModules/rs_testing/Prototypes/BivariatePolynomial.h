#ifndef RAPT_BIVARIATEPOLYNOMIAL_H
#define RAPT_BIVARIATEPOLYNOMIAL_H

/** Class for representing bivariate polynomials, i.e. polynomials in two variables x and y. They 
are represented by an MxN matrix A of coefficients which is supposed to be sandwiched between two 
vectors of powers of x and y as in X^T * A * Y, where X,Y denote the vectors constructed from 
powers of x and y and ^T denotes transposition, making X^T a row vector. For example, a polynomial
that has degree 2 in x and degree 3 in y looks like:

  p(x,y) = |x^0 x^1 x^2| * |a00 a01 a02 a03| * |y^0|
                           |a10 a11 a12 a13|   |y^1|
                           |a20 a21 a22 a23|   |y^2|
                                               |y^3|

         =   a00*x^0*y^0 + a01*x^0*y^1 + a02*x^0*y^2 + a03*x^0*y^3
           + a10*x^1*y^0 + a11*x^1*y^1 + a12*x^1*y^2 + a13*x^1*y^3
           + a20*x^2*y^0 + a21*x^2*y^1 + a22*x^2*y^2 + a23*x^2*y^3

so the shape of the matrix is 3x4 with M=3, N=4. In general, it's (degX+1)x(degY+1) where degX, 
degY are the degrees with respect to x and y. */

template<class T>
class rsBivariatePolynomial
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime 

  rsBivariatePolynomial() {}

  rsBivariatePolynomial(int degreeX, int degreeY) : coeffs(degreeX+1, degreeY+1) {}

  rsBivariatePolynomial(int degreeX, int degreeY, std::initializer_list<T> l) 
    : coeffs(degreeX+1, degreeY+1, l) {}


  //-----------------------------------------------------------------------------------------------
  // \name Evaluation (High Level)

  /** Evaluates the polynomial at the given (x,y). */
  T evaluate(T x, T y) const;

  /** Evaluates the polynomial for a given x. The result is a univariate polynomial in y. */
  rsPolynomial<T> evaluateX(T x) const;

  /** Evaluates the polynomial for a given y. The result is a univariate polynomial in x. */
  rsPolynomial<T> evaluateY(T y) const;

  /** Takes a univariate polynomial x = x(y) as input and replaces each occurrence of x in this
  bivariate polynomial p(x,y) by the polynomial expression x(y). This results in a univariate
  polynomial in y. */
  rsPolynomial<T> evaluateX(const rsPolynomial<T>& x) const;

  /** Takes a univariate polynomial y = y(x) as input and replaces each occurrence of y in this
  bivariate polynomial p(x,y) by the polynomial expression y(x). This results in a univariate
  polynomial in x. */
  rsPolynomial<T> evaluateY(const rsPolynomial<T>& y) const;

  //-----------------------------------------------------------------------------------------------
  // \name Evaluation (Low Level)

  /** Evaluates the polynomial for a given x. The result is a univariate polynomial in y whose 
  coefficients are stored in py. */
  void evaluateX(T x, T* py) const;
  // make function static, taking an rsMatrixView parameter, just like derivativeX

  /** Evaluates the polynomial for a given y. The result is a univariate polynomial in x whose 
  coefficients are stored in px. */
  void evaluateY(T y, T* px) const;

  /** Used internally by evaluateX(const rsPolynomial<T>& x). Let M,N be the degrees in x and y of
  this bivariate polynomial an K be the degree of the passed univariate polynomial. Then, the 
  degree of the result py will be given by K*M + N. The workspace must have a size of K*M+1.  */
  void evaluateX(const T* px, int xDeg, T* py, T* workspace) const;

  /** Used internally by evaluateY(const rsPolynomial<T>& y). Let M,N be the degrees in x and y of
  this bivariate polynomial an K be the degree of the passed univariate polynomial. Then, the 
  degree of the result px will be given by K*N + M. The workspace must have a size of K*N+1.  */
  void evaluateY(const T* py, int yDeg, T* px, T* workspace) const;


  //-----------------------------------------------------------------------------------------------
  // \name Arithmetic

  /** Computes the coefficients of a bivariate polynomial r that is given as the product of two
  univariate polynomial p and q that are functions of x and y alone, respectively, such that 
  r(x,y) = p(x) * q(y). */
  static rsBivariatePolynomial<T> multiply(const rsPolynomial<T>& p, const rsPolynomial<T>& q);

  static void multiply(const T* p, int pDeg, const T* q, int qDeg, rsMatrixView<T>& r);

  static void weightedSum(const rsMatrixView<T>& p, T wp, const rsMatrixView<T>& q, T wq, 
    rsMatrixView<T>& r);


  /** Computes the coefficients of a bivariate polynomial r that is given as the composition of 
  a linear combination of x and y and a univariate polynomial p: r(x,y) = p(a*x + b*y). */
  static rsBivariatePolynomial<T> composeWithLinear(const rsPolynomial<T>& p, T a, T b);

  static rsBivariatePolynomial<T> composeWithLinearOld(const rsPolynomial<T>& p, T a, T b);

  /** Given a bivariate polynomial p(x,y) and two univariate polynomials x(t), y(t), this function 
  computes p(x(t),y(t)) which is a univariate polynomial in t.  */
  static rsPolynomial<T> compose(const rsBivariatePolynomial<T>& p,
    const rsPolynomial<T>& x, const rsPolynomial<T>& y);



  /** Multiplies this bivariate polynomial with a univariate polynomial in y only and returns the 
  result which is again a bivariate polynomial. This amounts to convolving each row of our 
  coefficient matrix with the coefficient array of p. */
  rsBivariatePolynomial<T> multiplyY(const rsPolynomial<T>& p) const;


  void negate() { coeffs.negate(); }

  void scale(T s) { coeffs.scale(s); }

  // todo: make a similar multiplyX method - this needs to convolve the columns with p, so we will
  // need a convolution routine with strides


  /** Given a complex valued bivariate polynomial, this function splits it into two real-valued 
  bivariate polynomials...tbc... */
  static void splitRealImag(const rsBivariatePolynomial<std::complex<T>>& p,
    rsBivariatePolynomial<T>& pRe, rsBivariatePolynomial<T>& pIm);


  template<class T2>
  rsBivariatePolynomial<T2> convert(T2 dummy) const
  {
    rsBivariatePolynomial<T2> p(getDegreeX(), getDegreeY());
    for(int m = 0; m <= getDegreeX(); m++)
      for(int n = 0; n <= getDegreeY(); n++)
        p.coeff(m, n) = T2(coeffs(m, n));
    return p;
  }


  /** Given a complex polynomial (or more generally, a complex function), the associated Polya 
  vector field is the complex conjugate of the vector field that would result from just intepreting
  real and imaginary parts of the function's output as x- and y-coordinates of a 2D vector field. 
  That just means, the y-part is negated, i.e. fx(x,y) = Re(p(z)), fy(x,y) = -Im(p(z)) where 
  z = x + i*y. The reason for this seemingly unnatural negation is that the resulting vector field 
  will be conservative when the original function is analytic. This is a desirable property for 
  vector fields and it does not hold true without the negation (-> verify that).  */
  static void polyaVectorField(const rsPolynomial<std::complex<T>>& p,
    rsBivariatePolynomial<T>& px, rsBivariatePolynomial<T>& py);


  //-----------------------------------------------------------------------------------------------
  // \name Factory

  /** Creates the zero polynomial. */
  static rsBivariatePolynomial<T> zero() { return rsBivariatePolynomial<T>(0, 0, { 0 }); };


  //-----------------------------------------------------------------------------------------------
  // \name Calculus

  static void derivativeX(const rsMatrixView<T>& c, rsMatrixView<T>& d);
  rsBivariatePolynomial<T> derivativeX() const;

  static void derivativeY(const rsMatrixView<T>& c, rsMatrixView<T>& d);
  rsBivariatePolynomial<T> derivativeY() const;


  static void integralX(const rsMatrixView<T>& p, rsMatrixView<T>& pi, T c = T(0));
  rsBivariatePolynomial<T> integralX(T c = T(0)) const;
  // should the integration constant be a univariate polynomial in y instead of a fixed value?
  // ...yes, i think so ...or maybe get rid of the paramneter alltogether and just leave the matrix
  // entries corresponding to the terms not involving x as is

  static void integralY(const rsMatrixView<T>& p, rsMatrixView<T>& pi, T c = T(0));
  rsBivariatePolynomial<T> integralY(T c = T(0)) const;


  /** Computes the definite integral of the polynomial with respect to x with the integration 
  limits a,b. The result is a univariate polynomial in y whose coefficients are stored in py. */
  //static void integralX(T a, T b, T* py);

  /** Computes the definite integral of the polynomial with respect to y with the integration 
  limits a,b. The result is a univariate polynomial in x whose coefficients are stored in px. */
  //static void integralY(T a, T b, T* px);

  /** Computes the definite integral of the polynomial with respect to x with the integration 
  limits a,b. The result is a univariate polynomial in y. The types Ta, Tb can both be 
  independently either T (for a constant integration limit) or rsPolynomial<T> (for an integration
  limit that is a univariate polynomial in the variable that is not integrated over, here y). */
  template<class Ta, class Tb>
  rsPolynomial<T> integralX(Ta a, Tb b) const;
  // todo: use const references (important when Ta, Tb are polynomials)
  // needs tests for when a and/or b are polynomials
  // if it's called with integer parameters a, b the compiler actually generates a version with
  // Ta,Tb = int even if T = double...hmm...that may be a bit bloatsome...maybe client code should
  // make sure, it does actually call it with double arguments in this case

  /** Computes the definite integral of the polynomial with respect to y with the integration 
  limits a,b. The result is a univariate polynomial in x. */
  template<class Ta, class Tb>
  rsPolynomial<T> integralY(Ta a, Tb b) const;
  // needs tests for when a and/or b are polynomials

  /** Given two bivariate polynomials px(x,y), py(x,y) that together constitute a 2D vector field 
  and assuming that this vector field is conservative (implying that a potential exists), this 
  function computes the potential. The potential of a 2D vector field given by fx(x,y), fy(x,y) is
  a single bivariate function P(x,y) whose partial derivatives with respect to x and y give the 
  original two functions fx, fy. Of course, one function gives in general less information than 
  two, so this works only for special kinds of vector fields, namely conservative vector fields. 
  The caller must ensure that px, py actually satisfy this condition - if they don't, the returned
  function is meaningless. 
  See: https://mathinsight.org/conservative_vector_field_find_potential  */
  static rsBivariatePolynomial<T> getPotential(
    const rsBivariatePolynomial<T>& px, const rsBivariatePolynomial<T>& py);
  // be consistent with regard to using get - either we should rename integralX to getIntegralX 
  // etc. or rename getPotenteial to potential - choose the variant that is consistent with 
  // rsPolynomial and rsNumericDifferentiator
  // -the convention is that the potential's *NEGATIVE* gradient should give the original 
  //  functions back - here we take just the gradient - change that...we may need to ripple the 
  //  negation through to the getPolyaPotential
  // -document how it can be easisly checked, if a potential exists - i think, the condition is 
  //  px_x == py_y or similar

  /** The Polya vector field of an analytic complex function (such as a polynomial) is 
  conservative, so a potential exists for such a Polya vector field. This function computes that 
  potential for a given (complex) Polynomial. The result is a bivariate real polynomial whose 
  partial derivatives with respect to x and y give the real part and the negative imaginary part 
  of the original complex polynomial. */
  static rsBivariatePolynomial<T> getPolyaPotential(const rsPolynomial<std::complex<T>>& p);


  rsBivariatePolynomial<T> getHarmonicConjugate() const;


  rsBivariatePolynomial<T> getLaplacian() const
  { return derivativeX().derivativeX() + derivativeY().derivativeY(); }
  // rename to laplacian

  // implement vector Laplacian:
  // https://en.wikipedia.org/wiki/Laplace_operator#Two_dimensions


  static rsBivariatePolynomial<T> divergence2D(
    const rsBivariatePolynomial<T>& fx, const rsBivariatePolynomial<T>& fy)
  { return fx.derivativeX() + fy.derivativeY(); }
  // needs test, get rid of the 2D qualifier - it's clear that we need 2D divergence when we deal 
  // with bivariate polynomials

  static rsBivariatePolynomial<T> curl2D(
    const rsBivariatePolynomial<T>& fx, const rsBivariatePolynomial<T>& fy)
  { return fy.derivativeX() - fx.derivativeY(); }
  // needs test

  // https://en.wikipedia.org/wiki/Vector_calculus_identities#Divergence






  /** Computes the double integral of the polynomial over the given rectangle. This function 
  performs the integration over x first and then the integration over y. */
  T doubleIntegralXY(T x0, T x1, T y0, T y1) const;

  /** Computes the double integral of the polynomial over the given rectangle. This function 
  performs the integration over y first and then the integration over x. */
  T doubleIntegralYX(T x0, T x1, T y0, T y1) const;

  /** Given a scalar field p(x,y) and a parametric curve x(t), y(t) and start- and end-values a,b 
  for the parameter t, this function computes the path integral (a.k.a. line integral, curve 
  integral or contour integral) of the (polynomial) scalar field along the given (polynomial) path.
  Currently, this is implemented only for linear paths, i.e. the degrees of x(t),y(t) must be at 
  most 1. Path integrals over scalar fields can be thought of as computing the area of a curtain. 
  Imagine the path on the xy-plane being projected up to the surface landscape of p(x,y) and a 
  curtain hanging down from this projected line to the xy-plane. Of course, the "landscape" may
  also lie below the xy-plane, in which case the result would become negative. So it's more like
  a signed area difference. see: https://en.wikipedia.org/wiki/Line_integral */
  static T pathIntegral(const rsBivariatePolynomial<T>& p, 
    const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b);

  /** Given a vector field p(x,y), q(x,y) and a parametric curve x(t), y(t) and start- and 
  end-values a,b for the parameter t, this function computes the path integral (a.k.a. line 
  integral, curve integral or contour integral) of the (polynomial) vector field along the given
  (polynomial) path. Path integrals over vector fields can be thought of as computing the work that
  a force field does on a particle that moves through the field along the given path. It integrates
  over the scalar product of the vector field's direction vectors with the curve's "velocity" 
  vectors, i.e. over the component of the vector field that points along the curve's direction.
  see: https://en.wikipedia.org/wiki/Line_integral */
  static T pathIntegral(const rsBivariatePolynomial<T>& p, const rsBivariatePolynomial<T>& q,
    const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b);
  // aka path integral of 2nd kind or vector path integral?
  // -maybe rename to lineIntegral or contourIntegral ...or maybe flowIntegral
  // http://ndp.jct.ac.il/tutorials/infitut2/node61.html

  /** ....
  It integrates over the 2D cross product (which is a scalar) of the vector field's direction 
  vectors with the curve's "velocity" vectors, i.e. over the component of the vector field that 
  points perpendicular to the curve's direction.
  */
  static T fluxIntegral(const rsBivariatePolynomial<T>& p, const rsBivariatePolynomial<T>& q,
    const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b);

  /** Path integral over a vector field around the closed rectangular loop going along the 4 line 
  segments: (x0,y0) -> (x1,y0) -> (x1,y1) -> (x0,y1) -> (x0,y0). If x1 > x0 and y1 > y0, then we go
  right, up, left, down. By Green's theorem, this should be equal to the double integral over the 
  enclosed rectangle of the curl of the vector field. */
  static T loopIntegral(const rsBivariatePolynomial<T>& p, const rsBivariatePolynomial<T>& q,
    T x0, T x1, T y0, T y1);
  // aka circulation? maybe rename to circulationIntegral or loopFlow, flowAroundLoop

  /** Like loopIntegral but using the flux instead of the flow. The flux is the component of the 
  vector field is perpendicular to the curve. The integral measures, how much of a fluid flows
  out of the given rectangle. By Gauss' theorem, this should be equal to the double integral over 
  the enclosed rectangle of the divergence of the vector field. */
  static T outfluxIntegral(const rsBivariatePolynomial<T>& p, const rsBivariatePolynomial<T>& q,
    T x0, T x1, T y0, T y1);
  // i made up this name - figure out, if there is a more proper name for that, maybe loopFlux,
  // fluxThroughLoop

  // todo: implement flux integral around a rectangular loop

  // todo:
  // -implement scalar path integrals - they are possible only numerically, due to the square-root 
  //  in the integrand (for the speed) ...except when x(t), y(t) are both linear - in this case, the 
  //  speed is a constant number and we can just scale everything by it
  // -implement flux and circulation integrals, if possible - compare results to double integral
  //  (Gauss' and Green's theorem - the vector field must satisfy some conditions for them to hold)



  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void initialize(int degX, int degY)
  {
    coeffs.setShape(degX+1, degY+1);
    coeffs.setToZero();
  }
  // todo: write a setDegrees function that takes over the content of the old matrix and fills up
  // with zeros



  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  int getDegreeX() const { return coeffs.getNumRows()-1; }

  int getDegreeY() const { return coeffs.getNumColumns()-1; }

  /** Implements a sort of weak equality comparison of this polynomial with a second polynomial q. 
  It has a tolerance for the comparisons of individual coefficients and in the second polynomial q
  may have a differently shaped coefficient matrix, if the coeffs that have no corresponding partner 
  in the respective other matrix are (close to) zero. So, the "overlapping" coeffs must be close to 
  their partner and the non-overlapping coeffs must be close to zero. */
  bool isCloseTo(const rsBivariatePolynomial<T>& q, T tol = T(0)) const
  {

    int M = rsMax(coeffs.getNumRows(),    q.coeffs.getNumRows());
    int N = rsMax(coeffs.getNumColumns(), q.coeffs.getNumColumns());
    for(int m = 0; m < M; m++) {
      for(int n = 0; n < N; n++) {
        T d = coeffs.getElementPadded(m, n) - q.coeffs.getElementPadded(m, n);
        //if(rsAbs(d) > tol)
        if(rsGreaterAbs(d, tol))
          return false;     }}
    return true;
  }
  // move out of class

  bool isHarmonic(T tol = T(0)) const { return getLaplacian().isCloseTo(zero(), tol); }

  static bool areHarmonicConjugates(
    const rsBivariatePolynomial<T>& u, const rsBivariatePolynomial<T>& v, T tol = T(0));


  //-----------------------------------------------------------------------------------------------
  // \name Operators

  bool operator==(const rsBivariatePolynomial<T>& rhs) const { return coeffs == rhs.coeffs; }


  rsBivariatePolynomial<T> operator+(const rsBivariatePolynomial<T>& q) const 
  { 
    int M = rsMax(coeffs.getNumRows(),    q.coeffs.getNumRows());
    int N = rsMax(coeffs.getNumColumns(), q.coeffs.getNumColumns());
    rsBivariatePolynomial<T> r(M-1, N-1);
    weightedSum(coeffs, T(1), q.coeffs, T(1), r.coeffs);
    return r;
  }

  rsBivariatePolynomial<T> operator-(const rsBivariatePolynomial<T>& q) const 
  { 
    int M = rsMax(coeffs.getNumRows(),    q.coeffs.getNumRows());
    int N = rsMax(coeffs.getNumColumns(), q.coeffs.getNumColumns());
    rsBivariatePolynomial<T> r(M-1, N-1);
    weightedSum(coeffs, T(1), q.coeffs, T(-1), r.coeffs);
    return r;
  }

  rsBivariatePolynomial<T> operator*(const rsBivariatePolynomial<T>& q) const
  {
    rsBivariatePolynomial<T> r(getDegreeX()+q.getDegreeX(), getDegreeY()+q.getDegreeY());
    rsMatrixView<T>::convolve(coeffs, q.coeffs, &r.coeffs);
    return r;
  }

  //rsBivariatePolynomial<T> operator*(T a, const rsBivariatePolynomial<T>& q) const;

  /** Unary minus. */
  rsBivariatePolynomial<T> operator-() const
  {
    rsBivariatePolynomial<T> r = *this;
    r.negate();
    return r;
  }

  /** Evaluation. */
  T operator()(T x, T y) const { return evaluate(x, y); }


  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Read and write access to the (i,j)th coefficient. */
  T& coeff(int i, int j) { return coeffs(i, j); }

  /** Read access to the (i,j)th coefficient. */
  const T& coeff(int i, int j) const { return coeffs(i, j); }

  T getCoeffPadded(int i, int j, T padding = T(0)) const 
  { return coeffs.getElementPadded(i, j, padding); }


protected:

  rsMatrix<T> coeffs;

};
// maybe implement a function getHarmonicConjugate which returns 2*x*y when p = x^2 - y^2, etc.
// ...are they obtained by a rotation? does every polynomial have such a conjugate?

template<class T>
rsBivariatePolynomial<T> operator*(const T& a, const rsBivariatePolynomial<T>& q)
{
  rsBivariatePolynomial<T> r = q; r.scale(a); return r;
}


#endif