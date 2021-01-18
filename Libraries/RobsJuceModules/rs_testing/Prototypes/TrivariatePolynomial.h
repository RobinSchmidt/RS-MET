#ifndef RAPT_TRIVARIATEPOLYNOMIAL_H
#define RAPT_TRIVARIATEPOLYNOMIAL_H

template<class T>
class rsTrivariatePolynomial
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime 

  rsTrivariatePolynomial() {}

  rsTrivariatePolynomial(int degreeX, int degreeY, int degreeZ)
  {
    setDegrees(degreeX, degreeY, degreeZ);
  }

  rsTrivariatePolynomial(int degreeX, int degreeY, int degreeZ, std::initializer_list<T> l)
  {
    setDegrees(degreeX, degreeY, degreeZ);
    std::vector<T> vl(l);
    coeffs.setData(vl);  
  }


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setDegrees(int degreeX, int degreeY, int degreeZ)
  {
    std::vector<int> shape({degreeX+1, degreeY+1, degreeZ+1});
    coeffs.setShape(shape);
    coeffs.setToZero();  // todo: take over old data
  }

  void fillRandomly(T min = T(0), T max = T(1), int seed = 0, bool roundToInt = false)
  {
    coeffs.fillRandomly(min, max, seed, roundToInt);
  }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  int getDegreeX() const { return coeffs.getExtent(0)-1; }
  int getDegreeY() const { return coeffs.getExtent(1)-1; }
  int getDegreeZ() const { return coeffs.getExtent(2)-1; }

  bool isZero(T tolerance) const { return coeffs.isAllZeros(tolerance); }

  bool isCloseTo(const rsTrivariatePolynomial<T>& p, T tol) const;

  /** Returns true, iff the vector field represented by the 3 trivariate functions 
  fx(x,y,z), fy(x,y,z), fz(x,y,z) has a vector potential, i.e. a vector field whose curl is given
  by fx, fy, fz. The necessarry and sufficient condition for such a vector potential to exist is 
  that the divergence of fx, fy, fz must be zero. @see vectorPotential. */
  static bool hasVectorPotential(const rsTrivariatePolynomial<T>& fx, 
    const rsTrivariatePolynomial<T>& fy, const rsTrivariatePolynomial<T>& fz, T tolerance = T(0))
  { return divergence(fx, fy, fz).isZero(tolerance); }

  /** Checks whether F,G,H is a vector potential for f,g,h. */
  static bool isVectorPotential(const rsTrivariatePolynomial<T>& F, 
    const rsTrivariatePolynomial<T>& G, const rsTrivariatePolynomial<T>& H, 
    const rsTrivariatePolynomial<T>& f, const rsTrivariatePolynomial<T>& g, 
    const rsTrivariatePolynomial<T>& h, T tolerance = T(0));



  //-----------------------------------------------------------------------------------------------
  // \name Evaluation

  T evaluate(T x, T y, T z) const;

  rsBivariatePolynomial<T> evaluateX(T x) const;

  // todo: partial evaluation for x,y,z (returning bivriate polynomials), xy,xz,yz (returning 
  // univariate polynomials)


  //-----------------------------------------------------------------------------------------------
  // \name Arithmetic

  static void weightedSum(const rsTrivariatePolynomial<T>& p, T wp,
    const rsTrivariatePolynomial<T>& q, T wq, rsTrivariatePolynomial<T>& r)
  { 
    int L = rsMax(p.getDegreeX(), q.getDegreeX());
    int M = rsMax(p.getDegreeY(), q.getDegreeY());
    int N = rsMax(p.getDegreeZ(), q.getDegreeZ());
    r.setDegrees(L, M, N);
    rsMultiArray<T>::weightedSum(p.coeffs, wp, q.coeffs, wq, r.coeffs); 
  }

  rsTrivariatePolynomial<T> operator+(const rsTrivariatePolynomial<T>& p) const
  { rsTrivariatePolynomial<T> r; weightedSum(*this, T(1), p, T(1), r); return r; }

  rsTrivariatePolynomial<T> operator-(const rsTrivariatePolynomial<T>& p) const
  { rsTrivariatePolynomial<T> r; weightedSum(*this, T(1), p, T(-1), r); return r; }

  rsTrivariatePolynomial<T> operator*(const rsTrivariatePolynomial<T>& p) const;

  static rsPolynomial<T> compose(const rsTrivariatePolynomial<T>& p,
    const rsPolynomial<T>& x, const rsPolynomial<T>& y, const rsPolynomial<T>& z);

  static rsBivariatePolynomial<T> compose(const rsTrivariatePolynomial<T>& p,
    const rsBivariatePolynomial<T>& x, const rsBivariatePolynomial<T>& y,
    const rsBivariatePolynomial<T>& z);

  bool operator==(const rsTrivariatePolynomial<T>& p) const
  { return coeffs == p.coeffs; }

  void negate() { coeffs.negate(); }


  //-----------------------------------------------------------------------------------------------
  // \name Calculus

  // todo: use pointers for output parameters

  /** Returns the derivative with respect to x. */
  rsTrivariatePolynomial<T> derivativeX() const;
  static void derivativeX(const rsMultiArray<T>& c, rsMultiArray<T>& d);

  /** Returns the derivative with respect to y. */
  rsTrivariatePolynomial<T> derivativeY() const;
  static void derivativeY(const rsMultiArray<T>& c, rsMultiArray<T>& d);

  /** Returns the derivative with respect to z. */
  rsTrivariatePolynomial<T> derivativeZ() const;
  static void derivativeZ(const rsMultiArray<T>& c, rsMultiArray<T>& d);

  /** Given a scalar field f(x,y,z), this function computes the 3 partial derivatives with respect 
  to x,y,z and assigns them to f_x,f_y,f_z. */
  static void gradient(const rsTrivariatePolynomial<T>& f, rsTrivariatePolynomial<T>& f_x, 
    rsTrivariatePolynomial<T>& f_y, rsTrivariatePolynomial<T>& f_z)
  { f_x = f.derivativeX(); f_y = f.derivativeY(); f_z = f.derivativeZ(); }

  static rsTrivariatePolynomial<T> divergence(const rsTrivariatePolynomial<T>& fx, 
    const rsTrivariatePolynomial<T>& fy, const rsTrivariatePolynomial<T>& fz)
  { return fx.derivativeX() + fy.derivativeY() + fz.derivativeZ(); }

  static void curl(const rsTrivariatePolynomial<T>& fx, const rsTrivariatePolynomial<T>& fy, 
    const rsTrivariatePolynomial<T>& fz, rsTrivariatePolynomial<T>& cx, 
    rsTrivariatePolynomial<T>& cy, rsTrivariatePolynomial<T>& cz);

  rsTrivariatePolynomial<T> laplacian()
  {
    return derivativeX().derivativeX() + derivativeY().derivativeY() + derivativeZ().derivativeZ();
  }
  // todo: implement computing the 2nd partial derivatives in one go -> optimization


  static void integralX(const rsMultiArray<T>& a, rsMultiArray<T>& ai, T c = T(0));
  rsTrivariatePolynomial<T> integralX(T c = T(0)) const;

  static void integralY(const rsMultiArray<T>& a, rsMultiArray<T>& ai, T c = T(0));
  rsTrivariatePolynomial<T> integralY(T c = T(0)) const;

  static void integralZ(const rsMultiArray<T>& a, rsMultiArray<T>& ai, T c = T(0));
  rsTrivariatePolynomial<T> integralZ(T c = T(0)) const;


  template<class Ta, class Tb>
  rsBivariatePolynomial<T> integralX(Ta a, Tb b) const;

  // todo: integralY(Ta a, Tb b), integralZ(Ta a, Tb b) -> copy,paste,edit




  /** Computes the triple integral of the polynomial over the given cuboid. This function 
  performs the integration over x first, then over y, then over z. */
  T tripleIntegralXYZ(T x0, T x1, T y0, T y1, T z0, T z1) const
  { return integralX(x0, x1).doubleIntegralXY(y0, y1, z0, z1); }

  // maybe have a function where only the innermost limits x0,x1 are constant, the middle limits 
  // y0,y1 are univariate polynomials in x and the outermost limits z0, z1 are bivariate polynomials
  // in x,y


  /** Given a vector field represented by the 3 trivariate functions f(x,y,z), g(x,y,z), h(x,y,z),
  this function computes a scalar potential P(x,y,z) for the given vector field. A scalar
  potential (often just called unqualified "potential") is a scalar field whose gradient gives the 
  original vector field. Note that sometimes (especially in physics) the convention is used that 
  the negative gradient instead of the gradient itself is used. This convention is not adopted 
  here - if you want to adopt it, you need to add the minus yourself. The necessarry and sufficient
  condition for a scalar potential to exist is that the given vector field must have zero curl. The
  function assumes this condition to hold for f,g,h. If it doesn't hold, the computed result is 
  meaningless. A scalar potential is unique up to some constant.
  https://en.wikipedia.org/wiki/Scalar_potential  */ 
  static rsTrivariatePolynomial<T> scalarPotential(const rsTrivariatePolynomial<T>& f, 
    const rsTrivariatePolynomial<T>& g, const rsTrivariatePolynomial<T>& h);

  /** Given a vector field represented by the 3 trivariate functions f(x,y,z), g(x,y,z), h(x,y,z),
  this function computes a vector potential F(x,y,z), G(x,y,z), H(x,y,z) for the given vector 
  field. A vector potential is another vector field whose curl gives the original vector field. 
  The necessarry and sufficient condition for a vector potential to exist is that the given vector
  field must have zero divergence. The function assumes this condition to hold for f,g,h. If it 
  doesn't hold, the computed result is meaningless. A vector potential is only unqiue up to some
  vector field whose curl is zero. That means, given a vector potential for f,g,h, you can add any
  curl-free (i.e. conservative) vector field to it and it will still be a vector potential for 
  f,g,h. This leaves a lot of freedom of choice how to define F,G,H. Some of the degrees of freedom
  are used here to set H(x,y,z) to be identically zero. 
  https://en.wikipedia.org/wiki/Vector_potential  */
  static void vectorPotential(const rsTrivariatePolynomial<T>& f, 
    const rsTrivariatePolynomial<T>& g, const rsTrivariatePolynomial<T>& h, 
    rsTrivariatePolynomial<T>& F, rsTrivariatePolynomial<T>& G, rsTrivariatePolynomial<T>& H);

  /** Computes the path integral of a vector field fx(x,y,z), fy(x,y,z), fz(x,y,z) over a path that
  is given parametrically by 3 univariate polynomials x(t), y(t), z(t) where the parameter t runs 
  from a to b. */
  static T pathIntegral(const rsTrivariatePolynomial<T>& fx, const rsTrivariatePolynomial<T>& fy, 
    const rsTrivariatePolynomial<T>& fz, const rsPolynomial<T>& x, const rsPolynomial<T>& y, 
    const rsPolynomial<T>& z, T a, T b);

  /** Computes the path integral of a vector field fx(x,y,z), fy(x,y,z), fz(x,y,z) over a path that
  is given by linear segments connecting the vertices in the given array of 3D points. */
  static T pathIntegral(const rsTrivariatePolynomial<T>& fx, const rsTrivariatePolynomial<T>& fy, 
    const rsTrivariatePolynomial<T>& fz, const std::vector<rsVector3D<T>>& path);

  /** Computes the flux of a vector field given by 3 functions fx(x,y,z), fy(x,y,z), fz(x,y,z)
  through a parametric surface patch given by x(u,v), y(u,v), z(u,v) where u and v run from u0 to
  u1 and v0 to v1 respectively. If the vector field describes a fluid velocity, the flux integral 
  measures, how much of the fluid flows through the given surface patch per unit time. */
  static T fluxIntegral(const rsTrivariatePolynomial<T>& fx, const rsTrivariatePolynomial<T>& fy, 
    const rsTrivariatePolynomial<T>& fz, const rsBivariatePolynomial<T>& x, 
    const rsBivariatePolynomial<T>& y, const rsBivariatePolynomial<T>& z, T u0, T u1, T v0, T v1);

  /** Flux through a triangular patch between vertices P0, P1, P2. */
  static T fluxIntegral(const rsTrivariatePolynomial<T>& fx, const rsTrivariatePolynomial<T>& fy, 
    const rsTrivariatePolynomial<T>& fz, const rsVector3D<T>& P0, const rsVector3D<T>& P1, 
    const rsVector3D<T>& P2);

  /** Computes the flux of a vector field coming out of a cuboid bounded by the given coordinates. 
  By Gauss theorem, this should be equal to the triple integral of the divergence and i think, 
  computing it that way is more efficient and accurate. This function is mostly for the sake of 
  completeness and proof of concept. */
  static T outfluxIntegral(const rsTrivariatePolynomial<T>& fx, 
    const rsTrivariatePolynomial<T>& fy, const rsTrivariatePolynomial<T>& fz,
    T x0, T x1, T y0, T y1, T z0, T z1);
  // maybe rename to outfluxIntegralDirect, implement an outfluxIntegralViaDivergence and then
  // outfluxIntegral becomes just an alias to the better algorithm ...maybe it could dispatch, if
  // the choice of best algorithm is not always the same, so it becomes a meta-algorithm but then
  // we would have to figure out under which conditions which algo is better...but probably, the 
  // divergence-based algo is always better...i think...maybe



  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Read and write access to the (i,j)th coefficient. */
  T& coeff(int i, int j, int k) { return coeffs(i, j, k); }

  /** Read access to the (i,j)th coefficient. */
  //const T& coeff(int i, int j) const { return coeffs(i, j); }

  T getCoeffPadded(int i, int j, int k, T padding = T(0)) const 
  { return coeffs.getElementPadded3D(i, j, k, padding); }


protected:


  rsMultiArray<T> coeffs;

};


#endif