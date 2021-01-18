
template<class T>
T rsBivariatePolynomial<T>::evaluate(T x, T y) const
{
  T xm(1), yn(1), r(0);  // x^m, y^n, result
  for(int m = 0; m < coeffs.getNumRows(); m++) {
    yn = T(1);
    for(int n = 0; n < coeffs.getNumColumns(); n++) {
      r += coeffs(m, n) * xm * yn;
      yn *= y;  }
    xm *= x; }
  return r;
}
// todo: try find an algo that works like Horner's rule 

template<class T>
void rsBivariatePolynomial<T>::evaluateX(T x, T* py) const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  rsArrayTools::fillWithZeros(py, N);
  T xm(1);
  for(int m = 0; m < M; m++) {
    for(int n = 0; n < N; n++)
      py[n] += coeffs(m, n) * xm;
    xm *= x; }
}

template<class T>
rsPolynomial<T> rsBivariatePolynomial<T>::evaluateX(T x) const
{
  rsPolynomial<T> py(getDegreeY());
  evaluateX(x, py.getCoeffPointer());
  return py;
}

template<class T>
void rsBivariatePolynomial<T>::evaluateY(T y, T* px) const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  rsArrayTools::fillWithZeros(px, M);
  T yn(1);
  for(int n = 0; n < N; n++) {
    for(int m = 0; m < M; m++)
      px[m] += coeffs(m, n) * yn;
    yn *= y; }
}

template<class T>
rsPolynomial<T> rsBivariatePolynomial<T>::evaluateY(T y) const
{
  rsPolynomial<T> px(getDegreeX());
  evaluateY(y, px.getCoeffPointer());
  return px;
}

template<class T>
void rsBivariatePolynomial<T>::evaluateX(const T* x, int xDeg, T* r, T* xm) const
{
  using AT = rsArrayTools;
  int M = getDegreeX();
  int N = getDegreeY();
  int K = xDeg;
  AT::fillWithZeros(xm, K*M+1); 
  xm[0] = T(1);
  int Km = 1; 
  for(int n = 0; n <= N; n++)
    r[n] = coeffs(0, n);
  for(int m = 1; m <= M; m++) {
    AT::convolveInPlace(xm, Km, x, K+1);
    Km += K;
    for(int i = 0; i < Km; i++)
      for(int n = 0; n <= N; n++)
        r[i+n] += coeffs(m, n) * xm[i]; }
}

template<class T>
rsPolynomial<T> rsBivariatePolynomial<T>::evaluateX(const rsPolynomial<T>& px) const
{
  int M = getDegreeX();
  int N = getDegreeY();
  int K = px.getDegree();
  int L = K*M + N;
  rsPolynomial<T> r(L);
  std::vector<T> xm(K*M+1);
  evaluateX(px.getCoeffPointerConst(), K, r.getCoeffPointer(), &xm[0]);
  return r;
}

template<class T>
void rsBivariatePolynomial<T>::evaluateY(const T* y, int yDeg, T* r, T* yn) const
{
  using AT = rsArrayTools;
  int M = getDegreeX();
  int N = getDegreeY();
  int K = yDeg;
  AT::fillWithZeros(yn, K*N+1); 
  yn[0] = T(1);                           // initially, y^n the constant 1, i.e. y^0
  int Kn = 1;                             // current effective length of yn
  for(int m = 0; m <= M; m++)
    r[m] = coeffs(m, 0);                  // copy coeffs from the left column
  for(int n = 1; n <= N; n++) {
    AT::convolveInPlace(yn, Kn, y, K+1);  // multiply y^n by y again to get the next power
    Kn += K;                              // length of yn has increased by K
    for(int i = 0; i < Kn; i++)
      for(int m = 0; m <= M; m++)
        r[i+m] += coeffs(m, n) * yn[i]; } // accumulate
}

template<class T>
rsPolynomial<T> rsBivariatePolynomial<T>::evaluateY(const rsPolynomial<T>& py) const
{
  int M = getDegreeX();
  int N = getDegreeY();
  int K = py.getDegree();
  int L = K*N + M;            // degree of result
  rsPolynomial<T> r(L);       // result
  std::vector<T> yn(K*N+1);   // workspace, holds coeff-array of powers of y(x)
  evaluateY(py.getCoeffPointerConst(), K, r.getCoeffPointer(), &yn[0]);
  return r;
}

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::multiply(
  const rsPolynomial<T>& p, const rsPolynomial<T>& q)
{
  rsBivariatePolynomial<T> r(p.getDegree(), q.getDegree());
  multiply(p.getCoeffPointerConst(), p.getDegree(), q.getCoeffPointerConst(), q.getDegree(),
    r.coeffs);
  return r;
}

template<class T>
void rsBivariatePolynomial<T>::multiply(const T* p, int pDeg, const T* q, int qDeg, 
  rsMatrixView<T>& r)
{
  rsAssert(r.hasShape(pDeg+1, qDeg+1));
  for(int m = 0; m <= pDeg; m++)
    for(int n = 0; n <= qDeg; n++)
      r(m, n) = p[m] * q[n];
}

template<class T>
void rsBivariatePolynomial<T>::weightedSum(
  const rsMatrixView<T>& p, T wp, const rsMatrixView<T>& q, T wq, rsMatrixView<T>& r)
{
  rsMatrixView<T>::weightedSum(p, wp, q, wq, r);
  /*
  int M = rsMax(p.getNumRows(),    q.getNumRows());
  int N = rsMax(p.getNumColumns(), q.getNumColumns());
  rsAssert(r.hasShape(M, N));
  for(int m = 0; m < M; m++)
    for(int n = 0; n < N; n++)
      r(m, n) = wp * p.getElementPadded(m, n) + wq * q.getElementPadded(m, n);
      */
}
// todo: move to rsMatrixView

template<class T>
void rsBivariatePolynomial<T>::derivativeX(const rsMatrixView<T>& c, rsMatrixView<T>& d)
{
  int M = c.getNumRows();
  int N = c.getNumColumns();
  rsAssert(d.hasShape(M-1, N));
  for(int m = 1; m < M; m++) {
    T s(m);
    for(int n = 0; n < N; n++)
      d(m-1, n) = s * c(m, n); }
}

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::derivativeX() const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  if(M == 1)
    return rsBivariatePolynomial<T>::zero();
  rsBivariatePolynomial<T> q;
  q.coeffs.setShape(M-1, N);
  derivativeX(coeffs, q.coeffs);
  return q;
}

template<class T>
void rsBivariatePolynomial<T>::derivativeY(const rsMatrixView<T>& c, rsMatrixView<T>& d)
{
  int M = c.getNumRows();
  int N = c.getNumColumns();
  rsAssert(d.hasShape(M, N-1));
  for(int n = 1; n < N; n++) {
    T s(n);
    for(int m = 0; m < M; m++)
      d(m, n-1) = s * c(m, n); }
}

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::derivativeY() const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  if(N == 1)
    return rsBivariatePolynomial<T>::zero();
  rsBivariatePolynomial<T> q;
  q.coeffs.setShape(M, N-1);
  derivativeY(coeffs, q.coeffs);
  return q;
}

template<class T>
void rsBivariatePolynomial<T>::integralX(const rsMatrixView<T>& p, rsMatrixView<T>& pi, T c)
{
  int M = p.getNumRows();
  int N = p.getNumColumns();
  rsAssert(pi.hasShape(M+1, N));
  for(int n = 0; n < N; n++)   // i think, this is wrong: just set pi(0,0) to c and the rest to 0
    pi(0, n) = c; 
  for(int m = 1; m <= M; m++) {
    T s = T(1) / T(m);
    for(int n = 0; n < N; n++)
      pi(m, n) = s * p(m-1, n); }
}

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::integralX(T c) const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  rsBivariatePolynomial<T> q;
  q.coeffs.setShape(M+1, N);
  integralX(coeffs, q.coeffs, c);
  return q;
}

template<class T>
void rsBivariatePolynomial<T>::integralY(const rsMatrixView<T>& p, rsMatrixView<T>& pi, T c)
{
  int M = p.getNumRows();
  int N = p.getNumColumns();
  rsAssert(pi.hasShape(M, N+1));
  for(int m = 0; m < M; m++)   // i think, this is wrong: just set pi(0,0) to c and the rest to 0
    pi(m, 0) = c; 
  for(int n = 1; n <= N; n++) {
    T s = T(1) / T(n);
    for(int m = 0; m < M; m++)
      pi(m, n) = s * p(m, n-1); }
}

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::integralY(T c) const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  rsBivariatePolynomial<T> q;
  q.coeffs.setShape(M, N+1);
  integralY(coeffs, q.coeffs, c);
  return q;
}

template<class T>
template<class Ta, class Tb>
rsPolynomial<T> rsBivariatePolynomial<T>::integralX(Ta a, Tb b) const
{
  rsBivariatePolynomial<T> P = integralX();
  rsPolynomial<T> Pb = P.evaluateX(b);
  rsPolynomial<T> Pa = P.evaluateX(a);
  return Pb - Pa;
}
// todo: make workspace-based version(s)

template<class T>
template<class Ta, class Tb>
rsPolynomial<T> rsBivariatePolynomial<T>::integralY(Ta a, Tb b) const
{
  rsBivariatePolynomial<T> P = integralY();
  rsPolynomial<T> Pb = P.evaluateY(b);
  rsPolynomial<T> Pa = P.evaluateY(a);
  return Pb - Pa;
}
// todo: make workspace-based version(s)

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::getPotential(
  const rsBivariatePolynomial<T>& px, const rsBivariatePolynomial<T>& py)
{
  //rsAssert(px.derivativeY() == py.derivativeX(), "px, py are not a potential field");
  // we need a weaker notion of equality here: allow different formal shapes of the coeff matrices, 
  // allow tolerance for equality of coefficients, maybe have a function 
  // isPotentialField(px, py, tol) that calls px.isCloseTo(py, tol)

  rsBivariatePolynomial<T> Px, Px_y, gyp, gy;
  Px   = px.integralX();    // integrate px with respect to x
  Px_y = Px.derivativeY();  // differentiate the result with respect to y
  gyp  = py - Px_y;         // g'(y): derivative of integration "constant" g(y)..
  gy   = gyp.integralY();   // ..which is still a function of y
  return Px + gy;           // Px + gy is the desired potential function P(x,y)
}
// -maybe optimize: gyp has nonzero coeffs only for terms that are free of any x
// -maybe implement different algorithms, integrating py with respect to y first, etc. - they 
//  should all give the same result up to roundoff error


template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::getPolyaPotential(
  const rsPolynomial<std::complex<T>>& p)
{
  rsBivariatePolynomial<T> px, py;
  polyaVectorField(p, px, py);
  return getPotential(px, py);       // (px, py) is a potential field -> compute its potential
}

template<class T> 
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::getHarmonicConjugate() const
{
  rsAssert(isHarmonic());  // needs tolerance
  using BiPoly = rsBivariatePolynomial<T>;
  BiPoly u_dx_iy = derivativeX().integralY();
  BiPoly u_dy = derivativeY(); 
  for(int m = 0; m <= u_dy.getDegreeX(); m++)
    for(int n = 1; n <= u_dy.getDegreeY(); n++)
      u_dy.coeff(m, n) = T(0);
  return u_dx_iy - u_dy.integralX();
}
// -needs tests 

template<class T> 
T rsBivariatePolynomial<T>::doubleIntegralXY(T x0, T x1, T y0, T y1) const
{
  rsPolynomial<T>& ix = integralX(x0, x1);  // still a function of y
  return ix.definiteIntegral(y0, y1);       // just a number
}

template<class T> 
T rsBivariatePolynomial<T>::doubleIntegralYX(T x0, T x1, T y0, T y1) const
{
  rsPolynomial<T>& iy = integralY(y0, y1);  // still a function of x
  return iy.definiteIntegral(x0, x1);       // just a number
}
// make function names consistent maybe rename definiteIntegral to integral...but this could lead
// to ambiguities with the function that evaluates the indefinite integral at some x0 with 
// integration constant c

template<class T> 
T rsBivariatePolynomial<T>::pathIntegral(const rsBivariatePolynomial<T>& p, 
  const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b)
{
  rsAssert(x.getDegree() <= 1 && y.getDegree() <= 1, "Only implemented for linear paths");
  // ToDo: lift this restriction by switching to numerical integration, if any degree is larger.
  // When x(t) or y(t) is nonlinear, then xp or yp would still be functions of t and the sqrt would
  // appear in the integrand and we would have to evaluate x,y at t within the integral.

  using Poly   = rsPolynomial<T>;
  using BiPoly = rsBivariatePolynomial<T>;
  Poly pt = BiPoly::compose(p, x, y);     // p(t)
  Poly xp = x.derivative();               // x'(t) - optimize: it's just the coeff for x^1
  Poly yp = y.derivative();               // y'(t)   ...but take care: the degree may be zero
  T dx = xp(T(0));                        // dx/dt, constant bcs x is linear
  T dy = yp(T(0));                        // dy/dt, constant bcs y is linear
  T ds = sqrt(dx*dx + dy*dy);             // arc length (the "arc" is just a line element)
  return ds * pt.definiteIntegral(a, b);
}
// pt has 2 zero coeffs at the end in the test - check, why that happens - should we expect that or
// is that a bug?

template<class T> 
T rsBivariatePolynomial<T>::pathIntegral(
  const rsBivariatePolynomial<T>& u, const rsBivariatePolynomial<T>& v,
  const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b)
{
  using Poly   = rsPolynomial<T>;
  using BiPoly = rsBivariatePolynomial<T>;
  Poly ut = BiPoly::compose(u, x, y);  // u(t) - todo: support syntax: u.compose(x, y)
  Poly vt = BiPoly::compose(v, x, y);  // v(t)
  Poly xp = x.derivative();            // x'(t)
  Poly yp = y.derivative();            // y'(t)
  Poly P  = ut * xp + vt * yp;         // the scalar product in the integrand
  return P.definiteIntegral(a, b);
}
// todo: maybe have a function that leaves a,b, unspecified - this should return a single 
// univariate polynomial into which the limits can be inserted later, i.e. just return
// P.integral or p.indefiniteIntegral

template<class T> 
T rsBivariatePolynomial<T>::fluxIntegral(
  const rsBivariatePolynomial<T>& p, const rsBivariatePolynomial<T>& q,
  const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b)
{
  using Poly   = rsPolynomial<T>;
  using BiPoly = rsBivariatePolynomial<T>;
  Poly pt = BiPoly::compose(p, x, y);  // p(t)
  Poly qt = BiPoly::compose(q, x, y);  // q(t)
  Poly xp = x.derivative();            // x'(t)
  Poly yp = y.derivative();            // y'(t)
  Poly P  = pt * yp - qt * xp;
  return P.definiteIntegral(a, b);
}
// needs tests - this is on shaky grounds - especially with respect to computing and 
// (not) normalizing the normal vector (which is taken to be equal to (nx, ny) = (yp, -xp) here...
// but actually should be normalized...right?). ...in general, the function is very similar to
// pathIntegral with P = ut * xp + vt * yp  replaced by  P = pt * yp - qt * xp

template<class T> 
T rsBivariatePolynomial<T>::loopIntegral(const rsBivariatePolynomial<T>& p,
  const rsBivariatePolynomial<T>& q, T x0, T x1, T y0, T y1)
{
  using P = rsPolynomial<T>;
  P xt, yt;
  T r = T(0);
  xt = P({x0, x1-x0}); yt = P({y0}); r += pathIntegral(p, q, xt, yt, T(0), T(1));  // rightward
  xt = P({x1}); yt = P({y0, y1-y0}); r += pathIntegral(p, q, xt, yt, T(0), T(1));  // upward
  xt = P({x1, x0-x1}); yt = P({y1}); r += pathIntegral(p, q, xt, yt, T(0), T(1));  // leftward
  xt = P({x0}); yt = P({y1, y0-y1}); r += pathIntegral(p, q, xt, yt, T(0), T(1));  // downward
  return r;
}
// optimize: avoid creating so many temporary univariate polynomials xt, yt - reuse the existing 
// objects by re-assigning the coeffs

template<class T> 
T rsBivariatePolynomial<T>::outfluxIntegral(const rsBivariatePolynomial<T>& p, 
  const rsBivariatePolynomial<T>& q, T x0, T x1, T y0, T y1)
{
  using P = rsPolynomial<T>;
  P xt, yt;
  T r = T(0);
  xt = P({x0, x1-x0}); yt = P({y0}); r += fluxIntegral(p, q, xt, yt, T(0), T(1));  // rightward
  xt = P({x1}); yt = P({y0, y1-y0}); r += fluxIntegral(p, q, xt, yt, T(0), T(1));  // upward
  xt = P({x1, x0-x1}); yt = P({y1}); r += fluxIntegral(p, q, xt, yt, T(0), T(1));  // leftward
  xt = P({x0}); yt = P({y1, y0-y1}); r += fluxIntegral(p, q, xt, yt, T(0), T(1));  // downward
  return r;
}
// maybe get rid of the duplication by using a pointer to fluxIntegral/pathIntegral

template<class T>
bool rsBivariatePolynomial<T>::areHarmonicConjugates(
  const rsBivariatePolynomial<T>& u, const rsBivariatePolynomial<T>& v, T tol)
{
  // We must have: u_x == v_y and u_y == -v_x:
  using BiPoly = rsBivariatePolynomial<T>;
  bool r = true;
  BiPoly ux = u.derivativeX();
  BiPoly uy = u.derivativeY();
  BiPoly vx = v.derivativeX(); vx.negate();
  BiPoly vy = v.derivativeY();
  r &= ux.isCloseTo(vy, tol);
  r &= uy.isCloseTo(vx, tol);
  return r;
}
// try to optimize: avoid creating the temporary BiPolys. Instead, compare coeffs in u, v directly.

// optimized version of composeWithLinear
template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::composeWithLinear(
  const rsPolynomial<T>& p, T a, T b)
{
  int N = p.getDegree();
  rsBivariatePolynomial<T> r(N, N);
  r.coeffs.setToZero();
  const T* c = p.getCoeffPointerConst();
  int N1 = N+1;
  std::vector<T> wrk(3*N1);             // workspace
  T* B  = &wrk[0*N1];                   // binomial coeffs
  T* an = &wrk[1*N1];                   // powers of a
  T* bn = &wrk[2*N1];                   // powers of b
  an[0] = T(1);                         // a^0
  bn[0] = T(1);                         // b^0
  for(int n = 1; n <= N; n++) {
    an[n] = a * an[n-1];                // a^n
    bn[n] = b * bn[n-1]; }              // b^n
  for(int n = 0; n <= N; n++) {
    rsNextPascalTriangleLine(B, B, n);  // B[k] is now B(n,k) the binomial coeff "n-choose-k"
    for(int k = 0; k <= n; k++)
      r.coeffs(k, n-k) += c[n] * B[k] * an[k] * bn[n-k]; } // += c[n] * B(n,k) * a^k * b^(n-k)
  return r;
}
// todo:
// -let the workspace be passed by the user
// -operate on a pre-allocated rsMatrixView
// -the polynomial p should be passed as raw coefficient array
// -maybe generalize to composeWithAffine(..., T a, T b, T c) that computes the composition with 
//  (a*x + b*y + c) - i think, we need trinomial coeffients and Pascal's pyramid for this:
//  https://en.wikipedia.org/wiki/Multinomial_theorem, 
//  https://en.wikipedia.org/wiki/Pascal%27s_pyramid

// old version of composeWithLinear, just for reference
template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::composeWithLinearOld(
  const rsPolynomial<T>& p, T a, T b)
{
  int N = p.getDegree();
  rsBivariatePolynomial<T> r(N, N);
  r.coeffs.setToZero();
  const T* c = p.getCoeffPointerConst();
  for(int n = 0; n <= N; n++)
    for(int k = 0; k <= n; k++)
      r.coeffs(k, n-k) += c[n] * (T) rsBinomialCoefficient(n, k) * pow(a, k) * pow(b, n-k);
  return r;
}
// This works, but is very inefficient: "Shlemiel the painter" strikes again in 
// rsBinomialCoefficient and pow is expensive! On the other hand, it does not allocate memory.
// ....maybe keep both versions...or well.. it actually does allocate for the returned object
// 

template<class T> 
rsPolynomial<T> rsBivariatePolynomial<T>::compose(const rsBivariatePolynomial<T>& p,
  const rsPolynomial<T>& x, const rsPolynomial<T>& y)
{
  //   p(x,y) = \sum_m \sum_n a_{mn} x^m y^n
  // where: 
  //   x = x(t) = \sum_i b_i t^i
  //   y = y(t) = \sum_j c_j t^j
  // so the univariate p(t) is:
  //   p(t) = \sum_m \sum_n \left(  a_{mn} * (\sum_i b_i t^i)^m * (\sum_j c_j t^j)^n  \right)
  // so, the algo needs:
  //   -successive powers of the b,c arrays (iterated convolutions of the arrays with themselves)
  //   -the products of all possible combinations of those powers (more convolutions)
  //   -multiply these products by a coeff a_{mn} and accumulate the result into the output coeff
  //    array
  using Poly = rsPolynomial<T>;
  using AT   = rsArrayTools;
  int M = p.getDegreeX();
  int N = p.getDegreeY();
  int I = x.getDegree();
  int J = y.getDegree();
  int L = I+M + J+N;              // degree of result
  int strideX = M*(I+1)-1;        // number of columns in matrix xp (powers of x)
  int strideY = N*(J+1)-1;        // same for yp
  rsMatrix<T> xp(M+1, strideX);   // powers of x
  rsMatrix<T> yp(N+1, strideY);   // powers of y
  std::vector<T> tmp(L+1);        // holds coeffs for (x(t))^m * (y(t))^n as function of t
  Poly pt(L);                     // target univariate polynomial p(t)
  Poly::powers(x.getCoeffPointerConst(), I, xp.getDataPointer(), M, strideX);
  Poly::powers(y.getCoeffPointerConst(), J, yp.getDataPointer(), N, strideY);
  for(int m = 0; m <= M; m++) {
    for(int n = 0; n <= N; n++) {
      T*  xm  = xp.getRowPointer(m);            // coeffs for (x(t))^m
      T*  yn  = yp.getRowPointer(n);            // coeffs for (y(t))^n
      int lxm = m*I+1;                          // effective length of xm
      int lyn = n*J+1;                          // effective length of xm
      AT::convolve(xm, lxm, yn, lyn, &tmp[0]);
      for(int i = 0; i < lxm+lyn-1; i++)
        pt[i] += p.coeff(m, n) * tmp[i];  }}    // accumulation
  return pt;
}
// what about the composition where the inner polynomial is bivariate and the outer univariate?


template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::multiplyY(const rsPolynomial<T>& polyY) const
{
  int degX = getDegreeX();
  int degY = getDegreeY();
  int degP = polyY.getDegree();
  rsBivariatePolynomial<T> r(degX, degY+degP);
  const T* h = polyY.getCoeffPointerConst();     // "impulse response"
  for(int i = 0; i <= degX; i++) {
    const T* x = coeffs.getRowPointerConst(i);   // input
    T* y = r.coeffs.getRowPointer(i);            // output
    rsArrayTools::convolve(x, degY+1, h, degP+1, y); }
  return r;
}
// this assumes row-major storage of matrices

template<class T>
void rsBivariatePolynomial<T>::splitRealImag(const rsBivariatePolynomial<std::complex<T>>& p,
  rsBivariatePolynomial<T>& pRe, rsBivariatePolynomial<T>& pIm)
{
  int m = p.getDegreeX();
  int n = p.getDegreeY();
  pRe.initialize(m, n);
  pIm.initialize(m, n);
  for(int i = 0; i <= m; i++) {
    for(int j = 0; j <= n; j++) {
      pRe.coeff(i, j) = p.coeff(i, j).real();
      pIm.coeff(i, j) = p.coeff(i, j).imag(); }}
}

template<class T>
void rsBivariatePolynomial<T>::polyaVectorField(const rsPolynomial<std::complex<T>>& p,
  rsBivariatePolynomial<T>& px, rsBivariatePolynomial<T>& py)
{
  using Complex = std::complex<T>;
  using BiPolyC = rsBivariatePolynomial<Complex>;
  Complex one(1, 0), im(0, 1);
  BiPolyC bp = BiPolyC::composeWithLinear(p, one, im); // bp(x, y) = p(x + i*y) = p(z)
  splitRealImag(bp, px, py);                           // extract real and imaginary parts
  py.negate();                                         // apply complex conjugation
}

// ToDo:
// -compute Hessian and its determinant and maybe (squared) eigenvalues and -vectors
// -root finding: find points (x,y) for which p(x,y) = 0. i think, these are in general not 
//  isolated points but rather curves - for example, when p(x,y) = x^2 + y^2 - 1, the unit circle
//  gives the set of zeros and there we have y^2 = sqrt(1- x^2)...not sure how to deal with this
// -maybe implement constrained optimization via Lagarange multipliers - find extremum of p(x,y)
//  subject to c(x,y) = 0 where c is the constraint - the Lagrange function will be a trivariate 
//  polynomial
// -implement differentiation of implictit function (Kaprfinger, pg. 560)
// -implement differential operators in polar-coordinates (needs bivariate rational function for
//  the 1/r factor)
// -can we compute path integrals along circular or elliptic arcs? their parametric description 
//  involves sin/cos but their implicit description is polynomial. the explicit description as
//  y = f(x) involves the sqrt - but maybe we can somehow evade this by integrating not with 
//  respect to x but with respect to x^2 - maybe this can help:
//  https://en.wikipedia.org/wiki/Riemann%E2%80%93Stieltjes_integral
//  -figure out first, if the result is actually a polynomial (when one or both limits are left as
//   free parameters)...well...nope...i think, that can't be the case - arc lengths are in general
//   transcendental numbers

