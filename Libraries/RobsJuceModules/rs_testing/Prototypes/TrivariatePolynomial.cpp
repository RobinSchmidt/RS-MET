
// inquiry:

template<class T>
bool rsTrivariatePolynomial<T>::isCloseTo(const rsTrivariatePolynomial<T>& p, T tol) const
{
  int L = rsMax(getDegreeX(), p.getDegreeX());
  int M = rsMax(getDegreeY(), p.getDegreeY());
  int N = rsMax(getDegreeZ(), p.getDegreeZ());
  for(int i = 0; i <= L; i++)
    for(int j = 0; j <= M; j++)
      for(int k = 0; k <= N; k++)
        if( rsAbs( getCoeffPadded(i,j,k) - p.getCoeffPadded(i,j,k)) > tol )
          return false;
  return true;
}

template<class T>
bool rsTrivariatePolynomial<T>::isVectorPotential(const rsTrivariatePolynomial<T>& F,
  const rsTrivariatePolynomial<T>& G, const rsTrivariatePolynomial<T>& H,
  const rsTrivariatePolynomial<T>& f, const rsTrivariatePolynomial<T>& g,
  const rsTrivariatePolynomial<T>& h, T tol)
{
  rsTrivariatePolynomial<T> cx, cy, cz; // components of curl of F,G,H
  curl(F, G, H, cx, cy, cz);
  return f.isCloseTo(cx, tol) && g.isCloseTo(cy, tol) && h.isCloseTo(cz, tol);
}

// evaluation:

template<class T>
T rsTrivariatePolynomial<T>::evaluate(T x, T y, T z) const
{
  T xl(1), ym(1), zn(1), r(0);  // x^l, y^m, z^n, result
  for(int l = 0; l < coeffs.getExtent(0); l++) {
    ym = T(1);
    for(int m = 0; m < coeffs.getExtent(1); m++) {
      zn = T(1);
      for(int n = 0; n < coeffs.getExtent(2); n++) {
        r += coeffs(l, m, n) * xl * ym * zn;
        zn *= z; }
      ym *= y; }
    xl *= x; }
  return r;
}
// have functions where x,y,z are uni- or bivariate polynomials - the result is then again a 
// polynomial (of the same kind)

template<class T>
rsBivariatePolynomial<T> rsTrivariatePolynomial<T>::evaluateX(T x) const
{
  int L = getDegreeX();
  int M = getDegreeY();
  int N = getDegreeZ();
  rsBivariatePolynomial<T> p_yz(M, N);
  T xl(1);   // x^l
  for(int l = 0; l <= L; l++) {
    for(int m = 0; m <= M; m++) {
      for(int n = 0; n <= N; n++) {
        p_yz.coeff(m, n) += coeffs(l, m, n) * xl; }}
    xl *= x; }
  return p_yz;
}

// arithmetic:

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::operator*(
  const rsTrivariatePolynomial<T>& p) const
{
  rsTrivariatePolynomial<T> r;
  rsMultiArray<T>::convolve3D(coeffs, p.coeffs, r.coeffs);
  return r;
}

template<class T> 
rsPolynomial<T> rsTrivariatePolynomial<T>::compose(const rsTrivariatePolynomial<T>& p,
  const rsPolynomial<T>& x, const rsPolynomial<T>& y, const rsPolynomial<T>& z)
{
  using Poly = rsPolynomial<T>;
  Poly p_t;            // p(t) = p(x(t), y(t), z(t))
  //Poly one({1});       // constant polynomial one(t) = 1
  Poly one(0); one[0] = T(1);
  Poly xl, ym, zn;     // (x(t))^l, (y(t))^m, (z(t))^n
  xl = one;
  for(int l = 0; l < p.coeffs.getExtent(0); l++) {
    ym = one;
    for(int m = 0; m < p.coeffs.getExtent(1); m++) {
      zn = one;
      for(int n = 0; n < p.coeffs.getExtent(2); n++) {
        p_t = p_t + p.coeffs(l, m, n) * xl * ym * zn;
        zn = zn * z; }
      ym = ym * y; }
    xl = xl * x; }
  return p_t;
}
// try to avoid the duplication with the function below by templatizing - Poly and BiPoly both need
// a factory function to create the one-polynomial or maybe more generally a constant polynomial
// ...but maybe we later optimize in which case the code will look different in both cases

template<class T> 
rsBivariatePolynomial<T> rsTrivariatePolynomial<T>::compose(const rsTrivariatePolynomial<T>& p,
  const rsBivariatePolynomial<T>& x, const rsBivariatePolynomial<T>& y,
  const rsBivariatePolynomial<T>& z)
{
  //   p(x,y,z) = \sum_l \sum_m \sum_n a_{lmn} x^l y^m z^n
  // where: 
  //   x = x(u,v) = \sum_h \sum_i b_{hi} u^h v^i
  //   y = y(u,v) = \sum_j \sum_k c_{jk} u^j v^k
  //   z = z(u,v) = \sum_q \sum_r d_{qr} u^q v^r
  // so the bivariate p(u,v) is:
  //   p(u,v) = \sum_l \sum_m \sum_n  
  //            \left(   a_{lmn} 
  //                   * (\sum_h \sum_i b_{hi} u^h v^i)^l
  //                   * (\sum_j \sum_k c_{jk} u^j v^k)^m
  //                   * (\sum_q \sum_r d_{qr} u^q v^r)^n
  //            \right)
  // so, the algo needs:
  //   -successive powers of the b,c,d matrices (iterated 2D convolutions of the matrices with 
  //    themselves)
  //   -the products of all possible combinations of those powers (more 2D convolutions)
  //   -multiply these products by a coeff a_{lmn} and accumulate the result into the output coeff
  //    array
  //   -that's complicated - let's write an inefficient prototype first that just works like 
  //    evaluation but with bivariate polynomials instead of numbers


  using BiPoly  = rsBivariatePolynomial<T>;
  BiPoly p_uv;            // p(u,v) = p(x(u,v), y(u,v), z(u,v))
  BiPoly one(0, 0, {1});  // constant bivariate polynomial one(u,v) = 1
  BiPoly xl, ym, zn;      // (x(u,v))^l, (y(u,v))^m, (z(u,v))^n

  // This is quite inefficient (lots of temporary objects are created where we could potentially 
  // work in place) - but it's readable - may be optimized later:
  xl = one;
  for(int l = 0; l < p.coeffs.getExtent(0); l++) {
    ym = one;
    for(int m = 0; m < p.coeffs.getExtent(1); m++) {
      zn = one;
      for(int n = 0; n < p.coeffs.getExtent(2); n++) {
        p_uv = p_uv + p.coeffs(l, m, n) * xl * ym * zn;
        zn = zn * z; }
      ym = ym * y; }
    xl = xl * x; }
  return p_uv;
}

// calculus:

template<class T>
void rsTrivariatePolynomial<T>::derivativeX(const rsMultiArray<T>& c, rsMultiArray<T>& d)
{
  int L = c.getExtent(0);
  int M = c.getExtent(1);
  int N = c.getExtent(2);
  rsAssert(d.hasShape({ L-1, M, N }));
  for(int l = 1; l < L; l++) {
    T s(l);
    for(int m = 0; m < M; m++)
      for(int n = 0; n < N; n++)
        d(l-1, m, n) = s * c(l, m, n); }
}
// maybe don't require d to have the proper shape - instead, set up the shape of d

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::derivativeX() const
{
  rsTrivariatePolynomial<T> q(getDegreeX()-1, getDegreeY(), getDegreeZ());
  derivativeX(coeffs, q.coeffs);
  return q;
}

template<class T>
void rsTrivariatePolynomial<T>::derivativeY(const rsMultiArray<T>& c, rsMultiArray<T>& d)
{
  int L = c.getExtent(0);
  int M = c.getExtent(1);
  int N = c.getExtent(2);
  rsAssert(d.hasShape({ L, M-1, N }));
  for(int m = 1; m < M; m++) {
    T s(m);
    for(int l = 0; l < L; l++)
      for(int n = 0; n < N; n++)
        d(l, m-1, n) = s * c(l, m, n); }
}

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::derivativeY() const
{
  rsTrivariatePolynomial<T> q(getDegreeX(), getDegreeY()-1, getDegreeZ());
  derivativeY(coeffs, q.coeffs);
  return q;
}

template<class T>
void rsTrivariatePolynomial<T>::derivativeZ(const rsMultiArray<T>& c, rsMultiArray<T>& d)
{
  int L = c.getExtent(0);
  int M = c.getExtent(1);
  int N = c.getExtent(2);
  rsAssert(d.hasShape({ L, M, N-1 }));
  for(int n = 1; n < N; n++) {
    T s(n);
    for(int m = 0; m < M; m++)
      for(int l = 0; l < L; l++)
        d(l, m, n-1) = s * c(l, m, n); }
}

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::derivativeZ() const
{
  rsTrivariatePolynomial<T> q(getDegreeX(), getDegreeY(), getDegreeZ()-1);
  derivativeZ(coeffs, q.coeffs);
  return q;
}

template<class T>
void rsTrivariatePolynomial<T>::curl(const rsTrivariatePolynomial<T>& fx, 
  const rsTrivariatePolynomial<T>& fy, const rsTrivariatePolynomial<T>& fz, 
  rsTrivariatePolynomial<T>& cx, rsTrivariatePolynomial<T>& cy, rsTrivariatePolynomial<T>& cz)
{
  cx = fz.derivativeY() - fy.derivativeZ();
  cy = fx.derivativeZ() - fz.derivativeX();
  cz = fy.derivativeX() - fx.derivativeY();
}

template<class T>
void rsTrivariatePolynomial<T>::integralX(const rsMultiArray<T>& a, rsMultiArray<T>& ai, T c)
{
  int L = a.getExtent(0);
  int M = a.getExtent(1);
  int N = a.getExtent(2);
  rsAssert(ai.hasShape({ L+1, M, N }));

  ai(0, 0, 0) = c; 
  for(int m = 1; m < M; m++)
    for(int n = 1; n < N; n++)
      ai(0, m, n) = 0; 
  // i think, more generally, the integration "constant" c could be a bivariate polynomial in y,z 
  // and we would do: ai(0, m, n) = c(m, n) for (m,n) = (0,0)...(M-1,N-1)

  for(int l = 1; l <= L; l++) {
    T s = T(1) / T(l);
    for(int m = 0; m < M; m++)
      for(int n = 0; n < N; n++)
        ai(l, m, n) = s * a(l-1, m, n); }
}

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::integralX(T c) const
{
  rsTrivariatePolynomial<T> q(getDegreeX()+1, getDegreeY(), getDegreeZ());
  integralX(coeffs, q.coeffs, c);
  return q;
}
// needs test

template<class T>
void rsTrivariatePolynomial<T>::integralY(const rsMultiArray<T>& a, rsMultiArray<T>& ai, T c)
{
  int L = a.getExtent(0);
  int M = a.getExtent(1);
  int N = a.getExtent(2);
  rsAssert(ai.hasShape({ L, M+1, N }));
  ai(0, 0, 0) = c; 
  for(int l = 1; l < L; l++)
    for(int n = 1; n < N; n++)
      ai(l, 0, n) = 0; 
  for(int m = 1; m <= M; m++) {
    T s = T(1) / T(m);
    for(int l = 0; l < L; l++)
      for(int n = 0; n < N; n++)
        ai(l, m, n) = s * a(l, m-1, n); }
}

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::integralY(T c) const
{
  rsTrivariatePolynomial<T> q(getDegreeX(), getDegreeY()+1, getDegreeZ());
  integralY(coeffs, q.coeffs, c);
  return q;
}
// needs test

template<class T>
void rsTrivariatePolynomial<T>::integralZ(const rsMultiArray<T>& a, rsMultiArray<T>& ai, T c)
{
  int L = a.getExtent(0);
  int M = a.getExtent(1);
  int N = a.getExtent(2);
  rsAssert(ai.hasShape({ L, M, N+1 }));
  ai(0, 0, 0) = c; 
  for(int l = 1; l < L; l++)
    for(int m = 1; m < M; m++)
      ai(l, m, 0) = 0; 
  for(int n = 1; n <= N; n++) {
    T s = T(1) / T(n);
    for(int l = 0; l < L; l++)
      for(int m = 0; m < M; m++)
        ai(l, m, n) = s * a(l, m, n-1); }
}

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::integralZ(T c) const
{
  rsTrivariatePolynomial<T> q(getDegreeX(), getDegreeY(), getDegreeZ()+1);
  integralZ(coeffs, q.coeffs, c);
  return q;
}

template<class T>
template<class Ta, class Tb>
rsBivariatePolynomial<T> rsTrivariatePolynomial<T>::integralX(Ta a, Tb b) const
{
  rsTrivariatePolynomial<T> P = integralX();
  rsBivariatePolynomial<T> Pb = P.evaluateX(b);
  rsBivariatePolynomial<T> Pa = P.evaluateX(a);
  return Pb - Pa;
}

template<class T> 
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::scalarPotential(
  const rsTrivariatePolynomial<T>& px, const rsTrivariatePolynomial<T>& py, 
  const rsTrivariatePolynomial<T>& pz)
{
  rsTrivariatePolynomial<T> Px, Px_y, gyz_y, gyz, Pxy, Pxy_z, hz_z, hz, Pxyz;
  Px    = px.integralX();      // integrate px with respect to x

  Px_y  = Px.derivativeY();    // differentiate the result with respect to y
  gyz_y = py - Px_y;           // d g(y,z) / dy
  gyz   = gyz_y.integralY();   // g(y,z)
  Pxy   = Px + gyz;

  Pxy_z = Pxy.derivativeZ();
  hz_z  = pz - Pxy_z;          // d h(z) / dz
  hz    = hz_z.integralZ();    // h(z)
  Pxyz  = Pxy + hz;

  return Pxyz;
}
// maybe use rsAssert(hasScalarPotential(px, py, pz)); ..but with tolerance

template<class T> 
void rsTrivariatePolynomial<T>::vectorPotential(const rsTrivariatePolynomial<T>& f,
  const rsTrivariatePolynomial<T>& g, const rsTrivariatePolynomial<T>& h,
  rsTrivariatePolynomial<T>& F, rsTrivariatePolynomial<T>& G, rsTrivariatePolynomial<T>& H)
{
  //rsAssert(hasVectorPotential(f, g, h));  // we need a tolerance

  using TP = rsTrivariatePolynomial<T>;
  TP fz   = f.integralZ();
  TP gz   = g.integralZ();
  TP fz_x = fz.derivativeX();
  TP gz_y = gz.derivativeY();
  TP a_x  = h + fz_x + gz_y;
  TP a    = a_x.integralX();

  F = gz;           // F(x,y,z) =  gz + b(x,y) where b(x,y) = 0
  G = a - fz;       // G(x,y,z) = -fz + a(x,y)
  H = TP(0,0,0);    // H(x,y,z) =  0
}
// -see VectorPotentials.txt in the Notes folder for derivation of the algo
// -maybe get rid of the parameter H - client code may implicitly assume it to be zero
// -maybe give the user more choices with respect to the free choices that can be made, such as
//  setting H(x,y,z) = 0, b(x,y) = 0, see:
//  https://en.wikipedia.org/wiki/Gauge_fixing
//  https://en.wikipedia.org/wiki/Magnetic_vector_potential#Gauge_choices
//  https://www.reed.edu/physics/faculty/wheeler/documents/Electrodynamics/Class%20Notes/Chapter%204.pdf


template<class T> 
T rsTrivariatePolynomial<T>::pathIntegral(
  const rsTrivariatePolynomial<T>& fx, 
  const rsTrivariatePolynomial<T>& fy,
  const rsTrivariatePolynomial<T>& fz, 
  const rsPolynomial<T>& x, 
  const rsPolynomial<T>& y, 
  const rsPolynomial<T>& z, 
  T a, T b)
{
  using Poly    = rsPolynomial<T>;
  using TriPoly = rsTrivariatePolynomial<T>;
  Poly fxt = TriPoly::compose(fx, x, y, z);   // fx(t) = fx(x(t),y(t),z(t))
  Poly fyt = TriPoly::compose(fy, x, y, z);   // fy(t) = fy(x(t),y(t),z(t))
  Poly fzt = TriPoly::compose(fz, x, y, z);   // fz(t) = fz(x(t),y(t),z(t))
  Poly xp = x.derivative();                   // x'(t)
  Poly yp = y.derivative();                   // y'(t)
  Poly zp = z.derivative();                   // z'(t)
  Poly P  = fxt * xp + fyt * yp + fzt * zp;   // the scalar product in the integrand
  return P.definiteIntegral(a, b);
}

template<class T> 
T rsTrivariatePolynomial<T>::pathIntegral(
  const rsTrivariatePolynomial<T>& fx, 
  const rsTrivariatePolynomial<T>& fy,
  const rsTrivariatePolynomial<T>& fz,
  const std::vector<rsVector3D<T>>& path)
{
  T result = T(0);
  rsPolynomial<T> x(1), y(1), z(1);
  for(size_t i = 1; i < path.size(); i++) {
    x[0] = path[i-1].x; x[1] = path[i].x - x[0];
    y[0] = path[i-1].y; y[1] = path[i].y - y[0];
    z[0] = path[i-1].z; z[1] = path[i].z - z[0];
    result += pathIntegral(fx, fy, fz, x, y, z, T(0), T(1)); }
  return result;
}
// needs more tests

template<class T> 
T rsTrivariatePolynomial<T>::fluxIntegral(
  const rsTrivariatePolynomial<T>& fx,
  const rsTrivariatePolynomial<T>& fy,
  const rsTrivariatePolynomial<T>& fz,
  const rsBivariatePolynomial<T>& x,
  const rsBivariatePolynomial<T>& y,
  const rsBivariatePolynomial<T>& z,
  T u0, T u1, T v0, T v1)
{
  using BiPoly  = rsBivariatePolynomial<T>;
  using TriPoly = rsTrivariatePolynomial<T>;

  // vector field on the surface:
  BiPoly gx = TriPoly::compose(fx, x, y, z);  // gx(u,v) = fx(x(u,v), y(u,v), z(u,v))
  BiPoly gy = TriPoly::compose(fy, x, y, z);  // gy(u,v) = fy(x(u,v), y(u,v), z(u,v))
  BiPoly gz = TriPoly::compose(fz, x, y, z);  // gz(u,v) = fz(x(u,v), y(u,v), z(u,v))

  // partial derivatives of gx,gy,gz with respect to u,v:
  BiPoly ax = x.derivativeX();                // ax := d x(u,v) / du
  BiPoly ay = y.derivativeX();                // ay := d y(u,v) / du
  BiPoly az = z.derivativeX();                // az := d z(u,v) / du
  BiPoly bx = x.derivativeY();                // bx := d x(u,v) / dv
  BiPoly by = y.derivativeY();                // by := d y(u,v) / dv
  BiPoly bz = z.derivativeY();                // bz := d z(u,v) / dv 

  // components of the cross-product (Bärwolff, pg 355):
  BiPoly cx = ay*bz - az*by;
  BiPoly cy = az*bx - ax*bz;
  BiPoly cz = ax*by - ay*bx;

  // differential flux element and total flux through surface (Bärwollf, pg 600):
  BiPoly df = gx*cx + gy*cy + gz*cz;
  return df.doubleIntegralXY(u0, u1, v0, v1);
}

template<class T> 
T rsTrivariatePolynomial<T>::fluxIntegral(
  const rsTrivariatePolynomial<T>& fx,
  const rsTrivariatePolynomial<T>& fy,
  const rsTrivariatePolynomial<T>& fz,
  const rsVector3D<T>& P0, const rsVector3D<T>& P1, const rsVector3D<T>& P2)
{
  rsBivariatePolynomial<T> x(1,1), y(1,1), z(1,1);
  x.coeff(0, 0) = P0.x;
  y.coeff(0, 0) = P0.y;
  z.coeff(0, 0) = P0.z;

  x.coeff(1, 0) = P1.x - P0.x;
  y.coeff(1, 0) = P1.y - P0.y;
  z.coeff(1, 0) = P1.z - P0.z;

  x.coeff(1, 1) = P2.x - P1.x;
  y.coeff(1, 1) = P2.y - P1.y;
  z.coeff(1, 1) = P2.z - P1.z;

  return fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // Parametrization was derived by considering:
  //   p1(u) = P0 + u*(P1-P0), p2(u) = P0 + u*(P2-P0)
  // which leads to:
  //   p(u,v) = (1-v)*p1(u) + v*p2(u)  
  //          = P0 + (P1-P0)*u + (P2-P1)*u*v
}

template<class T> 
T rsTrivariatePolynomial<T>::outfluxIntegral(const rsTrivariatePolynomial<T>& fx,
  const rsTrivariatePolynomial<T>& fy, const rsTrivariatePolynomial<T>& fz,
  T x0, T x1, T y0, T y1, T z0, T z1)
{
  using BiPoly = rsBivariatePolynomial<T>;

  T dx = x1 - x0;
  T dy = y1 - y0;
  T dz = z1 - z0;
  BiPoly x, y, z;

  // flux through surface patch where z = z0:
  x = BiPoly(1, 1, { x0,   0, dx, 0 });
  y = BiPoly(1, 1, { y1, -dy,  0, 0 });
  z = BiPoly(0, 0, { z0             });
  T fz0 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // flux through surface patch where z = z1:
  //x = BiPoly(1, 1, { x0,  0, dx, 0  }); // has not changed
  y = BiPoly(1, 1, { y0, dy,  0, 0  });
  z = BiPoly(0, 0, { z1             });
  T fz1 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // flux through surface patch where y = y0:
  //x = BiPoly(1, 1, { x0, 0,  dx, 0 }); // has not changed
  y = BiPoly(0, 0, { y0            });
  z = BiPoly(1, 1, { z0,  dz, 0, 0 });
  T fy0 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // flux through surface patch where y = y1:
  //x = BiPoly(1, 1, { x0, 0,  dx, 0 }); // has not changed
  y = BiPoly(0, 0, { y1            });
  z = BiPoly(1, 1, { z1, -dz, 0, 0 });
  T fy1 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // flux through surface patch where x = x0:
  x = BiPoly(0, 0, { x0           });
  y = BiPoly(1, 1, { y0,  0, dy, 0 });
  z = BiPoly(1, 1, { z1, -dz, 0, 0 });
  T fx0 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // flux through surface patch where x = x1:
  x = BiPoly(0, 0, { x1           });
  //y = BiPoly(1, 1, { y0, 0, dy, 0 }); // has not changed
  z = BiPoly(1, 1, { z0, dz, 0, 0 });
  T fx1 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  return fx0 + fx1 + fy0 + fy1 + fz0 + fz1;
}
// -optimize 
// -create the BiPoly objects only once and re-assign the coeffs -> less allocation
// -remove the fz0, etc variables - directly add intermediate results into an accumulator

// ToDo:
// -path and surface integrals over scalar fields, restricted to paths and surfaces parametrized by 
//  (bi)linear functions, such that the factor for the differential path- or surface element 
//  becomes a constant (only in this case, we can evade the annoying square-root in the integrand). 
//  ..but what would be the interpretation of path integrals over scalar fields in 3D - does it 
//  even make sense? if not, implement only surface integrals. i think, it would somehow be the 
//  amount of 3D space between the (x,y,z)-space and the swept out line - whatever that means - by 
//  analogy with "curtain" interpretation of 2D path integrals over 2D scalar fields. the 
//  interpretation of surface integrals is: if the scalar field gives a mass-density per unit area,
//  the integral computes the total mass of a thin shell...where we are restricted to planar 
//  "shells" due to the pesky square-root problem - which is perhaps not very interesting, but for
//  completeness
// -transformed (volume) integrals (using the determinant of the Jacobian matrix) - needs perhaps 
//  also be restricted to linear transformations because of the inverse transformations that occurs
//  in the adjustment of the integration limits - (i think) we would then be able to integrate over 
//  general parallelepipeds instead of just axis-aligned cuboids
//  https://www.youtube.com/watch?v=geJ-36mnZ1I
// -constrained optimization via Lagrange multipliers (maybe - but how?)
// -Legendre transform (maybe) - nope, not possible because it involves inverse functions
// -for vector potentials:
//  -let user choose a divergence
//   -we need a function to create a gradient field with a given divergence, so we need a field 
//    whose Laplacian is zero and whose divergence is as given
//   -then, compute the actual divergence of the simple vector potential, obtain the difference 
//    with the desired divergence and add the gradient field that has as divergence this difference
//   -if divergence is zero, the vector potential itself has a vector potential
//  -can we choose the vector potential such that its Laplacian is zero? this would imply that we 
//   can find a scalar potential for it. ...probably not - because the curl of a gradient is always 
//   zero, so if we take the gradient of a scalar field (to get the vector potential) and then take 
//   the curl of that to get our vector field, we would always get zero