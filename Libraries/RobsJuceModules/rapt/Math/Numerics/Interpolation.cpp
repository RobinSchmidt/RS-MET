template<class Tx, class Ty>
void resampleNonUniformLinear(const Tx* xIn, const Ty* yIn, int inLength, 
  const Tx* xOut, Ty* yOut, int outLength)
{
  // code copied and adapted from from rsTimeWarper::invertMonotonousWarpMap
  int i = 0;
  for(int n = 0; n < outLength; n++) {
    Tx x = xOut[n];
    while(xIn[i+1] < x && i < inLength-2)
      i++;
    yOut[n] = rsInterpolateLinear(xIn[i], xIn[i+1], yIn[i], yIn[i+1], x);
  }
}
// optimize: recompute line-equation coeffs only when necessarry - avoid 
// repeated recomputation of still valid line-coeffs in case of upsampling

template<class Tx, class Ty>
void rsNaturalCubicSpline(const Tx *xIn, const Ty *yIn, int N, const Tx *xOut, Ty *yOut, 
  int Ni, Ty scaleRhs)
{
  if(N <= 3) {
    resampleNonUniformLinear(xIn, yIn, N, xOut, yOut, Ni);
    return;
    // preliminary...maybe we should use a line when N==2 and a quadratic when N==3
    // ...what about N==1 - maybe just fill the output array with a constant?
  }


  // Equations from Taschenbuch der Mathematik, (Bronstein u.a., 5.Auflage, Verlag Harri Deutsch), 
  // pages 955/956

  // polynomial coeffs per segment, there are N-1 segments but for a and c, we need one dummy 
  // coeff for the formulas:
  std::vector<Ty> a(N), b(N-1), c(N), d(N-1); 

  // establish array of a-coefficients - a[i] is actually the same as y[i] ...get rid
  for(int i = 0; i < N; i++) 
    a[i] = yIn[i];

  // establish the differences between the x-values:
  std::vector<Tx> h(N-1);
  for(size_t i = 0; i < h.size(); i++) 
    h[i] = xIn[i+1] - xIn[i];

  c[0] = c[N-1] = 0; // from the "natural" boundary conditions f''(0) = f''(N-1) = 0

  // establish the (N-2)x(N-2) tridiagonal matrix for computing the remaining c-coefficients 
  // (Eq. 19.240):
  std::vector<Ty> md(N-2), uld(N-3);    // arrays for main diagonal and upper/lower diagonal
  for(size_t i = 0; i < md.size(); i++)
    md[i] = 2*(h[i]+h[i+1]);
  for(size_t i = 0; i < uld.size(); i++)
    uld[i] = h[i+1];

  // establish the right-hand side for the tridiagonal system (Eq. 19.239):
  std::vector<Ty> rhs(N-2);
  Ty scl = scaleRhs * Ty(3);
  for(size_t i = 0; i < rhs.size(); i++)
    rhs[i] = scl * ((a[i+2]-a[i+1])/h[i+1] - (a[i+1]-a[i])/h[i]);
  // maybe introduce a scale factor - i think, when the rhs is 0, we get a linear interpolant

  // solve the tridiagonal system:
  rsLinearAlgebra::rsSolveTridiagonalSystem(&uld[0], &md[0], &uld[0], &rhs[0], &c[1], N-2);

  // compute b-coefficients by Eq. 19.238:
  for(size_t i = 0; i < b.size(); i++)
    b[i] = (a[i+1]-a[i])/h[i] - (2*c[i]+c[i+1])*h[i]/Ty(3);

  // compute d-coefficients by Eq. 19.237:
  for(size_t i = 0; i < d.size(); i++)
    d[i] = (c[i+1]-c[i])/(Ty(3)*h[i]);


  // OK, we have our polynomial coeffs for the segments - now do the actual interpolation:
  int i = 0;  // index into input data arrays (and polynomial coeffs)
  //int j = 0;  // index into output data arrays
  for(int j = 0; j < Ni; j++) {
    Tx x = xOut[j];
    while(xIn[i+1] < x && i < N-2)
      i++;
    yOut[j] = rsPolynomial<Ty>::evaluateCubic(x-xIn[i], a[i], b[i], c[i], d[i]); // Eq. 19.235
  }

  // todo: maybe refactor to separate the coefficient computation from the actual interpolation
  // the polynomial coefficients themselves could be very useful as basis for numerical 
  // integration - just add up the definite integrals over the spline segments between the
  // integration limits - maybe numeric differentiation could also be improved by using the 
  // spline - we could just evaluate the spline's derivative at the datapoints
}

template<class T>
void rsCubicSplineCoeffsFourPoints(T *a, T *y)
{
  T s[2];
  s[0] = 0.5*(y[1]-y[-1]); // left slope 
  s[1] = 0.5*(y[2]-y[0]);  // right slope
  rsPolynomial<T>::rsCubicCoeffsTwoPointsAndDerivatives(a, y, s);
}

template<class T>
void fitCubicWithDerivative(T x1, T x2, T y1, T y2, T yd1,
  T yd2, T *a3, T *a2, T *a1, T *a0)
{
  *a3 = -( x2*(yd2+yd1) + x1*(-yd2-yd1) - 2*y2 + 2*y1 );
  *a2 = x1*(x2*(yd2-yd1)-3*y2) + (x2*x2)*(yd2+2*yd1) + (x1*x1)*(-2*yd2-yd1)
       - 3*x2*y2 + (3*x2+3*x1)*y1;
  *a1 = x1*((x2*x2)*(2*yd2+yd1)-6*x2*y2) - (x1*x1*x1)*yd2 + (x1*x1)*x2*(-yd2-2*yd1)
       + (x2*x2*x2)*yd1 + 6*x1*x2*y1;
  *a1 = -*a1;
  *a0 = (x1*x1)*((x2*x2)*(yd2-yd1)-3*x2*y2) + (x1*x1*x1)*(y2-x2*yd2) + x1*(x2*x2*x2)*yd1
       + (3*x1*(x2*x2)-(x2*x2*x2))*y1;

  T scaler = 1.0 / ( -(x2*x2*x2) + 3*x1*(x2*x2) - 3*x1*x1*x2 + (x1*x1*x1) );
  *a3 *= scaler;
  *a2 *= scaler;
  *a1 *= scaler;
  *a0 *= scaler;

  // maybe precompute x1*x1*x1, x2*x2*x2 for optimization
}

template<class T>
void fitCubicWithDerivativeFixedX(T y0, T y1, T yd0, T yd1, T *a3, 
  T *a2, T *a1, T *a0)
{
  *a2 = -yd1-2*yd0+3*y1-3*y0;
  *a3 = yd1+yd0-2*y1+2*y0;
  *a1 = yd0;
  *a0 = y0;

  // maybe factor out the 3 in the 1st line and the 2 in the second for optimization - but make 
  // performance tests to check
}

template<class T>
void fitQuinticWithDerivativesFixedX(T y0, T y1, T yd0, T yd1, T ydd0, 
  T ydd1, T *a5, T *a4, T *a3, T *a2, T *a1, T *a0)
{
  *a0 = y0;
  *a1 = yd0;
  *a2 = ydd0/2;
  *a3 = (ydd1-8*yd1+20*y1-6*(*a2)-12*(*a1)-20*(*a0))/2;
  *a4 = -ydd1+7*yd1-15*y1+3*(*a2)+8*(*a1)+15*(*a0);
  *a5 = (ydd1-6*yd1+12*y1-2*(*a2)-6*(*a1)-12*(*a0))/2;
}

template<class T>
void getHermiteCoeffsM(const T *y0, const T *y1, T *a, int M)
{
  int n, i, j;

  // compute a[0],...,a[M]:
  int fn = 1;
  for(n = 0; n <= M; n++)
  {
    a[n] = y0[n] / fn;   // a[n] = y0[n] / factorial(n); 
    fn *= (n+1);
  }

  // establish right hand side (vector k):
  T *k = new T[M+1];
  for(n = 0; n <= M; n++)
  {
    k[n] = y1[n];
    for(i = n; i <= M; i++)
      k[n] -= rsProduct(i-n+1, i) * a[i];
  }

  // establish matrix A:
  rsMatrixDbl A(M+1, M+1);
  for(i = 1; i <= M+1; i++)
  {
    for(j = 1; j <= M+1; j++)
      A.set(i-1, j-1, rsProduct(M+j-i+2, M+j));
  }

  // solve the linear system and cleanup:
  rsLinearAlgebra::rsSolveLinearSystem(A.getDataPointer(), &a[M+1], k, M+1);
  delete[] k;
}

template<class T>
void getHermiteCoeffs1(const T *y0, const T *y1, T *a)
{
  a[0] = y0[0];     // a0 = y(0)
  a[1] = y0[1];     // a1 = y'(0)

  T k1 = y1[0] - a[1] - a[0];
  T k2 = y1[1] - a[1];

  a[2] = 3*k1 - k2;
  a[3] = k2 - 2*k1;
}
// ToDo: make a variant that takes x0, x1 as additional parameters - it should call this function
// internally and then convert the coeffs by considering the nested polynomial:
// p(x) = a0 + a1*(a*x+b) + a2*(a*x+b)^2 + a3*(a*x+b)^3 where a,b can be computed from x0,x1

template<class T>
void getHermiteCoeffs2(const T *y0, const T *y1, T *a)
{
  a[0] = y0[0];     // a0 = y(0)
  a[1] = y0[1];     // a1 = y'(0)
  a[2] = y0[2]/2;   // a2 = y''(0)/2

  T k1 = y1[0] - a[2]  - a[1] - a[0];
  T k2 = y1[1] - y0[2] - a[1];
  T k3 = y1[2] - y0[2];

  a[3] =  (k3-8*k2+20*k1)/2;
  a[4] =  -k3+7*k2-15*k1;
  a[5] =  (k3-6*k2+12*k1)/2;
}

template<class T>
void getHermiteCoeffs3(const T *y0, const T *y1, T *a)
{
  a[0] = y0[0];     // a0 = y(0)
  a[1] = y0[1];     // a1 = y'(0)
  a[2] = y0[2]/2;   // a2 = y''(0)/2
  a[3] = y0[3]/6;   // a3 = y'''(0)/6

  // can be streamlined, using: 2*a2 = y''(0), 6*a3 = y'''(0)
  T k1 = y1[0] -   a[3] -   a[2] - a[1] - a[0];
  T k2 = y1[1] - 3*a[3] - 2*a[2] - a[1];
  T k3 = y1[2] - 6*a[3] - 2*a[2];
  T k4 = y1[3] - 6*a[3];

  // maybe reorder summands so as to get rid of the unary minus(es):
  a[4] =  (-k4+15*k3-90*k2+210*k1)/6;
  a[5] = -(-k4+14*k3-78*k2+168*k1)/2;
  a[6] =  (-k4+13*k3-68*k2+140*k1)/2;
  a[7] = -(-k4+12*k3-60*k2+120*k1)/6;

  // optimize: replace divisions by multiplications
}

template<class T>
T getDelayedSampleLinear(T d, T *y)
{
  return (1.0-d)*y[0] + d*y[-1];
}

template<class T>
T getDelayedSampleAsymmetricHermite1(T d, T *y, T shape)
{
  T y0[2], y1[2], a[4];

  // desired signal values at endpoints:
  y0[0] = y[-1];
  y1[0] = y[0];

  // desired derivatives at endpoints:
  y0[1] = shape * (y[-1] - y[-2]);
  y1[1] = shape * (y[0]  - y[-1]);

  // compute polynomial coeffs:
  getHermiteCoeffs1(y0, y1, a);

  // evaluate:
  T x = 1.0 - d;
  return a[0] + a[1]*x + a[2]*x*x + a[3]*x*x*x;  // optimize this
}

template<class T>
T getDelayedSampleAsymmetricHermiteM(T d, T *y, int M, T shape)
{
  int N = 2*M+1;
  int i, j;

  /*
  T *y0   = (T*) alloca((M+1)*sizeof(T));
  T *y1   = (T*) alloca((M+1)*sizeof(T));
  T *yTmp = (T*) alloca((M+2)*sizeof(T));
  T *a    = (T*) alloca((N+1)*sizeof(T));
  */
  T *y0   = new T[M+1];
  T *y1   = new T[M+1];
  T *yTmp = new T[M+2];
  T *a    = new T[N+1];

  // create desired signal values and derivatives at the endpoints:
  y0[0] = y[-1];
  y1[0] = y[0];
  for(i = 0; i < M+2; i++)
    yTmp[i] = y[-i];
  for(j = 1; j <= M; j++)
  {
    for(i = 0; i < M+1; i++)
      yTmp[i] = yTmp[i] - yTmp[i+1];
    y1[j] = shape * yTmp[0];
    y0[j] = shape * yTmp[1];
    shape *= shape;
  }

  // compute polynomial coeffs:
  getHermiteCoeffsM(y0, y1, a, M);

  // evaluate:
  T x = 1.0 - d;
  T result = rsPolynomial<T>::evaluate(x, a, N);

  // cleanup and return result:
  delete[] y0;
  delete[] y1;
  delete[] yTmp;
  delete[] a;
  return result;
}

template<class T>
void fitCubicThroughFourPoints(T x0, T y0, T x1, T y1, T x2,
                               T y2, T x3, T y3, T *a, T *b,
                               T *c, T *d)
{
  // powers of the input value x:
  T x02 = x0*x0;    // x0^2
  T x12 = x1*x1;    // x1^2
  T x22 = x2*x2;    // x2^2
  T x32 = x3*x3;    // x3^2
  T x03 = x0*x02;   // x0^3
  T x13 = x1*x12;   // x1^3
  T x23 = x2*x22;   // x2^3
  T x33 = x3*x32;   // x3^3

  // some common subexpressions:
  T k1  = (x13*(y3-y2)-x23*y3+x33*y2+(x23-x33)*y1);
  T k2  = (x23*y3-x33*y2);
  T k3  = (x2*x33-x23*x3);
  T k4  = (x32*y2-x22*y3);
  T k5  = (x12*(x33-x23)-x22*x33+x23*x32+x13*(x22-x32));
  T k6  = (x2*y3+x1*(y2-y3)-x3*y2+(x3-x2)*y1);
  T k7  = (x22*x33-x23*x32);
  T k8  = (x2*x32-x22*x3);
  T k9  = (x3*y2-x2*y3);
  T k10 = (x22*y3-x32*y2);
  T k11 = (x23*x3-x2*x33);

  // a scaler that applies to all coefficients:
  T scaler = 1.0 / (x0*k5+x1*k7+x02*(x2*x33+x1*(x23-x33)+x13*(x3-x2)-x23*x3)+x12*k11+
    x03*(x1*(x32-x22)-x2*x32+x22*x3+x12*(x2-x3))+x13*k8);

  // the coefficients themselves:
  *a = (x0*(x12*(y3-y2)-x22*y3+x32*y2+(x22-x32)*y1)+x1*k10+x02*k6+x12*k9+k8*y1+
        (x1*(x32-x22)-x2*x32+x22*x3+x12*(x2-x3))*y0) * scaler;

  *b = -(x0*k1+x1*k2+x03*k6+x13*k9+k3*y1+(x1*(x33-x23)-x2*x33+x23*x3+x13*(x2-x3))*y0) * scaler;

  *c = (x02*k1+x12*k2+x03*(x22*y3+x12*(y2-y3)-x32*y2+(x32-x22)*y1)+x13*k4+k7*y1+k5*y0) * scaler;

  *d = -(x0*(x12*k2+x13*k4+k7*y1)+x02*(x1*(x33*y2-x23*y3)+x13*(x2*y3-x3*y2)+k11*y1)+x03*
         (x1*k10+x12*k9+k8*y1)+(x1*(x23*x32-x22*x33)+x12*k3+x13*(x22*x3-x2*x32))*y0) * scaler;

  // ...this seems all very complicated - can this be simplified? ...see formula for Lagrange 
  // polynomial in the handout "MUS421/EE367B Lecture 4" by Julius Smith ...this product formula 
  // with h_delta on page 13 -> generalize this function to compute coefficients for n-th order 
  // Lagrange interpolator
}

template<class T> 
T rsInterpolateWrapped(T x0, T x1, T t, T xMin, T xMax)
{
  RAPT::rsAssert(xMin <= x0 && x0 <= xMax);
  RAPT::rsAssert(xMin <= x1 && x1 <= xMax);

  T r  = xMax - xMin;  // range
  T du = x1 + r - x0;  // upper difference
  T dm = x1     - x0;  // middle difference
  T dl = x1 - r - x0;  // lower difference
  T au = rsAbs(du);
  T am = rsAbs(dm);
  T al = rsAbs(dl);
  T x;

  // add an appropriate fraction of the delta that has the smallest absolute difference to x0:
  if(au < am && au < al)
    x = x0 + t*du;
  else if(am < al)
    x = x0 + t*dm;
  else
    x = x0 + t*dl;

  // re-wrap result into allowed interval:
  if(x > xMax)
    x -= r;
  else if(x < xMin)
    x += r;
  return x;
}
// check if it is correct (maybe by a unit test?)...and if it can be simplified

template<class T>
void cubicSplineArcCoeffs2D(T x1, T dx1, T y1, T dy1, T x2, T dx2, T y2, T dy2, T* a, T* b)
{
  // Compute coeffs of the two polynomials:
  // x(t) = a0 + a1*t + a2*t^2 + a3*t^3
  // y(t) = b0 + b1*t + b2*t^2 + b3*t^3
  T z0[2], z1[2]; // y0, y1 inputs in getHermiteCoeffs1
  z0[0] = x1; z0[1] = dx1; z1[0] = x2; z1[1] = dx2; getHermiteCoeffs1(z0, z1, a);
  z0[0] = y1; z0[1] = dy1; z1[0] = y2; z1[1] = dy2; getHermiteCoeffs1(z0, z1, b);
}


template<class T>
void quadraticLineCoeffs2D(T x0, T y0, T x1, T y1, T* a, T* b)
{
  a[0] = x0;
  b[0] = y0;
  a[1] = x1-x0;
  b[1] = y1-y0;
  a[2] = 0;
  b[2] = 0;
}
// rename to lineCoeffs2D, don't access a[2], b[2], add to header

template<class T>
bool quadraticSplineArcCoeffs2D(T x0, T dx0, T y0, T dy0, T x1, T dx1, T y1, T dy1, T* a, T* b)
{
  a[0] = x0;
  b[0] = y0;
  T den, k; 
  T tol = T(1.e-6); // use (some multiple of) epsilon of T

  // catch 0/0 cases:
  if(  (abs(dx0) < tol && abs(dy0) < tol) 
    || (abs(dx1) < tol && abs(dy1) < tol) )
  {
    return false;
  }

  if(abs(dx0) > abs(dy0)) {
    T s0 = dy0/dx0;
    if(abs(dx1) > abs(dy1)) {  // compute a,b coeffs from s0, s1
      T s1 = dy1/dx1;
      den  = s0 - s1;
      if(abs(den) < tol)  
        return false;
      k    = T(1)/den;
      a[1] = 2*(s1*x0 - s1*x1 - y0 + y1)*k;
      a[2] = -((s0 + s1)*x0 - (s0 + s1)*x1 - 2*y0 + 2*y1)*k;
      b[1] = 2*(s0*s1*x0 - s0*s1*x1 - s0*y0 + s0*y1)*k;
      b[2] = -(2*s0*s1*x0 - 2*s0*s1*x1 - (s0 + s1)*y0 + (s0 + s1)*y1)*k;
    }
    else {                     // compute a,b coeffs from s0, r1
      T r1 = dx1/dy1;
      den  = (r1*s0 - 1);
      if(abs(den) < tol) 
        return false;
      k    = T(1)/den;
      a[1] = -2*(r1*y0 - r1*y1 - x0 + x1)*k;
      a[2] = -((r1*s0 + 1)*x0 - (r1*s0 + 1)*x1 - 2*r1*y0 + 2*r1*y1)*k;
      b[1] = -2*(r1*s0*y0 - r1*s0*y1 - s0*x0 + s0*x1)*k;
      b[2] = -(2*s0*x0 - 2*s0*x1 - (r1*s0 + 1)*y0 + (r1*s0 + 1)*y1)*k;
    }
  }
  else {
    T r0 = dx0/dy0;
    if(abs(dx1) > abs(dy1)) {  // compute a,b coeffs from r0, s1
      T s1 = dy1/dx1;
      den  = (r0*s1 - 1);
      if(abs(den) < tol) 
        return false;
      k    = T(1)/den;
      a[1] = -2*(r0*s1*x0 - r0*s1*x1 - r0*y0 + r0*y1)*k;
      a[2] = ((r0*s1 + 1)*x0 - (r0*s1 + 1)*x1 - 2*r0*y0 + 2*r0*y1)*k;
      b[1] = -2*(s1*x0 - s1*x1 - y0 + y1)*k;
      b[2] = (2*s1*x0 - 2*s1*x1 - (r0*s1 + 1)*y0 + (r0*s1 + 1)*y1)*k;
    }
    else {                     // compute a,b coeffs from r0, r1
      T r1 = dx1/dy1;
      den  = (r0 - r1);
      if(abs(den) < tol) 
        return false;
      k    = T(1)/den;
      a[1] = 2*(r0*r1*y0 - r0*r1*y1 - r0*x0 + r0*x1)*k;
      a[2] = -(2*r0*r1*y0 - 2*r0*r1*y1 - (r0 + r1)*x0 + (r0 + r1)*x1)*k;
      b[1] = 2*(r1*y0 - r1*y1 - x0 + x1)*k;
      b[2] = -((r0 + r1)*y0 - (r0 + r1)*y1 - 2*x0 + 2*x1)*k;
    }
  }
  return true;
}
/*
// optimized computation for s0,s1 cae in function above:
TCor dx, dy, ss, k, s1dx;
dx   = x1-x0;
dy   = y1-y0;
k    = 1/(s0-s1);
ss   = s0+s1;     // slope sum
s1dx = s1*dx;     // rename to k1
a[1] = 2*(    dy-s1dx) *k;
b[1] = 2*(s0*(dy-s1dx))*k;  // use k2 = (dy-s1dx)
a[2] = (ss*dx - 2*dy)*k;
b[2] = (2*s0*s1dx - ss*dy)*k;
// make optimized versions for all cases...
*/

template<class T>
void quadraticOrCubicSplineArcCoeffs2D(T x0, T dx0, T y0, T dy0, T x1, T dx1, T y1, T dy1, 
  T* a, T* b)
{
  bool success = quadraticSplineArcCoeffs2D(x0, dx0, y0, dy0, x1, dx1, y1, dy1, a, b);
  if(!success)
    cubicSplineArcCoeffs2D(x0, dx0, y0, dy0, x1, dx1, y1, dy1, a, b);
    //quadraticLineCoeffs2D(x0, y0, x1, y1, a, b); // preliminary - use cubic
  else {
    a[3] = 0;
    b[3] = 0;
  }
}
// it still sometimes seem to produce linear segments (try (2,3) lissjous figure with 11 samples)


template<class T>
void cubicSplineArcLength2D(T *a, T *b, T *t, T* s, int N)
{
  // The arc-length s(t) between 0 and t of the 2D cubic spline defined by the two polynomials:
  //   x(t) = a0 + a1*t + a2*t^2 + a3*t^3
  //   y(t) = b0 + b1*t + b2*t^2 + b3*t^3
  // is given by the definite integral from 0 to t over the integrand:
  //   c(t) = sqrt( (dx/dt)^2 + (dy/dt)^2 )
  // where the term inside the square-root is a fourth degree polynomial (the derivative of a cubic
  // is a quadratic, squaring that gives a quartic and adding two quartics gives still a quartic). 
  // We evaluate the integrand at the N values t[n] and perform a numeric integration over these 
  // integrand values.

  // Find coeffs for quartic polynomial under the square-root in the integrand:
  typedef rsPolynomial<T> PL;
  T c[5], d[5];                       // coeffs of:
  PL::derivative(a, c, 3);            // c is dx/dt (a is x(t))
  PL::derivative(b, d, 3);            // d is dy/dt (b is y(t))
  PL::multiply(c, 2, c, 2, c);        // c is (dx/dt)^2
  PL::multiply(d, 2, d, 2, d);        // d is (dy/dt)^2
  rsArrayTools::add(c, d, c, 5);      // c is (dx/dt)^2 + (dy/dt)^2
  // The coeffs of our desired quartic are now in our c-array.

  // Evaluate the integrand at the given t-values and perform numeric integration:
  for(int n = 0; n < N; n++)
    s[n] = sqrt(rsMax(T(0), PL::evaluate(t[n], c, 4))); // write and use optimized evaluateQuartic
  rsNumericIntegral(t, s, s, N); // integration works in place (we use s for integrand and integral)
}
// todo: make a 3D version (the only difference is that we have 3 polynomials x(t),y(t),z(t) that 
// we have to take derivates of, square and add...or maybe make an N-dimensional version - just 
// compute one (squared) derivative per dimension and accumulate the resulting quartics - the 
// result will always be just a 1D quartic, regardless of the number of dimensions of the space
// maybe make a function that can also handle higher order splines...at least quartics because
// we want to try a quartic interpolant - maybe instead of having a,b arrays of scalar 
// coefficients, the coefficients could be vectors (2D, 3D, nD)...hmmm...but at the end, we need 
// the c-array to be an array of scalars - i think the scalar coeff is just the sum of the 
// "vector-coeff" elements...figure out, if that can be made to work...
