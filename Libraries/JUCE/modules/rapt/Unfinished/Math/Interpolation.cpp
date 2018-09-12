template<class T>
void resampleLinear(T* xIn, T* yIn, int inLength, T* xOut, T* yOut, int outLength)
{
  // code copied and adapted from from rsTimeWarper::invertMonotonousWarpMap
  int iOld = 0;
  for(int n = 0; n < outLength; n++) {
    T x = xOut[n];
    for(int i = iOld; i < inLength-1; i++) {
      if( xIn[i] <= x && xIn[i+1] >= x )      // find i, such that xIn[i] <= x and xIn[i+1] >= x
        break;
    }
    yOut[n] = rsInterpolateLinear(xIn[i], xIn[i+1], yIn[i], yIn[i+1], x);
    iOld = i;
  }
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
void getHermiteCoeffsM(T *y0, T *y1, T *a, int M)
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
void getHermiteCoeffs1(T *y0, T *y1, T *a)
{
  a[0] = y0[0];     // a0 = y(0)
  a[1] = y0[1];     // a1 = y'(0)

  T k1 = y1[0] - a[1] - a[0];
  T k2 = y1[1] - a[1];

  a[2] = 3*k1 - k2;
  a[3] = k2 - 2*k1;
}

template<class T>
void getHermiteCoeffs2(T *y0, T *y1, T *a)
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
void getHermiteCoeffs3(T *y0, T *y1, T *a)
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
  T result = rsPolynomial<T>::evaluatePolynomialAt(x, a, N);

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
