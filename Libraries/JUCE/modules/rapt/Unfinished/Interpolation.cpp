using namespace RSLib;

void RSLib::rsCubicSplineCoeffsFourPoints(double *a, double *y)
{
  double s[2];
  s[0] = 0.5*(y[1]-y[-1]); // left slope 
  s[1] = 0.5*(y[2]-y[0]);  // right slope
  rsCubicCoeffsTwoPointsAndDerivatives(a, y, s);
}

void RSLib::fitCubicWithDerivative(double x1, double x2, double y1, double y2, double yd1,
  double yd2, double *a3, double *a2, double *a1, double *a0)
{
  *a3 = -( x2*(yd2+yd1) + x1*(-yd2-yd1) - 2*y2 + 2*y1 );
  *a2 = x1*(x2*(yd2-yd1)-3*y2) + (x2*x2)*(yd2+2*yd1) + (x1*x1)*(-2*yd2-yd1)
       - 3*x2*y2 + (3*x2+3*x1)*y1;
  *a1 = x1*((x2*x2)*(2*yd2+yd1)-6*x2*y2) - (x1*x1*x1)*yd2 + (x1*x1)*x2*(-yd2-2*yd1)
       + (x2*x2*x2)*yd1 + 6*x1*x2*y1;
  *a1 = -*a1;
  *a0 = (x1*x1)*((x2*x2)*(yd2-yd1)-3*x2*y2) + (x1*x1*x1)*(y2-x2*yd2) + x1*(x2*x2*x2)*yd1
       + (3*x1*(x2*x2)-(x2*x2*x2))*y1;

  double scaler = 1.0 / ( -(x2*x2*x2) + 3*x1*(x2*x2) - 3*x1*x1*x2 + (x1*x1*x1) );
  *a3 *= scaler;
  *a2 *= scaler;
  *a1 *= scaler;
  *a0 *= scaler;

  // maybe precompute x1*x1*x1, x2*x2*x2 for optimization
}

void RSLib::fitCubicWithDerivativeFixedX(double y0, double y1, double yd0, double yd1, double *a3, 
  double *a2, double *a1, double *a0)
{
  *a2 = -yd1-2*yd0+3*y1-3*y0;
  *a3 = yd1+yd0-2*y1+2*y0;
  *a1 = yd0;
  *a0 = y0;

  // maybe factor out the 3 in the 1st line and the 2 in the second for optimization - but make 
  // performance tests to check
}

void RSLib::fitQuinticWithDerivativesFixedX(double y0, double y1, double yd0, double yd1, double ydd0, double ydd1, double *a5, double *a4,
                                            double *a3, double *a2, double *a1, double *a0)
{
  *a0 = y0;
  *a1 = yd0;
  *a2 = ydd0/2;
  *a3 = (ydd1-8*yd1+20*y1-6*(*a2)-12*(*a1)-20*(*a0))/2;
  *a4 = -ydd1+7*yd1-15*y1+3*(*a2)+8*(*a1)+15*(*a0);
  *a5 = (ydd1-6*yd1+12*y1-2*(*a2)-6*(*a1)-12*(*a0))/2;
}

void RSLib::getHermiteCoeffsM(double *y0, double *y1, double *a, int M)
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
  double *k = new double[M+1];
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
  RSLib::rsSolveLinearSystem(A.getDataPointer(), &a[M+1], k, M+1);
  delete[] k;
}

void RSLib::getHermiteCoeffs1(double *y0, double *y1, double *a)
{
  a[0] = y0[0];     // a0 = y(0)
  a[1] = y0[1];     // a1 = y'(0)

  double k1 = y1[0] - a[1] - a[0];
  double k2 = y1[1] - a[1];

  a[2] = 3*k1 - k2;
  a[3] = k2 - 2*k1;
}

void RSLib::getHermiteCoeffs2(double *y0, double *y1, double *a)
{
  a[0] = y0[0];     // a0 = y(0)
  a[1] = y0[1];     // a1 = y'(0)
  a[2] = y0[2]/2;   // a2 = y''(0)/2

  double k1 = y1[0] - a[2]  - a[1] - a[0];
  double k2 = y1[1] - y0[2] - a[1];
  double k3 = y1[2] - y0[2];

  a[3] =  (k3-8*k2+20*k1)/2;
  a[4] =  -k3+7*k2-15*k1;
  a[5] =  (k3-6*k2+12*k1)/2;
}

void RSLib::getHermiteCoeffs3(double *y0, double *y1, double *a)
{
  a[0] = y0[0];     // a0 = y(0)
  a[1] = y0[1];     // a1 = y'(0)
  a[2] = y0[2]/2;   // a2 = y''(0)/2
  a[3] = y0[3]/6;   // a3 = y'''(0)/6

  // can be streamlined, using: 2*a2 = y''(0), 6*a3 = y'''(0)
  double k1 = y1[0] -   a[3] -   a[2] - a[1] - a[0];
  double k2 = y1[1] - 3*a[3] - 2*a[2] - a[1];
  double k3 = y1[2] - 6*a[3] - 2*a[2];
  double k4 = y1[3] - 6*a[3];

  // maybe reorder summands so as to get rid of the unary minus(es):
  a[4] =  (-k4+15*k3-90*k2+210*k1)/6;
  a[5] = -(-k4+14*k3-78*k2+168*k1)/2;
  a[6] =  (-k4+13*k3-68*k2+140*k1)/2;
  a[7] = -(-k4+12*k3-60*k2+120*k1)/6;
}

double RSLib::getDelayedSampleLinear(double d, double *y)
{
  return (1.0-d)*y[0] + d*y[-1];
}

double RSLib::getDelayedSampleAsymmetricHermite1(double d, double *y, double shape)
{
  double y0[2], y1[2], a[4];

  // desired signal values at endpoints:
  y0[0] = y[-1];
  y1[0] = y[0];

  // desired derivatives at endpoints:
  y0[1] = shape * (y[-1] - y[-2]);
  y1[1] = shape * (y[0]  - y[-1]);

  // compute polynomial coeffs:
  getHermiteCoeffs1(y0, y1, a);

  // evaluate:
  double x = 1.0 - d;
  return a[0] + a[1]*x + a[2]*x*x + a[3]*x*x*x;  // optimize this
}

double RSLib::getDelayedSampleAsymmetricHermiteM(double d, double *y, int M, double shape)
{
  int N = 2*M+1;
  int i, j;

  /*
  double *y0   = (double*) alloca((M+1)*sizeof(double));
  double *y1   = (double*) alloca((M+1)*sizeof(double));
  double *yTmp = (double*) alloca((M+2)*sizeof(double));
  double *a    = (double*) alloca((N+1)*sizeof(double));
  */
  double *y0   = new double[M+1];
  double *y1   = new double[M+1];
  double *yTmp = new double[M+2];
  double *a    = new double[N+1];

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
  double x = 1.0 - d;
  double result = RSLib::evaluatePolynomialAt(x, a, N);

  // cleanup and return result:
  delete[] y0;
  delete[] y1;
  delete[] yTmp;
  delete[] a;
  return result;
}

void RSLib::fitCubicThroughFourPoints(double x0, double y0, double x1, double y1, double x2,
                                      double y2, double x3, double y3, double *a, double *b,
                                      double *c, double *d)
{
  // powers of the input value x:
  double x02 = x0*x0;    // x0^2
  double x12 = x1*x1;    // x1^2
  double x22 = x2*x2;    // x2^2
  double x32 = x3*x3;    // x3^2
  double x03 = x0*x02;   // x0^3
  double x13 = x1*x12;   // x1^3
  double x23 = x2*x22;   // x2^3
  double x33 = x3*x32;   // x3^3

  // some common subexpressions:
  double k1  = (x13*(y3-y2)-x23*y3+x33*y2+(x23-x33)*y1);
  double k2  = (x23*y3-x33*y2);
  double k3  = (x2*x33-x23*x3);
  double k4  = (x32*y2-x22*y3);
  double k5  = (x12*(x33-x23)-x22*x33+x23*x32+x13*(x22-x32));
  double k6  = (x2*y3+x1*(y2-y3)-x3*y2+(x3-x2)*y1);
  double k7  = (x22*x33-x23*x32);
  double k8  = (x2*x32-x22*x3);
  double k9  = (x3*y2-x2*y3);
  double k10 = (x22*y3-x32*y2);
  double k11 = (x23*x3-x2*x33);

  // a scaler that applies to all coefficients:
  double scaler = 1.0 / (x0*k5+x1*k7+x02*(x2*x33+x1*(x23-x33)+x13*(x3-x2)-x23*x3)+x12*k11+
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
