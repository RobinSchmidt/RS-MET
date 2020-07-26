#include "Prototypes.h"

#include "FilterDesign/PoleZeroPrototype.cpp"
#include "FilterDesign/PoleZeroMapper.cpp"
#include "FilterDesign/PoleZeroDesignerAnalog.cpp"
#include "FilterDesign/PoleZeroDesignerDigital.cpp"
#include "FilterDesign/ComplementaryFilters.cpp"

#include "ModalAnalyzer.cpp"
#include "Probability.cpp"
#include "Projection3Dto2D.cpp"
#include "Polygon.cpp"
#include "Drawing.cpp"
#include "QuantumSystems.cpp"
#include "Relativity.cpp"
#include "SineParameterEstimator.cpp"

//#include "SinusoidalModeling.cpp" // moved to rapt

using namespace RAPT;

std::vector<double> solvePentaDiagonalSystem(
  std::vector<double>& M, std::vector<double>& L,
  std::vector<double>& D,
  std::vector<double>& U, std::vector<double>& V,
  std::vector<double>& B)
{
  int N = (int) D.size();
  rsAssert(B.size() == N);
  rsAssert(M.size() == N-2);
  rsAssert(L.size() == N-1);
  rsAssert(V.size() == N-2);
  rsAssert(U.size() == N-1);

  // Gaussian elimination without pivot-search - we just always use D[i] as pivot element:
  int i;
  double k;
  std::vector<double> piv(N-1);  // pivot elements for inspection only
  for(i = 0; i < N-2; i++) {
    piv[i]  = D[i];              // pivot element
    k       = L[i]/D[i];
    //L[i]   -= k*D[i];          // should give 0
    D[i+1] -= k*U[i];
    B[i+1] -= k*B[i];    
    U[i+1] -= k*V[i];
    k       = M[i]/D[i];
    //M[i]   -= k*D[i];          // should give 0
    L[i+1] -= k*U[i];
    D[i+2] -= k*V[i];
    B[i+2] -= k*B[i];
    //int dummy = 0;
  }
  piv[i]  = D[i];
  k       = L[i]/D[i];     // a final partial step outside the loop
  //L[i]   -= k*D[i];        // should give 0
  D[i+1] -= k*U[i];
  B[i+1] -= k*B[i];

  // we may get rid of the steps where 0 should result we could directly assign the values or maybe
  // just leave the step out (i think, the values are not needed in subsequent steps)
  // maybe use += (define k with negative sign)
  // maybe split the function into 2: elimination and backsubstitution (but have a convenience 
  // function that calls both)
  // maybe return the minimum absolute pivot that was encountered

  // Gaussian elimination is done - now do the backsubstitution to find the solution vector:
  std::vector<double> x(N);
  x[N-1] =  B[N-1]                  / D[N-1];
  x[N-2] = (B[N-2] - U[N-2]*x[N-1]) / D[N-2];
  for(i = N-3; i >= 0; i--)
    x[i] = (B[i] - U[i]*x[i+1] - V[i]*x[i+2]) / D[i];

  return x;
}

std::vector<double> pentaDiagMatVecMul(
  std::vector<double>& m, std::vector<double>& l,
  std::vector<double>& d,
  std::vector<double>& u, std::vector<double>& v,
  std::vector<double>& x)
{
  int N = (int) d.size();
  int i;
  std::vector<double> c(N);
  c[0] = d[0]*x[0] + u[0]*x[1] + v[0]*x[2];
  c[1] = l[0]*x[0] + d[1]*x[1] + u[1]*x[2] + v[1]*x[3];
  for(i = 2; i < N-2; i++)
    c[i] = m[i-2]*x[i-2] + l[i-1]*x[i-1] + d[i]*x[i] + u[i]*x[i+1] + v[i]*x[i+2];
  c[i] = m[i-2]*x[i-2] + l[i-1]*x[i-1] + d[i]*x[i] + u[i]*x[i+1];
  i++;
  c[i] = m[i-2]*x[i-2] + l[i-1]*x[i-1] + d[i]*x[i];
  return c;
}

std::vector<double> rsMinSqrDifFixSum(const std::vector<double>& s, 
  const std::vector<double>& w)
{
  int Ns = (int) s.size();  // number of sums
  int Nv = Ns + 1;          // number of values
  int Nm = Nv + Ns;         // number of linear equations, matrix size
  typedef std::vector<double> Vec;

  // handle special case for just one given sum:
  if(Ns == 1) {
    Vec v(2);
    v[0] = v[1] = 0.5*s[0];
    return v;
  }

  // establish the diagonals for the coefficient matrix:
  Vec d0(Nm), d1(Nm-1), d2(Nm-2);
  int i;
  d0[0] = 2;
  for(i = 1; i < Nm; i += 2) {
    d0[i]   = 0;
    d0[i+1] = 4; }
  d0[Nm-1] = 2;
  for(i = 0; i < Nm-1; i++)
    d1[i] = 1;
  for(i = 0; i < Nm-3; i += 2) {
    d2[i]   = -2;
    d2[i+1] =  0; }
  d2[i] = -2;

  // apply error weights, if desired:
  if(!w.empty())
  {
    rsAssert((int)w.size() == Ns); // weight-vector, if given, must have same size as sum-vector

    // modify outer sub- and superdiagonal:
    for(i = 0; i < Ns; i++)
      d2[2*i] *= w[i];

    // modify main diagonal:
    d0[0]    *= w[0];
    d0[Nm-1] *= w[Ns-1];
    for(i = 1; i < Ns; i++)
      d0[2*i] *= 0.5*(w[i-1]+w[i]);
  }

  // establish right-hand side vector:
  Vec b(Nm);
  int j = 0;
  for(i = 0; i <= Nm-2; i += 2) {
    b[i]   = 0;
    b[i+1] = s[j];
    j++;
  }
  b[Nm-1] = 0;

  // use temporaries, because things get messed up in the solver:
  Vec bt = b, l2 = d2, l1 = d1, d = d0, u1 = d1, u2 = d2;
  Vec x = solvePentaDiagonalSystem(l2, l1, d, u1, u2, bt);

  // verification - b2 should be equal to b, if the solver is legit:
  Vec b2 = pentaDiagMatVecMul(d2, d1, d0, d1, d2, x);

  // extract output array v (in x, the outputs are interleaved with the Lagrange multipliers):
  Vec v(Nv);
  for(i = 0; i < Nv; i++)
    v[i] = x[2*i];
  return v;
}

std::vector<double> rsMinSqrCrvFixSum(const std::vector<double>& s, const std::vector<double>& /*w*/)
{
  int Ns = (int) s.size();  // number of sums
  int Nv = Ns + 1;          // number of values
  int Nm = Nv + Ns;         // number of linear equations, matrix size
  typedef std::vector<double> Vec;

  // todo: catch special cases when Nm is small


  // create solver object and populate the coefficient matrix:
  rsBandDiagonalSolver<double> solver(Nm, 4, 4);
  solver.initMatrixWithZeros();
  int i;

  // main diagonal:
  solver.setDiagonalElement(  0,    0,  2);
  solver.setDiagonalElement(  0,    2, 10);
  for(i = 4; i < Nm-3; i += 2)
    solver.setDiagonalElement(0, i,    12);
  solver.setDiagonalElement(  0, Nm-3, 10);
  solver.setDiagonalElement(  0, Nm-1,  2);

  // first upper and lower diagonals (all ones):
  for(i = 0; i < Nm-1; i++) {
    solver.setDiagonalElement(+1, i, 1);  // upper
    solver.setDiagonalElement(-1, i, 1);  // lower
  }

  // second upper and lower diagonals: -4,0,-8,0,-8,...,-8,0,-8,0,-4
  solver.setDiagonalElement(  +2,    0, -4);
  solver.setDiagonalElement(  -2,    0, -4);
  for(i = 2; i < Nm-3; i += 2) {
    solver.setDiagonalElement(+2,    i, -8);
    solver.setDiagonalElement(-2,    i, -8);
  }
  solver.setDiagonalElement(  +2, Nm-3, -4);
  solver.setDiagonalElement(  -2, Nm-3, -4);

  // third upper and lower diagonals are all zeros and 3rd upper and lower diagonals are filled
  // with 2,0,2,0,2,0,...,0,2:
  for(i = 0; i < Nm-4; i += 2) {
    solver.setDiagonalElement(+4, i, 2);
    solver.setDiagonalElement(-4, i, 2);
  }

  // establish right-hand side vector:
  Vec b(Nm);
  int j = 0;
  for(i = 0; i <= Nm-2; i += 2) {
    b[i]   = 0;
    b[i+1] = s[j];
    j++;
  }
  b[Nm-1] = 0;

  // solve the system of equations:
  Vec x(Nm);
  solver.solve(&x[0], &b[0], 1);

  // extract output array v (in x, the outputs are interleaved with the Lagrange multipliers):
  Vec v(Nv);
  for(i = 0; i < Nv; i++)
    v[i] = x[2*i];
  return v;
}


double signalValueViaSincAt(double *x, int N, double t, double sincLength, double stretch,
  //FunctionPointer3DoublesToDouble windowFunction, 
  double (*windowFunction)(double,double,double),
  double windowParameter, bool normalizeDC)
{
  int L  = (int) floor(sincLength);
  int ti = (int) floor(t);         // integer part of t
  if( ti < 0 || ti >= N )
    return 0.0;
  double tf = t - ti;              // fractional part of t    
  double y  = 0.0;                 // output value
  double b;                        // weight for a sample
  double bs = 0.0;                 // sum of the weights
  int mMin  = -rsMin(L/2-1, ti);   // minimum shift
  int mMax  = +rsMin(L/2, N-ti-1); // maximum shift
  int    nr;                       // sample read position
  double ts;                       // time instant with shift
  for(int m = mMin; m <= mMax; m++)
  {
    ts = m-tf;
    b  = rsNormalizedSinc(ts/stretch) * windowFunction(ts, sincLength, windowParameter);
    bs += b;
    nr = ti+m;
    y += b * x[nr]; 
  }
  if( normalizeDC == true )
    return y/bs;
  else
    return y/stretch;
}

void halpernU(double *a, int K)
{
  RAPT::rsArrayTools::fillWithZeros(a, K+1);
  int k, m;
  int f1, f2, f3, f4; // factorials
  if( rsIsEven(K) )
  {
    k = K/2;
    for(m = 0; m <= k; m++)
    {
      // a[2*m] = (-1)^(k-m) * (m+k)! / ( (k-m)! * (m!)^2 ):
      f1 = rsFactorial(m+k);
      f2 = rsFactorial(k-m);
      f3 = rsFactorial(m);
      a[2*m] = pow(-1.0, k-m) * f1 / (f2*f3*f3);
    }
  }
  else
  {
    k = (K-1)/2;
    for(m = 0; m <= k; m++)
    {
      // a[2*(k-m)+1] = (-1)^m * (2*k+1-m)! / (m! * (k+1-m)! * (k-m)!):
      f1 = rsFactorial(2*k+1-m);
      f2 = rsFactorial(m);
      f3 = rsFactorial(k+1-m);
      f4 = rsFactorial(k-m);
      a[2*(k-m)+1] = pow(-1.0, m) * f1 / (f2*f3*f4);
    }
  }

  // we leave the normalization factors out, in order to keep all numbers integers. the squared
  // factors (which are also integers) will be applied, after the coefficient array was convolved
  // with itself
}
void halpernT2(double *c, int N)
{
  double u[30];  
  double s[30];
  double t[30];

  // create the U-polynomial:
  halpernU(u, N-1);

  // square it:
  RAPT::rsArrayTools::convolve(u, N, u, N, s);

  // apply squared scale factors:
  int k;
  if( rsIsEven(N) )
  {
    k = (N-2)/2;
    RAPT::rsArrayTools::scale(s, 2*N-1, 4*(k+1));
  }
  else
  {
    k = (N-1)/2;
    RAPT::rsArrayTools::scale(s, 2*N-1, 4*k+2);
  }

  // multiply it by x:
  RAPT::rsArrayTools::rightShift(s, 2*N, 1);

  // integrate it from 0 to w:
  double a[1] = { 0 };    // lower integration limit (as polynomial)
  double b[2] = { 0, 1 }; // upper integration limit (as polynomial)
  rsPolynomial<double>::integrateWithPolynomialLimits(s, 2*N-1, a, 0, b, 1, t);
    // this is actually the same as just taking the antiderivative of s

  // copy to output:
  RAPT::rsArrayTools::copy(t, c, 2*N+1);
}

void papoulisL2(double *v, int N)
{
  double P[20];  // array for Legendre polynomials
  int k;
  if( rsIsOdd(N) )
  {
    k = (N-1)/2;
    RAPT::rsArrayTools::fillWithZeros(v, k+1);
    for(int r = 0; r <= k; r++)
    {
      // add the weighted Legendre polynomial of order r to our v polynomial:
      rsPolynomial<double>::legendrePolynomial(P, r);
      rsPolynomial<double>::weightedSum(v, r, 1.0, P, r, 2*r+1.0, v);
    }

    // square it:
    RAPT::rsArrayTools::convolve(v, k+1, v, k+1, v);

    // account for leaving out the division by (k+1)*sqrt(2) in the weights in the accumulation:
    RAPT::rsArrayTools::scale(v, 2*k+1, 1.0/(2*(k+1)*(k+1)));

    // integrate from -1 to 2*w^2-1
    double a[1] = { -1 };
    double b[3] = { -1, 0, 2};
    rsPolynomial<double>::integrateWithPolynomialLimits(v, 2*k, a, 0, b, 2, v);
  }
  else
  {
    k = (N-2)/2;

    // generate Legendre Polynomial of order k+1, store in P:
    rsPolynomial<double>::legendrePolynomial(P, k+1);

    // take the derivative, store in v:
    rsPolynomial<double>::derivative(P, v, k+1);  // v has order k, length k+1

    // square it:
    RAPT::rsArrayTools::convolve(v, k+1, v, k+1, v);   // has order 2*k, length 2*k+1

    // multiply by (x+1):
    double xp1[2] = { 1, 1};
    RAPT::rsArrayTools::convolve(v, 2*k+1, xp1, 2, v); // has order 2*k+1, length 2*k+2=N

    // integrate from -1 to 2*w^2-1
    double a[1] = { -1 };
    double b[3] = { -1, 0, 2};
    rsPolynomial<double>::integrateWithPolynomialLimits(v, 2*k+1, a, 0, b, 2, v); // order: 2*N, length 2*N+1

    // scale, such that L^2(1) = 1:
    double K = 1.0 / rsPolynomial<double>::evaluate(1.0, v, 2*N);
    RAPT::rsArrayTools::scale(v, 2*N+1, K); 
  }
}

void rsDampedSineFilterAnalysis(double b0, double b1, double a1, double a2, double* w, double* A,
  double* d, double* p)
{
  rsAssert(0.25*a1*a1-a2 < 0.0, "no damped sine filter, poles not complex conjugate");
  double P, cw;
  P  = sqrt(a2);
  cw = -0.5*a1/P;
  *d = -1.0/log(P);
  *w = acos(cw);
  *p = atan2(sin(*w), b1/(P*b0)+cw);
  if(rsAbs(b0) > rsAbs(b1))
    * A = b0/sin(*p);
  else
    *A = b1/(P*sin(*w-*p));
  if(*A < 0.0)
  {
    *A  = -*A;
    *p += PI;
  }
}

void rsDampedSineFilterAnalysis2(double b0, double b1, double a1, double a2, double* w, double* A,
  double* d, double* p)
{
  rsAssert(0.25*a1*a1-a2 < 0.0, "no damped sine filter, poles not complex conjugate");
  typedef std::complex<double> Complex;
  Complex j(0.0, 1.0);                 // imaginary unit
  double P = sqrt(a2);                 // pole radius
  *w = acos(-0.5*a1/P);                // pole angle
  Complex q = P * exp(j * *w);         // pole location
  Complex r = (b1+b0*q)/(2*q.imag());  // residue location
  *d = -1.0/log(P);                    // normalized decay time constant
  *A = 2*abs(r);                       // amplitude
  *p = arg(r);                         // start phase...
  if(*p < 0.0)                         // ...in interval 0...2pi instead of -pi...pi
    * p += 2*PI;
  // Remark: There are actually two mathematical errors in this sequence of assignments which 
  // conveniently cancel each other and streamline the implementation, that's why I left them in.
  // Actually, it should be r = (b1+b0*q)/(2.0*j*q.im) and *p = arg(r) + 0.5*PI, so we have 
  // first missed a division by j (corresponding to a rotation by -pi/2) in the computation of the
  // residue r and that's why we later don't need to add pi/2 to the startphase value ;-)
}

void rsDampedSineFilterResidueAndPole(double b0, double b1, double a1, double a2,
  std::complex<double>* r, std::complex<double>* p)
{
  rsAssert(0.25*a1*a1-a2 < 0.0, "no damped sine filter, poles not complex conjugate");
  double P = sqrt(a2);                        // pole radius
  double w = acos(-0.5*a1/P);                 // pole angle
  std::complex<double> j(0.0, 1.0);           // imaginary unit
  *p = P * exp(j * w);                        // pole location
  *r = (b1 + b0 * *p) / (2. * j * p->imag()); // residue location
}


double cheby_poly(int n, double x) // Chebyshev polyomial T_n(x), does NOT work for x < -1
{
  if(fabs(x) <= 1) return cos( n*acos( x));
  else             return cosh(n*acosh(x)); 
}
void cheby_win(double *out, int N, double atten)
{
  // Prototype implementation with O(N^2) scaling of the computational cost
  int nn, i;
  double M, n, sum = 0; // max=0;
  double tg = pow(10,atten/20);         // 1/r term [2], 10^gamma [2]
  double x0 = cosh((1.0/(N-1))*acosh(tg));
  M = (N-1)/2;
  if(N%2==0) M = M + 0.5;               // handle even length windows 
  for(nn=0; nn<(N/2+1); nn++){
    n = nn-M;
    sum = 0;
    for(i=1; i<=M; i++){
      sum += cheby_poly(N-1,x0*cos(PI*i/N))*cos(2.0*n*PI*i/N); }
    out[nn] = tg + 2*sum;
    out[N-nn-1] = out[nn]; }
  RAPT::rsArrayTools::normalizeMean(out, N);
}

/*

The code above was adapted from here:
http://practicalcryptography.com/miscellaneous/machine-learning/implementing-dolph-chebyshev-window/

For more info about doing an FFT based implementation, see:

https://www.dsprelated.com/freebooks/sasp/Dolph_Chebyshev_Window.html
https://ccrma.stanford.edu/~jos/sasp/Dolph_Chebyshev_Window.html
https://en.wikipedia.org/wiki/Window_function#Dolph%E2%80%93Chebyshev_window
https://octave.sourceforge.io/signal/function/chebwin.html

This paper is good: The Dolph–Chebyshev Window: A Simple Optimal Filter (Lynch)
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.77.5098&rep=rep1&type=pdf

https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.signal.chebwin.html
https://github.com/scipy/scipy/blob/v0.19.0/scipy/signal/windows.py#L1293-L1416

This code can be entered directly into sage:


# Plot the window and its frequency response:

N = 21 # window length

from scipy import signal
from scipy.fftpack import fft, fftshift
import matplotlib.pyplot as plt
import numpy as np

window = signal.chebwin(N, at=100)
plt.plot(window)
plt.title("Dolph-Chebyshev window (100 dB)")
plt.ylabel("Amplitude")
plt.xlabel("Sample")
plt.show()

plt.figure()
A = fft(window, 2048) / (len(window)/2.0)
freq = np.linspace(-0.5, 0.5, len(A))
response = 20 * np.log10(np.abs(fftshift(A / abs(A).max())))
plt.plot(freq, response)
plt.axis([-0.5, 0.5, -120, 0])
plt.title("Frequency response of the Dolph-Chebyshev window (100 dB)")
plt.ylabel("Normalized magnitude [dB]")
plt.xlabel("Normalized frequency [cycles per sample]")
plt.show()


The relevant fragment of the implementation (with added comments) of signal.chebwin in scipy looks 
like this:

order = M - 1.0
beta = np.cosh(1.0 / order * np.arccosh(10 ** (np.abs(at) / 20.)))
k = np.r_[0:M] * 1.0
x = beta * np.cos(np.pi * k / M)
p = np.zeros(x.shape)

# when x > 1, use cosh/acosh:
p[x > 1] = np.cosh(order * np.arccosh(x[x > 1]))

# when x < -1, use cosh/acosh or 0, depending on order being even or odd:
p[x < -1] = (1 - 2 * (order % 2)) * np.cosh(order * np.arccosh(-x[x < -1]))

# when |x| <= 1, use cos/acos:
p[np.abs(x) <= 1] = np.cos(order * np.arccos(x[np.abs(x) <= 1]))

# Appropriate IDFT and filling up depending on even/odd M
if M % 2:                                             # if length M is odd
    w = np.real(fftpack.fft(p))                       #     do the FFT: w = fft(p) - why no ifft?
    n = (M + 1) // 2                                  #     compute shift amount
    w = w[:n]                                         #     ääähh - what does this do?
    w = np.concatenate((w[n - 1:0:-1], w))            #     apply circular shift by n?
else:                                                 # else (i.e. M is even)
    p = p * np.exp(1.j * np.pi / M * np.r_[0:M])      #     additional modulation of p required?
    w = np.real(fftpack.fft(p))                       #     w = fft(p)
    n = M // 2 + 1                                    #     compute shift amount
    w = np.concatenate((w[n - 1:0:-1], w[1:n]))       #     apply circular shift
*/




//=================================================================================================

template<class TSig, class TPar>
void rsStateVectorFilter<TSig, TPar>::setPoles(CRPar p1re, CRPar p1im, CRPar p2re, CRPar p2im)
{
  xx = p1re;
  yx = p1im;
  yy = p2re;
  xy = p2im;
  makePolesDistinct();  // avoid singular case of equal poles by fudging
}

template<class TSig, class TPar>
void rsStateVectorFilter<TSig, TPar>::setImpulseResponseStart(TPar h[3])
{
  TPar A[2][2];
  TPar c[2];
  A[0][0] = xx + xy;
  A[0][1] = yx + yy;
  A[1][0] = xx*A[0][0] + xy*A[0][1];
  A[1][1] = yx*A[0][0] + yy*A[0][1];
  rsLinearAlgebra::rsSolveLinearSystem2x2(A, c, &h[1]);
  cx = c[0];
  cy = c[1];
  ci = h[0] - cx - cy;

  // The first 3 impulse response samples of this filter can be calculated as:
  // h[0] = cx + cy + ci
  // h[1] = cx*(xx + xy) + cy*(yx + yy)
  // h[2] = cx*(xx*(xx+xy) + xy*(yx+yy)) + cy*(yx*(xx+xy) + yy*(yx+yy))
  // and this must match what is given to this function. From this, we can set up a system of
  // 3 linear equations for the 3 mixing coefficients and solve it. It can actually be solved as a 
  // 2x2 system for cx,cy and the 3rd equation for ci is then trivially ci = h[0]-cx-cy.

  // Problem: The A-matrix may become singular. This happens when we have two equal real poles as
  // happens in a cookbook lowpass with q = 0.5.
  // ...then, the filter can't be expressed as a parallel connection of two first order
  // filters - what to do in this case? ...
  // maybe work with a fudged determinant, like D_f = sign(D) * max(abs(D), fudgeFactor) where
  // fudgefactor is some small nonzero number
  // rewrite the equation system as:
  // |a b| * |x1| = |y1|
  // |c d|   |x2|   |y2|
  // and use D = a*d - b*c as preliminary determinant and then call D = fudgeDet(D) or something
  // ...very hacky/dirty/ugly but it may work (the filter may be slightly misadjusted in these
  // problematic cases but otherwise work well). don't call rsSolveLinearSystem2x2, write out the
  // equations here

  // or maybe use a fudgePoles function called at the end of setPoles to make sure, they are
  // somewhat distinct - yep, that's what's currently being done
}

template<class TSig, class TPar>
void rsStateVectorFilter<TSig, TPar>::setupFromBiquad(
  CRPar b0, CRPar b1, CRPar b2, CRPar a1, CRPar a2)
{
  // compute and set up poles from a1, a2:
  TPar d = TPar(0.25)*a1*a1 - a2;  // discriminant, term inside sqrt
  TPar p = -a1*TPar(0.5);          // -p/2 in p-q-formula, term before +/- sqrt
  if(d >= 0) {                     // d >= 0: we have 2 real poles
    TPar sq = rsSqrt(d);
    setPoles(p+sq, 0, p-sq, 0);
  }
  else {                           // d < 0: we have complex conjugate poles
    TPar sq = rsSqrt(-d);
    setPoles(p, sq, p, -sq);
  }

  // compute first 3 samples of biquad impulse response and set up mixing coeffs:
  TPar h[3]; 
  h[0] = b0;
  h[1] = b1 - a1*h[0];
  h[2] = b2 - a1*h[1] - a2*h[0];
  setImpulseResponseStart(h);
}

template<class TSig, class TPar>
void rsStateVectorFilter<TSig, TPar>::makePolesDistinct()
{
  TPar minDelta = TPar(0.001); // actually  2 * minimum pole-distance
    // with 0.001, impulse response of a cookbook lowpass (1kHz@44.1kHz, q=0.5) is visually 
    // indistiguishable from the correct, desired one and the mixing coeffs are around 9.
    // With 0.0001, mixing coeffs raise to about 90. Maybe make more formal tests, 
    // how  the upper bound of the mixing coeffs behaves as function of minDelta and choose 
    // something reasonable. The c-coeffs should ideally stay of the same order of magnitude as in
    // the naturally nonsingular cases (or at least not be vastly higher)
    // maybe plot error between desired and approximated impulse response

  TPar delta;
  if(xy == 0) { // two real poles
    delta = xx - yy;
    if(rsAbs(delta) < minDelta) {
      TPar avg = TPar(0.5)*(xx+yy);
      TPar d2  = TPar(0.5)*minDelta;
      if(xx >= yy) {
        xx = avg + d2;
        yy = avg - d2;
      }
      else {
        xx = avg - d2;
        yy = avg + d2;
      }
    }
  }
  else {
    delta = xy + yx; // + because they have opposite signs
    if(rsAbs(delta) < minDelta) {
      TPar d2  = TPar(0.5)*minDelta;
      if(xy > 0) {
        xy = +d2;
        yx = -d2;
      }
      else {
        xy = -d2;
        yx = +d2;
      }
    }
  }

  // todo: optimize: avoid the makePolesDistinct when it can be predicted that they are
  // distinct anyway (after computing sqrt(d) or sqrt(-d)) ...but not in prototype
  // but what if the shifting of poles makes an originally stable filter unstable? maybe a more
  // sophisticated approach is needed - and unit tests

}

// explicit instantiations:
template class rsStateVectorFilter<double, double>;
template class rsStateVectorFilter<rsFloat64x2, double>;
template class rsStateVectorFilter<float, float>;
template class rsStateVectorFilter<rsFloat32x4, float>;
template class rsStateVectorFilter<rsFloat32x4, rsFloat32x4>;

// can the mixing coeffs be computed more simply - from the b-coeffs of the biquad?

// todo: make functions setupFromPolesAndZeros - i think, the update matrix is (re,-im;im,re)
// setPoles(p1re, p1im, p2re, p2im)
// setImpRespStart(h[3])
// maybe refactor setupFromBiquad into pole (update matrix) and zero (mixing coeffs) computation
// in order to bypass the complicated pole computations when possible
// maybe also factor out a setImpulseResponse function which takes the desied first 3 samples
// of the impulse response - this may be useful when we already know the poles but want to match
// a filter that is not necessarily a biquad (for example, a state-variable filter)

// setupFromParameters - for various parametrizations
// like setupFromFrqDecAmpPhs (for the "modal" filter - take care to take into account the 
// injection into both parts - shifts phase by 45° and increases amplitude by sqrt(2))
// setFreqAndReso(freq, reso) 
// freq controls pole angle, reso is pole radius - maybe allow values > 1 and use saturation
// circle-saturation/distortion or (soft)clipping

// is it possible to find simpler formulas? it looks all very complicated - maybe from equating the
// transfer functions of biquad and svf? i think, the TF of this filter is:
// X(z) = 1 / (1 - xx * z^-1 * X(z) - xy * z^-1 * Y(z))
// Y(z) = 1 / (1 - yx * z^-1 * X(z) - yy * z^-1 * Y(z))
// H(z) = ci + cx*X(z) + cy*Y(z)
// maybe look here for the transfer function of a state-space filter:
// https://ccrma.stanford.edu/courses/250a-fall-2003/hiqfilters.pdf

// todo: make unit tests that compare outputs of biquads for various settings all combinations 
// of real/complex poles/zeros, also with unused higher order coeffs (like b2=0, etc)
// against the outputs of the state vector (and also state variable) implementation
// ->also good to have in place when we later change a-coeff sign convention

// maybe try a nonlinearity: multiply x and y by 1/(1 + x^2 + y^2) after state update
// saturates(with some foldover)/contracts state vector without changing its angle
// maybe apply this factor also to "in" because it would be weird to pass the input through
// undistorted ...but might be interesting to explore

//=================================================================================================
/*
void rsModalFilterFloatSSE2::setParametersTwoEnvs(
  double w, double A, double p, 
  double att1, double att2, double attB,
  double dec1, double dec2, double decB)
{
  //double A = E;  // preliminary
  // todo: compute energy, optionally divide A by sqrt(energy), optimize computations (some 
  // values are computed 4 times (cw,sw,cp,sp) in the 4 calls below - but do this in production 
  // code):

  RAPT::rsDampedSineFilter(w, (1-attB)*A, att1, p, &b0[0], &b1[0], &a1[0], &a2[0]);
  RAPT::rsDampedSineFilter(w,    attB *A, att2, p, &b0[1], &b1[1], &a1[1], &a2[1]);
  RAPT::rsDampedSineFilter(w, (1-decB)*A, dec1, p, &b0[2], &b1[2], &a1[2], &a2[2]);
  RAPT::rsDampedSineFilter(w,    decB *A, dec2, p, &b0[3], &b1[3], &a1[3], &a2[3]);

  // filters 0 and 1 need a minus sign to make the final getSum work:
  b0[0] = -b0[0], b1[0] = -b1[0];
  b0[1] = -b0[1], b1[1] = -b1[1];
}

void rsModalFilterFloatSSE2::setParameters(double w, double A, double p, 
  double att, double dec, double dw, double dp, double b, double attScl, double decScl)
{
  // 1st attack/decay filter pair:
  double c = 1-b;
  RAPT::rsDampedSineFilter(w-0.5*dw, c*A, att/attScl, p-0.5*dp, &b0[0], &b1[0], &a1[0], &a2[0]);
  RAPT::rsDampedSineFilter(w-0.5*dw, c*A, dec/decScl, p-0.5*dp, &b0[1], &b1[1], &a1[1], &a2[1]);

  // 2nd attack/decay filter pair:
  RAPT::rsDampedSineFilter(w+0.5*dw, b*A, att*attScl, p+0.5*dp, &b0[2], &b1[2], &a1[2], &a2[2]);
  RAPT::rsDampedSineFilter(w+0.5*dw, b*A, dec*decScl, p+0.5*dp, &b0[3], &b1[3], &a1[3], &a2[3]);

  // filters with indices 0 and 2 are for the attacks and need to be subtracted in final sum:
  b0[0] = -b0[0], b1[0] = -b1[0];
  b0[2] = -b0[2], b1[2] = -b1[2];

  // maybe dw, dp should also be scaled by b(lend) and not just by 0.5
}
*/


// Idea have a ModalSynth class that lets the user insert different types of mode filters
// simple decaying sines, attack/decay-sines, 4-env-sines, nonlinear modes (perhaps with amplitude
// dependent frequency - they should have a second "sidechain" input where we feed back the total
// summed output - so the nonlinear effects may depend on the total output value


//=================================================================================================


#undef small // this is #defined in some silly windows header - WTF?

template<class T>
void rsMovingQuantileFilterCore<T>::setLengthAndReadPosition(int newLength, int newPosition)
{
  int C = getCapacity();
  rsAssert(newLength   <= C, "Length cannot exceed capacity");
  rsAssert(newLength   >= 2, "If L < 2, we get access violations");
  rsAssert(newPosition >= 1, "If p < 1, we get access violations");

  if( modulatable )
    modulateLengthAndReadPosition(newLength, newPosition);
  else {
    L = newLength;
    p = newPosition;
    dblHp.setData(&small[0], p, C, &large[0], L-p, C);
    reset(); }
}

template<class T>
void rsMovingQuantileFilterCore<T>::modulateLengthAndReadPosition(int newLength, int newPosition)
{
  int C = getCapacity();

  //rsError("not yet implemented");

  int nSo = p;         // old number of small samples
  int nLo = L-p;       // old number of large samples
  L = newLength;
  p = newPosition;
  int nS  = p;         // new number of small samples
  int nL  = L-p;       // new number of large samples

  if(nS+nL == nSo+nLo)
  {
    // total length stays the same - we may need to re-distribute data from the min-heap to the
    // max-heap or vice versa
    int nSi = nSo;                  // current number of small samples
    int nLi = nLo;                  // current number of large samples
    if(nS > nSo) {                  // i think, the if/else is redundant with the whiles 
      while(nSi < nS) {             // number of small samples grows, redistribute from large to small
        moveFirstLargeToSmall(nSo); // move one sample from large heap to small heap
        nSi++; nLi--; }}
    else if(nS < nSo) {
      while(nLi < nL)  {            // number of large samples grows, redistribute from small to large
        moveFirstSmallToLarge(nSo); // move one sample from small heap to large heap
        nSi--; nLi++; }}
    // else: nothing to do at all


    //dblHp.setData(&small[0], p, C, &large[0], L-p, C);
    // we should not need that, should we? yes - seems to make no difference

  }
  else if(nS+nL < nSo+nLo)
  {
    // total length shrinks - we must discard some of the old datapoints in addition to possibly
    // redistribute data

    // preliminary - copied from non-modulatable code:
    dblHp.setData(&small[0], p, C, &large[0], L-p, C); 
    reset();
    int dummy = 0;
  }
  else
  {
    // total length grows - we must fill up missing data in addition to possibly redistribute 
    // data

    // preliminary - copied from non-modulatable code:
    dblHp.setData(&small[0], p, C, &large[0], L-p, C); 
    reset();
    int dummy = 0;
  }

  // In order to increase the length on the fly, we may need to insert some older samples which 
  // we do have stored anywhere here yet - to do so, we may have to use a circular buffer of old 
  // input values that has a length equal to C...maybe that should be done in a subclass - this 
  // buffer could at the same time by used to implement a highpass version - or maybe we store 
  // here a pointer to an rsRingBuffer which is null by default and if it's non-null, we use it 
  // here to fill up the missing past samples and if it's null, we just fill up with zeros or 
  // with the mean of the currently stored old values or maybe use the previous output sample

  // We may need to distinguish between a couple of cases:
  // case 1: length does not change - we need to only redistribute values between the heaps
  //  subcase 1.1: quantile decreases: redistribute from small to large
  //  subcase 1.2: quantile increases: redistribute from large to small
  // case 2: length decreases: we must discard some values (and possibly redistribute)
  // case 3: length increases: we must fill up the heaps (and possibly redistribute), ideally we
  //   should fill up with actual older signal values (to have them available, we may need and 
  //   additional ringbuffer to just store past values), if they are not available, we could use
  //   zeros, of the most recent value or the current quantile - maybe implement variations

  // Can the handling be unified? Or maybe we can reduce it to two cases for each of the two 
  // sub-heaps, i.e. nS grows or shrinks, nL grows or shrinks

}

template<class T>
void rsMovingQuantileFilterCore<T>::moveFirstLargeToSmall(int nSo)
{
  //dblHp.moveFirstLargeToSmall();
  //return;


  //// test - if uncommented, uncomment the reset at the bottom of the function, too:
  //auto swapNodes2 = [&](Node& a, Node& b) { this->swapNodes2(a, b); };
  //dblHp.setSwapFunction(swapNodes2);

  bool ok = isStateConsistent();

  //rsAssert(isStateConsistent(), "state inconsistent");

  // peek the node that we want to move - mainly to figure out, what bufIdx it has stored, so we 
  // can adapt the stored key-value at that bufIdx (the key will change, due to the move):
  Node n1 = dblHp.getSmallestLargeValue();    // this is the node we want to move
  int  bi = n1.bufIdx;                        // buffer index of to be moved node

  // a sanity check:
  int  bc = buf[bi];  // buffer content at index - is the node's key
  //rsAssert(bc < 0, "we expect an index into the large heap");
  // make a function isKeyInLargeHeap(int k), so we can swap between rsDoubleHeap/2

  // store the new key at appropriate position in buffer - the new key is preliminarily the current
  // number of elements in the small heap because heap-insertion first inserts at the end and then 
  // let's the element float up:
  int oldKey1 = buf[bi];
  int newKey1 = dblHp.getNumSmallValues();
  buf[bi] = newKey1;
  ok = isStateConsistent();
  
  /*
  // very experimental - i think, this code should be used only without the 
  // dblHp.decrementLargeIndexOffset() call below, maybe it can even be used with the i-nS 
  // convention:
  Node n2  = dblHp.getLastLargeValue();    // this is the node that gets its key invalidated (right?)
  int  bi2 = n2.bufIdx;
  int oldKey2 = buf[bi2];
  int newKey2 = dblHp.indexToKey(dblHp.getNumSmallValues()+1); // does that make sense?
  //int newKey2 = dblHp.indexToKey(dblHp.getNumSmallValues());
  buf[bi2] = newKey2;   // we must to re-assign the key to something - but what?
  ok = isStateConsistent();
  // something like toKey(dblHp.getNumLargeValues()-2) ..bcs the old key is 
  // toKey(dblHp.getNumLargeValues()-1)...but there already exists some other node with that key
  // but we have just freed the key for the 0th node in the large heap, because we just re-keyed
  // the node with that key - maybe it should get that? so would that be
  // toKey( dblHp.getNumSmallValues()+1 )
  // hmm..maybe that has to be done *after* the extract/insert because this operations does swaps
  // whose effect depends on the content of buf - but i think, we need to figure out which node will 
  // get invalidated before re-shuffling - afterwards, it will be hard to find
  */
  


  // move the first node of the large heap into the small heap:
  Node n = dblHp.large.extractFirst();  // shuffles large heap and buf
  //ok = isStateConsistent();// triggers key/index out of range but not itself
  //fixInconsistentBufferKeys(n);   // test
  int  i = dblHp.small.insert(n);       // shuffles small heap and buf
  ok = isStateConsistent();
  rsAssert(i  == 0);                    // it should end up at the front
  //ok &= isNodeConsistent(n);
  // bcs extractFirst already shuffles buf, we need to first peek the first element of the large
  // heap to obtain its stored bufIndex there - after the call to extract, it's invalid

  // all large-heap indices (not the keys!) get decremented by one, bcs. they move one index 
  // position forward - this is realized as a global offset stored in the dblHp object:
  //dblHp.decrementLargeIndexOffset();


  //// test - i have the hunch that we need some sort of dummy-swap between extract and insert:
  //Node nodeLarge = dblHp.large.extractFirst(); // creates a copy
  //Node nodeSmall = dblHp.small[0];
  //swapNodes(nodeSmall, nodeLarge);  // a swap with dummy nodes to update the buffer
  //dblHp.small.insert(nodeLarge);
  //// nope! that doesn't work either!



  makeBufferConsistent();
  ok = isStateConsistent();

  int dummy = 0;

  // We did the move by hand. Alternatively, we could have called  dblHp.moveFirstLargeToSmall().
  // This should perhaps be done, when it all works and then the friend declaration in
  // rsDoubleHeap should be removed (we access a protected member in the by-hand move).


  // ...that should be it - but it doesn't seem to work. :-( i suspect, that in the last operation,
  // the actual move, the keys in the buf are nor correctly re-shuffled. Maybe the swapNodes call
  // that occurs during the node extraction is wrong? Maybe the dblHp.decrement.. call should be 
  // done after movement? Should that make any difference, i.e. is the offset used during 
  // extract/insert? Nope!

  // perhaps, the exctract/insert should also be realized by a sort of swap?
  // -in the large heap, swap the first node with the last
  // -insert a dummy node at the end of the small heap that is a copy of the first small node
  // -swap them
  // -update heap sizes
  // -do the floating business
  // -or maybe between extract and insert, perform a swap with dummy nodes that are copies
  //  of the first small and first large node

  // with decrementing the 1,3 filter switch to a 2,2 filter, we don't get index errors anymore, 
  // but the output signal is wrong. also, when switching again to 3,1, we get index errors again


  // Some more sanity checks:

  /*
  //Node n2 = dblHp.atKey(i);                         // should be same as n1
  Node n2 = dblHp.atIndex(i);                         // should be same as n1
  rsAssert(n1 == n2);
  int bi2 = n2.bufIdx; 
  int bc2 = buf[bi2];
  rsAssert(bc2 == i);
  */

  //rsAssert(isIndexPermutation(&buf[0], L));
  // doe not work anymore when we use keys

  rsAssert(isStateConsistent(), "state inconsistent");

  rsAssert(dblHp.isDoubleHeap());

  // the buffer buf stores heap-indices - but the heap-index of a node changes when we move
  // it from the large to the small heap (it gets smaller by the old nS value)
  // i think, we need to subtract nSo from buf[bi]
  //buf[bi] -= nSo;  
  // is that correct? or should we do that before calling  dblHp.moveFirstLargeToSmall(); does
  // that make a difference? or will it mess up the swapping...it's a mess!!

  // aftermoving the node over, we must change the stored key in the buf to make it now refer to
  // a positive index

  // but no - the strategy withe the key doesn't work either because all large buffer keys become 
  // invalid after the move - actually it seems, we should decrement them all by one, but: that's
  // actually the same as incrementing an subtractive offset - which is preciely what we did before
  // with the i-nS convention - but maybe we need both: a flag (per key) to indicate the large heap 
  // and a global offset for the large heap keys - here, the global offset gest decremented


  // maybe to streamline it, copy the code from moveFirstLargeToSmall() to here and maybe
  // simplify, if possible


  //auto swapNodes = [&](Node& a, Node& b) { this->swapNodes(a, b); };
  //dblHp.setSwapFunction(swapNodes);
}
// this does not work yet - i think, we may have to do something with the stored bufIndex
// values, too
// no - the problem is that we use nS as offset in rsDoubleHeap: whenever an index i >= nS, we
// interpret that index as i-nS into the large heap. but when nS changes, everything gets messed
// up. We need to adopt some convention that is independent from nS, for example: whenever the 1st
// bit (of an unsigned int) is one, mask that away and use the remaining bits as actual index.
// i think, we can do this in an rsDoubleHeap2 subclass

template<class T>
void rsMovingQuantileFilterCore<T>::moveFirstSmallToLarge(int nSo)
{
  dblHp.moveFirstSmallToLarge();

  // this seems to magically work already, although it shouldn't - we don't update the buffer
  // content and don't increment the offset yet...more tests needed...
}

/*
template<class T>
int rsMovingQuantileFilterCore<T>::getNodeKey(rsMovingQuantileFilterCore<T>::Node n)
{
  for(int i = 0; i < dblHp.getSize(); i++)
  {
    Node n2 = dblHp.atIndex(i);
    if(n2 == n)
      return dblHp.indexToKey(i);
  }
  return 666; // indicate failure
}

template<class T>
void rsMovingQuantileFilterCore<T>::fixInconsistentBufferKeys(
  rsMovingQuantileFilterCore<T>::Node n)
{
  Node n2;
  int i, k;
  while(1)
  {
    i  = n.bufIdx;         // buffer index stored in node n
    k  = buf[i];           // node key stored in buffer at i
    n2 = dblHp.atKey(k);   // retrieve the node at key k
    if(n2 == n)
    {
      break; // we found a node that is consistent
    }
    else
    {
      // node is inconsistent - fix the key stored in the buffer and continue with the new node:
      buf[i] = k; // getKey(n)...but figuring out a node's key is O(n) 
      n = n2;
    }
  }
  int dummy = 0;
}
*/

// modifies the values in the buffer to make them consistent
template<class T>
void rsMovingQuantileFilterCore<T>::makeBufferConsistent()
{
  for(int i = 0; i < dblHp.getNumValues(); i++)
  {
    int  k  = dblHp.indexToKey(i);
    Node n  = dblHp.atIndex(i);
    int  bi = n.bufIdx;
    buf[bi] = k;
  }
}

// maybe hae a similar function that establishes consistency by modifying the stored bufIdx fields
// in the nodes in the dblHp

template<class T>
bool rsMovingQuantileFilterCore<T>::isStateConsistent()
{
  bool r = true;

  // Check that all nodes in the double-heap have the correct back-link to their index in the 
  // delay-buffer. This back-link is needed to update the delay-buffer when nodes in the heap are
  // swapped during floatUp/Down:
  for(size_t i = 0; i < L; i++)
    r &= isBufferSlotConsistent((int)i);

  // Check that each key occurs in the buf exactly once:
  std::vector<int> tmp(L);
  for(size_t i = 0; i < L; i++)
    tmp[i] = dblHp.keyToIndex(buf[i]);
  r &= isIndexPermutation(&tmp[0], (int) tmp.size());

  // todo: 
  // -check that the sum of the heap sizes matches the buffer size
  // -check that each buffer index occurs in the heap exactly once

  return r;
}

template<class T>
bool rsMovingQuantileFilterCore<T>::isNodeConsistent(
  const rsMovingQuantileFilterCore<T>::Node& n)
{
  bool r = true;

  for(size_t i = 0; i < buf.size(); i++) {
    int  k  = buf[i];         // key of node in i-th buffer slot
    Node n2 = dblHp.atKey(k); // retrieve the node
    if(n2 == n)
      r &= n.bufIdx == i; }

  return r;
}

template<class T>
bool rsMovingQuantileFilterCore<T>::isBufferSlotConsistent(int i)
{
  int  k = buf[i];         // key of node in i-th buffer slot
  Node n = dblHp.atKey(k); // retrieve the node
  bool result = n.bufIdx == i;
  return result;
}


//=================================================================================================

rsGroupString rsGroupString::inverse() const
{
  size_t len = s.size();
  std::vector<unsigned int> r(len);
  for(int i = 0; i < len; i++)
    r[i] = s[len-1-i];
  return r;
}

rsGroupString rsGroupString::operator+(const rsGroupString &rhs) const
{
  std::vector<unsigned int> r = s;
  for(int i = 0; i < rhs.s.size(); i++) {
    if(r.size() > 0 && r[r.size()-1] != rhs.s[i])
      r.push_back(rhs.s[i]);
    else if(r.size() > 0)     // avoid popping on empty vector
      r.pop_back();
  }
  return r;
  // todo: check, if this implementation actually does the right thing.
  // a more naive implementation would be: just concatenate and then call removeDoublets

  // as a generalization, we could disallow certain pairings - at the moment, the pairings aa, bb, 
  // cc, are disallowed - we could define a bijective function from the character-set to itself 
  // that defines an "annihilator" for each character. at the moment, each character is it own
  // annihilator - maybe use pairings that are unnatural for an actual language
}



rsGroupString2::rsGroupString2(const char* inStr)
{
  size_t len = strlen(inStr);
  s.resize(len);
  for(size_t i = 0; i < len; i++)
    s[i] = inStr[i] - 97;
  rsAssert(checkCharacters());
}

rsGroupString2::rsGroupString2(const rsGroupString& gs)
{
  //s = gs.s;  // compiler error: cannot access protected member - why?
  s = gs.get();
  rsAssert(checkCharacters());
}

std::string rsGroupString2::toString() const
{
  std::string outStr;
  size_t len = s.size();
  outStr.resize(len);
  for(size_t i = 0; i < len; i++)
    outStr[i] = (char)(s[i] + 97);
  rsAssert(checkCharacters());
  return outStr;
}

bool rsGroupString2::checkCharacters() const
{
  bool r = true;
  for(int i = 0; i < s.size(); i++)
    r &= s[i] >= 0 && s[i] <= 25;
  return r;
}

// Maybe we need to let 'a' map to 1 instead of mapping to 0 - otherwise "" and "a" would both be 
// neutral with respect to multiplication if we define multiplication via modular addition over
// the characters...or maybe the string "a" is the neutral element of multiplication? there 
// wouldn't be ambiguity because "aa", "aaa" etc. are not allowed anyway - but what about 
// "a" * ""? would that be "a" - should it be ""? ...see field axioms - no, must be "", see below
// A*0 = 0 in every ring
// or maybe we can at least get a ring if a field is not possible?


// ringe (weitz):
// https://www.youtube.com/watch?v=tpSFUXPyy0c&index=38&list=PLb0zKSynM2PD4-kkRIAWFdnFivbhEgfeO
// gruppen:
// https://www.youtube.com/watch?v=L_G3Lbo5z60&list=PLb0zKSynM2PD4-kkRIAWFdnFivbhEgfeO&index=21
// https://www.youtube.com/watch?v=fYArC336jd4&index=26&list=PLb0zKSynM2PD4-kkRIAWFdnFivbhEgfeO
// körper:
// https://www.youtube.com/watch?v=wnUdh5iMoz4&list=PLb0zKSynM2PD4-kkRIAWFdnFivbhEgfeO&index=45
// ah - das hilft auch bei der frage: soll "a" oder "" neutrales element der multiplikation sein
// (see 3:10 - "die null von hier oben kann man hier nicht mit rein nehmen")
// @9:05: bei multiplikation mit null ("") muß null rauskommen 
// soo: multiply: 
// -the result has the size (num chars) of the larger of the two inputs
// -each character is given by (modular) addition of the chars of the two inputs (assuming to add
//  0, for the part where there are no chars anymore in one of the inputs
//  abc * ab = bdc ...does that work out? what about
//  abd * ab = bdd = bd ? but bd / ab = ab != abc ...or maybe don't apply the 
//  "no-doublets" rule in multiplication? ...or maybe there should not actually be a strict 
//  "no-doublets" rule at all and we just leave the addition with that push-or-pop rule as 
//  implemented? ...this would be even better because it doesn't restrict our allowed strings
// -check distributive A*(B+C)=A*B+A*C, (A+B)*C = A*C+B*C (check, if the 2nd must hold or only
//  in abelian groups) and associative law...
// -maybe it would be better to reverse roles of add/mul...because this mul (if it works) is
//  commutative...but then the distributive law would have to work the other way around...but
//  commutativity of addition is more basic - it works for matrices as well (where mul is not
//  commutative)