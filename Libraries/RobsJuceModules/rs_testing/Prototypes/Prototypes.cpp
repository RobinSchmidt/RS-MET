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
#include "OdeSolver.cpp"
#include "AdditiveSynthEngine.cpp"

//#include "SinusoidalModeling.cpp" // moved to rapt

using namespace RAPT;

std::vector<double> splineSlopes(const std::vector<double>& x, const std::vector<double>& y,
  bool prescribe2ndDeriv, double ypStart, double ypEnd)
{
  int N = (int) x.size();
  rsAssert((int)y.size() == N);

  std::vector<double> s(N);       // the slopes


  return s;
}


void rsSinCos1(double x, double* s, double* c)
{
  // it seems like this code is slower than std::cos/sin
  // taken from: http://lab.polygonal.de/2007/07/18/fast-and-accurate-sinecosine-approximation/
  // low precision version:

  // always wrap input angle to -PI..PI:
  if (x < -3.14159265)
    x += 6.28318531;
  else
    if (x >  3.14159265)
      x -= 6.28318531;

  // compute sine:
  if (x < 0)
    *s = 1.27323954 * x + 0.405284735 * x * x;
  else
    *s = 1.27323954 * x - 0.405284735 * x * x;

  // compute cosine: sin(x + PI/2) = cos(x)
  x += 1.57079632;
  if (x >  3.14159265)
    x -= 6.28318531;
  if(x < 0)
    *c = 1.27323954 * x + 0.405284735 * x * x;
  else
    *c = 1.27323954 * x - 0.405284735 * x * x;
}

// approximations for arguments rdeuced to the range -0.5*pi..0.5*pi
double rsCosApproxReduced1(double x)
{
  constexpr double a2 = -4.0 / (PI*PI);
  return 1.0 + a2 * x*x;
}
//var("x a0 a2")
//f(x) = a0 + a2*x^2
//eq1  = f(0)  == 1         # 1st requirement
//eq2  = f(pi) == 0         # 2nd requirement
//solve([eq1,eq2],[a0,a2])

double rsSinApproxReduced1(double x)
{
  constexpr double pi2 = PI*PI;
  constexpr double a1 =  3.0 /  PI;
  constexpr double a3 = -4.0 / (PI*pi2);
  double x2 = x*x;
  return a1*x + a3*x*x2;
}
//var("x a1 a3")
//f(x)  = a1*x + a3*x^3
//f1(x) = diff(f(x), x, 1)   # 1st derivative
//eq1 = f( pi/2) == 1        # 1st requirement
//eq2 = f1(pi/2) == 0        # 2nd requirement
//solve([eq1,eq2],[a1,a3])
//
// [[a1 == 3/pi, a3 == -4/pi^3]]


// Controls only the 2nd derivative:
double rsCosApproxReduced2(double x)
{
  constexpr double pi2 = PI*PI;
  constexpr double a2 = -24.0 / (5.0 * pi2);
  constexpr double a4 = +16.0 / (5.0 * pi2*pi2);
  double x2 = x*x;
  return 1.0 + a2*x2 + a4*x2*x2;
}
//var("x a0 a2 a4")
//f(x)  = a0 + a2*x^2 + a4*x^4
//f2(x) = diff(f(x), x, 2)   # 2nd derivative
//eq1 = f(0)     == 1        # 1st requirement
//eq2 = f(pi/2)  == 0        # 2nd requirement
//eq3 = f2(pi/2) == 0        # 3rd requirement
//solve([eq1,eq2,eq3],[a0,a2,a4])

// Controls the 1st derivative at pi/2:
double rsCosApproxReduced2_all(double x)
{
  constexpr double pi2 = PI*PI;
  constexpr double a2 =    (PI-8) / (pi2);
  constexpr double a4 = -4*(PI-4) / (pi2*pi2);
  double x2 = x*x;
  return 1.0 + a2*x2 + a4*x2*x2;
}
//var("x a0 a2 a4")
//f(x)  = a0 + a2*x^2 + a4*x^4
//f1(x) = diff(f(x), x, 1)    # 1st derivative
//eq1 = f(0)     ==  1        # 1st requirement
//eq2 = f(pi/2)  ==  0        # 2nd requirement
//eq3 = f1(pi/2) == -1        # 3rd requirement
//solve([eq1,eq2,eq3],[a0,a2,a4])

// Controls 2nd and 4th derivative at pi/2:
double rsCosApproxReduced3(double x)
{
  constexpr double pi2 = PI*PI;
  constexpr double pi4 = pi2*pi2;
  constexpr double a2  = -300.0 / (61.0 * pi2);
  constexpr double a4  = +240.0 / (61.0 * pi2*pi2);
  constexpr double a6  = -64.0  / (61.0 * pi2*pi4);
  double x2 = x*x;
  double x4 = x2*x2;
  return 1.0 + a2*x2 + a4*x2*x2 + a6*x2*x4;
}
//var("x a0 a2 a4 a6")
//f(x)  = a0 + a2*x^2 + a4*x^4 + a6*x^6
//f2(x) = diff(f(x), x, 2)  # 2nd derivative
//f4(x) = diff(f(x), x, 4)  # 4th derivative
//eq1 = f(0)     == 1       # 1st requirement
//eq2 = f(pi/2)  == 0       # 2nd requirement
//eq3 = f2(pi/2) == 0       # 3rd requirement
//eq4 = f4(pi/2) == 0       # 4th requirement
//solve([eq1,eq2,eq3,eq4],[a0,a2,a4,a6])

// Controls the 1st and 2nd derivative at pi/2:
double rsCosApproxReduced3_all(double x)
{
  constexpr double pi2 = PI*PI;
  constexpr double pi4 = pi2*pi2;
  constexpr double a2 =  3*(3*PI-16) / (4*pi2);
  constexpr double a4 = -2*(7*PI-24) / (pi2*pi2);
  constexpr double a6 =  4*(5*PI-16) / (pi2*pi4);
  double x2 = x*x;
  double x4 = x2*x2;
  return 1.0 + a2*x2 + a4*x4 + a6*x2*x4;
}
//var("x a0 a2 a4 a6")
//f(x)  = a0 + a2*x^2 + a4*x^4 + a6*x^6
//f1(x) = diff(f(x), x, 1)    # 1st derivative
//f2(x) = diff(f(x), x, 2)    # 2nd derivative
//eq1 = f(0)     ==  1        # 1st requirement
//eq2 = f(pi/2)  ==  0        # 2nd requirement
//eq3 = f1(pi/2) == -1        # 3rd requirement
//eq4 = f2(pi/2) ==  0        # 4th requirement
//solve([eq1,eq2,eq3,eq4],[a0,a2,a4,a6])

// Controls the 2nd, 4th and 6th derivative at pi/2:
double rsCosApproxReduced4(double x)
{
  constexpr double pi2 = PI*PI;
  constexpr double pi4 = pi2*pi2;
  constexpr double a2  = -6832.0 / (1385.0 * pi2);
  constexpr double a4  = +1120.0 / (277.0  * pi2*pi2);
  constexpr double a6  = -1792.0 / (1385.0 * pi2*pi4);
  constexpr double a8  = +256.0  / (1385.0 * pi4*pi4);
  double x2 = x*x;
  double x4 = x2*x2;
  return 1.0 + a2*x2 + a4*x2*x2 + a6*x2*x4 + a8*x4*x4;
}
//var("x a0 a2 a4 a6 a8")
//f(x)  = a0 + a2*x^2 + a4*x^4 + a6*x^6 + a8*x^8
//f2(x) = diff(f(x), x, 2)  # 2nd derivative
//f4(x) = diff(f(x), x, 4)  # 4th derivative
//f6(x) = diff(f(x), x, 6)  # 6th derivative
//eq1 = f(0)     == 1       # 1st requirement
//eq2 = f(pi/2)  == 0       # 2nd requirement
//eq3 = f2(pi/2) == 0       # 3rd requirement
//eq4 = f4(pi/2) == 0       # 4th requirement
//eq5 = f6(pi/2) == 0       # 5th requirement
//solve([eq1,eq2,eq3,eq4,eq5],[a0,a2,a4,a6,a8])

// Controls the 1st, 2nd and 3rd derivative at pi/2:
double rsCosApproxReduced4_all(double x)
{
  constexpr double pi2 = PI*PI;
  constexpr double pi3 = PI*pi2;
  constexpr double pi4 = pi2*pi2;
  constexpr double a2 =  ( 1./96)*(348*PI - pi3 - 1536)/pi2;
  constexpr double a4 =  (-1./8 )*(244*PI - pi3 -  768)/pi4;
  constexpr double a6 =  ( 1./2 )*(172*PI - pi3 -  512)/(pi2*pi4);
  constexpr double a8 =  (-2./3 )*(132*PI - pi3 -  384)/(pi4*pi4);
  double x2 = x*x;
  double x4 = x2*x2;
  return 1.0 + a2*x2 + a4*x4 + a6*x2*x4 + a8*x4*x4;
}
// [[a0 == 1, 
// a2 == 1/96*(348*pi - pi^3 - 1536)/pi^2,
// a4 == -1/8*(244*pi - pi^3 - 768)/pi^4,
// a6 == 1/2*(172*pi - pi^3 - 512)/pi^6, 
// a8 == -2/3*(132*pi - pi^3 - 384)/pi^8]]

//var("x a0 a2 a4 a6 a8")
//f(x)  = a0 + a2*x^2 + a4*x^4 + a6*x^6 + a8*x^8
//f1(x) = diff(f(x), x, 1)    # 1st derivative
//f2(x) = diff(f(x), x, 2)    # 2nd derivative
//f3(x) = diff(f(x), x, 3)    # 3rd derivative
//eq1 = f(0)     ==  1        # 1st requirement
//eq2 = f(pi/2)  ==  0        # 2nd requirement
//eq3 = f1(pi/2) == -1        # 3rd requirement
//eq4 = f2(pi/2) ==  0        # 4th requirement
//eq5 = f3(pi/2) ==  1        # 5th requirement
//solve([eq1,eq2,eq3,eq4,eq5],[a0,a2,a4,a6,a8])


double rsCosApproxReduced5(double x)
{
  constexpr double pi2 = PI*PI;
  constexpr double pi4 = pi2*pi2;
  constexpr double pi6 = pi2*pi4;
  constexpr double a2  = -249300 / (50521 * pi2    );
  constexpr double a4  = +204960 / (50521 * pi2*pi2);
  constexpr double a6  = -67200  / (50521 * pi6    );
  constexpr double a8  = +11520  / (50521 * pi4*pi4);
  constexpr double a10 = -1024   / (50521 * pi4*pi6);
  double x2 = x*x;
  double x4 = x2*x2;
  double x8 = x4*x4;
  return 1.0 + a2*x2 + a4*x2*x2 + a6*x2*x4 + a8*x8 + a10*x2*x8;
}
//var("x a0 a2 a4 a6 a8 a10")
//f(x)  = a0 + a2*x^2 + a4*x^4 + a6*x^6 + a8*x^8 + a10*x^10
//f2(x) = diff(f(x), x, 2)  # 2nd derivative
//f4(x) = diff(f(x), x, 4)  # 4th derivative
//f6(x) = diff(f(x), x, 6)  # 6th derivative
//f8(x) = diff(f(x), x, 8)  # 8th derivative
//eq1 = f(0)     == 1       # 1st requirement
//eq2 = f(pi/2)  == 0       # 2nd requirement
//eq3 = f2(pi/2) == 0       # 3rd requirement
//eq4 = f4(pi/2) == 0       # 4th requirement
//eq5 = f6(pi/2) == 0       # 5th requirement
//eq6 = f8(pi/2) == 0       # 6th requirement
//solve([eq1,eq2,eq3,eq4,eq5,eq6],[a0,a2,a4,a6,a8,a10])

// It seems like controlling only even derivatives for the cosine gives a larger maximum error but
// the shape at the junction seems to be smoother compared to controlling even and odd derivatives.
// Yet another option would be to control function values at various points, like all zeros and/or
// maxima and minima. Maybe try some more options. The goal is not necessarily to minimize the 
// maximum error but rather to minimize the audible artifacts when using such an approximation as
// sine oscillator, so maybe as error measure we should look at the total harmonic distortion in 
// the spectrum and optimize that. Maybe try it with a frequency that visits the wrap-around points
// often and closely.
// https://www.desmos.com/calculator/6hvksqmtgl

// Approximation for arguments reduced to the range -1.5*pi..1.5*pi using internally the functions
// for the range -0.5pi..0.5pi and some range reduction that can (hopefully) be imlemented in a 
// branchless way:
double rsCos2(double x)
{
  double xa = rsAbs(x);

  // This implements the range reduction from +-1.5pi to +-0.5pi. The double(xa > 0.5*PI) can later
  // be done with a branchless simd instruction that will contains ones where the condition is true
  // and zeros where the condition is false. The computation of xr is conceptually a linear 
  // interpolation with coeffs 0 and 1 between the original and shifted argument:
  double o  = double(xa > 0.5*PI);     // |x| is outside range -> reflect to inside
  double xr = xa - o*PI;               // (1-o)*xa + o*(xa-PI)

  // Call the function for the range-reduced cosine (for the range -pi/2...pi/2):
  //double y = rsCosApproxReduced1(xr);     // error below 0.1     -> SNR > 20 dB

  //double y = rsCosApproxReduced2(xr);     // error below 0.01    -> SNR > 40 dB
  //double y = rsCosApproxReduced2_all(xr);     // max error less than above but less smooth

  //double y = rsCosApproxReduced3(xr);       // error below 0.001    -> SNR > 60 dB
  //double y = rsCosApproxReduced3_all(xr); // error below 0.0001   -> SNR > 80 dB

  //double y = rsCosApproxReduced4(xr);       // error below 0.0001  -> SNR > 80 dB
  double y = rsCosApproxReduced4_all(xr); // error below 0.00001  -> SNR > 100 dB

  //double y = rsCosApproxReduced5(xr);       // error below 0.00001  -> SNR > 100 dB

  // Final adjustment to take into account the possible range reduction. This is conceptually a 
  // linear interpolation with coeffs 0 and 1 between original and negated output value:
  y = y - 2*o*y;                       // (1-o)*y + o*(-y)
  return y;

  // For the approximations controlling only even derivatives, the error has qualitatively a shape
  // like cos(x) - cos(3x). Maybe we can create the cos(3x) term by a multiple angle formula and 
  // subtract an appropriately scaled amount. But for the triple angle formula, we also need an 
  // appoximate of the sine...soo that would call for approximating sin/cos jointly

  // ToDo: 
  // -implement approximations of various qualities, maybe with SNR between 60 and 160 dB, the 4
  //  derivative version attains 100 dB...but maybe that can be improved further by numeric 
  //  optimization of the coeffs
  // -check if errors are below expected tolerances
  // -Check, if the odd derivatives that are not controlled are nevertheless matched, i.e. their 
  //  values may be wrong at the junction, but maybe at least they match, so the resulting function
  //  is smooth up to that order. If that's the case, it means, by controlling only even 
  //  derivatives we may achieve a higher degree of smoothness at the expense of larger maximum 
  //  error for the same number of coefficients.
  // -Maybe find a minimax optimized set of coeffs that push the error below the 32bit floating 
  //  point rounding noise. Maybe try to use special polynomials constructed from multiplying
  //  simpler polynomials, like: (a0 + a2*x^2 + a4*x^4) * (b0 + b2*x^2 + b4*x^4) 
  //   = a0*b0 + a0*b2*x^2 + a0*b4*x^4...
  //  this gives an 8th degree polynomial with only 5 multiplications and 4 additions instead
  //  of 8 of each. The price is that the coeffs are coupled. We have less degrees of freedom to
  //  optimize and therefore potentially an error that is greater than for a general optimized
  //  8th degree polynomial

}

void rsSinCos2(double x, double* s, double* c)
{
  double xa = rsAbs(x);
  double xs = rsSign(x);

  double o  = double(xa > 0.5*PI);     // |x| is outside range -> reflect to inside
  double xr = (1.0-o)*xa + o*(PI-xa);
  *s = sin(xa);

  double n  = double(x < 0.0);         // x is negative sine needs to be multiplied by -1
  double fs = 1 - 2*n;                 // factor conceptually: f = (1-n)*1 + n*(-1)
  *s *= fs;

  xr  = (1.0-o)*xa + o*(xa-PI);
  *c  = cos(xa);  // or do we need xr?

  // ToDo:
  // -replace the call to sin/cos with a polynomial approximation (use odd polynomial for sin, even
  //  for cos)
  // -simplify, templatize, vectorize, the <,> operations shall be implemented by some branchless
  //  SIMD functions
  // -factor out the range reduction...maybe templatize on FSin, FCos or use function pointers
  //  for the range-reduced sin/cos functions
  // -factor out the pure cosine - it doesn't even need to compute the sign, so it could be even 
  //  more efficient, the sine can then be obtained by a simple shift - but careful: the shift will 
  //  change the valid domain
  // -implement an efficient sine-oscillator based on that, should also provide a SIMD version, 
  //  generating 4 float sines in parallel -> use it in the sinusoidal modeling framework
}

void rsSinCosApprox4(double x, double* s, double* c)
{
  double xa = rsAbs(x);
  double xs = rsSign(x);

  double o  = double(xa > 0.5*PI);     // |x| is outside range -> reflect to inside
  double xr = (1.0-o)*xa + o*(PI-xa);
  *s = rsCosApproxReduced4(xa);
  // This is still WRONG!!! We need a function rsSinApproxReduced4! but for benchmark, this is ok
  // as a stand in - an actual sine approximation should have similar complexity. But perhaps, 
  // it's better tor approximate sine and cosine jointly

  double n  = double(x < 0.0);         // x is negative sine needs to be multiplied by -1
  double fs = 1 - 2*n;                 // factor conceptually: f = (1-n)*1 + n*(-1)
  *s *= fs;

  xr  = (1.0-o)*xa + o*(xa-PI);
  *c  = rsCosApproxReduced4(xa);  // or do we need xr?
}
// ToDo: 
// -Implement approximations to the exponential function in the interval 1..2 using the same 
//  approach (matching value and derivatives at both endpoints). That makes for a smooth 
//  approximation, when we splice together such segments (scaled by 2^(float-exponent)).
//  Maybe try Chebychev approximation, too (i think, this finds an interpolation polynomial where
//  the evaluation points are given by the roots of a Chebychev polynomial)


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

void rsCircularShift(int* a, int N, int k)
{
  rsAssert(k > 0 && k < N, "other cases are not handled yet");
  using AT = rsArrayTools;
  AT::reverse(a, N);
  AT::reverse(a, k);
  AT::reverse(&a[k], N-k);
}
// if it works, templatize and replace the implementation in rsArrayTools, but keep the old version
// somewhere else...or maybe turn the old version into a workspace-based implementation...and then
// do benchmarks, which one ist faster
// todo: handle cases, where k >= N, k < 0, k <= -N, ... i think, currently, it only works for
// 0 < k < N

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

//-------------------------------------------------------------------------------------------------

const rsUint64 rsModularIntegerNTT::roots[15] =
{
  4179340454199820288,
  3360066027580426122,
  3324705732702508476,
  4093416561646622555,
  4129893269131444668,
  4086220048833014884,
  4075462463776626479,
  4128470768322469725,
  4161514758139463127,
  4178097261067721820,
  4176218026832679610,
  4178374021926307362,
  4177450540047517585,
  4179324170293557359,
  4179136626770643812
};

const rsUint64 rsModularIntegerNTT::rootsInv[15] =
{
  4179340454199820288,
  819274426619394167,
  1324460247237108937,
  2371771477745282871,
  411660242129521289,
  1632617699591430033,
  145924037567102543,
  2873393256140181485,
  3938187127825272781,
  4043479916509059316,
  2571746520142147840,
  3191651436916646396,
  2791086557696422406,
  2638202074896170702,
  3121307387385914037
};

const rsUint64 rsModularIntegerNTT::lengthsInv[15] =
{
  2089670227099910145,
  3134505340649865217,
  3656922897424842753,
  3918131675812331521,
  4048736065006075905,
  4114038259602948097,
  4146689356901384193,
  4163014905550602241,
  4171177679875211265,
  4175259067037515777,
  4177299760618668033,
  3191651436916646396,
  4178830280804532225,
  4179085367502176257,
  4179212910850998273
};

// The modulus 4179340454199820289 was found on this page:
// https://www.programmersought.com/article/8906403181/ 
// and the roots were generated with the sage code:
//
// p = 4179340454199820289 
// Zp = Integers(p)
// root_list = Zp.zeta(2^10, all=False);  # replace the 2^10 by 2^11, etc
// root_list
//
// one at a time, replacing the 2^10 by 2^1, 2^2, ... up to 2^15. There may be even higher order
// roots, but SageMathCell was not able to produce them, probably due to exceeded computation time.
// The result for 2^15 already took about 2 minutes. Once we have the roots themselvese, their 
// inverses and the inverses of the corresponding array lengths can be found one at a time by:
//
// p  = 4179340454199820289      # modulus
// Zp = Integers(p)
// r  = Zp(4179340454199820288)  # current root
// N  = Zp(2^1)                  # current length (exponent = array index + 2)
// ri = 1/r                      # inverse root
// Ni = 1/N                      # inverse length
// r, ri, Ni
//
// and replacing the 4179340454199820288 with the current root and the 1 in 2^1 with the array 
// index of the current root +1. ToDo: find more roots, maybe using a local installation of sage
// or some other math software.

//-------------------------------------------------------------------------------------------------

// The numbers in these tables were created using SageMath with the following script:
//
// k  = 1                        # exponent for length N = 2^k
// p  = 3221225473               # modulus
// Zp = Integers(p)
// r  = Zp.zeta(2^k, all=False); # currentl root
// N  = Zp(2^k)                  # current length (exponent = array index + 2)
// ri = 1/r                      # inverse root
// Ni = 1/N                      # inverse length
// k, N, r, ri, Ni
//
// and running it for one value of k at a time, k is the array index plus 1, r is the root, ri its 
// inverse and Ni the inverse of the length. k,N are only printed out for reference - it's r,ri,Ni
// that are relevant. For values of k higher than 16, the script just takes too long, but there 
// should be even higher order roots - but we need another software to find them. Maybe it's 
// feasible to write a little exhaustive search routing in C++. We have to check up to 2^32 
// numbers per root, so that's 4 billion tests per root which may indeed be feasible to do. 
// ...done - the code is in numberTheoreticTrafoModuli and the results are here:

const rsUint64 rsModularIntegerNTT_64::roots[numRoots] = 
{
  3221225472, 2207278994, 2607818977, 2831384513, 3154294145, 3003263009, 3017814531, 3220790800,
  3207219178, 3220512096, 3217240014, 3221031619, 3220623453, 3220981849, 3221224630, 3221209783,
  3221200472, 3221213459, 3221225125, 3221217292, 3221224362, 3221223061, 3221225387, 3221225384,
  3221225462, 3221225464, 3221225470, 3221225465, 3221225459, 3221225460 
};

const rsUint64 rsModularIntegerNTT_64::rootsInv[numRoots] = 
{
  3221225472, 1013946479, 2190011530, 1626607911,  122509875,  885572692,  511341355, 1442890367,
   261910979, 1771302162,  203364025, 1775783335, 2302174028,  149793384, 1448214062,  931465900,
  1378887045, 1340345109, 2230791204,   66542856, 1954190791, 3062300999, 1985173838,  904838616,
   878516038, 2505397590, 1073741824,  402653184, 2070787804, 1486719449 
};

const rsUint64 rsModularIntegerNTT_64::lengthsInv[numRoots] =
{
  1610612737, 2415919105, 2818572289, 3019898881, 3120562177, 3170893825, 3196059649, 3208642561,
  3214934017, 3218079745, 3219652609, 3220439041, 3220832257, 3221028865, 3221127169, 3221176321,
  3221200897, 3221213185, 3221219329, 3221222401, 3221223937, 3221224705, 3221225089, 3221225281,
  3221225377, 3221225425, 3221225449, 3221225461, 3221225467, 3221225470
};

std::vector<int> rsConvolveNTT(const std::vector<int>& x, const std::vector<int>& h)
{
  // Shorthands for convenience:
  using AT = rsArrayTools;
  using LT = rsLinearTransforms;
  using MI = rsModularIntegerNTT_64;

  // Figure out lengths:
  int Nx = (int) x.size();         // length of input signal
  int Nh = (int) h.size();         // length of impulse response
  int Ny = Nx+Nh-1;                // length of output signal 
  int N  = 2*rsNextPowerOfTwo(Ny); // length of NTT buffer
  int k  = (int)rsLog2(N) - 1;     // index of twiddle factor

  // Prepare NTT buffers:
  std::vector<MI> X(N), H(N);
  AT::convert(&x[0], &X[0], Nx); AT::fillWithZeros(&X[Nx], N-Nx);
  AT::convert(&h[0], &H[0], Nh); AT::fillWithZeros(&H[Nh], N-Nh);

  // Do NTT-based convolution:
  MI WN = MI(MI::rootsInv[k]);        // twiddle factor for forward NTT
  LT::fourierRadix2DIF(&X[0], N, WN); // transform X
  LT::fourierRadix2DIF(&H[0], N, WN); // transform H
  for(int i = 0; i < N; i++)          // spectral multiplication, 
    X[i] = X[i] * H[i];               // ...result goes to X
  WN = MI(MI::roots[k]);              // twiddle factor for inverse NTT
  LT::fourierRadix2DIF(&X[0], N, WN); // unscaled inverse NTT
  MI S(MI::lengthsInv[k]);            // scaler = 1/N
  for(int i = 0; i < N; i++)
    X[i] = S * X[i];                  // this magic scaling turns garbage into result

  // Convert result to output:
  std::vector<int> y(Ny);
  AT::convert(&X[0], &y[0], Ny);
  return y;
}
// For production code, we don't really need to scale the whole X[i] by S. It's enough to do it up
// to Ny and we also can write the result values directly into y. It's really amazing, how the X 
// array looks like total garbage and the final (modular) multiplication by the number S magically
// turns the garbage into our desired convolution result. ...feels almost like dealing with 
// cryptography - the NTT totally is some sort of crypto-FFT! :-) Maybe we can encrypt sequences
// by modular convolution and decrypt by deconvoltion...but wait...modular convolution is not 
// really what we are doing here. We actually compute a regular integer convolution by means of 
// going through modular arithmetic. But the end result is just regular convolved integer array.
// But the whole thing only works for short sequences of small numbers because the overflow 
// problems kick in really soon. That's a bit sad because it means, we can't use it for convolution
// reverb :-( Maybe with a 128 bit integer type, it could become practical? -> try it with 
// rsBigInt.

//-------------------------------------------------------------------------------------------------

template<class TSig, class TPar>
rsRationalFunction<TPar> rsLadderTest<TSig, TPar>::getTransferFunction(bool withGain)
{
  //rsAssert(B1 == 0 || B1 == B0, "Does not yet work in the general case");
  // the formulas below work only for zeroless or zero at z=-1

  using T  = TPar;
  using RF = RAPT::rsRationalFunction<T>;
  if(isBilinear())             // zero at z = -1
  {
    // Compute some intermediate variables:
    T a2   = a*a;          // a^2
    T a3   = a2*a;         // a^3
    T a4   = a2*a2;        // a^4
    T A    = a+1;          // a+1 == 2*b
    T A2   = A*A;          // A^2
    T A3   = A2*A;         // A^3
    T A4   = A2*A2;        // A^4
    T A4k  = A4*k;
    T C0   = 16*c[0];      // from here, they are relevant only for numerator
    T C1   =  8*c[1]*A;
    T C2   =  4*c[2]*A2;
    T C3   =  2*c[3]*A3;
    T C4   =    c[4]*A4;
    T ap3  = a   + 3;
    T a3p1 = a*3 + 1;

    // Create denominator:
    std::vector<T> D(6);
    D[0] =              A4k;
    D[1] = 4*( 4*a4 +   A4k);
    D[2] =    64*a3 + 6*A4k;
    D[3] = 4*(24*a2 +   A4k);
    D[4] =    64*a  +   A4k;
    D[5] =              16;

    // Create numerator:
    std::vector<T> N(6);
    N[0] =   T(0);
    N[1] =   C0*a4 + C1*a3     + C2*a2         + C3*a    + C4;
    N[2] = 4*C0*a3 + C1*ap3*a2 + C2*A*2*a      + C3*a3p1 + C4*4;
    N[3] = 6*C0*a2 + C1*A*3*a  + C2*(a2+4*a+1) + C3*A*3  + C4*6;
    N[4] = 4*C0*a  + C1*a3p1   + C2*A*2        + C3*ap3  + C4*4;
    N[5] =   C0    + C1        + C2            + C3      + C4;
    if(withGain) 
      rsScale(N, g);

    // Create and return rational function object:
    return RF(N, D);
  }
  else if(B1 == T(0))  // no zero
  {
    // Compute some intermediate variables:
    T b  = T(1)+a;
    T b2 = b*b;        // b^2
    T b4 = b2*b2;      // b^4
    T a2 = a*a;        // a^2
    T d0 = c[0];       // c0 * b^0  rename to C0 as above
    T d1 = c[1]*b;     // c1 * b^1
    T d2 = c[2]*b2;    // c2 * b^2
    T d3 = c[3]*b2*b;  // c3 * b^3
    T d4 = c[4]*b4;    // c4 * b^4

    // Compute coefficients of the numerator (note the binomial coeffs in the columns:
    // (1),(1,1),(1,2,1),(1,3,3,1),(1,4,6,4,1)):
    std::vector<T> N(5);
    N[4] = (d4 + d3 +      d2 +      d1 +      d0);
    N[3] = (     d3 + T(2)*d2 + T(3)*d1 + T(4)*d0) * a;
    N[2] = (               d2 + T(3)*d1 + T(6)*d0) * a2;
    N[1] = (                         d1 + T(4)*d0) * a2*a;
    N[0] = (                                   d0) * a2*a2;
    if(withGain) 
      rsScale(N, g);

    // Create rational function object and return it (note again the binomial coeffs 1,4,6,4,1) 
    // with the additional b4*k added in, so the denominator seems to be (a+1)^4 + b4*k*z^3
    return RF(N, { a2*a2, T(4)*a*a2, T(6)*a2, T(4)*a + b4*k, T(1) });
  }
  else  // the general case
  {
    // Compute some intermediate variables:
    T u    = B1;
    T u2   = u*u;
    T u3   = u2*u;
    T u4   = u2*u2;
    T a2   = a*a;          // a^2
    T a3   = a2*a;         // a^3
    T a4   = a2*a2;        // a^4
    T A    = a+1;
    T A2   = A*A;          // A^2
    T A3   = A2*A;         // A^3
    T A4   = A2*A2;        // A^4
    T A4k  = A4*k;
    T U1   = u-1;
    T U2   = U1*U1;        // (u-1)^2
    T U3   = U2*U1;        // (u-1)^3
    T U4   = U2*U2;        // (u-1)^4

    // Create denominator:
    std::vector<T> D(6);
    D[0] =     A4k*u4;
    D[1] = -  (A4k*(4*u4 - 4*u3)              -   a4);
    D[2] =  2*(A4k*(3*u4 - 6*u3 + 3*u2)       + 2*a3);
    D[3] = -2*(A4k*(2*u4 - 6*u3 + 6*u2 - 2*u) - 3*a2);
    D[4] =    (A4k*U4                         + 4*a );
    D[5] =  1;

    // Create numerator:
    std::vector<T> N(6);
    N[0] = T(0);

    // G0 contributions:
    N[1] += c[0]*a4;
    N[2] += c[0]*4*a3;
    N[3] += c[0]*6*a2;
    N[4] += c[0]*4*a;
    N[5] += c[0]*1;

    // G1 contributions:
    N[1] += c[1]*(   A *      u            * a3);
    N[2] += c[1]*(  -A * (  a*u - a - 3*u) * a2);
    N[3] += c[1]*(-3*A * (  a*u - a   - u) * a );
    N[4] += c[1]*(  -A * (3*a*u - 3*a - u)     );
    N[5] += c[1]*(  -A * (u - 1)               );

    // G2 contributions:
    N[1] += c[2]*(   A2 *  a2*u2);
    N[2] += c[2]*(-2*A2 * (a*u - a - u)*a*u);
    N[3] += c[2]*(   A2 * (a2*u2 - 2*a2*u - 4*a*u2 + a2 + 4*a*u + u2));
    N[4] += c[2]*( 2*A2 * (a*u - a - u)*U1 );
    N[5] += c[2]*(   A2 * U2);

    // G3 contributions:
    N[1] += c[3]*(  A3 *            a             * u3);
    N[2] += c[3]*( -A3 * (3*a*u - 3*a -   u)      * u2);
    N[3] += c[3]*(3*A3 * (  a*u -   a -   u) * U1 * u );
    N[4] += c[3]*( -A3 * (  a*u -   a - 3*u) * U2     );
    N[5] += c[3]*( -A3                       * U3     );

    // G4 contributions:
    N[1] +=    c[4]*A4*u4;
    N[2] += -4*c[4]*A4*u3*U1;
    N[3] +=  6*c[4]*A4*u2*U2;
    N[4] += -4*c[4]*A4*u*U3;
    N[5] +=    c[4]*A4*U4;

    if(withGain) 
      rsScale(N, g);
    return RF(N, D);
  }
}
// The expressions were obtained with sage via:
//
//   var("a z k u")              # u is shorthand for B1
//   b1 = (a+1)*u
//   b0 = (a+1)*(1-u)
//   G(z) = (b0+b1/z)/(1+a/z)
//   H    = G^4 / (1 + k*G^4/z)  # replace numerator by G^0,G^1,G^2,G^3 to select the stage
//   Hs = H.simplify_rational()
//   #Hs
//   Hs.factor()
//   #maxima.optimize(Hs)        # It's still a mess! -> optimize manually!
//
// and manual simplification of the results, also with the help of sage. Especially, 
// .factor() and .expand() turned out to be helpful. Also, The binomial theorem could be used a 
// lot. For each stage, such an expression was derived. They have all the same denominator. The 
// different numerators must be all weighted by the respective c-coefficient and accumulated.
//   Hs.factor()

// ToDo:
// Maybe make a function getCoeffsBA(TPar* b, TPar* a) which just fills the arrays (of length 5).
// Beware that the order must be reversed for a polynomial in z^-1 (maybe make it optional to be
// in this form)
// Check, if the patterns with the binomial coeffs extend to an N-stage ladder. If so, write down
// the general form of the transfer function in canonical form and implement a function that 
// computes numerator and denominator coeffs for the general case.

template<class TSig, class TPar>
rsRationalFunction<TPar> rsLadderTest<TSig, TPar>::getTransferFunctionOld()
{
  using T  = TPar;
  using RF = RAPT::rsRationalFunction<T>;
  T tol = 1024 * RS_EPS(T);  // ad hoc

  T B0 = T(1) - B1;
  T b0 = B0 * (1+a);
  T b1 = B1 * (1+a);

  RF G1( { b1, b0 }, { a, 1 }, tol); // G1(z) = (b0 + b1/z) / (1 + a/z) = (b1 + b0*z) / (a + 1*z)
  RF one({ 1      }, { 1    }, tol); // 1 = 1 / 1
  RF z  ({ 0,  1  }, { 1    }, tol); // z = (0 + 1*z) / 1
  RF G2 = G1*G1;                     // G1^2
  RF G3 = G1*G2;                     // G1^3
  RF G4 = G2*G2;                     // G1^4
  RF H  = g * (c[0]*one + c[1]*G1 + c[2]*G2 + c[3]*G3 + c[4]*G4) / (one + k * G4 / z); // H(z)
  return H;
}

template class rsLadderTest<double, double>; // explicit instantiation


//=================================================================================================

rsGroupString rsGroupString::inverse() const
{
  size_t len = s.size();
  std::vector<unsigned int> r(len);
  for(size_t i = 0; i < len; i++)
    r[i] = s[len-1-i];
  return r;
}

rsGroupString rsGroupString::operator+(const rsGroupString &rhs) const
{
  std::vector<unsigned int> r = s;
  for(size_t i = 0; i < rhs.s.size(); i++) {
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
  for(size_t i = 0; i < s.size(); i++)
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