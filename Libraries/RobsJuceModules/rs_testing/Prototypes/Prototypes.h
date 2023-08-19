#ifndef RS_PROTOTYPES_H
#define RS_PROTOTYPES_H

// todo: move all prototyes into rs_testing juce module - and maybe much of the other code, too

//#include "rapt/rapt.h"
#include "rosic/rosic.h"
using namespace RAPT;

// new implementation of classic IIR filter design:
#include "FilterDesign/PoleZeroPrototype.h"
#include "FilterDesign/PoleZeroMapper.h"
#include "FilterDesign/PoleZeroDesignerAnalog.h"
#include "FilterDesign/PoleZeroDesignerDigital.h"
#include "FilterDesign/ComplementaryFilters.h"

#include "SafeArray.h"
#include "ModalAnalyzer.h"
#include "ParticleBouncer.h"
#include "Probability.h"
#include "Projection3Dto2D.h"
#include "Polygon.h"
#include "PixelClassifier.h"
#include "Drawing.h"
#include "QuantumSystems.h"
#include "Relativity.h"
#include "SimdVector.h"
#include "SineParameterEstimator.h"
#include "OdeSolver.h"
#include "MeshStuff.h"
#include "DampedSine.h"
#include "AdditiveSynthEngine.h"
#include "MultiplicativeSynth.h"
#include "RenderScriptTools.h"
#include "FractalRenderer.h"

/** This file contains prototypical implementations of algorithms. These prototypes are not meant
to be used for production code but are useful for a more readable proof-of-concept (because of lack
of optimizations), for tweaking an algorithm's internal parameters which might not be even exposed
in the production-code versions, and to create reference output for the unit-tests for production
code. */

/*
moved to rapt:
static constexpr int allBits = -1;                                      // all bits are 1
static constexpr int allBitsButFirst = std::numeric_limits<int>::max(); // only 1st bit is 1
static constexpr int firstBitOnly = allBits ^ allBitsButFirst;          // only 1st bit is 0

// for unsiged int types, the bit twiddling is different:
//static size_t allBits = std::numeric_limits<size_t>::max();
//static size_t firstBitOnly = allBits - (allBits >> 1);
//static size_t allBitsButFirst= allBits ^ firstBitOnly;
*/

/** 
From pg 125 here: https://www.cs.otago.ac.nz/graphics/Geoff/tartini/papers/Philip_McLeod_PhD.pdf
@param x Input
@param y Output
@param n Size of Symmetric Hanning Window
@param len The size of x and y  */
template<class T>
void complexMovingAverage(const T* x, T* y, int N, int L);



/** Returns the index of the element in the array A (of length N) that is the best match to the 
given x according to some function "isCloser" that determines which of two given values is closer 
to a reference value. If there are multiple best matches in A with the same (minimal) distance to 
x, it will return the index of either the first or last of them, depending on whether the isCloser 
function uses a < or <= comparison for determining closeness (use <, if you want the first). The 
isCloser function should take 3 arguments a,b,r of type T and return true, iff a is closer to r 
than b (r is the reference value). */
template<class T, class F>
int findBestMatch(T* A, int N, const T& x, const F& isCloser)
{
  T   minVal = A[0];
  int minIdx = 0;
  for(int i = 0; i < N; i++)  // maybe we cant start at i = 1
  {
    if(isCloser(A[i], minVal, x)) 
    {
      minVal = A[i];
      minIdx = i;
    }
  }
  return minIdx;
}
// todo: 
// -move into rsArrayTools ...but maybe re-order the parameters of isCloser to a,r,b, maybe rename
//  isCloser to isBetterMatch
// -write unit test

/** Converts a Taylor approximation with coeffs given in t into a Pade approximation with numerator 
coeffs in p and denominator coeffs in q. It is assumed that p and q already have the correct sizes. 
The degrees of the polynomials p and q must add up to the degree of t: deg(p) + deg(q) = deg(t). 
The degree of a polynomial represented by its coefficient array is always size() - 1. */
template<class T>
void rsTaylorToPade(const std::vector<T>& t, std::vector<T>& p, std::vector<T>& q);

/** Under construction.
Solves the linear system A*x = b in an "optimal" way regardless of the question whether the system 
is over-, under- or critically determined. If it's critically determined and the matrix is regular, 
this just boils down to finding the unique exact solution. In the overdetermined case, it finds a 
least-squares approximation to a solution. In the underdetermined case, it finds the minimum-norm
solution among the infinitely many. 
...but what if a critically determined system is singular? ...tbc... */
template<class T>
void solveOptimal(rsMatrix<T>& A, rsMatrix<T>& X, rsMatrix<T>& B);

/** Given abscissa and ordinate values in x and y, this function computes the slopes at the 
x-locations, i.e. values for the 1st derivative of y, and returns them as vector which is of the 
same length as x and y. If these slopes are later used in cubic Hermite interpolation, the Hermite 
interpolant will actually become a spline interpolant, i.e. an interpolant that is 2nd order smooth
at the junctions. If you just use numerical differentiation to obtain such slope values, 2nd order
smoothness will in general not happen. Conversely, the computed slopes themselves can also be seen
as another way of obtaining numerical estimates of the 1st derivative. If prescribe2ndDeriv is 
true, the final 2 parameters prescribe values for the 2nd derivative at the endpoints (defaulting 
to zero, leading to what is called a "natural" cubic spline). If it is false, these values will be 
interpreted as prescriptions for the 1st derivative at the endpoints (leading to what is called a 
"complete" cubic spline. ..tbc... 

this needs unit-tests and experiments
*/
template<class T>
std::vector<T> splineSlopes(const std::vector<T>& x, const std::vector<T>& y, 
  bool prescribe2ndDeriv = true, T ypStart = T(0), T ypEnd = T(0));
// maybe the production version of this should go into rsNumericalDifferentiator
// todo: also implement periodic boundary conditions

// some experimental sin/cos approximations:
double rsCos2(double x);
void rsSinCos1(double x, double* s, double* c);
void rsSinCos2(double x, double* s, double* c);
void rsSinCosApprox4(double x, double* s, double* c);


/** Computes the coefficients for an Adams-Bashforth multistep method for numerically solving an 
ordinary differential equation (ODE). Directly implements the formula given here:
  https://en.wikipedia.org/wiki/Linear_multistep_method#Adams%E2%80%93Bashforth_methods
without any attempt to optimize the efficiency. I recommend to use rsFraction<int> for the template
parameter T. Maybe the algo can be turned into an O(N) algo by not creating the polynomial p from 
scratch leaving out the i=j factor each time but instead construction a "master" polynomial and 
dividing out the i=j factor in each iteration. Oh, and the factorials could be computed more 
efficiently, too. */
template<class T>
std::vector<T> coeffsAdamsBashforth(int order)
{
  using Poly = rsPolynomial<T>;
  int s = order;
  std::vector<T> b(s);
  T sign = T(1);                        // for the (-1)^j
  for(int j = 0; j < s; j++) {
    Poly p({ T(1) });
    for(int i = 0; i < s; i++) 
      if(i != j)
        p = p * Poly({T(i), T(1)});
    T d = p.definiteIntegral(T(0), T(1));
    b[s-j-1] = (sign*d) / (rsFactorial(j) * rsFactorial(s-j-1));
    sign *= T(-1);  }
  return b;
}

/** Similar to coeffsAdamsBashforth. See:
  https://en.wikipedia.org/wiki/Linear_multistep_method#Adams%E2%80%93Moulton_methods  */
template<class T>
std::vector<T> coeffsAdamsMoulton(int order)
{
  using Poly = rsPolynomial<T>;
  int s = order;
  std::vector<T> b(s+1);
  T sign = T(1);
  for(int j = 0; j <= s; j++) {
    Poly p({ T(1) });
    for(int i = 0; i <= s; i++) 
      if(i != j)
        p = p * Poly({T(i-1), T(1)}); 
    T d = p.definiteIntegral(T(0), T(1));
    b[s-j] = (sign*d) / (rsFactorial(j) * rsFactorial(s-j));
    sign *= T(-1); }
  return b;
}
// These are now superseded by the static functions in class rsOdeCoeffs. But maybe we should keep 
// them as prototypes anyway and then optimize the functions in rsOdeCoeffs


template<class T>
T newton(const std::function<T(T)>& f, const std::function<T(T)>& fp, T x0, T yt = T(0))
{
  T tol = 1.e-10;  // preliminary
  T x   = x0;
  int its = 0;
  while(true)
  {
    T y = f(x);    // should that be y = f(x) - yt; ?
    if(rsAbs(y) <= tol)
      break;
    x -= y / fp(x);
    its++;
  }
  return x;
}
// todo: use yt for the target value, provide a maxIts parameter, move to rsRootFinder, maybe 
// detect divergence



template<class T>
void weightedSum(const T* x1, int N1, T w1, const T* x2, int N2, T w2, T* y, int Ny)
{
  int n = 0;

  while(n < rsMin(N1, N2, Ny)) {
    y[n] = w1 * x1[n] + w2 * x2[n];
    n++; }

  if(n == Ny)
    return;

  // handle overhanging part in x1 or x2:
  if(N1 > N2) {
    while(n < rsMin(N1, Ny)) {
      y[n] = w1 * x1[n];
      n++; }}
  else if(N2 > N1) {
    while(n < rsMin(N2, Ny)) {
      y[n] = w2 * x2[n];
      n++; }}

  // handle trailing zeros in result, if necessarry:
  while(n < Ny) {
    y[n] = T(0);
    n++; }
}
// todo: document, write test, move into rsArrayTools


void solveTriDiagGauss(const std::vector<double>& lowerDiag, std::vector<double>& mainDiag, 
  const std::vector<double>& upperDiag, std::vector<double>& x, std::vector<double>& b);

void solveTriDiagThomas(const std::vector<double>& lowerDiag, const std::vector<double>& mainDiag, 
  std::vector<double>& upperDiag, std::vector<double>& x, std::vector<double>& b);

void solveWrappedTriDiag(const std::vector<double>& lowerDiag, std::vector<double>& mainDiag, 
  const std::vector<double>& upperDiag, std::vector<double>& x, std::vector<double>& b);

/** Solves a pentadiagonal linear system of equations with given diagonals and right-hand side
using a simple algorithm without pivot-search. lowerDiag1 is the one directly below the main
diagonal, lowerDiag2 the one below lowerDiag1 - and similarly for upperDiag1/upperDiag2. In the
process of the computations, the right hand side vector is destroyed. the same is true for mainDiag
and the two inner sub/superdiagonals lowerDiag1, upperDiag1. Note also that you can't use the same
array for lowerDiag1 and upperDiag1, even if your matrix is symmetric.

..What about lowerDiag2/upperDiag2? are these preserved and may these point to the same vector?
It's probably safest to assume that everything may get messed up and all arrays should be
distinct. */
std::vector<double> solvePentaDiagonalSystem(
  std::vector<double>& lowerDiag2, std::vector<double>& lowerDiag1,
  std::vector<double>& mainDiag,
  std::vector<double>& upperDiag1, std::vector<double>& upperDiag2,
  std::vector<double>& righHandSide);

/** Multiplies a pentadiagonal matrix with a vector...  */
std::vector<double> pentaDiagMatVecMul(
  std::vector<double>& lowerDiag2, std::vector<double>& lowerDiag1,
  std::vector<double>& mainDiag,
  std::vector<double>& upperDiag1, std::vector<double>& upperDiag2,
  std::vector<double>& input);




/** Minimizes the sum of squared differences between adjacent array elements under the constraint
that the sums of adjacent array elements must be equal to given values. The input array s is the
length N-1 array of the desired sums, the output array v is the length N value array, such that
v[i] + v[i+1] = s[i] and sum_i (v[i+1] - v[i])^2 = min. You may optionally pass an array of
weights for the squared differences in the cost function - if you do, the w array must have the
same length as s, if you don't, unit weights will be used for each squared difference. With
weights, we will minimize sum_i w[i] * (v[i+1] - v[i])^2 subject to the (same) constraints that
v[i] + v[i+1] = s[i] for all i = 0,..,N-2 */
std::vector<double> rsMinSqrDifFixSum(
  const std::vector<double>& s,
  const std::vector<double>& w = std::vector<double>() );

std::vector<double> rsMinSqrCrvFixSum(
  const std::vector<double>& s,
  const std::vector<double>& w = std::vector<double>() );



/** Prototype for rsResampler::signalValueViaSincAt(). It provides as additional parameters for
tweaking:
-pointer to a window-function
-parameter for the window (if applicable)
-switch for normalizing the output by the sum of the tap weights
*/
double signalValueViaSincAt(double *x, int N, double t, double sincLength, double stretch,
  //FunctionPointer3DoublesToDouble windowFunction = rsExactBlackmanWindow,
  double (*windowFunction)(double,double,double) = &RAPT::rsWindowFunction::exactBlackman,
  double windowParameter = 0.0, bool normalizeDC = true);

/** Generates polynomial coefficients of the polynomial used in Halpern filters. It's the T^2(w)
polynomial in Eq. 8.18 in Paarmann: Design and Analysis of Analog Filters. */
void halpernT2(double *c, int N);

/** Generates polynomial coefficients of the polynomial used in Papoulis filters. It's the L^2(w)
polynomial in Eq. 8.14 in Paarmann: Design and Analysis of Analog Filters */
void papoulisL2(double *c, int N);

/** Retrieves damped sine filter design parameters from its coefficients. See
@see rsDampedSineFilter for meaning of parameters. The phase p is returned in the interval
0...2pi. */
void rsDampedSineFilterAnalysis(double b0, double b1, double a1, double a2, double* w, double* A,
  double* d, double* p);

void rsDampedSineFilterAnalysis2(double b0, double b1, double a1, double a2, double* w, double* A,
  double* d, double* p); // other algorithm for the same thing

void rsDampedSineFilterResidueAndPole(double b0, double b1, double a1, double a2,
  std::complex<double>* residue, std::complex<double>* pole);


/** Calculates a chebyshev window of size N, store coeffs in out as in Antoniou
  -out should be array of size N
  -atten is the required sidelobe attenuation (e.g. if you want -60dB atten, use '60')
Dolph-Chebychev window generation code from here:
http://practicalcryptography.com/miscellaneous/machine-learning/implementing-dolph-chebyshev-window/
not recommended for production use because the complexity is O(N^2) - instead use an iFFT
approach
References:
[1] Lyons, R., "Understanding Digital Signal Processing", Prentice Hall, 2004.
[2] Antoniou, A., "Digital Filters", McGraw-Hill, 2000.  */
void cheby_win(double *out, int N, double atten);


void rsCircularShift(int* a, int N, int k);
// circular shift without additional memory (using 3 reversals) - needs test

//=================================================================================================

/** Linear fractional interpolation is an interpolation method constructed around the so called 
linear fractional transformations. These are functions of the form:

          a*x + b
  f(x) = ---------
          c*x + d

It is suitable only for strictly monotonic data. The interpolant that interpolates a segment 
between two data points will also be monotonic. Moreover, the interpolant will be easily 
invertible and the inverse interpolating function will be of the same kind and can easily be 
obtained from the forward interpolating function. The motivation for such a scheme lies in the 
desire to have an interpolation method that is both invertible and first order smooth, i.e. has 
matching derivatives at the nodes. Linear interpolation is easily invertible but not smooth. with
polynomial interpolation, it's the other way around.

...TBC...  

*/

template<class T>
class rsLinearFractionalInterpolator
{
  // Maybe have two template parameters for x and y. But maybe that would count as speculative
  // generality code smell?

public:

  /** Given two length N arrays x, y with x-axis values and corresponding y-axis values, this
  function fills the array yi with values corresponding to the xi by linear fractional ("linfrac")
  interpolation (or extrapolation, if necessary) of the input data. The xi and yi arrays are of 
  length Ni. With the extrapolateLinearly switch, you can choose whether or not the linfrac 
  method should also be used for extrapolation or if for extrapolation, we should 
  switch to linear. Extrapolating linearly is the default because it is the safe way. Linfrac
  extrapolation may shoot off to infinity due to a potentially present pole of the linfrac function
  in the extrapolated region so it should be used with great caution, if ever.  */
  static void interpolate(const T* x, const T* y, T* s, int N, const T* xi, T* yi, int Ni, 
    bool extrapolateLinearly = true, T shape = 0.5);
  // ToDo: 
  // -Document the shape parameter or maybe remove it or allow the user to pass a whole array
  //  of shape parameters to set it separately for each segment
  // -Rescale the shape parameter such that the neutral value is 0 instead of 0.5

  /** Implements the basic linear fractional transformation f(x) on which everything else is based.
  The function maps the unit interval [0,1] to itself with f(0) = 0 and f(1) = 1 and with given 
  slope s at the origin such that f'(0) = s. It produces a concave shape (like a saturation curve)
  for s > 1 and a convex shape (like one side of a bowl or parabola) for s < 1. For s = 1, it 
  produces the identity function. The slope at x,y = 1,1 will be given by 1/s. The slope s must be 
  in the interval (0,inf) excluding the boundaries. */
  static T simpleMap(T value, T slopeAt0);
  // It's actually the same function as implemented in rsBiRationalMap_01 in RealFunctions.h just 
  // parametrized differently. There, we use a parameter p in the range (-1,+1) and here we use the
  // slope at the origin in the range (0,inf) as patameter. The conversion formulas are:
  //   s = (1+p)/(1-p), p = (s-1)/(s+1)
  // Maybe rename that function there to rsLinFracMap


  /** Implements a symmetrized version of the (simple, basic, prototypical) linear fractional map. 
  It uses two appropriately scaled and shifted versions of the original map in the interval 0..0.5
  and 0.5..1 to produce a shape that is symmetric around x,y = 0.5,0.5. The shape looks like a 
  sigmoid (s-shaped) for s < 1 and like a saddle (like x^3) for s > 1. */
  static T symmetricMap(T value, T slopeAt0);

  /** Implements a function composed of three maps with slopes s1,s2,s3 where the middle one is the
  symmetrized variant. The end result is a function that sadwiches a sigmoid/saddle shape between 
  two convex/concave shapes. This is the final shape that we want. */
  static T tripleMap(T value, T slope1, T slope2, T slope3);

  /** Computes the partial slopes s1,s2,s3 for the three individual maps in tripleMap such that the
  composed/sandwiched map produces the desired slope targetSlopeAt0 at the origin x,y = 0,0 and the 
  desired slope targetSlopeAt1 at the end point x,y = 1,1 of the unit interval. */
  static void computeSlopes(T targetSlopeAt0, T targetSlopeAt1, 
    T* mapSlope1, T* mapSlope2, T* mapSlope3, T shape = T(0.5));
  // ToDo: 
  // -Rename to calcSlopes for consistency
  // -Document the shape parameter. And/or maybe get rid of it by commenting it out, so we can 
  //  easily add it back later

  /** Given the slope s1 of the first map at zero, this function computes the split point at which
  we should switch between the two partial (composed) maps. */
  static T getSplitPoint(T s1) { return simpleMap(T(0.5), T(1)/s1); }
  // We need to apply the first map's inverse to 0.5. We want to figure out at which input x the
  // first map produces the output 0.5. The inverse map is obtained by using the reciprocal slope.
  // We need 0.5 because that's the value at which the 2nd map switches between its two halves.


  static void calcComposedCoeffsLeft(T s1, T s2, T s3, T* a, T* b, T* c, T* d);

  static void calcComposedCoeffsRight(T s1, T s2, T s3, T* a, T* b, T* c, T* d);

  // Make a function 
  //static void calcSplitPointAndCoeffs(T targetSlopeAt0, T targetSlopeAt1, 
  //  T* splitPoint, T* aL, T* bL, T* cL, T* dL, T* aR, T* bR, T* cR, T* dR);


  /** Computes the a,b,c,d parameters for the combined triple map, i.e. the map which results from
  applying the 3 maps one after another. That composed map is not the same for all input values. 
  That's why you need to provide the input value in the first parameter, too. This is because the
  inner, symmetrized map contains a switch between two maps at 0.5. This switch translates to a 
  switch somewhere in the composed map (whereever the first map maps to 0.5). */
  static void composeTripleMapIntoOne(T value, T slope1, T slope2, T slope3, 
    T* a, T* b, T* c, T* d);

  /** Takes a normalized input value from the unit interval [0..1] and produces the cooresponding 
  normalized output value also in [0..1]. You need to specify the dsired normalized slopes at the
  origin xy = 0,0 and the end point x,y = 1,1. */
  static T getNormalizedY(T normalizedX, T slopeAt0, T slopeAt1, T shape);

  // ToDo:
  // -Order the functions here in a more top-down rather than bottom-up fashion because the 
  //  higher-level functions are more likely to be important for the user. For example,
  //  getNormalizedY is currently the most high level function - it should appear at or near the
  //  top. They are ordered bottom-up because that's the order in which they have been written.

};

template<class T>
T rsLinearFractionalInterpolator<T>::simpleMap(T x, T s)
{
  return s*x / ((s-1)*x + 1);
}

template<class T>
T rsLinearFractionalInterpolator<T>::symmetricMap(T x, T s)
{
  if(x < T(0.5))
    return T(0.5) * simpleMap(T(2)*x, s);
  else
    return T(0.5) * (simpleMap((T(2)*x-T(1)), T(1)/s) + T(1));
}

template<class T>
T rsLinearFractionalInterpolator<T>::tripleMap(T x, T s1, T s2, T s3)
{
  x = simpleMap(   x, s1);  // pre convexity or concavity
  x = symmetricMap(x, s2);  // sigmoidity or saddleness
  x = simpleMap(   x, s3);  // post convexity or concavity
  return x;
}

template<class T>
void rsLinearFractionalInterpolator<T>::computeSlopes(
  T slopeAt0, T slopeAt1, T* s1, T* s2, T* s3, T shape = T(0.5))
{
  // Compute slope at zero for middle map (controlling sigmoidity vs saddleness) and slope at 
  // zero for combined outer maps s13 = s1*s3 (controlling convexity vs concavity):
  *s2 = sqrt(slopeAt0 * slopeAt1);
  T s13 = slopeAt0 / *s2;            // == *s2 / slopeAt1

  // Compute slopes at zero for first and last map according to our shape parameter:
  if(shape == 0.5)
    *s1 = *s3 = sqrt(s13);           // Optimized special case for symmetric shape
  else {
    *s1 = pow(s13, shape);           // General case for user controlled shape
    *s3 = s13 / *s1; }

  // Question:
  // Could there be any numercial reasons for prefering s13 = slopeAt0 / *s2 vs using
  // s13 = *s2 / slopeAt1 or vice versa? Maybe depending on which of the two slopes is greater?
  // Maybe try it with some crappy fixed point format for T.
}

template<class T>
void rsLinearFractionalInterpolator<T>::calcComposedCoeffsLeft(
  T s1, T s2, T s3, T* a, T* b, T* c, T* d)
{
  *a = s1*s2*s3;
  *b = 0;
  *c = s1*s2*s3 + s1*s2 - s1 - 1;
  *d = 1; 
}

template<class T>
void rsLinearFractionalInterpolator<T>::calcComposedCoeffsRight(
  T s1, T s2, T s3, T* a, T* b, T* c, T* d)
{
  *a = s1*s3 - s2*s3 + s3;
  *b = s2*s3 - s3;
  *c = (s1*s3 - s2*s3 - s2 + s3);
  *d = s2*s3 + s2 - s3; 
}

template<class T>
void rsLinearFractionalInterpolator<T>::composeTripleMapIntoOne(T x, T s1, T s2, T s3,
  T* a, T* b, T* c, T* d)
{
  // Figure out the split point where we need to switch between the two maps. 
  T xs = getSplitPoint(s1);

  // Compute coeffs for the composed map depending on whether the input value x is before or after
  // the split point:
  if(x <= xs)
    calcComposedCoeffsLeft(s1, s2, s3, a, b, c, d);
  else
    calcComposedCoeffsRight(s1, s2, s3, a, b, c, d);
}

template<class T>
T rsLinearFractionalInterpolator<T>::getNormalizedY(T x, T slopeAt0, T slopeAt1, T shape)
{
  // Compute the slopes for the 3 individual maps that compose to the final overall map:
  T s1, s2, s3;
  computeSlopes(slopeAt0, slopeAt1, &s1, &s2, &s3, shape);

  // Compute the coefficients for the combined map which is of the more general form
  // f(x) = (a*x + b) / (c*x + d). 
  T a, b, c, d;
  composeTripleMapIntoOne(x, s1, s2, s3, &a, &b, &c, &d);
  // This step is where the switch between the two linear fractional maps happens. It is the 
  // switch between two maps that allows us to prescribe 2 derivatives and have actually even a 
  // degree of freedom left to control the shape.

  // Compute the output of the composed map:
  T y = (a*x + b) / (c*x + d); 
  //T error = y - tripleMap(x, s1, s2, s3);  // Compare to reference computation for debug
  return y;

  // I'm not sure, if it's worth it for interpolating a single data point to produce the coeffs 
  // for the composed map vs just literally applying all 3 maps in sequence. Obtaining the 
  // composed map makes more sense when the map is used on many points within a segment such as 
  // in upsampling an array of data. Then, we would for each segment compute the coeffs just 
  // twice (once for the sub-segment below the split-point and once for the sub-segment above it).
}

template<class T>
void rsLinearFractionalInterpolator<T>::interpolate(
  const T* x, const T* y, T* s, int N, const T* xi, T* yi, int Ni, 
  bool extrapolateLinearly, T shape)
{  
  // The code below follows closely rsInterpolateLinear.

  int n = 0;       // Index into input data
  int i = 0;       // Index into interpolated data
  T dx, dy, dxr;   // Segment length (or width), height and reciprocal of length.
  T a, b, c, d;    // Coeffs for the linfrac y = (a*x + b) / (c*x + d)

  // Linear front extrapolation, if necessary:
  if(extrapolateLinearly) {
    while(xi[i] < x[0]) {
      yi[i] = y[0] + s[0] * (xi[i] - x[0]);
      i++; }}
  // It is necessary, if the new x-values for the output data start before the first x-value of 
  // the input data and the caller has requested linear extrapolation. Unlike for the tail 
  // extrapolation below, we need no else-branch here because if we don't want to extrapolate 
  // linearly, a linear fractional extrapolation for the front will naturally occur in the main 
  // loop.

  // Main Loop over the input datapoints:
  while(n < N-1)
  {
    dx  = x[n+1] - x[n];                // Length of current segment
    dy  = y[n+1] - y[n];                // Height of current segment
    dxr = 1 / dx;                       // Reciprocal of the length dx

    // Retrieve and normalize the slopes:
    T slopeScale = dx/dy;
    T slopeAt0   = slopeScale * s[n];
    T slopeAt1   = slopeScale * s[n+1];

    // Do a sanity check:
    rsAssert(slopeAt0 > 0 && slopeAt1 > 0, "Data is not strictly monotonic");
    // If this happens, it means that you have passed in data that is not strictly monotonically
    // increasing or decreasing and/or one of the dervative values is not consistent with the
    // data. In such a case, this interpolation scheme will fail badly to the point of producing
    // NaNs.
    // ToDo: 
    // -Figure out, if the > 0 test is actually strict enough. Maybe extremely small nonzero values
    //  are also already problematic and it should be > eps or something? Maybe then do something 
    //  like slopeAt0 = max(slopAt0, eps); slopeAt1 = max(slopeAt1, eps); to avoid producing NaNs. 
    //  ...but the "eps" is just an ad hoc guess. Figure out, what value is suitable.

    // Compute values for the initial slopes (at x,y = 0,0) for the 3 normalized linear 
    // fractional maps:
    T s1, s2, s3;
    computeSlopes(slopeAt0, slopeAt1, &s1, &s2, &s3, shape);

    // Compute the split point, i.e. the point at which we have to switch between the two maps
    // for the two sub-segments:
    T xSplit = getSplitPoint(s1);       // Normalized split point in [0,1]
    xSplit = xSplit * dx + x[n];        // Denormalize for comparisons with the actual xi[i]

    // Calculate the a, b, c, d coeffs for the linear fractional map y = (a*x + b) / (c*x + d)
    // for the left sub-segment and interpolate it:
    calcComposedCoeffsLeft(s1, s2, s3, &a, &b, &c, &d);
    while(xi[i] <= xSplit && i < Ni) {
      T xn  = dxr * (xi[i] - x[n]);     // Normalized x in [0,1]
      T yn  = a*xn / (c*xn + d);        // Normalized y in [0,1], b == 0  -> was optimized out
      yi[i] = y[n] + dy*yn;             // Denormalize y and write to output
      i++; }

    // Calculate the coeffs for the right sub-segment and interpolate it:
    calcComposedCoeffsRight(s1, s2, s3, &a, &b, &c, &d);
    while(xi[i] < x[n+1] && i < Ni) {
      T xn = dxr * (xi[i] - x[n]);
      T yn = (a*xn + b) / (c*xn + d);   // Here, b != 0  -> use the full formula
      yi[i]   = y[n] + dy*yn; 
      i++; }

    n++;
  }

  // Tail extrapolation, if necessarry, i.e. if the new x-values for the output data go beyond
  // those of the input data:
  if(extrapolateLinearly) {
    while(i < Ni) {
      yi[i] = y[N-1] + s[N-1] * (xi[i] - x[N-1]);
      i++; }}
  else {
    while(i < Ni) {
      T xn = dxr * (xi[i] - x[N-2]);
      T yn = (a*xn + b) / (c*xn + d);  // Use the last computed a,b,c,d coeffs
      yi[i]   = y[N-2] + dy*yn;
      i++; }}
}
// ToDo:
// -Maybe include some sanity checks for the assumed montonicity of the data. We should have
//  either (y[n+1] > y[n] and s[n] > 0) or (y[n+1] < y[n] and s[n] < 0)  for all n. Otherwise,
//  the method will fail. Maybe in such a case, fall back to linear interpolation or just output 
//  all zeros. The rationale to not fall back to linear is that when giving a clearly invalid 
//  result in a production environment, the bug in the higher level code that gave rise to the
//  invalid input data will be detected rather than just potentially go unnoticed and cause a 
//  quality degradation of some higher level algorithm. In a test/debug environment, we may 
//  trigger an assert.
// -Maybe write some unit tests that cover extreme situations, like N = 0,1,2,3, nonomonotonic
//  data, target slopes of zero, etc.




/*
Notes:

Linear interpolation is easily invertible but it's not first order smooth. Polynomial 
interpolation schemes can be made even smoother than first order but are in general not 
invertible - at least not easily - and even if it is invertible (because one has somehow 
restricted the interpolant to be monotonic), the inverse interpolating function will be of a 
completely different kind, i.e. come from a different class of functions and has to be computed 
with a different (more complicated) algorithm. For example, if you interpolate data given in 
arrays x[n], y[n] using piecewise cubic polynomials to produce a continuous function y = f(x), 
the (exact) inversion of the produced interpolant y = f(x), will not be a piecewise cubic but 
instead be defined by using a root-finder applied to the piecewise cubic y = f(x), so the 
inverse interpolant will be something like a piecewise "cube-rooty" function (well, it's 
actually even more complicated than that but you get the point).

Linear fractional interpolation solves this problem by using as interpolants between the 
segments functions of the general form y = f(x) = (a*x + b) / (c*x + d). These functions
are also known as "linear fractional transformations" because they consist of a fraction or
quotient of two linear functions. One nice feature of the linear fractional transformations
(which we will abbreviate as "linfrac maps" or just "linfracs" in the following) is that they
form a group, implying that the inverses are also in the set and compositions of such 
functions yield another element from the set. One might hope that with the 4 tweakable 
parameters, one could satisfy 4 constraints to let the values and derivatives match some 
prescribed value at the ends of the intervals. However, unfortunately, that's not so simple.
When trying to solve the resulting system of equations, it turns out that it can only be 
solved, if the two slopes at the ends of the interval are reciprocals of one another so we
can't prescribe them independently. I think, it is because the linfrac has actually only 3 
degrees of freedom rather than 4 because scaling all coeffs by the same number gives the same
transformation. Maybe that's why fixing the endpoints and one derivative already locks 
everything in place. To allow ourselves more control, we will use two linfracs per segment 
that join together at an internal node somewhere within the segment. At this internal node, 
the two sub-segments will also join smoothly to first order such that the overall smoothness
of the interpolant is unchanged.  

The construction of the interpolant for a segment works as follows. We will assume that we
want to find a function f(x) that maps the unit interval monotonically and invertibly to 
itself. That means the function should produce f(0) = 0 and f(1) = 1. Additionally, it should
produce f'(0) = d0 and f'(1) = d1 for some pair of prescribed slope values at the interval 
boundaries. We use d0, d1 for "derivative" rather than s0, s1 for "slope" because we will use
s1 later for something else These slope values will ensure that the overall interpolant will 
have matching slopes at the segment boundaries. Where we get these d0, d1 values from doesn't 
matter here. It may make sense to produce them by numerical differentiation of the actual 
data, but that's a different topic. Here, we just assume them to be given. 

In addition to the regular linFrac, we will also make use of a symmetrized variant that 
consists of two appropriately scaled and shifted versions of the orginal map, one for the
domain 0..0.5 and another for the domain 0.5..1. They join smoothly at 0.5,0.5 and implement 
an inflection point there. The overall shape in 0..1 will be that of a sigmoid or a saddle.
The idea is now to combine a regular concave or convex linFrac with such a sigmoid/saddle and 
then with another convex/concave one. So all in all we use 3 linFracs applied in sequence:
to obtain a neesting of 3 maps:
  concave/convex  ->  sigmoid/saddle  ->  concave/convex
These 3 linFracs in succession will give use (more than) enough freedom to meet our 
requirements for the derivatives at 0 and 1 at the expense of having to introduce a switch.
Because at the end of the day, all 3 will combine into one single linFrac due to the 
composition rule - but it will have to be a different one for the first and last section of 
the interval 0..1. It's the fact that middle one is the symmetrized variant that gives us this
additional amount of control. The symmetrized variant lies not in our group of linFracs so it
gives us indeed soemthing new - namely, the ability to introduce inflection points. The end
result of this construction will be a piecewise linfrac made from two pieces. In a way, the 
symmetrized linfrac serves as a prototype for the more general piecewise linfracs.

OK - so we start from the assumption that we first use a regular linFrac in sequence with a
symmetrized one and another regular one. A symmetrized linFrac will have the same derivative 
as a regular one at the origin. the regular ones will have a reciprocated derivative at 1,1 
while the symmetrized one has the same derivative at 1,1 (after all, it's symmetric around 
0.5,0.5). So all in all, at the origin, the overall slope will be given by the product of 
the slopes of the 3 individual linFracs, which we will denote by s1,s2,s3. At 1,1 the overall 
slope will be given by s2/(s1*s3). The goal is now to produce the 3 slope parameters for our
3 linFracs from the 2 slope requirements d0, d1 that are given. So we have to solve the 
following system of equations
  d0 = s1*s2*s3
  d1 = s2/(s1*s3)
Let's define s13 = s1*s3 to get a system in 2 variables:
  d0 = s13*s2
  d1 = s2/(s13)
If we multiply both equations, we get
  d0*d1 = s2^2
so we should use 
  s2 = sqrt(d0*d1)
Now, we can use either of the two equations to compute s13 as:
  s13 = d0/s2 = s2/d1
What's left is to distribute the combined slope s13 to s1 and s3 (remember that s13 = s1*s3). 
The most natural way to do this is to just set s1 = s3 = sqrt(s13) but we could also introduce
a shape parameter in 0..1 and set s1 = s13^shape, s3 = s13/s1. For shape = 0.5, this reduces
to the sqrt formula. Im not sure, if such a shape parameter is useful, though. The shapes 
obtained from the sqrt formula look most symmetric and seem to make most sense. I think, for 
later inversion, the shape needs either to be 0.5 or we need to use 1 - shape for the inverse 
interpolation (I've not yet figured out the details of that).

Now that we have our 3 slopes s1,s2,s3, we can combine the three linFracs into one for each
section of the unit interval. Because the second, middle, symmetrized linFrac splits at 0.5, 
we need to figure out where the first map produces 0.5. This is easily done by evaluating
an inverse linFrac, i.e. one with slope 1/s1, at 0.5. This gives us our split point xs. The 
first combined map will be used for x = 0..xs and the second combined map will be used for
x = xs..1. Now we need to produce the coefficients for the two maps. The first map will also
be of the restriced kind with b = 0 and d = 1, but for the second, we'll need the fully 
flexible variant with all 4 coeffs a,b,c,d. To figure out the coeffs of the maps, we'll just
have to apply them in succession algebraically to a variable x and then read off the coeffs
of the result. This is best done with a computer algebra system. The formulas for the coeffs 
for the combined maps were found with the following Sage code. We have two version because we 
need to switch between them:

 1st part:

 var("s1 s2 s3")
 y1 = s1*x / ((s1-1)*x + 1)     # 1st map
 y2 = 2*y1                      # scale 
 y3 = s2*y2 / ((s2-1)*y2 + 1)   # 2nd map 
 y4 = y3/2                      # scale
 y5 = s3*y4 / ((s3-1)*y4 + 1)   # 3rd map
 num = numerator(y5)
 den = denominator(y5)
 num.expand().collect(x), den.expand().collect(x)


 2nd part:

 var("s1 s2 s3")
 y1 = s1*x / ((s1-1)*x + 1)              # 1st map
 y2 = 2*y1-1                             # scale 
 y3 = ((1/s2)*y2) / (((1/s2)-1)*y2 + 1)  # 2nd map, parentheses around ((1/s2)*y2) important!  
 y4 = (y3+1)/2                           # scale
 y5 = s3*y4 / ((s3-1)*y4 + 1)            # 3rd map
 num = numerator(y5)
 den = denominator(y5)
 num.expand().collect(x), den.expand().collect(x)


-I think, maps of the special form of simpleMap constitute a subgroup of the group of the more 
 general maps with all 4 coeffs. Their composition is given by just multiplying their slopes 
 together. This is also implied by the chain rule of differentiation. The special form results 
 the constraints f(0) = 0, f(1) = 1 leading to b = 0, c = a-1. With a given slope s, we use 
 s = a/d to fix/normalize d = 1 and use s for a directly leading to the especially simple form. 

-Combining the 3 maps into one possible because of the fact that the composition of linear 
 fractional transformations gives yet another linear fractional transformation. The affine 
 transforms before and after the middle map do not destroy this because they are also special 
 cases of the linfrac maps. That means, we can compose the whole composition of the 3 maps into 
 a switch between two single maps. 

*/


//=================================================================================================

/** Computes approximate powers recursively.

see: https://www.kvraudio.com/forum/viewtopic.php?f=33&t=567563
*/

template<class T>
class rsPowerIterator
{

public:

  rsPowerIterator(T power, T epsilon)
  {
    setup(power, epsilon);
    reset();
  }


  void setup(T power, T epsilon)
  {
    p  = power;
    e  = epsilon;
    y0 = pow(e, p);
    z0 = T(1)/e;
    s  = T(1) / (pow(T(1)+e, p) - y0);  // scaler
  }

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Returns one output value at a time (and updates internal state variables for next call).
  Is dx = segmentLength / sampleRate?  */
  inline T getValue(T dx)
  {
    x += dx; 

    z -= dx*z*z; 
    //z = z*(T(1)-dx*z);  // algebraically equivalent to z -= dx*z*z;

    y += dx*z*p*y;
    //y += p*y*dx / (x+e);  // not using z, more accurate than y += dx*z*p*y;

    return s*(y-y0);
  }


  void reset()
  {
    x = T(0);
    y = y0;
    z = z0;
  }

protected:

  T x, y, z;  // state variables
  T y0, z0;   // initial values for y,z, x0 = 0
  T s;        // precomputed scaler
  T e, p;     // user parameters: "epsilon" and the power/exponent

};

// Idea: We want to compute an iterative approximation of y = x^p where x = 0..1 and p, the power, 
// is a parameter. We actually compute y = (x+e)^p where e ("epsilon") is a small number, such as 
// 0.0001. The idea is to initialize by x = 0, y = e^p and then iteratively increment x by dx and y
// by dy. To compute dy, consider 
//   dy/dx = (d/dx) (x+e)^p = p * (x+e)^(p-1) = p * (x+e)^p / (x+e) = p * y / (x+e)
// so our increment for y can be computed from y and dx as:
//   dy = p*y*dx / (x+e)
// to get rid of the division, we define:
//   z = 1 / (x+e) = (x+e)^-1
// Taking the derivative of z with respect to x:
//   dz/dx = -(x+e)^-2 = -z^2
// gives a formula for an increment for z in terms of z and dx:
//   dz = -z^2 * dx
// so the full algorithm is:
//   initialize: x = 0; y = e^p; z = 1/e;
//   iterate:    x += dx; z += -z*z*dx; y += z*p*y*dx;
// As output, we use a shifted and scaled version of y:
//   output = (y-e^p) / ((1+e)^p - e^p)
// where e^p is our y[0]


//=================================================================================================

/** Simulates the dynamics of a rotating rigid body around its three pricipal axes of intertia. If
they are all different, when it initially rotates around the axis of the middle/intermediate moment
of inertia with some small perturbation of having a rotational component around any of the other
two axes, the rotation axis periodically flips over. This is known as the "tennis racket effect"
because it also occurs when throwing up a tennis racket in a particular way. It is due to the
rotation around the intermediate axis being an unstable equilibrium of the dynamic equations that
describe the rotation. Rotation around any of the other two principal axes (those with maximum and
minimum moment of inertia) are stable equlibria.


// see:
https://en.wikipedia.org/wiki/Tennis_racket_theorem
https://en.wikipedia.org/wiki/Euler%27s_equations_(rigid_body_dynamics)
https://en.wikipedia.org/wiki/Moment_of_inertia
https://en.wikipedia.org/wiki/Poinsot%27s_ellipsoid
https://en.wikipedia.org/wiki/Polhode
https://www.youtube.com/watch?v=1VPfZ_XzisU
https://arxiv.org/pdf/1606.08237.pdf

*/

template<class T>
class rsTennisRacket
{


public:

  /** \name Setup */

  /** Sets the ratio of the moments of inertia. This determines the frequency of the flips....
  todo: have a function setFrequency */
  void setInertiaRatio(T newRatio)
  {
    // maybe wrap into if(newRatio != ratio):
    ratio = newRatio;
    I1 = ratio;
    // I2 = 1;  // always
    I3 = T(1) / ratio;
  }

  /** Sets the current state consisting of the angular velocities along the 3 principal axes. */
  void setState(T w1, T w2, T w3)
  {
    this->w1 = w1;
    this->w2 = w2;
    this->w3 = w3;
  }

  /** Sets the step size for the numerical integration scheme */
  void setStepSize(T newSize)
  {
    h = newSize;
  }


  /** \name Inquiry */

  T getW1() const { return w1; }
  T getW2() const { return w2; }
  T getW3() const { return w3; }



  /** \name Processing */

  void updateState(T M1, T M2, T M3)
  {
    // compute angular acceleration vector:
    T a1 = (M1 - (I3 - I2)*w2*w3) / I1;
    T a2 = (M2 - (I1 - I3)*w3*w1) / I2;
    T a3 = (M3 - (I2 - I1)*w1*w2) / I3;
    // formula from https://en.wikipedia.org/wiki/Euler%27s_equations_(rigid_body_dynamics)

    // todo: add damping terms...


    // update angular velocities:
    w1 += h*a1;
    w2 += h*a2;
    w3 += h*a3;

    // optionally renormalize rotational energy:
    if(normalizeEnergy) {
      T E = (I1*w1*w1 + I2*w2*w2 + I3*w3*w3); // is actually twice the energy
      T s = sqrt(T(1)/E);
      w1 *= s; w2 *= s; w3 *= s;
    }
    // maybe factor this out - client code may call it after calling setState to ensure an
    // energy normalized state - maybe include a target energy
  }


  T getSample(T in)
  {
    updateState(0, in, 0); // todo: use injection vector
    return w2;             // todo: use output vector
  }

protected:

  // user parameters:
  T ratio = 2;


  // algo parameters:
  T I1 = 2, I2 = 1, I3 = T(0.5); // moments of inertia along principal axes
  T h = T(0.01);                 // step size
  bool normalizeEnergy = true;

  // state:
  T w1 = 0, w2 = 1, w3 = 0;      // angular velocities along principal axes

};




//=================================================================================================

/** A second order (2 poles, 2 zeros) filter implementation, whose internal state is represented as
a 2D vector, i.e. a point in the xy-plane. The state is updated by multiplying the current state
vector by a matrix (and adding the input value to both components). The output is computed as a
linear combination of the state-vector's coordinates and the input. The state update matrix will
have one of these two general forms:

  [p1 0 ]     or:     r * [cos(w)  -sin(w)]
  [0  p2]                 [sin(w)   cos(w)]

where in the first case, p1 and p2 are the two real poles and the x and y states decay
exponentially and independently from each other when the input is switched off. In the second case,
the numbers r*cos(w), r*sin(w) are the real and imaginary parts of a pair of complex conjugate
poles and we will see a spiraling/decaying rotation of the states when there's no input (anymore).
The filter structure can realize almost any biquad transfer function - the singular problematic
case is when there are two equal real poles - in this case, they will be made artificially distinct
by fudging them a bit. The effect of that fudging on the transfer function will be miniscule. The
advantage of that filter structure is that it (presumably) responds well to modulation. 

\todo:
-Try to find a real solution to the edge case of two equal poles that doesn't require fudging. 
 Maybe a matrix of the form:

   [p1 0 ]
   [1  p2]

 could be used to model a series connection of the two one-poles? But then the problem arises when
 to "switch" (or maybe "fade") from parallel to serial mode - or if at all. Maybe in the case of 
 two real poles, we should always use this kind of matrix? But what would be the implications of 
 doing so with regard to state invalidation? */

template<class TSig, class TPar>
class rsStateVectorFilter
{
  typedef const TSig& CRSig;
  typedef const TPar& CRPar;

public:

  /** Sets up the filter coefficients to simulate a biquad filter with given coeffs. */
  void setupFromBiquad(CRPar b0, CRPar b1, CRPar b2, CRPar a1, CRPar a2);

  /** Sets up the two poles of this filter. You need to pass real and imaginary parts of both
  poles separately. If there are two real poles, the imaginary parts p1im, p2im should both be zero
  and if there's a complex pair, the imaginary parts should be negatives of each other, i.e p2im
  should be -p1im. The poles determine the coefficients in the state update matrix. */
  void setPoles(CRPar p1re, CRPar p1im, CRPar p2re, CRPar p2im);

  /** Assuming the poles are already fixed, this function computes the mixing coefficients such
  that the first 3 samples of the impulse response will equal what you pass to this function. This
  is used to compute the mixing coefficients after the poles have been determined. */
  void setImpulseResponseStart(TPar h[3]);

  // maybe make a setZeros function, too


  /** Produces one output sample at a time. */
  inline TSig getSample(CRSig in)
  {
    updateState(in);
    return cx*x + cy*y + ci*in;
  }

  /** Resets the filter state. */
  void reset()
  {
    x = y = TSig(0);
  }

protected:

  /** Used internally in getSample to update the filter state. */
  inline void updateState(CRSig in)
  {
    TSig t = x;             // temporary
    x = xx*x + xy*y + in;   // update x
    y = yx*t + yy*y + in;   // update y
  }

  /** This is a function to fudge with the poles in cases where they are (almost) equal. Such a
  case cannot be represented exactly by this filter structure (a singular matrix in the mixing
  coefficient calculation would occur), so we use distinct poles close to the originally desired
  poles. The effect is a slight misadjustment of the filter in these particular cases. */
  void makePolesDistinct();
   // maybe return a bool to inform, if the poles were modified, maybe also return a bool from
   // setPoles in order to be able to make client code aware of the fudging

  TPar xx = 0, xy = 0, yx = 0, yy = 0;  // matrix coeffs
  TPar cx = 0, cy = 0, ci = 1;          // mixing coeffs
  TSig x  = 0, y  = 0;                  // state vector

};

//=================================================================================================

/** Another go at a general ordinary differential equation solver with a probably more convenient
interface than the old one (which required subclassing to define a concrete ODE (system)).

References
(1) Numerical Recipies
(2) Mathematik, Ahrens et al, 4.Aufl.

-maybe rename to rsExplicitInitialValueSolver (and provide also an implicit one)
-maybe factor out a solver that doesn't carry around x and where f only depends on y - or maybe
 subsume systems that depend explicitly on x by incorporating an identity function into the vector
 of functions f(y), i.e. f(y1, y2, y3, ...) = (y1, f2(y1,y2,y3..), f3(y1,y2,y3..), ...), the
 derivative of y1 is always 1, so x += h translates to y1 += 1*h in the new notation. that would
 simplify interface and implementation but requires more understanding from the user and does not
 allow to have a different datatype for x
-maybe move the state variables to a subclass (rsMultiStepInitialValueSolver or something)
*/

template<class Tx, class Ty>
class rsInitialValueSolver
{

public:


  /** \name Setup */

  void setStepSize(Tx newSize) { h = newSize; }




  /** \name Evaluation */

  Ty getSampleForwardEuler()
  {
    f0 = f(x, y);
    x += h;
    y += h*f0;
    return y;
  }

  Ty getSampleRungeKutta4()
  {
    Ty k1 = f(x,         y);
    Ty k2 = f(x + 0.5*h, y + 0.5*h*k1);
    Ty k3 = f(x + 0.5*h, y + 0.5*h*k2);
    Ty k4 = f(x +     h, y +     h*k3);
    x += h;
    y += (h/6.) * (k1 + 2*k2 + 2*k3 + k4); // optimize away division
    return y;
  }


  Ty getSampleAdamsBashforth2()
  {
    f0 = f(x, y);  // new evaluation
    x += h;
    y += 0.5*h*(3*f0 - f1);
    f1 = f0;
    return y;
  }
  // todo: make sure f0, f1, are initialized correctly (do this by doing an RK step)
  // see (2) page 481


  Ty getSampleAdamsBashforth4()
  {
    f0 = f(x, y);  // new evaluation
    x += h;
    y += (h/24.) * (55*f0 - 59*f1 + 37*f2 - 9*f3); // optimize away division
    f3 = f2;
    f2 = f1;
    f1 = f0;
    return y;
  }

  // rsCubicSpline<Tx, Ty> getSolution(Tx x0, Tx x1, double accuracy);
  // should produce an object of class rsCubicSpline that represents the solution



protected:

  Tx x = 0;
  Ty y = 0;
  Tx h = 1;  // step size

  std::function<Ty(Tx,Ty)> f; // this is the "f" in y'(x) = f(x,y)

  // state variables for multistep methods:
  Ty f0, f1, f2, f3, f4; // f(x[n],y[n]), f(x[n-1],y[n-1]), ...

};

//=================================================================================================

// -maybe make also a class rsBinarySearchTree where the left child is <= and the right child is >=
//  the parent
// -pehaps we could implement a general function: needsSwap(int parentIndex, int childIndex) - in
//  the case of a search tree, it would first have to figure out, if the childIndex is a left or
//  right child and order the arguments of the comparison according to that
// -in the old RSLib codebase, i did some sort of linked-tree - maybe that could be dragged in as
//  rsLinkedTree or rsDynamicTree or something like that - all the different trees could be in
//  a file Trees.h/cpp
// -i think, currently, the order of the children of a node is undefined - both children are <=
//  the parent, but maybe we could impose additional useful structure, if we also have
//  left <= right - would that be easy to implement and/or useful?
// -maybe implement also a MinMaxHeap https://en.wikipedia.org/wiki/Min-max_heap

/** Works similar to rsBinaryHeap just that the property satisfied by the nodes is different.
Here, it holds that:

   A[left(i)]  <= A[i]    and
   A[right(i)] >= A[i]

As a reminder, for a (max)heap, the property was A[left(i)] <= A[i] and A[right(i)] <= A[i], i.e.
both children must have values less-or-equal than the parent. Here, the left node must be
less-or-equal and the right node must be greater-or-equal. */

template<class T>
class rsBinarySearchTree : public rsBinaryTree<T>
{

public:

  using rsBinaryTree<T>::rsBinaryTree;  // inherit constructors


  bool isSearchTree(int i = 0) const
  {
    if(i >= this->size)
      return true;
    bool result = true;
    int l = this->left(i);
    int r = this->right(i);
    if(l < this->size) result &= !this->less(this->data[i], this->data[l]) && isSearchTree(l);
    if(r < this->size) result &= !this->less(this->data[r], this->data[i]) && isSearchTree(r);
    return result;
  }
  // actually, this test is not enough - it says yes to the "pseudotree"
  // 50,20,80,10,60,30,90 - it satisfies the property at each node with respect to direct
  // parents/children but not for the full subtrees - however - maybe for maintaining the
  // property, it's enough to check that? can replace create pseudotrees? ..oh - yes - that seems
  // to be the case! damn! maybe the whole idea of using the same structure as for the binary heap
  // does not work out as i hoped... :-(

  void buildSearchTree()
  {
    for(int i = this->size-1; i >= 0; i--)
    {
      //floatUp(i);
      //floatIntoPlace(i);
      //floatDown(i);
      fixNode(i);
    }
    // is it vital to do this in reverse order?
    // i think, we can't just let the nodes i float into place - we must call some sort of
    // fixNode function that either swaps p,l or p,r or l,r, or nothing
  }


  void fixNode(int i)
  {
    //floatIntoPlace(i);  // does not work

    int l = this->left(i);
    if(l >= this->size) return;       // node is leaf
    int r = this->right(i);
    if(r >= this->size) {             // node has only left child - make sure, the data at the
      if(this->less(this->data[i], this->data[l]))      // parent i is not less than at left child l - if it is: swap
        this->swap(this->data[i], this->data[l]);
      return; }

    // ok - we have a node that has both children - we must figure out which of the 3 nodes i,l,r
    // is the middle element m and if m is not already i, then do a swap. we must also make sure
    // that data[l] <= data[r]. we have 3 possible swapping scenarios: swap(i,l), swap(i,r),
    // swap(l,r) ...plus, of course, the no-swap scenario

    //int m = i;                           // index of mid element, preliminary set to i
    //if(less(data[i], data[l])) m = l;
    //if(less(data[r], data[m])) m = r;

    // fix order of children:
    if(this->less(this->data[r], this->data[l]))
      this->swap(this->data[l], this->data[r]);

    // fix order of i,l:
    if(this->less(this->data[i], this->data[l]))
      this->swap(this->data[l], this->data[i]);

    // fix order of i,r:
    if(this->less(this->data[r], this->data[i]))
      this->swap(this->data[r], this->data[i]);

    // try to do this only with 2 comparisons and one swap

  }



  int replace(int i, const T& x)
  {
    // data[i] = x; return floatIntoPlace(i); // this is what the heap needs

    // we insert it (preliminarily) either at i or at the sibling of i:
    int p = this->parent(i);
    if(this->isLeft(i))
      if(this->less(this->data[p], x))
        i = this->right(p);   // right sibling
    else
      if(this->less(x, this->data[p]))
        i = this->left(p);    // left sibling
    this->data[i] = x;
    return floatIntoPlace(i);
  }


protected:

  int floatIntoPlace(int i) { return floatUp(floatDown(i)); }

  int floatUp(int i)
  {
    while(i > 0) {
      int p = this->parent(i);
      if(this->isLeft(i)) {
        if(this->less(this->data[p], this->data[i]))
          this->swap(this->data[i], this->data[p]);
        else
          return i; }
      else {
        if(this->less(this->data[i], this->data[p]))      // i and p are swapped compared to left nodes
          this->swap(this->data[i], this->data[p]);
        else
          return i;   }
      i = p; }
    return i;
  }

  int floatDown1(int i)
  {
    int l = this->left(i);
    if(l >= this->size) return i;
    if(this->less(this->data[i], this->data[l])) {
      this->swap(this->data[i], this->data[l]);
      return floatDown(l); }
    int r = this->right(i);
    if(r >= this->size) return i;
    if(this->less(this->data[r], this->data[i])) {
      this->swap(this->data[i], this->data[r]);
      return floatDown(r); }
    return i;
  }
  // recursive implementation - not working


  int floatDown2(int i)
  {
    while(i < this->size-1)
    {
      int l = this->left(i);
      if(l >= this->size) return i;              // node i is a leaf
      int m = i;                           // index of mid element, preliminary set to i
      if(less(this->data[i], this->data[l])) m = l;    // left child is bigger than i, so mid index goes to l
      int r = this->right(i);

      /*
      if(r >= size)
      {
        //rsError("not yet implemented");

        if(m != i)
        {
          swap(data[m], data[i]);
          i = m;
        }
        continue;
      }
      else
      {
        if(less(data[r], data[m]))
          m = r;
      }
      */

      if(r < this->size && less(this->data[r], this->data[m])) m = r;

      if(m != i) {
        swap(this->data[m], this->data[i]);
        i = m;  }
      else
        return i;

      /*
      if(r < size && less(data[b], data[r])) b = r;

      if(b != r) {              // difference to heap: there it was: if(b != i)
        swap(data[i], data[b]);
        i = b;  }
      else
        return i;
        */


      //// no - this is wrong: we need to consider 3 cases:
      ////   swap l with i, swap r with i, swap r with l
      //// see  the heap implementation
      //if(l < size && less(data[i], data[l])) {
      //  swap(data[i], data[l]);
      //  i = l; }
      //else if(r < size && less(data[r], data[i])) {
      //  swap(data[i], data[r]);
      //  i = r; }
      //else
      //  return i;


    }
    return i;

    // hmm..i think, we have to figure out the index of the biggest and smallest element b,s among
    // i,l,r. if s==l and b==r, the node i should stay where it is
  }

  int floatDown(int i)
  {
    return floatDown1(i);
    //return floatDown2(i);
  }

};

//=================================================================================================

/** Extends rsQuantileFilter by a second core allowing to do more interesting stuff such as forming
linear combinations of lower an upper quantiles (such as min and max), etc. Using a second core
instead of just using two rsQuantileFilter objects is more economical because the delayline can be
shared between the cores. */

template<class T>
class rsDualQuantileFilter : public rsQuantileFilter<T>
{

public:

  rsDualQuantileFilter()
  {
    allocateResources();
    this->core.setDelayBuffer(&this->delayLine);
    this->core2.setDelayBuffer(&this->delayLine);
    this->dirty = true;
  }


  // setters for the parameters of the second core
  void setFrequency2(   T newFrequency) { setLength2(T(1) / newFrequency); }
  void setLength2(      T newLength)    { length2   = newLength;    this->dirty = true; }
  void setQuantile2(    T newQuantile1) { quantile2 = newQuantile1; this->dirty = true; }
  void setLowpassGain2( T newGain)      { loGain2   = newGain; }
  void setHighpassGain2(T newGain)      { hiGain2   = newGain; }

  void setCore2Equal()
  {
    length2   = this->length;
    quantile2 = this->quantile;
    loGain2   = this->loGain;
    hiGain2   = this->hiGain;
    delayScl2 = this->delayScl;
    this->dirty = true;
  }

  void setCore2Complementary()
  {
    setCore2Equal();
    quantile2 = T(1) - this->quantile;
    this->dirty = true;
  }


  // maybe have a setCore2Complementary which uses the same settings for core2 as core2 except the
  // quantile, for which we use: quantile2 = 1-quantile

  T getSample(T x)
  {
    this->delayLine.getSample(x);
    if(this->dirty)
      updateInternals();
    T yL1 = this->core.getSample(x);
    T yH1 = this->delayLine[this->delayScl*this->delay] - yL1;
    T yL2 = core2.getSample(x);
    T yH2 = this->delayLine[delayScl2*delay2] - yL2;
    return this->loGain * yL1 + this->hiGain * yH1 + loGain2 * yL2 + hiGain2 * yH2;
  }

  void reset() { this->core.reset(); core2.reset(); this->delayLine.reset(); }

  virtual void updateInternals() override
  {
    // rsError("has to be updated");
    // 
    /*
    // compute internal and set up core parameters:
    int L, p; T w;

    this->convertParameters(this->length, this->quantile, this->sampleRate, &L, &p, &w, &this->delay);
    this->core.setLengthAndReadPosition(L, p);
    this->core.setRightWeight(w);

    this->convertParameters(length2, quantile2, this->sampleRate, &L, &p, &w, &delay2);
    core2.setLengthAndReadPosition(L, p);
    core2.setRightWeight(w);
    */

    // update inherited baseclass members:
    T L = this->length * this->sampleRate;
    this->core.setLengthAndQuantile(L, this->quantile);
    this->delay = T(0.5)*(L-1);
    //this->delay = T(1.0) * (L-1) * this->quantile;  // experimental

    // update new subclass members:
    L = length2 * this->sampleRate;
    core2.setLengthAndQuantile(L, quantile2);
    delay2 = T(0.5)*(L-1);
    //delay2 = T(1.0) * (L-1) * quantile2;   // experimental

    this->dirty = false;
  }


protected:

  virtual void allocateResources() override
  {
    int mL = (int) ceil(this->maxLength * this->sampleRate);
    this->core.setMaxLength(mL);
    core2.setMaxLength(mL);
    this->delayLine.setCapacity(mL);
  }


  // the 2nd core and its set of parameters:
  rsQuantileFilterCore2<T> core2;
  T length2   = 0.0;
  T quantile2 = 0.5;
  T loGain2   = 0.0;
  T hiGain2   = 0.0;
  T delayScl2 = 1.0;
  T delay2    = 0.0;

};

//=================================================================================================

/** Extends rsQuantileFilter by producing a pseudo-resonance signal formed by... */

template<class T>
class rsQuantileFilterResonant : public rsQuantileFilter<T>
{

public:

  using Base = rsQuantileFilter<T>;

  rsQuantileFilterResonant()
  {
    //resoHighpass.setMode(rsOnePoleFilter<T,T>::HIGHPASS_MZT);
    resoHighpass.setMode(rsStateVariableFilter<T,T>::HIGHPASS);
    allocateResources();
    Base::dirty = true;
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the sample rate at which this filter should operate. May re-allocate memory. */
  void setSampleRate(T newSampleRate)
  {
    this->sampleRate = newSampleRate;
    allocateResources();
    resoHighpass.setSampleRate(this->sampleRate);
    Base::dirty = true;
  }

  /** Sets the maximum length (in seconds) for this filter. May re-allocate memory. */
  void setMaxLength(T newMaxLength)
  {
    this->maxLength = newMaxLength;
    allocateResources();
    Base::dirty = true;
  }

  /** Sets sample rate and maximum length at the same time. May re-allocate memory. This may avoid
  some re-allocations compared to using setSampleRate followed by setMaxLength (or vice versa), 
  depending on the situation - so, if possible, it's recommended to set both at the same time. */
  void setSampleRateAndMaxLength(T newSampleRate, T newMaxLength)
  {
    this->sampleRate = newSampleRate;
    this->maxLength  = newMaxLength;
    allocateResources();
    Base::dirty = true;
  }

  void setResonanceMix(T newMix)
  {
    resoMix = newMix;
  }
  // maybe use a resoGain instead

  /*
  enum class ResoMode
  {


  };
  */


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Produces one output sample from a given input sample. */
  T getSample(T x)
  {
    this->delayLine.getSample(x);
    if(this->dirty)
      updateInternals();
    T yL = this->core.getSample(x);                       // lowpass part
    T yD = this->delayLine[this->delayScl * this->delay]; // delayed input
    T yH = yD - yL;                                       // highpass part
    T yF = this->loGain * yL + this->hiGain * yH;         // non-resonant filtered output
    T yR = getResonance(x, yL, yD);                       // (pseudo) resonance
    return (T(1)-resoMix) * yF + resoMix * yR;            // crossfade between non-reso and reso
  }
  // maybe factor out a function to produce lowpass and highpass getSampleLoHi or something at the
  // same time - client code may find that useful - or maybe getOutputs to be consistent with
  // rsStateVariableFilter

  /** Resets the filter into its initial state. */
  void reset()
  {
    Base:reset();
    minCore.reset();
    maxCore.reset();
    bandpass.reset();
    resoHighpass.reset();
  }

  /** Updates the internal algorithm parameters and embedded objects according to the user
  parameters. This is called in getSample, if the state is dirty but sometimes it may be
  convenient to call it from client code, too. */
  virtual void updateInternals()
  {
    double L = this->length * this->sampleRate;
    this->core.setLengthAndQuantile(L, this->quantile);

    minCore.setLengthAndQuantile(L, T(0));
    maxCore.setLengthAndQuantile(L, T(1));
    //minCore.setLengthAndQuantile(100, T(0));
    //maxCore.setLengthAndQuantile(100, T(1));
    //minCore.setLengthAndQuantile(rsMax(L, 10.0), T(0));
    //maxCore.setLengthAndQuantile(rsMax(L, 10.0), T(1));


    this->delay = T(0.5)*(L-1);

    T f = this->getFrequency();

    //T w = T(2*PI)*sampleRate/length; // == 2*PI*sampleRate*frequency
    T w = T(2*PI) * f  / this->sampleRate; // == 2*PI*frequency/sampleRate
    bandpass.setFrequencyAndAbsoluteBandwidth(w, T(0.00001));  // preliminary
    // todo: let the resonance frequency have an adjustable pitch-offset with respect to core 
    // filters...maybe the min/max cores could have an offset as well
    // and maybe we should use a power rule to generalize from constant absolute bandwidth (0) and
    // constant relative bandwidth (1) to anything between and beyond...and we generally need a 
    // reasonable user-parameter to control the bandwidth

    //T gain = bandpass.getMagnitudeAt(w);  // excluding the output gain
    //bandpass.setOutputGain(T(1)/gain);
    // for test/debug 
    // -gain tends to get huge for small bandwidths (seems to be in reciprocal 
    //  relationship)! 
    // -todo: plot impulse response

    //T fHp = rsMax(T(0), 0.9*f - T(1000));
    T fHp = rsMax(T(0), 0.5*f);
    //resoHighpass.setCutoff(fHp);
    resoHighpass.setFrequency(fHp);


    this->dirty = false;
  }



protected:

  /** Computes the pseudo-resonance signal. Inputs are the filter input x, the "lowpass" quantile
  filter output yL and the delayed input yD (rename to xD) which are used to generate the 
  resonance. */
  T getResonance(T x, T yL, T yD)
  {
    T min  = minCore.getSample(x);
    T max  = maxCore.getSample(x);

    // this computation should be factored out and we should provide various different ways to
    // compute wMin/wMax that can bew switched by the user using a setResonanceMode function:
    T wMin, wMax;
    wMin = 0;                              // preliminary.. maybe these should be members and 
    wMax = 0;                              // computed in updateinternals

    //if(yD >= yL) wMax = 1;
    //else         wMin = 1;

    // this is not a very good way to compute wMin/wMax - it should use some condition that is 
    // likely to cause a switch once per cycle, i.e. on the avarage switches once within a length
    // of the filter buffer. this condition here switches far too often - try several things
    // ...oh - but maybe that's just because the input is white noise - if it's somewhat lowpassed
    // noise, i may be better

    computeMinMaxWeights(&wMin, &wMax, x, yL, yD, min, max);
    T yR = wMin*min + wMax*max;

    // experimental - highpass and amplify the resonance to counteract resonance loss at high 
    // cutoff frequencies:
    yR  = resoHighpass.getSample(yR);
    T a = T(1);    // preliminary - should later depend on cutoff in some way

    // experimental:
    T p = 2*this->getFrequency()/this->sampleRate;  // 0..1 as f goes 0..fs
    //p = 1-p; p = p*p; a = 1 / (0.1 + p);
    //p = cos(0.5*PI*p); a = 1 / (0.1 + p);
    p = 1-sin(0.5*PI*p); a = 1 / (0.1 + p);
    // todo: figure out the best shape


    return a*yR;
  }

  void computeMinMaxWeights(T* wMin, T* wMax, T x, T yL, T yD, T yMin, T yMax)
  {
    T yB = bandpass.getSample(yD);  // or maybe feed x

    T thresh = T(0.0) * yL;
    //T thresh = T(0.2) * rsMax(rsAbs(yMin), rsAbs(yMax));
    // make user-adjustable, determines pulse-width, should perhaps be scaled by yL...and maybe we
    // should normalize the amplitude of the bandpass with respect to the bandwidth (we probably 
    // want a constant peak gain bandpass - maybe the RBJ bandpass is good for that)
    // ...it should take into account the magnitude of the bandpass output which may be a function
    // of the bandwidth - we shoould normalize everything such that thresholds between -1..+1 make
    // sense and have the same effect regardless of input volume and bandwidth setting

    T a = T(1);  // scaler to control the resonance amplitude
    T f = this->getFrequency();
    //a = 10 / length;
    //a = T(1) + T(4)*rsClip(f / T(20000), T(0), T(1)); 
    // yes, we need a factor, but we also need a highpass, otherwise we just amplify noise
    // maybe the highpass freq should be something like max(0, cutoff-2000) such that the highpass
    // sets in above 2kHz cutoff an is tuned 2kHz below the cutoff...of maybe use 
    // max(0, 0.5*cutoff-2000) or something - tweak the factor and offset to taste

    if(yB >= thresh) { *wMin = T(0); *wMax = a;    }
    else             { *wMin = a;    *wMax = T(0); }
    // maybe instead of hard-switching, we can make a sort of soft-switch?


    /*
    if(yD >= yL) { *wMin = T(0); *wMax = T(1); }
    else         { *wMin = T(1); *wMax = T(0); }
    */
  }

  virtual void allocateResources() override
  {
    Base::allocateResources();
    int mL = this->getMaxRequiredLengthInSamples();
    minCore.setMaxLength(mL);
    maxCore.setMaxLength(mL);
  }

  // filter cores to extract min and max:
  rsQuantileFilterCore2<T> minCore, maxCore;
  // can this be done more efficiently? i could perhaps use the rsMovingMinMaxFilter class but 
  // then the length would not be modulatable and i don't think it's easy to make it modulatable

  T resoMix = T(0);

  rsTwoPoleFilter<T, T> bandpass;
  // the two-pole is actually not really a bandpass but a resonator - maybe an RBJ bandpass is 
  // better?

  //rsOnePoleFilter<T, T> resoHighpass;

  rsStateVariableFilter<T, T> resoHighpass;

};



//=================================================================================================

/** This is a naive implementation of (the core of) a moving quantile filter and meant only for
producing test outputs to compare the production version of the rsQuantileFilterCore against.
It's horribly inefficient - the cost per sample is O(N*log(N)) whereas the production version
should run in O(log(N)). ..maybe move to unit tests.. */

template<class T>
class rsQuantileFilterNaive
{

public:

  rsQuantileFilterNaive(int numSmaller = 20, int numLarger = 20)
    : buf(numSmaller+numLarger)
  {
    setLengths(numSmaller, numLarger);
  }

  void setLengths(int numSmaller, int numLarger)
  {
    nS = numSmaller;
    nL = numLarger;
    updateBufferLengths();
  }

  int getLength() const { return nS + nL; }


  void prepareSortedDelayBuffer(T x)
  {
    T y = buf.getSample(x);
    buf.copyTo(&tmp[0], true);          // O(N)
    rsHeapSort(&tmp[0], getLength());   // O(N*log(N))
  }

  T getSample(T x)
  {
    prepareSortedDelayBuffer(x);
    return tmp[nS];
    // pehaps, we should use interpolation here - provide functions setLength, setQuantile
  }

  T getSampleMedian(T x)
  {
    prepareSortedDelayBuffer(x);
    int L = getLength();
    T p = 0.5 * (L-1);
    int i = (int) floor(p);
    T   f = p - i;
    //return f*tmp[i] + (1-f)*tmp[i+1]; // this is wrong, but matches rsQuantileFilter
    return (1-f)*tmp[i] + f*tmp[i+1]; // this is correct and matches rsArrayTools::median
    // verify, if we use f and 1-f correctly
  }


  void reset()
  {
    buf.reset();
    //for(int i = 0; i < getLength(); i++)
    //  buf[i] = 0;
  }


protected:

  void updateBufferLengths()
  {
    buf.setLength((size_t)getLength());
    tmp.resize(   (size_t)getLength());
    //reset();
  }

  rsDelayBuffer<T> buf;   // circular buffer
  std::vector<T> tmp;
  int nS = 0, nL = 0;    // maybe use size_t

};


//=================================================================================================

/** This is a naive implementation of a "filter" that minimizes the distance between subsequent 
samples by re-ordering them. It always keeps a buffer of N samples from the past and when a new 
sample comes in, it replaces one of the buffered samples with the new one and the replaced sample 
from the buffer is used as output. Which one is replaced and used as output is one that is closest 
to the previous output sample. In the special case that the current input is closer to the previous
output than all buffered samples, it will be used as output and the buffer will left alone. */

template<class T>
class rsDistanceFilterNaive
{

public:

  void setLength(int newLength)
  {
    buf.resize(newLength);
    rsFill(buf, T(0));
  }
  // O(N)

  T getSample(T x)
  {
    T   dMin = rsAbs(y - x);
    int iMin = -1;
    for(int i = 0; i < (int) buf.size(); i++) {
      T d = rsAbs(y - buf[i]);
      if(d < dMin) {
        dMin = d;
        iMin = i; }}

    if(iMin != -1) {
      y = buf[iMin];
      buf[iMin] = x; }
    else
      y = x;

    return y;
  }
  // O(N)

  void reset()
  {
    rsFill(buf, T(0));
    y = T(0);
  }
  // O(N)

protected:

  std::vector<T> buf;     // buffer of stored samples
  T y = T(0);             // old output

};
// ..not yet tested
// todo:
// -implement this more efficiently using a heap: 
//  -we need a heap-search procedure which should run in log(N) time): int heap.findClosest(T x)
//  -this should be used to determine the sample to be replaced, replacement itself will also take
//   O(log(N)), so the overall complexity of getSample will be O(log(N))
//  ...oh - no, maybe a heap will not work and we need a binary search tree instead
// -use it to filter uniform white noise - it should preserve the uniform amplitude distribtion
//  while still imposing a correlation (that was actually the question that lead to the idea: how 
//  can we produce correlated noise with uniform amplitude distribution - because regular filters
//  tend to gaussianize it)
// -the naive version could generalize to other distance measures such as 2D Euclidean distance
//  for stereo signal - unfortunatley, i don't see, how this could be made efficient using a heap 
//  because there's no natural ordering for 2D vectors...but maybe one could be invented...but it 
//  should be related to distance between pairs of vectors

// ...how else could we impose correlations (and therefore, some sort of not-flat frequency 
// spectrum) without modifying the amplitude distribution? the amplitude distribution can be 
// modified by waveshaping and correlations can be introduced by filtering - so maybe a filtering
// process (like 2 or 3 point MA) can be used to spectrally shape the noise, which also turns a 
// uniform distribtuion into triangular or parabolic and then that can be followed by waveshaper 
// that counteracts this change (square? cube? sqrt?, cbrt?) ...or maybe it should be the 
// (inverse of?) the integral of the resulting distribution

//=================================================================================================

template<class T>
class rsNoiseGeneratorTriModal
{

public:

  rsNoiseGeneratorTriModal()
  {
    selectorLowpass.setMode(selectorLowpass.LOWPASS_IIT);
    selectorLowpass.setSampleRate(1);
    selectorLowpass.setCutoff(0.1);
    selector.setRange(-1, +1);
    ng1.setRange(-1.0, -0.3);
    ng2.setRange(-0.3, +0.3);
    ng3.setRange(+0.3, +1.0);
    setOrder(7);
  }

  void setOrder(int newOrder)
  {
    ng1.setOrder(newOrder);
    ng2.setOrder(newOrder);
    ng3.setOrder(newOrder);
  }

  inline T getSample()
  {
    T s = selector.getSample();
    s = selectorLowpass.getSample(s);
    if(s < thresh1)
      return ng1.getSample();
    if(s < thresh2)
      return ng2.getSample();
    return ng3.getSample();
  }

  rsOnePoleFilter<T, T> selectorLowpass;
  // the selector is lowpassed such that successive samples tend to be selected from the same
  // distribution


protected:

  rsNoiseGenerator<T> selector;
  rsNoiseGenerator2<T> ng1, ng2, ng3;
  T thresh1 = -0.2, thresh2 = +0.2;

};

//=================================================================================================
// continued fraction stuff:

/** A class for generating the (integer) continued fraction expansion coefficients of a given
(floating point) number. You pass the number to the constructor and after that, you can pull out
successive cofficients via getNext(). */

template<class TInt, class TFloat>
class rsContinuedFractionGenerator
{

public:

  rsContinuedFractionGenerator(TFloat number) : state(number) {}

  TInt getNext()
  {
    TFloat flr = floor(state);
    state = TFloat(1) / (state - flr);
    return (TInt) flr;
  }
  // maybe rename to getNextCoeff ...or whatever these numbers are called


protected:

  TFloat state;

};
// https://www.youtube.com/watch?v=CaasbfdJdJg


template<class T>
rsFraction<T> rsContinuedFractionConvergent(T* a, int N)
{
  T p0 = 0, p1 = 1, p2 = 0;
  T q0 = 1, q1 = 0, q2 = 1;
  for(int i = 0; i < N; ++i) {
    p2 = a[i]*p1 + p0; p0 = p1; p1 = p2;
    q2 = a[i]*q1 + q0; q0 = q1; q1 = q2; }
  return rsFraction<T>(p2, q2);
}
// algorithm adapted from cfcv.c (by Hollos) - i don't really know, why it works
// this can actually be done directly using the generator, without the need for explicitly
// computing and storing the array a
// maybe move into class rsContinuedFractionGenerator...maybe as static method
// ..but i don't think that this continued fraction stuff should go into rsFraction - it's stuff
// on top of it
// -note the the convegents are not eqaul to the best approximants. there's some additional stuff
//  that needs to be done - maybe implement that in a function rsRationalApproximant:
//  https://en.wikipedia.org/wiki/Continued_fraction#Best_rational_approximations
//
// can we somehow figure out, how many of the CFE coeffs are correct without knowing the correct
// CFE? maybe by converting the convergents back to double and only add more coeffs as long as
// the back-converted number actually gets closer to the original number?
// -maybe rename to rsContinuedToRegularFraction - it just converts an array of (simple) continued
//  fraction coeffs to a normal, regular fraction - it doesn't really matter whether or not the
//  coefficient array is supposed to converge to some irrational number or if it's just another
//  representation of a fraction
// -maybe also implement a generalization that does not assume a *simple* continued fraction (i.e.
//  one where all numerators are 1 and there are no minusses)
// https://en.wikipedia.org/wiki/Continued_fraction
// https://en.wikipedia.org/wiki/Generalized_continued_fraction

// returns the (simple) continued fraction coeffs of the given rational number
template<class T>
std::vector<T> rsContinuedFraction(rsFraction<T> x)
{
  std::vector<T> c;
  T p = x.getNumerator();
  T q = x.getDenominator();
  T a = p/q, b;
  c.push_back(a);
  while(p > q*a) {
    b = p - q*a; p = q; q = b; a = p/q;
    c.push_back(a); }
  return c;
}
// algo adapted from from cfrat.c
// i think, this is some variation of the Euclidean algorithm

// Wikipedia says: "Even-numbered convergents are smaller than the original number, while
// odd-numbered ones are larger."
// https://en.wikipedia.org/wiki/Continued_fraction#Infinite_continued_fractions_and_convergents
// ...verify this. Knowing this could be useful

// This is interesting for doing computations with continued fractions:
// https://srossd.com/posts/2020-09-18-gosper-1/

// other fun stuff that can be done with fractions:
// https://en.wikipedia.org/wiki/Egyptian_fraction


//=================================================================================================

template<class T>
void fitQuadratic(T x1, T y1, T x2, T y2, T x3, T y3, T* a0, T* a1, T* a2)
{
  T k1 = y1 / ((x1-x2)*(x1-x3));
  T k2 = y2 / ((x2-x1)*(x2-x3));
  T k3 = y3 / ((x3-x1)*(x3-x2));
  T b1 = -k1*(x2+x3);
  T b2 = -k2*(x1+x3);
  T b3 = -k3*(x1+x2);
  T c1 = k1*x2*x3;
  T c2 = k2*x1*x3;
  T c3 = k3*x1*x2;
  *a2  = k1 + k2 + k3;  // coeff for x^2
  *a1  = b1 + b2 + b3;  // coeff for x^1
  *a0  = c1 + c2 + c3;  // coeff for x^0

  // Formulas were derived from setting up 3 polynomials in product form, where each has zeros at
  // all but one of the datapoints, say xi, and to have value yi at xi and then adding them up
  // (idea due to Lagrange):
  //   p1(x) = k1*(x-x2)*(x-x3)       p1 has zeros at at x2,x3
  //   p2(x) = k2*(x-x1)*(x-x3)       p2 has zeros at at x1,x3
  //   p3(x) = k3*(x-x1)*(x-x2)       p3 has zeros at at x1,x2
  // Require:
  //   p1(x1) = y1, p2(x2) = y2, p3(x3) = y3
  // Solve these for the ki, i.e. k1,k2,k3. For example, k1 = y1 / ((x1-x2)*(x1-x3)). Plug, for
  // example, k1 back into the p1 equation and multiply it out to obtain its coeffs - do the same
  // for p2 and p3 and then obtain the final polynomial coeffs by adding the corresponding  coeffs
  // of each of the partial polynomials.

  // operations: add: 9, sub: 6, mul: 12, div: 3, neg: 3, ass: 12, tmp: 9
}

//=================================================================================================

/** Sawtooth oscillator based on a PolyBlep. Convenience class that also demonstrates, how the 
bleps are supposed to be used in a quite minimal setting. */

template<class T>
class rsBlepSawOsc : protected rsBlepReadyOsc<T>
{

public:

  using Base = rsBlepReadyOsc<T>;

  void setPhaseIncrement(T newInc) { Base::setPhaseIncrement(newInc); }

  T getSample()
  {
    T stepDelay, stepAmp;
    T y = Base::getSampleSaw(&stepDelay, &stepAmp); // todo: maybe switch waveform later
    if(stepAmp != 0.0)
      blep.prepareForStep(stepDelay, stepAmp);
    return blep.getSample(y);
  }

  void reset(T start = T(0))
  {
    Base::resetPhase(start);
    blep.reset();
  }

protected:

  rsPolyBlep2<T, T> blep;

};

//=================================================================================================

/** A class for representing and performing computations with sparse matrices which are matrices in 
which most elements are zero. This implementation uses a std::vector of "elements" where each 
element stores its row- and column-indices and the actual value. 

The class is meant to be used mostly in situations where computing matrix-vector products is common
and random access of matrix elements is uncommon. In a regular (dense) matrix implementation of 
shape MxN, the matrix-vector product is an O(M*N) operation and element access is O(1). Here, 
computing the matrix-vector product is an O(K) operation (where K is the number of nonzero entries) 
and element random access is O(log(K)) for reading an element. This is because the elements are 
stored as sorted array and random access requires binary search. In the product, we can just 
iterate over all the stored elements. 

ToDo:
-Document complexity of element write access. I think, it should also be O(log(K)) when we 
 overwrite an existing element and O(K) when we insert a new element (or remove one - removal 
 should be triggered on write when the caller writes a zero). */

template<class T>
class rsSparseMatrix
{


public:



  // todo: provide a (subset of) the same interface as rsMatrix(View)...maybe it's a good idea to 
  // keep the interface small to make it easier to provide different implementations with different
  // storage modes with the same interface that can be benchmarked against each other.

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  rsSparseMatrix() {}

  rsSparseMatrix(int numRows, int numColumns)
  {
    rsAssert(numRows >= 1 && numColumns >= 1); 
    this->numRows = numRows; 
    this->numCols = numColumns;
  }

  /** Copy constructor. */
  rsSparseMatrix(const rsSparseMatrix<T>& A)
  {
    numRows  = A.numRows;
    numCols  = A.numCols;
    elements = A.elements;
  }

  /** Creates a rsSparseMatrix from a regular (dense) rsMatrix. You can optionally set a relative 
  tolerance for the absolute value below which elements in A will be considered zero and not be 
  stored in th sparse representation. The tolerance is relative to the maximum absolute value in 
  A and it is zero by default. */
  static rsSparseMatrix<T> fromDense(const rsMatrix<T>& A, T tolRel = T(0));

  /** Creates a regular (dense) rsMatrix from a rsSparseMatrix. */
  static rsMatrix<T> toDense(const rsSparseMatrix<T>& A);
  // ToDo: explain rationale for making this a static function taking A as argument rather than
  // a member that can be called like "rsMatrix<T> D = A.toDense();". I think, it's for making
  // the call consistent with "fromDense()" which cannot be called like that unless we move it 
  // into rsMatrix which would couple rsMatrix to rsSparseMatrix which is undesirable. rsMatrix
  // should not depend on rsSparseMatrix. Having rsSparseMatrix dependent on rsMatrix is 
  // acceptable, though. I think of rsMatrix as the more fundamental class and rsSparseMatrix as
  // something on top of it.


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of nonzero elements in this matrix. */
  int getNumElements() const { return (int) elements.size(); }
  // Maybe rename to getSize() to match rsMatrix ...but maybe it's intentional to not match it 
  // because for a sparse matrix, the notion of "size" is a bit ambiguous because the number of 
  // stored elements is not predetermined?.

  int getNumRows() const { return numRows; }

  int getNumColumns() const { return numCols; }

  /** Returns true, iff this matrix has the given shape. */
  bool hasShape(int numRows, int numCols) const
  { return this->numRows == numRows && this->numCols == numCols; }

  /** Returns true, iff B has the same shape as this matrix. */
  bool hasSameShapeAs(const rsSparseMatrix<T>& B) const
  { return hasShape(B.numRows, B.numCols); }



  bool isValidIndexPair(int i, int j) const 
  { return i >= 0 && i < numRows && j >= 0 && j < numCols; }


  /** Returns the diagonal part D of this matrix A, such that A = D + N (where N is the 
  non-diagonal part). Diagonal part means that only the diagonal elements are included, the others
  are set to zero, i.e. left out. */
  rsSparseMatrix<T> getDiagonalPart() const
  {
    rsSparseMatrix<T> D(numRows, numCols);
    for(size_t k = 0; k < elements.size(); k++)
      if(elements[k].i == elements[k].j)
        D.set(elements[k].i, elements[k].j, elements[k].value);
    return D;
  }
  // maybe this should return a vector instead of a sparse matrix, maybe make a function that 
  // splits A into diagonal and non-diagonal at once: splitOutDiagonalPart(Vec& D, Mat& N)

  /** Returns the diagonal part N of this matrix A, such that A = D + N (where D is the 
  diagonal part). */
  rsSparseMatrix<T> getNonDiagonalPart() const
  {
    rsSparseMatrix<T> N(numRows, numCols);
    for(size_t k = 0; k < elements.size(); k++)
      if(elements[k].i != elements[k].j)
        N.set(elements[k].i, elements[k].j, elements[k].value);
    return N;
  }

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Reserves storage space for the given number of elements. */
  void reserve(int numElements) { elements.reserve(numElements); }


  /** Fast insertion of an element. It just appends it to the end without any regard to checking, 
  whether or not such an element already exists (that caller should be sure, that it doesn't) and
  without taking care of inserting it at the right position. It's much faster to build a matrix 
  using this method as opposed to using set() but it must be used with care because you may end up
  with an inconsistent state of the object. To restore consistency, you may have to call 
  sortElements (not yet implemented) afterwards to make sure that it's correctly sorted because 
  otherwise, random access of elements will not work. If all you want is to compute matrix-vector
  products, you actually do not need random access, but it's nevertheless good practice to keep the
  object in a consistent state anyway...tbc... */
  void appendFastAndUnsafe(int i, int j, T val)
  {
    Element e(i, j, T(val));
    elements.push_back(e);
  }
  // todo:
  // -implement sortElements
  // -maybe rename to insertFastUnsafe or appendFastAndUnsafe
  // -add a method isConsistent() which can be used by client code to check object consistency 
  //  after using fast and usafe methods to build a matrix

  void setShape(int newNumRows, int newNumColumns)
  {
    numRows = newNumRows;
    numCols = newNumColumns;

    // Ensure that all elements have i < numRows, j < numCols. Throw away elements that are out of
    // the new range:
    for(size_t k = 0; k < elements.size(); k++) {
      if(elements[k].i >= numRows || elements[k].j >= numRows) {
        rsRemove(elements, k);
        k--;  }}
  }

  /** Sets the element at position (i,j) to the given new value. This may lead to insertion of a 
  new element (if there's no element yet at i,j) or removal of existing elements (if val is 
  zero). */
  void set(int i, int j, T val);
  // todo: 
  // -element removal needs tests
  // -maybe provide an add method that accumulates into an already existing weight instead of
  //  overwriting it

  /** Like set but instead of setting the entry to the new value, it adds the new value to what is 
  already there. This may also lead to addition or removal. */
  void add(int i, int j, T val);


  /** Sets the matrix to an all zeros matrix. This does not change the shape. */
  void setToZero() { elements.clear(); };

  // ToDo: setToIdentity, setDiagonalValues(T value), setDiagonalValues(T* values)
  // ...just like in rsMatrix

  void sortElements()
  {
    std::sort(elements.begin(), elements.end());
    //RAPT::rsHeapSort(&elements[0], getNumElements()); // linker error
  }

  /** Transposes this matrix in place. */
  void transpose()
  {
    int tmp = numRows;
    numRows = numCols;
    numCols = tmp;
    for(size_t k = 0; k < elements.size(); k++) {
      tmp = elements[k].i;
      elements[k].i = elements[k].j;
      elements[k].j = tmp;  }
    sortElements();
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Accessors. Element access via these is slow, so they should probably be only used, when 
  a matrix is built once and for all as a precomputation. When the matrix is used later e.g. in an
  iterative linear solver, you will probably want to use more efficient functions like product. */

  /** Read access. */
  T operator()(int i, int j) 
  { 
    Element e(i, j, T(0));
    if(elements.empty()) 
      return T(0);
    size_t k = (size_t) rsArrayTools::findSplitIndex(&elements[0], getNumElements(), e);
    if(k >= elements.size() || e < elements[k])
      return T(0);
    else
      return elements[k].value;
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Computations */

  /** Computes the matrix-vector product y = A*x where x must be numCols long and y must be numRows 
  long. The complexity is O(N+K) where N is the number of rows and K is the number of nonzero 
  elements. */
  template<class Tx, class Ty>
  void product(const Tx* x, Ty* y) const
  {
    rsAssert((void*)x != (void*)y, "Can't be used in place");
    for(int j = 0; j < numRows; j++)
      y[j] = Ty(0);
    for(size_t k = 0; k < elements.size(); k++)
      y[elements[k].i] += Ty(elements[k].value * x[elements[k].j]);
  }
  //void product(const T* x, T* y) const
  //{
  //  rsAssert(x != y, "Can't be used in place");
  //  for(int j = 0; j < numRows; j++)
  //    y[j] = T(0);
  //  for(size_t k = 0; k < elements.size(); k++)
  //    y[elements[k].i] += elements[k].value * x[elements[k].j];
  //}
  // -maybe include optional strideX, strideY parameters - or maybe implement a separate function 
  //  with strides
  // -how can we implement the product with the transposed matrix? would it be
  //    y[elements[k].j] += elements[k].value * x[elements[k].i];


  /** Returns a transposed version of this matrix. */
  rsSparseMatrix<T> getTranspose() const
  {
    rsSparseMatrix<T> t(*this);
    t.transpose();
    return t;
  }








  /** Computes the matrix-vector product y = A*x and returns the maximum absolute value of the 
  change in y before and after the call. This is supposed to be used inside vector iterations of
  the form yNew = A * yOld, where y takes the role of yNew and x the role of yOld. This implies 
  that x and y must have the same dimensionality which in turn means the matrix A must be square. 
  The return value can be used in a convergence test. If x and y are the same, i.e. the function is
  applied in place, the returned y vector will actually not be the matrix-vector product A*x 
  because in the computations of the values with higher index, it will use already updated values 
  with lower index. That may seem undesirable from a mathematical cleanliness point of view, but in
  practice, that may actually speed up convergence (see Jacobi- vs Gauss-Seidel iteration - it's a 
  similar situation here). So, it can be used in place in such a context - just don't expect the 
  computed y to be the exact matrix-vector product then. (todo: test, if it really improves 
  convergence - i just assume it because of the similarity to Gauss-Seidel) */
  T iterateProduct(const T* x, T* y) const;
  // may not be needed


  // todo: maybe implement matrix-operators +,-,*. Such operations should result in another 
  // (possibly less) sparse matrix. A naive algorithm for addition can use element access and set. 
  // I think, this will have complexity O(N*M*(log(K1)+log(K2))) where N,M are number of rows and 
  // columns and K1, K2 are the numbers of nozero elements in operands 1 and 2. I think, it's 
  // possible to do in O(K1 + K2) or maybe O(K1*log(K1) + K2*log(K2))...or something. Maybe the 
  // naive algo should be part of the test-suite but not the class itself.


  /** Multiplies a matrix with a std::vector to give another vector: y = A * x. */
  std::vector<T> operator*(const std::vector<T>& x) const;


  /** unary minus. Negates the matrix. */
  rsSparseMatrix<T> operator-() const;

  /** Adds two sparse matrices. */
  rsSparseMatrix<T> operator+(const rsSparseMatrix<T>& x) const;

  /** Subtracts two sparse matrices. */
  rsSparseMatrix<T> operator-(const rsSparseMatrix<T>& x) const;

  /** Multiplies two sparse matrices. */
  rsSparseMatrix<T> operator*(const rsSparseMatrix<T>& x) const;



  /** Compares two sparse matrices for equality */
  bool operator==(const rsSparseMatrix<T>& A) const
  { return numRows == A.numRows && numCols == A.numCols && elements == A.elements; }

  //-----------------------------------------------------------------------------------------------
  // maybe move into class rsIterativeLinearAlgebra - hmm - but the implementation really relies
  // on the way the data is stored...so maybe not

  /** Given the diagonal part D and non-diagonal part N of a matrix A, such that A = D + N, this 
  function solves the linear system A*x = (D+N)*x = b via Gauss-Seidel iteration and returns the 
  number of iterations taken. The resulting solution vector will be written into x. Whatever is 
  stored in the vector x before the call will be taken as initial guess. For the iteration to 
  converge, the matrix A must be strictly diagonally dominant. That means the diagonal element on 
  each row must have larger absolute value than the sum of the absolute values of the off-diagonal
  elements in the same row. 
  https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
  https://en.wikipedia.org/wiki/Diagonally_dominant_matrix   */
  static int solveGaussSeidel(const rsSparseMatrix<T>& D, const rsSparseMatrix<T>& N, T* x, 
    const T* b, T tol);

  /** Convenience function that takes the matrix A and splits it internally into diagonal and 
  non-diagonal parts. This is slow, so it should be used only in testing. In production code, the 
  splitting can typically be done once and for all as pre-processing step. */
  static int solveGaussSeidel(const rsSparseMatrix<T>& A, T* x, const T* b, T tol)
  { return solveGaussSeidel(A.getDiagonalPart(), A.getNonDiagonalPart(), x, b, tol); }
  // maybe make it non-static, to be called like A.solveGaussSeidel


  /**
  https://en.wikipedia.org/wiki/Successive_over-relaxation
  needs a workspace of size N, becomes Gauss-Seidel for w = 1 and Jacobi for w = 0.  */
  static int solveSOR(const rsSparseMatrix<T>& D, const rsSparseMatrix<T>& N, T* x, 
    const T* b, T tol, T* workspace, T w);

  static int solveSOR(const rsSparseMatrix<T>& A, T* x, const T* b, T tol, T* workspace, T w)
  { return solveSOR(A.getDiagonalPart(), A.getNonDiagonalPart(), x, b, tol, workspace, w); }


protected:

  /** Given a matrix position with row- and column-indices i,j, this function either returns the 
  flat index k at which the element is found in our elements array or, if it's not found, the index
  at which that element should be inserted. */
  /*
  int flatIndex(int i, int j)
  {

    return 0;  // preliminary
    // todo: do a binary search for element at position i,j, return its the flat index k if 
    // present, otherwise return the place where this element should be inserted
  }
  */

  struct Element
  {
    int i, j;   // row and column index
    T value;

    Element(int row, int col, T val) { i = row; j = col; value = val; }

    /** The less-than operator compares indices, i.e. a < b, iff a is supposed to be stored before 
    b in our elements array. The actual stored value plays no role in this comparison. This may 
    seem weird but it is just what is needed for the binary search (which uses the operator) that 
    is used in element insertion and random access. 
    ToDo: 
    Try to find a more elegant solution or at least make that operater inaccessible from client 
    code. Client code would probably expect a value-based comparison and be totally confused
    by this definition. Apply the "principle of least astonishment". Maybe use a binary search 
    routine that takes the comparison function as parameter such that we do not have to rely on the
    < operator when calling it. */
    bool operator<(const Element& b) const
    {
      if(i   < b.i) return true;
      if(b.i <   i) return false;
      if(j   < b.j) return true;   // i == b.i, compare with respect to j
      return false;
    }

    /** Compares two elements for equality. */
    bool operator==(const Element& b) const { return i == b.i && j == b.j && value == b.value; }

  };


  int numRows = 0, numCols = 0;  // number of rows and columns - maybe try to get rid
  
  std::vector<Element> elements;

};
// ToDo: 
// -Make another implementation (with the same interface) that stores rows. This saves one 
//  integer of storage space per entry because the row index is given implicitly. Maybe make a 
//  column-wise version, too - but that's less useful because with row-wise storage, it's more 
//  convenient and efficient to execute matrix-vector multiplications which is the most important 
//  operation in iterative linear solvers, which are the main application of sparse matrices.
// -Maybe templatize also on the integer type used for indices i,j. Users may want to use short
//  integers (16 bit) to save even more storage space, especially when T is float because then 
//  one entry is exactly 2*16+32 = 64 bits long). Maybe use TIdx, TVal for the two template 
//  parameters.
// -Make a class rsSparseTensor

template<class T>
rsSparseMatrix<T> rsSparseMatrix<T>::fromDense(const rsMatrix<T>& A, T tolRel)
{
  int M = A.getNumRows();
  int N = A.getNumColumns();
  rsSparseMatrix<T> B(M, N);
  T tol = T(0);
  T maxA = A.getAbsoluteMaximum();
  if(maxA > T(0))
    tol = tolRel / maxA;
  for(int i = 0; i < M; i++)
    for(int j = 0; j < N; j++)
      if(rsAbs(A(i, j)) > tol)   // Must use "> tol", not ">= tol" to handle tol = 0
        B.set(i, j, A(i, j));    // correctly 
  return B;
}

template<class T>
rsMatrix<T> rsSparseMatrix<T>::toDense(const rsSparseMatrix<T>& A)
{
  rsMatrix<T> B(A.getNumRows(), A.getNumColumns());
  for(int k = 0; k < A.getNumElements(); k++) {
    Element el = A.elements[k];
    B(el.i, el.j) = el.value;   }
  return B;
}

template<class T>
void rsSparseMatrix<T>::set(int i, int j, T val) 
{ 
  rsAssert(isValidIndexPair(i, j), "Index out of range");
  Element e(i, j, T(val));
  if(elements.empty() && val != T(0)) {
    elements.push_back(e);
    return;  }
  size_t k = (size_t) rsArrayTools::findSplitIndex(&elements[0], getNumElements(), e);
  if((k >= elements.size() || e < elements[k]) && val != T(0))
    rsInsert(elements, e, k);
  else {
    if(val != T(0))   // use a tolerance! Should be a member, defaulting to zero
      elements[k] = e;
    else
      rsRemove(elements, k); }
}

template<class T>
void rsSparseMatrix<T>::add(int i, int j, T val) 
{ 
  rsAssert(isValidIndexPair(i, j), "Index out of range");
  Element e(i, j, T(val));
  if(elements.empty() && val != T(0)) {
    elements.push_back(e);
    return;  }
  size_t k = (size_t) rsArrayTools::findSplitIndex(&elements[0], getNumElements(), e);
  if((k >= elements.size() || e < elements[k]) && val != T(0))
    rsInsert(elements, e, k);
  else 
  {
    e.value += elements[k].value;
    if(e.value != T(0))  // use a tolerance!
      elements[k] = e;
    else
      rsRemove(elements, k);
  }
  // There's a lot of code duplication from the set() function. Try to get rid of that! The 
  // difference is only in the "else" branch at the bottom. 
}

template<class T>
std::vector<T> rsSparseMatrix<T>::operator*(const std::vector<T>& x) const
{
  rsAssert((int) x.size() == this->numCols, "Vector incompatible for left multiply by matrix");
  std::vector<T> y(this->numRows);
  product(&x[0], &y[0]);
  return y;
}

template<class T>
rsSparseMatrix<T> rsSparseMatrix<T>::operator-() const
{
  rsSparseMatrix<T> R(*this);                  // Difference matrix, initialize as copy of this
  for(size_t k = 0; k < elements.size(); k++)
    R.elements[k].value = -R.elements[k].value;
  return R;
}
// ToDo: Factor out the loop into a negate() function that negates an existing matrix in 
// place. See implementation of rsMatrix. It also has such a function. 

template<class T>
rsSparseMatrix<T> rsSparseMatrix<T>::operator+(const rsSparseMatrix<T>& B) const
{
  rsAssert(B.hasSameShapeAs(*this), "Matrices incompatible for addition");
  rsSparseMatrix<T> S(*this);                  // Sum matrix, initialize as copy of this
  for(size_t k = 0; k < B.elements.size(); k++)
    S.add(B.elements[k].i, B.elements[k].j, B.elements[k].value);
  return S;
}

template<class T>
rsSparseMatrix<T> rsSparseMatrix<T>::operator-(const rsSparseMatrix<T>& B) const
{
  rsAssert(B.hasSameShapeAs(*this), "Matrices incompatible for subtraction");
  rsSparseMatrix<T> D(*this);                  // Difference matrix, initialize as copy of this
  for(size_t k = 0; k < B.elements.size(); k++)
    D.add(B.elements[k].i, B.elements[k].j, -B.elements[k].value);
  return D;
}

template<class T>
rsSparseMatrix<T> rsSparseMatrix<T>::operator*(const rsSparseMatrix<T>& B) const
{
  rsAssert(numCols == B.numRows, "Matrices incompatible for multiplication");
  rsSparseMatrix<T> P(numRows, B.numCols);     // Product matrix, initialize to zero matrix
  for(size_t n = 0; n < elements.size(); n++)
    for(size_t k = 0; k < B.elements.size(); k++)
      if(elements[n].j == B.elements[k].i)
        P.add(elements[n].i, B.elements[k].j, elements[n].value * B.elements[k].value);
  return P;
}

template<class T>
T rsSparseMatrix<T>::iterateProduct(const T* x, T* y) const
{
  rsAssert(numRows == numCols, "Can be used only for square matrices");
  size_t k = 0;
  T dMax = T(0);
  while(k < elements.size())
  {
    size_t i = elements[k].i;
    T yi = T(0);
    while(k < elements.size() && elements[k].i == i) {
      yi += elements[k].value * x[elements[k].j];
      k++; }
    T dyi = y[i] - yi;
    dMax = rsMax(dMax, rsAbs(dyi));
    y[i] = yi; 
  }
  return dMax;
}
// this function may not be that useful after all, maybe remove it (or put it into some sort of 
// code-attic)

template<class T>
int rsSparseMatrix<T>::solveGaussSeidel(
  const rsSparseMatrix<T>& D, const rsSparseMatrix<T>& C, T* x, const T* b, T tol)
{
  size_t N = (size_t) D.numRows;  // Why size_t and not int?
  rsAssert(D.numCols == N);       // The matrices D,C must be square and have the same shape
  rsAssert(C.numRows == N);
  rsAssert(C.numCols == N);
  int numIts = 0;
  while(true) {

    // Perform one Gauss-Seidel step and record the maximum change in the value in any of the
    // elements in vector x:
    size_t k = 0;
    T dMax = T(0);
    for(size_t i = 0; i < N; i++) {
      T xi = b[i];  // xi will become the new, updated value for x[i]
      while(k < C.elements.size() && C.elements[k].i == i)  {
        xi -= C.elements[k].value * x[C.elements[k].j];
        k++;
      }
      xi /= D.elements[i].value;
      T dxi = x[i] - xi;
      dMax = rsMax(dMax, rsAbs(dxi));
      x[i] = xi; 
    }

    // Increment iteration counter and check convergence criterion:
    numIts++;
    if(dMax <= tol)
      break;
  }

  return numIts;
}
// todo: can the internal step be expressed using C.product(x, x) - this would be nice for 
// generalizing the algo to other implementations of sparse matrices...but it may be hard to keep
// track of the dMax

template<class T>
int rsSparseMatrix<T>::solveSOR(const rsSparseMatrix<T>& D, const rsSparseMatrix<T>& C, T* x,
  const T* b, T tol, T* wrk, T w)
{
  size_t N = (size_t) D.numRows;
  int numIts = 0;
  while(true)
  {
    // Perform one iteration of SOR and record the maximum change in the value in any of the
    // elements in vector x:
    for(size_t i = 0; i < N; i++)
      wrk[i] = x[i];

    size_t k = 0;
    T dMax = T(0);
    for(size_t i = 0; i < N; i++)
    {
      T xi = b[i];  // xi will become the new, updated value for x[i]
      while(k < C.elements.size() && C.elements[k].i == i) 
      {
        T xOld = wrk[C.elements[k].j];
        T xNew = x[  C.elements[k].j];
        xi -= C.elements[k].value * (w * xNew + (T(1)-w) * xOld);
        k++;
      }
      xi /= D.elements[i].value;
      T dxi = x[i] - xi;
      dMax = rsMax(dMax, rsAbs(dxi));
      x[i] = xi;
    }

    // Increment iteration counter and check convergence criterion:
    numIts++;
    if(dMax <= tol)
      break;
  }

  return numIts;
}

//=================================================================================================

/** A class that implements iterative algorithms for numerical linear algebra. It is written in
such a way that the same code can be used with dense or sparse matrices or any other kind of 
special matrix, as long as some appropriate low-level functions are added to this class for the
specific matrix-types. These low-level functions include things like retrieving the shape and for
computing matrix-vector or matrix-matrix products. These low-level special purpose implementations 
(one for each special matrix class) will then be used inside the actual computational algorithms 
which themselves will be the same for all the different matrix classes. If you want to use the 
algorithms for some new special matrix class, you will need to add just a small amount of 
boilerplate code for these operations - either as static member functions here in the class (if 
you are me) or as free functions (if you are a client and don't want to hack my library). */

class rsIterativeLinearAlgebra
{

public:



  //-----------------------------------------------------------------------------------------------
  // \name Low Level Interface


  /** Solves the linear system of equations A*x = b by means of the conjugate gradient (CG) method.
  In its original form, the CG method converges only for symmetric positive definite (SPD) matrices 
  A. However, if the "leastSquares" flag is true, the function applies the CG method to the 
  modified linear system A^T*A*x = A^T*b in which case the matrix A^T*A is indeed symmetric and at 
  least positive semi-definite (todo: figure out, if that is enough to ensure convergence). This 
  modified system can be obtained by either simply pre-multiplying the original system by A^T or by
  considering the least-squares problem dot(r,r) = min where r = b - A*x is the residual. I call 
  this the least squares conjugate gradient (LSCG) method but this is not standard terminology, as 
  far as i know. Theoretically, if N is the size of the matrix (which must be square), the 
  algorithm should converge in at most N steps. However, due to roundoff error, that may not be the 
  case and the numerical problems get worse when using LSCG instead of regular CG (because the 
  condition number of the matrix gets squared for LSCG). For LSCG, i.e. if leastSquares is true, 
  the workspace must be of size 4*N. For regular CG, a size of 3*N is enough. To make the 
  implementation suitable also for sparse matrices, the A^T * A matrix is not explicitly created
  (it may not be sparse even if A is sparse). Instead there will be one call to the regular product
  function and another to the transProduct function per iteration (in case of LSCG only). With the 
  optional shift parameter, you can have the algorithm work with a shifted matrix A' = A + shift*I 
  where I is the identity matrix. @see shift() */
  template<class T, class TMat>
  static int solveViaCG(const TMat& A, T* x, const T* b, T* workspace, T tol, int maxIts, 
    bool leastSquares = false, const T& shift = T(0));





  template<class T, class TMat>
  static int largestEigenValueAndVector(const TMat& A, T* val, T* vec, T tol, T* workspace);
  // maybe rename to vonMisesIteration or eigenViaVonMises/eigenViaPowerIteration
  // largestEigenpair

  template<class T, class TMat>
  static int eigenspace(const TMat& A, T* vals, T* vecs, T tol, T* workspace, int maxIts);
  // rename to eigenViaPower or eigensystemViaEPI (extended power iteration)
  // each eigenvector is found in turn from the largest to the smallest via a variation of the von 
  // Mises iteration in which the projection of the iterates onto the already found eigenspace is
  // subtracted from the iterates
  // returns the maximum number of iterations, i.e. the number of iteration for the eigenvector 
  // that took the most iterations...hmm - or maybe it should return the sum of all iterations
  // -maybe it should take an additional parameter to specify, how many eigenpairs should be found
  //  and also a maximum number of iterations
  // mayb rename to eigensystem
  // https://reference.wolfram.com/language/ref/Eigensystem.html
  // sizes: A: NxN, vals: N, vecs: N*N, workspace: N



  /** Returns true, iff the vector y is a scalar multiple of the vector x (up to some tolerance). 
  Both vectors must be of length N. The "factor" parameter will get the scale factor assigned if y
  is indeed a multiple of x or 0, if not. I deliberately chose 0 and not nan for this case for two
  reasons: (1) the function should be able to work with number types T that have no concept of nan 
  (fractions, modular integers, etc. (2) the function is mainly intended for dealing with 
  eigenvectors and they are defined to be nonzero, so y cannot be 0*x and still count as 
  eigenvector, so the zero is actually "free" and can be used to encode that condition. */
  template<class T>
  static bool isScalarMultiple(const T* x, const T* y, int N, T tol, T* factor);
  // todo: (decide and) document if tolerance is absolute or relative, maybe rename to tolR, if
  // relative
  // maybe rename to rsIsMultiple move into rsArrayTools (and maybe have an alias isScalarMultiple
  // here)


  // todo: implement:
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method
  // https://en.wikipedia.org/wiki/Biconjugate_gradient_method
  // https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method

protected:

  /** If leastSquares is false, this function just computes the regular matrix-vector product 
  y = A * x. If leastSquares is true, it computes the product y = A^T * A * x that occurs in the 
  least squares problem |A*x-b|^2 = min that is related to the problem A*x-b = 0. Only in the 
  second case the workspace is needed (it then must have the same size as x and y). */
  template<class T, class TMat>
  static void productLS(const TMat& A, const T* x, T* y, 
    bool leastSquares = false, T* workspace = nullptr, const T& shift = T(0));

  /** This is used to modify the result of a matrix-vector multiplication y = A*x by adding to each
  y[i] and additional s*x[i]. This has the same effect as if the matrix A would have been modified
  to A' = A + s*I where I is the identity matrix. The matrix A' has the same eigenvectors as A but
  all eigenvalues are shifted by s. Such shifts play a role in computations of eigenvectors: if 
  you know some eigenvalue s_k, to find the corresponding eigenvector v_k, you need to solve the 
  linear system (A - s_k*I) * v_k = 0 (verify this!). */
  template<class T>
  static void shift(const T& s, const T* x, T* y, int N)
  {
    if(s != T(0))
      for(int i = 0; i < N; i++)
        y[i] += s * x[i];
  }
  // needs tests



  // Specializations of some low-level functions (this is boilerplate):

  // Specializations for rsMatrix:
  template<class T> static int numRows(const rsMatrix<T>& A) { return A.getNumRows(); }
  template<class T> static int numColumns(const rsMatrix<T>& A) { return A.getNumColumns(); }
  template<class T> static void product(const rsMatrix<T>& A, const T* x, T* y) { A.product(x, y); }
  template<class T> static void transProduct(const rsMatrix<T>& A, const T* x, T* y) { A.transProduct(x, y); }

  // Specializations for rsSparseMatrix:
  template<class T> static int numRows(const rsSparseMatrix<T>& A) { return A.getNumRows(); }
  template<class T> static int numColumns(const rsSparseMatrix<T>& A) { return A.getNumColumns(); }
  template<class T> static void product(const rsSparseMatrix<T>& A, const T* x, T* y) { A.product(x, y); }


  template<class TMat>
  static bool isSquare(const TMat& A) { return numRows(A) == numColumns(A); }

};

template<class T, class TMat>
void rsIterativeLinearAlgebra::productLS(
  const TMat& A, const T* x, T* y, bool ls, T* w, const T& s)
{
  rsAssert(isSquare(A), "Matrix A must be square");
  int N = numRows(A);
  if(ls) {
    product(     A, x, w   );
    shift(       s, x, w, N);
    transProduct(A, w, y   ); 
    shift(       s, w, y, N); }
  else {
    product(     A, x, y   );
    shift(       s, x, y, N); }
}

template<class T, class TMat>
int rsIterativeLinearAlgebra::solveViaCG(const TMat& A, T* x, const T* b, 
  T* wrk, T tol, int maxIts, bool ls, const T& s)
{
  rsAssert(isSquare(A), "Matrix A must be square");
  using AT = rsArrayTools;            // shortcut for convenience
  int N = numRows(A);                 // number of inputs and outputs, size of x and b
  T* r = &wrk[0];                     // residual
  T* p = &wrk[N];                     // current direction for update
  T* t = &wrk[2*N];                   // temporary for product A*p
  T* w = nullptr;                     // workspace for least-squares product
  T rho0;
  if(ls) {                            // LSCG algorithm was requested
    w = &wrk[3*N];                    // additional workspace of size N required
    productLS(   A, x, t, ls, w, s);  // t = A^T * A * x (with optional shift)
    transProduct(A, b, p);            // p = A^T * b
    shift(       s, b, p, N);         // optional shift
    AT::subtract(p, t, r, N);         // r = p-t = A^T * b - A^T * A * x
    rho0 = AT::sumOfSquares(p, N);    // rho0 = dot(p, p)
  }
  else {                              // regular CG algorithm was requested
    product(A, x, t);                 // t = A*x
    shift(  s, x, t, N);              // optional shift
    AT::subtract(b, t, r, N);         // r = b-t = b - A*x
    rho0 = AT::sumOfSquares(b, N);    // rho0 = dot(b, b)
  }
  AT::copy(r, p, N);                  // p = r
  T rho = AT::sumOfSquares(r, N);     // rho = dot(r, r)
  T a;                                // current step size for update 
  T rhoP;                             // rho of previous iteration
  T rhoR;                             // ratio rho/rhoP
  T rhoC = tol*tol*rho0;              // convergence test threshold for rho
  int k;
  for(k = 0; k < maxIts; k++) {
    if(rho <= rhoC) return k;             // convergence test
    productLS(A, p, t, ls, w, s);         // t = A*p (CG) or t = A^T*A*p (LSCG) with optional shift
    a = rho / AT::sumOfProducts(p, t, N); // a = rho / dot(p,t)
    for(int i = 0; i < N; i++) {
      x[i] += a*p[i];                     // x = x + a*p
      r[i] -= a*t[i]; }                   // r = r - a*t
    rhoP = rho;
    rho  = AT::sumOfSquares(r, N);        // rho = dot(r,r)
    rhoR = rho / rhoP;
    for(int i = 0; i < N; i++)
      p[i] = r[i] + rhoR*p[i]; }
  return k;
}
// References:
// -https://en.wikipedia.org/wiki/Conjugate_gradient_method
// -Conjugate gradient type methods for unsymmetric and inconsistent systems of linear equations
//
// ToDo:
// -Introduce scalar shift parameter. After each product of the form y = A*x, do y += shift*x 
//  (maybe conditionally: if(shift != 0)...). This has the same effect as adding the shift 
//  parameter to each diagonal element of the matrix A and is important in the context of 
//  eigenvalues. It shifts the spectrum of the matrix A. In productLS it must be applied after 
//  product *and* after transProduct (i think). It should call the shift member function
// -needs tests with singular systems (consistent and inconsistent)
// -try to incorporate preconditioning...but maybe that should be done in pre/post processing steps
//  outside this function
//
// oh - the extension to the least squares problem has already a name: CGNR (Conjugate Gradient for 
// Normal Residual):
// https://en.wikipedia.org/wiki/Conjugate_gradient_method#Conjugate_gradient_on_the_normal_equations
// see also: https://web.stanford.edu/group/SOL/software/lsqr/

template<class T, class TMat>
int rsIterativeLinearAlgebra::largestEigenValueAndVector(
  const TMat& A, T* val, T* vec, T tol, T* wrk)
{
  rsAssert(numRows(A) == numColumns(A), "Can be used only for square matrices");
  using AT = rsArrayTools;
  int N = numRows(A);
  T L = AT::euclideanNorm(vec, N);
  AT::scale(vec, N, T(1) / L);
  int numIts = 0;
  while(true) {
    product(A, vec, wrk);
    L = AT::euclideanNorm(wrk, N);
    AT::scale(wrk, N, T(1) / L);
    T dMax = AT::maxDeviation(vec, wrk, N);
    if(dMax <= tol) {
      *val = L;                          // that's only the absolute value...
      int i = AT::maxAbsIndex(vec, N);
      if(vec[i] * wrk[i] < T(0))         // ...this figures out the sign
        *val = -(*val);
      AT::copy(wrk, vec, N);             // return the last iterate, not the 2nd to last
      break; }
    AT::copy(wrk, vec, N);
    numIts++; }
  return numIts;
}
// This implements the von Mises vector iteration.
// https://en.wikipedia.org/wiki/Power_iteration

// ToDo:
// -try using the maximum norm instead of the Euclidean - may be numerically more precise due to 
//  involving less arithmetic (none, actually)
// -maybe rename to MisesIteration
// -maybe include a maxNumIterations parameter
// -maybe factor out the stuff below product, such that the compiled version of this code can be 
//  re-used for dense matrices (to reduce binary size)..maybe into a normalizeAndTest function
//  that returns a boolean, the copy inside the if can be dragged outside the while(true) loop
// -try to come up with an in-place iteration like Gauss-Seidel - but it must take care to not
//  use the updated values as is but scale them by the current stimate of the eigenvalue, because 
//  the next iterate is expected to be scaled by the eigenvalue in the iteration - so the summation
//  loop need to be split inot two halfs: one with the scaling (using the updated values) and one 
//  without (using the old values)
// -try to find all eigenvalues and -vectors by subtracting the projection onto the already found
//  eigenspace from each iterate - this should be done right after product(vec, wrk); i think

template<class T, class TMat>
int rsIterativeLinearAlgebra::eigenspace(
  const TMat& A, T* vals, T* vecs, T tol, T* wrk, int maxIts)
{
  // Algorithm:
  // It works like the von Mises vector iteration for one eigenvector at the time, from large to 
  // small. To avoid converging to the same eigenvector again, instead of just forming the 
  // matrix-vector product of the matrix with the previous iterate, we subtract from that product 
  // the projection of it onto the space spanned by the already previously found eigenvectors.

  rsAssert(numRows(A) == numColumns(A), "Can be used only for square matrices");
  using AT = rsArrayTools;

  T (*norm)(const T*, int) = &AT::euclideanNorm;
  //T (*norm)(const T*, int) = &AT::maxAbs;
  // i think, any norm should work, but maxAbs apparently doesn't - why?

  int N = numRows(A);
  int totalIts = 0;
  for(int n = 0; n < N; n++) {         // loop over the eigenvectors
    T* val = &vals[n];                 // location of n-th eigenvalue
    T* vec = &vecs[n*N];               // start location of n-th eigenvector
    int itsPerVec = 0;
    while(true) {
      itsPerVec++; 
      product(A, vec, wrk);
      for(int i = 0; i < n; i++) {
        T pi = T(0);
        for(int j = 0; j < N; j++) 
          pi += wrk[j] * vecs[i*N + j];   // compute projection coeff
        for(int j = 0; j < N; j++) 
          wrk[j] -= pi * vecs[i*N + j]; } // subtract projection
      bool done = isScalarMultiple(vec, wrk, N, tol, val);
      T L = norm(wrk, N);
      AT::scale(wrk, N, T(1) / L);        // what if L == 0?
      AT::copy(wrk, vec, N);
      if(done || itsPerVec >= maxIts) {
        totalIts += itsPerVec;
        break;  }
    }
  }
  return totalIts;

  // todo: drag the T L = norm(vec, N); AT::scale(vec, N, T(1) / L); into the iteration before
  // forming the product, remove them after the convergence test
  // 
}
// -needs a maxIts parameters
// -maybe get rid of the function that that computes only the largest - give this function another
//  parameter that determines, how many eigenvalues should be computed, and if it's 1, it just 
//  reduces to the function that computes the largest.
// -maybe have a boolean parameter to switch between finding eigenspace of A or A^-1 - in the 2nd
//  case, we may have to use solve(A, vec, wrk) or solve(A, wrk, vec) instead of 
//  product(A, vec, wrk) ...maybe we could also find eigenvalues of A^T, A^-T
// -Maybe try to improve convergence by using a matrix with suitably shifted eigenvalues. The 
//  shift should be such as to maximize the ratio between the largest and second-to-largest 
//  (remaining, absolute) eigenvalue. Maybe the trace could help: it's the sum of all eigenvalues, 
//  so we may divide by N to get the average and subtract that average to center the eigenvalues 
//  around zero. But for the 2nd eigenpair, we want to center the *remaining* eigenvalues around 
//  zero. Maybe subtract all the already found eigenvalues from the trace before dividing by N 
//  -> experiments needed. ...the determinant is the product of all eigenvalues btw - but i don't 
//  think that's helpful here (it's hard to compute anyway). For example, if we have a matrix with
//  eigenvalues 1,2,3, the convergence speed will be given by |3|/|2| = 1.5. That's the factor by 
//  which the contribution of the 3rd eigenvector grows faster than that of the 2nd. If we shift 
//  by -1, the shifted eigenvalues are 0,1,2 and the convergence speed will be |2|/|1| = 2 which is
//  better. But we don't know the ratio of largest to second largest which (i think) is needed to 
//  make an informed shift. Maybe we could estimate it in each iteration by computing in each 
//  iteration 2 vectors: one as usual and one with the projection onto the usual vector subtracted,
//  which is going to be our estimate for the 2nd largest eigenvector ...dunno, just brainstorming

template<class T>
bool rsIterativeLinearAlgebra::isScalarMultiple(const T* x, const T* y, int N, T tol, T* factor)
{
  *factor = T(0);

  // Skip initial section of zeros:
  auto isZero = [&](const T& number) { return rsAbs(number) <= tol; };
  int i = 0;
  while(isZero(x[i]) && isZero(y[i]))  
    i++;
  if(i == N || isZero(x[i]) || isZero(y[i])) 
    return false;
    // either both x and y are all zeros or we found an i for which x[i] is zero but y[i] is 
    // nonzero or vice versa

  // If we reach this point, we have found an index i for which both x[i] and y[i] are nonzero. We 
  // compute their ratio, defined by r := y[i] / x[i] such that y[i] = r * x[i]:
  T r = y[i] / x[i];

  // Now we iterate through the remaining array and check for each pair x[i], y[i], if they either
  // obey the same ratio or are both zero:
  i++;
  while(i < N) {
    // Check, if the ratio of x[i] and y[i] matches r unless both are zero:
    if(!(isZero(x[i]) && isZero(y[i]))) {
      T yt = r  * x[i];     // expected target value for y
      T d  = yt - y[i];     // difference between target and actual
      if(rsAbs(d) > tol)
        return false;    }  // the ratio of x[i] and y[i] did not match r
    i++; }

  // If we reach this point, y is indeed a scalar multiple of x and the scale factor is r:
  *factor = r;
  return true;
}





template<class T>
T rsDot(const std::vector<T>& x, const std::vector<T>& y)
{
  rsAssert(x.size() == y.size());
  T sum = T(0);
  for(size_t i = 0; i < x.size(); i++)
    sum += x[i] * y[i];
  return sum;
}
// move elsewhere

template<class T>
bool rsStaysFixed(const std::vector<T>& x, const std::vector<T>& dx, T s = T(1), T tolR = T(0))
{
  rsAssert(x.size() == dx.size());
  return rsArrayTools::staysFixed(&x[0], &dx[0], (int)x.size(), s, tolR);
}

/** Solves the linear system A * x = b iteratively via the conjugate gradient (CG) method. The 
matrix A must be symmetric and positive definite (SPD) for the algorithm to converge. */
template<class T>
int rsSolveCG(const rsMatrix<T>& A, std::vector<T>& x, const std::vector<T>& b, T tol, int maxIts)
{
  using Vec = std::vector<T>;
  Vec r = b - A*x;         // residuum?
  Vec p = r;               // current direction?
  T rho0 = rsDot(b, b);
  T rho  = rsDot(r, r);
  Vec t;
  T a, rhos;
  int k;
  for(k = 0; k < maxIts; k++)
  {
    if(sqrt(rho/rho0) <= tol)  // avoid sqrt: use rho/rho0 <= tol*tol
      return k;                
    // ToDo: we need a better stopping criterion. Maybe something based on 
    // rsStaysFixed(x, p, a, tol) ...

    t    = A*p;
    a    = rho / rsDot(p,t);
    x    = x + a*p;
    r    = r - a*t;
    rhos = rho;
    rho  = rsDot(r, r);
    p    = r + (rho/rhos)*p;
  }
  return k;
  //return maxIts+1;
}
// References: Numerical Linear Algebra and Matrix Factorizations (Tom Lyche), pg 286

/** Solves the linear system A * x = b iteratively via applying the conjugate gradient method to 
the related problem: A^T * A * x = A^T * b. This related problem can be derived by either simply
pre-multiplying both sides by A^T or by minimizing dot(A*x-b, A*x-b) = (A*x-b)^T * (A*x-b). The 
nice thing is that the matrix A^T * A is guarenteed to be positive (semi?)definite (which is 
required by CG for convergence), so this strategy extends the applicability of the conjugate 
gradient method to (almmost?) general matrices. I call this the least squares conjugate gradient
(LSCG) method, because it can be seen as minimizing the squared norm of the residual r = A*x-b. I 
think, if the problem has no solution, the algo will produce a least squares approximation to a 
solution and if it has multiple solutions, it will produce the minimum norm solution (tests 
needed!). */
template<class T>
int rsSolveLSCG(const rsMatrix<T>& A, std::vector<T>& x, const std::vector<T>& b,
  T tol, int maxIts)
{
  rsMatrix<T> AT = A.getTranspose();
  return rsSolveCG(AT*A, x, AT*b, tol, maxIts);
}
// -more tests needed, especially with singular (consistent and inconsistent) systems
// -figure out and document, if there can be systems, where this method breaks down
// -an efficient implementation that works also for sparse matrices is needed - this should not
//  explicitly create the matrix A^T * A, because it may not be sparse even if A is sparse

template<class T>
int rsSolveShiftedLSCG(const rsMatrix<T>& A, std::vector<T>& x, const std::vector<T>& b,
  T tol, int maxIts, T shift, bool leastSquares = true)
{
  rsAssert(A.isSquare());
  int N = A.getNumRows();
  using Mat = rsMatrix<T>;
  Mat As = A + Mat::diag(N, shift);
  if(leastSquares) return rsSolveLSCG(As, x, b, tol, maxIts);
  else             return rsSolveCG(  As, x, b, tol, maxIts);
}

template<class T>
int rsSolveRichardson(const rsMatrix<T>& A, std::vector<T>& x, const std::vector<T>& b, T alpha,
  T tolR, int maxIts)
{
  int its = 0;
  std::vector<T> dx;
  while(its < maxIts) {
    dx = alpha * (b - A*x);
    if(rsStaysFixed(x, dx, T(1), tolR)) 
      return its;  // x has converged
    x = x + dx;
    its++; }
  return its;
}
// References: Numerical Linear Algebra and Matrix Factorizations (Tom Lyche), pg 261
// -if A has positive eigenvalues s1 >= s2 >= sN > 0, the method converges, iff 0 < alpha < 2/s1.
// -The iteration can also be expressed as xNew = (I-alpha*A)*xOld + b. Define kappa = s1/sN, 
//  alpha_0 = 2 / (s1 + sN), then we have:
//    rho(I-alpha*A) > rho(I-alpha_0*A) = (kappa-1)/(kappa+1). 
//  where rho is the spectral radius (i guess)
// -I think, this algorithm corresponds to gradient descent algorithm in numerical optimization. 
//  Starting from that as the baseline, we can try various tweaks to try to improve the 
//  convergence speed and robustness, such as:
//  -select update rate alpha per iteration (...somehow)
//  -use momentum
//  -use an alpha vector, i.e. an update rate per coordinate, possibly adaptive
//  -use (adaptive) momentum per coordinate
//  -maybe we could detect if the update vector alternates sign between the iterations for a given
//   coordinate, and if so, increase the momentum and/or decrease update rate for that coordinate
//   ...and decrease/increase it, if it doesn't alternate?
//  -or: use a 2-point moving avergae for the dx - that should effectively remove sign 
//   alternations in the updates. maybe the filter coeffs (which should sum to 1) can also be made
//   adaptive. maybe based on the difference of successive values: if the difference is large, the 
//   coeff should be large (i.e. smoothing should be stronger) and if the difference is zero, the 
//   coeff should be zero. maybe based on some function f(x) that is zero for x=0 and approaches
//   0.5, if |x| approaches infinity
//  -detect, if the error err = b - A*x has increased with respect to previous iteration, if so,
//   maybe undo the step and also adapt momentum and update rate
// -maybe interpret the problem as minimization of dot(A*x-b, A*x-b) = (A*x-b)^T * (A*x-b)
//    = ((A*x)^T-b^T) * (A*x-b) 
//    = (A*x)^T * A*x  -  (A*x)^T * b  -  b^T * A*x  +  b^T * b
//    = x^T*A^T*A*x  -  x^T*A^T*b  -  b^T*A*x  +  b^T*b = min
//  define:
//    P := 2*(A^T*A), q = (A^T + A)*b
//  then solve:
//    P*x - q = 0   or   P*x = q


// see also:
// https://en.wikipedia.org/wiki/Generalized_minimal_residual_method (has matlab code)


//=================================================================================================

template<class T> 
bool isHarmonic2(const rsBivariatePolynomial<T>& u, T tol = T(0))
{
  for(int m = 2; m <= u.getDegreeX(); m++) {
    for(int n = 2; n <= u.getDegreeY(); n++) {
      T c_xx = m * (m-1) * u.coeff(m  , n-2);   // coeff for x^(m-2) * y^(n-2) in u_xx
      T c_yy = n * (n-1) * u.coeff(m-2, n  );   // coeff for x^(m-2) * y^(n-2) in u_yy
      if(rsAbs(c_xx + c_yy) > tol)
        return false;   }}
  return true;
}
// -needs more tests
// -if it works, maybe replace the old implementation in rsBivariatePolynomial - this here is more
//  efficient
// -make a function makeHarmonicY that assigns:
//    u.coeff(m-2, n) = -m * (m-1) * u.coeff(m, n-2) / (n * (n-1))
//  and a makeHarmonicX that assigns:
//    u.coeff(m, n-2) = -n * (n-1) * u.coeff(m-2, n) / (m * (m-1))
//  and maybe some sort of symmetric version that assigns both - the goal is always that 
//  c_xx + c_yy = 0  ->  m * (m-1) * u[m, n-2] + n * (n-1) * u[m-2, n] = 0
//  m * (m-1) * u[m, n-2] = -n * (n-1) * u[m-2, n]
//  so:
//    u[m, n-2] / u[m-2, n] = -(n*(n-1)) / (m*(m-1))
//  i think, this is the general symmetry condition that the matrix for a harmonic polynomial must
//  satisfy

template<class T> 
rsPolynomial<std::complex<T>> getComplexFromHarmonicUV(
  const rsBivariatePolynomial<T>& u, const rsBivariatePolynomial<T>& v)
{
  rsAssert(rsBivariatePolynomial<T>::areHarmonicConjugates(u, v));
  // fails: seems like ux == -vy where we should have ux == vy - there's some sign confusion: i 
  // think P_y is not the harmonic conjugate of of P_x, -P_y is

  /*
  int M = rsMax(u.getDegreeX(), v.getDegreeX());
  rsPolynomial<std::complex<T>> p(M);
  for(int i = 0; i <= M; i++)
  p[i] = std::complex<T>(u.getCoeffPadded(i, 0), v.getCoeffPadded(i, 0));
  return p;
  // seems like if v.degX > u.degX, the last coeff is zero anyway - more tests needed
  */

  int M = rsMin(u.getDegreeX(), v.getDegreeX());
  rsPolynomial<std::complex<T>> p(M);
  for(int i = 0; i <= M; i++)
    p[i] = std::complex<T>(u.coeff(i, 0), v.coeff(i, 0));
  return p;
}
// needs more tests - especially: why should the loop go only up to M and not max(M,N), where 
// N = v.getDegreeX(), using a zero-padded accessor like getCoeffPadded(i, 0) for safe 
// out-of-range access returning 0? ...or maybe the loop should go only up to min(M,N)


template<class T> 
void firstFundamentalForm(const rsBivariatePolynomial<T>& x, const rsBivariatePolynomial<T>& y,
  const rsBivariatePolynomial<T>& z, rsBivariatePolynomial<T>& E, rsBivariatePolynomial<T>& F,
  rsBivariatePolynomial<T>& G)
{
  using BiPoly = rsBivariatePolynomial<T>;
  BiPoly x_u = x.derivativeX();     // dx/du
  BiPoly x_v = x.derivativeY();     // dx/dv
  BiPoly y_u = y.derivativeX();     // dy/du
  BiPoly y_v = y.derivativeY();     // dy/dv
  BiPoly z_u = z.derivativeX();     // dz/du
  BiPoly z_v = z.derivativeY();     // dz/dv
  E = x_u*x_u + y_u*y_u + z_u*z_u;
  F = x_u*x_v + y_u*y_v + z_u*z_v;
  G = x_v*x_v + y_v*y_v + z_v*z_v;
}
// see Weitz - Differentialgeometrie, p.154
// -maybe make a function that takes as input an rsVector3D of rsBivariatePolynomial instead of 
//  x,y,z separately
// -are E,F,G the entries of J^T * J where J is the Jacobian and E*G - 2*F is its determinant?

// todo: 
// -consider parametric surfaces given by a triple of bivariate polynomials:
//  x(u,v), y(u,v), z(u,v) (note that u,v are here the independent variables, i.e. the inputs to
//  the 3 component polynomials, in the context of complex analysis, they were used for the real 
//  and imaginary part of the function - don't get confused by this!)
// -write a function to compute the first fundamental form coeffs E,F,G of the surface
// -functions for principal curvatures, mean curvature and Gaussian curvature - these should
//  all be bivariate polynomials again (if i'm not mistaken)
// -circulation ond flux through a curve (path integrals over curl and divergence)
// -maybe integral functions that take univariate polynomials for the integration limits


// experimental - doesn't seem to work:
template<class T> 
void vectorPotential2(const rsTrivariatePolynomial<T>& f, const rsTrivariatePolynomial<T>& g,
  const rsTrivariatePolynomial<T>& h, const rsTrivariatePolynomial<T>& d, 
  rsTrivariatePolynomial<T>& F, rsTrivariatePolynomial<T>& G, rsTrivariatePolynomial<T>& H)
{
  using TP = rsTrivariatePolynomial<T>;
  TP fz   = f.integralZ();
  TP gz   = g.integralZ();
  TP fz_x = fz.derivativeX();
  TP gz_y = gz.derivativeY();
  TP a_x  = h + fz_x + gz_y;
  TP a    = a_x.integralX();

  TP fz_y = fz.derivativeY();
  TP gz_x = gz.derivativeX();
  TP a_y  = a.derivativeY();
  TP D    = d + fz_y - gz_x;
  TP b_x  = D - a_y;
  TP b    = b_x.integralX();

  F = b + gz;       // F(x,y,z) =  gz + b(x,y)
  G = a - fz;       // G(x,y,z) = -fz + a(x,y)
  H = TP(0,0,0);    // H(x,y,z) =  0
}


template<class T>
void potentialToDivergence(const rsBivariatePolynomial<T>& P, rsBivariatePolynomial<T>& D)
{
  int M = P.getDegreeX(); 
  int N = P.getDegreeY();
  D.initialize(M, N);
  for(int m = 0; m <= M; m++) {
    for(int n = 0; n <= N; n++) {
      if(m <= M-2 && n <= N-2)
        D.coeff(m, n) = (m+1)*(m+2)*P.coeff(m+2, n) + (n+1)*(n+2)*P.coeff(m, n+2);
      else if(m <= M-2 && n > N-2)
        D.coeff(m, n) = (m+1)*(m+2)*P.coeff(m+2, n);
      else if(m > M-2 && n <= N-2)
        D.coeff(m, n) = (n+1)*(n+2)*P.coeff(m, n+2);
      else
        D.coeff(m, n) = T(0); }}
}

template<class T>
void divergenceToPotential1(const rsBivariatePolynomial<T>& D, rsBivariatePolynomial<T>& P)
{
  int M = D.getDegreeX();
  int N = D.getDegreeY();
  P.initialize(M+2, N+2);
  for(int m = 0; m <= M; m++) 
    for(int n = 0; n <= N; n++) 
      P.coeff(m+2, n) = (D.coeff(m, n) - (n+1)*(n+2)*P.coeff(m, n+2)) / ((m+1)*(m+2));
}
// something is still wrong - now it produces too many nonzero entries - i think, it works only 
// when the bottom-right 2x2 section of the divergence is zero?
// seems like we could fix it by zeroing the last two rows ..or well - we need to do something
// that creates zeros in the last two rows in the computed divergence

template<class T>
void divergenceToPotential4(const rsBivariatePolynomial<T>& D, rsBivariatePolynomial<T>& P)
{
  int M = D.getDegreeX();
  int N = D.getDegreeY();
  P.initialize(M+2, N+2);
  for(int m = 0; m <= M; m++) 
    for(int n = 0; n <= N; n++) 
      P.coeff(m+2, n) = (D.coeff(m, n) - (n+1)*(n+2)*P.coeff(m, n+2)) / ((m+1)*(m+2));

  //for(int m = M; m <= M+2; m++) 
  //  for(int n = 0; n <= N+2; n++)
  //    P.coeff(m, n) = T(0);
}


template<class T>
void divergenceToPotential2(const rsBivariatePolynomial<T>& D, rsBivariatePolynomial<T>& P)
{
  int M = D.getDegreeX();
  int N = D.getDegreeY();
  P.initialize(M+2, N+2);
  for(int n = 0; n <= N; n++)
    for(int m = 0; m <= M; m++)
      P.coeff(m, n+2) = (D.coeff(m, n) - (m+1)*(m+2)*P.coeff(m+2, n)) / ((n+1)*(n+2));
}
// note that we had to exchange the loops compared to divergenceToPotential1 

// experimental:
template<class T>
void divergenceToPotential3(const rsBivariatePolynomial<T>& D, rsBivariatePolynomial<T>& P)
{
  int M = D.getDegreeX();
  int N = D.getDegreeY();
  P.initialize(M+2, N+2);
  for(int n = 0; n <= N; n++)
    for(int m = 0; m <= M; m++)
      if(m <= M-2 /* || n <= N-2*/)
        P.coeff(m, n+2) = (D.coeff(m, n) - (m+1)*(m+2)*P.coeff(m+2, n)) / ((n+1)*(n+2));
}

//=================================================================================================

/** A class for representing modular integers with respect to the modulus 4179340454199820289 which
is suitable for radix-2 number theoretic transforms using 64 bit unsigned integers. 

...nope! That doesn't work! We get immediately overflow when multiplying a root with itself. I 
think, the modulus must be less than 2^32 to avoid this. 2^31 - 1 actually happens to be a Mersenne
prime ...maybe that should be used? Maybe that can simplify the modulo operation? But no, for this, 
not even a 4th root of unity exists. Using a radix of 3, we only get roots upt to 3^2 = 9.

Rename this to rsModularIntegerNTT_128 and replace rsUint64 by rsUint128...once we have such a 
thing. I think, this may work with 128 bit wide integers. Make a class rsModularIntegerNTT_64 using
3221225473 as modulus.  ...done...and this does indeed work and supports radix-2 NTTs up to length
N = 2^30 ....that should be more than enough in practice */

//typedef signed int rsInt32;

//typedef __uint128 rsUint128;  // gcc?
//typedef __int128_t rsUint128;   // msc
//typedef uint128_t rsUint128;   // msc


class rsModularIntegerNTT // rename to rsModularIntegerNTT_128
{

public:

  using ModInt = rsModularIntegerNTT;

  rsUint64 value;  // i think, we need a 128 bit wide integer type for this modulus

  rsModularIntegerNTT() {}
  rsModularIntegerNTT(rsUint64 x) : value(x) {}

  ModInt operator+(const ModInt& b)
  {
    return (value + b.value) % modulus;
  }

  ModInt operator*(const ModInt& b)
  {
    return (value * b.value) % modulus;
  }

  ModInt& operator+=(const ModInt& b) { *this = *this+b; return *this; }
  ModInt& operator*=(const ModInt& b) { *this = *this*b; return *this; }


  bool operator==(const ModInt& b) const { return value == b.value; }
  bool operator!=(const ModInt& b) const { return value != b.value; }


  // The magic numbers (definitions of the arrays are in .cpp file):
  static const rsUint64 modulus = 4179340454199820289;
  static const rsUint64 roots[15];      // N-th roots of unity for N = 2^(k+1), k is array index
  static const rsUint64 rootsInv[15];   // modular inverses of the roots
  static const rsUint64 lengthsInv[15]; // modular inverses of the lengths N

};

//-------------------------------------------------------------------------------------------------

class rsModularIntegerNTT_64
{

public:

  using ModInt = rsModularIntegerNTT_64;

  rsUint64 value;

  rsModularIntegerNTT_64() {}
  rsModularIntegerNTT_64(rsUint64 x) : value(x) {}

  ModInt operator+(const ModInt& b) { return (value + b.value) % modulus; }
  ModInt operator-(const ModInt& b) 
  { 
    if( b.value > value ) return modulus + value - b.value;
    else                  return           value - b.value;
    // Can we do this in a branchless way? If not, implement unary minus and express subtraction as 
    // addition of negative....maybe we can just do: 
    //   return (modulus + value - b.value) % modulus
    // that would indeed be branchless but again involves the (supposedly expensive) modulo 
    // operator...but we really need to measure, what is expensive and what is not
  }
  // Maybe instead of explicit mod, do something like subtracting the modulus if the result is 
  // larger that the modulus in a branch-free way, like:
  // rsUint64 tmp = value + b.value;
  // rsUint64 g   = tmp >= modulus;
  // return (1-g)*tmp + g*(tmp-modulus);
  // this replaces the % by a cmp, 2*, 2-, 1+, 1=
  // we have actually plenty of headroom for additions, so we may do a lot of them without using 
  // mod and then only reduce once ...maybe that can be done better with radix-4 or radix-8 NTTs, 
  // when we can arrange the operations in a way to avoid a lot of the % operations. For 
  // multiplications, we need to immediately reduce, so let's try to avoid them as much as 
  // possible - for example, in the FFT routine, we compute twiddle factors on the fly by 
  // multiplications - but we may use tables instead

  ModInt operator*(const ModInt& b) { return (value * b.value) % modulus; }

  ModInt& operator+=(const ModInt& b) { *this = *this+b; return *this; }
  ModInt& operator*=(const ModInt& b) { *this = *this*b; return *this; }


  bool operator==(const ModInt& b) const { return value == b.value; }
  bool operator!=(const ModInt& b) const { return value != b.value; }


  operator rsUint64() const { return value; }

  // The magic numbers (definitions of the arrays are in .cpp file):
  static const rsUint64 modulus  = 3221225473;
  static const int      numRoots = 30;         // 30 is the number of (2^k)th roots of unity that exist for the modulus
  static const rsUint64 roots[numRoots];       // N-th roots of unity for N = 2^(k+1), k is array index
  static const rsUint64 rootsInv[numRoots];    // modular inverses of the roots
  static const rsUint64 lengthsInv[numRoots];  // modular inverses of the lengths N
  // They actually all fit into a 32 bit integer, so it's a bit wasteful to store them as 64 bit. 
  // But the 64 bit format is how they are needed in the computations...so..I'm not sure what's the 
  // best way to store them...If conversion from uint32 to uint64 is free, then it would seem 
  // better to store them 32 bit format...we'll see.

};

/** NTT-convolution routine based on rsModularIntegerNTT_64. */
std::vector<int> rsConvolveNTT(const std::vector<int>& x, const std::vector<int>& h);



//=================================================================================================

/** A baseclass for representing colors. A color is always represented as a triple of numbers in 
the range 0..1 but what those numbers mean may differ depending on the color space which is 
determined by the subclass. For example, in subclass rsColorRGB they mean red, green, blue and in 
rsColorHSL they mean hue, saturation, lightness. This baseclass stores the 3 values generically in 
the x,y,z members (which are inherited from the baseclass rsVector3D) and provides static functions
to convert between the various color spaces.

References:
  https://en.wikipedia.org/wiki/HSL_and_HSV
  https://www.rapidtables.com/convert/color/index.html  */

template<class T>
class rsColor : public rsVector3D<T>
{

public:



  rsColor(T x = T(0), T y = T(0), T z = T(0)) : rsVector3D<T>(x, y, z) {}


  rsColor(const rsVector3D<T>& v) : rsVector3D<T>(v) {}


  /** Converts HSL (hue, saturation, lightness) to RGB (red, green, blue). */
  static void hsl2rgb(T H, T S, T L, T* R, T* G, T* B);

  /** Converts RGB (red, green, blue) to HSL (hue, saturation, lightness). */
  static void rgb2hsl(T R, T G, T B, T* H, T* S, T* L);


  // needs more tests:
  static void lab2xyz(T L, T a, T b, T* X, T* Y, T* Z);
  static void xyz2lab(T X, T Y, T Z, T* L, T* a, T* b);

  // under construction:
  static void ch2ab(T C, T h, T* a, T* b);
  static void ab2ch(T a, T b, T* C, T* h);
  // this is the cyclidrical version of L*a*b, the L value remains the same, so it's not included
  // in the parameters

  static void xyz2rgb(T X, T Y, T Z, T* R, T* G, T* B);


  static void lab2rgb(T L, T a, T b, T* R, T* G, T* B);
  //static void rgb2lab(T R, T G, T B, T* L, T* a, T* b);
  static void lch2rgb(T L, T C, T h, T* R, T* G, T* B);


  /** Converts an RGB triple to a C-string with the hexadecimal representation, either with the 
  sharp symbol # prepended or not. So the char array must have space for 6 or 7 characters, 
  and one more for the null termination, if so selected. So, with a length of 8 characters, you're
  always on the save side, but you may get away with one or two less, when one or both flags are 
  false. */
  static void rgb2hex(T R, T G, T B, char* hex, bool withSharp = true, 
    bool withNullTermination = true);

  /** Works the same as rgb2hex(T, T, T, char*, bool, bool), just using unsigned chracters in 
  0..255 to represent the RGB values. */
  static void rgb2hex(unsigned char R, unsigned char G, unsigned char B, 
    char* hex, bool withSharp = true, bool withNullTermination = true);

  /** Convenience function to convert directly from HSL to hexadecimal. @see rgb2hex */
  static void hsl2hex(T H, T S, T L, char* hex, bool withSharp = true, 
    bool withNullTermination = true);

  // https://en.wikipedia.org/wiki/Web_colors#Hex_triplet

};
// needs test: try roundtrips between rgb and hsl and some special cases (pure colors, black, 
// white, gray, etc.)
// ToDo: 
// -implement more conversions - we don't currently need them but for the sake of completeness
// -maybe have functions like setHue, setLightness, setSaturation, setRed, setGreen, setBlue
// -have a function rgb2hex (and maybe hex2rgb too)


template<class T>
void rsColor<T>::hsl2rgb(T H, T S, T L, T* R, T* G, T* B)
{
  H  *= T(360);                                    // we expect 0 <= H <= 1 
  T C = (T(1) - rsAbs(T(2)*L-T(1))) * S;
  T X = C * (T(1) - rsAbs(fmod(H/T(60), T(2)) - T(1)));
  T m = L - C/2;
  if(     H <  60) { *R = C; *G = X; *B = 0; }
  else if(H < 120) { *R = X; *G = C; *B = 0; }
  else if(H < 180) { *R = 0; *G = C; *B = X; }
  else if(H < 240) { *R = 0; *G = X; *B = C; }
  else if(H < 300) { *R = X; *G = 0; *B = C; }
  else             { *R = C; *G = 0; *B = X; }
  *R += m;
  *G += m;
  *B += m;

  // see: https://www.rapidtables.com/convert/color/hsl-to-rgb.html
}

template<class T>
void rsColor<T>::rgb2hsl(T R, T G, T B, T* H, T* S, T* L)
{
  auto wrap = [](T x, T min, T max)  // todo: use library function
  {
    T range = max-min;
    while(x > max) x -= range;
    while(x < min) x += range;
    return x;
  };
  T Cmax = rsMax(R, G, B);
  T Cmin = rsMin(R, G, B);
  T D    = Cmax - Cmin;                                 // delta
  if(        D == 0) *H = 0;
  //else if(Cmax == R) *H = T(60) * rsWrapToInterval((G-B)/D, T(0), T(6));
  else if(Cmax == R) *H = T(60) * wrap((G-B)/D, T(0), T(6));
  else if(Cmax == G) *H = T(60) * ((B-R)/D + T(2));
  else if(Cmax == B) *H = T(60) * ((R-G)/D + T(4));
  *H /= T(360);                                         // convert from 0..360 to 0..1
  *L  = (Cmax + Cmin) / T(2);
  if(D != 0) *S = D / (T(1) - rsAbs(T(2) * *L - T(1))); // do we need this if?
  else       *S = 0;

  /*
  // for debug:
  T r,g,b;
  hsl2rgb(*H, *S, *L, &r, &g, &b);
  T dr = rsAbs(R - r);
  T dg = rsAbs(G - g);
  T db = rsAbs(B - b);
  T tol = 1.e-6;
  rsAssert(dr <= tol && dg <= tol && db <= tol);
  */

  // see: https://www.rapidtables.com/convert/color/rgb-to-hsl.html
}


template<class T>
void rsColor<T>::lab2xyz(T L, T a, T b, T* X, T* Y, T* Z)
{
  // input: L, a, b (in 0..1 or in 0..255?)
  // reference white Xr, Yr, Zr taken to be 1, 1, 1 ...nope! that's wrong!
  // i think, L may be in 0...100 and a,b, in -128...+128?

  // Use 
  // Xr = 95.0489, Yr = 100, Zr = 108.8840  for D65, or
  // Xr = 96.4212, Yr = 100, Zr = 82.5188   for D50
  // https://en.wikipedia.org/wiki/CIELAB_color_space#From_CIEXYZ_to_CIELAB
  // https://www.nixsensor.com/free-color-converter/


  T e  = T(216./24389);       // epsilon, rounded to 0.008856 in CIE standard
  T k  = T(24389./27);        // kappa, rounded to 903.3 in CIE standard
  T ke = T(8);                // k*e == 216/27 = 8

  T fy = (L+T(16)) / T(116);
  T fx = fy + a / T(500);
  T fz = fy - b / T(200);
  T fx3 = fx*fx*fx;           // fx^3
  T fz3 = fz*fz*fz;           // fz^3
  if(fx3 > e) *X = fx3;
  else        *X = (T(116)*fx-T(16)) / k;
  if(L  > ke) *Y = fy*fy*fy;
  else        *Y = L / k;
  if(fz3 > e) *Z = fz3;
  else        *Z = (T(116)*fz-T(16)) / k;

  // see http://www.brucelindbloom.com/index.html?Equations.html
}

template<class T>
void rsColor<T>::xyz2lab(T X, T Y, T Z, T* L, T* a, T* b)
{
  static const T e = T(216./24389);       // epsilon, rounded to 0.008856 in CIE standard
  static const T k = T(24389./27);        // kappa, rounded to 903.3 in CIE standard
  // maybe they should be static class members

  auto toF = [](T v)
  {
    if(v > e) return cbrt(v);
    else      return (k*v+T(16)) / T(116);
  };
  T fx = toF(X);
  T fy = toF(Y);
  T fz = toF(Z);

  *L = (116*fy - 16);
  *a = (500*(fx - fy));
  *b = (200*(fy - fz));
}


template<class T>
void rsColor<T>::ch2ab(T C, T h, T* a, T* b)
{
  h  *= 2*PI;        // experimental
  *a  = C * cos(h);  // don't we need a factor of 2*pi or something?
  *b  = C * sin(h); 
  // https://de.wikipedia.org/wiki/LCh-Farbraum
}

template<class T>
void rsColor<T>::ab2ch(T a, T b, T* C, T* h)
{
  *C = rsSqrt(*a * *a + *b * *b);
  *h = atan2(*b, *a);
  // https://de.wikipedia.org/wiki/LCh-Farbraum
}

template<class T>
void rsColor<T>::xyz2rgb(T X, T Y, T Z, T* R, T* G, T* B)
{
  //rsError("not yet ready to use");
  *R = X * ( 3.2404542) + Y * (-1.5371385) + Z * (-0.4985314);
  *G = X * (-0.9692660) + Y * ( 1.8760108) + Z * ( 0.0415560);
  *B = X * ( 0.0556434) + Y * (-0.2040259) + Z * ( 1.0572252);
  // Numbers taken from the XYZ -> sRGB matrix with D65 from here:
  // http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html

  /*
  // doesn't work - taken from http://www.easyrgb.com/en/math.php#text8
  *R = X *  3.2406 + Y * -1.5372 + Z * -0.4986;
  *G = X * -0.9689 + Y *  1.8758 + Z *  0.0415;
  *B = X *  0.0557 + Y * -0.2040 + Z *  1.0570;

  if(*R > 0.0031308) *R = 1.055 * pow(*R, (1. / 2.4)) - 0.055;
  else               *R = 12.92 * *R;
  if(*G > 0.0031308) *G = 1.055 * pow(*G, (1. / 2.4)) - 0.055;
  else               *G = 12.92 * *G;
  if(*B > 0.0031308) *B = 1.055 * pow(*B, (1. / 2.4)) - 0.055;
  else               *B = 12.92 * *B;
  */

  // clean this up!
}
// todo: add D65 to the function name, also implement the D50 variant
// there are many variations, how such a comversion can be done with respect to two things:
// -which matrix is used
// -which companding functions is used
//  -if gamma-companding is used, you have to also specify a value for gamma

// http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html ...how can we turn the matrices
// there into double precision? maybe by taking M, computing inv(M) in double precision? that gives
// us a double-precision matrix for inv(M). maybe we can also use the given inv(M) and invert it to 
// get a double preicison version of M. and then maybe take averages of given and computed? and 
// maybe iterate that until both matrices converge? the goal is to get as closely as possible to a 
// perfect roundtrip in double precision whiel at the same time staying as closely as possible to 
// the original given matrices. with these rounded matrices, the roundtrip will introduce
// an error that's determined by the single-precision format, even if double precision is used. 
// maybe write a general function: rsRefineInversion(rsMatrix& M, rsMatrix& Mi) 

// Libraries for color conversions:
// https://docs.python.org/3/library/colorsys.html
// https://github.com/python/cpython/blob/3.10/Lib/colorsys.py
// https://python-colormath.readthedocs.io/en/latest/conversions.html


template<class T>
void rsColor<T>::lab2rgb(T L, T a, T b, T* R, T* G, T* B)
{
  T X, Y, Z;
  lab2xyz(L, a, b, &X, &Y, &Z);  
  // returns very small values, roundtrip lab -> xyz -> lab has been tested in a unit test


  //X *= 100; Y *= 100; Z *= 100;  // ad-hoc test - nah! produces nonsense

  xyz2rgb(X, Y, Z, R, G, B);     // not yet tested
}
// https://stackoverflow.com/questions/7880264/convert-lab-color-to-rgb
// http://www.easyrgb.com/en/math.php#text8
// http://www.brucelindbloom.com/


template<class T>
void rsColor<T>::lch2rgb(T L, T C, T h, T* R, T* G, T* B)
{
  T a, b;
  ch2ab(C, h, &a, &b);        // not yet tested
  lab2rgb(L, a, b, R, G, B);
}


/*
template<class T>
void rsColor<T>::rgb2lab(T R, T G, T B, T* L, T* a, T* b)
{

}
// https://de.wikipedia.org/wiki/Lab-Farbraum#Umrechnung_von_RGB_zu_Lab
*/


template<class T>
void rsColor<T>::rgb2hex(T R, T G, T B, char* hex, bool sharp, bool null)
{
  using uchar = unsigned char;
  uchar r = (uchar) round(T(255) * R);
  uchar g = (uchar) round(T(255) * G);
  uchar b = (uchar) round(T(255) * B);
  rgb2hex(r, g, b, hex, sharp, null);
}

template<class T>
void rsColor<T>::rgb2hex(unsigned char R, unsigned char G, unsigned char B, 
  char* hex, bool sharp, bool null)
{
  auto toHex = [](unsigned char c)
  {
    rsAssert(c >= 0 && c <= 15);
    if(c >= 10) return c + 55;          // 65: ASCII code of A and we need to subtract 10
    else        return c + 48;          // 48: ASCII code of 0
  };
  int s = 4;                            // shift
  int m = 15;                           // mask
  int i = 0;                            // index
  if(sharp) { hex[i] = '#'; i++; }
  hex[i] = toHex((R >> s) & m); i++;
  hex[i] = toHex( R       & m); i++;
  hex[i] = toHex((G >> s) & m); i++;
  hex[i] = toHex( G       & m); i++;
  hex[i] = toHex((B >> s) & m); i++;
  hex[i] = toHex( B       & m); i++;
  if(null) hex[i] = '\0';
}

template<class T>
void rsColor<T>::hsl2hex(T H, T S, T L, char* hex, bool sharp, bool null)
{
  T R, G, B;
  hsl2rgb(H, S, L, &R, &G, &B);
  rgb2hex(R, G, B, hex, sharp, null);
}

template<class T>
class rsColorRGB : public rsColor<T>
{
public:
  using rsColor<T>::rsColor;
  T getRed()   const { return this->x; }
  T getGreen() const { return this->y; }
  T getBlue()  const { return this->z; }
};

template<class T>
class rsColorHSL : public rsColor<T>
{
public:
  using rsColor<T>::rsColor;
  T getHue()        const { return this->x; }
  T getSaturation() const { return this->y; }
  T getLightness()  const { return this->z; }
};

template<class T>
class rsColorLCH : public rsColor<T>
{
public:
  using rsColor<T>::rsColor;
  T getLightness() const { return this->x; }
  T getChroma()    const { return this->y; }
  T getHue()       const { return this->z; }
};


// todo: implement:
// https://en.wikipedia.org/wiki/CIELUV
// https://en.wikipedia.org/wiki/CIELAB_color_space#CIELAB
// https://de.wikipedia.org/wiki/Lab-Farbraum
// the wikipedia article on HSL/HSV says:
// "The usual formulations of HSB and HLS are flawed with respect to the properties of color 
// vision. Now that users can choose colors visually, or choose colors related to other media 
// (such as PANTONE), or use perceptually-based systems like L*u*v* and L*a*b*, HSB and HLS should
// be abandoned."
// https://en.wikipedia.org/wiki/HSL_and_HSV#Disadvantages
// i think, in the L*u*v and L*a*b spaces, the we can compute the a,b values as cos(H), sin(H)?
// ...ah, yes, indeed:  
// https://en.wikipedia.org/wiki/CIELAB_color_space#Cylindrical_model
// https://en.wikipedia.org/wiki/CIELUV#Cylindrical_representation_(CIELCH)
// https://de.wikipedia.org/wiki/LCh-Farbraum
// https://sensing.konicaminolta.us/us/blog/understanding-the-cie-lch-color-space/
// we need CIE_HLC or CIE_HLC_uv ...or maybe (a cylidrical version) YUV:
// https://en.wikipedia.org/wiki/YUV
// http://www.niwa.nu/2013/05/understanding-yuv-values/
// ..it has no real yellow, but maybe that's actually good for grpahs in plots

// https://en.wikipedia.org/wiki/CIE_1931_color_space  for XYZ
// ""The CIE XYZ color space encompasses all color sensations that are visible to a person with 
// average eyesight. That is why CIE XYZ (Tristimulus values) is a device-invariant representation
// of color. It serves as a standard reference against which many other color spaces are defined.


//=================================================================================================

/** Subclass of rsLadderFilter with some additional functionality that shall not (yet) go into the
production version. */

template<class TSig, class TPar>
class rsLadderTest : public RAPT::rsLadderFilter<TSig, TPar>
{

public:

  /** Returns the transfer function as rsRationalFunction object. This is mainly useful for 
  research and development and not suitable for use actual products and a total no-go to use at 
  realtime. For plots in products, you should probably use getTransferFunctionAt or 
  getMagnitudeResponseAt. */
  rsRationalFunction<TPar> getTransferFunction(bool withGain = true);

  /** Old implementation of getTransferFunction, using rsRationalFunction's arithmetic instead of
  just assigning the coeffs via analytically derived formulas (as the new one does). It's less 
  efficient and less precise than the new one, but nicely demonstrates how rsRationalFunction can 
  be used for such computations. It will produce a function that is formally 8-pole, but features 
  pole/zero cancellations. */
  rsRationalFunction<TPar> getTransferFunctionOld();

};

//=================================================================================================

/** Just for testing, if a class compiles that defines comparison operators with same input 
parameters but different return types. */

/*
template<class T>
class rsOperatorTest
{

public:

  rsOperatorTest(T a) { v = a; }


  inline bool operator<(const rsOperatorTest& b) const 
  {
    return (v < b.v);
  }

  //inline T operator<(const rsOperatorTest& b) const 
  //{
  //  if(v < b.v) return T(1);
  //  else        return T(0);
  }
  // Nope - doesn't compile - gives error: "operator redefinition ...". I think, in the simd
  // classes, the implementation should return a bool which should implement a logical "and" of
  // all the scalar results. If we need a 0 or 1 of the given type for branchless code, a function
  // should be used (like compareLess, compareLessOrEqual, etc.)

protected:

  T v = T(0);

};
*/


//=================================================================================================
// the stuff below is just for playing around - maybe move code elsewhere, like the research-repo:

/** A class for representing a particular kind of string with which we can do some computations
just like with numbers. The set of all such strings forms a group (see group theory). The group
operation (which we call addition here) is to concatenate two strings and then delete all pairs of
equal characters, i.e. the string "aaab" would be reduced to "ab", one "aa" pair is deleted. The
inverse element to each string is obtained by reversing it. Adding "cba" to "abc" like abc+cba
results in abccba wich would subsequently be reduced to the empty string (the rule of deleting
equal pairs is used as often as applicable). The additive neutral element is the empty string. */

class rsGroupString
{

public:

  rsGroupString() {}

  rsGroupString(const std::vector<unsigned int>& initialString) { s = initialString; }

  rsGroupString(int length) { s.resize(length); }

  // define operator =,
  // maybe: < (lexicographical order), * (i still have to invent a suitable multiplication rule)

  bool operator==(const rsGroupString& t) const { return t.s == s; }


  /** Adds two GroupString objects. This addition is the group operation and is (conceptually)
  performed by concatenating two strings and then deleting all doublets (iteratively, as often as
  necessary to eliminate all of them). */
  rsGroupString operator+(const rsGroupString &rhs) const;

  /** Unary minus. Returns the additive inverse */
  rsGroupString operator-() const { return inverse(); }

  /** Binary subtraction by adding the additive inverse. */
  rsGroupString operator-(const rsGroupString &rhs) { return *this + (-rhs); }

  /** Read/write access to i-th character. */
  unsigned int& operator[](int i) { return s[i]; }

  unsigned int last() const
  {
    if(s.size() > 0) return s[s.size()-1];
    else             return 0; // is this reasonable? or should we use a special "error" signal
  }

  void append(unsigned int x) { s.push_back(x); }

  void removeLast()
  {
    if(s.size() > 0)
      s.pop_back();
  }

  /** Length of the string. */
  int length() const { return (int) s.size(); }



  /** Returns the (additive) inverse which is just the string in reversed order. */
  rsGroupString inverse() const;
   // maybe later (when we have multiplication), rename to additiveInverse and implement a
   // multiplicativeInverse, too
   // then the class should be renamed to fieldString

   // maybe let integers 0 and 1 be used and implement 1/s = s.multiplicativeInverse, etc.


  void resize(int newSize) { s.resize((int)newSize); }

  std::vector<unsigned int> get() const { return s; }

  // bool hasDoublets();
  // void removeDoublets();


protected:



  std::vector<unsigned int> s;  // we represent the characters as unsigned integers

  //int modulus = 26;
  // the modulus, we use the 26 lowercase letters, but that is tweakable...but we don't need that
  // yet but maybe later when we do operations on individual characters

};

/** Subclass of rsGroupString that lets use more conveniently work with strings over the alphabet
a,b,c,..,x,y,z. The class provides the conversions from/to std::string comparison operators etc.
But all these convenenience functions have nothing to do with the actual algebraic structure, which
is why they have been factored out to keep the baseclass pure. */

class rsGroupString2 : public rsGroupString
{

public:

  /** Convenience function meant to be used for strings over the aplhabet a,b,c,...,x,y,z. We
  represent 'a' as 0 and then count up to 'z' = 25. */
  rsGroupString2(const char* initialString);

  //rsGroupString2(const std::string& initialString);

  /** Constructor for conversion from baseclass object (?) */
  rsGroupString2(const rsGroupString& gs);





  /** Converts to a std::string. */
  std::string toString() const;

  /** Comparison operator (for some reason, we need to override it or the compiler balks). */
  bool operator==(const rsGroupString2& t) const  { return t.s == s; }

  //bool operator==(const std::string& str) const { return str == toString(); }  // obsolete?

  //rsGroupString operator+(const rsGroupString &rhs) const;

  //rsGroupString operator+(const std::string& &rhs) const
  //{

  //}



  /** Checks, if the passed unsigned integer corresponds to one of the allowed characters, i.e. is
  from the alphabet. */
  static bool isLowerCaseLetter(unsigned int c) { return c >= 97 && c <= 122; } // 'a' = 97, 'z' = 122

  /** Checks, if all characters are inthe valid range, i.e. inside our restricted alphabet. */
  bool checkCharacters() const;

};


#endif
