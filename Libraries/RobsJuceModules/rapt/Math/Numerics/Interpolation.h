#ifndef RAPT_INTERPOLATION_H
#define RAPT_INTERPOLATION_H

// \todo wrap into class rsInterpolation ...or maybe merge with Interpolator from the Unfinished
// folder

/** Given two arrays of input abscissa- and ordinate values xIn, yIn of length inLength and an
array of new abscissa values xOut, this function fills the array yOut with values that correspond 
to the xOut values by linearly interpolating the yIn array. Array xOut and yOut are of length 
outLength. xIn must be strictly monotonically increasing. If the xOut array contains values outside
the range of the xIn value, the will linearly extrapolate beyond the original range based on the 
two initila and/or final points. */
template<class Tx, class Ty>
void resampleNonUniformLinear(const Tx* xIn, const Ty* yIn, int inLength, 
  const Tx* xOut, Ty* yOut, int outLength);
// rename

/** Computes coefficients for a cubic polynomial that goes through the 4 points (-1,y[-1]),
(0,y[0]), (1,y[1]), (2,y[2]) that will have matched 1st derivatives when used on successive
positions in an y-array. [TEST THIS! it has not yet been tested] */
template<class T>
void rsCubicSplineCoeffsFourPoints(T *a, T *y);

/** Fits a cubic polynomial of the form:
\f[ f(x) = a3*x^3 + a2*x^2 + a1*x + a0  \f]
to two points (x1,y1), (x2,y2) and matches values of the derivative (given by yd1, yd2) at these
points.
\todo change order of a coeffs - but make sure that all client code using this is updated
accordingly. or maybe pass an array a[4] - this will force client code to be updated  */
template<class T>
void fitCubicWithDerivative(T x1, T x2, T y1, T y2, T yd1,
  T yd2, T *a3, T *a2, T *a1, T *a0);

/** Similar to fitCubicWithDerivative, but the x-coodinates of the two points are fixed at x0=0,
x1=1 such that we fit the points (0,y0), (1,y1) and match values of the derivative (given by yd0,
yd1) there. This simplifies the computation a lot compared to the general case. */
template<class T>
void fitCubicWithDerivativeFixedX(T y0, T y1, T yd0, T yd1,
  T *a3, T *a2, T *a1, T *a0);

/** Similar to fitCubicWithDerivativeFixedX, but fits a quintic (5th order) polynomial in order
to additionaly match desired values for the 2nd derivatives (given by ydd0, ydd1) at the sample
points. */
template<class T>
void fitQuinticWithDerivativesFixedX(T y0, T y1, T yd0, T yd1,
  T ydd0, T ydd1, T *a5, T *a4, T *a3, T *a2, T *a1,
  T *a0);

/** Computes coefficients for a polynomial that passes through the points (x0 = 0, y0[0]),
(x1 = 1, y1[0]) and in addition to passing through these points, it also matches a number "M" of
derivatives to prescribed values. These values should be passed in y0[1], y0[2], etc. for the
values at x0 = 0 and in y1[1], y1[2], etc. for the values at x1 = 1. The argument "a" is used to
return the computed 2*M+2 polynomial coefficients for the polynomial of order 2*M+1.
y0: (M+1)-element array containing y(0), y'(0), y''(0), y'''(0), etc.
y1: (M+1)-element array containing y(1), y'(1), y''(1), y'''(1), etc.
a:  (2*M+2)-element array for returning a0, a1, a2, a3, a4, a5, a6, a7, etc.  */
template<class T>
void getHermiteCoeffsM(const T *y0, const T *y1, T *a, int M);

/** Optimized version of getHermiteCoeffsM for the case M == 1. */
template<class T>
void getHermiteCoeffs1(const T *y0, const T *y1, T *a);

/** Optimized version of getHermiteCoeffsM for the case M == 2. */
template<class T>
void getHermiteCoeffs2(const T *y0, const T *y1, T *a);

/** Optimized version of getHermiteCoeffsM for the case M == 3. */
template<class T>
void getHermiteCoeffs3(const T *y0, const T *y1, T *a);

/** Computes a delayed sample with a fractional delay of "d" (0 <= d <= 1) behind y[0]. To
compute the output value, the function uses y[0], y[-1], y[-2], ..., y[-(M+1)] to obtain finite
difference approximations for a number "M" of derivatives and uses these as desired derivatives
for a Hermite interpolation. The resulting interpolating function (seen in the continuous time
domain) will pass through all the sample-points and has M continuous derivatives. The continuous
time impulse response of the interpolator is asymmetric due to the way in which the derivatives
are approximated (using only past values). The "shape" parameter controls, which value will
actually be used for  the desired derivative by multiplying the raw finite difference
approximation for the n-th derivative by shape^n.  */
template<class T>
T getDelayedSampleAsymmetricHermiteM(T d, T *y, int M, T shape = 1.0);

/** Optimized version of getDelayedSampleAsymmetricHermiteM for the case M == 1. */
template<class T>
T getDelayedSampleAsymmetricHermite1(T d, T *y, T shape = 1.0);

/** Computes a delayed sample with a fractional delay of "d" (0 <= d <= 1) behind y[0] by
linearly interpolating between y[0] and y[-1]. */
template<class T>
T getDelayedSampleLinear(T d, T *y);

// \todo make symmetric versions for Hermite interpolators (using symmetric finite differences),
// make similar functions for Lagrange interpolation and sinc-interpolation

/** Fits a cubic polynomial of the form:
\f[ f(x) = a*x^3 + b*x^2 + c*x + d  \f]
to four points (x0,y0), (x1,y1),(x2,y2), (x3,y3)  */
template<class T>
void fitCubicThroughFourPoints(T x0, T y0, T x1, T y1, T x2,
  T y2, T x3, T y3, T *a, T *b, T *c, T *d);

/** Given two x-axis values x1, x2 and the corresponding y-axis values y1, y2, this function
returns the y-value corresponding to x by linearly interpolating between x1, x2. If x is outside
the range of x1, x2, the function will extrapolate. */
template<class Tx, class Ty>
Ty rsInterpolateLinear(Tx x1, Tx x2, Ty y1, Ty y2, Tx x);
// maybe change interface to x1, y1, x2, y2 to make it consistent with other functions for example
// in rsLine2D
// but this is a change that would silently break client code
// ...maybe deprecate these free functions and replace them by a class rsInterpolator with a new 
// API
// hmm...but then rsInterpolateCubicHermite would also have to be changed

/** !!! NOT YET TESTED !!!
Linearly interpolates between two values x0, x1 where the values are supposed to be wrapping 
around and always be in xMin..xMax, for example -1..1, 0..1, 0..2*pi, -pi..pi, etc. t is the 
interpolation parameter between 0..1 such that x = (1-t)*x0 + t*x1 = x0 + t*(x1-x0). When the value
is restricted to a given range, we have two possible directions of the interpolation path - we can 
go forward or backward from x0 to x1 (with wrap around) as t traverses 0..1. This function chooses
the direction for which the absolute distance between x0 and x1 is shorter. */
template<class T> 
T rsInterpolateWrapped(T x0, T x1, T t, T xMin, T xMax);

/** Given two length N arrays x, y with x-axis values and corresponding y-axis values, this
function fills the array yi with values corresponding to the xi by linear interpolation
(or extrapolation, if necessary) of the input data. The xi and yi arrays are of length Ni. */
template<class T>
void rsInterpolateLinear(const T *x, const T *y, int N, const T *xi, T *yi, int Ni);

/** Given four x-axis values x1, x2, x3, x4 and the corresponding y-axis values y1, y2, y3, y4,
this function returns the y-value corresponding to x by cubic hermite interpolation. It is
assumed that the value of x is between x2 and x3 and the interpolation will be such that
successive interpolants on the same data set will have matching derivatives at the data points
(in addition to pass through the points). It can be used for interpolating data sets for which
the x-axis values are not equidistant. */
template<class T>
T rsInterpolateCubicHermite(T x1, T x2, T x3, T x4, T y1, T y2, T y3, T y4, T x);

/** Given two length N arrays x, y with x-axis values and corresponding y-axis values, and a
2D array of size MxN containing M derivative values at the sample points y', y'', this function
fills the array yi with values corresponding to the xi by spline interpolation (or extrapolation,
if necessary) of the input data. The xi and yi arrays are of length Ni. The number M gives the
number of derivatives that should match at the data points where two subsequent splines join and
the actual values of the m-th derivative at a data-point with index n is given by yd[m-1][n].

\todo have also an optional parameter **ydi where interpolated values of the derivative will be
stored - we can just use evaluatePolynomialAndDerivativesAt instead of evaluatePolynomialAt and
do some copying
*/
template<class Tx, class Ty>
void rsInterpolateSpline(const Tx *x, const Ty *y, Ty **yd, int N, int M, 
  const Tx *xi, Ty *yi, int Ni);
// todo: rename to rsInterpolateHermite - a spline is different - it doesn't prescribe derivatives
// and instead makes higher order deriavtives match at the nodes

/** Given two length N arrays x, y with x-axis values and corresponding y-axis values, this
function fills the array yi with values corresponding to the xi by spline interpolation
(or extrapolation, if necessary) of the input data. The xi and yi arrays are of length Ni. The
smoothness parameter gives the number of derivatives that should match at the data points where
two subsequent splines join. The higher this value, the smoother the resulting spline function
will be in terms of continuous derivatives, but it will also tend to oscillate more between the
datapoints for higher smoothness values. With a value of 1, which is the default, a cubic spline
will be used and the 1st derivative will match at the data points. Generally, a polynomial of
order 2*smoothness+1 will be used. */
template<class Tx, class Ty>
void rsInterpolateSpline(const Tx *x, const Ty *y, int N, const Tx *xi, Ty *yi, int Ni, 
  int smoothness = 1);
// -rename to rsInterpolateHermite
// -try to avoid the oscillatory behavior by using downscaled numerical derivatives - let the user 
//  pass a factor k between 0..1: 1 means no downscaling (take numerical derivatives as is), 0 
//  means all derivatives are taken to be zero - maybe the 1st derivative should be multiplied by 
//  k, the 2nd by k^2, etc. - this will happen naturally, when we just always multiply the 
//  estimates taken form previous numercila derivatives by k (factors accumulate)...but maybe it
//  may also make sense to have different factors for different derivatives


/** Given arrays....

The value scaleRhs is a factor by which the interpolant can be morphed between linear (=0) and the
actual, proper cubic spline (=1). It doesn't appear in the literature, but it may be a useful tweak
in some situations. */
template<class Tx, class Ty>
void rsNaturalCubicSpline(const Tx *x, const Ty *y, int N, const Tx *xi, Ty *yi, int Ni, Ty scaleRhs = Ty(1));
// rename to rsInterpolateNaturalCubic, or maybe resampleNaturalSpline...or 
// rsResampler::naturalSpline
// ...or maybe wrap into class rsInterpolator

// todo: implement a natural cubic spline interpolation as describen here:
// http://mathworld.wolfram.com/CubicSpline.html
// i think, this can use the method above by assigning derivative values according to the 
// tridiagonal system at the bottom instead of using a finite difference approximation
// other links:
// http://www.maths.nuigalway.ie/~niall/teaching/Archive/1617/MA378/2-2-CubicSplines.pdf

// http://mathonline.wikidot.com/natural-cubic-spline-function-interpolation
// http://mathonline.wikidot.com/natural-cubic-spline-function-interpolation-examples-1


/** Given two points (x0,y0), (x1,y1), this function computes the cubic polynomial coefficients for
a spline arc that has the parametric equations:
x(t) = a0 + a1*t + a2*t^2 + a3*t^3
y(t) = b0 + b1*t + b2*t^2 + b3*t^3
The desired derivatives of x and y with respect to the parameter t at the two points are passed in 
by (dx0,dy0) and (dx1,dy1) and the polynomial coefficients are written into the 
arrays a and b which are supposed to be of length 4. */
template<class T>
void cubicSplineArcCoeffs2D(T x0, T dx0, T y0, T dy0, T x1, T dx1, T y1, T dy1, T* a, T* b);
// make notation consisten with quadratic case
// use x0, x1 throughut the library instead of x1, x2 - they correspond to parameter-values 0 and 1
// and also counting from zero is customary in programming

/** Similar to cubicSplineArcCoeffs2D but tries to fit a quadratic spline between the two points. 
It doens't require dx/dt and dy/dt to have certian values at the joints but only their quotient 
(dx/dt)/(dy/dt) = dx/dy must be equal at the joints. This is less restrictive but still gives a 
curve that is 2nd order continuous in 2D (although both 1D function may be not). However, this may
sometimes fail due to division by zero conditions. In this case, the function does nothing and just 
returns false (otherwise true). */
template<class T>
bool quadraticSplineArcCoeffs2D(T x0, T dx0, T y0, T dy0, T x1, T dx1, T y1, T dy1, T* a, T* b);

/** Tries to fit a quadratic (via quadraticSplineArcCoeffs2D) and if this fails, falls back to
cubicSplineArcCoeffs2D. */
template<class T>
void quadraticOrCubicSplineArcCoeffs2D(T x0, T dx0, T y0, T dy0, T x1, T dx1, T y1, T dy1, 
  T* a, T* b);

/** Given the polynomial coeffcient arrays of a cubic spline in 2D described by
x(t) = a0 + a1*t + a2*t^2 + a3*t^3
y(t) = b0 + b1*t + b2*t^2 + b3*t^3
this function numercially approximates the arc length at the given values of t (passed in as input 
array) and stores them in the output array s. Both arrays are assumed to be of length N. */
template<class T>
void cubicSplineArcLength2D(T* a, T* b, T* t, T* s, int N);
// maybe this should go to the Geometry folder?



//===============================================================================================
// implementation of template-functions (move to cpp file):

template<class Tx, class Ty>
Ty rsInterpolateLinear(Tx x1, Tx x2, Ty y1, Ty y2, Tx x)
{
  if(x1 == x2)
    return Ty(0.5) * (y1+y2); // at a discontinuity, we return the average
  Ty a = (y2-y1) / (x2-x1);
  Ty b = y1 - a*x1;
  return a*x + b;
    // factor out computation of a and b
}

template<class T>
void rsInterpolateLinear(const T *x, const T *y, int N, const T *xi, T *yi, int Ni)
{
  int n = 0;  // index into input data
  int i = 0;  // index into interpolated data
  T a, b;     // parameters of the line y = a*x + b

  while(n < N-1)
  {
    a = (y[n+1]-y[n]) / (x[n+1]-x[n]);
    b = y[n] - a*x[n];
    while(xi[i] < x[n+1] && i < Ni)
    {
      yi[i] = a*xi[i] + b;
      i++;
    }
    n++;
  }

  while(i < Ni)
  {
    yi[i] = a*xi[i] + b;
    i++;
  }
}

template<class T>
T rsInterpolateCubicHermite(T x1, T x2, T x3, T x4, T y1, T y2, T y3, T y4, T x)
{
  T s, d1, d3, s2, s3, a[4];
  s  = T(1)/(x3-x2);
  d1 = s*(x2-x1);
  d3 = s*(x4-x2)-T(1);
  s2 = ((y2-y1)/d1 + (y3-y2)*d1) / (d1+T(1));
  s3 = ((y3-y2)*d3 + (y4-y3)/d3) / (d3+T(1));
  fitCubicWithDerivativeFixedX(y2, y3, s2, s3, &a[3], &a[2], &a[1], &a[0]);
  return rsPolynomial<T>::evaluate(s*(x-x2), a, 3);
  // maybe factor out a function that returns the polynomial coefficients (and s) because the
  // same set of coefficients may get used to interpolate at multiple values of x between the
  // same x2 and x3 and it is wasteful to recompute the coefficients each time
}

template<class Tx, class Ty>
void rsInterpolateSpline(const Tx *x, const Ty *y, Ty **yd, int N, int M, 
  const Tx *xi, Ty *yi, int Ni)
{
  int n = 0;              // index into input data
  int i = 0;              // index into interpolated data
  int m;                  // index of the derivative
  Tx scale, shift;        // scaler and shifter for the input value for the polynomial
  Ty *a  = new Ty[2*M+2]; // polynomial coefficients
  Ty *y0 = new Ty[M+1];   // y0 values and derivatives passed to getHermiteCoeffsM
  Ty *y1 = new Ty[M+1];   // y1 values and derivatives passed to getHermiteCoeffsM

  while(n < N-1)
  {
    // compute coeffs for spline and shift- and scale values:
    shift = x[n];                    // maybe rename to x0
    scale = x[n+1]-x[n];             // maybe rename to dx
    y0[0] = y[n];                    // y-values are taken as is
    y1[0] = y[n+1];
    for(m = 1; m <= M; m++)
    {
      y0[m] = yd[m-1][n]   * scale;  // y' and higher derivatives have to be scaled by dx
      y1[m] = yd[m-1][n+1] * scale;
    }
    getHermiteCoeffsM(y0, y1, a, M);

    // extra-/interpolate:
    scale = Tx(1) / scale;
    while(xi[i] < x[n+1] && i < Ni)
    {
      yi[i] = rsPolynomial<Ty>::evaluate(scale*(xi[i]-shift), a, 2*M+1);
      i++;
    }

    n++;
  }

  // extrapolate tail:
  while(i < Ni)
  {
    yi[i] = rsPolynomial<Ty>::evaluate(scale*(xi[i]-shift), a, 2*M+1);
    i++;
  }

  // cleanup:
  delete[] a;
  delete[] y0;
  delete[] y1;
}

template<class Tx, class Ty>
void rsInterpolateSpline(const Tx *x, const Ty *y, int N, const Tx *xi, Ty *yi, int Ni, int M)
{
  // compute numeric derivatives of y, to be used for the spline at data points:
  Ty **yd = nullptr;
  if(M > 0)
  {
    rsMatrixTools::allocateMatrix(yd, M, N);
    rsNumericDifferentiator<Ty>::derivative(x, y, yd[0], N);
    for(int m = 1; m < M; m++)
      rsNumericDifferentiator<Ty>::derivative(x, yd[m-1], yd[m], N);
  }

  // interpolate with these numeric derivatives and cleanup:
  rsInterpolateSpline(x, y, yd, N, M, xi, yi, Ni);
  if(M > 0)
    rsMatrixTools::deallocateMatrix(yd, M, N);
}

#endif
