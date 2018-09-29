#ifndef rosic_Interpolation_h
#define rosic_Interpolation_h

// commented functions are superseded by corresponding template functions in RAPT (and may eventually
// be deleted)

namespace rosic
{

  /** Fits a cubic polynomial of the form: 
  \f[ f(x) = a3*x^3 + a2*x^2 + a1*x + a0  \f]
  to two points (x1,y1), (x2,y2) and matches values of the derivative (given by yd1, yd2) at these 
  points. */
  void fitCubicWithDerivative(double x1, double x2, double y1, double y2, double yd1, 
    double yd2, double *a3, double *a2, double *a1, double *a0);

  /** Similar to fitCubicWithDerivative, but the x-coodinates of the two points are fixed at x0=0, x1=1 such that we fit the points
  (0,y0), (1,y1) and match values of the derivative (given by yd0, yd1) there. This simplifies the computation a lot compared to the 
  general case. */
  void fitCubicWithDerivativeFixedX(double y0, double y1, double yd0, double yd1, double *a3, double *a2, double *a1, double *a0);

  /** Similar to fitCubicWithDerivativeFixedX, but fits a quintic (5th order) polynomial in order to additionaly match desired values for 
  the 2nd derivatives (given by ydd0, ydd1) at the sample points. */
  void fitQuinticWithDerivativesFixedX(double y0, double y1, double yd0, double yd1, double ydd0, double ydd1, double *a5, double *a4,  
                                       double *a3, double *a2, double *a1, double *a0);

  /** Computes coefficients for a polynomial that passes through the points (x0 = 0, y0[0]), (x1 = 1, y1[0]) and in addition to passing
  through these points, it also matches a number "M" of derivatives to prescribed values. These values should be passed in y0[1], 
  y0[2], etc. for the values at x0 = 0 and in y1[1], y1[2], etc. for the values at x1 = 1. The argument "a" is used to return the computed 
  coefficients. 
  y0: (M+1)-element array containing y(0), y'(0), y''(0), y'''(0), etc.
  y1: (M+1)-element array containing y(1), y'(1), y''(1), y'''(1), etc.
  a:  (2*M+1)-element array for returning a0, a1, a2, a3, a4, a5, a6, a7, etc.  */
  //void getHermiteCoeffsM(double *y0, double *y1, double *a, int M);

  /** Optimized version of getHermiteCoeffsM for the case M == 1. */
  //void getHermiteCoeffs1(double *y0, double *y1, double *a);

  /** Optimized version of getHermiteCoeffsM for the case M == 2. */
  //void getHermiteCoeffs2(double *y0, double *y1, double *a);

  /** Optimized version of getHermiteCoeffsM for the case M == 3. */
  //void getHermiteCoeffs3(double *y0, double *y1, double *a);

  /** Computes a delayed sample with a fractional delay of "d" (0 <= d <= 1) behind y[0]. To compute the output value, the function uses 
  y[0], y[-1], y[-2], ..., y[-(M+1)] to obtain finite difference approximations for a number "M" of derivatives and uses these as desired 
  derivatives for a Hermite interpolation. The resulting interpolating function (seen in the continuous time domain) will pass through all 
  the sample-points and has M continuous derivatives. The continuous time impulse response of the interpolator is asymmetric due to the way
  in which the derivatives are approximated (using only past values). The "shape" parameter controls, which value will actually be used for 
  the desitred derivative by multiplying the raw finite difference approximation for the n-th derivative by shape^n.  */
  double getDelayedSampleAsymmetricHermiteM(double d, double *y, int M, double shape = 1.0);

  /** Optimized version of getDelayedSampleAsymmetricHermiteM for the case M == 1. */
  double getDelayedSampleAsymmetricHermite1(double d, double *y, double shape = 1.0);

  /** Computes a delayed sample with a fractional delay of "d" (0 <= d <= 1) behind y[0] by linearly interpolating between y[0] and 
  y[-1]. */
  double getDelayedSampleLinear(double d, double *y);

  // \todo make symmetric versions for Hermite interpolators (using symmetric finite differences), make similar functions for Lagrange 
  // interpolation and sinc-interpolation


  /** Fits a cubic polynomial of the form: 
  \f[ f(x) = a*x^3 + b*x^2 + c*x + d  \f]
  to four points (x0,y0), (x1,y1),(x2,y2), (x3,y3)  */
  void fitCubicThroughFourPoints(double x0, double y0, double x1, double y1, double x2, double y2,
    double x3, double y3, double *a, double *b, double *c, double *d);

}  // end namespace rosic

#endif // rosic_Polynom_h
