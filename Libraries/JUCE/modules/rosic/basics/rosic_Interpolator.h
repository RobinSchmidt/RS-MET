#ifndef rosic_Interpolator_h
#define rosic_Interpolator_h

// rosic-indcludes:
#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This class implements some interpolation formulas to reconstruct a continuous signal from samples,
  that is: it allows for reading out the signal at arbitrary non-integer time instants. Most of the
  member-functions which perform the interpolation are declared as static such that one does not
  need to actually create an instance of the class if one wants to use them. However, some other
  interpolation-methods (namely the allpass-methods) need to keep track of internal state
  variables - these are handled here as member variables, thus those interpolation-methods can not
  work as static methods. So if you want to use them, you will need to instantiate an object of
  the Interpolator class.

  References:
   -Olli Niemitalo: Polynomial Interpolators for High-Quality Resampling of Oversampled Audio
   -Jon Datorro: Effect Design Part 2 - Delay-Line Modulation and Chorus (JAES, Vol 45, No. 10)


   \todo: avoid recalculation of the warped allpas coefficient when it stays the same between
   successive calls

  */

  class Interpolator
  {

  public:

    /** An enumeration of the available interpolation methods. */
    enum interpolationMethods
    {
      TRUNCATED_INDEX = 0,
      NEAREST_NEIGHBOUR,
      LINEAR,
      CUBIC_ZERO_DERIVATIVE,
      ALLPASS,
      WARPED_ALLPASS,           // up to here, we can do delaymodulation downto zero

      HERMITE_4P_3O,
      PARABOLIC_4P_2O_2X,
      //.....


      WARPED_ALLPASS_APPROXIMATE
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Interpolator();

    /** Destructor. */
    ~Interpolator();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Chooses one of the interpolation methods. @see enum interpolationMethods. */
    void setInterpolationMethod(int newMethod);

    //---------------------------------------------------------------------------------------------
    // inquiry (get-, is-, etc. functions):

    /** Returns the currently chosen interpolation method. @see enum interpolationMethods. */
    int getInterpolationMethod() const { return interpolationMethod; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one interpolated output-sample. It will use  the interpolation-method that has
    been selected via the setInterpolationMethod()-function. */
    INLINE double getSample(double x, double *y);

    /** Calculates an output sample using allpass interpolation as described in Datorro's paper -
    note that the x-value used here as input argument corresponds to 1-frac in the papaer
    because the Interpolator class sees it from the viewpoint of a phase-accumulating oscillator
    whereas the paper deals with delaylines. */
    INLINE double getSampleAllpass(double x, double *y);

    /** Calculates an output sample using warped allpass interpolation as described in Datorro's
    paper - note that the x-value used here as input argument corresponds to 1-frac in the papaer
    because the Interpolator class sees it from the viewpoint of a phase-accumulating oscillator
    whereas the paper deals with delaylines. */
    INLINE double getSampleWarpedAllpass(double x, double *y);

    /** Returns the nearest neighbour sample - that is y[0] when x < 0.5 and y[1] otherwise. */
    static INLINE double getSampleNearestNeighbour(double x, double *y);

    /** Calclulates one output sample by linearly interpolating between y[0] and y[1] at
    position x. y should point to an array of at least two values, x should be a number
    between 0...1. */
    static INLINE double getSampleLinear(double x, double *y);

    /** Calclulates one output sample by linearly interpolating between y[0] and y[1] at
    position x. y should point to an array of at least two values, x should be a number between
    0...1. */
    static INLINE double getSampleLinear(double x, float *y);



    static INLINE double getSampleCubicZeroDerivative(double x, double *y);

    /** Calclulates one output sample by interpolating parabolically between 4 points. Warning:
    this interpolator has a zero at Nyquist-frequency which makes it not so suitable for
    non-oversampled audio. */
    static INLINE double getSampleParabolic4p2o2x(double x, double *y);

    /** Calclulates one output sample by interpolating between 4 points with a 3rd order Hermite
    polynomial at position x. The position x should be the fractional part of the actual
    time-instant, the y[1] should be the sample value at the integer part of that time-instant such
    that the *y-array contains two sample-values before and two samples values after the time
    instant, at which the signal is to be read out. */
    static INLINE double getSampleHermite4p3o(double x, double *y);

    /** Calclulates one output sample by interpolating between 4 points with a 2nd-order-osculating
    5th order polynomial. */
    static INLINE double getSample2ndOsculating4p5o(double x, double *y);

    /** Calclulates one output sample by interpolating between 4 points with a 3rd order polynomial,
    optimized for 2x oversampling. */
    static INLINE double getSampleOptimal4p3o2x(double x, double *y);

    /** Calclulates one output sample by interpolating between 6 points with a 5th order
    polynomial, optimized for 2x oversampling. The position x should be the fractional part of the
    actual time-instant, the y[2] should be the sample value at the integer part of that
    time-instant such that the *y-array contains 3 sample-values before and 3 samples values after
    the time instant, at which the signal is to be read out. */
    static INLINE double getSampleOptimal6p5o2x(double x, double *y);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state variable to zero (relevant only for allpass interpolators). */
    void reset();

    //=============================================================================================

  protected:

    double previousOutput;       // previous output sample (for the allpass interpolator)
    int    interpolationMethod;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double Interpolator::getSample(double x, double *y)
  {
    switch( interpolationMethod )
    {
    case NEAREST_NEIGHBOUR:     return getSampleNearestNeighbour(x, y);
    case LINEAR:                return getSampleLinear(x, y);
    case CUBIC_ZERO_DERIVATIVE: return getSampleCubicZeroDerivative(x, y);
    case HERMITE_4P_3O:         return getSampleHermite4p3o(x, y);
    case ALLPASS:               return getSampleAllpass(x, y);
    case WARPED_ALLPASS:        return getSampleWarpedAllpass(x, y);
    case PARABOLIC_4P_2O_2X:    return getSampleParabolic4p2o2x(x, y);

    default:                    return getSampleLinear(x, y);
    }
  }

  INLINE double Interpolator::getSampleAllpass(double x, double *y)
  {
    //double coeff   = 1.0-x;
    double coeff   = x;
    double result  = y[0] + coeff*y[1] - coeff*previousOutput;
    previousOutput = result;
    return result;
  }

  INLINE double Interpolator::getSampleWarpedAllpass(double x, double *y)
  {
    double coeff   = (1.0-x)/(1.0+x);
    //double coeff   = x / (2.0-x); 
    double result  = y[0] + coeff*y[1] - 0.999*coeff*previousOutput; 
      // multiplication of feedback with 0.999... avoids self-oscillation at fs/2
    previousOutput = result;
    return result;
  }

  INLINE double Interpolator::getSampleNearestNeighbour(double x, double *y)
  {
    if( x < 0.5 )
      return y[0];
    else
      return y[1];
  }

  INLINE double Interpolator::getSampleLinear(double x, double *y)
  {
    return y[0] + x*(y[1]-y[0]);
  }

  INLINE double Interpolator::getSampleLinear(double x, float *y)
  {
    return y[0] + x*(y[1]-y[0]);
  }

  INLINE double Interpolator::getSampleCubicZeroDerivative(double x, double *y)
  {
    double d = y[1]-y[0];
    return d*x*x*(3.0-2.0*x) + y[0];
  }

  INLINE double Interpolator::getSampleParabolic4p2o2x(double x, double *y)
  {
    //doubleA y1mym1, c0, c1, c2;

    // 4-point, 2nd-order parabolic 2x (x-form)
    double y1mym1 = y[2]-y[0];
    double c0 = 1/2.0*y[1] + 1/4.0*(y[0]+y[2]);
    double c1 = 1/2.0*y1mym1;
    double c2 = 1/4.0*(y[3]-y[1]-y1mym1);
    return (c2*x+c1)*x+c0;
  }

  INLINE double Interpolator::getSampleHermite4p3o(double x, double *y)
  {
    static doubleA c0, c1, c2, c3;

    // 4-point, 3rd-order Hermite (x-form)
    c0 = y[1];
    c1 = (1.0/2.0)*(y[2]-y[0]);
    c2 = (y[0] - (5.0/2.0)*y[1]) + (2.0*y[2] - (1.0/2.0)*y[3]);
    c3 = (1.0/2.0)*(y[3]-y[0]) + (3.0/2.0)*(y[1]-y[2]);
    return ((c3*x+c2)*x+c1)*x+c0;
  }

  INLINE double Interpolator::getSample2ndOsculating4p5o(double x, double *y)
  {
    static doubleA z, even1, odd1, even2, odd2, c0, c1, c2, c5;

    // 4-point, 5th-order 2nd-order-osculating (z-form)
    z     = x - 1/2.0;
    even1 = y[0]+y[3];
    odd1  = y[0]-y[3];
    even2 = y[1]+y[2];
    odd2  = y[1]-y[2];
    c0    = 9/16.0*even2 - 1/16.0*even1;
    c1    = 3/16.0*odd1 - 25/16.0*odd2;
    c2    = 1/4.0*(even1-even2);
    c5    = odd1 - 3*odd2;
    return (((c5*z*z-c5)*z+c2)*z+c1)*z+c0;
  }

  INLINE double Interpolator::getSampleOptimal4p3o2x(double x, double *y)
  {
    static doubleA z, even1, odd1, even2, odd2, c0, c1, c2, c3;

    // Optimal 2x (4-point, 3rd-order) (z-form)
    z     = x - 1/2.0;
    even1 = y[2]+y[1];
    odd1  = y[2]-y[1];
    even2 = y[3]+y[0];
    odd2  = y[3]-y[0];
    c0    = even1*0.45868970870461956 + even2*0.04131401926395584;
    c1    = odd1*0.48068024766578432 + odd2*0.17577925564495955;
    c2    = even1*-0.246185007019907091 + even2*0.24614027139700284;
    c3    = odd1*-0.36030925263849456 + odd2*0.10174985775982505;
    return ((c3*z+c2)*z+c1)*z+c0;
  }

  INLINE double Interpolator::getSampleOptimal6p5o2x(double x, double *y)
  {
    doubleA z,
      even1, odd1, even2, odd2, even3, odd3,
      c0, c1, c2, c3, c4, c5;

    // Optimal 2x (6-point, 5th-order) (z-form)
    z     = x - 1/2.0;
    even1 = y[3]+y[2];
    odd1  = y[3]-y[2];
    even2 = y[4]+y[1];
    odd2  = y[4]-y[1];
    even3 = y[5]+y[0];
    odd3  = y[5]-y[0];
    c0    = even1*0.40513396007145713 + even2*0.09251794438424393 + even3*0.00234806603570670;
    c1    = odd1*0.28342806338906690 +  odd2*0.21703277024054901 +  odd3*0.01309294748731515;
    c2    = even1*-0.191337682540351941 + even2*0.16187844487943592 + even3*0.02946017143111912;
    c3    = odd1*-0.16471626190554542 + odd2*-0.00154547203542499 + odd3*0.03399271444851909;
    c4    = even1*0.03845798729588149 + even2*-0.05712936104242644 + even3*0.01866750929921070;
    c5    = odd1*0.04317950185225609 + odd2*-0.01802814255926417 + odd3*0.00152170021558204;
    return ((((c5*z+c4)*z+c3)*z+c2)*z+c1)*z+c0;
  }

}  // end namespace rosic

#endif // rosic_Interpolator_h
