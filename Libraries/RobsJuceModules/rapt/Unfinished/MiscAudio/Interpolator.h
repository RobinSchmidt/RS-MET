#ifndef RAPT_INTERPOLATOR_H
#define RAPT_INTERPOLATOR_H

/** This class implements some interpolation formulas to reconstruct a continuous signal from 
samples, that is: it allows for reading out the signal at arbitrary non-integer time instants. 
Most of the member-functions which perform the interpolation are declared as static such that one 
does not need to actually create an instance of the class if one wants to use them. However, some 
other interpolation-methods (namely the allpass-methods) need to keep track of internal state
variables - these are handled here as member variables, thus those interpolation-methods can not
work as static methods. So if you want to use them, you will need to instantiate an object of
the Interpolator class.

References:
 -Olli Niemitalo: Polynomial Interpolators for High-Quality Resampling of Oversampled Audio
 -Jon Datorro: Effect Design Part 2 - Delay-Line Modulation and Chorus (JAES, Vol 45, No. 10)

\todo: avoid recalculation of the warped allpas coefficient when it stays the same between
  successive calls */

template<class T>  // maybe we will need to distinguish between position/index type and value type (TIdx or TPos, TVal)
class rsInterpolator
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


  /** \name Construction/Destruction */

  /** Constructor. */
  rsInterpolator()
  {
    previousOutput      = 0;
    interpolationMethod = LINEAR;
  }


  ///** Destructor. */
  //~rsInterpolator();


  /** \name Setup */

  /** Chooses one of the interpolation methods. @see enum interpolationMethods. */
  RS_INLINE void setInterpolationMethod(int newMethod)
  {
    interpolationMethod = newMethod;
  }


  /** \name Inquiry */

  /** Returns the currently chosen interpolation method. @see enum interpolationMethods. */
  int getInterpolationMethod() const { return interpolationMethod; }


  /** \name Audio Processing */

  /** Calculates one interpolated output-sample. It will use  the interpolation-method that has
  been selected via the setInterpolationMethod()-function. */
  RS_INLINE T getSample(T x, T *y);

  /** Calculates an output sample using allpass interpolation as described in Datorro's paper -
  note that the x-value used here as input argument corresponds to 1-frac in the papaer
  because the Interpolator class sees it from the viewpoint of a phase-accumulating oscillator
  whereas the paper deals with delaylines. */
  RS_INLINE T getSampleAllpass(T x, T *y);

  /** Calculates an output sample using warped allpass interpolation as described in Datorro's
  paper - note that the x-value used here as input argument corresponds to 1-frac in the papaer
  because the Interpolator class sees it from the viewpoint of a phase-accumulating oscillator
  whereas the paper deals with delaylines. */
  RS_INLINE T getSampleWarpedAllpass(T x, T *y);

  /** Returns the nearest neighbour sample - that is y[0] when x < 0.5 and y[1] otherwise. */
  static RS_INLINE T getSampleNearestNeighbour(T x, T *y);

  /** Calclulates one output sample by linearly interpolating between y[0] and y[1] at
  position x. y should point to an array of at least two values, x should be a number
  between 0...1. */
  static RS_INLINE T getSampleLinear(T x, T *y);





  static RS_INLINE T getSampleCubicZeroDerivative(T x, T *y);

  /** Calclulates one output sample by interpolating parabolically between 4 points. Warning:
  this interpolator has a zero at Nyquist-frequency which makes it not so suitable for
  non-oversampled audio. */
  static RS_INLINE T getSampleParabolic4p2o2x(T x, T *y);

  /** Calclulates one output sample by interpolating between 4 points with a 3rd order Hermite
  polynomial at position x. The position x should be the fractional part of the actual
  time-instant, the y[1] should be the sample value at the integer part of that time-instant such
  that the *y-array contains two sample-values before and two samples values after the time
  instant, at which the signal is to be read out. */
  static RS_INLINE T getSampleHermite4p3o(T x, T *y);

  /** Calclulates one output sample by interpolating between 4 points with a 2nd-order-osculating
  5th order polynomial. */
  static RS_INLINE T getSample2ndOsculating4p5o(T x, T *y);

  /** Calclulates one output sample by interpolating between 4 points with a 3rd order polynomial,
  optimized for 2x oversampling. */
  static RS_INLINE T getSampleOptimal4p3o2x(T x, T *y);

  /** Calclulates one output sample by interpolating between 6 points with a 5th order
  polynomial, optimized for 2x oversampling. The position x should be the fractional part of the
  actual time-instant, the y[2] should be the sample value at the integer part of that
  time-instant such that the *y-array contains 3 sample-values before and 3 samples values after
  the time instant, at which the signal is to be read out. */
  static RS_INLINE T getSampleOptimal6p5o2x(T x, T *y);


  /** \name Misc */

  /** Resets the internal state variable to zero (relevant only for allpass interpolators). */
  RS_INLINE void reset() { previousOutput = 0; }


protected:

  T previousOutput;       // previous output sample (for the allpass interpolator)
  int    interpolationMethod;

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class T>
RS_INLINE T rsInterpolator<T>::getSample(T x, T *y)
{
  switch(interpolationMethod)
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

template<class T>
RS_INLINE T rsInterpolator<T>::getSampleAllpass(T x, T *y)
{
  T result  = y[0] + x*y[1] - 0.999*x*previousOutput;
  previousOutput = result;
  return result;
  // multiplication of feedback with 0.999 avoids self-oscillation at fs/2
}

template<class T>
RS_INLINE T rsInterpolator<T>::getSampleWarpedAllpass(T x, T *y)
{
  T coeff   = (1.0-x)/(1.0+x);
  //T coeff   = x / (2.0-x); 
  T result  = y[0] + coeff*y[1] - 0.999*coeff*previousOutput;
    // multiplication of feedback with 0.999... avoids self-oscillation at fs/2
  previousOutput = result;
  return result;
}

template<class T>
RS_INLINE T rsInterpolator<T>::getSampleNearestNeighbour(T x, T *y)
{
  if(x < 0.5)
    return y[0];
  else
    return y[1];
}

template<class T>
RS_INLINE T rsInterpolator<T>::getSampleLinear(T x, T *y)
{
  return y[0] + x*(y[1]-y[0]);
}

template<class T>
RS_INLINE T rsInterpolator<T>::getSampleCubicZeroDerivative(T x, T *y)
{
  T d = y[1]-y[0];
  return d*x*x*(3.0-2.0*x) + y[0];
}

template<class T>
RS_INLINE T rsInterpolator<T>::getSampleParabolic4p2o2x(T x, T *y)
{
  //TA y1mym1, c0, c1, c2;

  // 4-point, 2nd-order parabolic 2x (x-form)
  T y1mym1 = y[2]-y[0];
  T c0 = 1/2.0*y[1] + 1/4.0*(y[0]+y[2]);
  T c1 = 1/2.0*y1mym1;
  T c2 = 1/4.0*(y[3]-y[1]-y1mym1);
  return (c2*x+c1)*x+c0;
}

template<class T>
RS_INLINE T rsInterpolator<T>::getSampleHermite4p3o(T x, T *y)
{
  static T c0, c1, c2, c3;

  // 4-point, 3rd-order Hermite (x-form)
  c0 = y[1];
  c1 = (1.0/2.0)*(y[2]-y[0]);
  c2 = (y[0] - (5.0/2.0)*y[1]) + (2.0*y[2] - (1.0/2.0)*y[3]);
  c3 = (1.0/2.0)*(y[3]-y[0]) + (3.0/2.0)*(y[1]-y[2]);
  return ((c3*x+c2)*x+c1)*x+c0;
}

template<class T>
RS_INLINE T rsInterpolator<T>::getSample2ndOsculating4p5o(T x, T *y)
{
  static T z, even1, odd1, even2, odd2, c0, c1, c2, c5;

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

template<class T>
RS_INLINE T rsInterpolator<T>::getSampleOptimal4p3o2x(T x, T *y)
{
  static T z, even1, odd1, even2, odd2, c0, c1, c2, c3;

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

template<class T>
RS_INLINE T rsInterpolator<T>::getSampleOptimal6p5o2x(T x, T *y)
{
  T z,
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


#endif
