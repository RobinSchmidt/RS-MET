#ifndef RS_BELLFUNCTIONS_H
#define RS_BELLFUNCTIONS_H

namespace RSLib
{

/** This class is a collection of bell shaped functions where we assume that the input is
nonnegative: x >= 0. If your input value may be negative, you can just pass the absolute value, 
such that the positive side will be mirrored at the y-axis. The functions don't do that 
themselves for possible effciency gains, when you already know that the input is nonnegative. */

class RSLib_API rsPositiveBellFunctions
{

public:

  /** Ramps down from 1 to 0 in a straight line such that f(x = 0) = 1, f(x >= 1) = 0. */
  static double linear(double x);

  /** Uses a cubic polynomial to smoothly ramp down from 1 to 0 such that  f(x = 0) = 1, 
  f(x >= 1) = 0. The polynomial between 0 <= x <= 1 is chosen such that f'(0) = 0 and f'(1) = 0. */
  static double cubic(double x);

  /** Like cubic, but uses a quintic (5th order) polynomial that additionally satisfies  
  f''(0) = 0 and f''(1) = 0. */
  static double quintic(double x);

  /** Like quintic, but uses a heptic (7th order) polynomial that additionally satisfies  
  f'''(0) = 0 and f'''(1) = 0. */
  static double heptic(double x);

};

//=================================================================================================

/** This is a class for representing and evaluating bell shaped functions with adjustable center,
width and (relative) length of an inserted flat top zone. Picture this as cutting the bell into two
halves and dragging them apart while squeezing them along the x-axis (such that the overall width
stays the same) and inserting a constant value in between. The relative width of this flat top zone
is controlled by a parameter between 0 and 1 where 0 means no flat top zone and 1 reduces the whole
bell to a rectangular pulse. The shape of the transition is are based on a prototype bell shaped 
function that you can pass as function pointer (suitable functions are those from the class 
rsPositiveBellFunctions, but you may also use others as well). This class here will then scale and 
shift the input according to your desired center value, width and flat-top length.  */

class RSLib_API rsParametricBellFunction
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsParametricBellFunction();


  /** \name Setup */

  /** Sets the center value for the function. */
  void setCenter(double newCenter);

  /** Sets the width value for the function. */
  void setWidth(double newWidth);

  /** Sets the relative width of the flat top. If 0, there will be no flat top at all. If 1, the
  whole bell will reduce to a rectangular pulse between center-width/2...center+width/2. */
  void setFlatTopWidth(double newWidth);

  /** Sets the prototype function to use. This prototype is supposed to represent a zero-centered
  bell shaped function. It is sufficient, if this function is defined for nonnegative input values
  only, because symmetrization is done internally in this class anyway. */
  void setPrototypeBell(double (*newFunction)(double));


  /** \name Function Evaluation */

  /** Returns an output value for given x. */
  RS_INLINE double getValue(double x);


protected:


  /** \name Data */

  // user parameters:
  double center, width, flat;  // center, width and relative length of flat zone
  double (*bell)(double);      // prototype bell function

  // internal coefficients:
  double a, b;

  // width, flat are actually redundant now

};

RS_INLINE double rsParametricBellFunction::getValue(double x)
{
  double tmp = rsAbs(a*(x-center));
  if(tmp < flat)
    return bell(0.0);
  else
    return bell(b*(tmp-flat));



  //double tmp = rsAbs(2*(x-center) / width);
  //if(tmp < flat)
  //  return bell(0.0);
  //else
  //  return bell((tmp-flat) / (1-flat));



  // \todo: precompute 2/width, 1/(1-flat) to avoid divisions

}

}

#endif
