#ifndef RAPT_BELLFUNCTIONS_H_INCLUDED
#define RAPT_BELLFUNCTIONS_H_INCLUDED

/** This class is a collection of bell shaped functions where we assume that the input is
nonnegative: x >= 0. If your input value may be negative, you can just pass the absolute value, 
such that the positive side will be mirrored at the y-axis. The functions don't do that 
themselves for possible effciency gains, when you already know that the input is nonnegative. */

template<class T>
class rsPositiveBellFunctions
{

public:

  /** Ramps down from 1 to 0 in a straight line such that f(x = 0) = 1, f(x >= 1) = 0. */
  static T linear(T x);

  /** Uses a cubic polynomial to smoothly ramp down from 1 to 0 such that  f(x = 0) = 1, 
  f(x >= 1) = 0. The polynomial between 0 <= x <= 1 is chosen such that f'(0) = 0 and f'(1) = 0. */
  static T cubic(T x);

  /** Like cubic, but uses a quintic (5th order) polynomial that additionally satisfies  
  f''(0) = 0 and f''(1) = 0. */
  static T quintic(T x);

  /** Like quintic, but uses a heptic (7th order) polynomial that additionally satisfies  
  f'''(0) = 0 and f'''(1) = 0. */
  static T heptic(T x);

  /** Implements the one-dimensional bump function which is both smooth and compactly supported. 
  Smooth means that is infinitely often differentiable. Compactly supported meas that it is nonzero
  only inside a finite interval (here in -1..+1). It's defined piecewise as f(x) = exp(-1/(1-x^2)) 
  for x in (-1,1) and 0 outside (-1,1). see:
  https://en.wikipedia.org/wiki/Bump_function 
  https://en.wikipedia.org/wiki/Distribution_(mathematics)#Basic_idea  */
  static T bump(T x);

  /** Generalization of the bump function that replaces x^2 by |x|^p containing a parameter p, i.e.
  f(x) = exp(-1/(1-x^2)) by f(x) = exp(-1/(1-|x|^p)) in the interval (-1,1). */
  static T bump(T x, T p);
  // ToDo: 
  // -Document what the parameter p does intuitively. It somehow controls the shape/transition. 
  // -Figure out and document, for what values of p the function f is still smooth at -1 and +1.

};

//=================================================================================================

/** This is a class for representing and evaluating bell shaped functions with adjustable center,
width and (relative) length of an inserted flat top zone. Picture this as cutting the bell into two
halves and dragging them apart while squeezing them along the x-axis (such that the overall width
stays the same) and inserting a constant value in between. The relative width of this flat top zone
is controlled by a parameter between 0 and 1 where 0 means no flat top zone and 1 reduces the whole
bell to a rectangular pulse. The shape of the transition is based on a prototype bell shaped 
function that you can pass as function pointer (suitable functions are those from the class 
rsPositiveBellFunctions, but you may also use others as well). This class here will then scale and 
shift the input according to your desired center value, width and flat-top length.  */

template<class T>
class rsParametricBellFunction
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsParametricBellFunction();


  /** \name Setup */

  /** Sets the center value for the function. */
  void setCenter(T newCenter);

  /** Sets the width value for the function. */
  void setWidth(T newWidth);

  /** Sets the relative width of the flat top. If 0, there will be no flat top at all. If 1, the
  whole bell will reduce to a rectangular pulse between center-width/2...center+width/2. */
  void setFlatTopWidth(T newWidth);

  /** Sets the prototype function to use. This prototype is supposed to represent a zero-centered
  bell shaped function. It is sufficient, if this function is defined for nonnegative input values
  only, because symmetrization is done internally in this class anyway. */
  void setPrototypeBell(T (*newFunction)(T));


  /** \name Function Evaluation */

  /** Returns an output value for given x. */
  inline T getValue(T x);


protected:


  /** \name Data */

  // user parameters:
  T center, flat;     // center and relative length of flat zone
  T (*bell)(T);       // prototype bell function

  // internal coefficients:
  T a, b;

};

template<class T>
inline T rsParametricBellFunction<T>::getValue(T x)
{
  T tmp = rsAbs(a*(x-center));
  //double tmp = fabs(a*(x-center));
  if(tmp <= flat)
    return bell(0.0);
  else
    return bell(b*(tmp-flat));
}

#endif
