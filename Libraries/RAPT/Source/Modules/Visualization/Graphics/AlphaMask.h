#ifndef RAPT_ALPHAMASK_H_INCLUDED
#define RAPT_ALPHAMASK_H_INCLUDED

/** This is a subclass of Image intended to be used as alpha mask for blending pixels in an 
existing image with new colors. Your pixel type should be a monochromatic type, i.e. a scalar
such as float or double. */

template<class TPix>  // pixel type
class AlphaMask : public ImageResizable<TPix>
{

public:


  /** \name Construction/Destruction */

  /** Constructor. */
  AlphaMask();


  /** \name Setup */

  /** Sets the size/diameter of the mask in pixels. */
  void setSize(double newSize);

  /** Sets the width of the transition from full alpha to zero as value between 0 and 1. 0 means a 
  hard transition, 1 a maximally soft transition. */
  void setTransitionWidth(double newWidth);

  //void setShape(int newShape);
    // circle, rectangle, etc


protected:

  void renderMask();

  /*double size, transitionWidth;*/

  rsParametricBellFunction<double> bell;
  // instead of the parameteric bell, use a 3rd order polynomial with adjustable slopes at start-
  // and endpoint.
  // f(x) = a0 + a1*x + a2*x^2 + a3*x^3 with f(0) = 1, f'(0) = s0, f(1) = 0, f'(1) = s1
  // solving it yields: a0 = 1, a1 = s0, a2 = -s1 - 2*s0 - 3, a3 = 2 + s0 + s1
  // maybe the user parameter should be the negative slope
  // we should figure out the condition for not having a local minimum or maximum between 0..1 and 
  // perhaps restrict the parameter range

};

#endif