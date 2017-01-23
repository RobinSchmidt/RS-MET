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

};

#endif