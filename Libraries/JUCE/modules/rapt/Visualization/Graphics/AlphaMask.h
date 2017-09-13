#ifndef RAPT_ALPHAMASK_H_INCLUDED
#define RAPT_ALPHAMASK_H_INCLUDED

/** This is a subclass of Image intended to be used as alpha mask for blending pixels in an 
existing image with new colors. Your pixel type should be a monochromatic type, i.e. a scalar
such as float or double. */

template<class TPix>  // pixel type
class AlphaMask : public rsImageResizable<TPix>
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

  /** Sets the slope of the transition at the inside of the flat-color zone. */
  void setInnerSlope(double newSlope);

  /** Sets the slope of the transition at the outside of the dot. */
  void setOuterSlope(double newSlope);

  //void setShape(int newShape);
    // circle, rectangle, etc

  /** Copies the user parameters which determine the shape (transition width, slopes, etc.) from
  some other mask. Basically, that's all user parameters except the size. It's intended to be used
  for creating a preview of a mask, where the preview does not necessarily need to be of the same 
  size as the actual mask used for painting. The other mask may use a different template parameter
  for the pixels. */
  template<class T>
  void copyShapeParametersFrom(const AlphaMask<T>& otherMask);

  /** \name Inquiry */

  inline double getTransitionWidth() const { return flat; }
  inline double getInnerSlope()      const { return slope0; }
  inline double getOuterSlope()      const { return slope1; }


  /** \name Misc */

  static double cubicBell(double x, double steepnessAt0, double steepnessAt1);
    // move to Math

  //static double rationalMap(double x, double a);

protected:

  /** Renders the mask. Called internally whenever a parameter chenges. */
  void renderMask();

  /** Given a (normalized) distance from the center, this function returns the alpha value that 
  should be rendered into the mask at this distance. */
  double getAlphaForDistance(double distance);

  // data:
  double flat;           // flat top radius
  double slope0, slope1; // (negative) slope at center (0) and border (1)
};

#endif