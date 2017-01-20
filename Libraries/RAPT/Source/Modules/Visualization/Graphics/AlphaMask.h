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

  //void setMaxPixelSize(int newMaxWidth, int newMaxHeight);

  void setSize(double newSize);

  //void setShape(int newShape);
    // circle, rectangle, etc


protected:

  void renderMask();

  double size;

};

#endif