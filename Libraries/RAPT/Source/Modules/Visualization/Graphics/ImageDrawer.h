#ifndef RAPT_IMAGEDRAWER_H_INCLUDED
#define RAPT_IMAGEDRAWER_H_INCLUDED

/** A baseclass for drawing on images. It consolidates the data and functionality that all drawers 
have in common, regardless of what they draw. Subclasses like LineDrawer will do the actual drawing
and likely define more data and methods. */

template<class TPix, class TWgt, class TCor>  // pixel, weight, coordinate types
class ImageDrawer
{

public:

  /** The blend modes for determining a new pixel color as function of its current color, an 
  incoming color and a blend amount. */
  enum blendModes
  {
    BLEND_LINEAR = 0,
    BLEND_ADD_CLIP,
    BLEND_ADD_SATURATE,
    //BLEND_MULTIPLY,

    NUM_BLEND_MODES
  };


  /** \name Construction/Destruction */

  /** Constructor. */
  ImageDrawer(Image<TPix> *imageToDrawOn);


  /** \name Setup */

  /** Sets the image taht we will draw on. */
  void setImageToDrawOn(Image<TPix> *imageToDrawOn);

  /** Selects one of the blend modes. The blend mode is the function that is used to compute a new
  color for a pixel from an incoming desired color, the pixel's old color and a weight between
  0 and 1 that determines how to mix the old and the new color. */
  void setBlendMode(int newMode);



  /** \name Blend functions */

  static void linearBlend(TPix &pixel, TPix color, TWgt blend);

  static void addAndClip(TPix &pixel, TPix color, TWgt blend);

  static void addAndSaturate(TPix &pixel, TPix color, TWgt blend);



protected:

  Image<TPix> *image;

  int blendMode;
  void (*blendFunction)(TPix& pixel, TPix color, TWgt weight);

};

//=================================================================================================

/** ...add the LineDrawer class here... */

#endif