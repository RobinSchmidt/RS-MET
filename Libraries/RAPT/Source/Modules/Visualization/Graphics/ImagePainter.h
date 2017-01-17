#ifndef RAPT_IMAGEPAINTER_H_INCLUDED
#define RAPT_IMAGEPAINTER_H_INCLUDED

/** A class for painting on an Image object. It is based on a second image that we use as "brush".
Whenever a "dot" is painted onto the image at a particular location, a copy of the brush image will
be blended into the target image at that location. The pixel types for the target image and the 
brush image can in general be different - for example the target image could use a RBGA type 
represented by a vector of 4 floats and the brush could use just a single float. The brush pixels
will be used as weights by which a new color will be multiplied when it is accumulated into the 
image. */

template<class TPix, class TWgt, class TCor>  // pixel, weight, coordinate types
class ImagePainter
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  ImagePainter(Image<TPix> *imageToPaintOn = nullptr, Image<TWgt> *brushToUse = nullptr);


  /** \name Setup */

  /** Sets the image on which we paint.  */
  void setImageToPaintOn(Image<TPix> *imageToPaintOn);

  /** Sets the brush that we use. It is basically a matrix of weights. */
  void setBrushToUse(Image<TWgt> *brushToUse);



  /** \name Inquiry */



  /** \name Painting */

  void paintDot(TCor x, TCor y, TPix color);



protected:

  // data members:

  Image<TPix> *image;
  Image<TWgt> *brush;

};

#endif