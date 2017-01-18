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
  ImagePainter(Image<TPix> *imageToPaintOn, Image<TWgt> *brushToUse);


  /** \name Setup */

  /** Sets the image on which we paint.  */
  void setImageToPaintOn(Image<TPix> *imageToPaintOn);

  /** Sets the brush that we use. It is basically a matrix of weights. */
  void setBrushToUse(Image<TWgt> *brushToUse);


  /** \name Inquiry */


  /** \name Painting */

  /** Function for painting a simple 3x3 dot at given integer position. */
  void paintDot3x3(int x, int y, TPix color, TWgt weightStraight = 0, TWgt weightDiagonal = 0);  

  /** Function for painting a simple 3x3 dot at given noninteger position. */
  void paintDot3x3(TCor x, TCor y, TPix color, TWgt weightStraight = 0, TWgt weightDiagonal = 0);

  /** Paints a dot at an integer position using out stored brush (which represents a prototype 
  dot). */
  void paintDot(int x, int y, TPix color);

  void paintDot(TCor x, TCor y, TPix color);


protected:

  /** Internal functions. */

  /** Accumulates the given value into the accumulator accu. We use a rather peculiar accumulation
  function here: newAccu = (oldAccu + value) / (1 + value). When accu starts out a zero and all 
  accumulated values are >= 0, this function will ensure that accu is always < 1 and it will go
  into saturation smoothly. */
  inline void accumulate(TPix &accu, TPix value)
  {
    accu = (accu + value) / (TPix(1) + value);
  }
  // rename to addAndSaturate

  // data members:

  Image<TPix> *image;
  Image<TWgt> *brush;
  //ImageBrush<TWgt, TCor> *brush;




  int wi, hi;    // image width and height
  int wb, hb;    // brush width and height

};

#endif