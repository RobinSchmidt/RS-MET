#ifndef RAPT_IMAGEPAINTER_H_INCLUDED
#define RAPT_IMAGEPAINTER_H_INCLUDED

/** A class for painting on an Image object. It is based on an "alpha-mask" that is used as 
prototype "dot". Whenever a dot is painted onto the image at a particular location, the mask will
be used to blend the existing colors at these pixels with a new target color.  */

template<class TPix, class TWgt, class TCor>  // pixel, weight, coordinate types
class ImagePainter
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  ImagePainter(Image<TPix> *imageToPaintOn, AlphaMask<TWgt> *maskToUse);


  /** \name Setup */

  /** Sets the image on which we paint.  */
  void setImageToPaintOn(Image<TPix> *imageToPaintOn);

  /** Sets the alpha mask that we use as prototye "dot". It is basically a matrix of weights. */
  void setAlphaMaskForDot(AlphaMask<TWgt> *maskToUse);


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
  AlphaMask<TWgt> *mask;

  int wi, hi;    // image width and height

};

#endif