#ifndef RAPT_IMAGEPAINTER_H_INCLUDED
#define RAPT_IMAGEPAINTER_H_INCLUDED

/**  */

template<class TPix, class TWgt, class TCor>  // pixel, weight, coordinate types
class ImagePainter
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  ImagePainter(Image<TPix> *imageToPaintOn = nullptr, Image<TWgt> *brushToUse = nullptr);




  /** \name Setup */

  /**  */
  void setImageToPaintOn(Image<TPix> *imageToPaintOn);


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