#ifndef RAPT_IMAGEPROCESSOR_H_INCLUDED
#define RAPT_IMAGEPROCESSOR_H_INCLUDED

/** A collection of (basic) image processing algorithms. */

template<class T>
class rsImageProcessor
{

public:


  static void gammaCorrection(rsImage<T>& img, T gamma);

  /** Inverts the brightness values of all pixels in the given image */
  static void invert(rsImage<T>& img);



  static void normalize(rsImage<T>& img);


  static void normalizeFast(rsImage<T>& img);



  /** Joint normalization of two images */
  static void normalizeJointly(rsImage<T>& img1, rsImage<T>& img2);





  /** Scales the image up by the given factor by simply repeating pixel values */
  static rsImage<T> scaleUp(const rsImage<T>& img, int scl);


  /** Shapes a ramp for 0 to 1 into a smooth sine curve. */
  static void sineShape(rsImage<T>& img);

  // move code from GraphicsExperiments to here



};


#endif