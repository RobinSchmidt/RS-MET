#ifndef RAPT_IMAGEPROCESSOR_H_INCLUDED
#define RAPT_IMAGEPROCESSOR_H_INCLUDED


// maybe rename class rsImageProcessor to rsImageEffect and/or rsImageTransform and/or rename the
// file to ImageProcessors.h/cpp (plural, because we have several class in these files)

/** A collection of (basic) image processing algorithms. */

template<class T>
class rsImageProcessor
{

public:

  /** Applies the exponent "gamma" to all pixel values: y = x^gamma. */
  static void gammaCorrection(rsImage<T>& img, T gamma);
  // todo: check conventions used in image processing (i think, IrfanView's gamma correction uses
  // the reciprocal)

  /** Inverts the brightness values of all pixels in the given image */
  static void invert(rsImage<T>& img);

  /** Normalizes the range of the pixel values to be 0..1. */
  static void normalize(rsImage<T>& img);

  /** Faster but less numerically precise implementation of normalization. */
  static void normalizeFast(rsImage<T>& img);

  /** Joint normalization of two images. */
  static void normalizeJointly(rsImage<T>& img1, rsImage<T>& img2);

  /** Scales the image up by the given factor by simply repeating pixel values. */
  static rsImage<T> scaleUp(const rsImage<T>& img, int scl);

  /** Shapes a ramp for 0 to 1 into a smooth sine curve. */
  static void sineShape(rsImage<T>& img);



  // todo:
  // move more code from GraphicsExperiments to here
  // static void filter(const rsImage<T>& input, const rsImage<T>& kernel);
  // todo: have options for left/right/top/bottom border handling - maybe make an enum
  // borderHandling with values: zero, repeat, mirror, periodic

};



//=================================================================================================

/** A class for contour plotting from an input image. The input image represents height-levels and 
the generated outputs are images with either the contour lines or contour fills. */

template<class TPix, class TLvl>
class rsImageContourPlotter
{

public:


  rsImage<TPix> getContourLines(const rsImage<TPix>& z, const std::vector<TLvl>& levels, 
    const std::vector<TPix>& colors, bool antiAlias);

  rsImage<TPix> getContourFills(const rsImage<TPix>& z, const std::vector<TLvl>& levels,
    const std::vector<TPix>& colors, bool antiAlias);


protected:

  // internal sub-routines (maybe make some of them public)



  static void drawContour(const rsImage<TLvl>& z, TLvl level, rsImage<TPix>& target, TPix color, 
    bool antiAlias);

  static void fillBetweenContours(const rsImage<TLvl>& z, TLvl lo, TLvl hi, rsImage<TPix>& target,
    TPix fillColor, bool antiAlias = false);


  /** Used in drawContour for anti-aliasing. */
  static void contourSubPixelPosition(TLvl z00, TLvl z01, TLvl z10, TLvl z11, TLvl c,
    TLvl* x, TLvl* y, TLvl* weight);
  // x,y may actually be another type (coordinates), and weight may be yet another

  /** Used in fillBetweenContours for anti-aliasing. */
  static TLvl contourPixelCoverage(TLvl z00, TLvl z01, TLvl z10, TLvl z11, TLvl c);

  /** Used in contourSubPixelPosition and contourPixelCoverage. */
  static int contourSegmentCoeffs(TLvl z00, TLvl z01, TLvl z10, TLvl z11, TLvl c,
    TLvl& x0, TLvl& y0, TLvl& x1, TLvl& y1);



  //rsImagePainter<TPix, TLvl, TLvl> painter;



};



#endif