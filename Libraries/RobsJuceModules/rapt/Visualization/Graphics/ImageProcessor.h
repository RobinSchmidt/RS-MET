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
  // the reciprocal - check what is the most widely used convention). maybe rename to just gamma

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

  // todo: flipLeftRight, flipUpDown, transpose





  // todo:
  // move more code from GraphicsExperiments to here
  // static void filter(const rsImage<T>& input, const rsImage<T>& kernel);
  // todo: have options for left/right/top/bottom border handling - maybe make an enum
  // borderHandling with values: zero, repeat, mirror, periodic

};



//=================================================================================================

/** A class for contour plotting from an input image. The input image represents height-levels and 
the generated outputs are images with either the contour lines or contour fills. */

template<class TPix, class TVal>  // pixel and value types (for coordinates, heights, etc.)
class rsImageContourPlotter
{

public:


  rsImageContourPlotter();


  //-----------------------------------------------------------------------------------------------
  // \name Contour Drawing

  /** For a given input image that is interpreted as containing values of a height function 
  z = f(x,y), this function returns an image that contains the contour lines for the given set of
  contour levels. You may also give the lines different colors...  */
  rsImage<TPix> getContourLines(const rsImage<TPix>& z, const std::vector<TVal>& levels, 
    const std::vector<TPix>& colors, bool antiAlias);

  /** Like getContourLines, but instead of drawing the contour-lines themselves, it fills the areas
  between the contour lines with given colors. */
  rsImage<TPix> getContourFills(const rsImage<TPix>& z, const std::vector<TVal>& levels,
    const std::vector<TPix>& colors, bool antiAlias);


  //-----------------------------------------------------------------------------------------------
  // \name Internal Subroutines

  void drawContour(const rsImage<TVal>& z, TVal level, rsImage<TPix>& target, TPix color, 
    bool antiAlias = true);

  void fillBetweenContours(const rsImage<TVal>& z, TVal lo, TVal hi, rsImage<TPix>& target,
    TPix fillColor, bool antiAlias = true);


  /** Used in drawContour for anti-aliasing. */
  void contourSubPixelPosition(TVal z00, TVal z01, TVal z10, TVal z11, TVal c,
    TVal* x, TVal* y, TVal* weight);
  // x,y may actually be another type (coordinates), and weight may be yet another
  // or maybe use som generic TVal for "value"

  /** Used in fillBetweenContours for anti-aliasing. */
  TVal contourPixelCoverage(TVal z00, TVal z01, TVal z10, TVal z11, TVal c);

  /** Used in contourSubPixelPosition and contourPixelCoverage.
  Given values at a pixel z00 and its neighbours, this function computes the coeffs of a 
  parametric line equation x(t) = x0 + t*(x1-x0), y(t) = y0 + t*(y1-y0). The line is supposed to 
  approximate a segment of a contour of a function/image with level given by c. This is how the 
  function values are distributed over the pixels: z00 is the value at the pixel under 
  investigation and z01, z10, z11 are bottom, right and diagonal neighbours (in that order):
     z00--z10
      |    |
     z01--z11  
  so the first index refers to the x-coordinate, the second to the y-coordinate and y increases 
  downward as is common in image processing when coordinates are expressed in pixel-coordinates. 
  The function also returns an integer indicating which branch we ended up in: 
    0: top-left, 1: top-right, 2: bottom-left, 3: bottom-right, 4: horizontalish, 5: verticalish 
  where "top-left" etc. means that the segment "cuts off" the top-left corner of the pixel. */
  int contourSegmentCoeffs(TVal z00, TVal z01, TVal z10, TVal z11, TVal c,
    TVal& x0, TVal& y0, TVal& x1, TVal& y1);



  rsImagePainter<TPix, TVal, TVal> painter;

};

#endif