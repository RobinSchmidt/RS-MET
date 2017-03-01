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

  /** Sets the image that we will draw on. */
  void setImageToDrawOn(Image<TPix> *imageToDrawOn);

  /** Sets the color with which we draw on the image. */
  inline void setColor(TPix newColor) { color = newColor; }

  /** Selects one of the blend modes. The blend mode is the function that is used to compute a new
  color for a pixel from an incoming desired color, the pixel's old color and a weight between
  0 and 1 that determines how to mix the old and the new color. */
  void setBlendMode(int newMode);
    // maybe make it possible to provide a function pointer to a custom blend function 


  /** \name Drawing */

  /** Blends the pixel color in the image at given coordinates with the color of this drawer 
  according to some weight. If the weight is 0, the pixel's color is unchanged, if it's 1, the new 
  color has the biggest impact.  What exactly happens depends on the blend-mode setting. */
  inline void plot(int x, int y, TWgt weight)
  {
    blendFunction((*image)(x, y), color, TPix(weight));
  }


protected:

  Image<TPix> *image;

  TPix color;
  int blendMode;
  void (*blendFunction)(TPix& pixel, TPix color, TWgt weight);


  /** \name Blend functions */

  static void linearBlend(   TPix &pixel, TPix color, TWgt blend);
  static void addAndClip(    TPix &pixel, TPix color, TWgt blend);
  static void addAndSaturate(TPix &pixel, TPix color, TWgt blend);

};

//=================================================================================================

/** A class for drawing straight lines. The lines do not necessarily have to have a solid color -
instead, you can choose one of the line profiles that lets the color vary in dependence on the
distance of a pixel from the ideal geometric line.  */

template<class TPix, class TWgt, class TCor>  // pixel, weight, coordinate types
class LineDrawer : public ImageDrawer<TPix, TWgt, TCor>
{

public:


  /** The blend modes for determining a new pixel color as function of its current color, an 
  incoming color and a blend amount. */
  enum lineProfiles
  {
    PROFILE_FLAT = 0,      // solid/flat color
    PROFILE_LINEAR,        // metallic
    PROFILE_PARABOLIC,     // plastic
    PROFILE_CUBIC,         // cloudy
    //PROFILE_PARAMETRIC,  // add later: define solid-width and 2 slopes (like in the dot)

    NUM_LINE_PROFILES
  };


  /** Constructor. */
  LineDrawer(Image<TPix> *imageToDrawOn);


  /** \name Setup */

  /** Selects one of the line profiles. */
  void setLineProfile(int newProfile);

  /** Sets the width of the line in pixels. */
  void setLineWidth(TCor newWidth);

  /** Sets the drawer up to draw half-circular end caps, if true. If false, the caps will be 
  rectangular. */
  void setRoundCaps(bool shouldBeRound);


  /** \name Drawing */

  /** Draws a line from (x0,y0) to (x1,y1). */
  void drawLine(TCor x0, TCor y0, TCor x1, TCor y1);

  /** Special line drawing function that is supposed to be used for drawing sequences of connected
  lines. After an initial call to drawLine or a previous call to lineTo, you can call this 
  function in order to avoid artifacts (phantom circles) at the line joints. */ 
  void lineTo(TCor x1, TCor y1);


  //void drawConnectedLine(TCor x0, TCor y0, TCor x1, TCor y1);



  // todo: have members for simplified 1-pixel line drawing: drawLineWu, drawLineBresenham


protected:

  bool roundCaps = true;
  TCor w2;                // lineWidth/2
  TCor x0 = 0, y0 = 0;    // start-point for lineTo function
  int  profileIndex;
  TWgt (*lineProfile)(TCor distance, TCor halfWidth);

  //vector<Line2D> lines;


  /** \name Line profile functions */

  static TWgt profileFlat(     TCor distance, TCor halfWidth);
  static TWgt profileLinear(   TCor distance, TCor halfWidth);
  static TWgt profileParabolic(TCor distance, TCor halfWidth);
  static TWgt profileCubic(    TCor distance, TCor halfWidth);

private:

  /** Convenience function to possibly plot a pixel with swapped x/y coordinates (needed for steep 
  lines). */
  inline void plot(int x, int y, TWgt weight, bool swapXY)
  {
    if(swapXY)
      ImageDrawer::plot(y, x, weight);
    else
      ImageDrawer::plot(x, y, weight);
  }

  /** Sets up the internal variables for the line drawing algorithm for the two given 
  endpoints. */
  void setupAlgorithmVariables(TCor x0, TCor y0, TCor x1, TCor y1);

  /** Draws the middle section between the end caps. */
  void drawMiddleSection();

  /** Draws the left end cap of the line. */
  inline void drawLeftCap() { drawCap(xs, xel); }

  /** Draws the right end cap of the line. */
  inline void drawRightCap() { drawCap(xsr, xe); }

  /** Draws either left opr right end cap, called internally from drawLeftCap and drawRightCap. The
  code is exactly the same except for different loop start and end indices, because we need to 
  check against both caps everytime - because any pixel may be part of both caps at the same time
  (occurs for very short slanted lines). */
  void drawCap(int start, int end);

  /** A special variant of the left cap drawing code to be used for joined lines (in lineTo). It 
  avoids drawing anything within the left end-circle. When used in round caps mode, it avoids the
  phantom circles that woul otherwise appear in the line joints. */
  //void drawLeftCapForJoint();
   // currently works only in round cap mode

  void drawCapForJoint(int start, int end, TCor x1, TCor y1);



  // internal variables for the actual drawing algorithm:
  TCor dx, dy, a, b, yf, dp, d, L, A, B, C0, C1, AxBy;
  TWgt sc; // scaler for color
  int xMax, yMax, xs, xe, xel, xsr, ys, ye, x, y, dvy;
  bool steep;

};

#endif