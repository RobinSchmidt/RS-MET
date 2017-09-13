#ifndef RAPT_IMAGEDRAWER_H_INCLUDED
#define RAPT_IMAGEDRAWER_H_INCLUDED

/** A baseclass for drawing on images. It consolidates the data and functionality that all drawers
have in common, regardless of what they draw. Subclasses like LineDrawer will do the actual drawing
and likely define more data and methods. */

template<class TPix, class TWgt, class TCor>  // pixel, weight, coordinate types
class rsImageDrawer
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
  rsImageDrawer(rsImage<TPix> *imageToDrawOn);


  /** \name Setup */

  /** Sets the image that we will draw on. */
  void setImageToDrawOn(rsImage<TPix> *imageToDrawOn);

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

  rsImage<TPix> *image;

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
class LineDrawer : public rsImageDrawer<TPix, TWgt, TCor>
{

public:


  /** The function to determine the brightness/weight as function of the perpendicular distance
  to the ideal geometric line. */
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
  LineDrawer(rsImage<TPix> *imageToDrawOn);


  /** \name Setup */

  /** Selects one of the line profiles. */
  void setLineProfile(int newProfile);

  /** Sets the width of the line in pixels. */
  void setLineWidth(TCor newWidth);

  /** Sets the drawer up to draw half-circular end caps, if true. If false, the caps will be
  rectangular. */
  void setRoundCaps(bool shouldBeRound);


  /** \name Drawing */

  /** Draws a line from (x0,y0) to (x1,y1). The joinableStart/End parameters scale the color weight
  inside the end-circles (in round cap mode) by 0.5 - this is supposed to lead to nice line
  joints, but that doesn't work yet - the joints still look funky :-(  */
  void drawLine(TCor x0, TCor y0, TCor x1, TCor y1, bool joinableStart = false,
    bool joinableEnd = false);

  /** Special line drawing function that is supposed to be used for drawing sequences of connected
  lines. After an initial call to initLine or a previous call to lineTo, you can call this
  function in order to avoid artifacts (phantom circles) at the line joints. The optional
  uniformColor parameter indicates that between successive calls to lineTo, the color will not be
  changed (by setColor). If the whole polyline has the same color, you can pass true and then a
  simpler, more efficient end cap handling code will be invoked. */
  void lineTo(TCor x, TCor y, bool uniformColor = false);

  /** Function to initialize our xOld, yOld members which are used for polyline drawing. Call this
  once with the start point of the polyline before repeatedly calling lineTo  */
  inline void initPolyLine(TCor x, TCor y)
  {
    xOld = x;
    yOld = y;
    // maybe it should draw a circle?
  }


  //void drawConnectedLine(TCor x0, TCor y0, TCor x1, TCor y1);



  // todo: have members for simplified 1-pixel line drawing: drawLineWu, drawLineBresenham


protected:

  bool roundCaps = true;
  TCor w2;                // lineWidth/2
  //TCor x0 = 0, y0 = 0;    // start-point for lineTo function
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
      rsImageDrawer<TPix, TWgt, TCor>::plot(y, x, weight);
    else
      rsImageDrawer<TPix, TWgt, TCor>::plot(x, y, weight);
  }

  /** Sets up the internal variables for the line drawing algorithm for the two given
  endpoints. */
  void setupAlgorithmVariables(TCor x0, TCor y0, TCor x1, TCor y1);

  /** Draws the middle section between the end caps. */
  void drawMiddleSection();

  /** Draws the left end cap of the line. */
  inline void drawLeftCap(bool joinable) { drawCap(xs, xel, joinable); }

  /** Draws the right end cap of the line. */
  inline void drawRightCap(bool joinable) { drawCap(xsr, xe, joinable); }

  /** Draws either left or right end cap, called internally from drawLeftCap and drawRightCap. The
  code is exactly the same except for different loop start and end indices, because we need to
  check against both caps everytime - because any pixel may be part of both caps at the same time
  (occurs for very short slanted lines). */
  void drawCap(int start, int end, bool joinable = false);

  /** A special variant of the cap drawing code to be used for joined lines (in lineTo) with
  uniform color. It avoids drawing anything within the left end-circle. When used in round caps
  mode, it avoids the phantom circles that woul otherwise appear in the line joints. xj, yj are the
  coordinates of the cap (which sits inside the joint) that shall not be drawn. */
  void drawCapForJointUniformColor(int start, int end, TCor xj, TCor yj);



  // internal variables for the actual drawing algorithm:
  TCor dx, dy, a, b, yf, dp, d, L, A, B, C0, C1, AxBy;
  TWgt sc; // scaler for color
  int xMax, yMax, xs, xe, xel, xsr, ys, ye, x, y, dvy;
  bool steep, back;


  TCor xOld = 0, yOld = 0;

};

#endif
