#ifndef RS_GRAPHICSRENDERER2D_H
#define RS_GRAPHICSRENDERER2D_H

namespace RSLib
{


  /**

  A plain old data class for holding the state of an rsGraphicsRenderer2D object.

  */

  struct RSLib_API rsGraphicsRenderer2DState
  {
    rsColorRGBA color;      // maybe have fillColor and outlineColor members
    double lineThickness;

    // joinStyle, capStyle, miterLimit, strokeStyle, fillStyle (solid, gradient, texture, image)
    // *font, textJustification

    rsGraphicsRenderer2DState()
    {
      lineThickness = 1.f;
    }
  };

  //===============================================================================================

  /**

  This is an abstract class that represents an drawing canvas which is an area with associated
  functions that perform drawing operations.

  \todo: text-rendering (vector- and pixel-fonts (pixel first))

  \todo: stroke-styles, join-styles, cap-styles, bezier-splines, fill-styles
  (or maybe use "type" instead of "style")
  -but first: develop window management further
  -then: make a test-application: SimpleDraw -> test for widgets, drawing and XML (SVG)
   ->export to bmp, png, jpeg, ppm
   ->maybe also used a test for TiledWindow

   \todo: have a member rsStack<rsGraphicsRenderer2DState> to store a state-stack, have 
   member-functions pushState/popState


  */

  class RSLib_API rsGraphicsRenderer2D
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    rsGraphicsRenderer2D();

    /** Destructor. */
    virtual ~rsGraphicsRenderer2D();


    /** \name Setup */

    /** Sets the color to be used for subsequent drawing actions. */
    virtual void setColor(const rsColorRGBA& newColor);

    /** Sets the line-thickness to be used for subsequent drawing of outlined shapes. */
    virtual void setLineThickness(const double newThickness);


    /** \name Drawing */

    /** Override this function to fill the whole area with the color that is chosen in the current
    state.color variable. */
    virtual void fillAll() = 0;

    /** Override this function in order draw a filled polygon with the current color. */
    virtual void drawFilledPolygon(const rsPolygon2D<double> &polygon) = 0;

    /** Override this function in order draw an outlined filled polygon with the current color and 
    current line thickness. */
    virtual void drawOutlinedPolygon(const rsPolygon2D<double> &polygon) = 0;

    //virtual void drawLineStrip(int numVertices, double xValues, double yValues);

    /** ...... */
    virtual void drawText(const rsString &text, int x, int y, int width, int height) = 0;
       // more parameters rsJustification::centered, &rsNormalFont (or maybe put them into the 
       // state)

    /** Override this function in order to create a color-gradient between the 4 corners of an 
    rectangle. */
    virtual void fillRectangleWithBilinearGradient(int x, int y, int w, int h,
                                                   const rsColorRGBA &topLeftColor,
                                                   const rsColorRGBA &topRightColor,
                                                   const rsColorRGBA &bottomLeftColor,
                                                   const rsColorRGBA &bottomRightColor) = 0;

    /** Draws a line between (x1,y1) and (x2,y2) using the current color and line-thickness. The
    baseclass implementation creates an rectangle that represents the line and then uses 
    fillPolygon to render it. You may override it in subclasses in order to use a better 
    approach. */
    virtual void drawLine(double x1, double y1, double x2, double y2);

    /** Draws an outlined rectangle with origin at (x,y) and given width and height. 
    Uses fillRectangle */
    virtual void drawOutlinedRectangle(double x, double y, double width, double height);

    /** Draws a filled rectangle with origin at (x,y) and given width and height. 
    Uses fillPolygon */
    virtual void drawFilledRectangle(double x, double y, double width, double height);

    /** This function draws polygon that inscribed in an ellipse centered at (centerX, centerY) 
    with given horizontal and vertical radii and the given number of vertices. In the limit when 
    numVertices goes to infinity, it becomes an ellipse with given horizontal and vertical 
    radii. The initialAngle parameter specifies the angle of the first vertex with respect to
    the horizontal. */
    virtual void drawFilledEllipseInscribedPolygon(int numVertices, double centerX, double centerY, 
      double horizontalRadius, double verticalRadius, double initialAngle = 0.0);

    /** Draws a filled ellipse with center (centerX,centerY) and given horizontal and vertical 
    radii. Uses fillEllipseInscribedPolygon. */
    virtual void drawFilledEllipse(double centerX, double centerY, 
                                   double horizontalRadius, double verticalRadius);

    /** Draws an outlined pixel rectangle with one pixel width. The coordinates passed to this 
    function are assumed to be pixel-coordinates which differ by 0.5 pixels from grid coordinates
    (using the convention that a pixel is centered on its coordinates. */
    virtual void drawOutlinedPixelRectangle(int x, int y, int w, int h);

    /** Fills the rectangle given in pixel-coordinates with the current color. */
    virtual void drawFilledPixelRectangle(int x, int y, int w, int h);

    /*
    \todo
    virtual void drawTriangle(rsPoint2D<double> a, rsPoint2D<double> b, rsPoint2D<double> c);
    virtual void drawPolygon(const rsPolygon2D<double> &polygon) = 0;
    virtual void drawEllipse(double center, double horizontalRadius, double verticalRadius);
    virtual void drawText(const rsString &text, double xLeft, yBaseline);
    virtual void drawImage(const rsImage& imageToDraw, double xLeft, double yTop); 
      // hmm maybe include clipping and/or scaling
    */

    inline void drawLine(const rsPoint2D<double> &p1, const rsPoint2D<double> &p2)
    {
      drawLine(p1.x, p1.y, p2.x, p2.y);
    }
    //inline void drawLine(rsLine2D<double> theLine, double thickness);

    inline void drawOutlinedRectangle(const rsRectangle2D<double> &r)
    {
      drawOutlinedRectangle(r.x, r.y, r.w, r.h);
    }

    inline void drawOutlinedRectangle(const rsRectangle2D<rsInt32> &r)
    {
      drawOutlinedRectangle((double) r.x, (double) r.y, (double) r.w, (double) r.h);
    }

    // \todo drawRectangle, drawEllipse

  protected:

    /** \name Data */

    rsGraphicsRenderer2DState state;

  };

  //===============================================================================================

  /**

  This is an concrete subclass that implements rsGraphicsRenderer2D using software-based rendering 
  on an rsImageRegionRGBA object.

  \todo implement drawing of outlined polygons - i think, AGGs agg_vcgen_stroke can be used for
  this ("vertex generator" to generate a sequence of output vertices from a sequence of input
  vertices?) - so quite possibly, we may invoke fillPolygon once again, but this time for a 
  generated polygon that traces the outlines

  \todo write own implementations of rasterization algorithms using classes: 
  class rsPixelSpan - members: rsUint32 xStart, length; rsUint8 *coverages;
  class rsScanLine - members: rsUint32 y; rsArray<PixelSpan> spans;

  */

  class RSLib_API rsGraphicsRenderer2DImage : public rsGraphicsRenderer2D
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. Takes a reference to an rsImage object on which all drawing operations will be 
    performed as parameter. */
    rsGraphicsRenderer2DImage(rsImageRegionRGBA &imageRegionToRenderOn);

    /** Destructor. */
    virtual ~rsGraphicsRenderer2DImage();


    /** \name Setup */

    /** Sets the image region that this renderer should render on. */
    void setImageRegionToRenderOn(const rsImageRegionRGBA &newRegion)
    {
      imageRegion = newRegion;
    }


    /** \name Inquiry */

    /** Returns the image region that this renderer renders on. */
    rsImageRegionRGBA getImageRegionToRenderOn() const
    {
      return imageRegion;
    }


    /** \name Drawing Functions */

    // moveTo, lineTo

    // mandatory overrides:
    virtual void fillAll();
    virtual void drawFilledPolygon(const rsPolygon2D<double> &polygon);
    virtual void drawOutlinedPolygon(const rsPolygon2D<double> &polygon);
    virtual void drawText(const rsString &text, int x, int y, int width, int height);
    virtual void fillRectangleWithBilinearGradient(int x, int y, int w, int h,
                                                   const rsColorRGBA &topLeftColor,
                                                   const rsColorRGBA &topRightColor,
                                                   const rsColorRGBA &bottomLeftColor,
                                                   const rsColorRGBA &bottomRightColor);

    // optional overrides:
    virtual void drawOutlinedPixelRectangle(int x, int y, int w, int h);
    virtual void drawFilledPixelRectangle(int x, int y, int w, int h);

    /** Draws a horizontal pixel line at pixel-height y between x1 and x2. */
    virtual void drawHorizontalPixelLine(int y, int x1, int x2);

    /** Draws a vertical pixel line at pixel x between y1 and y2. */
    virtual void drawVerticalPixelLine(int x, int y1, int y2);


    // new member functions (not present in baseclass):

    /** Draws an aliased line using the Bresenham algorithm. */
    void drawBresenhamLine(int x1, int y1, int x2, int y2);

    /** Draws a vertical line with color gradient between yTop and yBottom at x with the given top 
    and bottom colors. */
    void drawVerticalGradientLine(int x, int yTop, int yBottom,  
                                  const rsColorRGBA &topColor, const rsColorRGBA &bottomColor);

    /** Draws a horizontal line with color gradient between xLeft and xRight at y with the given 
    left and right colors. */
    void drawHorizontalGradientLine(int y, int xLeft, int xRight,  
                                    const rsColorRGBA &leftColor, const rsColorRGBA &rightColor);

    /** Fills the with origin at (x,y) and width/height (w,h) with the currently selected 
    fill-color using the passed grayscale image as alpha-mask. The area will get a color that is
    somewhere between the values that are currently there and the new values that would be 
    determined by the filling color. The exact pixel colors will be determined by the values
    in the passed alpha(opacity)-mask. The optional xOffset and yOffset parameters define 
    coordinate offsets into the alpha-mask. These can be used to clip objects from left and top
    edges (for example drawText uses these to clip the text when it would otherwise extend beyond
    the desired box to the left and/or above). */
    void fillRectangleUsingAlphaMask(int x, int y, int w, int h, const rsImageGray *mask, 
      int xOffset = 0, int yOffset = 0);

    // image processing functions:
    // blur/contrast/brightness/greyOut/etc.

  protected:

        
    /** \name Misc */

    /** Clips the values x, y, w, h into the valid ranges. \todo rename to clipToValidRange */
    void clipToValidRange(int &x, int &y, int &w, int &h);


    /** \name Data */

    rsImageRegionRGBA &imageRegion;

  };

  //===============================================================================================

  /**

  This is an concrete subclass that implements rsGraphicsRenderer2D using OpenGL commands.

  */

  class RSLib_API rsGraphicsRenderer2DOpenGL : public rsGraphicsRenderer2D
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    rsGraphicsRenderer2DOpenGL();

    /** Destructor. */
    virtual ~rsGraphicsRenderer2DOpenGL();


    /** \name Setup */

    virtual void setColor(const rsColorRGBA& newColor);


    /** \name Drawing */

    //virtual void drawLine(rsPoint2D<double> p1, rsPoint2D<double> p2);

  };

}

#endif
