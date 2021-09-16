#ifndef jura_GraphicsTools_h
#define jura_GraphicsTools_h

// Maybe wrap these functions into a class rsDrawer or rsImageTools at some point 

/** Copies the pixel-data from one image into another. For this to work, the images must be 
compatible (same dimensions and pixel-format). */
JUCE_API void copyImage(juce::Image *sourceImage, juce::Image *targetImage);

/** Given an image (assumed to be of RGBA type), this function returns the indices of the red, 
green, blue and alpha components, i.e. assigns ri,gi,bi,ai to some permutation of 0,1,2,3. 
\TODO: the implementation is preliminary and just returns a fixed order that seems to be valid on
a windows PC. this is likely to be changed for other platforms. */
JUCE_API void colorComponentIndices(juce::Image& image, int &ri, int &gi, int &bi, int &ai);

/** Converts the given data (assumed to represent a flattened representation of 2D data with width 
and height equal to the dimensions of the passed image) into pixel colors on the image using the 
given colormap. The data is assumed to be normalized to the range 0..1 (ends inclusive). 
WARNING: it is VERY important that the data is in this range indeed, otherwise you'll get access 
violations!!! */
JUCE_API void normalizedDataToImage(float *data, juce::Image &image, 
  const jura::ColorMap& colorMap);
// \todo: make a version of this function that deals with unnormalized data - we need to pass a min
// and max value and then use these internally to scale the data values during readout

/** Similar to version of this function with the colormap but uses grayscale values. It also 
assumes normalized data, but if it's not you won't get access violations but just a wraparound from 
white to black when the data value is > 1. */
JUCE_API void normalizedDataToImage(float *data, juce::Image &image);

/** Converts a data-matrix of float values into an image. The matrix values are supposed to be in 
the range 0..1 and their values will determine the pixel brightnesses. You can set a base color in 
terms of red/green/blue values. The output image will then use the color so defined multiplied by 
the matrix values for the pixels. */
JUCE_API void dataMatrixToImage(float **data, juce::Image &image, 
  uint8 red = 255, uint8 green = 255, uint8 blue = 255);

/** Takes a pointer to floating point data and converts it into an image. The float array is 
supposed to be of length 4*w*h where w,h are width and height of the image. Presumably the array 
represents a w-times-h matrix of data where each datapoint consists of 4 float numbers 
representing red, green, blue and alpha (although the alpha channel is not used here - we write the
RGB values with full opacity into the image). So it's a conversion function from an internal 
float-RGBA pixel format to the juce::Image pixel format. */
JUCE_API void dataToImageOpaqueFloat32x4(float *data, juce::Image &image);
// maybe rename into something more meaningful

/** Draws a text with a BitmapFont with the given style-settings and returns the x-coordinate where 
the drawn text ends (and subsequent text can be appended, if desired). */
JUCE_API int drawBitmapFontText(Graphics &g, int x, int y, const juce::String& textToDraw, 
  const BitmapFont* fontToUse, const Colour& colourToUse, int kerning = -1, 
  Justification justification = Justification::topLeft);

/** Returns a colour that is a weighted mix of two other colours. It will simply linearly
interpolates the ARGB values separately. */
JUCE_API Colour getMixedColour(const Colour colour1, const Colour colour2,
  double weight1 = 0.5, double weight2 = 0.5);

/** Fills an rectangle with a 4-point bilinear gradient. */
JUCE_API void fillRectWithBilinearGradient(
  Graphics &g,
  int x, int y, int w, int h,
  Colour topLeftColour     = Colours::red,
  Colour topRightColour    = Colours::green,
  Colour bottomLeftColour  = Colours::blue,
  Colour bottomRightColour = Colours::white);

/** Fills an rectangle with a 4-point bilinear gradient. */
JUCE_API void fillRectWithBilinearGradientSlow(
  Graphics &g,
  int x, int y, int w, int h,
  Colour topLeftColour     = Colours::red,
  Colour topRightColour    = Colours::green,
  Colour bottomLeftColour  = Colours::blue,
  Colour bottomRightColour = Colours::white);

/** Fills an rectangle with a 4-point bilinear gradient. */
JUCE_API void fillRectWithBilinearGradient(
  Graphics &g, juce::Rectangle<int> r,
  Colour topLeftColour     = Colours::red,
  Colour topRightColour    = Colours::green,
  Colour bottomLeftColour  = Colours::blue,
  Colour bottomRightColour = Colours::white);

/** Fills an rectangle with a 4-point bilinear gradient. */
JUCE_API void fillRectWithBilinearGrayScaleGradient(
  Graphics &g,
  int x, int y, int w, int h,
  float topLeftWhite     = 0.9,
  float topRightWhite    = 0.7,
  float bottomLeftWhite  = 0.7,
  float bottomRightWhite = 0.9);

/** Fills an rectangle with the default background. */
JUCE_API void fillRectWithDefaultBackground(Graphics &g, int x, int y, int w, int h);

/** Fills an rectangle with the default background. */
JUCE_API void fillRectWithDefaultBackground(Graphics &g, juce::Rectangle<int> r);

/** Draws a triangle based on 3 points and optionally fills it. */
JUCE_API void drawTriangle(Graphics &g, float x1, float y1, float x2, float y2,
  float x3, float y3, bool fill);

JUCE_API void drawBlockDiagramPlus(Graphics &g, float x, float y, float w, float h, float thickness);

/** Prolongs and/or shortens the line through x1,y1 and x2,y2 such that it fits into the rectangle
bounded by xMin, yMin, and xMax yMax. */
JUCE_API void fitLineToRectangle(double &x1, double &y1, double &x2, double &y2,
  double xMin, double yMin, double xMax, double yMax);

/** Clips a line off at boundaries xMin, yMin, and xMax yMax. */
JUCE_API void clipLineToRectangle(double &x1, double &y1, double &x2, double &y2,
  double xMin, double yMin, double xMax, double yMax);

//=================================================================================================
// functions for coordinate system drawing (they replace the old, clunky and inflexible member 
// functions of rsPlot) ...maybe put into a class CoordinateSystemDrawer - more
// convenient (the functions tend to have too many parameters):

// maybe move into rsPlotDrawer (maybe as as static functions)

/** Sets up the output range of the passed coordinate mapper to the bounds of the given 
component. */
JUCE_API void setupCoordinateMapper(RAPT::rsCoordinateMapper2D<double>& mapper, 
  const Component* cmp);

/** Sets up the output range of the passed coordinate mapper to the bounds of the given 
image. */
JUCE_API void setupCoordinateMapper(RAPT::rsCoordinateMapper2D<double>& mapper, 
  const Image* img);

/** Sets up the output range of the passed coordinate mapper to the bounds of the given 
svg xml element. The svg must have "width" and "height" attributes. */
JUCE_API void setupCoordinateMapper(RAPT::rsCoordinateMapper2D<double>& mapper, 
  const XmlElement* svg);

/** Specifically used in coordinate-system drawing - maybe make more general, rename. */
JUCE_API void drawBitmapText(Graphics &g, const juce::String &text, double x, double y,
  double w, double h, BitmapFont const* font, Justification justification, Colour color);

/** Draws equidistant horizontal grid-lines. The spacing must be given in model coordinates and 
will be interpreted as a factor in the case of a logarithmic y-axis. */
JUCE_API void drawHorizontalGrid(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, float thickness);

/** Draws tick marks for the y-axis with values next to them. Tick marks are drawn in whatever 
color the Graphics object is set up to, values are drawn in textColor. */
JUCE_API void drawAxisValuesY(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, double xPosition, juce::String (*yToString) (double y), Colour textColor);

/** Analog to drawHorizontalGrid. */
JUCE_API void drawVerticalGrid(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, float thickness);

/** Analog to drawAxisValuesY. */
JUCE_API void drawAxisValuesX(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, double yPosition, juce::String (*xToString) (double x), Colour textColor);

/** Draws concentric circles centered at the origin of given spacing of the radii. */
JUCE_API void drawRadialGrid(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, float thickness);

/** Draws straight lines that pass through the origin and extend over the visible range. The 
spacing of the angles should be given in degrees. */
JUCE_API void drawAngularGrid(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, float thickness);
 
/** Draws the x-axis for a coordinate system. */
JUCE_API void drawAxisX(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double yPosition, const juce::String& label, Colour labelColor);

/** Draws the y-axis for a coordinate system. */
JUCE_API void drawAxisY(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double xPosition, const juce::String& label, Colour labelColor);

//=================================================================================================
// corresponding functions that create the svg elements instead of drawing to a Graphics object

JUCE_API void addLineToSvgDrawing(XmlElement* svg, float x1, float y1, float x2, float y2,
  float thickness, Colour color, bool withArrowHead = false);

JUCE_API void addTextToSvgDrawing(XmlElement* svg, juce::String text, float x, float y,
  Justification justification, Colour color);

JUCE_API void drawHorizontalGrid(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, float thickness, Colour colour);

JUCE_API void drawAxisValuesY(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, double xPosition, juce::String (*yToString) (double y), Colour textColor);

JUCE_API void drawVerticalGrid(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, float thickness, Colour colour);

JUCE_API void drawAxisValuesX(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, double yPosition, juce::String (*xToString) (double x), Colour textColor);

JUCE_API void drawRadialGrid(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, float thickness, Colour color);

JUCE_API void drawAngularGrid(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double spacing, float thickness, Colour color);

JUCE_API void drawAxisX(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double yPosition, const juce::String& label, Colour color);

JUCE_API void drawAxisY(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper, 
  double xPosition, const juce::String& label, Colour color);


inline void setPixelRGB(uint8* p, uint8 r, uint8 g, uint8 b)
{
#if JUCE_MAC
  p[0] = r;
  p[1] = g;
  p[2] = b;
#else
  p[0] = b;       // blue comes first
  p[1] = g;       // green comes second
  p[2] = r;       // red comes third in memory
#endif

/*
//old:
#ifndef _MSC_VER   // #ifdef JUCE_LITTLE_ENDIAN seems to not solve the wrong-color thing
  p[0] = r;
  p[1] = g;
  p[2] = b;
#else
  p[0] = b;       // blue comes first,
  p[1] = g;       // green comes second,
  p[2] = r;       // red comes third in memory
#endif
  // This does not yet seem to work correctly on my new M1 mac
  // see juce::PixelRGB, PixelARGB, Colour
*/
}



#endif