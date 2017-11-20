#ifndef rojue_GraphicsTools_h
#define rojue_GraphicsTools_h

#include "../includesForRojue.h"
//#include "rojue_RLookAndFeel.h"
#include "rojue_BitmapFont.h"

namespace rojue
{

  /** Copies the pixel-data from one image into another. For this to work, the images must be compatible (same dimensions and 
  pixel-format). */
  void copyImage(juce::Image *sourceImage, juce::Image *targetImage);

  /** Draws a text with a BitmapFont with the given style-settings and returns the x-coordinate where the drawn text ends (and subsequent 
  text can be appended, if desired). */
  int drawBitmapFontText(Graphics &g, int x, int y, const juce::String& textToDraw, const BitmapFont* fontToUse, 
    const Colour& colourToUse, int kerning = -1, Justification justification = Justification::topLeft);

  /** Returns a colour that is a weighted mix of two other colours. It will simply linearly 
  interpolates the ARGB values separately. */
  Colour getMixedColour(const Colour colour1, const Colour colour2, 
    double weight1 = 0.5, double weight2 = 0.5);

  /** Fills an rectangle with a 4-point bilinear gradient. */
  void fillRectWithBilinearGradient(
    Graphics &g, 
    int x, int y, int w, int h, 
    Colour topLeftColour     = Colours::red, 
    Colour topRightColour    = Colours::green, 
    Colour bottomLeftColour  = Colours::blue, 
    Colour bottomRightColour = Colours::white);

  /** Fills an rectangle with a 4-point bilinear gradient. */
  void fillRectWithBilinearGradientSlow(
    Graphics &g, 
    int x, int y, int w, int h, 
    Colour topLeftColour     = Colours::red, 
    Colour topRightColour    = Colours::green, 
    Colour bottomLeftColour  = Colours::blue, 
    Colour bottomRightColour = Colours::white);

  /** Fills an rectangle with a 4-point bilinear gradient. */
  void fillRectWithBilinearGradient(
    Graphics &g, Rectangle<int> r, 
    Colour topLeftColour     = Colours::red, 
    Colour topRightColour    = Colours::green, 
    Colour bottomLeftColour  = Colours::blue, 
    Colour bottomRightColour = Colours::white);

  /** Fills an rectangle with a 4-point bilinear gradient. */
  void fillRectWithBilinearGrayScaleGradient(
    Graphics &g, 
    int x, int y, int w, int h, 
    float topLeftWhite     = 0.9, 
    float topRightWhite    = 0.7, 
    float bottomLeftWhite  = 0.7, 
    float bottomRightWhite = 0.9);

  /** Fills an rectangle with the default background. */
  void fillRectWithDefaultBackground(Graphics &g, int x, int y, int w, int h);

  /** Fills an rectangle with the default background. */
  void fillRectWithDefaultBackground(Graphics &g, Rectangle<int> r);

  /** Draws a triangle based on 3 points and optionally fills it. */
  void drawTriangle(Graphics &g, float x1, float y1, float x2, float y2, 
    float x3, float y3, bool fill);

  void drawBlockDiagramPlus(Graphics &g, float x, float y, float w, float h, float thickness);


  /** Prolongs and/or shortens the line through x1,y1 and x2,y2 such that it fits into the rectangle
  bounded by xMin, yMin, and xMax yMax. */
  void fitLineToRectangle(double &x1, double &y1, double &x2, double &y2, 
    double xMin, double yMin, double xMax, double yMax);

  /** Clips a line off at boundaries xMin, yMin, and xMax yMax. */
  void clipLineToRectangle(double &x1, double &y1, double &x2, double &y2, 
    double xMin, double yMin, double xMax, double yMax);



}

#endif