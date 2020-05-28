#ifndef jura_ColourizableBitmap_h
#define jura_ColourizableBitmap_h

/** This class ...needs to be documented...it's actually used as an alpha-mask to composit the 
shape stored in the alpha-mask with some existing juce::Image in a particular color - i.e. the 
desired shape (stored in the alpha-mask) is composited in the desired color with the existing 
image. This is used as basis for the pixel-fonts.

\todo try to get rid of inheriting from juce::Image - newer versions of juce disallow this by 
declaring the Image class final. I had to hack the juce library and remove this final declaration 
in my copy of the juce source tree - which is a dirty and hacky thing to do. ..but how? Maybe use 
rsMatrix or rsImage as baseclass in place of juce::Image and in renderImage, pass the juce::Image 
onto which we should render? Maybe during transition, provide two renderImage functions - one 
without parameters (the current one) and one with a reference to a juce::Image representing the 
image to render on. Then transition to using the second one everywhere, then remove the first, 
then replace the juce::Image baseclass with rsImage or something. */

class JUCE_API ColourizableBitmap : public juce::Image
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ColourizableBitmap(const int width, const int height, unsigned char *newOpacities);

  /** Destructor. */
  virtual ~ColourizableBitmap();

  //-----------------------------------------------------------------------------------------------
  // setup:

  virtual void setColour(const Colour& newColour);

  /** You must take care that the passed array has the correct size. */
  virtual void setPixelOpacities(unsigned char* newOpacities);


  virtual void setPixelOpacityAt(int x, int y, unsigned char newOpacity)
  {
    pixelOpacities[y*Image::getWidth()+x] = newOpacity;
  }

  virtual void clearOpacities();

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  virtual unsigned char getPixelOpacityAt(int x, int y)
  {
    return pixelOpacities[y*Image::getWidth()+x];
  }

  virtual float getFloatPixelOpacityAt(int x, int y)
  {
    return (1.f/255.f) * pixelOpacities[y*Image::getWidth()+x];
  }

  juce_UseDebuggingNewOperator;

protected:

  /** Renders the image at some particular Colour. */
  virtual void renderImage();

  unsigned char*  pixelOpacities;
  Colour currentColour;

};

#endif  
