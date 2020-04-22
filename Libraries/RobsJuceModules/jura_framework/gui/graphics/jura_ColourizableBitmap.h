#ifndef jura_ColourizableBitmap_h
#define jura_ColourizableBitmap_h

/** This class ...needs to be documented

\todo try to get rid of inheriting from juce::Image - newer versions of juce disallow this by 
declaring the image class final. I had to hack the juce library and remove this final declaration 
in my copy of the juce source tree - which is a dirty and hacky thing to do. */

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
