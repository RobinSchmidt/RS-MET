#include "rojue_ColourizableBitmap.h"
using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ColourizableBitmap::ColourizableBitmap(const int width, const int height, unsigned char *newOpacities) 
: Image(Image::ARGB, width, height, true)
{
  pixelOpacities = new unsigned char[width*height];
  clearOpacities();
  currentColour = Colours::white;
  setPixelOpacities(newOpacities);
}

ColourizableBitmap::~ColourizableBitmap()
{
  deleteAndZero(pixelOpacities);
}

//-------------------------------------------------------------------------------------------------
// setup:

void ColourizableBitmap::setColour(const juce::Colour &newColour)
{
  if( newColour == currentColour )
    return; // nothing to do
  else
  {
    currentColour = newColour;
    renderImage();
  }
}

void ColourizableBitmap::setPixelOpacities(unsigned char *newOpacities)
{
  int w = Image::getWidth();
  int h = Image::getHeight();
  //float opacity;
  for(int y=0; y<h; y++)
  {
    for(int x=0; x<w; x++)
    {
      //opacity = newOpacities[y*w+x];
      pixelOpacities[y*w+x] = newOpacities[y*w+x];
    }
  }
  renderImage();
}

void ColourizableBitmap::clearOpacities()
{
  for(int x=0; x<Image::getWidth(); x++)
  {
    for(int y=0; y<Image::getHeight(); y++)
    {
      setPixelOpacityAt(x, y, 0);
    }
  }
  Image::clear(Rectangle<int>(0, 0, Image::getWidth(), Image::getHeight()), Colours::transparentBlack);
}

void ColourizableBitmap::renderImage()
{
  // initialize the image with all pixels having the specified colour:
  Image::clear(Rectangle<int>(0, 0, Image::getWidth(), Image::getHeight()), currentColour);

  // multiply all the alphas with the value stored in our array of the opacities:
  //float opacity;
  for(int x=0; x<Image::getWidth(); x++)
  {
    for(int y=0; y<Image::getHeight(); y++)
    {
      //opacity = getPixelOpacityAt(x, y);
      Image::multiplyAlphaAt(x, y, getFloatPixelOpacityAt(x, y));
    }
  }
}

//-------------------------------------------------------------------------------------------------
// inquiry:

