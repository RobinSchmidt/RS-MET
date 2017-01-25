#ifndef IMAGEPAINTERTESTS_H_INCLUDED
#define IMAGEPAINTERTESTS_H_INCLUDED

#include "../../../../../JuceLibraryCode/JuceHeader.h"
using namespace jura;

/** A component to be used to paint on a RAPT::Image object and to display it on the screen. Mainly 
intended to test the RAPT::ImagePainter class */

class PaintCanvas : public Component
{

public:

  PaintCanvas();


  //virtual void mouse


protected:

  RAPT::Image<float> paintImage, dotMask;
  //RAPT::ImagePainter<float, float, float> painter;

  juce::Image displayImage;  // the juce image that is used for displaying the paintImage


};

//=================================================================================================
 
/** */

class PainterComponent : public Component
{

public:

  PainterComponent();

protected:

  PaintCanvas canvas;

  // widgets for setting up the brush/pen:
  RSlider sliderSize, sliderBlur, sliderBrightness;

};

#endif