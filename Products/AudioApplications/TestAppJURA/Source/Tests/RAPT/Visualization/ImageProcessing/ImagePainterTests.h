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

  // Component overrides:
  virtual void mouseDown(const MouseEvent &event) override;
  virtual void mouseDrag(const MouseEvent &event) override;
  virtual void paint(Graphics &g)	override;
  virtual void resized() override;

  // dot setup:
  void setDotSize(double newSize);
  void setDotBlur(double newBlur);
  void setDotBrightness(double newBrightness);


protected:

  /** Paints a dot at the given position. */
  void paintDot(int x, int y);

  // image/painting stuff:
  RAPT::ImageResizable<float> paintImage;           // image to paint on
  RAPT::AlphaMask<float> dotMask;                   // prototype dot
  juce::Image displayImage;                         // for displaying the paintImage
  RAPT::ImagePainter<float, float, float> painter;  // handles to actual painting

  // data:
  float brightness;  // brightness of inserted dots


  friend class PainterComponent;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PaintCanvas)
};

//=================================================================================================
 
/** A component to wrap a PaintCanvas and add widgets for setting up the brush size, blur, color, 
etc. 

\todo:
-have a preview of the dot below the widgets

*/

class PainterComponent : public Component, public RSliderListener
{

public:

  PainterComponent();
  
  virtual void paint(Graphics &g)	override;
  virtual void resized() override;

  // widget callbacks:
  virtual void rSliderValueChanged(RSlider* rSlider) override;

protected:

  PaintCanvas canvas;

  // widgets for setting up the brush/pen:
  RSlider sliderSize, sliderBlur, sliderBrightness;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PainterComponent)
};

#endif