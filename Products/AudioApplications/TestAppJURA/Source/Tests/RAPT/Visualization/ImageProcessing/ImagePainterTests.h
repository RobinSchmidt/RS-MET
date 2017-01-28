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
etc. */

class PainterComponent : public Component, public RSliderListener
{

public:

  PainterComponent();
  
  virtual void paint(Graphics &g)	override;
  virtual void resized() override;

  // widget callbacks:
  virtual void rSliderValueChanged(RSlider* rSlider) override;

protected:

  /** Updates the image for previewing the dot */
  void updatePreviewDot();

  PaintCanvas canvas;
  juce::Image previewDot; 

  // widgets for setting up the brush/pen:
  RSlider sliderSize, sliderBlur, sliderBrightness;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PainterComponent)
};

/** Runs a unit test for the RAPT::ImagePainter class and returns true when the test was passed and 
false if the test has failed. */
bool unitTestImagePainter();


////=================================================================================================
//
///** A unit test for the RAPT::ImagePainter class. */
//
//class PainterUnitTest : public UnitTest
//{
//
//public:
//
//  PainterUnitTest();
//  virtual void runTest() override;
//
//protected:
//
//  bool maskDot1x1();
//  bool maskDot5x5();
//
//  RAPT::ImageResizable<float> image;               // image to paint on
//  RAPT::AlphaMask<float> mask;                     // prototype dot
//  RAPT::ImagePainter<float, float, float> painter; // painter object 
//
//  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PainterUnitTest)
//};


#endif