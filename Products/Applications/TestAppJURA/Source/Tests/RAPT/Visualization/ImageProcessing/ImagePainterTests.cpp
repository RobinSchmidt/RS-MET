#include "ImagePainterTests.h"

PaintCanvas::PaintCanvas()
{
  brightness = 1.f;

  paintImage.setMaxSize(2000, 1000);
  dotMask.setMaxSize(200, 200);

  painter.setImageToPaintOn(&paintImage);
  painter.setAlphaMaskForDot(&dotMask);
}

// Component callbacks:

void PaintCanvas::mouseDown(const MouseEvent &e)
{
  paintDot(e.x, e.y);
}

void PaintCanvas::mouseDrag(const MouseEvent &e)
{
  paintDot(e.x, e.y);
}

void PaintCanvas::paint(Graphics &g)
{
  normalizedDataToImage(paintImage.getPixelPointer(0, 0), displayImage);  // maybe use colormap
  g.drawImage(displayImage, juce::Rectangle<float>(0.f, 0.f, (float) getWidth(), (float) getHeight()));
}

void PaintCanvas::resized()
{
  int w, h;
  w = getWidth();
  h = getHeight();
  paintImage.setSize(w, h);
  displayImage = juce::Image(juce::Image::ARGB, w, h, false);
}

// setup:

void PaintCanvas::setDotSize(double newSize)
{
  dotMask.setSize(newSize);
}

void PaintCanvas::setDotBlur(double newBlur)
{
  dotMask.setTransitionWidth(newBlur);
}

void PaintCanvas::setDotInnerSlope(double newSlope)
{
  dotMask.setInnerSlope(newSlope);
}

void PaintCanvas::setDotOuterSlope(double newSlope)
{
  dotMask.setOuterSlope(newSlope);
}

void PaintCanvas::setDotBrightness(double newBrightness)
{
  brightness = (float) newBrightness;
}

// misc:

void PaintCanvas::clearImage()
{
  paintImage.clear();
  repaint();
}

void PaintCanvas::paintDot(int x, int y)
{
  painter.paintDotViaMask(x, y, brightness);

  //painter.paintDotViaMask(x+0.25f, y+0.75f, brightness); // test - use anti-aliased function
  //painter.paintDot3x3(x, y, brightness, 0.5, 0.25); // test
  repaint();
}

//=================================================================================================

PainterComponent::PainterComponent()
{
  addAndMakeVisible(canvas);

  sliderSize.setSliderName("Size");
  sliderSize.addListener(this);
  //sliderSize.setRange(1.0, 160.0, 0.125, 50.0);
  sliderSize.setRange(1.0, 160.0, 1.0, 50.0);
  addAndMakeVisible(sliderSize);
  rSliderValueChanged(&sliderSize);

  sliderBlur.setSliderName("Blur");
  sliderBlur.addListener(this);
  addAndMakeVisible(sliderBlur);
  rSliderValueChanged(&sliderBlur);

  sliderInnerSlope.setSliderName("InnerSlope");
  sliderInnerSlope.addListener(this);
  sliderInnerSlope.setRange(0.0, 3.0, 0.0, 0.5);
  addAndMakeVisible(sliderInnerSlope);
  rSliderValueChanged(&sliderInnerSlope);

  sliderOuterSlope.setSliderName("OuterSlope");
  sliderOuterSlope.addListener(this);
  sliderOuterSlope.setRange(0.0, 3.0, 0.0, 1.0);
  addAndMakeVisible(sliderOuterSlope);
  rSliderValueChanged(&sliderOuterSlope);

  sliderBrightness.setSliderName("Brightness");
  sliderBrightness.addListener(this);
  addAndMakeVisible(sliderBrightness);
  rSliderValueChanged(&sliderBrightness);

  clearButton.setButtonText("Clear");
  clearButton.addRButtonListener(this);
  addAndMakeVisible(clearButton);


  //updatePreviewDot();
}

void PainterComponent::paint(Graphics &g)
{
  g.fillAll(Colour::greyLevel(0.25f));

  float x, y, w, h;
  x = (float)sliderBrightness.getX();
  y = (float)sliderBrightness.getBottom() + 8.f;
  w = (float)previewDot.getWidth();
  h = (float)previewDot.getHeight();
  g.drawImage(previewDot, juce::Rectangle<float>(x, y, w, h));
}

void PainterComponent::resized()
{
  int x, y, w, h, dy;

  int widgetZoneWidth = 180;
  int widgetHeight = 16;
  int margin = 4;

  // set up canvas:
  x = 0;
  y = 0;
  w = getWidth()  - x - widgetZoneWidth;
  h = getHeight() - y;

  canvas.setBounds(x, y, w, h);

  // arrange widgets:
  x  = w + margin;
  y  = margin;
  w  = widgetZoneWidth - 2*margin;
  h  = widgetHeight;
  dy = h+margin;

  clearButton.setBounds(     x, y, w, h); y += dy;
  sliderSize.setBounds(      x, y, w, h); y += dy;
  sliderBlur.setBounds(      x, y, w, h); y += dy;
  sliderInnerSlope.setBounds(x, y, w, h); y += dy;
  sliderOuterSlope.setBounds(x, y, w, h); y += dy;
  sliderBrightness.setBounds(x, y, w, h); y += dy;
}

void PainterComponent::rSliderValueChanged(RSlider* rSlider)
{
  if(rSlider == &sliderSize)
    canvas.setDotSize(sliderSize.getValue());
  else if(rSlider == &sliderBlur)
    canvas.setDotBlur(sliderBlur.getValue());
  else if(rSlider == &sliderInnerSlope)
    canvas.setDotInnerSlope(sliderInnerSlope.getValue());
  else if(rSlider == &sliderOuterSlope)
    canvas.setDotOuterSlope(sliderOuterSlope.getValue());
  else if(rSlider == &sliderBrightness)
    canvas.setDotBrightness(sliderBrightness.getValue());

  updatePreviewDot();

  // \todo: maybe we should use a baseclass to maintain a set of parameters instead of deriving 
  // from RSliderListener - but this would require some refactoring
}

void PainterComponent::rButtonClicked(RButton* button)
{
  if(button == &clearButton)
    canvas.clearImage();
}

void PainterComponent::updatePreviewDot()
{
  int w = canvas.dotMask.getWidth();
  int h = canvas.dotMask.getHeight();
  previewDot = juce::Image(juce::Image::ARGB, w, h, false);
  normalizedDataToImage(canvas.dotMask.getPixelPointer(0, 0), previewDot);
  repaint();
}

//=================================================================================================

//bool unitTestImagePainter()
//{
//
//  return true;
//}

//
//PainterUnitTest::PainterUnitTest() : UnitTest("ImagePainter") 
//{
//
//}
//
//void PainterUnitTest::runTest()
//{
//  beginTest("Alpha Mask Dot Painting");
//  //expect(maskDot1x1());
//  //expect(maskDot2x2());
//  //expect(maskDot3x3());
//  //expect(maskDot4x4());
//  expect(maskDot5x5());
//}
//
//bool PainterUnitTest::maskDot1x1()
//{
//  return true; // preliminary
//}
//bool PainterUnitTest::maskDot5x5()
//{
//  return true; // preliminary
//}

