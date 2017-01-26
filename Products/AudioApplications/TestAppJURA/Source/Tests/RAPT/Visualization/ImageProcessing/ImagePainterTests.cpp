#include "ImagePainterTests.h"

PaintCanvas::PaintCanvas()
{
  brightness = 1.f;

  paintImage.setMaxSize(2000, 1000);
  dotMask.setMaxSize(50, 50);

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
  g.drawImage(displayImage, Rectangle<float>(0.f, 0.f, (float) getWidth(), (float) getHeight()));
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

void PaintCanvas::setDotBrightness(double newBrightness)
{
  brightness = (float) newBrightness;
}

// misc:

void PaintCanvas::paintDot(int x, int y)
{
  painter.paintDotViaMask(x, y, brightness);
  repaint();
}

//=================================================================================================

PainterComponent::PainterComponent()
{
  addAndMakeVisible(canvas);


  sliderSize.setSliderName("Size");
  sliderSize.addListener(this);
  sliderSize.setRange(0.0, 50.0, 0.125, 15.0);
  addAndMakeVisible(sliderSize);

  sliderBlur.setSliderName("Blur");
  sliderBlur.addListener(this);
  addAndMakeVisible(sliderBlur);

  sliderBrightness.setSliderName("Brightness");
  sliderBrightness.addListener(this);
  addAndMakeVisible(sliderBrightness);

  updatePreviewDot();
}

void PainterComponent::paint(Graphics &g)
{
  g.fillAll(Colour::greyLevel(0.25f));

  float x, y, w, h;
  x = sliderBrightness.getX();
  y = sliderBrightness.getBottom() + 8;
  w = previewDot.getWidth();
  h = previewDot.getHeight();
  g.drawImage(previewDot, Rectangle<float>(x, y, w, h));
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

  sliderSize.setBounds(x, y, w, h); y += dy;
  sliderBlur.setBounds(x, y, w, h); y += dy;
  sliderBrightness.setBounds(x, y, w, h); y += dy;
}

void PainterComponent::rSliderValueChanged(RSlider* rSlider)
{
  if(rSlider == &sliderSize)
    canvas.setDotSize(sliderSize.getValue());
  else if(rSlider == &sliderBrightness)
    canvas.setDotBrightness(sliderBrightness.getValue());
  else if(rSlider == &sliderBlur)
    canvas.setDotBlur(sliderBlur.getValue());

  updatePreviewDot();

  // \todo: maybe we should use a baseclass to maintain a set of parameters instead of deriving 
  // from RSliderListener - but this would require some refactoring
}

void PainterComponent::updatePreviewDot()
{
  int w = canvas.dotMask.getWidth();
  int h = canvas.dotMask.getHeight();
  previewDot = juce::Image(juce::Image::ARGB, w, h, false);
  normalizedDataToImage(canvas.dotMask.getPixelPointer(0, 0), previewDot);
  repaint();
}