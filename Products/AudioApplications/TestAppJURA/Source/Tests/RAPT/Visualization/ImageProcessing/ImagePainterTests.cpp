#include "ImagePainterTests.h"

PaintCanvas::PaintCanvas()
{
  brightness = 1.f;

  paintImage.setMaxSize(2000, 1000);
  dotMask.setMaxSize(50, 50);

  painter.setImageToPaintOn(&paintImage);
  painter.setAlphaMaskForDot(&dotMask);
}

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
  normalizedDataToImage(paintImage.getPixelPointer(0, 0), displayImage);
  //normalizedDataToImage(paintImage.getPixelPointer(0, 0), displayImage, colorMap); // maybe later
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
  sliderSize.setRange(0.0, 50.0, 0.125, 15.0);
  addAndMakeVisible(sliderSize);

  sliderBlur.setSliderName("Blur");
  addAndMakeVisible(sliderBlur);

  sliderBrightness.setSliderName("Brightness");
  addAndMakeVisible(sliderBrightness);
}

void PainterComponent::paint(Graphics &g)
{
  g.fillAll(Colour::greyLevel(0.25f));
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
  // \todo: update the respective setting in the ImagePainter
}