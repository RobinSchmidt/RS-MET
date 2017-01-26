#include "ImagePainterTests.h"

PaintCanvas::PaintCanvas()
{

}

void PaintCanvas::mouseDown(const MouseEvent &event)
{

}

void PaintCanvas::mouseDrag(const MouseEvent &event)
{

}

void PaintCanvas::paint(Graphics &g)
{
  g.fillAll(Colours::black); // preliminary
}

void PaintCanvas::resized()
{

}

void PaintCanvas::paintDot(int x, int y)
{

}

//=================================================================================================

PainterComponent::PainterComponent()
{
  addAndMakeVisible(canvas);


  sliderSize.setSliderName("Size");
  sliderSize.setRange(0.0, 30.0, 0.125, 3.0);
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

}