//-------------------------------------------------------------------------------------------------
// construction/destruction:
/*
MeteringDisplay::MeteringDisplay(const juce::String& description) : RWidget(description)
{

}

MeteringDisplay::~MeteringDisplay()
{
  deleteAllChildren();  // can't we leave that to the baseclass?
}
*/

//-------------------------------------------------------------------------------------------------
// setup:

void MeteringDisplay::setMeterStyle(int newMeterStyle)
{
  jassert(newMeterStyle >= 0 && newMeterStyle < numMeterStyles); // unknown meter style
  if( newMeterStyle >= 0 && newMeterStyle < numMeterStyles )
    style = newMeterStyle;
}

void MeteringDisplay::setRange(float newMinimum, float newMaximum)
{
  jassert(newMinimum < newMaximum); // max must be > min, equality not allowed due to div-by-zero
  if(newMinimum < newMaximum) {
    minValue = newMinimum;
    maxValue = newMaximum; }
}

void MeteringDisplay::setReferenceValue(float newReferenceValue)
{
  referenceValue = newReferenceValue;
}

void MeteringDisplay::setCurrentValue(float newValue)
{
  if(newValue != currentValue) {
    currentValue = jlimit(minValue, maxValue, newValue);
    repaintOnMessageThread(); }
}

//-------------------------------------------------------------------------------------------------
// overrides:

void MeteringDisplay::paint(Graphics &g)
{
  //g.fillAll(Colours::black);
  float x = 0.f;
  float y = 0.f;
  float w = (float) getWidth();
  float h = (float) getHeight();
  float relVal = (currentValue-minValue) / (maxValue-minValue); // relative value
  relVal = jlimit(0.f, 1.f, relVal);


  auto createGradientFill = [&](bool vertical)
  {
    int x1, y1, x2, y2; 
    
    if(vertical)
    {
      x1 = 0;  y1 = getHeight();  // left-bottom -> green
      x2 = 0;  y2 = 0;            // left-top    -> magenta
      // When x1 = x2, we get a purely vertical gradient.
    }
    else
    {
      y1 = 0;  x1 = 0;
      y2 = 0;  x2 = getWidth();
      // When y1 = y2, we get a purely horizontal gradient.
    }

    ColourGradient gradient = ColourGradient(Colours::green, x1, y1, Colours::magenta, x2, y2, false);
    float s = (referenceValue-minValue) / (maxValue-minValue);
    gradient.addColour(     s, Colours::red);
    gradient.addColour(0.75*s, Colours::yellow);
    FillType fill(gradient);
    g.setFillType(fill);
    g.fillAll();
  };


  switch( style )
  {
  case levelMeterStyle:
    {
      // Fill the whole meter with the gradient (maybe we should optimize this by caching the 
      // gradient as image?):
      createGradientFill(isVertical());

      // Cover some some part of the gradient with the background color:
      g.setColour(getBackgroundColour());
      if(isVertical())
        g.fillRect(0.f,        0.f, w, h * (1.f - relVal));
      else
        g.fillRect(w * relVal, 0.f, w * (1.f - relVal), h);
    }
    break;
  case triangularPointerStyle:
    {
      // Fill the background:
      g.setColour(getBackgroundColour());
      g.fillRect(x, y, w, h);

      // Draw the indicator:
      if(isVertical())
      {
        float pos = h * relVal;
        g.setColour(getWeakHighlightColour());
        drawTriangle(g, x, h-pos-w/2.f, x, h-pos+w/2.f, w, h-pos, true);
        g.setColour(getStrongHighlightColour());
        g.drawLine(x, h-pos, w, h-pos);
      }
      else
      {
        float pos = w * (1.f - relVal);
        g.setColour(getWeakHighlightColour());
        drawTriangle(g, w-pos-h/2.f, y, w-pos+h/2.f, y, w-pos, h, true);
        g.setColour(getStrongHighlightColour());
        g.drawLine(w-pos, y, w-pos, h);
      }
    }
    break;
  case horizontalRatio:
  {
    // under construction....

    // Fill the background:
    g.setColour(getBackgroundColour());
    g.fillRect(x, y, w, h);

    // Draw the indicator bar:
    float right = w * relVal;
    g.setColour(getWeakHighlightColour());
    g.fillRect(x, y, right, h);

    // Maybe do these in a subclass MeteringDisplayWithText
    // Display name of indicator at left:

    // Display formatted numeric value at right:
    // ...
  }
  break;


  } // end of switch(style)
}

//=================================================================================================

MeteringDisplayWithText::MeteringDisplayWithText() 
{ 
  style = horizontalRatio;
}

void MeteringDisplayWithText::setMeasurementName(const juce::String& newName) 
{ 
  measurementName = newName;
  // Should we trigger a repaint here? This may be appropriate, if the name is supposed to be 
  // dynamically updated. We currently don't need such a thing anywhere, so for the moment, it's 
  // good enough as is.
}

void MeteringDisplayWithText::setStringConversion(juce::String (*f) (double val, double max))
{
  valueToString = f;
  // Similar considerations with regard to repainting apply as for the name. The string conversion 
  // might be more likely to change at runtime though - for example to switch between units 
  // (although such switches can be handled within the conversion function, too).
}

void MeteringDisplayWithText::paint(Graphics& g)
{
  MeteringDisplay::paint(g);
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  int m = 4;              // margin

  // Draw the name:
  //y = handleRectangle.getY() + handleRectangle.getHeight()/2 - font->getFontAscent()/2;
  y = (h - font->getFontAscent())/2;
  drawBitmapFontText(g, x+m, y, measurementName, font, getTextColour());

  // Draw the value:
  String valueString = valueToString(currentValue, maxValue);
  x = getWidth() - font->getTextPixelWidth(valueString, font->getDefaultKerning());
  drawBitmapFontText(g, x-m, y, valueString, font, getTextColour());
}

