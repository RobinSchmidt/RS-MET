//-------------------------------------------------------------------------------------------------
// construction/destruction:

MeteringDisplay::MeteringDisplay(const juce::String &componentName) : RWidget(componentName)
{

}

MeteringDisplay::~MeteringDisplay()
{
  deleteAllChildren();  // can't we leave that to the baseclass?
}

//-------------------------------------------------------------------------------------------------
// setup:

void MeteringDisplay::setMeterStyle(int newMeterStyle)
{
  // the parameter 'newMeterStyle' seems not to be one of the predefined styles
  jassert(newMeterStyle >= 0 && newMeterStyle < numMeterStyles);

  if( newMeterStyle >= 0 && newMeterStyle < numMeterStyles )
    style = newMeterStyle;
}

void MeteringDisplay::setRange(float newMinimum, float newMaximum)
{
  jassert(newMinimum < newMaximum);
  if(newMinimum < newMaximum)
  {
    minValue = newMinimum;
    maxValue = newMaximum;
  }
}

void MeteringDisplay::setReferenceValue(float newReferenceValue)
{
  referenceValue = newReferenceValue;
}

void MeteringDisplay::setCurrentValue(float newValue)
{
  currentValue = newValue;
  repaint();
}

//-------------------------------------------------------------------------------------------------
// overrides:

void MeteringDisplay::paint(Graphics &g)
{
  g.fillAll(Colours::black);

  float x = 0.f;
  float y = 0.f;
  float w = (float) getWidth();
  float h = (float) getHeight();

  float relVal = (currentValue-minValue) / (maxValue-minValue); // relative value

  switch( style )
  {
  case levelMeterStyle:
    {
      // Create and paint the gradient:
      ColourGradient gradient = ColourGradient(Colours::green, w,h, Colours::magenta, x,y, false);
      gradient.addColour(     (referenceValue-minValue) / (maxValue-minValue), Colours::red);
      gradient.addColour(0.75*(referenceValue-minValue) / (maxValue-minValue), Colours::yellow);
      FillType fill(gradient);
      g.setFillType(fill);
      g.fillAll();

      // Cover some some part of the gradient with the background color:
      float top = h * relVal;
      g.setColour(getBackgroundColour());
      g.fillRect(x, y, w, jlimit(x, h, h-top));
    }
    break;
  case triangularPointerStyle:
    {
      // Fill the background:
      g.setColour(getBackgroundColour());
      g.fillRect(x, y, w, h);

      // Draw the indicator:
      float pos = h * relVal;
      g.setColour(getWeakHighlightColour()); 
      drawTriangle(g, x, h-pos-w/2.f, x, h-pos+w/2.f, w, h-pos, true);
      g.setColour(getStrongHighlightColour());
      g.drawLine(x, h-pos, w, h-pos);
    }
    break;
  case horizontalRatio:
  {
    // under construction....
  }
  break;


  } // end of switch(style)
}
