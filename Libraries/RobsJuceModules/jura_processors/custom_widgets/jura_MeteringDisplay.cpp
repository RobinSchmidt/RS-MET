//-------------------------------------------------------------------------------------------------
// construction/destruction:

MeteringDisplay::MeteringDisplay(const juce::String &componentName) : RWidget(componentName)
{
  //style          = levelMeterStyle;
  //verticalMode   = true;
  //minValue       = -48.0;
  //maxValue       = +6.0;
  //currentValue   = +0.0;
  //referenceValue = 0.0;
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

  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  // Factor out into drawVertically(g)
  // calculate the positions of the reference- and the current value in component coordinates:
  //double referencePosition = h * (referenceValue-minValue) / (maxValue-minValue);

  float relVal = (currentValue-minValue) / (maxValue-minValue); // relative value
  
  switch( style )
  {
  case levelMeterStyle:
    {
      // Create and paint the gradient:
      ColourGradient gradient = ColourGradient(Colours::green, (float) w, (float) h,
        Colours::magenta, (float) x, (float) y, false);
      gradient.addColour(     (referenceValue-minValue) / (maxValue-minValue), Colours::red);
      gradient.addColour(0.75*(referenceValue-minValue) / (maxValue-minValue), Colours::yellow);
      FillType fill(gradient);
      g.setFillType(fill);
      g.fillAll();

      // draw a black rectangle over some part of the gradient:
      float top = h * relVal;
      g.setColour(Colours::black);
      g.fillRect((float) x, (float) y, (float) w,
        jlimit((float) x, (float) h, (float) (h-top)) );
    }
    break;
  case triangularPointerStyle:
    {
      float pos = h * relVal;

      g.setColour(Colours::blue.brighter(2.0f)); // use some color from the widget-colorscheme

      drawTriangle(g, (float) x, (float) (h-pos-w/2),
                      (float) x, (float) (h-pos+w/2),
                      (float) w, (float) (h-pos), true);

      g.setColour(Colours::white);  // ...here too
      g.drawLine((float) x, (float) (h-pos), (float) w, (float) (h-pos));
    }
    break;
  case horizontalRatio:
  {

  }
  break;


  } // end of switch(style)
}

