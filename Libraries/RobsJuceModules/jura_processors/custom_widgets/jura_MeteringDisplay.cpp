//-------------------------------------------------------------------------------------------------
// construction/destruction:

MeteringDisplay::MeteringDisplay(const juce::String &componentName) : RWidget(componentName)
{
  style          = levelMeterStyle;
  verticalMode   = true;
  minValue       = -48.0;
  maxValue       = +6.0;
  currentValue   = +0.0;
  referenceValue = 0.0;
}

MeteringDisplay::~MeteringDisplay()
{
  deleteAllChildren();
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

void MeteringDisplay::setRange(double newMinimum, double newMaximum)
{
  jassert(newMinimum < newMaximum);
  if(newMinimum < newMaximum)
  {
    minValue = newMinimum;
    maxValue = newMaximum;
  }
}

void MeteringDisplay::setReferenceValue(double newReferenceValue)
{
  referenceValue = newReferenceValue;
}

void MeteringDisplay::setCurrentValue(double newValue)
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

  // calculate the positions of the reference- and the current value in component coordinates:
  //double referencePosition = h * (referenceValue-minValue) / (maxValue-minValue);
  double currentPosition   = h * (currentValue-minValue) / (maxValue-minValue);

  switch( style )
  {
  case levelMeterStyle:
    {
      // create and paint the gradient:
      ColourGradient gradient = ColourGradient(Colours::green, (float) w, (float) h,
        Colours::magenta, (float) x, (float) y, false);
      gradient.addColour(     (referenceValue-minValue) / (maxValue-minValue), Colours::red);
      gradient.addColour(0.75*(referenceValue-minValue) / (maxValue-minValue), Colours::yellow);

      //GradientBrush  brush(gradient);
      //g.setBrush(&brush);
      FillType fill(gradient);
      g.setFillType(fill);

      g.fillAll();

      // draw a black rectangle over some part of the gradient:
      g.setColour(Colours::black);
      g.fillRect((float) x, (float) y, (float) w,
        jlimit((float) x, (float) h, (float) (h-currentPosition)) );
    }
    break;
  case triangularPointerStyle:
    {
      g.setColour(Colours::blue.brighter(2.0f));
      drawTriangle(g, (float) x, (float) (h-currentPosition-w/2),
                      (float) x, (float) (h-currentPosition+w/2),
                      (float) w, (float) (h-currentPosition), true);

      g.setColour(Colours::white);
      g.drawLine((float) x, (float) (h-currentPosition), (float) w, (float) (h-currentPosition));
    }
    break;
  } // end of switch(style)
}

void MeteringDisplay::setBounds(int x, int y, int w, int h)
{
  Component::setBounds(x, y, w, h);
}
