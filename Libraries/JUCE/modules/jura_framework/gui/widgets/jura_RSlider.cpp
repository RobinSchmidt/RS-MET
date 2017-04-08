//-------------------------------------------------------------------------------------------------
// construction/destruction:

RSlider::RSlider(const String& name)
{
  currentValue             = 0.5;
  defaultValue             = 0.5;
  minValue                 = 0.0;
  maxValue                 = 1.0;
  interval                 = 0.01;
  stringConversionFunction = &valueToString2;
  scaling                  = Parameter::LINEAR;
  layout                   = NAME_INSIDE;

  addAndMakeVisible( nameRectangle = new Component() );
  nameRectangle->addMouseListener(this, true);
  setWantsKeyboardFocus(true); // we need this for the mousewheel to work

  ParameterObserver::localAutomationSwitch = true; // do we need this?
}

RSlider::~RSlider()
{
  // remove ourselves as listener from the Parameter object, such that it does not try to notify a 
  // nonexistent listener:
  ParameterObserver::localAutomationSwitch = false;

  if( assignedParameter != NULL )
    assignedParameter->deRegisterParameterObserver(this);
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void RSlider::setSliderName(const String &newName)
{
  sliderName = newName;
}

void RSlider::setLayout(int newLayout)
{
  layout = newLayout;
  resized();
}

void RSlider::setRange(double newMin, double newMax, double newInt, double newDef, bool initToDefault)
{
  jassert(newMin <= newMax);
  jassert(newInt >= 0.0);

  minValue     = newMin;
  maxValue     = newMax;
  interval     = newInt;
  defaultValue = newDef;

  valueSanityCheck();

  // initialize to default-value or keep the current value inside the new range:
  if( initToDefault == true )
    setValue(defaultValue, false, false);
  else
    setValue(currentValue, false, false);

  repaint();
}

void RSlider::setScaling(int newScaling)
{
  if( newScaling >= Parameter::LINEAR && newScaling <= Parameter::LINEAR_BIPOLAR )
  {
    scaling = newScaling;
    valueSanityCheck();
    repaint();
  }

  //if( scaling == Parameter::EXPONENTIAL && minValue <= 0.0 )
  //  jassertfalse; // minValue must be strictly greater than 0 for exponential scaling
}

void RSlider::setValue(double newValue, const bool sendUpdateMessage, 
  const bool sendMessageSynchronously)
{
  newValue = constrainValue(newValue);
  valueSanityCheck();

  if(currentValue != newValue)
  {
    currentValue = newValue;

    // currently being tested:
    if( assignedParameter != NULL )
    {
      ParameterObserver::localAutomationSwitch = false;
      //assignedParameter->setValue(currentValue);  //old
      //assignedParameter->setValue(currentValue, false); 
      assignedParameter->setValue(currentValue, sendUpdateMessage, sendUpdateMessage); 
       // is the thing with the sendUpdateMessage correct?
      ParameterObserver::localAutomationSwitch = true;
    }
    repaint();

    if (sendUpdateMessage)
      triggerChangeMessage (sendMessageSynchronously);
  }
}

void RSlider::setStateFromString(const juce::String &stateString, bool sendChangeMessage)
{
  double value = stateString.getDoubleValue();
  setValue(value, sendChangeMessage);
}

void RSlider::setProportionalValue(double newProportionalValue, 
      const bool sendUpdateMessage, const bool sendMessageSynchronously)
{
  setValue(proportionOfLengthToValue(newProportionalValue), sendUpdateMessage, 
    sendMessageSynchronously);
}

void RSlider::setDefaultValue(double newDefaultValue)
{
  newDefaultValue = constrainValue(newDefaultValue);
  defaultValue    = newDefaultValue;
  valueSanityCheck();
  if( assignedParameter != NULL )
    assignedParameter->setDefaultValue(defaultValue);
  repaint();
}

void RSlider::setDefaultValues(juce::Array<double> newDefaultValues)
{
  defaultValues = newDefaultValues;
}

void RSlider::setToDefaultValue(const bool sendUpdateMessage, const bool sendMessageSynchronously)
{
  setValue(defaultValue, sendUpdateMessage, sendMessageSynchronously);
}

void RSlider::assignParameter(Parameter *parameterToAssign)
{
  RWidget::assignParameter(parameterToAssign);
  if( assignedParameter != NULL )
    parameterRangeChanged(assignedParameter);
}

void RSlider::parameterRangeChanged(Parameter* parameterThatHasChangedRange)
{
  if( parameterThatHasChangedRange == assignedParameter && assignedParameter != NULL )
  {
    currentValue  = assignedParameter->getValue();
    defaultValue  = assignedParameter->getDefaultValue();
    defaultValues = assignedParameter->getDefaultValues();  
    minValue      = assignedParameter->getMinValue();
    maxValue      = assignedParameter->getMaxValue();
    scaling       = assignedParameter->getScaling();
    interval      = assignedParameter->getInterval();
    setSliderName(assignedParameter->getName());
    valueSanityCheck();
    //updateWidgetFromAssignedParameter(false);
    repaint();
  }
}

void RSlider::copySettingsFrom(const RSlider* otherSlider)
{
  sliderName  = otherSlider->sliderName;
  description = otherSlider->description;
  stringConversionFunction = otherSlider->stringConversionFunction;
  setRange(otherSlider->minValue, otherSlider->maxValue, otherSlider->interval, 
    otherSlider->defaultValue, true);
  scaling = otherSlider->scaling;

  // ...maybe more to add  here later
}

void RSlider::updateWidgetFromAssignedParameter(bool sendChangeMessage)
{
  if( assignedParameter != NULL )
  {
    setRange(assignedParameter->getMinValue(), assignedParameter->getMaxValue(), 
             assignedParameter->getInterval(), assignedParameter->getDefaultValue(), false);
    setValue(assignedParameter->getValue(), sendChangeMessage, false);
  }
}

void RSlider::setStringConversionFunction(String (*newConversionFunction) 
                                          (double valueToBeConverted) )
{
  stringConversionFunction = newConversionFunction;
  repaint();
}

const String RSlider::getTextFromValue(double value) const
{
  return stringConversionFunction(value);
}

double RSlider::getValueFromText (const String& text) const
{
  String t (text.trimStart());

  while (t.startsWithChar ('+'))
    t = t.substring (1).trimStart();

  return t.initialSectionContainingOnly("0123456789.-").getDoubleValue();
}

double RSlider::proportionOfLengthToValue(double proportion) const
{
  // todo:  use linToLin, linToExp, expToLin
  switch( scaling )
  {
  case Parameter::LINEAR:         return minValue + (maxValue - minValue) * proportion;
  case Parameter::EXPONENTIAL:    return minValue * exp( proportion*(log(maxValue)-log(minValue)) );
  case Parameter::LINEAR_BIPOLAR: return minValue + (maxValue - minValue) * proportion;

  default:          return 0.0;
  }
}

double RSlider::valueToProportionOfLength(double value) const
{
  // catch cases where minValue == maxValue
  // todo:  use linToLin, linToExp, expToLin

  switch( scaling )
  {
  case Parameter::LINEAR:      
    {
      if( maxValue > minValue )
        return (value - minValue) / (maxValue - minValue);
      else
        return 0.0;
    }
  case Parameter::EXPONENTIAL: 
    {
      if( minValue > 0.0 )
        return jlimit(0.0, 1.0, log(value/minValue) / (log(maxValue)-log(minValue)) );
      else
        return 0.0;
    }
  case Parameter::LINEAR_BIPOLAR:      
    {
      if( maxValue > minValue )
        return (value - minValue) / (maxValue - minValue);
      else
        return 0.0;
    }
  default:          return 0.0;
  }
}

juce::String RSlider::getStateAsString() const
{
  return juce::String(getValue());
}

void RSlider::addListener(RSliderListener* listener) throw()
{
  jassert( listener != 0 );
  if( listener != 0 )
    listeners.add(listener);
}

void RSlider::removeListener(RSliderListener* listener) throw()
{
  listeners.removeFirstMatchingValue(listener);
}

//-------------------------------------------------------------------------------------------------
// inquiry:

const String& RSlider::getSliderName() const
{
  return sliderName;
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void RSlider::mouseDown(const MouseEvent& e)
{
  if( e.originalComponent != this )
    return; // ignore click on the label, when it's not 'inside' the actual slider

  if( isEnabled() )
  {
    if( e.mods.isCommandDown() )
    {
      //if( assignedParameter != NULL )
      //  setValue( assignedParameter->getDefaultValue() );
      //else
      //  setValue(defaultValue);   // these default-values should actually be in sync...
      setValue(defaultValue);
    }
    else if( e.mods.isLeftButtonDown() )
    {
      double tmpValue = proportionOfLengthToValue((double) e.x / (double) getWidth());
      setValue(constrainAndQuantizeValue(tmpValue), true, false);
    }
  }
}

void RSlider::mouseDrag(const MouseEvent& e)
{
  if( e.originalComponent != this )
    return; // ignore drag on the label, when it's not 'inside' the actual slider

  double x = e.getMouseDownX() + e.getDistanceFromDragStartX();
  double tmpValue;
  if( isEnabled() )
  {
    if( !e.mods.isRightButtonDown() && !e.mods.isCommandDown() )
    {
      tmpValue = proportionOfLengthToValue((double) x / (double) getWidth());
      setValue(constrainAndQuantizeValue(tmpValue), true, false);
    }
  }
}

void RSlider::mouseDoubleClick(const MouseEvent& e)
{
  if( isEnabled() )
    setValue(openModalNumberEntryField(), true, false);
}

void RSlider::mouseWheelMove(const MouseEvent &event, const MouseWheelDetails &wheel)
{
  double tmpValue;
  if( isEnabled() )
  {
    float s;
    if( wheel.deltaY >= 0.0 )
      s = 1.0;
    else
      s = -1.0;

    float scale = 1.0;
    if(ModifierKeys::getCurrentModifiers().isCtrlDown())
      scale = 0.0125f;

    if( interval > 0.0 )
    {
      tmpValue = getValue() + scale * s * interval;
      setValue(constrainAndQuantizeValue(tmpValue), true, false);
    }
    else
      setProportionalValue(getProportionalValue() + scale * 0.01 * wheel.deltaY, true, false);
  }
}

void RSlider::paint(Graphics& g)
{
  int w = getWidth();
  int h = getHeight();

  if( w < 1 || h < 1 )
    return;

  g.setColour(getBackgroundColour());
  g.fillRect(handleRectangle);
  double thumbLeftValue, thumbRightValue, thumbWidth;
  g.setColour(getHandleColour());
  if( scaling == Parameter::LINEAR_BIPOLAR )
  {
    if( getValue() > defaultValue )
    {
      thumbLeftValue  = w * valueToProportionOfLength(defaultValue);
      thumbRightValue = w * valueToProportionOfLength(getValue());
    }
    else
    {
      thumbLeftValue  = w * valueToProportionOfLength(getValue());
      thumbRightValue = ::ceil(w * valueToProportionOfLength(defaultValue));
    }
    thumbWidth = thumbRightValue - thumbLeftValue;
    g.fillRect((int) (handleRectangle.getX() + thumbLeftValue), handleRectangle.getY(), 
      (int) thumbWidth, h);

    //g.setColour(thumbColour.brighter(0.25));
    //float x = (float) (w*valueToProportionOfLength(defaultValue));
    //g.drawLine(x, (float) sliderRectangle.getY(), x, (float) h, 2.f);
  }
  else
  {
    thumbWidth = w * valueToProportionOfLength(getValue());
    g.fillRect(handleRectangle.getX(), handleRectangle.getY(), (int) thumbWidth, h);
  }

  // draw the default value indicator:
  //g.setColour(specialColour1);
  g.setColour(getMixedColour(getBackgroundColour(), getHandleColour(), 0.25, 0.75));
  float xf = (float) (w*valueToProportionOfLength(defaultValue));
  g.drawLine(xf, (float) handleRectangle.getY(), xf, (float) h, 2.f);

  // draw the outline:
  g.setColour(getOutlineColour());
  g.drawRect(handleRectangle, 2);

  // draw the value:
  String valueString = stringConversionFunction(currentValue);
  int x = getWidth() - font->getTextPixelWidth(valueString, font->getDefaultKerning());
  int y = handleRectangle.getY() + handleRectangle.getHeight()/2 - font->getFontAscent()/2;
  drawBitmapFontText(g, x-4, y, valueString, font, getTextColour());

  // draw the name:
  x             = 4;
  Colour colour = getTextColour();
  if( layout == NAME_ABOVE )
  {
    x = 0;
    y = handleRectangle.getY() - font->getFontAscent() - 2;
    //colour = getOutlineColour();
    ColourSchemeComponent *colorParent = getOutlyingColourSchemeComponent();
    if( colorParent != NULL )
      colour = colorParent->getEditorColourScheme().headline;
  }
  drawBitmapFontText(g, x, y, sliderName, font, colour);

  // gray out the slider if it's disabled:
  //if( !isEnabled() )
  //  g.fillAll(Colours::lightgrey.withAlpha(0.75f));
}

void RSlider::resized()
{
  int w = getWidth();
  int h = getHeight();

  if( w < 1 || h < 1 )
    return;

  switch( layout )
  {
  case NAME_ABOVE:
    {
      nameRectangle->setBounds (0, 0,   w, h/2);
      handleRectangle.setBounds(0, h/2, w, h/2);
    };
    break;

  //case NAME_LEFT:

  case NAME_INSIDE:
    {
      handleRectangle.setBounds(0, 0, w, h);
      nameRectangle->setBounds( 0, 0, 0, 0);

      //int w2 = (int) (nameBoxRelativeWidth * (double) w);
      //nameBox->setBounds       (0,    0, w2, h);
      //int x = nameBox->getRight();
      //w2 = (int) (valueBoxRelativeWidth * (double) w);
      //valueBox->setBounds      (w-w2, 0, w2, h);
    };
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// internal functions:

double RSlider::constrainValue(double value) const throw()
{
  if( value <= minValue || maxValue <= minValue )
    value = minValue;
  else if( value >= maxValue )
    value = maxValue;

  if( assignedParameter != NULL )
  {
    value = jmax(value, assignedParameter->getMinValue());
    value = jmin(value, assignedParameter->getMaxValue());
  }

  return value;
}

double RSlider::constrainAndQuantizeValue(double value) const throw()
{
  if( interval > 0 )
    value = minValue + interval * floor ((value - minValue) / interval + 0.5);
  return constrainValue(value);
}

void RSlider::valueSanityCheck()
{
  jassert( !( minValue <= 0.0 && scaling == Parameter::EXPONENTIAL ) );
  if( ( minValue <= 0.0 && scaling == Parameter::EXPONENTIAL ) )
    minValue = 0.1;

  jassert( maxValue >= minValue );
  if( maxValue <= minValue )
    maxValue = minValue+1.0;

  currentValue = constrainValue(currentValue);
  defaultValue = constrainValue(defaultValue);
}

void RSlider::handleAsyncUpdate()
{
  RWidget::handleAsyncUpdate();
  cancelPendingUpdate();
  for(int i = listeners.size(); --i >= 0;)
  {
    ((RSliderListener*) listeners.getUnchecked (i))->rSliderValueChanged(this);
    i = jmin(i, listeners.size());
  }
}

void RSlider::triggerChangeMessage (const bool synchronous)
{
  if(synchronous)
    handleAsyncUpdate();
  else
    triggerAsyncUpdate();
}
