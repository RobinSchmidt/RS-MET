// class ParameterObserver:

bool ParameterObserver::globalAutomationSwitch = true;
bool ParameterObserver::guiAutomationSwitch    = true;

ParameterObserver::ParameterObserver()
{

}

ParameterObserver::~ParameterObserver()
{

}

bool ParameterObserver::wantsAutomationNotification()
{
  bool result = localAutomationSwitch && globalAutomationSwitch;
  if( this->isGuiElement() )
    return result && guiAutomationSwitch;
  else
    return result;
}

//=================================================================================================
// class Parameter:

bool Parameter::storeDefaultValues = false;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Parameter::Parameter(const juce::String& newName, double newMin, double newMax,
  double newDefault, Scaling newScaling, double newInterval)
{
  jassert(!(newMin >= newMax));                           // invalid range
  jassert(!(newMin <= 0.0 && newScaling == EXPONENTIAL)); // exponential scaling requires strictly positive minimum value

  if(newScaling == STRING)
    newInterval = 1.0;      // multiple choice are parameters represented by integers

  mapper   = new rsParameterMapperLinear();
  scaling  = newScaling;    // delete soon - hmm...or maybe not? Not sure. Maybe we should keep it
  setRange(newMin, newMax);
  setScaling(newScaling);

  setName(newName);
  interval        = newInterval;
  defaultValue    = restrictValueToParameterRange(newDefault);
  value           = defaultValue;
  normalizedValue = valueToProportion(value);

  valueChangeCallbackFunction = [](double) {};
}

Parameter::Parameter(CriticalSection *criticalSectionToUse, const String& newName,
  double newMinValue, double newMaxValue, double newInterval, double newDefaultValue,
  Scaling newScaling)
{
  jassert(!(newMinValue >= newMaxValue));                      // that would result in a zero or negative range
  jassert(!(newMinValue <= 0.0 && newScaling == EXPONENTIAL)); // exponential scaling requires strictly positive minimum value

  mapper   = new rsParameterMapperLinear();
  scaling  = newScaling;  // delete soon
  setRange(newMinValue, newMaxValue);
  setScaling(newScaling);  // set this before restrictValueToParameterRange is called

  mutex         = criticalSectionToUse;
  name          = newName;
  interval      = newInterval;
  defaultValue  = restrictValueToParameterRange(newDefaultValue);
  value         = defaultValue;
  normalizedValue = valueToProportion(value);

  valueChangeCallbackFunction = [](double) {};
}

Parameter::~Parameter()
{
  ScopedPointerLock spl(mutex);

  clearValueChangeCallbacks();

  while( parameterObservers.size() > 0 )
  {
    ParameterObserver *o = parameterObservers[0];
    o->parameterWillBeDeleted(this);
    deRegisterParameterObserver(o);
  }
  // remark: we use a while-loop to account for the possibility that the observer de-registers
  // itself in the callback to parameterWillBeDeleted in which case the array-size shrinks
  // inside the iteration which would make a for-loop ...mmm...a bug

  delete mapper;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void Parameter::setValue(double newValue, bool sendNotification, bool callCallbacks)
{
  ScopedPointerLock spl(mutex);
  if(value == newValue)
    return;
  value = restrictValueToParameterRange(newValue);
  normalizedValue = valueToProportion(value);
  if( callCallbacks == true )
    callValueChangeCallbacks(value);
  if( sendNotification == true )
    notifyObservers();
}

void Parameter::setNormalizedValue(double newValue, bool sendNotification, bool callCallbacks)
{
  ScopedPointerLock spl(mutex);
  setValue(proportionToValue(newValue), sendNotification, callCallbacks);
}

void Parameter::setRangeAndValue(double newMin, double newMax, double newValue,
  bool sendNotification, bool callCallbacks)
{
  ScopedPointerLock spl(mutex);
  jassert(newMin <= newValue && newValue <= newMax); // inconsistent values

  mapper->setRange(newMin, newMax);
  value = restrictValueToParameterRange(newValue);
  normalizedValue = valueToProportion(value);

  if( callCallbacks == true )    callValueChangeCallbacks(value);
  if( sendNotification == true ) notifyObservers();
}

void Parameter::resetToDefaultValue(bool sendNotification, bool callCallbacks)
{
  ScopedPointerLock spl(mutex);
  setValue(defaultValue, sendNotification, callCallbacks);
}

void Parameter::setStringValue(const juce::String &newString, bool sendNotification,
  bool callCallbacks)
{
  ScopedPointerLock spl(mutex);
  for(int i = 0; i < stringValues.size(); i++)
  {
    if( stringValues[i] == newString )
    {
      setValue((double) i, sendNotification, callCallbacks);
      return;
    }
  }
  setValue(defaultValue, sendNotification, callCallbacks);
  // use the default option when the string was not found
}

void Parameter::setRange(double newMinValue, double newMaxValue)
{
  ScopedPointerLock spl(mutex);
  mapper->setRange(newMinValue, newMaxValue);
  normalizedValue = mapper->unmap(value);
  valueSanityCheck();
  for(int i=0; i < (int) parameterObservers.size(); i++)
  {
    if( parameterObservers[i]->wantsAutomationNotification() )
      parameterObservers[i]->parameterRangeChanged(this);
  }
}

void Parameter::setMinValue(double newMinValue)
{
  if( newMinValue <= getMaxValue() )
    setRange(newMinValue, getMaxValue());
  else
    setRange(getMaxValue(), getMaxValue());    // triggers jassert
}

void Parameter::setMaxValue(double newMaxValue)
{
  if( newMaxValue >= getMinValue() )
    setRange(getMinValue(), newMaxValue);
  else
    setRange(getMinValue(), getMinValue());    // triggers jassert
}

void Parameter::setDefaultValue(double newDefaultValue, bool setToDefault)
{
  ScopedPointerLock spl(mutex);
  defaultValue = restrictValueToParameterRange(newDefaultValue);
  if( setToDefault == true )
    setValue(defaultValue, true, true);
}

void Parameter::setScaling(Scaling newScaling)
{
  ScopedPointerLock spl(mutex);
  if( newScaling < IDENTITY || newScaling >= NUM_SCALINGS )
  {
    jassertfalse;
    return; // invalid scaling index
  }
  scaling = newScaling;

  // create a new mapper object, if necessarry (maybe factor out):
  if(scaling == EXPONENTIAL)
  {
    rsParameterMapperExponential* tmp = dynamic_cast<rsParameterMapperExponential*>(mapper);
    if(tmp == nullptr) // old mapper is the wrong kind, we need to create a new one
    {
      tmp = new rsParameterMapperExponential;
      tmp->setRange(mapper->getMin(), mapper->getMax());
      delete mapper;
      mapper = tmp;
    }
  }
  else if(scaling == IDENTITY)
  {
    rsParameterMapperIdentity* tmp = dynamic_cast<rsParameterMapperIdentity*>(mapper);
    if(tmp == nullptr) 
    {
      delete mapper;
      mapper = new rsParameterMapperIdentity;
    }
  }
  else
  {
    // default case, catches LINEAR, LINEAR_BIPOLAR, INTEGER, BOOLEAN, STRING, uses linear mapper
    rsParameterMapperLinear* tmp = dynamic_cast<rsParameterMapperLinear*>(mapper);
    if(tmp == nullptr) 
    {
      tmp = new rsParameterMapperLinear;
      tmp->setRange(mapper->getMin(), mapper->getMax());
      delete mapper;
      mapper = tmp; 
    }
  }
}

void Parameter::setScalingFromString(String newScalingString)
{
  ScopedPointerLock spl(mutex);
  if( newScalingString == String("Boolean") )
    setScaling(Scaling::BOOLEAN);
  else if( newScalingString == String("Integer") )
    setScaling(Scaling::INTEGER);
  else if( newScalingString == String("Linear") )
    setScaling(Scaling::LINEAR);
  else if( newScalingString == String("Exponential") )
    setScaling(Scaling::EXPONENTIAL);
  else if( newScalingString == String("LinearBipolar") )
    setScaling(Scaling::LINEAR_BIPOLAR);
  else
    setScaling(Scaling::LINEAR);
}

void Parameter::setMapper(rsParameterMapper* newMapper)
{
  scaling = Scaling::CUSTOM;
  if(mapper != newMapper) {
    delete mapper;
    mapper = newMapper;
    normalizedValue = valueToProportion(value); // make normalized and actual value consistent
  }
}

void Parameter::addStringValue(const String& valueToAdd)
{
  ScopedPointerLock spl(mutex);
  stringValues.addIfNotAlreadyThere(valueToAdd);
  mapper->setRange(0.0, jmax(double(stringValues.size()-1), DBL_MIN)); // max must be > min
  //minValue = 0.0;
  //maxValue = (double) (stringValues.size()-1);
}

void Parameter::addNumericStringValues(int min, int max, int step)
{
  for(int i = min; i <= max; i++)
    addStringValue(String(i));
}

//-------------------------------------------------------------------------------------------------
// inquiry:

double Parameter::valueToProportion(double value) const
{
  return mapper->unmap(value);

  // old:
  /*
  if(minValue >= maxValue)
    return 0.0;
  switch( scaling )
  {
  case Parameter::EXPONENTIAL:
  {
    if( minValue > 0.0 )
      return jlimit(0.0, 1.0, log(value/minValue) / (log(maxValue/minValue)) );
    else
      return 0.0;
  }
  default: return (value - minValue) / (maxValue - minValue); // LINEAR(_BIPOLAR)
  }
  */
}

double Parameter::proportionToValue(double prop) const
{
  return mapper->map(prop);

  // old:
  /*
  switch( scaling )
  {
  case EXPONENTIAL: return minValue * exp(prop*(log(maxValue/minValue)));
  default:          return minValue + (maxValue - minValue) * prop;
  }
  */
}

String Parameter::getStringValue() const
{
  ScopedPointerLock spl(mutex);
  if( this->isStringParameter() )
  {
    // int index = (int) getValue(); // doesn't work for subclass ParameterGridInterval because
    int index = (int) value;         // it wrangles the value, so we need this
    if( index >=0 && index < stringValues.size() )
      return stringValues[index];
    else
    {
      jassertfalse; // value does represent a valid index in the array of string values
      return String();
    }
  }
  else
    return String();
}

String Parameter::getOptionStringAtIndex(int index) const
{
  ScopedPointerLock spl(mutex);
  if( !isStringParameter() )
  {
    jassertfalse; // tyring to retrieve an option-string from a non-string parameter
    return String();
  }

  if( index >= 0 && index < stringValues.size() )
    return stringValues[index];
  else
  {
    jassertfalse; // tyring to retrieve an option-string with invalid index
    return String();
  }
}

String Parameter::getDefaultStringValue() const
{
  ScopedPointerLock spl(mutex);
  if( !isStringParameter() )
    return String();
  else
  {
    int index = (int) getDefaultValue();
    if( index >= 0 && index < stringValues.size() )
      return stringValues[index];
    else
      return String();
  }
}

String Parameter::getScalingString() const
{
  ScopedPointerLock spl(mutex);
  switch( scaling )
  {
  case BOOLEAN:        return String("Boolean");
  case INTEGER:        return String("Integer");
  case LINEAR:         return String("Linear");
  case EXPONENTIAL:    return String("Exponential");
  case LINEAR_BIPOLAR: return String("LinearBipolar");
  default:             return String("Linear");
  }
}

//-------------------------------------------------------------------------------------------------
// state recall:

void Parameter::saveToXml(XmlElement* xml) const
{
  if(shouldBeSavedAndRecalled())
  {
    if(!isCurrentValueDefaultValue() || storeDefaultValues == true)
    {
      if(isStringParameter())
        xml->setAttribute(getName(), getStringValue());
      else
        xml->setAttribute(getName(), juce::String(getValue()));
    }
  }
}

void Parameter::recallFromXml(const XmlElement& xml)
{
  if(shouldBeSavedAndRecalled()) 
  {
    if(isStringParameter())
      setStringValue(xml.getStringAttribute(getName(), getDefaultStringValue()), true, true);
    else
    {
      double dbg = xml.getDoubleAttribute(getName(), INF);
      setValue(xml.getDoubleAttribute(getName(), getDefaultValue()), true, true);
    }
  }
}

//-------------------------------------------------------------------------------------------------
// add/remove parameterObservers:

void Parameter::registerParameterObserver(ParameterObserver *observerToAdd)
{
  ScopedPointerLock spl(mutex);
  appendIfNotAlreadyThere(parameterObservers, observerToAdd);
  //parameterObservers.addIfNotAlreadyThere(observerToAdd);
}

void Parameter::deRegisterParameterObserver(ParameterObserver *observerToRemove)
{
  ScopedPointerLock spl(mutex);
  removeFirstOccurrence(parameterObservers, observerToRemove);
  //parameterObservers.removeFirstMatchingValue(observerToRemove);
}

//-------------------------------------------------------------------------------------------------
// others:

void Parameter::notifyObservers()
{
  ScopedPointerLock spl(mutex);
  for(int i = 0; i < (int) parameterObservers.size(); i++)
  {
    //ParameterObserver *observer = parameterObservers[i];
    if( parameterObservers[i]->wantsAutomationNotification() )
      parameterObservers[i]->parameterChanged(this);
  }
}
/*
void Parameter::notifyGuiObservers()
{
  ScopedPointerLock spl(mutex);
  for(int i = 0; i < (int) parameterObservers.size(); i++)
  {
    if( parameterObservers[i]->isGuiElement() && 
      parameterObservers[i]->wantsAutomationNotification() )
      parameterObservers[i]->parameterChanged(this);
  }
}

void Parameter::notifyNonGuiObservers()
{
  ScopedPointerLock spl(mutex);
  for(int i = 0; i < (int) parameterObservers.size(); i++)
  {
    if( !parameterObservers[i]->isGuiElement() && 
      parameterObservers[i]->wantsAutomationNotification() )
      parameterObservers[i]->parameterChanged(this);
  }
}
*/
void Parameter::callValueChangeCallbacks(double argument)
{
  ScopedPointerLock spl(mutex);
  if(valueChangeCallbackFunction)
	  valueChangeCallbackFunction(argument);
  if( valueChangeCallbackDouble != nullptr )
    valueChangeCallbackDouble->call(argument);
  if( valueChangeCallbackInt != nullptr )
    valueChangeCallbackInt->call(juce::roundToInt(argument));
  if( valueChangeCallbackBool != nullptr )
    valueChangeCallbackBool->call(argument >= 0.5);
}

double Parameter::restrictValueToParameterRange(double valueToRestrict)
{
  ScopedPointerLock spl(mutex);

  // quantize:
  if( interval > 0 )
    valueToRestrict = getMinValue() 
    + interval * floor((valueToRestrict - getMinValue()) / interval + 0.5);
  if( scaling == BOOLEAN )
    valueToRestrict = (double) (valueToRestrict >= 0.5);
  if( scaling == STRING || scaling == INTEGER )
    valueToRestrict = (double) round(valueToRestrict);

  // clip:
  if( valueToRestrict > getMaxValue() )
    return getMaxValue();
  if( valueToRestrict < getMinValue() )
    return getMinValue();

  return valueToRestrict;
}

void Parameter::valueSanityCheck()
{
  ScopedPointerLock spl(mutex);
  jassert( !( getMinValue() <= 0.0 && scaling == EXPONENTIAL ) );

  /*
  // todo: update this for use with new mapper object:
  if((getMinValue() <= 0.0 && scaling == EXPONENTIAL))
  {
    minValue = 0.1;
  }
  jassert( getMaxValue() > getMinValue() );  // maybe allow >=
  if(getMaxValue() <= getMinValue())
  {
    maxValue = getMinValue+1.0;
  }
  */
  // updated:
  jassert( getMaxValue() >= getMinValue() ); 


  value        = restrictValueToParameterRange(value);
  defaultValue = restrictValueToParameterRange(defaultValue);
}


/*
Ideas:

\todo: getStorageString, setFromStorageString, getDisplayString
-allows for storing string/enum parameters with different stings than what is displayed in the
 ComboBox - example:
 displayString == 24 dB/oct (if we are already in a submenu 'Lowpass')
 storageString == Lowpass24 or Lowpass>24
-include a description-string ...will be copied into the corresponding widget's description on
 assignment

\todo
-for performance, provide two versions of setValue - the regular one and a version
 setValueWithoutLocking - where the latter should be called only under the premise that the
 mutex-lock is already held by the caller

the set/getValue stuff is messed up when smoothing/meta-control/modulation comes into play
ensure the following:
-set/get(Normalized)Value: sets/returns unmodulated (normalized) target value in all of Parameter's
 subclasses
 -setNormalizedValue is called from sliders in mouseDrag, etc., getNormalizedValue is called in 
  paint, so all subclasses must make sure to return a value that is suitable for that purpose
 -getValue is also called in paint and used for the numeric readout - so subclasses must make
  sure that get getValue returns the mapped value without smoothing and modulation



*/


//=================================================================================================
// class ParameterSetHolder:

void ParameterSetHolder::registerParameterSetObserver(ParameterSetObserver *observerToRegister)
{
  appendIfNotAlreadyThere(parameterSetObservers, observerToRegister);
}

void ParameterSetHolder::deRegisterParameterSetObserver(ParameterSetObserver *observerToDeRegister)
{
  removeFirstOccurrence(parameterSetObservers, observerToDeRegister);
}

void ParameterSetHolder::sendParameterSetChangeNotification(
  ParameterSetHolder* parameterSetHolderThatHasChanged)
{
  for(int i = 0; i < size(parameterSetObservers); i++)
    parameterSetObservers[i]->parameterSetChanged(parameterSetHolderThatHasChanged);
}
