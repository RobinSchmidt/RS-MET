// class ParameterObserver:

bool ParameterObserver::globalAutomationSwitch = true;
bool ParameterObserver::guiAutomationSwitch    = true;

ParameterObserver::ParameterObserver()
{
  localAutomationSwitch  = false;
  isGuiElement           = false;
}

ParameterObserver::~ParameterObserver()
{

}

bool ParameterObserver::wantsAutomationNotification()
{
  bool result = localAutomationSwitch && globalAutomationSwitch;
  if( this->isGuiElement == true )
    return result && guiAutomationSwitch;
  else
    return result;
}

//=================================================================================================
// class Parameter:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Parameter::Parameter(const juce::String& newName, double newMin, double newMax,
  double newDefault, int newScaling, double newInterval)
{
  jassert(!(newMin >= newMax));                           // invalid range
  jassert(!(newMin <= 0.0 && newScaling == EXPONENTIAL)); // exponential scaling requires strictly positive minimum value

  name         = newName;
  minValue     = newMin;
  maxValue     = newMax;
  interval     = newInterval;
  scaling      = newScaling;
  defaultValue = restrictValueToParameterRange(newDefault);
  value        = defaultValue;
}

Parameter::Parameter(CriticalSection *criticalSectionToUse, const String& newName, 
  double newMinValue, double newMaxValue, double newInterval, double newDefaultValue, 
  int newScaling) 
{
  jassert(!(newMinValue >= newMaxValue));                      // that would result in a zero or negative range
  jassert(!(newMinValue <= 0.0 && newScaling == EXPONENTIAL)); // exponential scaling requires strictly positive minimum value

  mutex         = criticalSectionToUse;
  name          = newName;
  minValue      = newMinValue;
  maxValue      = newMaxValue;
  interval      = newInterval;
  scaling       = newScaling;  // set this before restrictValueToParameterRange is called
  defaultValue  = restrictValueToParameterRange(newDefaultValue);
  value         = defaultValue;
}

Parameter::~Parameter()
{
  ScopedPointerLock spl(mutex);

  clearValueChangeCallbacks();

  while( parameterObservers.size() > 0 )
  {
    ParameterObserver *o = parameterObservers[0];
    o->parameterIsGoingToBeDeleted(this);
    deRegisterParameterObserver(o);
  }
  // remark: we use a while-loop to account for the possibility that the observer de-registers 
  // itself in the callback to parameterIsGoingToBeDeleted in which case the array-size shrinks
  // inside the iteration which would make a for-loop ...mmm...a bug
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void Parameter::setValue(double newValue, bool sendNotification, bool callCallbacks)
{
  ScopedPointerLock spl(mutex);
  value = restrictValueToParameterRange(newValue);
  if( callCallbacks == true )
    callValueChangeCallbacks();
  if( sendNotification == true )
    notifyObservers();
}

void Parameter::setProportionalValue(double newProportionalValue,
  bool sendNotification, bool callCallbacks)
{
  ScopedPointerLock spl(mutex);
  setValue(proportionToValue(newProportionalValue), sendNotification, callCallbacks);
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
  minValue = newMinValue;
  maxValue = newMaxValue;
  valueSanityCheck();
  for(int i=0; i < (int) parameterObservers.size(); i++)
  {
    if( parameterObservers[i]->wantsAutomationNotification() )
      parameterObservers[i]->parameterRangeChanged(this);
  }
}

void Parameter::setMinValue(double newMinValue)
{
  if( newMinValue <= maxValue )
    setRange(newMinValue, maxValue);
  else
    setRange(maxValue, maxValue);
}

void Parameter::setMaxValue(double newMaxValue)
{
  if( newMaxValue >= minValue )
    setRange(minValue, newMaxValue);
  else
    setRange(minValue, minValue);
}

void Parameter::setDefaultValue(double newDefaultValue, bool setToDefault)
{
  ScopedPointerLock spl(mutex);
  defaultValue = restrictValueToParameterRange(newDefaultValue);
  if( setToDefault == true )
    setValue(defaultValue, true, true);
}

void Parameter::setScaling(int newScaling)
{
  ScopedPointerLock spl(mutex);
  if( newScaling < BOOLEAN || newScaling > EXPONENTIAL )
  {
    jassertfalse;
    return; // invalid scaling index
  }
  scaling = newScaling;
}

void Parameter::setScalingFromString(String newScalingString)
{
  ScopedPointerLock spl(mutex);
  if( newScalingString == String("Boolean") )
    setScaling(BOOLEAN);
  else if( newScalingString == String("Integer") )
    setScaling(INTEGER);
  else if( newScalingString == String("Linear") )
    setScaling(LINEAR);
  else if( newScalingString == String("Exponential") )
    setScaling(EXPONENTIAL);
  else if( newScalingString == String("LinearBipolar") )
    setScaling(LINEAR_BIPOLAR);
  else
    setScaling(LINEAR);
}

void Parameter::addStringValue(const String& valueToAdd)
{
  ScopedPointerLock spl(mutex);
  stringValues.addIfNotAlreadyThere(valueToAdd);
  minValue = 0.0;
  maxValue = (double) (stringValues.size()-1);
}

//-------------------------------------------------------------------------------------------------
// inquiry:

double Parameter::valueToProportion(double value)
{
  if(minValue >= maxValue)
    return 0.0;
  switch( scaling )
  {
  case Parameter::EXPONENTIAL: 
  {
    if( minValue > 0.0 )
      return jlimit(0.0, 1.0, log(value/minValue) / (log(maxValue)-log(minValue)) ); // optimize
    else
      return 0.0;
  }
  default: return (value - minValue) / (maxValue - minValue); // LINEAR(_BIPOLAR)
  }
}

double Parameter::proportionToValue(double prop)
{
  switch( scaling )
  {
  case Parameter::LINEAR:         return minValue + (maxValue - minValue) * prop;
  case Parameter::LINEAR_BIPOLAR: return minValue + (maxValue - minValue) * prop;  // is the same - get rid of duplication - use as default
  case Parameter::EXPONENTIAL:    return minValue * exp(prop*(log(maxValue)-log(minValue))); // optimize: log(a)-log(b) = log(a/b)
  default: return 0.0;
  }
  // maybe make valueToProportion/proportionToValue static methods and pass in min, max, scaling
}

String Parameter::getStringValue() const
{
  ScopedPointerLock spl(mutex);
  if( this->isStringParameter() )
  {
    int index = (int) getValue();
    if( index >=0 && index < stringValues.size() )
      return stringValues[index];
    else
    {
      jassertfalse; // value does represent a valid index in the array of string values
      return String::empty;
    }
  }
  else
    return String::empty;
}

String Parameter::getOptionStringAtIndex(int index) const
{
  ScopedPointerLock spl(mutex);
  if( !isStringParameter() )
  {
    jassertfalse; // tyring to retrieve an option-string from a non-string parameter
    return String::empty;
  }

  if( index >= 0 && index < stringValues.size() )
    return stringValues[index];
  else
  {
    jassertfalse; // tyring to retrieve an option-string with invalid index
    return String::empty;
  }
}

String Parameter::getDefaultStringValue() const
{
  ScopedPointerLock spl(mutex);
  if( !isStringParameter() )
    return String::empty;
  else
  {
    int index = (int) getDefaultValue();
    if( index >= 0 && index < stringValues.size() )
      return stringValues[index];
    else
      return String::empty;
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
    ParameterObserver *observer = parameterObservers[i];
    if( parameterObservers[i]->wantsAutomationNotification() )
      parameterObservers[i]->parameterChanged(this);
  }
}

void Parameter::callValueChangeCallbacks()
{
  ScopedPointerLock spl(mutex);
  if( valueChangeCallbackDouble != nullptr )
    valueChangeCallbackDouble->call(value);
  if( valueChangeCallbackInt != nullptr )
    valueChangeCallbackInt->call(juce::roundDoubleToInt(value));
  if( valueChangeCallbackBool != nullptr )
    valueChangeCallbackBool->call(value >= 0.5);
}

double Parameter::restrictValueToParameterRange(double valueToRestrict)
{
  ScopedPointerLock spl(mutex);
  if( scaling == BOOLEAN )
    valueToRestrict = (double) (valueToRestrict >= 0.5);
  if( scaling == STRING )
    valueToRestrict = (double) round(valueToRestrict);
  if( valueToRestrict > maxValue )
    return maxValue;
  if( valueToRestrict < minValue )
    return minValue;
  return valueToRestrict;
}

void Parameter::valueSanityCheck()
{
  ScopedPointerLock spl(mutex);
  jassert( !( minValue <= 0.0 && scaling == EXPONENTIAL ) );
  if( ( minValue <= 0.0 && scaling == EXPONENTIAL ) )
    minValue = 0.1;
  jassert( maxValue > minValue );
  if( maxValue <= minValue )
    maxValue = minValue+1.0;
  value                = restrictValueToParameterRange(value);
  defaultValue         = restrictValueToParameterRange(defaultValue);
}

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
  for(int i = 0; i < parameterSetObservers.size(); i++)
    parameterSetObservers[i]->parameterSetChanged(parameterSetHolderThatHasChanged);
}
