
//-------------------------------------------------------------------------------------------------
// construction/destruction:

AutomatableParameter::AutomatableParameter(CriticalSection *criticalSectionToUse, 
  const String& newName, double newMinValue, double newMaxValue, double newInterval, 
  double newDefaultValue, int newScaling, int newDefaultMidiController)                     
: Parameter(criticalSectionToUse, newName, newMinValue, newMaxValue, newInterval, 
  newDefaultValue, newScaling) 
{
  lowerAutomationLimit   = newMinValue;
  upperAutomationLimit   = newMaxValue;
  defaultMidiController  = newDefaultMidiController;
  assignedMidiController = newDefaultMidiController;
  learnMode              = false;
  updateAutomationFromParameter();
}

AutomatableParameter::~AutomatableParameter()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void AutomatableParameter::setAutomationValue(double newAutomationValue, bool notifyListeners, 
  bool callCallbacks)
{
  ScopedPointerLock spl(mutex);
  if( newAutomationValue < 0.0 || newAutomationValue > 1.0 )
    jassertfalse; // automation has to be in the range 0...1
  automationValue = newAutomationValue;
  updateParameterFromAutomation();
  setValue(value, notifyListeners, callCallbacks);
}

void AutomatableParameter::setLowerAutomationLimit(double newLimit)
{
  ScopedPointerLock spl(mutex);
  lowerAutomationLimit = restrictValueToParameterRange(newLimit);
  updateAutomationFromParameter();
}

void AutomatableParameter::setUpperAutomationLimit(double newLimit)
{
  ScopedPointerLock spl(mutex);
  upperAutomationLimit = restrictValueToParameterRange(newLimit);
  updateAutomationFromParameter();
}

void AutomatableParameter::setDefaultMidiController(int newDefaultControllerNumber)
{
  ScopedPointerLock spl(mutex);
  if( newDefaultControllerNumber < -1 || newDefaultControllerNumber > 127 )
  {
    jassertfalse; // invalid controller number
    return;
  }
  defaultMidiController = newDefaultControllerNumber;
}

void AutomatableParameter::assignMidiController(int controllerNumberToAssign)
{
  ScopedPointerLock spl(mutex);
  if( controllerNumberToAssign < -1 || controllerNumberToAssign > 127 )
  {
    jassertfalse; // invalid controller number
    return;
  }
  assignedMidiController = controllerNumberToAssign;
}

void AutomatableParameter::revertToDefaultMidiController()
{
  ScopedPointerLock spl(mutex);
  assignedMidiController = defaultMidiController;
}

void AutomatableParameter::setMidiController(int controllerNumber, float controllerValue)
{
  ScopedPointerLock spl(mutex);
  if( controllerNumber == assignedMidiController )
    setAutomationValue( (double) controllerValue / 127.0, true, true );
  else if( learnMode == true )
  {
    assignedMidiController = controllerNumber;
    learnMode              = false;
    setAutomationValue( (double) controllerValue / 127.0, true, true );
  }
}

void AutomatableParameter::switchIntoMidiLearnMode(bool shouldLearn)
{
  ScopedPointerLock spl(mutex);
  learnMode = shouldLearn;
}

void AutomatableParameter::revertToDefaults(bool resetValueAlso, bool sendNotification, 
  bool callCallbacks)
{
  ScopedPointerLock spl(mutex);
  assignedMidiController = defaultMidiController;
  lowerAutomationLimit   = minValue;
  upperAutomationLimit   = maxValue;
  if( resetValueAlso == true )
    setValue(defaultValue, sendNotification, callCallbacks);
}

//-------------------------------------------------------------------------------------------------
// inquiry:

bool AutomatableParameter::isInDefaultState() const
{
  ScopedPointerLock spl(mutex);
  if( assignedMidiController != defaultMidiController )
    return false;
  if( minValue != lowerAutomationLimit )
    return false;
  if( maxValue != upperAutomationLimit )
    return false;
  return true;
}

//-------------------------------------------------------------------------------------------------
// others:

double AutomatableParameter::restrictValueToAutomatableRange(double valueToRestrict)
{
  ScopedPointerLock spl(mutex);
  if( valueToRestrict > upperAutomationLimit )
    return upperAutomationLimit;
  if( valueToRestrict < lowerAutomationLimit )
    return lowerAutomationLimit;
  return valueToRestrict;
}

void AutomatableParameter::valueSanityCheck()
{
  ScopedPointerLock spl(mutex);
  Parameter::valueSanityCheck();
  lowerAutomationLimit = restrictValueToParameterRange(lowerAutomationLimit);
  upperAutomationLimit = restrictValueToParameterRange(upperAutomationLimit);
}

void AutomatableParameter::updateParameterFromAutomation()
{
  ScopedPointerLock spl(mutex);
  switch( scaling )
  {
  case LINEAR:
    {
      value = linToLin(
        automationValue, 0.0, 1.0, lowerAutomationLimit, upperAutomationLimit);
    }
    break;

  case EXPONENTIAL:
    {
      if( lowerAutomationLimit <= 0.0 )
        jassertfalse; // exponential scaling requires strictly positive minimum value
      value = linToExp(
        automationValue, 0.0, 1.0, lowerAutomationLimit, upperAutomationLimit);
    }
    break;

  default: // linear mapping by default
    {
      value = linToLin(
        automationValue, 0.0, 1.0, lowerAutomationLimit, upperAutomationLimit);
    }
  }
}

void AutomatableParameter::updateAutomationFromParameter()
{
  ScopedPointerLock spl(mutex);
  // use a value restricted to the automatable range instead of the full range value itself in 
  // order to have the automationValue in the range 0...1:
  double tmpValue = value;
  if( value < lowerAutomationLimit )
    tmpValue = lowerAutomationLimit;
  if( value > upperAutomationLimit )
    tmpValue = upperAutomationLimit;

  switch( scaling )
  {
  case LINEAR:
    {
      automationValue = linToLin(
        tmpValue, lowerAutomationLimit, upperAutomationLimit, 0.0, 1.0);
    }
    break;
  case EXPONENTIAL:
    {
      if( lowerAutomationLimit <= 0.0 )
        jassertfalse; // exponential scaling requires strictly positive minimum value

      automationValue = expToLin(
        tmpValue, lowerAutomationLimit, upperAutomationLimit, 0.0, 1.0);
    }
    break;
  default: // linear mapping by default
    {
      automationValue = linToLin(
        tmpValue, lowerAutomationLimit, upperAutomationLimit, 0.0, 1.0);
    }
  }
}

//=================================================================================================

MetaParameter::MetaParameter()
{

}

void MetaParameter::addParameter(AutomatableParameter* p)
{
  p->setAutomationValue(value);
  p->registerParameterObserver(this);
  appendIfNotAlreadyThere(params, p);
}

void MetaParameter::removeParameter(AutomatableParameter* p)
{
  p->deRegisterParameterObserver(this);
  removeFirstOccurrence(params, p);
}

void MetaParameter::setAutomationValue(double v)
{  
  value = v;
  localAutomationSwitch = false; // so we don't call ourselves recursively
  for(int i = 0; i < size(params); i++)  
    params[i]->setAutomationValue(value, true, true);
  localAutomationSwitch = true;
}

void MetaParameter::parameterChanged(Parameter* p)
{
  AutomatableParameter* ap = dynamic_cast<AutomatableParameter*>(p);
  value = ap->getAutomationValue();
  localAutomationSwitch = false; // so we don't call ourselves recursively
  for(int i = 0; i < size(params); i++)
  {
    if(params[i] != ap)
      params[i]->setAutomationValue(value, true, true);
  }
  localAutomationSwitch = true;
}

void MetaParameter::parameterIsGoingToBeDeleted(Parameter* p)
{
  AutomatableParameter* ap = dynamic_cast<AutomatableParameter*>(p);
  removeParameter(ap);
}