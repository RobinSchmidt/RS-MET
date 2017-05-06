//-------------------------------------------------------------------------------------------------
// construction/destruction:

PolyphonicInstrumentAudioModule::PolyphonicInstrumentAudioModule(CriticalSection *newPlugInLock, 
  rosic::PolyphonicInstrument *instrumentToWrap) : AudioModule(newPlugInLock)
{
  jassert(instrumentToWrap != NULL); // you must pass a valid rosic-object to the constructor
  underlyingRosicInstrument = instrumentToWrap;

  // create and initialize the automatable parameters:
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void PolyphonicInstrumentAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedLock scopedLock(*lock);
  AudioModule::parameterRangeChanged(parameterThatHasChanged);

  if( underlyingRosicInstrument == NULL )
    return;

  // find out the index in the vector of the parameter that has been changed:
  int parameterIndex = getIndexOfParameter(parameterThatHasChanged);

  // parameterIndex now contains the index in the array of the parameter that has changed now set 
  // up the signal processing:

  // \todo replace this dispatching switch statement by the new callback-system

  double value = parameterThatHasChanged->getValue();
  switch( parameterIndex )
  {
  case   0: underlyingRosicInstrument->setMasterLevel(        value);     break;
  case   1: underlyingRosicInstrument->setVoiceLevelByKey(    value);     break;
  case   2: underlyingRosicInstrument->setVoiceLevelByVel(    value);     break;
  case   3: underlyingRosicInstrument->setMasterLevelByVoices(value);     break;
  case   4: underlyingRosicInstrument->setMidSideRatio(       value);     break;
  case   5: underlyingRosicInstrument->setGlideMode(   0.0 != value);     break;
  case   6: underlyingRosicInstrument->setGlideTime(          value);     break;
  case   7: underlyingRosicInstrument->setMasterTuneA4(       value);     break;
  case   8: underlyingRosicInstrument->setPitchWheelRange(    value);     break;
  default:
    {
      // do nothing
    }

  } // end of switch( parameterIndex )

}

void PolyphonicInstrumentAudioModule::setStateFromXml(const XmlElement& xmlState, 
  const juce::String& stateName, bool markAsClean) 
{
  underlyingRosicInstrument->allNotesOff();
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
}

/*
XmlElement* PolyphonicInstrumentAudioModule::getStateAsXml(XmlElement* xmlElementToStartFrom)
{
  if( underlyingRosicInstrument == NULL )
    return NULL;

  XmlElement* xmlState;
  if( xmlElementToStartFrom == NULL )
    xmlState = new XmlElement(juce::String(T("InstrumentState"))); 
  else
    xmlState = xmlElementToStartFrom;

  // store the setting of the inherited AudioModule object:
  xmlState = AudioModule::getStateAsXml(xmlState);

  return xmlState;
}

void PolyphonicInstrumentAudioModule::setStateFromXml(const XmlElement &xmlState)
{
  if( underlyingRosicInstrument == NULL )
    return;

  // we want to preserve the clean state but when we set automatable parameters it will be set to 
  // dirty - so we remember if it was clean and restore it after setting the parameters:
  bool presetIsClean = !underlyingRosicInstrument->isPresetDirty();

  // restore the settings of the inherited AudioModule object:
  AudioModule::setStateFromXml(xmlState);

  if( presetIsClean )
    underlyingRosicInstrument->markPresetAsClean();
}
*/

//-------------------------------------------------------------------------------------------------
// internal functions:

void PolyphonicInstrumentAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particular parameter has changed.

  // this pointer will be used to temporarily store the addresses of the created 
  // Parameter-objects:
  AutomatableParameter* p;

  // #000:
  p = new AutomatableParameter(lock, "MasterLevel", -36.0, 12.0, 0.1, 0.0, Parameter::LINEAR, 7);
  addObservedParameter(p);

  // #001:
  p = new AutomatableParameter(lock, "VoiceLevelByKey", -24.0, 24.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #002:
  p = new AutomatableParameter(lock, "VoiceLevelByVel", 0.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #003:
  p = new AutomatableParameter(lock, "MasterLevelByVoices", 0.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #004:
  p = new AutomatableParameter(lock, "MidSideRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p);

  // #005:
  p = new AutomatableParameter(lock, "GlideSwitch", 0.0, 1.0, 0.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #006:
  p = new AutomatableParameter(lock, "GlideTime", 5.0, 2000.0, 0.0, 50.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #007:
  p = new AutomatableParameter(lock, "MasterTuneA4", 220.0, 880.0, 0.01, 440.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #008:
  p = new AutomatableParameter(lock, "PitchWheelRange", 0.0, 24.0, 0.1, 12.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}
