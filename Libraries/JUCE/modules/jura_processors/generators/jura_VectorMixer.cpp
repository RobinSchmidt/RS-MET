
//-------------------------------------------------------------------------------------------------
// construction/destruction:

VectorMixerAudioModule::VectorMixerAudioModule(CriticalSection *newPlugInLock, 
  rosic::VectorMixer *newVectorMixerToWrap)
: AudioModule(newPlugInLock)
{
  jassert( newVectorMixerToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  wrappedVectorMixer = newVectorMixerToWrap;
  moduleName = juce::String("VectorMixer");

  // create and initialize the automatable parameters:
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void VectorMixerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedVectorMixer == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case  0: wrappedVectorMixer->setX(value);  break;
  case  1: wrappedVectorMixer->setY(value);  break;
  } // end of switch( parameterIndex )
}

/*
//-------------------------------------------------------------------------------------------------
// state saving and recall:

XmlElement* VectorMixerAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  // store the parameters of the underlying core object:
  if( wrappedVectorMixer != NULL )
    xmlState = vectorPadStateToXml(wrappedVectorMixer, xmlState);

  return xmlState;
}

void VectorMixerAudioModule::setStateFromXml(const XmlElement& xmlState,
                                             const juce::String& stateName, bool markAsClean)
{
  // restore the inherited controller mappings:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  // restore the parameters of the underlying core object:
  if( wrappedVectorMixer != NULL )
    vectorPadStateFromXml(wrappedVectorMixer, xmlState);
}
*/

//-------------------------------------------------------------------------------------------------
// internal functions:

void VectorMixerAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  // this pointer will be used to temporarily store the addresses of the created 
  // AutomatableParameter-objects:
  AutomatableParameter* p;

  // #000:
  p = new AutomatableParameter(lock, "X", -1.0, 1.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #001:
  p = new AutomatableParameter(lock, "Y", -1.0, 1.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================