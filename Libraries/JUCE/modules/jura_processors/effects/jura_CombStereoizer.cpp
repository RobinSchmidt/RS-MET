#include "rosof_CombStereoizerAudioModule.h"
using namespace rosof;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

CombStereoizerAudioModule::CombStereoizerAudioModule(CriticalSection *newPlugInLock, rosic::CombStereoizer *stereoizerToWrap)
 : AudioModule(newPlugInLock)
{
  jassert(stereoizerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedCombStereoizer = stereoizerToWrap;
  moduleName        = juce::String(T("CombStereoizer"));
  setActiveDirectory(getApplicationDirectory() + juce::String(T("/CombStereoizerPresets")) );
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// state management:

XmlElement* CombStereoizerAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  if( wrappedCombStereoizer == NULL )
    return xmlState;

  // add attributes for the non-automatable parameters (the automatable ones are already taken care
  // of by AudioModule::getStateAsXml()):
  xmlState->setAttribute(T("ChannelSwitch"), wrappedCombStereoizer->isWetSignalChannelSwitched() );

  return xmlState;
}

void CombStereoizerAudioModule::setStateFromXml(const XmlElement& xmlState,
                                            const juce::String& stateName, bool markAsClean)
{
  // restore the inherited controller mappings:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  if( wrappedCombStereoizer == NULL )
    return;

  // restore the values of the non-automatable parameters (the automatable ones are already taken 
  // care of by automatableModuleStateFromXml():
  wrappedCombStereoizer->setSwitchWetLeftForRight(xmlState.getBoolAttribute(T("ChannelSwitch"), false));
}

//-------------------------------------------------------------------------------------------------
// automation:

void CombStereoizerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedCombStereoizer == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedCombStereoizer->setDryWetRatio(                  (float) value);   break;
  case   1: wrappedCombStereoizer->delayLine.setDelayInMilliseconds((float) value);   break;
  case   2: wrappedCombStereoizer->wetFilter.setLowpassCutoff(      (float) value);   break;
  case   3: wrappedCombStereoizer->wetFilter.setHighpassCutoff(     (float) value);   break;
  case   4: wrappedCombStereoizer->wetFilter.setAllpassFrequency(   (float) value);   break;
  } // end of switch( parameterIndex )
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void CombStereoizerAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  // this pointer will be used to temporarily store the addresses of the created 
  // Parameter-objects:
  AutomatableParameter* p;

  // #000:
  p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 0.5);
  addObservedParameter(p);

  // #001:
  p = new AutomatableParameter(plugInLock, "DelayInMilliseconds", 1.0, 100.0, 0.0, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #002:
  p = new AutomatableParameter(plugInLock, "WetLowpass", 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #003:
  p = new AutomatableParameter(plugInLock, "WetHighpass", 20.0, 20000.0, 0.0, 20.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #004:
  p = new AutomatableParameter(plugInLock, "WetAllpass", 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);


  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect 
  // the values the automatable parameters:
  for(int i=0; i < (int) observedParameters.size(); i++ )
    parameterChanged(observedParameters[i]);
}
