//-------------------------------------------------------------------------------------------------
// construction/destruction:

StereoDelayAudioModule::StereoDelayAudioModule(CriticalSection *newPlugInLock, 
  rosic::StereoDelay *stereoDelayToWrap)
: AudioModule(newPlugInLock)
{
  jassert(stereoDelayToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedStereoDelay = stereoDelayToWrap;
  moduleName = juce::String("StereoDelay");
  setActiveDirectory(getApplicationDirectory() + juce::String(("/StereoDelayPresets")) );
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// state management:

XmlElement* StereoDelayAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  if( wrappedStereoDelay == NULL )
    return xmlState;

  // add attributes for the non-automatable parameters (the automatable ones are already taken care
  // of by AudioModule::getStateAsXml()):
  xmlState->setAttribute(("DelayInBeatsL"), wrappedStereoDelay->getDelayInBeats(StereoDelay::LEFT)  );
  xmlState->setAttribute(("DelayInBeatsR"), wrappedStereoDelay->getDelayInBeats(StereoDelay::RIGHT) );

  return xmlState;
}

void StereoDelayAudioModule::setStateFromXml(const XmlElement& xmlState,
                                             const juce::String& stateName, bool markAsClean)
{
  // restore the inherited controller mappings:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  if( wrappedStereoDelay == NULL )
    return;

  // restore the values of the non-automatable parameters (the automatable ones are already taken 
  // care of by automatableModuleStateFromXml():
  wrappedStereoDelay->setDelayInBeats(
    xmlState.getDoubleAttribute(("DelayInBeatsL"),   1.0), StereoDelay::LEFT  );
  wrappedStereoDelay->setDelayInBeats(
    xmlState.getDoubleAttribute(("DelayInBeatsR"),   1.0), StereoDelay::RIGHT );
}

//-------------------------------------------------------------------------------------------------
// automation:

void StereoDelayAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedStereoDelay == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedStereoDelay->setDelayScale(                value, StereoDelay::LEFT);                       break;
  case   1: wrappedStereoDelay->setDelayScale(                value, StereoDelay::RIGHT);                      break;
  case   2: wrappedStereoDelay->setInjection(                 value, StereoDelay::LEFT,  StereoDelay::LEFT);   break;
  case   3: wrappedStereoDelay->setInjection(                 value, StereoDelay::RIGHT, StereoDelay::LEFT);   break;
  case   4: wrappedStereoDelay->setInjection(                 value, StereoDelay::LEFT,  StereoDelay::RIGHT);  break;
  case   5: wrappedStereoDelay->setInjection(                 value, StereoDelay::RIGHT, StereoDelay::RIGHT);  break;
  case   6: wrappedStereoDelay->setDiffusorTimeInMilliseconds(value, StereoDelay::LEFT);                       break;
  case   7: wrappedStereoDelay->setDiffusorTimeInMilliseconds(value, StereoDelay::RIGHT);                      break;
  case   8: wrappedStereoDelay->setDiffusorAmount(            value, StereoDelay::LEFT);                       break;
  case   9: wrappedStereoDelay->setDiffusorAmount(            value, StereoDelay::RIGHT);                      break;
  case  10: wrappedStereoDelay->setLowpassCutoff(             value, StereoDelay::LEFT);                       break;
  case  11: wrappedStereoDelay->setLowpassCutoff(             value, StereoDelay::RIGHT);                      break;
  case  12: wrappedStereoDelay->setHighpassCutoff(            value, StereoDelay::LEFT);                       break;
  case  13: wrappedStereoDelay->setHighpassCutoff(            value, StereoDelay::RIGHT);                      break;
  case  14: wrappedStereoDelay->setCutoffScale(               value);                                          break;
  case  15: wrappedStereoDelay->setFeedback(                  value, StereoDelay::LEFT,  StereoDelay::LEFT);   break;
  case  16: wrappedStereoDelay->setFeedback(                  value, StereoDelay::RIGHT, StereoDelay::RIGHT);  break;
  case  17: wrappedStereoDelay->setFeedback(                  value, StereoDelay::LEFT,  StereoDelay::RIGHT);  break;
  case  18: wrappedStereoDelay->setFeedback(                  value, StereoDelay::RIGHT, StereoDelay::LEFT);   break;
  case  19: wrappedStereoDelay->setOutputMix(                 value, StereoDelay::LEFT,  StereoDelay::LEFT);   break;
  case  20: wrappedStereoDelay->setOutputMix(                 value, StereoDelay::LEFT,  StereoDelay::RIGHT);  break;
  case  21: wrappedStereoDelay->setOutputMix(                 value, StereoDelay::RIGHT, StereoDelay::LEFT);   break;
  case  22: wrappedStereoDelay->setOutputMix(                 value, StereoDelay::RIGHT, StereoDelay::RIGHT);  break;
  case  23: wrappedStereoDelay->setWetDelayInBeats(           value, StereoDelay::LEFT);                       break;
  case  24: wrappedStereoDelay->setWetDelayInBeats(           value, StereoDelay::RIGHT);                      break;
  case  25: wrappedStereoDelay->setDryWet(                    value);                                          break;
  } // end of switch( parameterIndex )
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void StereoDelayAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  // this pointer will be used to temporarily store the addresses of the created 
  // Parameter-objects:
  AutomatableParameter* p;

  // #000:
  p = new AutomatableParameter(lock, "DelayScaleL", 0.5, 2.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #001:
  p = new AutomatableParameter(lock, "DelayScaleR", 0.5, 2.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #002:
  p = new AutomatableParameter(lock, "LeftInToLeftDelay", -100.0, 100.0, 0.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #003:
  p = new AutomatableParameter(lock, "RightInToLeftDelay", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #004:
  p = new AutomatableParameter(lock, "LeftInToRightDelay", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #005:
  p = new AutomatableParameter(lock, "RightInToRightDelay", -100.0, 100.0, 0.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #006:
  p = new AutomatableParameter(lock, "DiffusorTimeL", 1.0, 50.0, 0.0, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #007:
  p = new AutomatableParameter(lock, "DiffusorTimeR", 1.0, 50.0, 0.0, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #008:
  p = new AutomatableParameter(lock, "DiffusorAmountL", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #009:
  p = new AutomatableParameter(lock, "DiffusorAmountR", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #010:
  p = new AutomatableParameter(lock, "LowpassCutoffL", 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #011:
  p = new AutomatableParameter(lock, "LowpassCutoffR", 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #012:
  p = new AutomatableParameter(lock, "HighpassCutoffL", 20.0, 20000.0, 0.0, 20.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #013:
  p = new AutomatableParameter(lock, "HighpassCutoffR", 20.0, 20000.0, 0.0, 20.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #014:
  p = new AutomatableParameter(lock, "CutoffScale", 0.0625, 16.0, 0.0, 1.0, Parameter::EXPONENTIAL, 74);
  addObservedParameter(p);

  // #015:
  p = new AutomatableParameter(lock, "FeedbackLeftToLeft", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #016:
  p = new AutomatableParameter(lock, "FeedbackRightToRight", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #017:
  p = new AutomatableParameter(lock, "FeedbackLeftToRight", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #018:
  p = new AutomatableParameter(lock, "FeedbackRightToLeft", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #019:
  p = new AutomatableParameter(lock, "LeftDelayToLeftOut", -100.0, 100.0, 0.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #020:
  p = new AutomatableParameter(lock, "LeftDelayToRightOut", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #021:
  p = new AutomatableParameter(lock, "RightDelayToLeftOut", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #022:
  p = new AutomatableParameter(lock, "RightDelayToRightOut", -100.0, 100.0, 0.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #023:
  p = new AutomatableParameter(lock, "WetDelayInBeatsL", 0.0, 4.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #025:
  p = new AutomatableParameter(lock, "WetDelayInBeatsR", 0.0, 4.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #025:
  p = new AutomatableParameter(lock, "DryWet", 0.0, 100.0, 0.0, 50.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect 
  // the values the automatable parameters:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================