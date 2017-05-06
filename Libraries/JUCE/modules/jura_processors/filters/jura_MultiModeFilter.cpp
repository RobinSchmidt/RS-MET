
//-------------------------------------------------------------------------------------------------
// construction/destruction:

MultiModeFilterAudioModule::MultiModeFilterAudioModule(CriticalSection *newPlugInLock, rosic::MultiModeFilter *newMultiModeFilterToWrap)
 : AudioModule(newPlugInLock)
{
  jassert( newMultiModeFilterToWrap != NULL ); // you must pass a valid rosic-object
  wrappedMultiModeFilter = newMultiModeFilterToWrap;
  moduleName = juce::String("MultiModeFilter");

  // initialize the current directory for preset loading and saving:
  setActiveDirectory(getApplicationDirectory() + juce::String("/MultiModeFilterPresets") );

  // create and initialize the automatable parameters:
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void MultiModeFilterAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedMultiModeFilter == NULL )
    return;

  // depending on the index of the parameterThatHasChanged, we now set up the corresponding 
  // parameter in the underlying rosic object:
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case  0: wrappedMultiModeFilter->setFrequencyNominal(          value);     break;  
  case  1: wrappedMultiModeFilter->setFrequencyByKey  (          value);     break;
  case  2: wrappedMultiModeFilter->setFrequencyByVel  (          value);     break;
  case  3: wrappedMultiModeFilter->setResonance(                 value);     break;
  case  4: wrappedMultiModeFilter->setQ(                         value);     break;
  case  5: wrappedMultiModeFilter->setAllpassFreq(               value);     break; 
  case  6: wrappedMultiModeFilter->setMakeUp(                    value);     break; 
  case  7: wrappedMultiModeFilter->setDrive(                     value);     break;
  //case  8: wrappedMultiModeFilter->setDc(                        value);     break;
  case  9: wrappedMultiModeFilter->setGain(                      value);     break;
  case 10: wrappedMultiModeFilter->setMorph(                     value);     break;
  case 11: wrappedMultiModeFilter->setOrder(               (int) value);     break;

  default:
    {
      // do nothing
    }
  } // end of switch( parameterIndex )

  markStateAsDirty();
}

//-------------------------------------------------------------------------------------------------
// state saving and recall:

// temporary - move into get/setStateFrom/ToXml

int stringToFilterModeIndex(const juce::String &modeString)
{
  if( modeString == juce::String("Bypass") )  
    return rosic::MultiModeFilterParameters::BYPASS;
  else if( modeString == juce::String("Moogish Lowpass") )  
    return rosic::MultiModeFilterParameters::MOOGISH_LOWPASS;
  else if( modeString == juce::String("Lowpass 6 dB/oct") )  
    return rosic::MultiModeFilterParameters::LOWPASS_6;
  else if( modeString == juce::String("Lowpass 12 dB/oct") )  
    return rosic::MultiModeFilterParameters::LOWPASS_RBJ;
  else if( modeString == juce::String("Highpass 6 dB/oct") )  
    return rosic::MultiModeFilterParameters::HIGHPASS_6;
  else if( modeString == juce::String("Highpass 12 dB/oct") )  
    return rosic::MultiModeFilterParameters::HIGHPASS_RBJ;
  else if( modeString == juce::String("Bandpass 2*6 dB/oct") )  
    return rosic::MultiModeFilterParameters::BANDPASS_RBJ;
  else if( modeString == juce::String("Bandstop 2*6 dB/oct") )  
    return rosic::MultiModeFilterParameters::BANDREJECT_RBJ;
  else if( modeString == juce::String("Peak/Dip") )  
    return rosic::MultiModeFilterParameters::PEAK_OR_DIP_RBJ;
  else if( modeString == juce::String("Low Shelv 1st order") )  
    return rosic::MultiModeFilterParameters::LOW_SHELV_1ST;
  else if( modeString == juce::String("Low Shelv 2nd order") )  
    return rosic::MultiModeFilterParameters::LOW_SHELV_RBJ;
  else if( modeString == juce::String("High Shelv 1st order") )  
    return rosic::MultiModeFilterParameters::HIGH_SHELV_1ST;
  else if( modeString == juce::String("High Shelv 2nd order") )  
    return rosic::MultiModeFilterParameters::HIGH_SHELV_RBJ;
  else if( modeString == juce::String("Allpass 1st order") )  
    return rosic::MultiModeFilterParameters::ALLPASS_1ST;
  else if( modeString == juce::String("Allpass 2nd order") )  
    return rosic::MultiModeFilterParameters::ALLPASS_RBJ;

  else if( modeString == juce::String("Morph Low/Band/High") )  
    return rosic::MultiModeFilterParameters::MORPH_LP_BP_HP;

  else if( modeString == juce::String("Morph Low/Peak/High") )  
    return rosic::MultiModeFilterParameters::MORPH_LP_PK_HP;

  // some more else ifs to come...

  else                                                 
    return rosic::MultiModeFilterParameters::BYPASS;
}

const juce::String filterModeIndexToString(int modeIndex)
{
  if(modeIndex == rosic::MultiModeFilterParameters::BYPASS ) 
    return juce::String("Bypass");
  else if( modeIndex == rosic::MultiModeFilterParameters::MOOGISH_LOWPASS ) 
    return juce::String("Moogish Lowpass");

  else if( modeIndex == rosic::MultiModeFilterParameters::LOWPASS_6 ) 
    return juce::String("Lowpass 6 dB/oct");
  else if( modeIndex == rosic::MultiModeFilterParameters::LOWPASS_RBJ ) 
    return juce::String("Lowpass 12 dB/oct");
  else if( modeIndex == rosic::MultiModeFilterParameters::HIGHPASS_6 ) 
    return juce::String("Highpass 6 dB/oct");
  else if( modeIndex == rosic::MultiModeFilterParameters::HIGHPASS_RBJ ) 
    return juce::String("Highpass 12 dB/oct");
  else if( modeIndex == rosic::MultiModeFilterParameters::BANDPASS_RBJ ) 
    return juce::String("Bandpass 2*6 dB/oct");
  else if( modeIndex == rosic::MultiModeFilterParameters::BANDREJECT_RBJ ) 
    return juce::String("Bandstop 2*6 dB/oct");
  else if( modeIndex == rosic::MultiModeFilterParameters::PEAK_OR_DIP_RBJ ) 
    return juce::String("Peak/Dip");
  else if( modeIndex == rosic::MultiModeFilterParameters::LOW_SHELV_1ST ) 
    return juce::String("Low Shelv 1st order");
  else if( modeIndex == rosic::MultiModeFilterParameters::LOW_SHELV_RBJ ) 
    return juce::String("Low Shelv 2nd order");
  else if( modeIndex == rosic::MultiModeFilterParameters::HIGH_SHELV_1ST ) 
    return juce::String("High Shelv 1st order");
  else if( modeIndex == rosic::MultiModeFilterParameters::HIGH_SHELV_RBJ ) 
    return juce::String("High Shelv 2nd order");
  else if( modeIndex == rosic::MultiModeFilterParameters::ALLPASS_1ST ) 
    return juce::String("Allpass 1st order");
  else if( modeIndex == rosic::MultiModeFilterParameters::ALLPASS_RBJ ) 
    return juce::String("Allpass 2nd order");

  else if(modeIndex == rosic::MultiModeFilterParameters::MORPH_LP_BP_HP ) 
    return juce::String("Morph Low/Band/High");

  else if(modeIndex == rosic::MultiModeFilterParameters::MORPH_LP_PK_HP ) 
    return juce::String("Morph Low/Peak/High");

  // some more else ifs to come...

  else                                              
    return juce::String("Bypass");
}

XmlElement* multiModeFilterStateToXml(MultiModeFilter* filter, 
  XmlElement* xmlElementToStartFrom)
{
  // the XmlElement which stores all the releveant state-information:
  XmlElement* xmlState;
  if( xmlElementToStartFrom == NULL )
    xmlState = new XmlElement(juce::String("MultiModeFilterState")); 
  else
    xmlState = xmlElementToStartFrom;

  // convert the filter mode into a string and store the mode-string in the XmlElement:
  juce::String modeString = filterModeIndexToString(filter->getMode());
  xmlState->setAttribute("Mode", modeString);

  // now store the filter parameters - depending on the chosen mode, the set of parameters to be 
  // stored may differ. we distiguish the different cases by means of the modeString:
  if( filter->getMode() == MultiModeFilterParameters::BYPASS )
  {
    // no more parameters to store, when filter is in bypass mode
  }
  else
  {
    // store the cutoff/center frequency parameters and it's key and velocity dependence:
    xmlState->setAttribute("Frequency",      filter->getFrequencyNominal()              );
    xmlState->setAttribute("FrequencyByKey", filter->getFrequencyByKey()                );
    xmlState->setAttribute("FrequencyByVel", filter->getFrequencyByVel()                );
  }

  if( filter->getMode() == MultiModeFilterParameters::MOOGISH_LOWPASS )
  {
    // store the parameters which are specific to the Moogish Lowpass mode:
    xmlState->setAttribute("Resonance",      filter->getResonance()                     );
    xmlState->setAttribute("Drive",          filter->getDrive()                         );
    xmlState->setAttribute("Order",          filter->getOrder()                         );
    xmlState->setAttribute("PreAllpass",     filter->getAllpassFreq()                   );
    xmlState->setAttribute("MakeUp",         filter->getMakeUp()                        );
  }
  if( filter->currentModeSupportsGain() )
    xmlState->setAttribute("Gain",           filter->getGain()                          );
  if( filter->currentModeSupportsQ() )
    xmlState->setAttribute("Q",              filter->getQ()                             );
  if( filter->currentModeSupportsTwoStages() )
    xmlState->setAttribute("TwoStages",      filter->usesTwoStages()                    );

  if( filter->getMode() == MultiModeFilterParameters::MORPH_LP_PK_HP )
    xmlState->setAttribute("Morph",              filter->getMorph()                     );

  return xmlState;
}

bool multiModeFilterStateFromXml(MultiModeFilter* filter, const XmlElement &xmlState)
{
  bool success = true;

  // restore the settings from the XmlElement:
  int modeIndex = stringToFilterModeIndex( xmlState.getStringAttribute("Mode", "Bypass") );
  filter->setMode( modeIndex );
  filter->setFrequencyNominal(xmlState.getDoubleAttribute("Frequency",      1000.0) );
  filter->setFrequencyByKey(  xmlState.getDoubleAttribute("FrequencyByKey",    0.0) );
  filter->setFrequencyByVel(  xmlState.getDoubleAttribute("FrequencyByVel",    0.0) );
  filter->setResonance(       xmlState.getDoubleAttribute("Resonance",         0.0) );
  filter->setQ(               xmlState.getDoubleAttribute("Q",           sqrt(0.5)) );
  filter->setGain(            xmlState.getDoubleAttribute("Gain",              0.0) );
  filter->setMorph(           xmlState.getDoubleAttribute("Morph",             0.0) );
  filter->setDrive(           xmlState.getDoubleAttribute("Drive",             0.0) );
  filter->useTwoStages(       xmlState.getBoolAttribute  ("TwoStages",       false) );
  filter->setOrder(           xmlState.getIntAttribute   ("Order",               4) );
  filter->setAllpassFreq(     xmlState.getDoubleAttribute("PreAllpass",    20000.0) );
  filter->setMakeUp(          xmlState.getDoubleAttribute("MakeUp",          100.0) );

  return success = true;
}

XmlElement* MultiModeFilterAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  // store the parameters of the underlying core object:
  if( wrappedMultiModeFilter != NULL )
    xmlState = multiModeFilterStateToXml(wrappedMultiModeFilter, xmlState);

  return xmlState;
}

void MultiModeFilterAudioModule::setStateFromXml(const XmlElement& xmlState, 
                                                 const juce::String& stateName, bool markAsClean)
{
  // restore the inherited controller mappings:
  AudioModule::setStateFromXml(xmlState, stateName, false);

  // restore the parameters of the underlying core object:
  if( wrappedMultiModeFilter != NULL )
    multiModeFilterStateFromXml(wrappedMultiModeFilter, xmlState);

  if( markAsClean == true )
    markStateAsClean();
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void MultiModeFilterAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  // WARNING: if you change some slider's name here, make sure to also change it in  the editor's
  // calls to getSliderByName - otherwise the editor wil dereference a NULL pointer

  std::vector<double> defaultValues;

  // this pointer will be used to temporarily store the addresses of the created 
  // AutomatableParameter-objects:
  AutomatableParameter* p;

  // #00
  p = new AutomatableParameter(lock, "Frequency", 20.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL, 74);
  defaultValues.push_back(20.0);
  defaultValues.push_back(200.0);
  defaultValues.push_back(2000.0);
  defaultValues.push_back(20000.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
    // todo: define more meaningful default values here - for example tune the frequency to 
    // harmonics (assuming keytrack==100%)

  // #01
  p = new AutomatableParameter(lock, "FrequencyByKey", -200.0, 200.0, 0.1, 0.0, Parameter::LINEAR);
  defaultValues.clear(); 
  defaultValues.push_back(0.0);
  defaultValues.push_back(25.0);
  defaultValues.push_back(50.0);
  defaultValues.push_back(75.0);
  defaultValues.push_back(100.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #02
  p = new AutomatableParameter(lock, "FrequencyByVel", -200.0, 200.0, 0.1, 0.0,  Parameter::LINEAR);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #03
  p = new AutomatableParameter(lock, "Resonance", 0.0, 100.0, 0.1, 10.0, Parameter::LINEAR, 71);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #04
  p = new AutomatableParameter(lock, "Q", 0.5, 50.0, 0.001, sqrt(0.5), Parameter::EXPONENTIAL, 71);
  defaultValues.clear(); 
  defaultValues.push_back(0.5);
  defaultValues.push_back(sqrt(0.5));
  defaultValues.push_back(1.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #05
  p = new AutomatableParameter(lock, "PreAllpass", 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  defaultValues.clear(); 
  defaultValues.push_back(20.0);
  defaultValues.push_back(200.0);
  defaultValues.push_back(2000.0);
  defaultValues.push_back(20000.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #06
  p = new AutomatableParameter(lock, "MakeUp", 0.0, 100.0, 1.0, 0.0,  Parameter::LINEAR);
  defaultValues.clear(); 
  defaultValues.push_back(0.0);
  defaultValues.push_back(25.0);
  defaultValues.push_back(50.0);
  defaultValues.push_back(75.0);
  defaultValues.push_back(100.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #07
  p = new AutomatableParameter(lock, "Drive", -24.0, 24.0, 0.01, 0.0, Parameter::LINEAR);
  defaultValues.clear(); 
  defaultValues.push_back(-24.0);
  defaultValues.push_back(-18.0);
  defaultValues.push_back(-12.0);
  defaultValues.push_back(-6.0);
  defaultValues.push_back(-3.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(6.0);
  defaultValues.push_back(12.0);
  defaultValues.push_back(18.0);
  defaultValues.push_back(24.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #08
  p = new AutomatableParameter(lock, "Dc", -1.0, 1.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #09
  p = new AutomatableParameter(lock, "Gain", -60.0, 30.0, 0.01, 0.0, Parameter::LINEAR);
  defaultValues.clear(); 
  defaultValues.push_back(-24.0);
  defaultValues.push_back(-18.0);
  defaultValues.push_back(-12.0);
  defaultValues.push_back(-6.0);
  defaultValues.push_back(-3.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(6.0);
  defaultValues.push_back(12.0);
  defaultValues.push_back(18.0);
  defaultValues.push_back(24.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #10
  p = new AutomatableParameter(lock, "Morph", -0.99, 0.99, 0.01, 0.0, Parameter::LINEAR);
  defaultValues.clear(); 
  defaultValues.push_back(-0.99);
  defaultValues.push_back(0.5);
  defaultValues.push_back(0.99);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #11
  p = new AutomatableParameter(lock, "Order", 0.0, 4.0, 1.0, 4.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

