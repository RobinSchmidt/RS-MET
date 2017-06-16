
//-------------------------------------------------------------------------------------------------
// construction/destruction:

ModuluxuryAudioModule::ModuluxuryAudioModule(CriticalSection *newPlugInLock,
  rosic::Moduluxury *moduluxuryToWrap)
: AudioModule(newPlugInLock)
{
  //jassert(moduluxuryToWrap != NULL); // you must pass a valid rosic-object to the constructor

  if(moduluxuryToWrap != nullptr)
    wrappedModuluxury = moduluxuryToWrap;
  else
  {
    wrappedModuluxury = new rosic::Moduluxury;
    wrappedModuluxuryIsOwned = true;
  }

  moduleName = juce::String("Moduluxury");
  setActiveDirectory(getApplicationDirectory() + juce::String("/ModuluxuryPresets"));
  initializeAutomatableParameters();
}

ModuluxuryAudioModule::~ModuluxuryAudioModule()
{
  if(wrappedModuluxuryIsOwned)
    delete wrappedModuluxury;
}

//-------------------------------------------------------------------------------------------------
// automation:

void ModuluxuryAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedModuluxury == NULL )
    return;

  //double value = parameterThatHasChanged->getValue();

  /*
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedModuluxury->setDryWetRatio(           (float) value );   break;
  case   1: wrappedModuluxury->setLateReverbLevel(       (float) value );   break;
  case   2: wrappedModuluxury->fdn.setReferenceDelayTime(        value );   break;
  case   3: wrappedModuluxury->fdn.setInjectionVector(     (int) value );   break;
  case   4: wrappedModuluxury->fdn.setFeedbackMatrix(      (int) value );   break;
  case   5: wrappedModuluxury->fdn.setOutputVector(        (int) value );   break;
  case   6: wrappedModuluxury->fdn.setAllpassMode(         value>=0.5  );   break;
  } // end of switch( parameterIndex )
  */

  markStateAsDirty();
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void ModuluxuryAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.



  /*
  juce::Array<double> defaultValues;
  Parameter* p;

  p = new Parameter("DryWetRatio", 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new Parameter("LateLevel", -48.0, 6.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new Parameter("ReferenceDelayTime", 5.0, 100.0, 1.0, 50.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  p = new Parameter("InjectionVector", 0.0, 1.0, 1.0, 0.0, Parameter::STRING);
  p->addStringValue(juce::String(T("AllOnes")));
  addObservedParameter(p);

  p = new Parameter("FeedbackMatrix", 0.0, 5.0, 1.0, 4.0, Parameter::STRING);
  p->addStringValue(juce::String(T("Identity")));
  p->addStringValue(juce::String(T("MinusIdentity")));
  p->addStringValue(juce::String(T("SeriesConnection")));
  p->addStringValue(juce::String(T("SeriesWithFeedback")));
  p->addStringValue(juce::String(T("Hadamard")));
  p->addStringValue(juce::String(T("MinusHadamard")));
  addObservedParameter(p);

  p = new Parameter("OutputVector", 0.0, 5.0, 1.0, 4.0, Parameter::STRING);
  p->addStringValue(juce::String(T("AllOnes")));
  p->addStringValue(juce::String(T("Out01")));
  p->addStringValue(juce::String(T("Out02")));
  p->addStringValue(juce::String(T("Out03")));
  p->addStringValue(juce::String(T("Out04")));
  addObservedParameter(p);

  p = new Parameter("AllpassMode", 0.0, 1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  */

  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect
  // the values the automatable parameters:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}
