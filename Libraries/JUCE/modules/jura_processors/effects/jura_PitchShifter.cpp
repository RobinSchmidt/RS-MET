
//-------------------------------------------------------------------------------------------------
// construction/destruction:

PitchShifterAudioModule::PitchShifterAudioModule(CriticalSection *newPlugInLock, 
  rosic::PitchShifterGrainAdaptive *pitchShifterToWrap) : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*lock);
  jassert(pitchShifterToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedPitchShifter = pitchShifterToWrap;
  moduleName = juce::String("PitchShifter");
  setActiveDirectory(getApplicationDirectory() + juce::String("/PitchShifterPresets"));
  createStaticParameters();
}

//-------------------------------------------------------------------------------------------------
// automation and state management:

void PitchShifterAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*lock);

  std::vector<double> defaultValues;
  AutomatableParameter* p;

  p = new AutomatableParameter(lock, "DetuneCoarse", -48.0, 48.0, 0.1, 0.0, 
    Parameter::LINEAR_BIPOLAR, 74);
  defaultValues.clear();
  defaultValues.push_back(2.03910002); // f/f0 =  9/ 8 = 1.125
  defaultValues.push_back(3.15641287); // f/f0         = 1.2
  defaultValues.push_back(3.86313714); // f/f0 =  5/ 4 = 1.25
  //defaultValues.push_back(4.54213948); // f/f0 = 13/10 = 1.3
  defaultValues.push_back(4.98044999); // f/f0 =  4/ 3 = 1.333...
  defaultValues.push_back(7.01955001); // f/f0 =  3/ 2 = 1.5
  defaultValues.push_back(8.84358713); // f/f0 =  5/ 3 = 1.666...
  defaultValues.push_back(9.68825906); // f/f0 =  7/ 4 = 1.75
  defaultValues.push_back(12.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setDetuneCoarse);
    // p->setValueChangeCallback(wrappedPitchShifter, &PitchShifterGrainAdaptive::setDetuneCoarse); does not work because the MS compiler 
    // seems not to be able to infer, that the template argument is "PitchShifterGrainAdaptive" - it ambiguates it with the "PitchShifter" 
    // baseclass - that's why we must pass the template argument here explicitly

  p = new AutomatableParameter(lock, "DetuneFine", -200.0, 200.0, 0.1, 0.0, 
    Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setDetuneFine);

  p = new AutomatableParameter(lock, "GrainLengthInMilliseconds", 1.0, 2000.0, 0.001, 
    18.1818181818181818, Parameter::EXPONENTIAL);
     // a grain length of 18.18... will tune the amplitude modulation artifacts to +-55 Hz with 
     // respect to the carrier frequency
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setGrainLengthInMilliseconds);

  p = new AutomatableParameter(lock, "GrainLengthInPitchCycles", 0.25, 64.0, 0.01, 4.0, 
    Parameter::EXPONENTIAL);
  defaultValues.clear();
  defaultValues.push_back(2.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(4.0);
  defaultValues.push_back(6.0);
  defaultValues.push_back(8.0);
  defaultValues.push_back(12.0);
  defaultValues.push_back(16.0);
  defaultValues.push_back(24.0);
  defaultValues.push_back(32.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setGrainLengthInPitchCycles);

  p = new AutomatableParameter(lock, "GrainLengthInBeats", 0.125, 4.0, 0.001, 0.5, 
    Parameter::EXPONENTIAL);
  defaultValues.clear();
  defaultValues.push_back(0.125);
  defaultValues.push_back(0.25);
  defaultValues.push_back(0.5);
  defaultValues.push_back(0.75);
  defaultValues.push_back(1.0);
  defaultValues.push_back(2.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(4.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setGrainLengthInBeats);

  p = new AutomatableParameter(lock, "GrainLengthUnit", 0.0, 2.0, 1.0, 0.0, Parameter::STRING);
  p->addStringValue(juce::String("ms"));
  p->addStringValue(juce::String("cycles"));
  p->addStringValue(juce::String("beats"));
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setGrainLengthUnit);

  p = new AutomatableParameter(lock, "Feedback", -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setFeedback);

  p = new AutomatableParameter(lock, "DryWet", 0.0, 100.0, 0.1, 100.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setDryWet);

  p = new AutomatableParameter(lock, "AntiAlias", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setAntiAliasing);

  p = new AutomatableParameter(lock, "Reverse", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setReversePlayback);

  p = new AutomatableParameter(lock, "Invert", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setNegativePolarity);

  p = new AutomatableParameter(lock, "FormantPreserve", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setFormantPreserve);

  // make sure that the parameters are initially in sync with the audio engine:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameters[i]->resetToDefaultValue(true, true);
}

XmlElement PitchShifterAudioModule::convertXmlStateIfNecessary(const XmlElement& xmlState)
{
  ScopedLock scopedLock(*lock);

  // retrieve the patch format of the xml-file to enable different interpretations of the patch for 
  // backwards compatibility:
  int xmlPatchFormat = xmlState.getIntAttribute("PatchFormat", 0);
  if( xmlPatchFormat == 0 ) // this is an old preset
  {
    // we had formerly only one "GrainLength" parameter which was either interpreted as being in ms 
    // or in pitch-cycles depending on a boolean flag "GrainLengthAdaption" - now we have 
    // different parameters for the different units and a string-parameter "GrainLengthUnit" to 
    // select which value is to be used
    XmlElement convertedState = xmlState;
    double d = xmlState.getDoubleAttribute("GrainLength", 8.0);
    bool   b = xmlState.getBoolAttribute("GrainLengthAdaption", false);
    if( b == true )
    {
      convertedState.setAttribute("GrainLengthInPitchCycles", d);
      convertedState.setAttribute("GrainLengthUnit", "cycles");
    }
    else
    {
      convertedState.setAttribute("GrainLengthInMilliseconds", d);
      convertedState.setAttribute("GrainLengthUnit", "ms");
    }
    return convertedState;
  }
  else
    return xmlState;
}
