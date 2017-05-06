

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

EchoLabDelayLineAudioModule::EchoLabDelayLineAudioModule(CriticalSection *newPlugInLock, rosic::EchoLabDelayLine *echoLabDelayLineToWrap)
: AudioModule(newPlugInLock)
{
  jassert( echoLabDelayLineToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  wrappedEchoLabDelayLine = echoLabDelayLineToWrap;
  moduleName = juce::String("DelayLine");
  setActiveDirectory(getApplicationDirectory() + juce::String("/EchoLabDelayLinePresets") );
  initializeAutomatableParameters();

  inputEqualizerModule = new EqualizerAudioModule(lock, &echoLabDelayLineToWrap->inputEqualizer);
  inputEqualizerModule->setModuleName(juce::String("InputFilter"));
  inputEqualizerModule->setActiveDirectory(getApplicationDirectory() + juce::String("/EchoLabPresets/InputFilterPresets") );
  addChildAudioModule(inputEqualizerModule);

  feedbackEqualizerModule = new EqualizerAudioModule(lock, &echoLabDelayLineToWrap->feedbackEqualizer);
  feedbackEqualizerModule->setModuleName(juce::String("FeedbackFilter"));
  feedbackEqualizerModule->setActiveDirectory(getApplicationDirectory() + juce::String("/EchoLabPresets/FeedbackFilterPresets") );
  addChildAudioModule(feedbackEqualizerModule);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void EchoLabDelayLineAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedEchoLabDelayLine == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case 0: wrappedEchoLabDelayLine->setDelayTime(        value         );   break;
  case 1: wrappedEchoLabDelayLine->setGlobalGainFactor( value         );   break;
  case 2: wrappedEchoLabDelayLine->setFeedbackInPercent(value         );   break;
  case 3: wrappedEchoLabDelayLine->setPan(              value         );   break;
  case 4: wrappedEchoLabDelayLine->setPingPongMode(     value >= 0.5  );   break;
  case 5: wrappedEchoLabDelayLine->setMute(             value >= 0.5  );   break;
  } 

  //sendChangeMessage();
  markStateAsDirty();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:







//-----------------------------------------------------------------------------------------------------------------------------------------
// audio processing:

/*
void EchoLabDelayLineAudioModule::getSampleFrameStereo(double* inOutL, double* inOutR)
{ 
  *inOutL = 0.0;
  *inOutR = 0.0;
  jassertfalse; // this function is not supposed to be used
}

void EchoLabDelayLineAudioModule::processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages) 
{ 
  buffer.clear();
  jassertfalse; // this function is not supposed to be used
} 
*/

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void EchoLabDelayLineAudioModule::initializeAutomatableParameters()
{
  std::vector<double> defaultValues;
  AutomatableParameter* p;

  // #00:
  p = new AutomatableParameter(lock, "DelayTime", 0.0001, 4.25, 0.0001, 1.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(0.125);
  defaultValues.push_back(0.25);
  defaultValues.push_back(0.5);
  defaultValues.push_back(1.0);
  defaultValues.push_back(2.0);
  defaultValues.push_back(4.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  addObservedParameter( new AutomatableParameter(lock, "Amplitude",  -1.0,  1.0,  0.001, 0.5, Parameter::LINEAR) );  // #01
  addObservedParameter( new AutomatableParameter(lock, "Feedback",  -99.0, 99.0,  0.1,   0.0, Parameter::LINEAR) );  // #02
  addObservedParameter( new AutomatableParameter(lock, "Pan",        -1.0,  1.0,  0.01,  0.0, Parameter::LINEAR) );  // #03
  addObservedParameter( new AutomatableParameter(lock, "PingPong",    0.0,  1.0,  1.0,   0.0, Parameter::BOOLEAN));  // #04
  addObservedParameter( new AutomatableParameter(lock, "Mute",        0.0,  1.0,  1.0,   0.0, Parameter::BOOLEAN));  // #05

  //addObservedParameter( new AutomatableParameter(plugInLock, "Solo",        0.0,  1.0,  1.0,   0.0, Parameter::BOOLEAN));  // #05
    // nah - "Solo" is not a parameter - it's a GUI feature - when slo is on, alway the currently selected delayline will play solo

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

/*
void EchoLabDelayLineAudioModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, bool markAsClean)
{
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);


  int dummy = 0;

}
*/

void EchoLabDelayLineAudioModule::reset()
{
  //wrappedEchoLabDelayLine->reset();  // maybe we should indeed implement this rest-function...
}

