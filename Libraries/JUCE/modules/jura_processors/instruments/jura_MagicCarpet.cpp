//=================================================================================================
// class DelayPhaserAudioModule

DelayPhaserAudioModule::DelayPhaserAudioModule(CriticalSection *newPlugInLock, 
  rosic::DelayPhaser *newDelayPhaserToWrap) : AudioModule(newPlugInLock)
{
  jassert( newDelayPhaserToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedDelayPhaser = newDelayPhaserToWrap;
  moduleName = juce::String("DelayPhaser");
  setActiveDirectory(getApplicationDirectory() + juce::File::separatorString 
    + juce::String(("DelayPhaserPresets")) );
  initializeAutomatableParameters();

  phaser1Module = new PhaserAudioModule(lock, &wrappedDelayPhaser->phaser1);
  phaser1Module->setModuleName(juce::String("Phaser1"));
  addChildAudioModule(phaser1Module);

  delayModule = new PingPongEchoAudioModule(lock, &wrappedDelayPhaser->delay);
  delayModule->setModuleName(juce::String(("Delay")));
  addChildAudioModule(delayModule);

  phaser2Module = new PhaserAudioModule(lock, &wrappedDelayPhaser->phaser2);
  phaser2Module->setModuleName(juce::String(("Phaser2")));
  addChildAudioModule(phaser2Module);
}

void DelayPhaserAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedDelayPhaser == NULL )
    return;
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  //case  0: wrappedDelayPhaser->setDryWetRatio(      value); break;
  case  1: wrappedDelayPhaser->setFeedback1(   0.01*value); break;
  case  2: wrappedDelayPhaser->setFeedback2(   0.01*value); break;
  case  3: wrappedDelayPhaser->setFeedback3(   0.01*value); break;
  } 
}

void DelayPhaserAudioModule::initializeAutomatableParameters()
{
  AutomatableParameter* p;

  p = new AutomatableParameter(lock, "DryWetRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(lock, "FeedbackPhaser1Delay", -99.0, 99.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(lock, "FeedbackDelayPhaser2", -99.0, 99.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(lock, "FeedbackGlobal",       -99.0, 99.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================
// class MagicCarpetAudioModule

//-------------------------------------------------------------------------------------------------
// construction/destruction:

MagicCarpetAudioModule::MagicCarpetAudioModule(CriticalSection *newPlugInLock, rosic::MagicCarpet *magicCarpetToWrap)
: PolyphonicInstrumentAudioModule(newPlugInLock, magicCarpetToWrap)
{
  jassert(magicCarpetToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedMagicCarpet        = magicCarpetToWrap;
  underlyingRosicInstrument = magicCarpetToWrap;
  moduleName                = juce::String(("MagicCarpet"));

  // initialize the current directory for preset loading and saving:
  setActiveDirectory(getApplicationDirectory() + juce::String(("/MagicCarpetPresets")) );

  oscSectionModule = new VectorSamplePlayerAudioModule(lock, &wrappedMagicCarpet->voiceArray[0].oscSection);
  oscSectionModule->setModuleName(juce::String(("Oscillators")));
  addChildAudioModule(oscSectionModule);

  filterModule = new FourPoleFilterAudioModule(lock, &wrappedMagicCarpet->voiceArray[0].filter);
  filterModule->setModuleName(juce::String(("Filter")));
  addChildAudioModule(filterModule);

  filterEnvModule = new BreakpointModulatorAudioModule(lock, &wrappedMagicCarpet->voiceArray[0].filterEnv);
  filterEnvModule->setModuleName(juce::String(("FilterEnvelope")));
  addChildAudioModule(filterEnvModule);

  ampEnvModule = new BreakpointModulatorAudioModule(lock, &wrappedMagicCarpet->voiceArray[0].ampEnv);
  ampEnvModule->setModuleName(juce::String(("AmpEnvelope")));
  addChildAudioModule(ampEnvModule);

  equalizerModule = new EqualizerAudioModule(lock, &wrappedMagicCarpet->equalizer);
  equalizerModule->setModuleName(juce::String(("Equalizer")));
  addChildAudioModule(equalizerModule);

  delayPhaserModule = new DelayPhaserAudioModule(lock, &wrappedMagicCarpet->delayPhaser);
  delayPhaserModule->setModuleName(juce::String(("DelayPhaser")));
  addChildAudioModule(delayPhaserModule);
}