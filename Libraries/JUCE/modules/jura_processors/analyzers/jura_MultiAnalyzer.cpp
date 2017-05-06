

//=================================================================================================
// class OscilloscopeAudioModule:

OscilloscopeAudioModule::OscilloscopeAudioModule(CriticalSection *newPlugInLock, 
  rosic::SyncedWaveformDisplayBuffer *displayBufferToUse)
 : AudioModule(newPlugInLock)
{
  jassert(displayBufferToUse != NULL); // you must pass a valid rosic-object to the constructor
  waveformDisplayBuffer = displayBufferToUse;
  moduleName          = juce::String(T("Oscilloscope"));
  setActiveDirectory(getApplicationDirectory() + juce::String(T("/OscilloscopePresets")) );
  initializeAutomatableParameters();
  displayIsFrozen = false;
}

void OscilloscopeAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( waveformDisplayBuffer == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case 0: waveformDisplayBuffer->setSyncMode(         (int) value  );  break;
  //case 1: waveformDisplayBuffer->setMidSideMode(      value != 0.0 );  break;
  case 2: waveformDisplayBuffer->setTimeWindowLength( value        );  break;
  } 

  markStateAsDirty(); // this feature is de-activated because the state will be marked as dirty immediately after preset-load which 
                        // renders it meaningless - fix this!
}

void OscilloscopeAudioModule::initializeAutomatableParameters()
{
  Parameter *p = new Parameter(plugInLock, "SyncMode", 0.0, 1.0, 1.0, 1.0, Parameter::STRING);
  p->addStringValue(juce::String(T("Free Running")));
  p->addStringValue(juce::String(T("Lowpass Zeros")));
  //p->setValue(1.0);
  addObservedParameter(p);
  addObservedParameter( new Parameter(plugInLock, "MidSideMode",       0.0,   1.0, 1.0,    0.0,  Parameter::BOOLEAN) );
  addObservedParameter( new Parameter(plugInLock, "TimeWindowLength",  0.001, 1.5, 0.001,  0.1,  Parameter::LINEAR)  );
  addObservedParameter( new Parameter(plugInLock, "MinAmplitude",     -2.0,  +2.0, 0.1,   -1.5,  Parameter::LINEAR)  );
  addObservedParameter( new Parameter(plugInLock, "MaxAmplitude",     -2.0,  +2.0, 0.1,    1.5,  Parameter::LINEAR)  );
  addObservedParameter( new Parameter(plugInLock, "FrameRate",        10.0,  50.0, 1.0,   15.0,  Parameter::LINEAR)  );
  addObservedParameter( new Parameter(plugInLock, "Freeze",            0.0,   1.0, 1.0,    1.0,  Parameter::BOOLEAN) );

  for(int i=0; i < (int) observedParameters.size(); i++ )
    parameterChanged(observedParameters[i]);
}

//=================================================================================================
// class SpectrumAnalyzerAudioModule:

SpectrumAnalyzerAudioModule::SpectrumAnalyzerAudioModule(CriticalSection *newPlugInLock, rosic::SpectrumAnalyzer *spectrumAnalyzerToWrap)
 : AudioModule(newPlugInLock)
{
  jassert(spectrumAnalyzerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedSpectrumAnalyzer = spectrumAnalyzerToWrap;
  moduleName              = juce::String(T("SpectrumAnalyzer"));
  setActiveDirectory(getApplicationDirectory() + juce::String(T("/SpectrumAnalyzerPresets")) );
  initializeAutomatableParameters();
  displayIsFrozen  =  false;
}

void SpectrumAnalyzerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedSpectrumAnalyzer == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case 0: wrappedSpectrumAnalyzer->setBlockSize(   (int) pow(2.0, value+8)     );  break;
  case 1: wrappedSpectrumAnalyzer->setMidSideMode(                value != 0.0 );  break;

  } 

  markStateAsDirty(); // this feature is de-activated because the state will be marked as dirty immediately after preset-load which 
                        // renders it meaningless - fix this!
}

void SpectrumAnalyzerAudioModule::initializeAutomatableParameters()
{
  addObservedParameter( new ParameterPowersOfTwo(plugInLock, "FFTSize", 256.0, 32768.0, 0.0, 1024.0) );

  addObservedParameter( new Parameter(plugInLock, "MidSideMode",        0.0,       1.0, 1.0,     0.0,   Parameter::BOOLEAN)      );
  addObservedParameter( new Parameter(plugInLock, "LinearFrequency",    0.0,       1.0, 1.0,     0.0,   Parameter::BOOLEAN)      );
  addObservedParameter( new Parameter(plugInLock, "MinFrequency",      15.625, 32000.0, 0.0,    15.625, Parameter::EXPONENTIAL)  );
  addObservedParameter( new Parameter(plugInLock, "MaxFrequency",      15.625, 32000.0, 0.0, 32000.0,   Parameter::EXPONENTIAL)  );
  addObservedParameter( new Parameter(plugInLock, "MinLevel",        -100.0,      10.0, 0.0,  -100.0,   Parameter::LINEAR)       );
  addObservedParameter( new Parameter(plugInLock, "MaxLevel",        -100.0,      10.0, 0.0,    10.0,   Parameter::LINEAR)       );
  addObservedParameter( new Parameter(plugInLock, "FrameRate",        10.0,       50.0, 1.0,   15.0,    Parameter::LINEAR)       );
  addObservedParameter( new Parameter(plugInLock, "Freeze",            0.0,        1.0, 1.0,    1.0,    Parameter::BOOLEAN)      );

  for(int i=0; i < (int) observedParameters.size(); i++ )
    parameterChanged(observedParameters[i]);
}

//=================================================================================================
// class MultiAnalyzerAudioModule:

MultiAnalyzerAudioModule::MultiAnalyzerAudioModule(CriticalSection *newPlugInLock, 
                                                   OscilloscopeAudioModule *newOscilloscopeModule, 
                                                   SpectrumAnalyzerAudioModule *newSpectrumAnalyzerModule)
                                                    : AudioModule(newPlugInLock)
{
  oscilloscopeModule     = newOscilloscopeModule;
  spectrumAnalyzerModule = newSpectrumAnalyzerModule;

  jassert( oscilloscopeModule     != NULL );
  jassert( spectrumAnalyzerModule != NULL );

  addChildAudioModule(oscilloscopeModule);
  addChildAudioModule(spectrumAnalyzerModule);

  moduleName = juce::String(T("MultiAnalyzer"));
  setActiveDirectory(getApplicationDirectory() + juce::String(T("/MultiAnalyzerPresets")) );

  initializeAutomatableParameters();
}

void MultiAnalyzerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case 0: mode = (int) value;  break;
  } 

  markStateAsDirty(); // this feature is de-activated because the state will be marked as dirty immediately after preset-load due to some
  // async scrollBarMove callback - this renders the feature meaningless - fix this!
}

void MultiAnalyzerAudioModule::initializeAutomatableParameters()
{
  Parameter *p = new Parameter(plugInLock, "Mode", 0.0, 1.0, 1.0, 1.0, Parameter::STRING);
  p->addStringValue(juce::String(T("Oscilloscope")));
  p->addStringValue(juce::String(T("SpectrumAnalyzer")));
  p->setValue(1.0, false, false);
  addObservedParameter(p);

  //addObservedParameter( new Parameter("FrameRate",  10.0,  50.0, 1.0,  25.0, Parameter::LINEAR)  );

  for(int i=0; i < (int) observedParameters.size(); i++ )
    parameterChanged(observedParameters[i]);
}


