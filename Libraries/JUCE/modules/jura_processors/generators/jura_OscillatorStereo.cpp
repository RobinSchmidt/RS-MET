
//-------------------------------------------------------------------------------------------------
// construction/destruction:

OscillatorStereoAudioModule::OscillatorStereoAudioModule(CriticalSection *newPlugInLock,
  rosic::OscillatorStereo *newOscillatorStereoToWrap) : AudioModule(newPlugInLock)
{
  jassert( newOscillatorStereoToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  setModuleTypeName("OscillatorStereo");
  wrappedOscillatorStereo = newOscillatorStereoToWrap;
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void OscillatorStereoAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedOscillatorStereo == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case  0: wrappedOscillatorStereo->setLevel(                      value);  break;
  case  1: wrappedOscillatorStereo->setLevelByKey(                 value);  break;
  case  2: wrappedOscillatorStereo->setLevelByVel(                 value);  break;
  case  3: wrappedOscillatorStereo->setMidSide(                    value);  break;
  case  4: wrappedOscillatorStereo->setPan(                        value);  break;

  case  5: wrappedOscillatorStereo->setDetuneSemitones(            value);  break;
  case  6: wrappedOscillatorStereo->setDetuneHz(                   value);  break;
  case  7: wrappedOscillatorStereo->setStereoDetuneSemitones(      value);  break;
  case  8: wrappedOscillatorStereo->setStereoDetuneHz(             value);  break;
  case  9: wrappedOscillatorStereo->setPitchEnvelopeDepth(         value);  break;

  case 10: wrappedOscillatorStereo->setStartPhase(                 value);  break;
  case 11: wrappedOscillatorStereo->setFullWavePhaseWarp(          value);  break;
  case 12: wrappedOscillatorStereo->setHalfWavePhaseWarp(          value);  break;
  case 13: wrappedOscillatorStereo->waveTable->setCombHarmonic(    value);  break;
  case 14: wrappedOscillatorStereo->waveTable->setCombAmount(      value);  break;

  case 15: wrappedOscillatorStereo->setSpectralContrast(           value);  break;
  case 16: wrappedOscillatorStereo->setSpectralSlope(              value);  break;
  case 17: wrappedOscillatorStereo->setHighestHarmonicToKeep( (int)value);  break;
  case 18: wrappedOscillatorStereo->setLowestHarmonicToKeep(  (int)value);  break;
  case 19: wrappedOscillatorStereo->setEvenOddRatio(               value);  break;

  case 20: wrappedOscillatorStereo->waveTable->setPhaseScale(      value);  break;
  case 21: wrappedOscillatorStereo->waveTable->setPhaseShift(      value);  break;
  case 22: wrappedOscillatorStereo->setEvenOddPhaseShift(          value);  break;
  case 23: wrappedOscillatorStereo->setStereoPhaseShift(           value);  break;
  case 24: wrappedOscillatorStereo->setEvenOddStereoPhaseShift(    value);  break;
  } // end of switch( parameterIndex )
  markStateAsDirty();
}

//-------------------------------------------------------------------------------------------------
// state saving and recall:

XmlElement* oscillatorStereoStateToXml(OscillatorStereo* osc,
  XmlElement* xmlElementToStartFrom)
{
  // the XmlElement which stores all the releveant state-information:
  XmlElement* xmlState;
  if( xmlElementToStartFrom == NULL )
    xmlState = new XmlElement(juce::String("MultiModeFilterState"));  //wrong? should be OscillatorStereoState
  else
    xmlState = xmlElementToStartFrom;

  // store the settings in the XmlElement:
  xmlState->setAttribute("AudioFileRelativePath", juce::String(osc->waveTable->getSampleName() )    );
  xmlState->setAttribute("Mute",                    osc->isMuted()                              );
  xmlState->setAttribute("Level",                   osc->getLevel()                             );
  xmlState->setAttribute("LevelByKey",              osc->getLevelByKey()                        );
  xmlState->setAttribute("LevelByVel",              osc->getLevelByVel()                        );
  xmlState->setAttribute("Pan",                     osc->getPan()                               );
  xmlState->setAttribute("MidSideRatio",            osc->getMidSide()                           );
  xmlState->setAttribute("StartPhase",              osc->getStartPhase()                        );
  xmlState->setAttribute("StartPhaseByKey",         osc->getStartPhaseByKey()                   );
  xmlState->setAttribute("StartPhaseByVel",         osc->getStartPhaseByVel()                   );
  xmlState->setAttribute("TimeReverse",             osc->waveTable->isTimeReversed()            );
  xmlState->setAttribute("PolarityInvert",          osc->waveTable->isPolarityInverted()        );
  xmlState->setAttribute("FullWavePhaseWarp",       osc->waveTable->getFullWavePhaseWarp()      );
  xmlState->setAttribute("HalfWavePhaseWarp",       osc->waveTable->getHalfWavePhaseWarp()      );
  xmlState->setAttribute("CombHarmonic",            osc->waveTable->getCombHarmonic()           );
  xmlState->setAttribute("CombAmount",              osc->waveTable->getCombAmount()             );
  xmlState->setAttribute("Tune",                    osc->getDetuneSemitones()                   );
  xmlState->setAttribute("DetuneHz",                osc->getDetuneHz()                          );
  xmlState->setAttribute("StereoDetuneSemitones",   osc->getStereoDetuneSemitones()             );
  xmlState->setAttribute("StereoDetuneHz",          osc->getStereoDetuneHz()                    );
  xmlState->setAttribute("PitchModulationDepth",    osc->getPitchEnvelopeDepth()                );
  xmlState->setAttribute("SpectralContrast",        osc->waveTable->getSpectralContrast()       );
  xmlState->setAttribute("SpectralSlope",           osc->waveTable->getSpectralSlope()          );
  xmlState->setAttribute("HighestHarmonicToKeep",   osc->waveTable->getHighestHarmonicToKeep()  );
  xmlState->setAttribute("LowestHarmonicToKeep",    osc->waveTable->getLowestHarmonicToKeep()   );
  xmlState->setAttribute("PhaseScale",              osc->waveTable->getPhaseScale()             );
  xmlState->setAttribute("PhaseShift",              osc->waveTable->getPhaseShift()             );
  xmlState->setAttribute("EvenOddRatio",            osc->waveTable->getEvenOddRatio()           );
  xmlState->setAttribute("EvenOddPhaseShift",       osc->waveTable->getEvenOddPhaseShift()      );
  xmlState->setAttribute("StereoPhaseShift",        osc->waveTable->getStereoPhaseShift()       );
  xmlState->setAttribute("EvenOddStereoPhaseShift", osc->waveTable->getEvenOddStereoPhaseShift());

  return xmlState;
}

bool oscillatorStereoStateFromXml(OscillatorStereo* osc, const XmlElement &xmlState)
{
  bool success = true;

  // temporarily switch off the automatic re-rendering of the mip-map, to multiple rendering (one for each
  // parameter):
  bool oldAutoReRenderState = osc->waveTable->isMipMapAutoReRenderingActive();
  osc->waveTable->setAutomaticMipMapReRendering(false);

  // restore the settings from the XmlElement:
  osc->setMute(                     xmlState.getBoolAttribute(           "Mute",                  false));
  osc->setLevel(                    xmlState.getDoubleAttribute(         "Level",                   0.0));
  osc->setLevelByKey(               xmlState.getDoubleAttribute(         "LevelByKey",              0.0));
  osc->setLevelByVel(               xmlState.getDoubleAttribute(         "LevelByVel",              0.0));
  osc->setPan(                      xmlState.getDoubleAttribute(         "Pan",                     0.0));
  osc->setMidSide(                  xmlState.getDoubleAttribute(         "MidSideRatio",            0.5));
  osc->setStartPhase(               xmlState.getDoubleAttribute(         "StartPhase",              0.0));
  osc->setStartPhaseByKey(          xmlState.getDoubleAttribute(         "StartPhaseByKey",         0.0));
  osc->setStartPhaseByVel(          xmlState.getDoubleAttribute(         "StartPhaseByVel",         0.0));
  osc->setTimeReverse(              xmlState.getBoolAttribute(           "TimeReverse",           false));
  osc->setPolarityInversion(        xmlState.getBoolAttribute(           "PolarityInvert",        false));
  osc->setFullWavePhaseWarp(        xmlState.getDoubleAttribute(         "FullWavePhaseWarp",       0.0));
  osc->setHalfWavePhaseWarp(        xmlState.getDoubleAttribute(         "HalfWavePhaseWarp",       0.0));
  osc->waveTable->setCombHarmonic(  xmlState.getDoubleAttribute(         "CombHarmonic",            1.0));
  osc->waveTable->setCombAmount(    xmlState.getDoubleAttribute(         "CombAmount",              0.0));
  osc->setDetuneSemitones(          xmlState.getDoubleAttribute(         "Tune",                    0.0));
  osc->setDetuneHz(                 xmlState.getDoubleAttribute(         "DetuneHz",                0.0));
  osc->setStereoDetuneSemitones(    xmlState.getDoubleAttribute(         "StereoDetuneSemitones",   0.0));
  osc->setStereoDetuneHz(           xmlState.getDoubleAttribute(         "StereoDetuneHz",          0.0));
  osc->setPitchEnvelopeDepth(       xmlState.getDoubleAttribute(         "PitchModulationDepth",    0.0));
  osc->waveTable->setSpectralContrast( xmlState.getDoubleAttribute(      "SpectralContrast",        1.0));
  osc->waveTable->setSpectralSlope(    xmlState.getDoubleAttribute(      "SpectralSlope",           0.0));
  osc->waveTable->setHighestHarmonicToKeep(xmlState.getIntAttribute(     "HighestHarmonicToKeep",  1024));
  osc->waveTable->setLowestHarmonicToKeep( xmlState.getIntAttribute(     "LowestHarmonicToKeep",      1));
  osc->waveTable->setPhaseScale( xmlState.getDoubleAttribute(            "PhaseScale",              1.0));
  osc->waveTable->setPhaseShift( xmlState.getDoubleAttribute(            "PhaseShift",              0.0));
  osc->waveTable->setEvenOddRatio( xmlState.getDoubleAttribute(          "EvenOddRatio",            0.5));
  osc->waveTable->setEvenOddPhaseShift(xmlState.getDoubleAttribute(      "EvenOddPhaseShift",       0.0));
  osc->waveTable->setStereoPhaseShift( xmlState.getDoubleAttribute(      "StereoPhaseShift",        0.0));
  osc->waveTable->setEvenOddStereoPhaseShift(xmlState.getDoubleAttribute("EvenOddStereoPhaseShift", 0.0));

  // let the oscillator and all its slaves calculate their increment:
  //osc->setFrequency(1000.0);
  osc->calculateIncrementForAllSlaves();

  // load the audio-file into the wavetable for the oscillator:
  juce::String samplePath = xmlState.getStringAttribute("AudioFileRelativePath", juce::String::empty);
  AudioSampleBuffer* buffer = AudioFileManager::createAudioSampleBufferFromFile(samplePath, true);
  if( buffer != NULL )
  {
    // pass the actual audio data:
    float* channelPointers[2];
    channelPointers[0] = buffer->getWritePointer(0, 0);
    if( buffer->getNumChannels() >= 2 )
      channelPointers[1] = buffer->getWritePointer(1, 0);
    else
      channelPointers[1] = buffer->getWritePointer(0, 0);
    osc->waveTable->setWaveform(channelPointers, buffer->getNumSamples() );

    // pass the path as c-string:
    /*
    long  length    = samplePath.length();
    char* fileNameC = new char[length+1];
    samplePath.copyToBuffer(fileNameC, length);
    osc->waveTable->setSampleName(fileNameC);
    */
    char* fileNameC = toZeroTerminatedString(samplePath);
    osc->waveTable->setSampleName(fileNameC);
    delete[] fileNameC;

    delete buffer;
    success = true;
  }
  else
  {
    juce::String errorString = samplePath + juce::String(" !error!");

    /*
    long  length    = errorString.length();
    char* fileNameC = new char[length+1];
    errorString.copyToBuffer(fileNameC, length);
    */
    char* fileNameC = toZeroTerminatedString(errorString);
    osc->waveTable->setSampleName(fileNameC);
    delete[] fileNameC;

    osc->waveTable->fillWithAllZeros();
    success = false;
  }

  // let the (wavetable inside the) oscillator render the mip map and restore the old state of the
  // automatic re-rendering:
  osc->waveTable->renderMipMap();
  osc->waveTable->setAutomaticMipMapReRendering(oldAutoReRenderState);

  return success;
}

XmlElement* OscillatorStereoAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  // store the parameters of the underlying core object:
  if( wrappedOscillatorStereo != NULL )
    xmlState = oscillatorStereoStateToXml(wrappedOscillatorStereo, xmlState);

  return xmlState;
}

void OscillatorStereoAudioModule::setStateFromXml(const XmlElement& xmlState,
                                                 const juce::String& stateName, bool markAsClean)
{
  // restore the inherited controller mappings:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  // restore the parameters of the underlying core object:
  if( wrappedOscillatorStereo != NULL )
    oscillatorStereoStateFromXml(wrappedOscillatorStereo, xmlState);
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void OscillatorStereoAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  std::vector<double> defaultValues;

  // this pointer will be used to temporarily store the addresses of the created
  //AutomatableParameter* p;
  Parameter* p;

  // amplitude related parameters:

  // #000:
  p = new AutomatableParameter(lock, "Level",    -36.0, 12.0, 0.01, 0.0, Parameter::LINEAR);
  defaultValues.push_back(-3.01029996);   // compensation gain for  2 uncorrelated sources
  defaultValues.push_back(-4.77121255);   // compensation gain for  3 uncorrelated sources
  defaultValues.push_back(-6.02059991);   // compensation gain for  4 uncorrelated or 2 in-phase sources
  defaultValues.push_back(-9.03089987);   // compensation gain for  8 uncorrelated
  defaultValues.push_back(-12.0411998);   // compensation gain for 16 uncorrelated or 4 in-phase sources
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #001:
  p = new AutomatableParameter(lock, "LevelByKey", -24.0, 24.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #002:
  p = new AutomatableParameter(lock, "LevelByVel", -12.0, 12.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #003:
  p = new AutomatableParameter(lock, "MidSide",    0.0,  1.0, 0.01, 0.5, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(0.0);
  defaultValues.push_back(0.25);
  defaultValues.push_back(0.5);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #004:
  p = new AutomatableParameter(lock, "Pan",       -1.0,  1.0, 0.01, 0.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(-1.0);
  defaultValues.push_back(-0.75);
  defaultValues.push_back(-0.5);
  defaultValues.push_back(-0.25);
  defaultValues.push_back(0.0);
  defaultValues.push_back(0.25);
  defaultValues.push_back(0.5);
  defaultValues.push_back(0.75);
  defaultValues.push_back(1.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // tuning related parameters:

  // #005:
  p = new AutomatableParameter(lock, "Tune",     -36.0, 36.0, 0.01, 0.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(-24.0);      // f/f0 =  1/ 4 = 0.25
  defaultValues.push_back(-12.0);      // f/f0 =  1/ 2 = 0.50
  defaultValues.push_back(2.03910002); // f/f0 =  9/ 8 = 1.125
  defaultValues.push_back(3.15641287); // f/f0         = 1.2
  defaultValues.push_back(3.86313714); // f/f0 =  5/ 4 = 1.25
  //defaultValues.push_back(4.54213948); // f/f0 = 13/10 = 1.3
  defaultValues.push_back(4.98044999); // f/f0 =  4/ 3 = 1.333...
  defaultValues.push_back(7.01955001); // f/f0 =  3/ 2 = 1.5
  defaultValues.push_back(8.84358713); // f/f0 =  5/ 3 = 1.666...
  defaultValues.push_back(9.68825906); // f/f0 =  7/ 4 = 1.75
  defaultValues.push_back(12.0);       // f/f0 =  2/ 1 = 2.00
  defaultValues.push_back(24.0);       // f/f0 =  4/ 1 = 4.00
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #006:
  p = new AutomatableParameter(lock, "DetuneHz", -20.0, 20.0, 0.01, 0.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(-4.0);
  defaultValues.push_back(-3.0);
  defaultValues.push_back(-2.0);
  defaultValues.push_back(-1.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(1.0);
  defaultValues.push_back(2.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(4.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #007:
  p = new AutomatableParameter(lock, "StereoDetune", -1.0, 1.0, 0.01, 0.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(-0.2);
  defaultValues.push_back(-0.1);
  defaultValues.push_back(0.0);
  defaultValues.push_back(0.1);
  defaultValues.push_back(0.2);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #008:
  p = new AutomatableParameter(lock, "StereoDetuneHz", -10.0, 10.0, 0.01, 0.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(-4.0);
  defaultValues.push_back(-3.0);
  defaultValues.push_back(-2.0);
  defaultValues.push_back(-1.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(1.0);
  defaultValues.push_back(2.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(4.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #009:
  p = new AutomatableParameter(lock, "PitchModulationDepth", -8.0, 8.0, 0.01, 0.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(-1.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(1.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // time domain waveform related parameters:

  // #010:
  p = new Parameter(lock, "StartPhase", 0.0, 360.0, 1.0, 0.01, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(0.0);
  defaultValues.push_back(45.0);
  defaultValues.push_back(90.0);
  defaultValues.push_back(135.0);
  defaultValues.push_back(180.0);
  defaultValues.push_back(225.0);
  defaultValues.push_back(270.0);
  defaultValues.push_back(315.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #011:
  p = new Parameter(lock, "FullWaveWarp", -0.99, 0.99, 0.001, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #012:
  p = new Parameter(lock, "HalfWaveWarp", -0.99, 0.99, 0.001, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #013:
  /*
  p = new AutomatableParameter("CombOffset", 0.0, 360.0, 0.1, 0.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(0.0);
  defaultValues.push_back(45.0);
  defaultValues.push_back(90.0);
  defaultValues.push_back(135.0);
  defaultValues.push_back(180.0);
  defaultValues.push_back(225.0);
  defaultValues.push_back(270.0);
  defaultValues.push_back(315.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
  */

  // #014:
  //p = new Parameter("CombHarmonic", 0.1, 1024.0, 0.001, 0.0, Parameter::EXPONENTIAL);
  p = new Parameter(lock, "CombHarmonic", 1.0, 128.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #014:
  p = new Parameter(lock, "CombAmount", -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // magnitude spectrum related parameters:

  // #015:
  p = new Parameter(lock, "SpectralContrast", 0.25, 4.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #016:
  p = new Parameter(lock, "SpectralSlope", -6.0, 6.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #017:
  p = new Parameter(lock, "HighestHarmonic", 1.0, 1024.0, 1.0, 1024.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #018:
  p = new Parameter(lock, "LowestHarmonic", 1.0, 1024.0, 1.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #019:
  p = new Parameter(lock, "EvenOddRatio", 0.0, 1.0, 0.005, 0.5, Parameter::LINEAR);
  addObservedParameter(p);

  //-----------------------------------------------------------------------------------------------
  // phase spectrum related parameters:

  // #020:
  p = new Parameter(lock, "PhaseScale", -1.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #021:
  p = new Parameter(lock, "PhaseShift", -180.0, 180.0, 1.0, 0.0, Parameter::LINEAR);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #022:
  p = new Parameter(lock, "EvenOddPhaseShift", -180.0, 180.0, 1.0, 0.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(-90.0);
  defaultValues.push_back(-45.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(45.0);
  defaultValues.push_back(90.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #023:
  p = new Parameter(lock, "StereoPhaseShift", -180.0, 180.0, 1.0, 0.0, Parameter::LINEAR);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #024:
  p = new Parameter(lock, "EvenOddStereoPhaseShift", -180.0, 180.0, 1.0, 0.0, Parameter::LINEAR);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);


  // make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================

// construction/destruction:

OscillatorStereoEditorContextMenu::OscillatorStereoEditorContextMenu(
  OscillatorStereoAudioModule* newOscillatorModuleToEdit, Component* componentToAttachTo)
  //: ComponentMovementWatcher(componentToAttachTo)
{
  // init the pointer to the oscillatoror to be edited:
  //jassert( newOscillatorToEdit != NULL )
  oscillatorModuleToEdit = newOscillatorModuleToEdit;

  addWidget( ampHeadline = new RTextField("Amplitude:") );
  ampHeadline->setNoBackgroundAndOutline(true);
  ampHeadline->setDescription(juce::String("Manipulations of the amplitude"));

  addWidget( tuningHeadline = new RTextField( "Tuning:") );
  tuningHeadline->setNoBackgroundAndOutline(true);
  tuningHeadline->setDescription(juce::String("Manipulations of the tuning/detuning of the oscillator"));

  addWidget( timeHeadline = new RTextField("Time:") );
  timeHeadline->setNoBackgroundAndOutline(true);
  timeHeadline->setDescription(juce::String("Time domain manipulations of the waveform"));

  addWidget( magSpectrumHeadline = new RTextField("Magnitude Spectrum:") );
  magSpectrumHeadline->setNoBackgroundAndOutline(true);
  magSpectrumHeadline->setDescription(juce::String("Manipulations of the magnitude spectrum"));

  addWidget( phaseSpectrumHeadline = new RTextField("Phase Spectrum:") );
  phaseSpectrumHeadline->setNoBackgroundAndOutline(true);
  phaseSpectrumHeadline->setDescription(juce::String("Manipulations of the phase spectrum"));

  // sliders for amplitude related parameters:

  addWidget( levelSlider = new RSlider("LevelSlider") );
  levelSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("Level") );
  levelSlider->setDescription(juce::String("Output level of the oscillator"));
  levelSlider->setStringConversionFunction(&decibelsToStringWithUnit2);
  levelSlider->addListener(this); // to send out the change-message for display update

  addWidget( levelByKeySlider = new RSlider("LevelByKeySlider") );
  levelByKeySlider->assignParameter( oscillatorModuleToEdit->getParameterByName("LevelByKey") );
  levelByKeySlider->setSliderName(juce::String("Key"));
  levelByKeySlider->setDescription(juce::String("Key dependence of oscillator's output level"));
  levelByKeySlider->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( levelByVelSlider = new RSlider("LevelVelSlider") );
  levelByVelSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("LevelByVel") );
  levelByVelSlider->setSliderName(juce::String("Vel"));
  levelByVelSlider->setDescription(juce::String("Velocity dependence of oscillator's output level"));
  levelByVelSlider->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( midSideSlider = new RSlider("MidSideSlider") );
  midSideSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("MidSide") );
  midSideSlider->setSliderName(juce::String("Mid/Side"));
  midSideSlider->setDescription("Mid/side adjustment for stereo(ized) waveforms");
  midSideSlider->setStringConversionFunction(&ratioToString0);
  midSideSlider->addListener(this);

  addWidget( panSlider = new RSlider("PanSlider") );
  panSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("Pan") );
  panSlider->setSliderName(juce::String("Pan"));
  panSlider->setDescription("Panorama position of the oscillator");
  panSlider->setStringConversionFunction(&valueToString2);
  panSlider->addListener(this);

  // sliders for tuning related parameters:

  addWidget( tuneSlider = new TuningSlider("TuneSlider") );
  tuneSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("Tune") );
  tuneSlider->setDescription(juce::String("Tuning of the oscillator in semitones"));
  tuneSlider->setStringConversionFunction(&semitonesToStringWithUnit2);

  addWidget( detuneHzSlider = new RSlider("DetuneSlider") );
  detuneHzSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("DetuneHz") );
  detuneHzSlider->setSliderName(juce::String("Detune Hz"));
  detuneHzSlider->setDescription(juce::String("Detuning of the oscillator in Hz"));
  detuneHzSlider->setStringConversionFunction(&hertzToStringWithUnit2);

  addWidget( stereoDetuneSlider = new RSlider("StereoDetuneHzSlider") );
  stereoDetuneSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("StereoDetune") );
  stereoDetuneSlider->setSliderName(juce::String("Stereo Detune"));
  stereoDetuneSlider->setDescription(juce::String("Detuning between left and right channel in semitones"));
  stereoDetuneSlider->setStringConversionFunction(&semitonesToStringWithUnit2);

  addWidget( stereoDetuneHzSlider = new RSlider("StereoDetuneHzSlider") );
  stereoDetuneHzSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("StereoDetuneHz") );
  stereoDetuneHzSlider->setSliderName(juce::String("Stereo Detune Hz"));
  stereoDetuneHzSlider->setDescription(juce::String("Detuning between left and right channel in Hz"));
  stereoDetuneHzSlider->setStringConversionFunction(&hertzToStringWithUnit2);

  // sliders for time domain related parameters:

  addWidget( startPhaseSlider = new RSlider("StartPhaseSlider") );
  startPhaseSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("StartPhase") );
  startPhaseSlider->setSliderName(juce::String("Start Phase"));
  startPhaseSlider->setDescription("Start phase of the oscillator");
  startPhaseSlider->setStringConversionFunction(&degreesToStringWithUnit0);
  startPhaseSlider->addListener(this);

  addWidget( fullWavePhaseWarpSlider = new RSlider("FullWavePhaseWarpSlider") );
  fullWavePhaseWarpSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("FullWaveWarp") );
  fullWavePhaseWarpSlider->setSliderName(juce::String("Full Wave Warp"));
  fullWavePhaseWarpSlider->setDescription("Applies phase warping to the entire waveform");
  fullWavePhaseWarpSlider->setStringConversionFunction(&degreesToStringWithUnit0);
  fullWavePhaseWarpSlider->addListener(this);

  addWidget( halfWavePhaseWarpSlider = new RSlider("HalfWavePhaseWarpSlider") );
  halfWavePhaseWarpSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("HalfWaveWarp") );
  halfWavePhaseWarpSlider->setSliderName(juce::String("Half Wave Warp"));
  halfWavePhaseWarpSlider->setDescription("Applies phase warping both half cycles of the waveform");
  halfWavePhaseWarpSlider->setStringConversionFunction(&degreesToStringWithUnit0);
  halfWavePhaseWarpSlider->addListener(this);

  addWidget( combHarmonicSlider = new RSlider("CombHarmonicSlider") );
  combHarmonicSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("CombHarmonic") );
  combHarmonicSlider->setSliderName(juce::String("Comb Harmonic"));
  combHarmonicSlider->setDescription("Harmonic on which the comb filter acts");
  combHarmonicSlider->setStringConversionFunction(&degreesToStringWithUnit0);
  combHarmonicSlider->addListener(this);
  //combHarmonicSlider->setVisible(false); // not yet meaningfully implemented

  addWidget( combAmountSlider = new RSlider("CombAmountSlider") );
  combAmountSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("CombAmount") );
  combAmountSlider->setSliderName(juce::String("Comb Amount"));
  combAmountSlider->setDescription("Amount of comb filtering");
  combAmountSlider->setStringConversionFunction(percentToStringWithUnit1);
  combAmountSlider->addListener(this);
  //combAmountSlider->setVisible(false); // not yet meaningfully implemented

  addWidget( reverseButton = new RButton(juce::String("Reverse")) );
  reverseButton->addRButtonListener(this);
  reverseButton->setDescription(juce::String("Time reverses the oscillator's waveform"));
  reverseButton->setClickingTogglesState(true);

  addWidget( invertButton = new RButton(juce::String("Invert")) );
  invertButton->addRButtonListener(this);
  invertButton->setDescription(juce::String("Inverts polarity of the oscillator's ouput"));
  invertButton->setClickingTogglesState(true);

  // sliders for magnitude spectrum related parameters:

  addWidget( spectralContrastSlider = new RSlider("SpectralContrastSlider") );
  spectralContrastSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("SpectralContrast") );
  spectralContrastSlider->setSliderName(juce::String("Contrast"));
  spectralContrastSlider->setDescription("Spectral contrast acting as exponent on the harmonic's magnitude");
  spectralContrastSlider->setStringConversionFunction(&valueToString2);
  spectralContrastSlider->addListener(this);

  addWidget( spectralSlopeSlider = new RSlider("SpectralSlopeSlider") );
  spectralSlopeSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("SpectralSlope") );
  spectralSlopeSlider->setSliderName(juce::String("Slope"));
  spectralSlopeSlider->setDescription("Spectral slope applied to the waveform in dB/oct");
  spectralSlopeSlider->setStringConversionFunction(&decibelsPerOctaveToString2);
  spectralSlopeSlider->addListener(this);

  addWidget( highestHarmonicSlider = new RSlider("LowpassSlider") );
  highestHarmonicSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("HighestHarmonic") );
  highestHarmonicSlider->setSliderName(juce::String("Highest Harmonic"));
  highestHarmonicSlider->setDescription("Highest harmonic in the waveform");
  highestHarmonicSlider->setStringConversionFunction(&valueToString0);
  highestHarmonicSlider->addListener(this);

  addWidget( lowestHarmonicSlider = new RSlider("HighpassSlider") );
  lowestHarmonicSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("LowestHarmonic") );
  lowestHarmonicSlider->setSliderName(juce::String("Lowest Harmonic"));
  lowestHarmonicSlider->setDescription("Lowest harmonic in the waveform");
  lowestHarmonicSlider->setStringConversionFunction(&valueToString0);
  lowestHarmonicSlider->addListener(this);

  addWidget( evenOddSlider = new RSlider("EvenOddSlider") );
  evenOddSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("EvenOddRatio") );
  evenOddSlider->setSliderName(juce::String("Even/Odd"));
  evenOddSlider->setDescription("Ratio of even and odd harmonics");
  evenOddSlider->setStringConversionFunction(&ratioBothFullAtCenterToString0);
  evenOddSlider->addListener(this);

  // sliders for phase spectrum related parameters:

  addWidget( evenOddPhaseShiftSlider = new RSlider("EvenOddPhaseShiftSlider") );
  evenOddPhaseShiftSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("EvenOddPhaseShift") );
  evenOddPhaseShiftSlider->setSliderName(juce::String("Even/Odd Shift:"));
  evenOddPhaseShiftSlider->setDescription("Applies a phase shift between even and odd harmonics");
  evenOddPhaseShiftSlider->setStringConversionFunction(&degreesToStringWithUnit0);
  evenOddPhaseShiftSlider->addListener(this);

  addWidget( phaseScaleSlider = new RSlider("PhaseScaleSlider") );
  phaseScaleSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("PhaseScale") );
  phaseScaleSlider->setSliderName(juce::String("Scale:"));
  phaseScaleSlider->setDescription("Scales the phase of each harmonic ");
  phaseScaleSlider->setStringConversionFunction(&degreesToStringWithUnit0);
  phaseScaleSlider->addListener(this);

  addWidget( phaseShiftSlider = new RSlider("PhaseShiftSlider") );
  phaseShiftSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("PhaseShift") );
  phaseShiftSlider->setSliderName(juce::String("Shift:"));
  phaseShiftSlider->setDescription("Shifts the phase of each harmonic by a constant");
  phaseShiftSlider->setStringConversionFunction(&degreesToStringWithUnit0);
  phaseShiftSlider->addListener(this);

  addWidget( stereoPhaseShiftSlider = new RSlider("StereoPhaseShiftSlider") );
  stereoPhaseShiftSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("StereoPhaseShift") );
  stereoPhaseShiftSlider->setSliderName(juce::String("Stereo Shift:"));
  stereoPhaseShiftSlider->setDescription("Applies stereoization via phase-shifting of harmonics");
  stereoPhaseShiftSlider->setStringConversionFunction(&degreesToStringWithUnit0);
  stereoPhaseShiftSlider->addListener(this);

  addWidget( evenOddStereoPhaseShiftSlider = new RSlider("EvenOddPhaseShiftSlider") );
  evenOddStereoPhaseShiftSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("EvenOddStereoPhaseShift") );
  evenOddStereoPhaseShiftSlider->setSliderName(juce::String("Even/Odd Stereo Shift:"));
  evenOddStereoPhaseShiftSlider->setDescription("Phase shift between even/odd harmonics, applied with opposite signs to left/right channels");
  evenOddStereoPhaseShiftSlider->setStringConversionFunction(&degreesToStringWithUnit0);
  evenOddStereoPhaseShiftSlider->addListener(this);

  addWidget( closeButton = new RButton(RButton::CLOSE) );
  closeButton->setDescription(juce::String("Closes the oscillator context menu"));
  closeButton->setClickingTogglesState(false);
  // we don't listen to this button ourselves - this is the job of the outlying editor object

  updateWidgetsAccordingToState();
  setSize(180, 488);
}

OscillatorStereoEditorContextMenu::~OscillatorStereoEditorContextMenu()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void OscillatorStereoEditorContextMenu::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( oscillatorModuleToEdit == NULL )
    return;
  if( oscillatorModuleToEdit->wrappedOscillatorStereo == NULL )
    return;

  rosic::OscillatorStereo* o = oscillatorModuleToEdit->wrappedOscillatorStereo;

  if( buttonThatWasClicked == reverseButton )
    o->setTimeReverse( reverseButton->getToggleState() );
  if( buttonThatWasClicked == invertButton )
    o->setPolarityInversion( invertButton->getToggleState() );

  sendChangeMessage();
}

void OscillatorStereoEditorContextMenu::rSliderValueChanged(RSlider *sliderThatHasChanged)
{
  if( oscillatorModuleToEdit == NULL )
    return;
  if( oscillatorModuleToEdit->wrappedOscillatorStereo == NULL )
    return;
  //rosic::OscillatorStereo* o = oscillatorModuleToEdit->wrappedOscillatorStereo;
  sendChangeMessage();
}

void OscillatorStereoEditorContextMenu::updateWidgetsAccordingToState()
{
  if( oscillatorModuleToEdit == NULL )
    return;
  if( oscillatorModuleToEdit->wrappedOscillatorStereo == NULL )
    return;
  rosic::OscillatorStereo* o         = oscillatorModuleToEdit->wrappedOscillatorStereo;
  rosic::MipMappedWaveTableStereo* w = o->waveTable;

  // update the widgets:
  levelSlider->setValue(                  o->getLevel(),                   false);
  levelByKeySlider->setValue(             o->getLevelByKey(),              false);
  levelByVelSlider->setValue(             o->getLevelByVel(),              false);
  panSlider->setValue(                    o->getPan(),                     false);
  midSideSlider->setValue(                o->getMidSide(),                 false);
  startPhaseSlider->setValue(             o->getStartPhase(),              false);
  fullWavePhaseWarpSlider->setValue(      w->getFullWavePhaseWarp(),       false);
  halfWavePhaseWarpSlider->setValue(      w->getHalfWavePhaseWarp(),       false);
  combHarmonicSlider->setValue(           w->getCombHarmonic(),            false);
  combAmountSlider->setValue(             w->getCombAmount(),              false);
  tuneSlider->setValue(                   o->getDetuneSemitones(),         false);
  detuneHzSlider->setValue(               o->getDetuneHz(),                false);
  stereoDetuneSlider->setValue(           o->getStereoDetuneSemitones(),   false);
  stereoDetuneHzSlider->setValue(         o->getStereoDetuneHz(),          false);
  spectralContrastSlider->setValue(       w->getSpectralContrast(),        false);
  spectralSlopeSlider->setValue(          w->getSpectralSlope(),           false);
  highestHarmonicSlider->setValue(        w->getHighestHarmonicToKeep(),   false);
  lowestHarmonicSlider->setValue(         w->getLowestHarmonicToKeep(),    false);
  evenOddSlider->setValue(                w->getEvenOddRatio(),            false);
  phaseScaleSlider->setValue(             w->getPhaseScale(),              false);
  phaseShiftSlider->setValue(             w->getPhaseShift(),              false);
  stereoPhaseShiftSlider->setValue(       w->getStereoPhaseShift(),        false);
  evenOddPhaseShiftSlider->setValue(      w->getEvenOddPhaseShift(),       false);
  evenOddStereoPhaseShiftSlider->setValue(w->getEvenOddStereoPhaseShift(), false);
  reverseButton->setToggleState(          w->isTimeReversed(),             false);
  invertButton->setToggleState(           w->isPolarityInverted(),         false);
}

void OscillatorStereoEditorContextMenu::resized()
{
  Component::resized();
  int x  = 0;
  int y  = 0;
  int w  = getWidth();
  int w2 = w/2;
  //int h  = getHeight();

  closeButton->setBounds(w-16, 0, 16, 16);

  int sh   = 16;    // slider height
  int inc  = sh+4;
  //int inc2 = inc+8;

  ampHeadline->setBounds(x+4, y+4, w-8-16, sh);
  y += inc;
  levelSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  levelByKeySlider->setBounds(x+4,  y+4, w2-8, sh);
  levelByVelSlider->setBounds(w2+4, y+4, w2-8, sh);
  y += inc;
  midSideSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  panSlider->setBounds(x+4, y+4, w-8, sh);
  y += inc;

  tuningHeadline->setBounds(x+4, y+4, w-8, sh);
  y += inc;
  tuneSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  detuneHzSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  stereoDetuneSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  stereoDetuneHzSlider->setBounds(x+4, y+4, w-8, sh);
  y += inc;
  timeHeadline->setBounds(x+4, y+4, w-8, sh);
  y += inc;
  startPhaseSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  combHarmonicSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  combAmountSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  fullWavePhaseWarpSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  halfWavePhaseWarpSlider->setBounds(x+4, y+4, w-8, sh);

  //y += sh;
  //startPhaseByKeySlider->setBounds(x+4,  y+4, w2-8, sh);
  //startPhaseByVelSlider->setBounds(w2+4, y+4, w2-8, sh);
  y += inc;
  reverseButton->setBounds(x+4,  y+4, w2-8, sh);
  invertButton->setBounds( w2+4, y+4, w2-8, sh);
  y += inc;

  magSpectrumHeadline->setBounds(x+4, y+4, w-8, sh);
  y += inc;
  spectralContrastSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  spectralSlopeSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  highestHarmonicSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  lowestHarmonicSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  evenOddSlider->setBounds(x+4, y+4, w-8, sh);
  y += inc;


  phaseSpectrumHeadline->setBounds(x+4, y+4, w-8, sh);
  y += inc;
  phaseScaleSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  phaseShiftSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  evenOddPhaseShiftSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  stereoPhaseShiftSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  evenOddStereoPhaseShiftSlider->setBounds(x+4, y+4, w-8, sh);
}

/*
void OscillatorStereoEditorContextMenu::componentMovedOrResized(bool wasMoved, bool wasResized)
{
int dummy = 0;
}

void OscillatorStereoEditorContextMenu::componentPeerChanged()
{
int dummy = 0;
}
*/

//=================================================================================================
// class OscillatorStereoEditor:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

OscillatorStereoEditor::OscillatorStereoEditor(CriticalSection *newPlugInLock,
  OscillatorStereoAudioModule* newOscillatorStereoAudioModule)
  : AudioModuleEditor(newOscillatorStereoAudioModule)
{
  // init the pointer to the modulator to be edited to NULL:
  oscillatorToEdit = NULL;
  jassert( newOscillatorStereoAudioModule != NULL );



  addPlot( waveformDisplay = new WaveformDisplayOld() );
  waveformDisplay->setAutoReRendering(false);
  waveformDisplay->setDescription(juce::String("Click on the display to switch oscillator on/off"));
  waveformDisplay->setAxisPositionX(CoordinateSystemOld::INVISIBLE);
  waveformDisplay->setAxisPositionY(CoordinateSystemOld::INVISIBLE);
  waveformDisplay->setCurrentRangeY(-1.2, 1.2);
  waveformDisplay->setSampleRate(1.0);
  waveformDisplay->addMouseListener(this, true);
  waveformDisplay->setAutoReRendering(true);

  addPlot( emptyDisplay = new CoordinateSystemOld() );
  emptyDisplay->setAutoReRendering(false);
  emptyDisplay->setDescription(waveformDisplay->getDescription());
  emptyDisplay->setAxisPositionX(CoordinateSystemOld::INVISIBLE);
  emptyDisplay->setAxisPositionY(CoordinateSystemOld::INVISIBLE);
  emptyDisplay->setCaption(juce::String("Off"), CoordinateSystemOld::CENTER);
  emptyDisplay->addMouseListener(this, true);
  emptyDisplay->setAutoReRendering(true);

  contextMenu = new OscillatorStereoEditorContextMenu(newOscillatorStereoAudioModule, this);
  contextMenu->addChangeListener(this);
  contextMenu->setAlwaysOnTop(true);
  contextMenu->setOpaque(true);
  contextMenu->closeButton->addRButtonListener(this);
  addChildColourSchemeComponent(contextMenu, false, false);
  //contextMenu->setSize(200, 200);

  /*
  contextMenuViewport = new Viewport(juce::String(T("OscillatorStereoSliderViewport")));
  contextMenuViewport->setViewedComponent(contextMenu);
  contextMenuViewport->setScrollBarThickness(12);
  contextMenuViewport->setScrollBarsShown(true, false);
  addAndMakeVisible(contextMenuViewport);
  */

  addWidget( waveFileLabel = new RTextField() );
  waveFileLabel->setNoBackgroundAndOutline(true);
  waveFileLabel->setJustification(Justification::centredBottom);
  waveFileLabel->setDescription(juce::String("Name of the currently loaded single cycle audiofile"));

  addWidget( waveLoadButton = new RButton(juce::String("Load")) );
  waveLoadButton->addRButtonListener(this);
  waveLoadButton->setDescription(juce::String("Load a waveform"));
  waveLoadButton->setClickingTogglesState(false);

  addWidget( wavePlusButton = new RButton(RButton::ARROW_RIGHT) );
  wavePlusButton->addRButtonListener(this);
  wavePlusButton->setDescription(juce::String("Next waveform in current directory"));
  wavePlusButton->setClickingTogglesState(false);

  addWidget( waveMinusButton = new RButton(RButton::ARROW_LEFT) );
  waveMinusButton->addRButtonListener(this);
  waveMinusButton->setDescription(juce::String("Previous waveform in current directory"));
  waveMinusButton->setClickingTogglesState(false);

  addWidget( moreButton = new RButton(juce::String("More")) );
  moreButton->addRButtonListener(this);
  moreButton->setDescription(juce::String("Open/close context menu with more options"));
  moreButton->setClickingTogglesState(true);

  addWidget( levelSlider = new RSlider("VolumeSlider") );
  levelSlider->assignParameter( moduleToEdit->getParameterByName("Level") );
  levelSlider->setSliderName(juce::String("Level"));
  levelSlider->setDescription(juce::String("Output level of the oscillator"));
  levelSlider->setStringConversionFunction(&decibelsToStringWithUnit2);
  levelSlider->setLayout(RSlider::NAME_ABOVE);
  levelSlider->addListener(this); // only to update the plot

  addWidget( tuneSlider = new TuningSlider("TuneSlider") );
  tuneSlider->assignParameter( moduleToEdit->getParameterByName("Tune") );
  tuneSlider->setSliderName(juce::String("Tune"));
  tuneSlider->setDescription(juce::String("Tuning of the oscillator in semitones"));
  tuneSlider->setStringConversionFunction(&semitonesToStringWithUnit2);
  tuneSlider->setLayout(RSlider::NAME_ABOVE);

  addWidget( pitchModulationSlider = new RSlider("PitchModulationSlider") );
  pitchModulationSlider->assignParameter( moduleToEdit->getParameterByName("PitchModulationDepth") );
  pitchModulationSlider->setSliderName(juce::String("Mod"));
  pitchModulationSlider->setDescription("Modulation depth for the pitch modulator");
  pitchModulationSlider->setStringConversionFunction(&valueToString2);

  numSamplesInPlot    = 0;
  waveformBuffer      = NULL;
  waveformPointers    = new double*[2];
  waveformPointers[0] = NULL;
  waveformPointers[1] = NULL;

  isTopLevelEditor = false;
  setHeadlineStyle(NO_HEADLINE);

  // initialize the current directory for waveform loading and saving:
  juce::String appDir = getApplicationDirectory();
  AudioFileManager::setActiveDirectory(appDir + juce::String("/Samples/SingleCycle/Classic") );

  //setOscillatorToEdit(newOscillatorStereoAudioModule->wrappedOscillatorStereo);
  // this will also set up the widgets according to the state of the oscillator

  oscillatorToEdit = newOscillatorStereoAudioModule->wrappedOscillatorStereo; // get rid of this...
  updateWidgetsAccordingToState();

  // widget arrangement is optimized for this size:
  //setSize(232, 136);
}

OscillatorStereoEditor::~OscillatorStereoEditor()
{
  delete contextMenu; // this is not a child component -> must be deleted separately (later, when we use a viewport, we don't need this anymore)

  if( waveformBuffer != NULL )
    delete[] waveformBuffer;
  if( waveformPointers != NULL )
    delete[] waveformPointers;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

void OscillatorStereoEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( oscillatorToEdit == NULL )
    return;

  if( buttonThatWasClicked == wavePlusButton )
    AudioFileManager::loadNextFile();
  else if( buttonThatWasClicked == waveMinusButton )
    AudioFileManager::loadPreviousFile();
  else if( buttonThatWasClicked == waveLoadButton )
    AudioFileManager::openLoadingDialog();
  else if( buttonThatWasClicked == moreButton )
  {
    if( moreButton->getToggleState() == true )
    {
      int x = getScreenX() + getWidth();
      int y = getScreenY() - 102; // this is a kludgy thing here with the context menu positioning
      contextMenu->setTopLeftPosition(x, y);
      contextMenu->addToDesktop(
        ComponentPeer::windowHasDropShadow | ComponentPeer::windowIsTemporary);
      contextMenu->setVisible(true);
      contextMenu->toFront(true);

      //contextMenuViewport->setVisible(true);
    }
    else
    {
      contextMenu->setVisible(false);
      //contextMenuViewport->setVisible(false);
    }
    return;
  }
  else if( buttonThatWasClicked == contextMenu->closeButton )
  {
    moreButton->setToggleState(false, true);
    return;
  }


  moduleToEdit->markStateAsDirty();
  //sendChangeMessage();
}

void OscillatorStereoEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  /*
  updateWidgetsAccordingToState();
  sendChangeMessage();
  */
  if( oscillatorToEdit == NULL )
    return;

  if( objectThatHasChanged == contextMenu )
  {
    moduleToEdit->markStateAsDirty();
    updatePlot();
  }
  else
    AudioModuleEditor::changeListenerCallback(objectThatHasChanged);
}

void OscillatorStereoEditor::rSliderValueChanged(RSlider *sliderThatHasChanged)
{
  if( oscillatorToEdit == NULL )
    return;

  moduleToEdit->markStateAsDirty();
  contextMenu->updateWidgetsAccordingToState();

  // although certain sliders would not require an update, most will, so we do the update in any
  // case:
  updatePlot();
}

//-------------------------------------------------------------------------------------------------
//

void OscillatorStereoEditor::updateWidgetsAccordingToState()
{
  if( oscillatorToEdit == NULL )
    return;

  // update the widgets:
  levelSlider->setValue(oscillatorToEdit->getLevel(),                                      false);
  tuneSlider->setValue(oscillatorToEdit->getDetuneSemitones(),                             false);
  pitchModulationSlider->setValue(oscillatorToEdit->getPitchEnvelopeDepth(),               false);

  // update the waveform display plot:
  updatePlot();
  updateWidgetVisibility();

  // update tzhe widgets of the context menu, too:
  contextMenu->updateWidgetsAccordingToState();
}

void OscillatorStereoEditor::mouseDown(const MouseEvent &e)
{
  if( oscillatorToEdit == NULL )
    return;

  MouseEvent e2 = e.getEventRelativeTo(waveformDisplay);
  if( waveformDisplay->contains(Point<int>(e2.x, e2.y)) )
  {
    oscillatorToEdit->setMute( !oscillatorToEdit->isMuted() );
    updateWidgetVisibility();
  }
}

void OscillatorStereoEditor::resized()
{
  Editor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  x = 0;
  w = getWidth()-x;
  y += 16;
  x = getWidth()-20;
  wavePlusButton->setBounds( x, y, 16, 16);
  x -= 14;
  waveMinusButton->setBounds(x, y, 16, 16);
  x = waveMinusButton->getX()-40+2;
  waveLoadButton->setBounds( x, y, 40, 16);

  x = 0;
  y = waveLoadButton->getY();
  w = waveLoadButton->getX();
  h = getHeight() - y;
  waveformDisplay->setBounds(x, y, w-4, h);
  //waveformDisplay->setBounds(x, y, 160, h);  // for debug test
  emptyDisplay->setBounds(waveformDisplay->getBounds());
  waveFileLabel->setBounds(0, waveformDisplay->getY()-16+2, w-4, 16);


  /*
  contextMenuViewport->setBounds(waveformDisplay->getBounds());
  contextMenu->setSize(
  contextMenuViewport->getWidth()-contextMenuViewport->getScrollBarThickness(),
  contextMenu->getHeight());
  */

  x = waveformDisplay->getRight();
  w = getWidth()-x;
  y = waveLoadButton->getBottom()-4;
  levelSlider->setBounds(x+4, y+4, w-8, 32);
  y += 32;
  tuneSlider->setBounds(x+4, y+4, w-8, 32);
  y += 34;
  pitchModulationSlider->setBounds(x+4, y+4, w-8, 16);


  moreButton->setBounds(getWidth()-40-4, getHeight()-16-4, 40, 16);

  // (re-)allocate memory for the waveform buffer and initialize content of the array:
  if( 2*waveformDisplay->getWidth() != numSamplesInPlot )
  {
    numSamplesInPlot = jmax(1, 2*waveformDisplay->getWidth()); // 2x oversampling of the plot to avoid jaggedness
    if( waveformBuffer != NULL )
      delete[] waveformBuffer;
    waveformBuffer = new double[2*numSamplesInPlot];
    for(int n=0; n<2*numSamplesInPlot; n++)
      waveformBuffer[n] = 0.0;
    waveformPointers[0] = &(waveformBuffer[0]);
    waveformPointers[1] = &(waveformBuffer[numSamplesInPlot]);
    if( oscillatorToEdit != NULL )
      oscillatorToEdit->getWaveformForDisplay(waveformPointers, numSamplesInPlot);
    waveformDisplay->setWaveform(waveformPointers, numSamplesInPlot, 2);
  }

  // we must update the plot which has also a new size now:
  updatePlot();
}

/*
void OscillatorStereoEditor::moved()
{
//int x = getScreenX() + getWidth();
//int y = getScreenY();
//contextMenu->setTopLeftPosition(x, y);
contextMenu->setVisible(false);
}
*/

void OscillatorStereoEditor::updatePlot()
{
  if( oscillatorToEdit == NULL )
    return;

  // update the waveform display:
  //int numSamples = oscillatorToEdit->waveTable->getPrototypeNumSamples();
  if( numSamplesInPlot > 0 )
  {
    // retrieve the filename (the full path), cut it down to the filename only and display it in
    // the filename field:
    juce::String fileName = juce::String(oscillatorToEdit->waveTable->getSampleName());
    if( fileName.contains("\\") )
      fileName = fileName.fromLastOccurrenceOf("\\", false, false);
    if( fileName.contains("/") )
      fileName = fileName.fromLastOccurrenceOf("/", false, false);
    waveFileLabel->setText(fileName);


    // ToDo: include a setActiveDirectory method in the FileManager class....
    //AudioFileManager::setActiveDirectory...


    // update the waveform display:
    waveformDisplay->setAutoReRendering(false);
    waveformDisplay->setMaximumRangeX(0.0, numSamplesInPlot);
    waveformDisplay->setCurrentRangeX(0.0, numSamplesInPlot);
    oscillatorToEdit->getWaveformForDisplay(waveformPointers, numSamplesInPlot);
    waveformDisplay->setWaveform(waveformPointers, numSamplesInPlot, 2);
    waveformDisplay->updatePlotImage();
    waveformDisplay->setAutoReRendering(true);
  }

  waveformDisplay->setVisible( !oscillatorToEdit->isMuted() );
  emptyDisplay->setVisible(     oscillatorToEdit->isMuted() );
}

void OscillatorStereoEditor::updateWidgetVisibility()
{
  if( oscillatorToEdit == NULL )
    return;

  waveLoadButton->setVisible(        !oscillatorToEdit->isMuted() );
  wavePlusButton->setVisible(        !oscillatorToEdit->isMuted() );
  waveMinusButton->setVisible(       !oscillatorToEdit->isMuted() );
  moreButton->setVisible(            !oscillatorToEdit->isMuted() );
  levelSlider->setVisible(           !oscillatorToEdit->isMuted() );
  tuneSlider->setVisible(            !oscillatorToEdit->isMuted() );
  pitchModulationSlider->setVisible( !oscillatorToEdit->isMuted() );
  waveFileLabel->setVisible(         !oscillatorToEdit->isMuted() );
  waveformDisplay->setVisible(       !oscillatorToEdit->isMuted() );
  emptyDisplay->setVisible(           oscillatorToEdit->isMuted() );
}

bool OscillatorStereoEditor::setAudioData(AudioSampleBuffer* newBuffer,
  const juce::File& underlyingFile, bool markAsClean)
{
  if( newBuffer != NULL && newBuffer->getNumChannels() > 0 && newBuffer->getNumSamples() > 0 )
  {
    float* channelPointers[2];
    channelPointers[0] = newBuffer->getWritePointer(0, 0);
    if( newBuffer->getNumChannels() >= 2 )
      channelPointers[1] = newBuffer->getWritePointer(1, 0);
    else
      channelPointers[1] = newBuffer->getWritePointer(0, 0);
    oscillatorToEdit->waveTable->setWaveform(channelPointers, newBuffer->getNumSamples());
    juce::String relativePath = underlyingFile.getRelativePathFrom(rootDirectory);
    char* fileNameC           = toZeroTerminatedString(relativePath);
    oscillatorToEdit->waveTable->setSampleName(fileNameC);
    delete[] fileNameC;
    waveFileLabel->setText(underlyingFile.getFileName());
    updatePlot();
    return true;
  }
  return false;
}
