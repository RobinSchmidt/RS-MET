
//-------------------------------------------------------------------------------------------------
// construction/destruction:

OscillatorStereoAudioModule::OscillatorStereoAudioModule(CriticalSection *newPlugInLock,     
  rosic::OscillatorStereo *newOscillatorStereoToWrap) : AudioModule(newPlugInLock)
{
  jassert( newOscillatorStereoToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  wrappedOscillatorStereo = newOscillatorStereoToWrap;
  moduleName = juce::String("OscillatorStereo");
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

