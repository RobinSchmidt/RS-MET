
//-------------------------------------------------------------------------------------------------
// construction/destruction:

WaveOscModule::WaveOscModule(CriticalSection *newPlugInLock,
  rosic::OscillatorStereo *oscToWrap) : AudioModuleWithMidiIn(newPlugInLock)
{
  if(oscToWrap == nullptr) {
    wrappedOsc = new rosic::OscillatorStereo;
    waveTable  = new rosic::MipMappedWaveTableStereo;
    wrappedOsc->setWaveTableToUse(waveTable);
    wrappedOscIsOwned = true;
  }
  else
    wrappedOsc = oscToWrap;
  //setModuleTypeName("OscillatorStereo"); // old name
  setModuleTypeName("WaveOscillator");
  createParameters();
  loadDefaultWaveform();
}

WaveOscModule::~WaveOscModule()
{
  if(wrappedOscIsOwned) {
    delete wrappedOsc;
    delete waveTable;
  }
}

AudioModuleEditor* WaveOscModule::createEditor(int type)
{
  return new WaveOscEditor(lock, this); // get rid of passing the lock
}

//-------------------------------------------------------------------------------------------------
// state saving and recall:


/*
// These two functions are legacy code and shoudl go away - still using them is actually a bug:
XmlElement* oscillatorStereoStateToXml(OscillatorStereo* osc, XmlElement* xmlElementToStartFrom)
{
  // the XmlElement which stores all the releveant state-information:
  XmlElement* xmlState;
  if(xmlElementToStartFrom == nullptr)
    xmlState = new XmlElement("WaveOscillator");
    //xmlState = new XmlElement(juce::String("OscillatorStereo"));  //wrong? should be OscillatorStereoState
  else
    xmlState = xmlElementToStartFrom;

  juce::String samplePath = juce::String(osc->waveTable->getSampleName());
  if(!samplePath.isEmpty())   // new
    xmlState->setAttribute("AudioFileRelativePath", samplePath);
  return xmlState;
}
bool oscillatorStereoStateFromXml(OscillatorStereo* osc, const XmlElement &xmlState)
{
  // this gets called multiple times on startup - why?

  bool success = false;

  // Temporarily switch off the automatic re-rendering of the mip-map, to avoid multiple renderings
  // (one for each parameter)...but that should probably be done in WaveOscModule::getStateAsXml.
  // Doing it here worked in the legacy code:
  bool oldAutoReRenderState = osc->waveTable->isMipMapAutoReRenderingActive();
  osc->waveTable->setAutomaticMipMapReRendering(false);
  osc->waveTable->fillWithAllZeros();

  // let the oscillator and all its slaves calculate their increment:
  //osc->setFrequency(1000.0);
  osc->calculateIncrementForAllSlaves();

  // load the audio-file into the wavetable for the oscillator:
  juce::String relativePath = xmlState.getStringAttribute("AudioFileRelativePath", juce::String());
  juce::String absolutePath = getSupportDirectory() + File::getSeparatorString() + relativePath;

  // old presets were stored with backslashes - these don't work on mac unless we replace 
  // the backslashes by forward slashes:
  relativePath = relativePath.replaceCharacter('\\', '/');
  absolutePath = absolutePath.replaceCharacter('\\', '/');


  //juce::File audioFile(absolutePath);
  //if( !audioFile.existsAsFile() )
  //  return false;
  // ...Nope! That bypasses the warning message box triggered in 
  // AudioFileManager::createAudioSampleBufferFromFile when the file couldn't be loaded. But we 
  // want to see that warning! It's important to alert the user when waveform loading goes wrong.
  // Just silently skipping it is not acceptable.


  AudioSampleBuffer* buffer = AudioFileManager::createAudioSampleBufferFromFile(absolutePath, true);
  if( buffer != nullptr )
  {
    // pass the actual audio data:
    float* channelPointers[2];
    channelPointers[0] = buffer->getWritePointer(0, 0);
    if( buffer->getNumChannels() >= 2 )
      channelPointers[1] = buffer->getWritePointer(1, 0);
    else
      channelPointers[1] = buffer->getWritePointer(0, 0);
    osc->waveTable->setWaveform(channelPointers, buffer->getNumSamples() );
    delete buffer;
    success = true;
  }
  else
  {
    // ToDo: maybe init waveform with silence - call a function like:
    //osc->waveTable->initWaveform();
    // which would have to be written. It makes more sense to use silence than to keep whatever
    // waveform is currently loaded. That would create a false impression of the state.
    // ...what about fillWithAllZeros that is called above
  }

  // pass the path as c-string:
  char* fileNameC = toZeroTerminatedString(relativePath);
  //jassert(fileNameC[0] != 'C'); // for debug
  osc->waveTable->setSampleName(fileNameC);
  delete[] fileNameC;

  // let the (wavetable inside the) oscillator render the mip map and restore the old state of the
  // automatic re-rendering:
  osc->waveTable->renderMipMap();
  osc->waveTable->setAutomaticMipMapReRendering(oldAutoReRenderState);

  return success;
}
*/



XmlElement* WaveOscModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  // New: 
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  xmlState->setAttribute("AudioFileRelativePath", samplePathRelative);
  return xmlState;


  /*
  // Old:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  if( wrappedOsc != nullptr ) // that should actually never be the case
    xmlState = oscillatorStereoStateToXml(wrappedOsc, xmlState);
  return xmlState;
  */
}

void WaveOscModule::setStateFromXml(
  const XmlElement& xmlState, const juce::String& stateName, bool markAsClean)
{
  // New:
  jassert(wrappedOsc != nullptr && wrappedOsc->waveTable != nullptr);
  if( wrappedOsc != nullptr && wrappedOsc->waveTable != nullptr )
    wrappedOsc->waveTable->setAutomaticMipMapReRendering(false);
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  juce::String samplePath = xmlState.getStringAttribute("AudioFileRelativePath", "");
  loadWaveform(samplePath);
  if( wrappedOsc != nullptr && wrappedOsc->waveTable != nullptr )
    wrappedOsc->waveTable->setAutomaticMipMapReRendering(true);

  // Note:
  // The reason for switching the automatic re-rendering of the mip-map in the embedded 
  // wrappedOsc->waveTable temporarily off is that because otherwise our call to 
  // AudioModule::setStateFromXml would trigger multiple re-renderings - one for each recalled 
  // parameter that affects the mip-map. I think so, at least ...verify!
  //
  // ToDo:
  // Maybe if the "AudioFileRelativePath" attribute is not found, use some default file that
  // always exists like Silence.flac. That would avoid the popping up of a warning box that says 
  // that the audio file [with empty name] wasn't found. In this case, it's acceptable to go 
  // without warning because the xml did not pecify a sample file anyway. Or maybe have a special
  // function clearWaveform or setToEmptyWaveform or something - it shouldn't even need a file 
  // with silence for this


  /*
  // Old:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  if( wrappedOsc != nullptr )
    oscillatorStereoStateFromXml(wrappedOsc, xmlState);
    */
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void WaveOscModule::createParameters()
{
  typedef MetaControlledParameter Param;
  Param* p;

  typedef Parameter StaticParam;  // for-non-automatable parameters
  StaticParam* sp; 

  typedef rosic::OscillatorStereo OS;
  OS* os = wrappedOsc;

  typedef rosic::MipMappedWaveTableStereo WT;
  WT* wt = os->waveTable;

  std::vector<double> defaultValues;

  // amplitude related parameters:

  p = new Param("Mute", 0.0, 1.0, 0.0, Parameter::BOOLEAN, 1.0);
  p->setValueChangeCallback<OS>(os, &OS::setMute);
  addObservedParameter(p);

  p = new Param("Level",    -36.0, 12.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<OS>(os, &OS::setLevel);
  defaultValues.push_back(-3.01029996);   // compensation gain for  2 uncorrelated sources
  defaultValues.push_back(-4.77121255);   // compensation gain for  3 uncorrelated sources
  defaultValues.push_back(-6.02059991);   // compensation gain for  4 uncorrelated or 2 in-phase sources
  defaultValues.push_back(-9.03089987);   // compensation gain for  8 uncorrelated
  defaultValues.push_back(-12.0411998);   // compensation gain for 16 uncorrelated or 4 in-phase sources
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("LevelByKey", -24.0, 24.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<OS>(os, &OS::setLevelByKey);
  addObservedParameter(p);

  p = new Param("LevelByVel", -12.0, 12.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<OS>(os, &OS::setLevelByVel);
  addObservedParameter(p);

  p = new Param("MidSide",    0.0,  1.0, 0.5, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<OS>(os, &OS::setMidSide);
  defaultValues.clear();
  defaultValues.push_back(0.0);
  defaultValues.push_back(0.25);
  defaultValues.push_back(0.5);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("Pan",       -1.0,  1.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<OS>(os, &OS::setPan);
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

  p = new Param("Tune",     -36.0, 36.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<OS>(os, &OS::setDetuneSemitones);
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

  p = new Param("DetuneHz", -20.0, 20.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<OS>(os, &OS::setDetuneHz);
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

  p = new Param("StereoDetune", -1.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<OS>(os, &OS::setStereoDetuneSemitones);
  defaultValues.clear();
  defaultValues.push_back(-0.2);
  defaultValues.push_back(-0.1);
  defaultValues.push_back(0.0);
  defaultValues.push_back(0.1);
  defaultValues.push_back(0.2);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("StereoDetuneHz", -10.0, 10.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<OS>(os, &OS::setStereoDetuneHz);
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

  p = new Param("PitchModulationDepth", -8.0, 8.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<OS>(os, &OS::setPitchEnvelopeDepth);
  defaultValues.clear();
  defaultValues.push_back(-1.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(1.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // time domain waveform related parameters:

  p = new Param("StartPhase", 0.0, 360.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<OS>(os, &OS::setStartPhase);
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

  sp = new StaticParam("FullWaveWarp", -0.99, 0.99, 0.0, Parameter::LINEAR, 0.001); 
  sp->setValueChangeCallback<OS>(os, &OS::setFullWavePhaseWarp);
  addObservedParameter(sp);

  sp = new StaticParam("HalfWaveWarp", -0.99, 0.99, 0.0, Parameter::LINEAR, 0.001);
  sp->setValueChangeCallback<OS>(os, &OS::setHalfWavePhaseWarp);
  addObservedParameter(sp);

  sp = new StaticParam("TimeReverse", 0.0, 1.0, 0.0, Parameter::BOOLEAN, 1.0);
  sp->setValueChangeCallback<OS>(os, &OS::setTimeReverse);
  addObservedParameter(sp);

  sp = new StaticParam("PolarityInvert", 0.0, 1.0, 0.0, Parameter::BOOLEAN, 1.0);
  sp->setValueChangeCallback<OS>(os, &OS::setPolarityInversion);
  addObservedParameter(sp);


  // #013:
  /*
  sp = new StaticParam("CombOffset", 0.0, 360.0, 0.0, Parameter::LINEAR, 0.1);
  defaultValues.clear();
  defaultValues.push_back(0.0);
  defaultValues.push_back(45.0);
  defaultValues.push_back(90.0);
  defaultValues.push_back(135.0);
  defaultValues.push_back(180.0);
  defaultValues.push_back(225.0);
  defaultValues.push_back(270.0);
  defaultValues.push_back(315.0);
  sp->setDefaultValues(defaultValues);
  addObservedParameter(sp);
  */

  sp = new StaticParam("CombHarmonic", 1.0, 128.0, 1.0, Parameter::EXPONENTIAL, 0.01);
  sp->setValueChangeCallback<WT>(wt, &WT::setCombHarmonic);
  addObservedParameter(sp);

  sp = new StaticParam("CombAmount", -100.0, 100.0, 0.0, Parameter::LINEAR, 0.1);
  sp->setValueChangeCallback<WT>(wt, &WT::setCombAmount);
  addObservedParameter(sp);

  // magnitude spectrum related parameters:

  sp = new StaticParam("SpectralContrast", 0.25, 4.0, 1.0, Parameter::EXPONENTIAL, 0.01);
  sp->setValueChangeCallback<WT>(wt, &WT::setSpectralContrast);
  addObservedParameter(sp);

  sp = new StaticParam("SpectralSlope", -6.0, 6.0, 0.0, Parameter::LINEAR, 0.01);
  sp->setValueChangeCallback<WT>(wt, &WT::setSpectralSlope);
  addObservedParameter(sp);

  sp = new StaticParam("HighestHarmonic", 1.0, 1024.0, 1024.0, Parameter::EXPONENTIAL, 1.0);
  sp->setValueChangeCallback<WT>(wt, &WT::setHighestHarmonicToKeep);
  addObservedParameter(sp);

  sp = new StaticParam("LowestHarmonic", 1.0, 1024.0, 1.0, Parameter::EXPONENTIAL, 1.0);
  sp->setValueChangeCallback<WT>(wt, &WT::setLowestHarmonicToKeep);
  addObservedParameter(sp);

  sp = new StaticParam("EvenOddRatio", 0.0, 1.0, 0.5, Parameter::LINEAR, 0.005);
  sp->setValueChangeCallback<WT>(wt, &WT::setEvenOddRatio);
  addObservedParameter(sp);

  //-----------------------------------------------------------------------------------------------
  // phase spectrum related parameters:

  sp = new StaticParam("PhaseScale", -1.0, 1.0, 1.0, Parameter::LINEAR, 0.01);
  sp->setValueChangeCallback<WT>(wt, &WT::setPhaseScale);
  sp->setDefaultValues(defaultValues); // ?
  addObservedParameter(sp);

  sp = new StaticParam("PhaseShift", -180.0, 180.0, 0.0, Parameter::LINEAR, 1.0);
  sp->setValueChangeCallback<WT>(wt, &WT::setPhaseShift);
  sp->setDefaultValues(defaultValues);
  addObservedParameter(sp);

  sp = new StaticParam("EvenOddPhaseShift", -180.0, 180.0, 0.0, Parameter::LINEAR, 1.0);
  sp->setValueChangeCallback<WT>(wt, &WT::setEvenOddPhaseShift);
  defaultValues.clear();
  defaultValues.push_back(-90.0);
  defaultValues.push_back(-45.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(45.0);
  defaultValues.push_back(90.0);
  sp->setDefaultValues(defaultValues);
  addObservedParameter(sp);

  sp = new StaticParam("StereoPhaseShift", -180.0, 180.0, 0.0, Parameter::LINEAR, 1.0);
  sp->setValueChangeCallback<WT>(wt, &WT::setStereoPhaseShift);
  sp->setDefaultValues(defaultValues);
  addObservedParameter(sp);

  sp = new StaticParam("EvenOddStereoPhaseShift", -180.0, 180.0, 0.0, Parameter::LINEAR, 1.0);
  sp->setValueChangeCallback<WT>(wt, &WT::setEvenOddStereoPhaseShift);
  sp->setDefaultValues(defaultValues);
  addObservedParameter(sp);
}

bool WaveOscModule::setWaveform(AudioSampleBuffer* buffer, const juce::File& waveFile)
{
  jassert(wrappedOsc != nullptr && wrappedOsc->waveTable != nullptr);
  if(wrappedOsc == nullptr || wrappedOsc->waveTable == nullptr)
    return false;

  // Set up the waveform in the DSP core from the given buffer or initialize as empty in case of 
  // nullptr:
  if( buffer != nullptr )
  {
    juce::String relPath = waveFile.getRelativePathFrom(getSupportDirectory());
    float* channelPointers[2];
    channelPointers[0] = buffer->getWritePointer(0, 0);
    if( buffer->getNumChannels() >= 2 )
      channelPointers[1] = buffer->getWritePointer(1, 0);
    else
      channelPointers[1] = buffer->getWritePointer(0, 0);
    wrappedOsc->waveTable->setWaveform(channelPointers, buffer->getNumSamples());
    wrappedOsc->waveTable->renderMipMap(); // Shouldn't this happen automatically?
    wrappedOsc->waveTable->setSampleName(relPath.toStdString().c_str()); // Redundant...
    samplePathRelative = relPath;                                        // ...this line should be enough
    return true;
  }
  else
  {
    wrappedOsc->waveTable->fillWithAllZeros();
    wrappedOsc->waveTable->setSampleName("");  // Redundant...
    samplePathRelative = "";                   // ...this line should be enough
    return false;
  }

  // ToDo: 
  // -Store the sample name only directly here. The core DSP class should not be concerned with 
  //  this. The fact that we store the sample name there has only historical reasons.
  // -What if the buffer has the wrong number of channels? can this happen?
  // -Get rid of storing the sample path redundantly in wrappedOsc
  // -Figure out if we need the call to renderMipMap
  // -Why do we acquire write-pointers in buffer->getWritePointer. Maybe because setWaveform isn't
  //  const-correct? -> make it so and acquire red pointers.
}

bool WaveOscModule::loadWaveform(const String& relativePath)
{
  //jassert(wrappedOsc != nullptr && wrappedOsc->waveTable != nullptr);
  //if(wrappedOsc == nullptr || wrappedOsc->waveTable == nullptr)
  //  return false;

  // Old presets were stored with backslashes - these don't work on mac unless we replace 
  // the backslashes by forward slashes:
  juce::String relPath = relativePath.replaceCharacter('\\', '/');
  juce::String absPath = getSupportDirectory() + File::getSeparatorString() + relPath;
  juce::File waveFile(absPath);
  AudioSampleBuffer* buffer = AudioFileManager::createAudioSampleBufferFromFile(waveFile, true);
  bool success = setWaveform(buffer, waveFile);
  delete buffer;
  return success;

  //return setWaveform(buffer, waveFile);

  // ToDo:
  // -Try to avoid usage of raw pointer and manual deletion for the buffer. Let 
  //  AudioFileManager::createAudioSampleBufferFromFile return some sort of smart pointer.
}

void WaveOscModule::loadDefaultWaveform()
{
  //loadWaveform(wrappedOsc, "/Samples/SingleCycle/Classic/Saw.flac");
  loadWaveform("/Samples/SingleCycle/Classic/Sine.flac");
  // I think, the saw is too aggressive to be fired up immediately when the osc is plugged in. The
  // sine is better as default (= initial) waveform. Maybe we should even reduce the volume 
  // initially, i.e. set the Level parameter to soething like -10 initially (but leave the deafult
  // at zero nonetheless).
}

//=================================================================================================

// construction/destruction:

WaveOscEditorContextMenu::WaveOscEditorContextMenu(
  WaveOscModule* newOscillatorModuleToEdit, Component* componentToAttachTo)
  //: ComponentMovementWatcher(componentToAttachTo)
{
  jassert(newOscillatorModuleToEdit != nullptr);
  oscillatorModuleToEdit = newOscillatorModuleToEdit;
  createWidgets();
  //updateWidgetsAccordingToState();
  setSize(180, 488);
}

WaveOscEditorContextMenu::~WaveOscEditorContextMenu()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void WaveOscEditorContextMenu::rButtonClicked(RButton *buttonThatWasClicked)
{
  /* not needed anymore
  if( oscillatorModuleToEdit == NULL )
    return;
  if( oscillatorModuleToEdit->wrappedOscillatorStereo == NULL )
    return;

  rosic::OscillatorStereo* o = oscillatorModuleToEdit->wrappedOscillatorStereo;

  if( buttonThatWasClicked == reverseButton )
    o->setTimeReverse( reverseButton->getToggleState() );
  if( buttonThatWasClicked == invertButton )
    o->setPolarityInversion( invertButton->getToggleState() );
    */

  sendChangeMessage();
}

void WaveOscEditorContextMenu::rSliderValueChanged(RSlider *sliderThatHasChanged)
{
  /*
  if( oscillatorModuleToEdit == NULL )
    return;
  if( oscillatorModuleToEdit->wrappedOscillatorStereo == NULL )
    return;
  //rosic::OscillatorStereo* o = oscillatorModuleToEdit->wrappedOscillatorStereo;
  */
  sendChangeMessage();
}
// ToDo: Document why we send change messages. Who picks them up? What gets updated in response?

void WaveOscEditorContextMenu::createWidgets()
{
  typedef rsAutomatableSlider Sld;  Sld* s;
  typedef rsAutomatableButton Btn;  Btn* b;
  typedef RTextField        Txf;  Txf* t;

  addWidget( ampHeadline = t = new Txf("Amplitude:") );
  t->setNoBackgroundAndOutline(true);
  t->setDescription("Manipulations of the amplitude");

  addWidget( tuningHeadline = t = new Txf("Tuning:") );
  t->setNoBackgroundAndOutline(true);
  t->setDescription("Manipulations of the tuning/detuning of the oscillator");

  addWidget( timeHeadline = t = new Txf("Time:") );
  t->setNoBackgroundAndOutline(true);
  t->setDescription("Time domain manipulations of the waveform");

  addWidget( magSpectrumHeadline = t = new Txf("Magnitude Spectrum:") );
  t->setNoBackgroundAndOutline(true);
  t->setDescription("Manipulations of the magnitude spectrum");

  addWidget( phaseSpectrumHeadline = t = new Txf("Phase Spectrum:") );
  t->setNoBackgroundAndOutline(true);
  t->setDescription("Manipulations of the phase spectrum");

  // sliders for amplitude related parameters:

  addWidget( levelSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("Level") );
  s->setDescription("Output level of the oscillator");
  s->setStringConversionFunction(&decibelsToStringWithUnit2);
  s->addListener(this); // to send out the change-message for display update

  addWidget( levelByKeySlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("LevelByKey") );
  s->setSliderName("Key");
  s->setDescription("Key dependence of oscillator's output level");
  s->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( levelByVelSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("LevelByVel") );
  s->setSliderName("Vel");
  s->setDescription("Velocity dependence of oscillator's output level");
  s->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( midSideSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("MidSide") );
  s->setSliderName("Mid/Side");
  s->setDescription("Mid/side adjustment for stereo(ized) waveforms");
  s->setStringConversionFunction(&ratioToString0);
  s->addListener(this);

  addWidget( panSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("Pan") );
  s->setSliderName("Pan");
  s->setDescription("Panorama position of the oscillator");
  s->setStringConversionFunction(&valueToString2);
  s->addListener(this);

  // sliders for tuning related parameters:

  addWidget( tuneSlider = new TuningSlider("TuneSlider") );
  tuneSlider->assignParameter( oscillatorModuleToEdit->getParameterByName("Tune") );
  tuneSlider->setDescription("Tuning of the oscillator in semitones");
  tuneSlider->setStringConversionFunction(&semitonesToStringWithUnit2);

  addWidget( detuneHzSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("DetuneHz") );
  s->setSliderName("Detune Hz");
  s->setDescription("Detuning of the oscillator in Hz");
  s->setStringConversionFunction(&hertzToStringWithUnit2);

  addWidget( stereoDetuneSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("StereoDetune") );
  s->setSliderName("Stereo Detune");
  s->setDescription("Detuning between left and right channel in semitones");
  s->setStringConversionFunction(&semitonesToStringWithUnit2);

  addWidget( stereoDetuneHzSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("StereoDetuneHz") );
  s->setSliderName("Stereo Detune Hz");
  s->setDescription("Detuning between left and right channel in Hz");
  s->setStringConversionFunction(&hertzToStringWithUnit2);

  // sliders for time domain related parameters:

  addWidget( startPhaseSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("StartPhase") );
  s->setSliderName("Start Phase");
  s->setDescription("Start phase of the oscillator");
  s->setStringConversionFunction(&degreesToStringWithUnit0);
  s->addListener(this);

  addWidget( fullWavePhaseWarpSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("FullWaveWarp") );
  s->setSliderName("Full Wave Warp");
  s->setDescription("Applies phase warping to the entire waveform");
  s->setStringConversionFunction(&degreesToStringWithUnit0);
  s->addListener(this);

  addWidget( halfWavePhaseWarpSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("HalfWaveWarp") );
  s->setSliderName("Half Wave Warp");
  s->setDescription("Applies phase warping both half cycles of the waveform");
  s->setStringConversionFunction(&degreesToStringWithUnit0);
  s->addListener(this);

  addWidget( combHarmonicSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("CombHarmonic") );
  s->setSliderName("Comb Harmonic");
  s->setDescription("Harmonic on which the comb filter acts");
  //s->setStringConversionFunction(&degreesToStringWithUnit0);
  s->setStringConversionFunction(&valueToString2);
  s->addListener(this);
  //s->setVisible(false); // not yet meaningfully implemented

  addWidget( combAmountSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("CombAmount") );
  s->setSliderName("Comb Amount");
  s->setDescription("Amount of comb filtering");
  s->setStringConversionFunction(percentToStringWithUnit1);
  s->addListener(this);
  //s->setVisible(false); // not yet meaningfully implemented

  addWidget( reverseButton = b = new Btn("Reverse") );
  b->assignParameter( oscillatorModuleToEdit->getParameterByName("TimeReverse") );
  b->addRButtonListener(this);
  b->setDescription("Time reverses the oscillator's waveform");
  b->setClickingTogglesState(true);

  addWidget( invertButton = b = new Btn("Invert") );
  b->assignParameter( oscillatorModuleToEdit->getParameterByName("PolarityInvert") );
  b->addRButtonListener(this);
  b->setDescription("Inverts polarity of the oscillator's ouput");
  b->setClickingTogglesState(true);

  // sliders for magnitude spectrum related parameters:

  addWidget( spectralContrastSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("SpectralContrast") );
  s->setSliderName("Contrast");
  s->setDescription("Spectral contrast acting as exponent on the harmonic's magnitude");
  s->setStringConversionFunction(&valueToString2);
  s->addListener(this);

  addWidget( spectralSlopeSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("SpectralSlope") );
  s->setSliderName("Slope");
  s->setDescription("Spectral slope applied to the waveform in dB/oct");
  s->setStringConversionFunction(&decibelsPerOctaveToString2);
  s->addListener(this);

  addWidget( highestHarmonicSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("HighestHarmonic") );
  s->setSliderName("Highest Harmonic");
  s->setDescription("Highest harmonic in the waveform");
  s->setStringConversionFunction(&valueToString0);
  s->addListener(this);

  addWidget( lowestHarmonicSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("LowestHarmonic") );
  s->setSliderName("Lowest Harmonic");
  s->setDescription("Lowest harmonic in the waveform");
  s->setStringConversionFunction(&valueToString0);
  s->addListener(this);

  addWidget( evenOddSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("EvenOddRatio") );
  s->setSliderName("Even/Odd");
  s->setDescription("Ratio of even and odd harmonics");
  s->setStringConversionFunction(&ratioBothFullAtCenterToString0);
  s->addListener(this);

  // sliders for phase spectrum related parameters:

  addWidget( evenOddPhaseShiftSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("EvenOddPhaseShift") );
  s->setSliderName("Even/Odd Shift:");
  s->setDescription("Applies a phase shift between even and odd harmonics");
  s->setStringConversionFunction(&degreesToStringWithUnit0);
  s->addListener(this);

  addWidget( phaseScaleSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("PhaseScale") );
  s->setSliderName("Scale:");
  s->setDescription("Scales the phase of each harmonic ");
  //s->setStringConversionFunction(&degreesToStringWithUnit0);
  s->setStringConversionFunction(&valueToString2);
  s->addListener(this);

  addWidget( phaseShiftSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("PhaseShift") );
  s->setSliderName("Shift:");
  s->setDescription("Shifts the phase of each harmonic by a constant");
  s->setStringConversionFunction(&degreesToStringWithUnit0);
  s->addListener(this);

  addWidget( stereoPhaseShiftSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("StereoPhaseShift") );
  s->setSliderName("Stereo Shift:");
  s->setDescription("Applies stereoization via phase-shifting of harmonics");
  s->setStringConversionFunction(&degreesToStringWithUnit0);
  s->addListener(this);

  addWidget( evenOddStereoPhaseShiftSlider = s = new Sld );
  s->assignParameter( oscillatorModuleToEdit->getParameterByName("EvenOddStereoPhaseShift") );
  //s->setSliderName("Even/Odd Stereo Shift:");
  s->setSliderName("Ev/Od Ster Shft:");
  s->setDescription("Phase shift between even/odd harmonics, applied with opposite signs to left/right channels");
  s->setStringConversionFunction(&degreesToStringWithUnit0);
  s->addListener(this);

  addWidget( closeButton = new RButton(RButton::CLOSE) );
  closeButton->setDescription("Closes the oscillator context menu");
  closeButton->setClickingTogglesState(false);
  // we don't listen to this button ourselves - this is the job of the outlying editor object
}

/*
// may be obsolete - verify and delete
void WaveOscEditorContextMenu::updateWidgetsAccordingToState()
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
*/

void WaveOscEditorContextMenu::resized()
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

  ampHeadline->setBounds(x+4, y+4, w-8-16, sh);     y += inc;
  levelSlider->setBounds(x+4, y+4, w-8, sh);        y += sh-2;
  levelByKeySlider->setBounds(x+4,  y+4, w2-8, sh);
  levelByVelSlider->setBounds(w2+4, y+4, w2-8, sh); y += inc;
  midSideSlider->setBounds(x+4, y+4, w-8, sh);      y += sh-2;
  panSlider->setBounds(x+4, y+4, w-8, sh);          y += inc;

  tuningHeadline->setBounds(x+4, y+4, w-8, sh);          y += inc;
  tuneSlider->setBounds(x+4, y+4, w-8, sh);              y += sh-2;
  detuneHzSlider->setBounds(x+4, y+4, w-8, sh);          y += sh-2;
  stereoDetuneSlider->setBounds(x+4, y+4, w-8, sh);      y += sh-2;
  stereoDetuneHzSlider->setBounds(x+4, y+4, w-8, sh);    y += inc;
  timeHeadline->setBounds(x+4, y+4, w-8, sh);            y += inc;
  startPhaseSlider->setBounds(x+4, y+4, w-8, sh);        y += sh-2;
  combHarmonicSlider->setBounds(x+4, y+4, w-8, sh);      y += sh-2;
  combAmountSlider->setBounds(x+4, y+4, w-8, sh);        y += sh-2;
  fullWavePhaseWarpSlider->setBounds(x+4, y+4, w-8, sh); y += sh-2;
  halfWavePhaseWarpSlider->setBounds(x+4, y+4, w-8, sh); y += inc;
  reverseButton->setBounds(x+4,  y+4, w2-8, sh);
  invertButton->setBounds( w2+4, y+4, w2-8, sh);         y += inc;

  magSpectrumHeadline->setBounds(x+4, y+4, w-8, sh);     y += inc;
  spectralContrastSlider->setBounds(x+4, y+4, w-8, sh);  y += sh-2;
  spectralSlopeSlider->setBounds(x+4, y+4, w-8, sh);     y += sh-2;
  highestHarmonicSlider->setBounds(x+4, y+4, w-8, sh);   y += sh-2;
  lowestHarmonicSlider->setBounds(x+4, y+4, w-8, sh);    y += sh-2;
  evenOddSlider->setBounds(x+4, y+4, w-8, sh);           y += inc;
  phaseSpectrumHeadline->setBounds(x+4, y+4, w-8, sh);   y += inc;
  phaseScaleSlider->setBounds(x+4, y+4, w-8, sh);        y += sh-2;
  phaseShiftSlider->setBounds(x+4, y+4, w-8, sh);        y += sh-2;
  evenOddPhaseShiftSlider->setBounds(x+4, y+4, w-8, sh); y += sh-2;
  stereoPhaseShiftSlider->setBounds(x+4, y+4, w-8, sh);  y += sh-2;
  evenOddStereoPhaseShiftSlider->setBounds(x+4, y+4, w-8, sh);
}

/*
void WaveOscEditorContextMenu::componentMovedOrResized(bool wasMoved, bool wasResized)
{
int dummy = 0;
}

void WaveOscEditorContextMenu::componentPeerChanged()
{
int dummy = 0;
}
*/

//=================================================================================================
// class WaveOscEditor:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WaveOscEditor::WaveOscEditor(CriticalSection *newPlugInLock,
  WaveOscModule* newWaveOscModule) : SampleBasedAudioModuleEditor(newWaveOscModule)
{
  // init the pointer to the modulator to be edited to NULL:
  jassert( newWaveOscModule != nullptr );
  oscModule = newWaveOscModule;

  createWidgets();

  // Overwrite the descriptions of the inherited sample-load widgets:
  sampleFileLabel->setDescription("Name of the currently loaded waveform");
  sampleLoadButton->setDescription("Load a waveform");
  samplePlusButton->setDescription("Next waveform in current directory");
  sampleMinusButton->setDescription("Previous waveform in current directory");

  // Maybe create this only when needed and init to nullptr as in lazy initialization:
  contextMenu = new WaveOscEditorContextMenu(newWaveOscModule, this);
  contextMenu->addChangeListener(this);
  contextMenu->setAlwaysOnTop(true);
  contextMenu->setOpaque(true);
  contextMenu->closeButton->addRButtonListener(this);
  addChildColourSchemeComponent(contextMenu, false, false);
  //contextMenu->setSize(200, 200);

  numSamplesInPlot    = 0;
  waveformBuffer      = nullptr;
  waveformPointers    = new double*[2];
  waveformPointers[0] = nullptr;
  waveformPointers[1] = nullptr;

  isTopLevelEditor = false;
  setHeadlineStyle(NO_HEADLINE);

  // Initialize the current directory for waveform loading:
  AudioFileManager::setActiveDirectory(getSupportDirectory() + "/Samples/SingleCycle/Classic");




  updateWidgetsAccordingToState();

  setSize(232, 136);  // Widget arrangement is optimized for this size
}

WaveOscEditor::~WaveOscEditor()
{
  delete contextMenu; 
  // This is not a child component -> must be deleted separately (later, when we use a viewport,
  // we don't need this anymore). In the constructor, we actually call 
  //   addChildColourSchemeComponent(contextMenu, false, false);
  // but the first "false" parameter says that it should not be added as child component. But then
  // who is the parent? ...it's actually a nullptr.

  if( waveformBuffer != nullptr )
    delete[] waveformBuffer;
  if( waveformPointers != nullptr )
    delete[] waveformPointers;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

void WaveOscEditor::rButtonClicked(RButton *b)
{
  if( oscModule == nullptr )
    return;

  if( b == moreButton )
  {
    if( moreButton->getToggleState() == true )
    {
      int x = getScreenX() + getWidth();
      int y = getScreenY() - 102; // this is a kludgy thing here with the context menu positioning

      // Move this code into a function WaveOscEditorContextMenu::showAt(x, y):

      contextMenu->setTopLeftPosition(x, y);
      contextMenu->addToDesktop(
        ComponentPeer::windowHasDropShadow | ComponentPeer::windowIsTemporary);

      // Experimental since 2024/03/01 to try to fix the problem that the menu is behind the main 
      // GUI when ToolChain is used as plugin:
      /*
      juce::ComponentPeer* contextPeer = contextMenu->getPeer();
      if(contextPeer)
      {
        contextPeer->toFront(false);  // false: do not take keyboard focus
        // ...nope - this also doesn't help. 
        // There are also methods like ComponentPeer::setAlwaysOnTop etc. - Check them out, too.
      }
      */
      // We don't have such a problem in standalone mode. I think, the problem might be related to
      // the fact that plugin GUIs are themselves opened in "always on top" mode by the host 
      // whereas in standalone mode this is not the case? Maybe check the Surge code, how they 
      // handle their context menus. Or maybe look into the code of juce::PopupMenu
      // juce::PopupMenu::showWithOptionalCallback - there, it creates a pointer to a window that
      // is never deleted - that looks like a memory leak to me. This problem also does not occur 
      // in Bitwig - which is good.

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
  else if( b == contextMenu->closeButton )
  {
    moreButton->setToggleState(false, true);
    return;
  }
  else
    SampleBasedAudioModuleEditor::rButtonClicked(b); // it was one of the inherited buttons


  moduleToEdit->markStateAsDirty();
  //sendChangeMessage();
}

void WaveOscEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  /*
  updateWidgetsAccordingToState();
  sendChangeMessage();
  */
  if( oscModule == nullptr )
    return;

  if( objectThatHasChanged == contextMenu )
  {
    moduleToEdit->markStateAsDirty();
    updatePlot();
  }
  else
    AudioModuleEditor::changeListenerCallback(objectThatHasChanged);
}

void WaveOscEditor::rSliderValueChanged(RSlider *sliderThatHasChanged)
{
  if( oscModule == nullptr )
    return;

  moduleToEdit->markStateAsDirty();
  contextMenu->updateWidgetsAccordingToState();

  // although certain sliders would not require an update, most will, so we do the update in any
  // case:
  updatePlot();
}

//-------------------------------------------------------------------------------------------------
//

void WaveOscEditor::updateWidgetsAccordingToState()
{
  // Check, if this is needed - if not, remove:
  if( oscModule == nullptr || oscModule->wrappedOsc == nullptr )
    return;


  // Update the waveform display plot:
  updatePlot();
  updateWidgetVisibility();

  // Update the widgets of the context menu, too:
  contextMenu->updateWidgetsAccordingToState();
}

void WaveOscEditor::mouseDown(const MouseEvent &e)
{
  // Check, if this is needed - if not, remove:
  if( oscModule == nullptr || oscModule->wrappedOsc == nullptr )
    return;

  if(containsPoint(waveformDisplay, e.x, e.y))
  {
    // Old, buggy:
    //oscModule->wrappedOsc->setMute( !oscModule->wrappedOsc->isMuted() );
    // !!!BUG!!! Directly accessing the DSP core object bypasses "Mute" parameter.

    // New:
    jura::Parameter* p = oscModule->getParameterByName("Mute");
    jassert(p);
    bool muted = p->getValue();
    p->setValue(!muted, true, true);


    updateWidgetVisibility();
  }
}

void WaveOscEditor::resized()
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
  samplePlusButton->setBounds( x, y, 16, 16); x -= 14;
  sampleMinusButton->setBounds(x, y, 16, 16); x = sampleMinusButton->getX()-40+2;
  sampleLoadButton->setBounds( x, y, 40, 16);

  x = 0;
  y = sampleLoadButton->getY();
  w = sampleLoadButton->getX();
  h = getHeight() - y;
  waveformDisplay->setBounds(x, y, w-4, h);
  //waveformDisplay->setBounds(x, y, 160, h);  // for debug test
  emptyDisplay->setBounds(waveformDisplay->getBounds());
  sampleFileLabel->setBounds(0, waveformDisplay->getY()-16+2, w-4, 16);

  /*
  contextMenuViewport->setBounds(waveformDisplay->getBounds());
  contextMenu->setSize(
  contextMenuViewport->getWidth()-contextMenuViewport->getScrollBarThickness(),
  contextMenu->getHeight());
  */

  x = waveformDisplay->getRight();
  w = getWidth()-x;
  y = sampleLoadButton->getBottom()-4;
  levelSlider->setBounds(x+4, y+4, w-8, 32); y += 32;
  tuneSlider->setBounds(x+4, y+4, w-8, 32);  y += 34;
  pitchModulationSlider->setBounds(x+4, y+4, w-8, 16);

  moreButton->setBounds(getWidth()-40-4, getHeight()-16-4, 40, 16);

  // (re-)allocate memory for the waveform buffer and initialize content of the array:
  if( 2*waveformDisplay->getWidth() != numSamplesInPlot )
  {
    numSamplesInPlot = jmax(1, 2*waveformDisplay->getWidth()); // 2x oversampling of the plot to avoid jaggedness
    if( waveformBuffer != nullptr )
      delete[] waveformBuffer;
    waveformBuffer = new double[2*numSamplesInPlot];
    for(int n=0; n<2*numSamplesInPlot; n++)
      waveformBuffer[n] = 0.0;
    waveformPointers[0] = &(waveformBuffer[0]);
    waveformPointers[1] = &(waveformBuffer[numSamplesInPlot]);
    if(oscModule != nullptr && oscModule->wrappedOsc != nullptr)
    {
      oscModule->wrappedOsc->getWaveformForDisplay(waveformPointers, numSamplesInPlot);
    }
    waveformDisplay->setWaveform(waveformPointers, numSamplesInPlot, 2);
  }

  // we must update the plot which has also a new size now:
  updatePlot();
}

/*
void WaveOscEditor::moved()
{
//int x = getScreenX() + getWidth();
//int y = getScreenY();
//contextMenu->setTopLeftPosition(x, y);
contextMenu->setVisible(false);
}
*/

void WaveOscEditor::createWidgets()
{
  addPlot( waveformDisplay = new rsWaveformPlot() );
  waveformDisplay->setAutoReRendering(false);
  waveformDisplay->setDescription(juce::String("Click on the display to switch oscillator on/off"));
  waveformDisplay->setAxisPositionX(rsPlotSettings::INVISIBLE);
  waveformDisplay->setAxisPositionY(rsPlotSettings::INVISIBLE);
  waveformDisplay->setCurrentRangeY(-1.2, 1.2);
  waveformDisplay->setSampleRate(1.0);
  waveformDisplay->addMouseListener(this, true);
  waveformDisplay->setAutoReRendering(true);

  addPlot( emptyDisplay = new rsPlot() );
  emptyDisplay->setAutoReRendering(false);
  emptyDisplay->setDescription(waveformDisplay->getDescription());
  emptyDisplay->setAxisPositionX(rsPlotSettings::INVISIBLE);
  emptyDisplay->setAxisPositionY(rsPlotSettings::INVISIBLE);
  emptyDisplay->setCaption(juce::String("Off"), rsPlotSettings::CENTER);
  emptyDisplay->addMouseListener(this, true);
  emptyDisplay->setAutoReRendering(true);

  /*
  contextMenuViewport = new Viewport(juce::String(T("OscillatorStereoSliderViewport")));
  contextMenuViewport->setViewedComponent(contextMenu);
  contextMenuViewport->setScrollBarThickness(12);
  contextMenuViewport->setScrollBarsShown(true, false);
  addAndMakeVisible(contextMenuViewport);
  */

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
}

void WaveOscEditor::updatePlot()
{
  if( oscModule == nullptr || oscModule->wrappedOsc == nullptr )
    return;

  // update the waveform display:
  //int numSamples = oscillatorToEdit->waveTable->getPrototypeNumSamples();
  if( numSamplesInPlot > 0 )
  {
    // retrieve the filename (the full path), cut it down to the filename only and display it in
    // the filename field:
    juce::String fileName = juce::String(oscModule->wrappedOsc->waveTable->getSampleName());
    if( fileName.contains("\\") )
      fileName = fileName.fromLastOccurrenceOf("\\", false, false);
    if( fileName.contains("/") )
      fileName = fileName.fromLastOccurrenceOf("/", false, false);
    sampleFileLabel->setText(fileName);

    // ToDo: include a setActiveDirectory method in the FileManager class....
    //AudioFileManager::setActiveDirectory...

    // update the waveform display:
    waveformDisplay->setAutoReRendering(false);
    waveformDisplay->setMaximumRangeX(0.0, numSamplesInPlot);
    waveformDisplay->setCurrentRangeX(0.0, numSamplesInPlot);
    oscModule->wrappedOsc->getWaveformForDisplay(waveformPointers, numSamplesInPlot);
    waveformDisplay->setWaveform(waveformPointers, numSamplesInPlot, 2);
    waveformDisplay->updatePlotImage();
    waveformDisplay->setAutoReRendering(true);
  }

  waveformDisplay->setVisible( !oscModule->wrappedOsc->isMuted() );
  emptyDisplay->setVisible(     oscModule->wrappedOsc->isMuted() );
}

void WaveOscEditor::updateWidgetVisibility()
{
  if( oscModule == nullptr || oscModule->wrappedOsc == nullptr )
    return;
  bool muted = oscModule->wrappedOsc->isMuted();
  sampleLoadButton->setVisible(     !muted);
  samplePlusButton->setVisible(     !muted);
  sampleMinusButton->setVisible(    !muted);
  sampleFileLabel->setVisible(      !muted);
  moreButton->setVisible(           !muted);
  levelSlider->setVisible(          !muted);
  tuneSlider->setVisible(           !muted);
  pitchModulationSlider->setVisible(!muted);
  waveformDisplay->setVisible(      !muted);
  emptyDisplay->setVisible(          muted);
}

bool WaveOscEditor::setAudioData(AudioSampleBuffer* newBuffer,
  const juce::File& underlyingFile, bool markAsClean)
{
  // New - causes access violations when clicking on the "load Next" button:
  jassert(oscModule != nullptr);
  if(oscModule == nullptr)
    return false;
  bool success = oscModule->setWaveform(newBuffer, underlyingFile);
  updatePlot();
  return success;


  /*
  // old:
  if( oscModule == nullptr || oscModule->wrappedOsc == nullptr )
    return false;

  // This needs to be rewritten

  if( newBuffer != nullptr && newBuffer->getNumChannels() > 0 && newBuffer->getNumSamples() > 0 )
  {
    float* channelPointers[2];
    channelPointers[0] = newBuffer->getWritePointer(0, 0);
    if( newBuffer->getNumChannels() >= 2 )
      channelPointers[1] = newBuffer->getWritePointer(1, 0);
    else
      channelPointers[1] = newBuffer->getWritePointer(0, 0);
    oscModule->wrappedOsc->waveTable->setWaveform(channelPointers, newBuffer->getNumSamples());

    //juce::String relativePath = underlyingFile.getRelativePathFrom(rootDirectory);
    // this seems wrong - we must use the support-folder as root directory

    juce::String relativePath = underlyingFile.getRelativePathFrom(getSupportDirectory());

    //juce::String root = rootDirectory;
    //juce::String support = getSupportDirectory();

    char* fileNameC = toZeroTerminatedString(relativePath);
    oscModule->wrappedOsc->waveTable->setSampleName(fileNameC);
    delete[] fileNameC;
    sampleFileLabel->setText(underlyingFile.getFileName());
    updatePlot();
    return true;
  }
  return false;
  */
}


/*

ToDo:

-Maybe have a cutoff slider that applies to the slope: below the cutoff, the specrum is left alone,
 above it, the slope kicks in - or vice versa depending on whether the slope is upward or downward
-Maybe have a second slope and cutoff slider such that we may generate bandpass spectra
-Allow steeper slopes - that makes sense in the context of having different oscs generate different
 bandpass spectra.

*/