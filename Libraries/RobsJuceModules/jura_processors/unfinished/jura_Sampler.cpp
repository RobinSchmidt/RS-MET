
SamplerModule::SamplerModule(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse) 
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("Sampler");
  setModuleName("Sampler");
}

void SamplerModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef SamplerModule SM;
  typedef Parameter Param;
  Param* p;

  p = new Param("Gain", -48.0, 12.0, 0.1, Parameter::LINEAR); 
  addObservedParameter(p);
  p->setValueChangeCallback<SM>(this, &SM::setGain);
}

AudioModuleEditor* SamplerModule::createEditor(int type)
{
  return new jura::SamplerEditor(this);
}

void SamplerModule::setSampleRate(double newSampleRate)
{ 
  engine.setSampleRate(newSampleRate); 
}

void SamplerModule::setGain(double newGain)
{
  //engine.setGlobalGain(newGain);
}

void SamplerModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  bool markAsClean)
{
  // Recall the global playback parameters such as gain, max num layers, max polyphony, resampling 
  // quality, etc.:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  // The actual instrument definition is loaded from an sfz file that is defined in the xml (just 
  // like sample-files are defined in the xml for the wavetable oscillator):
  juce::String jSfzPath = xmlState.getStringAttribute("InstrumentFile", juce::String());
  std::string  sSfzPath = jSfzPath.toStdString();
  int rc = engine.loadFromSFZ(sSfzPath.c_str());
  if(rc == ReturnCode::fileLoadError)
  {
    //RAPT::rsError("File load error"); // preliminary, for debug
    showWarningBox("SFZ Load Error", "File " + jSfzPath + " could not be loaded.");
    // The .sfz file that was specified in the .xml file was not found
  }
  // todo:
  // -catch other errors, such as sfzParseError, unknownOpcode, sampleLoadError, etc.

  // Notes:
  // -it currently only works, if the .sfz file and its required .wav files reside in the project
  //  directory, i.e.:
  //    RS-MET\Products\AudioPlugins\ToolChain\Builds\VisualStudio2019  
  //  when launching ToolChain from Visual Studio for debugging or in the same directory where 
  //  ToolChain.exe resides, i.e.:
  //    RS-MET\Products\AudioPlugins\ToolChain\Builds\VisualStudio2019\x64\Debug\Standalone Plugin
  //  when launching the .exe directly
  // -when using RAPT::rsError("File load error"); directly starting the exe doesn't work. I 
  //  guess it triggers the error on startup (because it tries to load a file with empty path)
  //  and the debug break makes the app exit immediately?). 

  // ToDo:
  // -Support placing the .sfz and .wav files in different folders. Maybe the path to the .sfz 
  //  file should be specified relatively to the .xml file. The paths to the .wav files should be
  //  either relatively to the .sfz file or relatively to some global sample directory that can be
  //  specified in the xml and/or sfz ...maybe the latter is better but requires us to introduce a 
  //  new opcode to the sfz spec...maybe sample_folder or sample_directory. options should be
  // -retrieve the .sfz and sample folder ...the latter either from the xml or from sfz (not sure
  //  yet, what's best)  
  // -set up the folders in the engine before caling engine.loadFromSFZ
}

XmlElement* SamplerModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  juce::String sfzPath = sfzFile.getRelativePathFrom(getPresetDirectory());
  xmlState->setAttribute("InstrumentFile", sfzPath);
  //xmlState->setAttribute("SampleDirectory", sampleDir);
  return xmlState;
}

void SamplerModule::noteOn(int key, int vel)
{
  Event ev(Event::Type::noteOn, (float)key, (float)vel);
  engine.handleMusicalEvent(ev);
}

void SamplerModule::noteOff(int key)
{
  Event ev(Event::Type::noteOff, (float)key, 0.f);
  engine.handleMusicalEvent(ev);
}

//void SamplerModule::handleMidiMessage(MidiMessage message) 
//{
//
//}

void SamplerModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
    engine.processFrame(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}

void SamplerModule::processStereoFrame(double *dL, double *dR)
{
  engine.processFrame(dL, dR);
}

void SamplerModule::reset()
{ 
  engine.reset();
}

//=================================================================================================

SamplerEditor::SamplerEditor(SamplerModule* samplerToEdit) : AudioModuleEditor(samplerToEdit)
{
  ScopedLock scopedLock(*lock);
  createWidgets();
  setSize(400, 200);
}

void SamplerEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;

  // ...
}

void SamplerEditor::createWidgets()
{
  addWidgetSet(sfzFileLoader = new FileSelectionBox);
  // todo: connect to the FileManager

  addWidget(maxNumLayersSlider = new RDraggableNumber);
  // todo: connect maxNumLayers parameter in module


  addWidget(numLayersLabel = new RTextField);

  addWidget(numLayersField = new RTextField);
}


/*
Bugs:

ToDo:
-add some widgets to the gui that show: 
 -numActiveLayers, numActiveVoices, 
 -loaded sfz file (maybe it should be a load/save widget set)
 -output level meter
 -sliders for MaxLayers, ...
 -menus for: Resampling, SignalFlow (default sfz (setting override), accumulative DSP)
-write some more complex .sfz files with multiple samples, regions, groups, etc. 
-make use of subdirectories for the samples and test, if that works
 -Pluck2.xml loads without errors

-override handleMidiEvent to pass the events to the engine



Maybe the xml presets for the SamplerModule should contain the filename for an sfz file that sould 
be loaded plus some global settings such as the maximum number of layers, the resampling algo, etc. 
Don't use the xml preset to store the actual sfz opcodes. Preset loading would become a nested 
process the outer level is the xml and the inner level is the sfz. The GUI should provide an 
additional set of load/save widgets for the sfz.
*/