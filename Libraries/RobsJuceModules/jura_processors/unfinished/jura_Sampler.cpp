
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
  samplerModule = samplerToEdit;

  // initialize the current directory for sfz loading:
  //FileManager::setActiveDirectory(getSupportDirectory() + "/SamplerInstruments");

  FileManager::setActiveDirectory(getApplicationDirectory());  // preliminary
  createWidgets();
  setSize(400, 200);
  //startTimer(20);  // in ms
  startTimerHz(50);   // in Hz, i.e. fps (frames per second)
}

bool SamplerEditor::loadFile(const juce::File& fileToLoad)
{
  RAPT::rsError("not yet implemented");
  return false;
}

bool SamplerEditor::saveToFile(const juce::File& fileToSaveTo)
{
  RAPT::rsError("not yet implemented");
  return false;
}

void SamplerEditor::timerCallback()
{
  jassert(samplerModule != nullptr);

  int num = samplerModule->getNumActiveLayers();
  numLayersField->setText(juce::String(num));

  // see: TrackMeterModuleEditor::timerCallback
}

void SamplerEditor::resized()
{
  ScopedLock scopedLock(*lock);  // do we actually nee the lock here?
  AudioModuleEditor::resized();
  int x = 2;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;

  // Status info widgets
  numLayersLabel->setBounds(x, y, 48, 16);
  x = numLayersLabel->getRight()+2;
  numLayersField->setBounds(x, y, 32, 16);
  x = numLayersField->getRight()+2;
  //maxNumLayersSlider->setBounds(x, y, 32, 16);

  // Instrument definition widgets:
  y += 20;
  x  = 2;
  //instrumentLabel->setBounds(x, y, 68, 16);
  x = instrumentLabel->getRight()+2;
  //sfzFileLoader->setBounds(x, y, w-x-4, 16);
}

void SamplerEditor::createWidgets()
{
  addWidget(instrumentLabel = new RTextField);
  instrumentLabel->setText("Instrument:");
  instrumentLabel->setDescription("Currently loaded sfz instrument");
  instrumentLabel->setDescriptionField(infoField);

  //addWidgetSet(sfzFileLoader = new FileSelectionBox("", this) );
  // causes a crash on destruction

  addWidget(numLayersLabel = new RTextField);
  numLayersLabel->setText("Layers:");
  numLayersLabel->setDescription("Number of currently playing layers");
  numLayersLabel->setDescriptionField(infoField);

  addWidget(numLayersField = new RTextField);
  numLayersField->setText("0");
  numLayersField->setDescription("Number of currently playing layers");
  numLayersField->setDescriptionField(infoField);

  addWidget(maxNumLayersSlider = new RDraggableNumber);
  //maxNumLayersSlider->assignParameter( samplerModule->getParameterByName("MaxNumLayers") );
  maxNumLayersSlider->setDescription("Number of available layers. Drag up/down to adjust.");
  maxNumLayersSlider->setDescriptionField(infoField);
}


/*
Bugs:
-it seems, sometimes the noteOff is not received or handeld properly - the layer doesn't 
 immediately stop playing but instead (i think) plays until the end of the sample. It happens when
 triggering at least 5 notes simultaneously 
 -has it to do with the voiceManager in ToolChain - that's plausible, because it's currently set to
  numVoices = maxNumVoices = 4
 -workaround:
  -for the moment, i have called setMaxNumVoices(16); in the constructor of rsVoiceManager (it 
   formerly was called with 4). when trying to call it in the constructor of ToolChain, we hit
   an assert on noteOn
 -fix: 
  -either bypass the voiceManager and treat the sampler as monophonic module (mixing voices 
   internally
  -or somehow route layers belonging to different notes to different voice outputs
  -or make the behavior switchable

ToDo:
-implement better sfz parsing, allowing tokens be separated by any combinations of whitspaces, 
 newlines and (maybe) tabs
-add some widgets to the gui that show: 
 -numActiveLayers, numActiveVoices, cpu-load, memory occupation
 -loaded sfz file (maybe it should be a load/save widget set)
 -output level meter
 -sliders for MaxLayers, ...
 -menus for: Resampling, SignalFlow (default sfz (setting override), accumulative DSP), mode (DFD, 
  RAM 32 Bit, RAM 16 Bit)
 -have an editor for the sfz. 
  -left is a sort of tree-view browser showing the groups and regions
  -for the selected region, show lokey/hiKey, lovel/hivel and a text editor with the opcodes
   ...or maybe some sort of "spreadsheet" like editor - 1 line per opcode, opcodes are sorted by 
   type, righ-click - insert with a menu of available opcodes
   -maybe show the region's lokey/hikey lovel/hivel settings as rectangle in an k*(128x128) square
    (k >= 1 to it a reasonable size)...or maybe 512x256 or something
  -maybe show a keyboard that shows the currently held notes (maybe the releasing one should be 
   shown also, but in fainter color - maybe the color could fade with the release time)
-write some more complex .sfz files with multiple samples, regions, groups, etc. 
-make use of subdirectories for the samples and test, if that works
 -Pluck2.xml loads without errors

-maybe implement multitimbrality:
 -allow multiple sfz files to be loaded simultaneously that respond to different midi channels
 -have one engine for each channel
 -let the engines use a shared sample pool (that will complicate the delete process from the sample
  pool because then, we should not always delete samples that are not used in some sfz - we must 
  check all loaded sfzs)



Maybe the xml presets for the SamplerModule should contain the filename for an sfz file that sould 
be loaded plus some global settings such as the maximum number of layers, the resampling algo, etc. 
Don't use the xml preset to store the actual sfz opcodes. Preset loading would become a nested 
process the outer level is the xml and the inner level is the sfz. The GUI should provide an 
additional set of load/save widgets for the sfz.
*/