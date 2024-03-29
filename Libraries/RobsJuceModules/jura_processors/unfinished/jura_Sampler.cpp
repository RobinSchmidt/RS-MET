
SamplerModule::SamplerModule(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse) 
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("Sampler");
  setModuleName("Sampler");
  setupDirectories();
  createParameters();
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

  p = new Param("BusMode", 0.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<SM>(this, &SM::setBusMode);

  // ToDo: InterpolationMethod (Linear, LagrangeCubic, HermiteCubic, Sinc), SincLength (2-512), 
  // Oversample (1-16), MaxLayers (16-8192), MaxFilters (16-8192), MaxEqualizers, MaxWaveshapers,

}

void SamplerModule::setupDirectories()
{
  ScopedLock scopedLock(*lock);
  juce::String sep = File::getSeparatorString();
  sfzRootDir = jura::getSupportDirectory() + sep + "SFZ" + sep;

  // Sanity check:
  juce::File sfzDirAsFile(sfzRootDir);
  if(!sfzDirAsFile.exists())
    showWarningBox("Error", "SFZ directory: " + sfzRootDir + " does not exist.");
  else if(!sfzDirAsFile.isDirectory())
    showWarningBox("Error", "SFZ directory: " + sfzRootDir + " is not a directory.");

  // Tell the engine, where to find the sfz files:
  bool ok = engine.setSfzRootDir(sfzRootDir.toStdString().c_str());
  if(!ok)
    showWarningBox("Error", "SFZ directory: " + sfzRootDir + " could not be set.");

  // todo: 
  // -assign sampleRootDir, then allow the samples to be located either in the sampleRootDir
  //  *or* as relative path with respect to the .sfz file. Maybe the sampler engine should search 
  //  for the files in the following directories (in  that order):
  //  -subdirectory of the .sfz file
  //  -user defined sample directory
  //  -global default sample directory (defined app-wide)
}

bool SamplerModule::doesSfzFileExist(const juce::String& relativePath)
{
  ScopedLock scopedLock(*lock);
  juce::String path = sfzRootDir + File::getSeparatorString() + relativePath;
  juce::File sfzFile(path);
  return sfzFile.existsAsFile();
}

void SamplerModule::setBusMode(bool shouldAccumulate)
{
  engine.setBusMode(shouldAccumulate);
}

AudioModuleEditor* SamplerModule::createEditor(int type)
{
  return new jura::SamplerEditor(this);
}

void SamplerModule::setSampleRate(double newSampleRate)
{ 
  ScopedLock scopedLock(*lock);
  engine.setSampleRate(newSampleRate); 
}

void SamplerModule::setGain(double newGain)
{
  //engine.setGlobalGain(newGain);
}

void SamplerModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  bool markAsClean)
{
  ScopedLock scopedLock(*lock);

  File appDir = File(getApplicationDirectory());
  appDir.setAsCurrentWorkingDirectory();
  // This is a bit quick-and dirty. On windows, the current working directory is the application
  // directory anyway, unless the app is started from the debugger, in which case it's the 
  // project directory. On the mac, the current working directory is root, i.e. "/", regardless
  // how the app is started. For commandline apps on mac started from the debugger, it's /build
  // inside the project directory. ToDo: find a cleaner way to handle directories. Maybe instead of
  // setting the current working directory, do something like engine.setCurrentSfzRootDirectory()
  // which should set the directory, relative to which the passed paths are interpreted when
  // loadFromSfz is called. hmm...when we set the sfz root director to the xml directory, we may 
  // get porblems when setStateFromXml is called in a total recall situation - in this case, there
  // is no xml file, but we still need some way to know where to look for the sfz file. Maybe the 
  // xml itself should specify the sfz root folder in an attribute? If the sfz is not found, open
  // a dialog that allows the user to locate it. Store the location somewhere and in the next call
  // to getStateAsXml, write it into the xml, such that next time, the updated path is used.

  // Recall the global playback parameters such as gain, max num layers, max polyphony, resampling 
  // quality, etc.:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  // The actual instrument definition is loaded from an sfz file that is defined in the xml (just 
  // like sample-files are defined in the xml for the wavetable oscillator):
  juce::String jSfzPath = xmlState.getStringAttribute("InstrumentFile", juce::String());
  if(jSfzPath.isEmpty()) {
    engine.clearInstrument();
    return; }

  // Check, if such an instrument file exists in the folder where we expect it:
  if(!doesSfzFileExist(jSfzPath))
  {
    showWarningBox("SFZ Load Error", "File " + jSfzPath + " does not exist.");
    // ToDo: Maybe have a warning box with a more comprehensive error message that makes some 
    // suggestions why this could have happened. In this case, it could be that the sfzRootDir
    // does not exist on the machine
  }

  std::string sSfzPath = jSfzPath.toStdString();
  int rc = engine.loadFromSFZ(sSfzPath.c_str());
  if(rc == ReturnCode::fileLoadError)
  {
    //RAPT::rsError("File load error"); // preliminary, for debug
    showWarningBox("SFZ Load Error", "File " + jSfzPath + " could not be loaded.");
    // The .sfz file that was specified in the .xml file was not found
  }
  // todo:
  // -catch other errors, such as sfzParseError, unknownOpcode, sampleLoadError, etc.
  // -figure out, hwat exactly went wrong and show a more specific error message



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
  //  ...i think, this should be fixed now because we now load the sfz files from the proper 
  //  content directory

  // ToDo:
  // -Support placing the .sfz and .wav files in different folders. Maybe the path to the .sfz 
  //  file should be specified relatively to the .xml file. The paths to the .wav files should be
  //  either relatively to the .sfz file or relatively to some global sample directory that can be
  //  specified in the xml and/or sfz ...maybe the latter is better but requires us to introduce a 
  //  new opcode to the sfz spec...maybe sample_folder or sample_directory. options should be
  // -retrieve the .sfz and sample folder ...the latter either from the xml or from sfz (not sure
  //  yet, what's best)  
  // -set up the folders in the engine before caling engine.loadFromSFZ
  // ToDo (maybe): 
  // -extract the directory from the path
  // -set the current working directory to that directory
  // -use as sSfzPath below only the part of the path with the directory stripped off
}

XmlElement* SamplerModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  /*
  // old:
  juce::String sfzPath = sfzFile.getRelativePathFrom(getPresetDirectory());
  xmlState->setAttribute("InstrumentFile", sfzPath);
  //xmlState->setAttribute("SampleDirectory", sampleDir);
  */

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
  {
    float fL, fR;
    engine.processFrame(&fL, &fR);
    inOutBuffer[0][n] = fL;
    inOutBuffer[1][n] = fR;
  }
}

void SamplerModule::processStereoFrame(double *dL, double *dR)
{
  float fL, fR;
  engine.processFrame(&fL, &fR);
  *dL = fL;
  *dR = fR;
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

  layersMeter->setCurrentValue(num);
  //numLayersField->setText(juce::String(num));

  // see: TrackMeterModuleEditor::timerCallback
}

void SamplerEditor::resized()
{
  ScopedLock scopedLock(*lock);  // do we actually nee the lock here?
  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;
  int m = 4;      // margin

  // Status info widgets
  layersMeter->setBounds(x+m, y, w/2-2*m, 16);

  //numLayersLabel->setBounds(x, y, 48, 16);
  //x = numLayersLabel->getRight();
  //numLayersField->setBounds(x, y, 32, 16);

  //x = numLayersField->getRight();
  //numLayersOfLabel->setBounds(x, y, 16, 16);
  //x = numLayersOfLabel->getRight();
  //maxNumLayersSlider->setBounds(x, y, 32, 16);

  // Instrument definition widgets:
  y += 20;
  x  = 2;
  //instrumentLabel->setBounds(x, y, 68, 16);
  //x = instrumentLabel->getRight()+2;
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

  layersMeter = new MeteringDisplayWithText();
  layersMeter->setMeasurementName("Layers");
  layersMeter->setDescription("Number of currently playing layers");
  layersMeter->setDescriptionField(infoField);
  //layersMeter->setMeterStyle(MeteringDisplay::horizontalRatio);  // is the default anyway
  layersMeter->setRange(0, 16);  // should change dynamically on xml load
  //layersMeter->setReferenceValue(0.0);
  //layersMeter->setCurrentValue(0.0);
  addWidget(layersMeter);

  // todo: add an overload warning lamp and a slider to adjust the maximum. maybe the slider should
  // be a Draggable number and doubleas warning lamp - it should become red, when the number was 
  // exceeded

  /*
  addWidget(maxNumLayersSlider = new RDraggableNumber);
  //maxNumLayersSlider->assignParameter( samplerModule->getParameterByName("MaxNumLayers") );
  maxNumLayersSlider->setValue(16);  // preliminary
  maxNumLayersSlider->setDescription("Number of available layers. Drag up/down to adjust.");
  maxNumLayersSlider->setDescriptionField(infoField);
  */

  // ToDo: have display widgets for monitoring the system activity and load: CPU and RAM usage, the 
  // latter for samples and DSP objects seperately, show how many DSP resources are available and 
  // how they are used - like with horizontal meters like:

  //  CPU:         #############         47%           # CPU load
  //  Sample-RAM:  ####               270MB/4GB        # RAM used for samples
  //  Process-RAM: ##                 140MB/4GB        # RAM used for DSPs
  //
  //  Layers:      #####               352/2048        # active region players
  //  Voices:      ######               23/64          # midi active notes
  //
  //  Filters:     #############       467/1024        # used filters
  //  Equalizers:  #######             234/1024        # used equalizers
  //  Amplifiers:
  //  Envelopes:
  //  LowFreqOscs:
  //  etc.
  //
  // Maybe add threads, disk-load, etc. later when we implement such features. Maybe have (sticky) 
  // warning lamps for overloads. maybe have a total-ram display, too? maybe that should inquire
  // the ram-suage of the process from the OS like int the system monitor? maybe we can also 
  // display information about the smoothness/speikeyness of the CPU load?
}


/*
Bugs:
-FilterBlip.xml behaves weird in the low keyrange at low velocities - that's strange because the 
 patch doesn't specify any veltrack stuff. It does have amp-keytrack, though. Also, multiple hits
 seem to get louder with each hit - check if there's maybe some remnant filter state that 
 accumulates over the key-hits
-when switching the preset while holding a note, it crashes (at least on mac)
 -> maybe we need to acquire the lock in setStateFromXml (done -> needs tests)
 -> it's apparently an access violation due to removing samples from the pool in the gui-thread
    while playing them in the audio thread.
-it seems, sometimes the noteOff is not received or handled properly - the layer doesn't 
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
-Not sure if this should be considered a bug: when we run out of RegionPlayers on a noteOn that
 triggers multiple layers, the engine plays as many layers/regions as it can for the given voice.
 Maybe we should not play the voice at all in such a case? Maybe on noteOn, we should make sure 
 that all layers are triggered succesfully or not play the note. actually, in 
 rsSamplerEngine::handleNoteOn there's code that should do the rollback already - mayb it'S bugy?

ToDo:
-On xml-load, update the widgets to reflect the new maximum settings for layers etc.
 -loading an xml may also have to update the max-number of group players - in this case, i think,
  the total number of groups in the patch is a strict upper limit for how many players will ever be
  needed - the saem cannot be said about the relation of regions and layers because a single region
  may produce more than one playback layer, if we allow notes to be retriggered while the old ones
  are still releasing (like in single-shot mode where the sample is always played fully, 
  overlapping with itself - regions may overlap with themselves in playback)
-When an xml file does not specify and sfz file, just load the new engine parameters and keep
 the old sfz file
-Maybe also show a list of the loaded samples. maybe in the left section, a widget to load sfz
 files, belwo it a list of samples showing also the ram usage for them. maybe have also a sort of
 TreeView that represents the SFZ and may allows to edit it

-GUI:
 -show, which .sfz file is loaded and give the user a text editor to edit it
 -show some data about the loaded patch: number of samples, regions, groups, filters, equalizers, 
  waveshapers etc.
 -maybe display also, how many filters, eqs, etc. are allocated and how many are currently in use
 -maybe show a warning, if the number of dsp objects was exceeded...or somehow prevent that from
  ever happening
 -Have 2 views: one for performing and one for editing the patch. the edit view may have largei'sh
  sample editor view (which might even be expanded into its own window). There, we can actually do
  destructive edits on samples with undo/save/etc. But maybe such edits could also be 
  non-destructive? the sampler just stores all actions and when the patch is loaded, re-applies 
  them? ...in this context, we may also develop a GUI for the sinusoidal models.
 
  
 
  
-optimize by using block-wise processing instead of sample-by-sample
-maybe instead of opening a message box on load error, show the error message in the load/save
 widget set 
 -maybe the module should keep some sort of atomic integer "status" variable for that purpose that
  is inquired in the timerCallback
 -to later show a level meter, we perhaps should refactor the editor of TrackMeter such that some 
  intermediate class does not derive from Timer - we don't want the embedded editor have manag its 
  own timer callbacks - that would be redundant - we can do all the required updates here
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
additional set of load/save widgets for the sfz and maybe an sfz editor.

On the GUI, we could also display the number of samples, regions, groups, masters, curves, etc. 
like sfizz: https://github.com/sfztools/sfizz/  ...what are masters and curves? we could also 
display the number of DSP modules...if that makes sense - but i think it doesn't because it depends
on the number of played notes...but maybe we can somehow compute the maximum number that could ever 
be needed? that could also help with allocating the objects...like MaxNumFilters, 
MaxNumWaveShapers, etc. we could also display the number of modulation connections. maybe we could 
do: Regions 3/128, Filters 5/256, etc., RAM usage, CPU usage, output level meters

Another nice project that may help to figure out, how to build and deploy a plugin is Surge:
https://github.com/surge-synthesizer/surge


An editor GUI could be organized as follows:
-Display the instrument on the left as a tree-view. The user can browse the groups and regions.
 -Groups can open/close (+/-) regions
 -Regions can open/close their parameter sets
 -When clicking on a parameter, a slider appears to adjust it
 -the slider can be dragged to a pane to remain in the view persistently
-Display the regions layed out in a rectangle with keyboard at the bottom like in NI Kontakt
 
 
 
 
 
 


*/
