SfzPlayer::SfzPlayer()
{
  setupDirectories();
  wildcardPatterns = String("*.sfz");  // To show only .sfz files in the dialog box?
  defaultExtension = String(".sfz");   // To append this as default extension when saving?
  updateFileList();
}

bool SfzPlayer::loadFile(const juce::File& fileToLoad)
{
  bool ok = loadFile(fileToLoad.getRelativePathFrom(sfzRootDir));             // This may fail,
  if(!ok) {                                                       // ...and in case of failure, 
    bool ok2 = setupFromSfzString(lastValidSfz);   // ...we revert to the last known valid sfz
    jassert(ok2); }
    // This should actually never fail because lastValidSfz is supposed to be a valid sfz. But what
    // if the user has deleted some required sample files in the meantime? Maybe in that case, we
    // should reset lastValidSfz and reset the engine, too

  return ok;
}

bool SfzPlayer::saveToFile(const juce::File& fileToSaveTo)
{
  showWarningBox("Error", "SfzPlayer::saveToFile not yet implemented");
  return false;
  // ToDo: we somehow need to connect this to the string shown in the code-editor. Maybe we need to
  // keep (a copy of) that string here in this class as member, maybe a tmpSfz juce::String
  // I think, we can use the lastValid Sfz member for this
}

bool SfzPlayer::loadFile(const juce::String& relativePath)
{
  //juce::String path = sfzRootDir + File::getSeparatorString() + relativePath;

  // Check, if such an sfz file exists in the folder where we expect it:
  juce::String path = sfzRootDir + relativePath;
  juce::File sfzFile(path);
  if(!sfzFile.existsAsFile())
  {
    showWarningBox("SFZ Load Error", "File " + path + " does not exist.");
    return false; // The .sfz file that was specified in the relativePath file was not found
    // ToDo: Maybe have a warning box with a more comprehensive error message that makes some 
    // suggestions why this could have happened. In this case, it could be that the sfzRootDir
    // does not exist on the machine
  }
  // ToDo: I think, this error handling here is redundant with the error handling in the call
  // Engine::loadFromSFZ below - verify that and if so, maybe get rid of it here. But before 
  // getting rid, we should make sure that the subsequent error reporting code is a bit more 
  // specific than it currently is. File not existent is a different condition than the general
  // "unable to load file" condition which may have lots of other reasons.

  // OK, the file exists. Let the engine try to load it:
  std::string sSfzPath = relativePath.toStdString();
  int rc = Engine::loadFromSFZ(sSfzPath.c_str());
  if(rc == ReturnCode::fileLoadError)
  {
    showWarningBox("SFZ Load Error", "File " + path + " could not be loaded.");
    return false;
  }
  // todo:
  // -catch other errors, such as sfzParseError, unknownOpcode, sampleLoadError, etc.
  // -figure out, what exactly went wrong and show a more specific error message

  // OK, the engine did apparently load the new sfz file successully. Now we need to store the
  // content of the file here in a member here too:
  lastValidSfz = sfzFile.loadFileAsString();
  return true;

  // Old Notes (obsolete):
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
}

bool SfzPlayer::setupFromSfzString(const juce::String& newSfz)
{
  bool ok = Engine::setFromSFZ(newSfz.toStdString());
  if(!ok) {
    bool ok2 = setupFromSfzString(lastValidSfz); 
    jassert(ok2);        // This should never happen, see comments in SfzPlayer::loadFile.
    return false; }
  lastValidSfz = newSfz;
  return true;
}

void SfzPlayer::setupDirectories()
{
  juce::String sep = File::getSeparatorString();
  sfzRootDir = jura::getSupportDirectory() + sep + "SFZ" + sep;
  FileManager::setActiveDirectory(jura::getSupportDirectory() + "/SFZ");  // use sfzRootDir
  // Maybe using SamplerPatches is more aesthetically pleasing than SFZ because of consistency 
  // with the existing folder structure

  // Sanity check:
  juce::File sfzDirAsFile(sfzRootDir);
  if(!sfzDirAsFile.exists())
    showWarningBox("Error", "SFZ directory: " + sfzRootDir + " does not exist.");
  else if(!sfzDirAsFile.isDirectory())
    showWarningBox("Error", "SFZ directory: " + sfzRootDir + " is not a directory.");

  // Tell the engine, where to find the sfz files:
  bool ok = Engine::setSfzRootDir(sfzRootDir.toStdString().c_str());
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

//=================================================================================================

SamplerModule::SamplerModule(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse) 
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("Sampler");
  setModuleName("Sampler");
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
  juce::String sfzPath = xmlState.getStringAttribute("InstrumentFile", juce::String());
  if(sfzPath.isEmpty())  
    engine.clearInstrument();    // no sfz file was specified in the xml
  else
    engine.loadFile(sfzPath);
}

XmlElement* SamplerModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  /*
  // old:
  juce::String sfzPath = sfzFile.getRelativePathFrom(getPresetDirectory());
  xmlState->setAttribute("InstrumentFile", sfzPath);
  //xmlState->setAttribute("SampleDirectory", sampleDir);
  // ToDo: Maybe Retrieve the filename of the currently loaded .sfz file and store it in the xml
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

SamplerEditor::SamplerEditor(SamplerModule* samplerToEdit) 
  : AudioModuleEditor(samplerToEdit)
  , samplerModule(samplerToEdit)
  , sfzEditor(sfzDoc, nullptr)
  //, sfzFileManager(samplerToEdit)
{
  ScopedLock scopedLock(*lock);
  //samplerModule = samplerToEdit;

  createWidgets();
  setSize(400, 200);
  //startTimer(20);  // in ms
  startTimerHz(50);   // in Hz, i.e. fps (frames per second) for the metering widgets
  sfzDoc.addListener(this);
  parseButton->addRButtonListener(this);

  samplerModule->engine.addFileManagerListener(this); 
  // We want to receive activeFileChanged callbacks when the currently active sfz file changes.
}

SamplerEditor::~SamplerEditor()
{
  samplerModule->engine.removeFileManagerListener(this);
}

void SamplerEditor::timerCallback()
{
  jassert(samplerModule != nullptr);

  int num = samplerModule->getNumActiveLayers();

  layersMeter->setCurrentValue((float)num);
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

  sfzFileLoader->setBounds(x, y, 300, 16);
  y += 16;
  sfzEditor.setBounds(x, y, w-8, getHeight()-y);  // preliminary, uses almost the full available space
  y -= 16;
  w  = 48;
  x  = sfzEditor.getRight() - w;
  parseButton->setBounds(x, y, w, 16);
}

void SamplerEditor::codeDocumentTextInserted(const String& newText, int insertIndex)
{
  setCodeIsDirty();
}

void SamplerEditor::codeDocumentTextDeleted(int startIndex, int endIndex)
{
  setCodeIsDirty();
}

void SamplerEditor::rButtonClicked(RButton* b)
{
  if(b == parseButton)
    parseCurrentEditorContent();
  else
    AudioModuleEditor::rButtonClicked(b);
}

void SamplerEditor::activeFileChanged(FileManager* fileMan)
{

  // There are actually two FileManager objects that we could listen to: the SamplerModule object
  // which manages the preset .xml file of the whole sampler and the SfzPlayer object embedded in
  // the former which manages the currently loaded .sfz instrument definition. We currently listen
  // only to the latter. Changes in the former are of no interest because a preset change due to 
  // the user loading another xml will usually imply also a change of the sfz file, i.e. loading
  // a new xml will trigger loading a new sfz. Which sfz should be loaded is defined in the xml 
  // preset. The xml-presets are on a higher level than the sfz instruments. Think of the loaded
  // sfz instrument more like a waveform sample that is loaded into an oscillator. We treat it in
  // a similar way.
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


  // The SFZ editor:
  sfzFileLoader = new jura::FileSelectionBox("", &samplerModule->engine);
  addWidgetSet(sfzFileLoader);
  sfzFileLoader->setDescription("Current SFZ file");

  addWidget(parseButton = new jura::RClickButton("Parse"));  // Maybe use a right-arrow ("Play")
  parseButton->setDescription("Parse the current sfz document");

  addAndMakeVisible(sfzEditor);
  // ToDo: set up the description of the editor...but it's not a subclass of RWidget...hmmm..

}

void SamplerEditor::setCodeIsParsed(bool isParsed)
{
  codeIsParsed = isParsed;
  //parseButton->setEnabled(!codeIsParsed); // doesn't seem to have any effect
  parseButton->setVisible(!codeIsParsed);   // works but is a bit drastic
}

void SamplerEditor::setCodeIsSaved(bool isSaved)
{
  // ToDo: make an atserisk appear in the sfz file box, if is not saved, disappea, if it is saved
}

void SamplerEditor::parseCurrentEditorContent()
{
  bool ok = samplerModule->setupFromSfzString(sfzEditor.getDocument().getAllContent());
  setCodeIsParsed(ok);

  // -try to parse the code
  // -setCodeIsParsed(true) ...maybe only in case of success?
}

void SamplerEditor::saveCurrentEditorContent()
{

}

//=================================================================================================

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
 -Make sfz editor functional:
  -Show the sfz content in the editor after loading an xml or sfz file
  -We may need to dirtify the xml-state when the user loads a new sfz
  -We need to implement saving of the sfz.
  -Maybe we should have an update button that tries to parse the current editor content. If it 
   fails, it would be nice to know where...give the user some sort of error indicator in the file. 
   Maybe have some sort of sfzValidator class.
 -show some data about the loaded patch: number of samples, regions, groups, filters, equalizers, 
  waveshapers etc.
 -maybe display also, how many filters, eqs, etc. are allocated and how many are currently in use
 -In the left section, show a tree-view of the instrument with levels <global>, <group>, <region>
  where <global> may be implicit, i.e. the top-level without an opening "folder" widget
  -The leaves are the individual opcodes for the pamaters. When the user clicks on one, it gets 
   selected an widgets for adjusting it appear in the main view (maybe it the middle)
  -If the selected leaf is a "sample" opcode, show a little waveform view. Maybe it could allow the
   use to load a new sample.
  -Maybe mutliple leaves can be selected simultaneously.
  -Maybe the selection should be stored in the .xml file, such that the right widgets can be made
   to appear on preset load
  -Maybe there could be some section with performance parameters, maybe use MIDI-controllers for 
   that. What about polyphonic and non-destructive modulation, i.e. the new features of CLAP? We 
   should definitely consider to support them in some way.
 -maybe show a warning, if the number of dsp objects was exceeded...or somehow prevent that from
  ever happening
 -Have 2 views: one for performing and one for editing the patch. the edit view may have largei'sh
  sample editor view (which might even be expanded into its own window). There, we can actually do
  destructive edits on samples with undo/save/etc. But maybe such edits could also be 
  non-destructive? the sampler just stores all actions and when the patch is loaded, re-applies 
  them? ...in this context, we may also develop a GUI for the sinusoidal models.
 -Maybe we should use a Mediator pattern similar to the Liberty GUI to coordinate TreeView and
  CodeEditor.
 
 
  
 
-Allow optional voice retriggering and voice stealing. At the moment, Each note starts a new voice.
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
