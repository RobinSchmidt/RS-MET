SfzPlayer::SfzPlayer()
{
  setupDirectories();
  wildcardPatterns = String("*.sfz");  // To show only .sfz files in the dialog box?
  defaultExtension = String(".sfz");   // To append this as default extension when saving?
  updateFileList();

  // ToDo: maybe support this opcode:
  //   https://sfzformat.com/opcodes/default_path
  // ...but that actually be handled by the sfz-engine
}

bool SfzPlayer::loadFile(const juce::File& fileToLoad)
{
  if(!fileToLoad.existsAsFile()) {
    showWarningBox("SFZ Load Error", "File " + fileToLoad.getFullPathName() + " does not exist.");
    return false; }
    // ToDo: Maybe have a warning box with a more comprehensive error message that makes some 
    // suggestions why this could have happened. In this case, it could be that the sfzRootDir
    // does not exist on the machine 
  else {
    juce::String sfzString = fileToLoad.loadFileAsString();
    bool ok = setupFromSfzString(sfzString, true);
    if(ok)
      FileManager::setActiveFile(fileToLoad);  // makes the sfz-file box reflect the new file
    return ok; }
}

bool SfzPlayer::saveToFile(const juce::File& fileToSaveTo)
{
  juce::Result res = fileToSaveTo.create();
  if(!res.wasOk()) {
    showWarningBox("Error", "File could not be created");
    return false; }
  else {
    bool ok = fileToSaveTo.replaceWithText(lastValidSfz);
    if(!ok) {
      showWarningBox("Error", "File could not be written");
      return false; }
    return true; }
}

bool SfzPlayer::loadFile(const juce::String& relativePath)
{
  return loadFile(juce::File(sfzRootDir + relativePath));
}

bool SfzPlayer::setupFromSfzString(const juce::String& newSfz, bool stringComesFromFile)
{
  bool ok = Engine::setFromSFZ(newSfz.toStdString());  // This may fail due to malformed sfz, etc.
  if(!ok) {
    // In case of failure to set up the engine from the new sfz string, we restore the old patch.
    // Restoring should never fail because lastValidSfz is always supposed to be a valid sfz. But 
    // what if the user has deleted some required sample files in the meantime? Maybe in that case,
    // we should reset lastValidSfz and reset the engine, too?
    bool restored = Engine::setFromSFZ(lastValidSfz.toStdString());
    jassert(restored);
    return false; }
  else {
    lastValidSfz = newSfz; // If all went well, the newSfz becomes the lastValidSfz for next time
    markFileAsClean(stringComesFromFile); // Controls "dirty" asterisk in sfz file-widget
    return true; }
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
  sfzPlayer.addFileManagerListener(this); 
}

void SamplerModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  // I'm not yet sure, if it makes sense to have the sampler engine have parameters in the same 
  // sense that the other AudioModule subclasses have. Maybe we need none of these...


  /*
  // Hmm maybe that's not the right way to go about it:
  // Create a fixed number of sfzParameter objects that can be dynamically connected to SFZ 
  // opcodes.
  using SM = SamplerModule;
  using OP = SfzOpcodeParameter;
  for(int i = 0; i < numOpcodeParams; i++)
  {
    OP* p =  new OP("OpcodeParam" + juce::String(i), 0.0, 1.0, 0.01, jura::Parameter::LINEAR);
    addObservedParameter(p);              // This call will register ouselves as observer for p
  }
  // Or: Maybe create the parameters dynamically - create always one parameter for each opcode and
  // show the sliders in the TreeView...maybe make a TreeView, whose leafs are the sliders.
  */


  /*
  // Code currently not used.

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
  */

}

void SamplerModule::setBusMode(bool shouldAccumulate)
{
  sfzPlayer.setBusMode(shouldAccumulate);
}

AudioModuleEditor* SamplerModule::createEditor(int type)
{
  return new jura::SamplerEditor(this);
}

void SamplerModule::setSampleRate(double newSampleRate)
{ 
  ScopedLock scopedLock(*lock);
  sfzPlayer.setSampleRate(newSampleRate); 
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
    sfzPlayer.clearInstrument();    // no sfz file was specified in the xml
  else
    sfzPlayer.loadFile(sfzPath);
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

void SamplerModule::activeFileChanged(FileManager* fileMan)
{
  jassert(fileMan == &sfzPlayer);
  markStateAsDirty();
  // The .xml preset file may not reflect the currently loaded .sfz file anymore. It's like when 
  // you load a new waveform for an oscillator.
}

void SamplerModule::parameterChanged(Parameter* param)
{
  // ToDo:
  // -Try to dyncmically cast param into an SfzOpcodeParamter
  // -If successful:
  //  -Update the corresponding setting in the instrument definition. To do this, we need:
  //   -Figure out where the setting is in the instrument defintion datastructure. This is given by
  //    the groupIndex/regionIndex/opcode data
  //   -Call sfzPlayer->setRegionSetting or setGroupSetting or setInstrumentSetting depending on 
  //    whether its a global, group or instrument setting.
  // -If unsuccessful, fall back to baseclass implementation

  int dummy = 0;
}

void SamplerModule::noteOn(int key, int vel)
{
  Event ev(Event::Type::noteOn, (float)key, (float)vel);
  sfzPlayer.handleMusicalEvent(ev);
}

void SamplerModule::noteOff(int key)
{
  Event ev(Event::Type::noteOff, (float)key, 0.f);
  sfzPlayer.handleMusicalEvent(ev);
}

void SamplerModule::setMidiController(int idx, float val)
{
  Event ev(Event::Type::controlChange, (float)idx, val);
  sfzPlayer.handleMusicalEvent(ev);
  midiCtrlDirty = true;  // ToDo: later in C++20, use std::atomic_flag::test()

  // ToDo:
  // Maybe have a controllersDirty field in the class and set it to true whenenever a controller is
  // received. Receiving a controller may bring the state of the engine out of sync with what is 
  // defined in the sfz code via the set_ccN opcodes. Maybe when we want to save a patch, we don't
  // want to save it with the old, out-of-sync values. Maybe we should respond to CC changes in a
  // similar way as to changes of low-level params in the TreeView: update the code
}

/*
void SamplerModule::handleMidiMessage(MidiMessage message) 
{
  int dummy = 0;
}
*/

void SamplerModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
  {
    float fL, fR;
    sfzPlayer.processFrame(&fL, &fR);
    inOutBuffer[0][n] = fL;
    inOutBuffer[1][n] = fR;
  }
}

void SamplerModule::processStereoFrame(double *dL, double *dR)
{
  float fL, fR;
  sfzPlayer.processFrame(&fL, &fR);
  *dL = fL;
  *dR = fR;
}

void SamplerModule::reset()
{ 
  sfzPlayer.reset();
}

//=================================================================================================

void SamplerInterfaceComponent::handleMediatorNotification(MediatedColleague* originator, 
  int messageCode, rsMessageData* messageData)
{
  PatchChangeInfo *info = dynamic_cast<PatchChangeInfo*>(messageData);
  if(info != nullptr)
    handlePatchUpdate(*info);
  else
    jassertfalse; // messageData must be of type jura::PatchChangeInfo
}

//=================================================================================================

SfzNodeData::SfzNodeData()
{

}  

SfzNodeData::SfzNodeData(const SfzNodeData& other)
{
  type = other.type;
  data = other.data;
  groupIndex  = other.groupIndex;
  regionIndex = other.regionIndex;
}

SfzNodeData SfzNodeData::createEmptyNode()
{
  SfzNodeData node;
  return node;
}

SfzNodeData SfzNodeData::createGroupNode(int groupIndex)
{
  SfzNodeData node;
  node.type = Type::group;
  node.groupIndex = groupIndex;
  return node;
}

SfzNodeData SfzNodeData::createRegionNode(int groupIndex, int regionIndex)
{
  SfzNodeData node;
  node.type = Type::region;
  node.groupIndex = groupIndex;
  node.regionIndex = regionIndex;
  return node;
}

SfzNodeData SfzNodeData::createPlaybackSettingNode(int groupIndex, int regionIndex,
  const rosic::Sampler::PlaybackSetting& playbackSetting)
{
  SfzNodeData node;
  node.type = Type::playbackSetting;
  node.groupIndex = groupIndex;
  node.regionIndex = regionIndex;
  node.data.playbackSetting = playbackSetting;
  return node;
}

SfzNodeData SfzNodeData::createModulationRoutingNode(int groupIndex, int regionIndex,
  const rosic::Sampler::ModulationRouting& modulationRouting)
{
  SfzNodeData node;
  node.type = Type::modulationRouting;
  node.groupIndex = groupIndex;
  node.regionIndex = regionIndex;
  node.data.modRouting = modulationRouting;
  return node;
}

bool SfzNodeData::isOpcodeNode() const
{
  return type == Type::modulationRouting || type == Type::playbackSetting;
}

SfzNodeData::OpcodeFormat SfzNodeData::getOpcodeFormat() const
{
  using namespace rosic::Sampler;
  if(type == Type::modulationRouting) {
    return OpcodeFormat::Float; }
  if(type == Type::playbackSetting) {
    Opcode op = data.playbackSetting.getOpcode();
    SfzCodeBook* cb = SfzCodeBook::getInstance();
    return cb->getOpcodeFormat(op); }
  return OpcodeFormat::Unknown;
}

//=================================================================================================

SfzOpcodeWidgetSet::SfzOpcodeWidgetSet()
{
  createWidgets();
}

bool SfzOpcodeWidgetSet::wantsExponentialSlider(rosic::Sampler::Opcode op) const
{
  return false;
  // Preliminary. ToDo: either apply a heuristic based on the range of the opcode's values or 
  // decide it in a case-statement based on the opcode. Maybe later we also want scalings other 
  // than linear or exponential?

  using namespace rosic::Sampler;
  SfzCodeBook* cb = SfzCodeBook::getInstance();
  float minVal = cb->opcodeMinValue(op);
  float maxVal = cb->opcodeMaxValue(op);
  if(minVal > 0.f && maxVal > 0.f && maxVal >= 50.f*minVal)
    return true;
    // This criterion to switch between linear and exponential mode is ad hoc and heuristic. It 
    // doesn't seem to work though, because for cutoff, the minVal is actually 0. Maybe we need 
    // some expWithOffset characteristic for that. Maybe the offset should be maxVal / rangeFactor 
    // where rangeFactor is 20000/20 = 1000 in the case of frequencies, i.e. the maximum over the 
    // minimum meaningful value. Or maybe something based on sinh?
}

void SfzOpcodeWidgetSet::resized()
{
  int w = getWidth();
  int h = getHeight();
  slider->setBounds(0, 0, w, h);
  comboBox->setBounds(0, 0, w, h);
  textField->setBounds(0, 0, w, h);
}

void SfzOpcodeWidgetSet::setSfzNodeToEdit(const SfzNodeData& nodeData)
{
  // ToDo: intead of taking a rosic::Sampler::PlaybackSetting parameter, take the full SfzNodeData
  // object. It has the playbackSetting as member and contains some more information - in 
  // particular, it also contians the info, whether it's a regular playback setting or a modulation
  // setting, etc.

  using namespace rosic::Sampler;
  using OF = OpcodeFormat;
  using WM = WidgetMode;


  if(!nodeData.isOpcodeNode())
  {
    patchChangeInfo.reset();  // We clean up what has become meaningless
    setWidgetMode(WM::none);
    return;
  }

  int groupIndex = nodeData.getGroupIndex();
  int regionIndex = nodeData.getRegionIndex();

  rosic::Sampler::PlaybackSetting setting = nodeData.getPlaybackSetting();
  // todo: retrieve that only after we had made sure that the data at the node actually is a
  // PlaybackSetting by putting it in an if-block. It could also be a ModulationRouting or a 
  // <region> tag

  SfzCodeBook* cb = SfzCodeBook::getInstance();

  // Figure out what kind of opcode we are dealing with and set up the widgetMode accordingly:
  Opcode op = setting.getOpcode();
  int idx = setting.getIndex();
  std::string opStr = cb->opcodeToString(op, idx);
  OF fmt = cb->getOpcodeFormat(op);

  // Retrieve value, range and default value:
  float val    = setting.getValue();
  float minVal = jmin(cb->opcodeMinValue(op), val); // jmin/jmax because values in the code may go
  float maxVal = jmax(cb->opcodeMaxValue(op), val); // beyond the nominal range in SFZ spec
  float defVal = cb->opcodeDefaultValue(op, idx);

  // Update our info about what we actually edit:
  patchChangeInfo.type = PatchChangeType::opcodeValueChanged;
  patchChangeInfo.groupIndex = groupIndex;
  patchChangeInfo.regionIndex = regionIndex;
  patchChangeInfo.oldSetting = setting;
  patchChangeInfo.newValue = val; 

  // Display the appropriate widget and set it up:
  if(fmt == OF::Float || fmt == OF::Integer || fmt == OF::Natural)
  {
    double sliderInterval = 0.0;
    if(fmt != OF::Float )
      sliderInterval = 1.0;
    slider->setRange(minVal, maxVal, sliderInterval, defVal, false);
    slider->setValue(val, false);
    slider->setSliderName(opStr);
    setWidgetMode(WM::slider);
    if(wantsExponentialSlider(op))
      slider->setScaling(jura::Parameter::scalings::EXPONENTIAL);
    else
      slider->setScaling(jura::Parameter::scalings::LINEAR);
  }
  else if(fmt == OF::String)
  {
    // Maybe we need to distinguish between choice-parameters and free text. In case of the former,
    // we want acombobox, in case of the latter, a text-edit-fiel
  }
  else
  {
    setWidgetMode(WM::none);
  }

  // ToDo: 
  // -Maybe have a different quantization interval depending on the parameter
  // -Maybe also have linear or exponential scaling depending on parameter

  int dummy = 0;
}

void SfzOpcodeWidgetSet::handlePatchUpdate(const PatchChangeInfo& info)
{
  updateWidgetContent(info);
  patchChangeInfo.oldSetting.setValue(info.newValue);
}

void SfzOpcodeWidgetSet::rSliderValueChanged(RSlider* s)
{
  patchChangeInfo.newValue = s->getValue(); // Store the desired new value
  notifyMediator(0, &patchChangeInfo);      // Notify colleague objects (TreeView, CodeEditor, ...)
}

void SfzOpcodeWidgetSet::rComboBoxChanged(RComboBox* cb)
{
  RAPT::rsError("Not yet implemented");
}

void SfzOpcodeWidgetSet::textChanged(RTextEntryField* tf)
{
  RAPT::rsError("Not yet implemented");
}

void SfzOpcodeWidgetSet::setWidgetMode(WidgetMode newMode)
{
  if(newMode != mode)
  {
    mode = newMode;
    updateVisibilities();
  }

  // May be obsolete - verify if it still applies. This comment is older, from before factoring out
  // SfzOpcodeWidgetSet from SfzOpcodeEditor. 
  repaint();
  // Without calling repaint() here, sometimes the slider doesn't update correctly when selecting a
  // new parameter in the TreeView. For example, in teh patch NoiseWhistle.sfz, selecting first 
  // cutoff and then resonance, the resonance slider is correctly displayed only when calling 
  // repaint here
  // ToDo: 
  // -Maybe use something like repaintOnMessageThread() or repaintOnMessageThread(this) 
  //  instead.
}

void SfzOpcodeWidgetSet::createWidgets()
{
  addWidget(slider = new jura::RSlider(), true, false);
  slider->addListener(this);

  addWidget(comboBox = new jura::RComboBox(), true, false);
  comboBox->registerComboBoxObserver(this);

  addWidget(textField = new jura::RTextEntryField(), true, false);
  textField->registerTextEntryFieldObserver(this);

  // May add descriptions to these widgets, too - but maybe these descriptions should also change
  // dynamically? We'll see
}

void SfzOpcodeWidgetSet::updateWidgetContent(const PatchChangeInfo& info)
{
  if(mode == WidgetMode::slider) {
    if(slider->getValue() != info.newValue)
      slider->setValue(info.newValue, false); }
}

void SfzOpcodeWidgetSet::updateVisibilities()
{
  slider->setVisible(false);
  comboBox->setVisible(false);
  textField->setVisible(false);
  switch(mode)
  {
  case WidgetMode::slider:  slider->setVisible(true);    break;
  case WidgetMode::chooser: comboBox->setVisible(true);  break;
  case WidgetMode::text:    textField->setVisible(true); break;
  }
}

//=================================================================================================

SfzTreeView::SfzTreeView()
{
  rootNode.setDeleteChildNodesOnDestruction(true); // root node manages lifetime of child nodes
  rootNode.setNodeText("<global>");
  RTreeView::setRootNode(&rootNode);
  RTreeView::setDrawRootNode(false);
  createWidgets();
}

void SfzTreeView::buildTreeFromSfz(const rosic::Sampler::SfzInstrument& sfz)
{
  //return; 
  // For debug - when uncommenting this, the GUI becomes much more responsive. Try to optimize the
  // GUI such that it remains responsive when we actually run the code of this function

  // Perhaps we should first check, if the current tree alreday is in sync with the given sfz and 
  // if so, avoid re-building the tree from scratch. Maybe such situation could occur when the user
  // has just edited some comments in the code editor. Such edits don't change the structure of the
  // instrument. Maybe if the structure matches and just some numerical opcode values differ, we
  // can also avoid rebuilding the tree and just update the tree nodes accordingly.

  // Some shorthands:
  using namespace rosic::Sampler;
  using Node     = SfzTreeViewNode;
  using Level    = SfzInstrument::HierarchyLevel;
  using Global   = SfzInstrument::Global;
  using Group    = SfzInstrument::Group;
  using Region   = SfzInstrument::Region;

  using OC = Opcode;
  using PS = PlaybackSetting;
  using MR = ModulationRouting;

  using Settings = std::vector<PS>;
  using Routings = std::vector<MR>;
  // Maybe get rid of some of them - we seem to use them only one time so defining these aliases 
  // doesn't seem to pull its weight.


  SfzCodeBook* cb = SfzCodeBook::getInstance();

  auto addSettingNode = [&](Node* parent, const PS& setting, int groupIndex, int regionIndex)
  {
    Node* child = Node::createPlaybackSettingNode(
      cb->settingToString(setting), groupIndex, regionIndex, setting);
    parent->addChildNode(child);
  };

  auto addModRoutingNode = [&](Node* parent, const MR& routing, int groupIndex, int regionIndex)
  {
    Node* child = Node::createModulationRoutingNode(
      cb->modRoutingToString(routing), groupIndex, regionIndex, routing);
    parent->addChildNode(child);
  };

  // Helper function to add the leaf nodes:
  auto addOpcodeChildNodes = [&](const Level* lvl, Node* node, int groupIndex, int regionIndex)
  {
    int gi = groupIndex, ri = regionIndex;

    // Add nodes for certain special opcode settings. We display the sample opcode only if a sample
    // is defined at this level and show the key/vel settings only when they are not at their 
    // defaults:
    std::string s; int i;
    s = lvl->getSamplePath(); 
    if(s !=  "")
    {
      PS ps(OC::Sample);
      node->addChildNode(Node::createPlaybackSettingNode(("sample=" + s), gi, ri, ps)); 
    }
     
    i = lvl->getLoKey(); if(i !=   0) addSettingNode(node, PS(OC::LoKey, i), gi, ri);
    i = lvl->getHiKey(); if(i != 127) addSettingNode(node, PS(OC::HiKey, i), gi, ri);
    i = lvl->getLoVel(); if(i !=   0) addSettingNode(node, PS(OC::LoVel, i), gi, ri);
    i = lvl->getHiVel(); if(i != 127) addSettingNode(node, PS(OC::HiVel, i), gi, ri);
    // Values are shown like 45.0000 - we need a better float-to-string conversion

    // Add nodes for the general opcode settings:
    const Settings& settings = lvl->getSettings();
    for(int i = 0; i < settings.size(); i++)
      addSettingNode(node, settings[i], groupIndex, regionIndex);

    // Add nodes for modulation routings:
    const Routings& routings = lvl->getModRoutings();
    for(int i = 0; i < routings.size(); i++)
      addModRoutingNode(node, routings[i], groupIndex, regionIndex);

    // ToDo: 
    // -Later, we may want to allow more than one sample per level. Then we need to display the 
    //  sample index, too. Maybe we should be able to treat the sample-opcode just like all the
    //  others. It's unelegant to have these exceptions.
  };

  // Build the tree with 3 hierarchy levels:
  rootNode.deleteChildNodesRecursively();
  const Global* global = sfz.getGlobal();
  addOpcodeChildNodes(global, &rootNode, -1, -1);
  for(int gi = 0; gi < sfz.getNumGroups(); gi++)
  {
    Node* groupNode = Node::createGroupNode("<group>", gi);
    const Group* group = sfz.getGroup(gi);
    addOpcodeChildNodes(group, groupNode, gi, -1);
    for(int ri = 0; ri < sfz.getNumRegions(gi); ri++)
    {
      Node* regionNode = Node::createRegionNode("<region>", gi, ri);
      const Region* region = sfz.getRegion(gi, ri);
      addOpcodeChildNodes(region, regionNode, gi, ri);
      groupNode->addChildNode(regionNode);
    }
    rootNode.addChildNode(groupNode);
  }

  repaintOnMessageThread();

  // -The update of the GUI is quite slow/laggy. I don't think that it's a general problem with the
  //  TreeView because it's much more responsive in the module-selector in ToolChain. Figure out 
  //  what makes it so unresponsive in this context. Do we get bombarded with lots of redundant
  //  callbacks or something? Or is it because we allocate the nodes on the heap and they are more 
  //  scattered in memory than in the other case? Maybe we should pre-allocate a pool of nodes?
  //  Use Visual Studio's profiler...maybe try also CodeBlocks and XCode profiling
  //  ...btw: this was already the case before SfzTreeViewNode got its additional datafields
  //  commenting out the repaintOnMessageThread calls does not seem to help. Even the patch-loading
  //  and editor update takes unreasobaly long. Try to tmporarily get rid of the whole tree-view.
  //  OK - commenting out the line addWidget(sfzTree) in createWidgets, does indeed make patch 
  //  loading much more responsive. Try making the TreeView smaller
  //  -> done -> hmm...it helps a little bit but not that much
  //  returning early from RTreeView::paint or RTreeView::drawNode helps. It is apparently the 
  //  painting that bogs down the machine. Maybe because we alway draw all nodes, even if they
  //  are not within the visible area?
  //  ...ooookayy...after adding the optimization to draw only the visible nodes, it got a lot
  //  better. ...still a bit sluggish, though...although...actually, now it seems to feel quite OK,
  //  at least in a release build
  // -Scrollbars do not correctly appear/disappear
  // -The TreeView seems to update/repaint itself only on mouseOver after loading a new patch. It 
  //  should update immediately. -> fixed by calling repaintOnMessageThread() at the end
  // -The instruments with 1 regions show 2 region nodes - Test with FilterBlip.sfz. the sft does
  //  indeed have two regions - the first contains all the settings, the 2nd only the sample. Oddly
  //  enough, the FilterBlipe patch is nevertheless playable. It actually shouldn't be...i think.
  // -Maybe the function should also get the sfz-string passed along and identify and store the
  //  code-locations in the nodes - but finding them is a nontrivial task.
  // -Maybe trying to edit the patch via the tree is a too tall order. Maybe use the tree only for
  //  display purposes and for selecting opcodes to interact with via sliders but leave the patch 
  //  editing to the code editor. Maybe put the editor to the left and the tree to the right to 
  //  clarify the dirction of information flow (editor -> tree). The "Parse" button may somehow be 
  //  associated with this "flow"...clikcing it, makes the information "flow" from code editor to 
  //  tree - so maybe place it somewhere at the intersection. Maybe indicate invalidation/
  //  dirtification of the tree after code edits

  // ToDo:
  // -Maybe show the <global> section explicitly
  // -Show the <control> section
}
// needs test

void SfzTreeView::clearTree()
{
  RTreeView::setRootNode(nullptr);
  repaintOnMessageThread();
}
// needs test

void SfzTreeView::handlePatchUpdate(const PatchChangeInfo& info)
{
  SfzTreeViewNode* node = findNode(info);
  jassert(node != nullptr);
  if(node != nullptr)
  {
    // Update the content of the node:
    rosic::Sampler::PlaybackSetting newSetting = info.oldSetting;
    newSetting.setValue(info.newValue);
    rosic::Sampler::SfzCodeBook* cb = rosic::Sampler::SfzCodeBook::getInstance();
    std::string newText = cb->settingToString(newSetting);
    node->setNodeText(newText);

    node->data.data.playbackSetting = newSetting;
    // I'm not sure, if it's a good idea to have a plybackSetting member in the node. Isn't that 
    // redundant with the patchChangeInfo.oldSetting field in SfzOpcodeEditor? Figure out! But 
    // maybe the we shoul rather get rid of that "oldSetting" field. Not sure yet what's best.
    // Also, the node may also represent a mod-routing rather than a playbackSetting, so if we keep
    // it, we need to dispatch here.

    repaintOnMessageThread();
    // Calling repaint here makes the opcode slider a bit unresponsive but not calling it will 
    // leave the tree node dirty, i.e. it will continue to show an outdated value which will update
    // itself only when we hover with the mouse over the tree-view.
  }
}

SfzTreeViewNode* SfzTreeView::findNode(const PatchChangeInfo& info)
{
  int gi = info.groupIndex;
  int ri = info.regionIndex;
  rosic::Sampler::Opcode op = info.oldSetting.getOpcode();
  int idx = info.oldSetting.getIndex();

  // Find the node that corresponds to the hierarchy level that contains the opcode that needs to
  // be updated:
  RTreeViewNode* node = &rootNode;
  if(gi >= 0)
    node = getGroupNode(node, gi);
  if(ri >= 0)
    node = getRegionNode(node, ri);

  // Find the node for the opcode itself and return it as SfzTreeViewNode:
  node = getOpcodeNode(node, info);
  return dynamic_cast<SfzTreeViewNode*>(node);

  /*
  // Find the node for the opcode itself:
  node = getOpcodeNode(node, info);

  SfzTreeViewNode* sfzNode = dynamic_cast<SfzTreeViewNode*>(node);
  jassert(sfzNode != nullptr);
  return sfzNode;
  */

  // ToDo:
  // -Clean up jasserts - some are redundant with those in the caller
  // -Maybe streamline the code - get rid of some local assignments
}

jura::RTreeViewNode* SfzTreeView::getGroupNode(jura::RTreeViewNode* parent, int groupIndex)
{
  jassert(parent != nullptr);
  return parent->findDirectChildByText("<group>", groupIndex);
}

jura::RTreeViewNode* SfzTreeView::getRegionNode(jura::RTreeViewNode* parent, int regionIndex)
{
  jassert(parent != nullptr);
  return parent->findDirectChildByText("<region>", regionIndex);
}

jura::RTreeViewNode* SfzTreeView::getOpcodeNode(jura::RTreeViewNode* parent, 
  const PatchChangeInfo& info)
{
  jassert(parent != nullptr);
  rosic::Sampler::SfzCodeBook* cb = rosic::Sampler::SfzCodeBook::getInstance();
  std::string str = cb->settingToString(info.oldSetting);
  return parent->findDirectChildByText(str, 0);
  // We pass 0 to the index parameter of findDirectChildByText because the string already uniquely
  // identifies the node. The index (and more) is already baked into the search-string.
}

void SfzTreeView::mouseMove(const MouseEvent& e)
{
  if(e.mods.isCtrlDown())
  {
    int y = e.getPosition().y;
    //RTreeViewNode* node = getNodeAtY(y);
    SfzTreeViewNode* node = dynamic_cast<SfzTreeViewNode*> (getNodeAtY(y));
    if(node != nullptr)
      showOverlayWidget(node, y);
  }
  else
  {
    hideOverlayWidgets();
    jura::RTreeView::mouseMove(e);
  }
}

void SfzTreeView::mouseExit(const MouseEvent& e)
{
  //MouseEvent e2 = e.getEventRelativeTo(this); // seems to be the same as e
  int margin = 4;
  if(e.x < margin || e.x > getWidth()-margin)
    hideOverlayWidgets();
  if(e.y < margin || e.y > getHeight()-margin)
    hideOverlayWidgets();

  jura::RTreeView::mouseExit(e);

  // Just calling hideOverlayWidgets() no matter what prevents the widgets from appearing in the 
  // first place. I guess mouseExit events are spawned also when the overlay widget appears because
  // then it is under the mouse. The intention is to make the overlay widgets disappear whenever 
  // the mouse exits the tree-view - but mousing over the overlay-widget should not count as 
  // exiting the tree-view. So we check the mouse-coordinates before calling hideOverlayWidgets
}

void SfzTreeView::setMediator(Mediator* newMediator)
{
  SamplerInterfaceComponent::setMediator(newMediator);
  overlayWidgets->setMediator(newMediator);
}

void SfzTreeView::setColourScheme(const WidgetColourScheme& newColourScheme)
{
  RTreeView::setColourScheme(newColourScheme);
  jassert(overlayWidgets != nullptr);
  overlayWidgets->setColourScheme(newColourScheme);
  // What if this function is called before overlayWidgets is assigned? Should we test for
  // nullptr first
}

void SfzTreeView::createWidgets()
{ 
  addChildComponent(overlayWidgets = new SfzOpcodeWidgetSet());
}

void SfzTreeView::hideOverlayWidgets()
{
  overlayWidgets->setVisible(false);
}

void SfzTreeView::showOverlayWidget(SfzTreeViewNode* node, int y)
{
  if(node->isOpcodeNode())  // ToDo: Check node against nullptr before deref or assert non-null
  {
    y -= overlayWidgets->getHeight() / 2;
    overlayWidgets->setBounds(16, y, getWidth()-32, 16);
    overlayWidgets->setVisible(true);
    overlayWidgets->setSfzNodeToEdit(node->getNodeData()); 
    int dummy = 0;
    //RAPT::rsError("Not yet finished and tested");
  }
  else
  {
    overlayWidgets->setSfzNodeToEdit(SfzNodeData::createEmptyNode());
    //overlayWidgets->setVisible(false);
  }
}

//=================================================================================================

SfzOpcodeEditor::SfzOpcodeEditor()
{
  createWidgets();
}

void SfzOpcodeEditor::handlePatchUpdate(const PatchChangeInfo& info)
{
  //patchChangeInfo.oldSetting.setValue(info.newValue);  // old

  int dummy = 0;
  // Maybe we can leave this function empty? I don't really see what we would need to do here. At 
  // the moment only the widget needs to be updated which takes care of itself because only the 
  // widget triggers such patch changes. But perhaps later, we could visualize the effect of the 
  // parameter somehow? Maybe a plot of a freq-response in case of a cutoff or resonance parameter?
  // We'll see....
  // Or: Maybe if such a patch-update event later can be triggered by something else than the 
  // opcode-widget, we may have to take care to update the widget here manually. But at the moment, 
  // there seems to be nothing to do here. Of course, we can also change values in tne code but 
  // that actually triggers a full re-parse at the moment (but that may change, too)
}

void SfzOpcodeEditor::setMediator(Mediator* newMediator)
{
  SamplerInterfaceComponent::setMediator(newMediator);
  opcodeWidgets->setMediator(newMediator);
}

void SfzOpcodeEditor::resized()
{
  int m = 4; // margin;
  int x = m;
  int y = m;
  int w = getWidth() - x - m;
  int h = 16;
  opcodeField->setBounds(x, y, w, h); // todo: Maybe center it in the editor
  y += h + m;
  helpField->setBounds(x, y, w, h); 
  y += h + 2*m;
  opcodeWidgets->setBounds(x, y, w, h);

  // Maybe the opcodeField should just say e.g. volume= and the slider should be right next to it
  // showing and adjusting the value. The help text can appear below
}

void SfzOpcodeEditor::createWidgets()
{
  // The widgets that are always present:
  addWidget(opcodeField = new jura::RTextField());
  opcodeField->setText("Opcode"); // preliminary...maybe get rid
  opcodeField->setDescription("Name of currently active opcode");

  addWidget(helpField = new jura::RTextField());
  helpField->setText("ToDo: Opcode description goes here");   // also preliminary
  helpField->setDescription("Short description of the opcode");

  addChildColourSchemeComponent(opcodeWidgets = new SfzOpcodeWidgetSet());
  // Maybe add descriptions to these widgets, too - but maybe these descriptions should also 
  // change dynamically? We'll see...
}

void SfzOpcodeEditor::updateVisibilities()
{
  opcodeWidgets->updateVisibilities();  // Needed? Verify!
}

//=================================================================================================

void SfzCodeEditor::handlePatchUpdate(const PatchChangeInfo& info)
{
  // ToDo - Maybe:
  // These patch-update callbacks can potentially come in quite frequently, for example when they
  // are generated by tweaking a slider. It may make sense to not immediately do this rather costly
  // text processing here but rather just set some sort of dirty flag and schedule the update for 
  // later. Or maybe just enqeue the info into a fifo buffer to be handled later where we may also
  // consolidate several updates of the same parameter into a single update. Maybe SfzCodeEditor 
  // could derive from juce::Timer, set the callback frequency to something like 0.5 seconds and do 
  // the text updates inside the timer callback. A similar strategy could also be applied to handle
  // updating the set_ccN opcodes in response to receiving midi cc messages.

  // Find the location in the code that is affected by the change:
  int startPos, endPos;
  findCodeSegment(info, &startPos, &endPos);
  if(startPos == -1 || endPos == -1) {
    //RAPT::rsError("No appropriate code segment found"); 
    return; }
    // I'm not yet sure, if we should trigger an error in such a situation. Maybe we should expect
    // that to happen in normal operation? We'll see....

  // Apply the required change to the code document:
  juce::CodeDocument& doc = CodeEditorComponent::getDocument();


  juce::String newValueString = juce::String(info.newValue);
  // ToDo: let info have a function getNewValueString() that dispatches between string-values and 
  // numeric values depending on the opcode...or maybe have that dispatcher functionality somewhere
  // else - but currently, we just assume a numeric value which is sometimes not the case...

  // doc.replaceSection(startPos, endPos, newValueString);  // that was wrong
  doc.replaceSection(startPos, endPos+1, newValueString);   // fixed 2023/05/21
  // But: maybe the off-by-one error is not here but inside findCodeSegment (inside one of the 
  // functions it calls)? Check that! Or maybe it's a different interpretation of what "endPos"
  // means? findCodeSegment includes the endPos and replaceSection excludes it ot vice versa?
  // Figure that out! I think my API uses the closed interval [start, end] but juce uses the 
  // half-open interval [start, end)?

  // Note:
  // The change we make to the code should not trigger a re-parsing, at least not automatically 
  // (but that doesn't really happen anyway). Otherwise, we would reparse on slider movement which 
  // is obviously totally silly and impractical. If the "Parse" button appears, that might be 
  // tolerable, although it would be better if it doesn't.

  // Done:
  // -Find the location in the code that is affected by the change, write a function 
  //  findCodeToModify. It should perhaps return a pair of indices or an index and a length for 
  //  where the relevant opcode (or region, etc.) substring starts and ends.
  // -Change the code there. Replace the substring with the appropriate new string. Don't reparse.
  //  Otherwise, we would reparse on slider movement which is obviously totally silly and 
  //  impractical

  // ToDo:
  // -We need to update not only the sfz-code but also the sfz-datastructure. Currently, when 
  //  changing a value via a slider, the value in the tree is not updated. When we then select 
  //  another opcode in the tree and the the previous opcode again, the slider will start out at
  //  the previous value. But that should be done inside the callbacks of the other colleagues, not
  //  here. (check, if that's still the case)
  // -Actually, it would be better to find the code-segment (or at least its start), when the user
  //  selects a new node in the tree - not on every slider-movement. The starting position does not
  //  change. The length may, depending on the text-formatting of floating point numbers and also
  //  when we are dealing with a choice opcode. But even the length of the segment is a thing, we 
  //  may keep track of without repeatedly figuring it out again and again. 
  // -Maybe include a test that the updated document indeed reflects the settings of the engine.
  //  Retrieve the sfz-data, re-parse the code, check if there is a match. This should be done only
  //  in debug situations.
}

void SfzCodeEditor::handleMidiUpdate(const rosic::Sampler::rsMusicalEvent<float>& ev)
{
  int startPos, endPos;
  findCodeSegment(ev, &startPos, &endPos);
  if(startPos == -1 || endPos == -1) {
    //RAPT::rsError("No appropriate code segment found"); 
    return; }

  // Apply the required change to the code document:
  juce::CodeDocument& doc = CodeEditorComponent::getDocument();
  juce::String newValueString = juce::String(ev.getValue2());   // maybe use custom function
  doc.replaceSection(startPos, endPos+1, newValueString);

  // Possible optimization idea:
  // Maybe we should not immediately update the code but rather schedule an update. Doing it 
  // immediately may be a bit costly because this function gets called from the timerCallback in
  // the SamplerEditor class and this timer has a rather high frequency (like 50 Hz). Maybe the 
  // code editor should use a less frequent update. Maybe every second or half or so. I guess text 
  // replacements are a rather costly operation that we don't want to do at 50 Hz rate. We could
  // derive the SfzCodeEditor also from juce::Timer and set a "dirty" flag here and in our own 
  // timerCallback, we do the update if the flag is dirty and then clear the flag.
}

void SfzCodeEditor::findCodeSegment(const PatchChangeInfo& info, int* startPos, int* endPos)
{
  // Maybe streamline these calls into a chain like getDocument().getallContent().toStdString().
  // But maybe we actually need to hold the reference to the doc for calling replaceSection() on 
  // it. We'll see..
  juce::CodeDocument& doc = CodeEditorComponent::getDocument();
  juce::String jStr = doc.getAllContent();
  std::string  code = jStr.toStdString();

  rosic::Sampler::SfzCodeBook::getInstance()->findOpcodeValueString(code, info.groupIndex, 
    info.regionIndex, info.oldSetting.getOpcode(), info.oldSetting.getIndex(), startPos, endPos);

  // ToDo:
  // -Try to avoid the conversion from juce::String to std::string by letting findOpcodeValueString
  //  operate on const char*
  // -Make sure, this function is called only when the user selects a new opcode from the tree, not
  //  on every slider movement as we do now.
  // -Make it so that we immediately hear the effect of the new setting, not only on a new noteOn.
  //  To achieve this, we must decouple the setup function from the noteOn, i.e. not only set up
  //  the modules on note on and then leave them be but allow to call setup during the module is
  //  running.
}

void SfzCodeEditor::findCodeSegment(const rosic::Sampler::rsMusicalEvent<float>& ev, 
  int* startPos, int* endPos)
{
  // This is the same silly conversion as in the other findCodeSegment function. Try to get rid of
  // that.
  juce::CodeDocument& doc = CodeEditorComponent::getDocument();
  juce::String jStr = doc.getAllContent();
  std::string  code = jStr.toStdString();

  auto  cb = rosic::Sampler::SfzCodeBook::getInstance();
  using Ev = rosic::Sampler::rsMusicalEvent<float>;
  using OC = rosic::Sampler::Opcode;
  if(ev.getType() == Ev::Type::controlChange)
    cb->findMidiControllerValueString(code, (int) ev.getValue1(), startPos, endPos);
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
  setSize(900, 540);  // width should be divisible by 3 - ToDo: document why
  // Maybe choose an initial size such that lines with a 100 character limit fit exactly into the 
  // code editor window. That's my preferred line-width in code. Hmm - that's a bit wide. Maybe
  // use a limit of 80 in the sfz files. Also, maybe remove the line-numbering from the code editor
  // to get more space for the actual code

  initOldControllersCache(); 
  // It's perhaps a good idea to do that before starting the timer as the timer callbacks may 
  // trigger actions look into this cached array of old controller values. Maybe we should also
  // do this in setVisible()? Maybe we should also start the timer when the GUI becomes visible and
  // stop it when it becomes invisible?

  //startTimer(20);  // in ms
  startTimerHz(50);   // in Hz, i.e. fps (frames per second) for the metering widgets
  // ToDo: maybe stop the timer when in edit mode - only in play mode, the metering widgets are 
  // shown


  sfzDoc.addListener(this);
  parseButton->addRButtonListener(this);

  samplerModule->sfzPlayer.addFileManagerListener(this); 
  // We want to receive activeFileChanged callbacks when the currently active sfz file changes.


  //setCodeIsClean();  //

  connectGuiElementsToMediator();
}

SamplerEditor::~SamplerEditor()
{
  ScopedLock scopedLock(*lock);
  samplerModule->sfzPlayer.removeFileManagerListener(this);
}

void SamplerEditor::timerCallback()
{
  jassert(samplerModule != nullptr);

  // Update the metering display for the active layers:
  int num = samplerModule->getNumActiveLayers();
  layersMeter->setCurrentValue((float)num);
  // see also: TrackMeterModuleEditor::timerCallback

  // Update the sliders for the MIDI controllers, if necessary:
  if(samplerModule->isMidiControlStateDirty()) {
    for(int i = 0; i < numCtrlSliders; i++) {
      float val = samplerModule->getMidiControllerCurrentValue(i);
      if(val != oldCtrlValues[i]) {
        using Event = rosic::Sampler::rsMusicalEvent<float>;
        Event ev(Event::Type::controlChange, (float)i, val);
        sfzEditor.handleMidiUpdate(ev);
        if(ctrlSliders[i]->isVisible())
          ctrlSliders[i]->setValue(val);           // seems costly - see below
        oldCtrlValues[i] = val;  }}                // prepare for the next cc event
    samplerModule->setMidiControlStateClean(); }

  // ToDo:
  // According to some test with a release build and a system CPU monitor, the 
  // ctrlSliders[i]->setValue updates seem to be quite costly! There's a big difference in CPU load
  // when this update is commented out. Why is that? It's just a basic slider! Check the drawing 
  // code! Maybe it can be optimized. Or maybe there's something else going on in setValue?
}

void SamplerEditor::resized()
{
  ScopedLock scopedLock(*lock);  // do we actually nee the lock here?

  // Let the baseclass take care of positioning the xml preset widget set, then reduce its width to 
  // make space for the buttons that switch between the GUI pages:
  AudioModuleEditor::resized();
  int x  = stateWidgetSet->getX();
  int y  = stateWidgetSet->getY();
  int h  = stateWidgetSet->getHeight();
  int m  = 4;                              // margin
  int bw = 48;                             // button width
  int w  = getWidth() - x - 2*(bw+m) - m;  // new width for preset section
  w = jmin(w, 360);
  stateWidgetSet->setBounds(x, y, w, h);
  // Actually, that's still a bit too wide. But maybe later we'll get more pages such that we need
  // more buttons. Then we can reduce the width further to make more space. Maybe then the 
  // 2*(bw + m) should be replaced by numButtons*(bw + m). We could actually reduce the width 
  // regardless of whether or not we need more page buttons and maybe just leave some empty space. 
  // We'll see....

  // Place the sfz load/save widgets to the extreme left of the editor, directly on top of the code 
  // editor. Next to these widgets, we place the GUI page buttons. The y-position should leave a 
  // comfortable vertical space between xml loader and sfz loader:
  y = getPresetSectionBottom() + 2*m;
  x = 0;
  w = 360;
  sfzFileLoader->setBounds(x, y, w, h);
  x = sfzFileLoader->getRight() + 2*m; 
  playButton->setBounds(x, y, bw, h);
  x += bw-2;
  editButton->setBounds(x, y, bw, h);

  // Lay out widgets for "Play" page:
  // ToDo: 
  // -Split the page into 3 regions: left/center/right
  // -Into the left region goes some information about the patch (author, comments, number of 
  //  samples, regions, groups, etc.)
  // -Into the center region go realtime performance controllers like sliders, vector-pads etc. 
  //  that can be assigned to midi controllers.
  // -Into right section go metering devices like RAM and CPU usage, maybe a scope, output levels,
  //  stereo width etc.
  // -Hmm...or maybe it's better to put the performance controls left, patch info centered? 
  x = 2 * getWidth() / 3 + m;
  y = sfzFileLoader->getBottom() + m;
  w = getWidth() / 3 - 2*m;
  h = 20;
  layersMeter->setBounds(x, y, w, h); 

  // Lay out widgets for "Edit" page:
  // ToDo:
  // -Split the page into 3 sections: top-left, top-right and bottom. the top-sections should take
  //  roughly 3/4 of the vertical space and the bottom section the remaining 1/4 (or maybe use even
  //  4/5, 1/5). And the top-left should take roughly 2/3 to 3/4 of the horizontal space
  // -The top-left section is the sfz code editor
  // -The top-right section is the sfz structure tree-view
  // -Maybe the code editor should use the full height. We may not need that much space for the
  //  opcode edit widgets...although, for a sample preview, it's rather nice to have the full width
  //  available.
  // -The bottom section dynamically shows appropriate edit or display widgets for the node that is
  //  selected in the tree-view, like a slider for continuously adjustable opcodes, a sample-view 
  //  for sample opcodes (maybe it should allow to load samples using a file-browser)
  // -It would be nice, if we could sync the code-editor and the tree-view, i.e. highlight the code
  //  corresponding to the selected node and maybe even modify it via the widgets.
  x = 0;
  y = sfzFileLoader->getBottom();
  w = 2 * getWidth() / 3;
  //h = 3 * (getHeight() - y) / 4;    // 3/4 of remeining height
  h = getHeight() - y;
  sfzEditor.setBounds(x, y, w, h);
  x = sfzEditor.getRight();
  h = 3 * (getHeight() - y) / 4;      // 3/4 of remeining height
  w = getWidth() - x;
  sfzTree->setBounds(x, y, w, h);
  x = x - bw;
  y = sfzFileLoader->getY();
  w = bw; 
  h = 16;
  parseButton->setBounds(x, y, w, h);
  y = sfzFileLoader->getY();
  w = structureField->getWidthToFitText();
  x = (sfzTree->getX() + sfzTree->getRight() - w) / 2;  // centered above tree-view
  structureField->setBounds(x, y, w, h);
  x = sfzTree->getX();
  y = sfzTree->getBottom() - 2;
  w = getWidth()  - x;
  h = getHeight() - y;
  opcodeEditor->setBounds(x, y, w, h);


  // ToDo:
  // -The free space in the top-right corner can perhaps be used to display some status into like 
  //  number of playing notes, last received event or maybe a level meter. The most important of 
  //  the metering widgets which should be shown permanently. This can be done for other modules
  //  with big editors, too. Maybe make a baseclass AudioModuleEditorWithMetering. Maybe the stereo
  //  meters shoudl both start int the middle and expand to left and right.
  // 

  /*
  // Under construction: display the status of the code editor content (Parsed/Malformed/Edited)
  x = sfzFileLoader->getRight() + 24;
  sfzStatusField->setLabelWidth(64);
  sfzStatusField->setBounds(x, y, 150, 16);
  // or maybe move to the right, directly next to (left of) the parseButton.
  // Maybe don't have an extra field for that - reuse the tree. fill it only, if the instrument is
  // successfully parsed, otherwise show an empty tree and a message that tells the user that the
  // sfz could not be parsed...ideally with some more specific error message.
  */
}

void SamplerEditor::treeNodeClicked(RTreeView* treeView, RTreeViewNode* node, 
  const MouseEvent& mouseEvent,int clickPosition)
{
  // ToDo:
  // -Put the node (exclusively) in selected state
  // -Tell the OpcodeEditor about what kind of node is selected to make it update its widgets

  SfzTreeViewNode* sfzNode = dynamic_cast<SfzTreeViewNode*>(node);
  jassert(sfzNode);  // You should only add SfzTreeViewNodes into the tree!


  int gi = sfzNode->data.groupIndex;
  int ri = sfzNode->data.regionIndex;
  //using TP = SfzTreeViewNode::Data::Type;
  using TP = SfzNodeData::Type;
  switch(sfzNode->data.type)
  {
  case TP::playbackSetting:
  {
    //rosic::Sampler::PlaybackSetting ps = sfzNode->data.data.playbackSetting;
    // That's kinda ugly, especially the data.data part. Try to do better!
    //opcodeEditor->setSettingToEdit(gi, ri, ps);  // old

    opcodeEditor->setSfzNodeToEdit(sfzNode->data);  // new


    int dummy = 0;
  } break;
  case TP::modulationRouting:
  {

  } break;
  case TP::group:
  {

  } break;
  case TP::region:
  {

  } break;
  default:
  {

  }


  }

  int dummy = 0;
}

void SamplerEditor::treeNodeChanged(RTreeView* treeView, RTreeViewNode* node)
{

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
    parseCodeEditorContent();
  else if(b == playButton || b == editButton)
    updateVisibilities();
  else
    AudioModuleEditor::rButtonClicked(b);
}

void SamplerEditor::rSliderValueChanged(RSlider* sld)
{
  // Update MIDI controller value in the engine for the control sliders. I'm not so sure, if it's
  // a good idea to directly send a midi cc event to the engine on the GUI thread. Normally, MIDI 
  // is received in the audio thread. On the other hand, I can't really see why that should be 
  // problematic. So let's just do it:
  for(int i = 0; i < numCtrlSliders; i++) {
    if(sld == ctrlSliders[i]) {
      samplerModule->setMidiController(i, (float) sld->getValue()); }}
}

void SamplerEditor::activeFileChanged(FileManager* fileMan)
{
  // BUG:
  // This callback is invoked twice when loading a file via the load button. When loading a file
  // via forward/backward buttons, it's only invoked once, as it should be. Figure out why it gets
  // invoked twice in the former case and fix that!

  jassert(fileMan == &samplerModule->sfzPlayer);
  // There are actually two FileManager objects that we could listen to: the SamplerModule object
  // itself which manages the preset .xml file of the whole sampler and the SfzPlayer object 
  // embedded in the former which manages the currently loaded .sfz instrument definition. We 
  // currently listen only to the latter. Changes in the former are of no interest because a preset
  // change due to the user loading another xml will usually imply also a change of the sfz file, 
  // i.e. loading a new xml will trigger loading a new sfz. Which sfz should be loaded is defined 
  // in the xml preset. The xml-presets are on a higher level than the sfz instruments. Think of 
  // the loaded sfz instrument more like a waveform sample that is loaded into an oscillator. We 
  // treat it in a similar way.

  juce::String currentSfz = samplerModule->sfzPlayer.getCurrentSfz();
  sfzDoc.replaceAllContent(currentSfz);
  setCodeIsClean();
  // Needed because sfzDoc.replaceAllContent will trigger calls to setCodeIsDirty() but when we 
  // receive a callback from the samplerModule->engine here, then the code shown in the editor
  // actually is actually clean (i.e. parsed and saved).
  // ...verify this...

  // Initialize the old controller values from the engine:
  initOldControllersCache();

  
  //oldCtrlValues[numCtrlSliders];


  //samplerModule->markStateAsDirty();  // Dirtify the xml load/save widget
  // but should we do this in the module? i think so? I think, we should not reach through into the
  // embedded SfzPlayer but instead call samplerModule
  // ...nnnah - the samplerModule should isefl listen to the sfzPlayer and dirtify itself in the
  // callback

  updateTreeView();
  updateVisibilities();
}

void SamplerEditor::handlePatchUpdate(const PatchChangeInfo& info)
{
  rosic::Sampler::PlaybackSetting os = info.oldSetting;
  int gi = info.groupIndex;
  int ri = info.regionIndex;

  if(gi == -1)         // It's a global setting
    samplerModule->setInstrumentSetting(os.getOpcode(), info.newValue, os.getIndex());
  else if(ri == -1)    // It's a group setting
    samplerModule->setGroupSetting(gi, os.getOpcode(), info.newValue, os.getIndex());
  else                 // It's a region setting
    samplerModule->setRegionSetting(gi, ri, os.getOpcode(), info.newValue, os.getIndex());

  int dummy = 0;

  // ToDo:
  // -Inspect the return values of samplerModule->set... They return a success or error code. Maybe
  //  we should assert here that it returns the code for success. We don't expect it to fail here. 
  //  That would be a bug.
  // -Instead of checking gi, ri against -1 here, let info have member functions isGroupSetting, 
  //  isGlobalSetting, isRegionSetting. The intention is that the class PatchChangeInfo should be
  //  responsible for encoding and retrieving the information, which kind of setting is being 
  //  modified. We are basically looking into an implementation detail of PatchChangeInfo here 
  //  which is no good style (although perhaps more efficient - but we are not in a performance
  //  critical code section here).
  //  -isRegionSetting returns true, iff gi != 0 && ri != 0
  //  -isGlobalSetting returns true, iff gi == -1
  //  -isGroupSetting return true, iff gi >= 0 && ri == -1
}

void SamplerEditor::createWidgets()
{
  /*
  addWidget(instrumentLabel = new RTextField);
  instrumentLabel->setText("Instrument:");
  instrumentLabel->setDescription("Currently loaded sfz instrument");
  instrumentLabel->setDescriptionField(infoField);
  */


  addWidget(playButton = new jura::RRadioButton("Play"));
  playButton->setDescription("Switch GUI to play/perform mode");
  playButton->addToRadioButtonGroup(&guiPageButtonGroup);
  playButton->setToggleState(true, false);
  playButton->addRButtonListener(this);

  addWidget(editButton = new jura::RRadioButton("Edit"));
  editButton->setDescription("Switch GUI to SFZ editing mode");
  editButton->addToRadioButtonGroup(&guiPageButtonGroup);
  editButton->addRButtonListener(this);


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


  // SFZ file loading widgets:
  sfzFileLoader = new jura::FileSelectionBox("", &samplerModule->sfzPlayer);
  addWidgetSet(sfzFileLoader);
  //sfzFileLoader->setDescription("Code of Current SFZ instrument");  // ??? Do better!
  //sfzFileLoader->addFileManagerListener(this); 


  // The SFZ TreeView and adjacent widgets:
  juce::String desc = "Structure of current SFZ patch";
  addWidget(structureField = new RTextField("Patch Structure"));
  structureField->setDescription(desc);
  addWidget(sfzTree = new SfzTreeView);
  sfzTree->setDescription(desc);
  sfzTree->registerTreeViewObserver(this);



  //addWidget(sfzStatusField = new jura::RLabeledTextEntryField);
  //sfzStatusField->setLabelText("SFZ Status:");
  //sfzStatusField->setEntryEditable(false);

  // The SFZ CodeEditor:
  addAndMakeVisible(sfzEditor);
  sfzEditor.setLineNumbersShown(false);
  // ToDo: set up the description of the editor...but it's not a subclass of RWidget...hmmm... 
  // maybe we can make a subclass of juce::CodeEditorComponent and jura::DescribedItem

  addWidget(parseButton = new jura::RClickButton("Parse"));  // Maybe use a right-arrow ("Play")
  parseButton->setDescription("Parse the current content of the code editor as sfz");
  //parseButton->addRButtonListener(this); // Needed? It initially worked without that. ...but why?
  // ToDo: Ctrl-P should also trigger parsing. Write that into the description when it's 
  // implemented.
 
  addChildEditor(opcodeEditor = new jura::SfzOpcodeEditor);
  opcodeEditor->setDescription("Editor for the currently selected opcode in the tree");

  // The sliders for the MIDI controllers:
  for(int i = 0; i < 128; i++){
    jura::RSlider* sld = new jura::RSlider("MIDI CC " + juce::String(i));
    sld->setRange(0.0, 127.0, 1.0, 0.0);
    sld->setStringConversionFunction(&valueToString0);
    sld->addListener(this);
    addWidget(sld);
    ctrlSliders[i] = sld; }


  // The order in which we call add...() for the widgets determines, which widgets are in the 
  // foreground when they overlap.

  // ToDo: 
  // -set up better descriptions for the xml and sfz loader widget sets. Maybe:
  //  Load .xml preset from file, Load .sfz instrument from file, etc.
}

void SamplerEditor::connectGuiElementsToMediator()
{
  // ToDo: 
  // -We need to override setMediator in the sfzTree opcodeEditor to allow it to connect not only 
  //  itself but also its embedded SfzOpcodeWidgetSet

  this->setMediator(&guiMediator);
  sfzTree->setMediator(&guiMediator);       // needs override (now)
  sfzEditor.setMediator(&guiMediator);
  opcodeEditor->setMediator(&guiMediator);  // needs override (later)
  // It's important that opcodeEditor is registered last because it holds the patchChangeInfo which
  // must be updated last (or at least, after the sfzTree)

  updateVisibilities();
}

void SamplerEditor::setCodeIsParsed(bool isParsed)
{
  codeIsParsed = isParsed;
  //parseButton->setEnabled(!codeIsParsed); // doesn't seem to have any effect
  parseButton->setVisible(!codeIsParsed);   // works but is a bit drastic, graying out would be better
  // ToDo: 
  // -Figure out why graying out the button doesn't work and fix it
}

void SamplerEditor::setCodeIsSaved(bool isSaved)
{
  // ToDo: make an atserisk appear in the sfz file box, if is not saved, disappea, if it is saved
}

void SamplerEditor::parseCodeEditorContent()
{
  bool ok = samplerModule->setupFromSfzString(sfzEditor.getDocument().getAllContent(), false);
  if(!ok)
  {
    showWarningBox("Error", "Code editor content could not be parsed as SFZ.");
    // ToDo: maybe add some persistent warning sign to the GUI. Like a warning lamp or something. 
    // Maybe a red frame around the code editor. Or maybe an info/status field showing
    // Parsed/Unparsed/Malformed, Parsed/Edited/Malformed
    // ...maybe the frame around the code editor should be green if parsed...but maybe not. It may
    // destroy the color-scheme
  }

  setCodeIsParsed(ok);  // ...maybe we should do this only in case of success?
  updateTreeView();
}

/*
void SamplerEditor::saveCodeEditorContent()
{

}
// hmm - maybe this is not needed?
*/

void SamplerEditor::updateTreeView()
{
  if(codeIsParsed)
    sfzTree->buildTreeFromSfz(samplerModule->sfzPlayer.getInstrumentData());
  else
    sfzTree->clearTree();
  // I'm not yet sure, what the desired behavior should be in case of having unparsed code in the
  // code editor. It means that the current editor content failed to parse. The player will 
  // actually have reverted to the last parsable version though, so maybe we should display that 
  // in the TreeView. To do so, we would just have to get rid of the "if" and just always take the
  // first branch. We'll see...
}

void SamplerEditor::updateVisibilities()
{
  if(playButton->getToggleState())
  {
    makePlayWidgetsVisible(true);
    makeEditWidgetsVisible(false);
  }
  else
  {
    makePlayWidgetsVisible(false);
    makeEditWidgetsVisible(true);
  }

  // ToDo: 
  // -Maybe stop the timer when in edit mode - only in play mode, the metering widgets are 
  //  shown. Maybe the function should be renamed to switchGuiPage...but maybe we need it for the
  //  MIDI status and level meters anyway on all screens.
  // -Maybe rename the "Edit" button to "Code" and have have another "Edit" button that shows more
  //  user oriented editing widgets. We want to make a convenient sfz authoring tool.
}

void SamplerEditor::makePlayWidgetsVisible(bool visible)
{
  layersMeter->setVisible(visible);


  // The sliders for the MIDI controllers:
  if(visible)
  {
    int x = 16;
    int y = sfzFileLoader->getBottom() + 16;
    int w = 240;  // slider width
    int h = 20;   // slider height
    for(int i = 0; i < 128; i++)
    {
      juce::String ctrlLabel = samplerModule->getMidiControllerLabel(i);
      int value = samplerModule->getMidiControllerCurrentValue(i);
      if(ctrlLabel != "")
      {
        ctrlSliders[i]->setBounds(x, y, w, h);
        y += h-2;
        ctrlSliders[i]->setSliderName(ctrlLabel + " - CC #" + juce::String(i));
        ctrlSliders[i]->setValue(value);
        ctrlSliders[i]->setVisible(true);
      }
      else
        ctrlSliders[i]->setVisible(false);
    }
  }
  else
  {
    for(int i = 0; i < 128; i++)
      ctrlSliders[i]->setVisible(false);
  }


  
  // ToDo:
  // -Maybe make only those MIDI control sliders visible that have a name/label assigned. This will 
  //  require to call setBounds on these sliders, too
  // -Maybe setting the value should not depend on what's written in the sfz but what the last
  //  received MIDI CC value was?
}

void SamplerEditor::makeEditWidgetsVisible(bool visible)
{
  sfzTree->setVisible(visible);
  parseButton->setVisible(visible);
  sfzEditor.setVisible(visible);
  structureField->setVisible(visible);
  opcodeEditor->setVisible(visible);
}

void SamplerEditor::initOldControllersCache()
{
  for(int i = 0; i < numCtrlSliders; i++)
    oldCtrlValues[i] = samplerModule->getMidiControllerCurrentValue(i);
}

//=================================================================================================

/*

ToDo:

-When saving, it saves the most recently parsed document rather than the current content of the 
 editor - that is wrong!
-When we click on e.g. volume to make it appear in the opcode editor, the manipulate e.g. 
 pitch_keycenter via the overlay widgets, the volume slider updates itself to show the new
 keycenter value - which makes no sense. Maybe before updating the widget content, we should check,
 if the patchChange applies to the shown widget
-Keep the selcted TreeNode highlighted as long as the slider for it is visible. Also highlight 
 the relevant section of the code in the editor.
-The slider needs exponential characteristic for certain parameters. Maybe to start, just use a 
 heuristic: if min and max are both strictly greater than 0 and max >= k*min for some constant k
 like k=50, then use exp-mode. Maybe make a function wantsExponentialSlider(Opcode op). There's 
 some preliminary code in SfzOpcodeEditor::setSettingToEdit - refine that. See also the commented
 SfzOpcodeEditor::wantsExponentialSlider function. We may need a more flexible way for dealing 
 with the mapping, so we may need to use a jura::Parameter object to whcih we can assign a
 ParameterMapper.
-The currently active parameter should stay selected in the tree as long as we can tweak it with 
 the slider -> better visual feedback for the user
-It would be nice if the change would be audibel immeidately wehn movin the slider -> we need to 
 decouple the "setup" function form "noteOn" in the engine
-We canot yet properly edit mod-routing parameters
-Maybe SamplerModule should override parameterChanged. There, it should figure out, if the changed
 parameter was an opcode-parameter (maybe via dynamic_cast) and if so, update the corresponding
 sfz-opcode in the engine's currently loaded instrument. We also somehow make the editor aware of 
 it to let it update its contents.
-Status: Parsed/Unparsed/Malformed or e.g.:   Parsed OK, Unsaved,  Unparsed, Malformed etc.
 or: Parsed: OK/Error/Not yet, or Parsed: Yes/No/Error, Saved: Yes/No
-Try what happens when the editor content is not a valid sfz. The old instruemnt should be retained
 (or restored) and maybe we should somehow indicate a parsing error to the user
-Make Save button disappear (or grayed out) when the lastValidSfz in the SfzPlayer is not in sync 
 with the code editor content. We always save the lastValidSfz but the user will assume that we 
 save the code editor's content so we need to make sure that they are in sync when saving is 
 invoked. This is needed because the SfzPlayer that is out FileManager has no access to the 
 editor's content. We always save what we hear, not what we see - so to speak. That may be 
 confusing to the user when these two things are not the same (i.e. when the code editor content 
 has not yet been parsed). Maybe we should always parse? But no, that would try to parse 
 excessively often, namely after each character deletion/insertion. Moreover, it may often fail to
 parse because duringg editing, it is likely to encounter an intermediate malformed sfz. Maybe
 we should try to parse when the editor loses keyboard focus?
-add activeFileBecameDirty callback to FileManagerListener for subclasses to optionally override
-implement it in FileSelectionBox
-> this should make the astersik appear in the sfz-widget after parsing the editor content
But: the widget will not save the editor content but instead the most recently parsed *valid*
content...that may actually be a plus: it disallows saving invalid sfz files - but it's not what a
user would expect. How should we deal with that? Maybe the save button should disappear when the
editor content is not in sync with the lastvalidSfz? That may be good solution
-Total recall: Maybe we should have an xml-tag that defines whether the sfz text is referenced, 
 i.e. resides in its own file, or is textually embedded in the xml. Maybe as 
 SfzLocation="EmbeddedInThisXml"|"ReferencedFromFile". Or maybe we could have either a tag
 SfzFile or a tag SfzText. If both are present (which perhaps should actually not occurr), the 
 embedded text takes precedence over the referenced file

Bugs:
-For the modulation connections, we don't get any sliders to appear
 -write some sort of unit test for GUI interactions in the TestAppJURA and then use that to find 
  and fix this bug
-[fixed] Saw_1forAllFilterEnv.sfz: change cutoff -> assert! We fail to find the correct code segment in
 sfz document. -> make a unit test with an sfz patch content similar to the content of that patch.
 maybe it's because the patch start witha comment? Or is it because we now use a backslash? Try 
 some variations of the patch. Yes! I made a copy where I stripped the initial comment and with 
 that patch, it works! In the call to
   findOpcodeValueString(const std::string& code, int groupIndex, int regionIndex, ...
 is 1 - i think, it should be zero - also, the patch-structure in the tree shown an empty
 region - something is wrong in the creation of the patch tree. I think, this is the reason.
 ToDo: after creating the tree, compare it to the actual sfz structure and put an assert like
 jassert( treeViewMatchesPatch(...) );
 ...looks like the problem is not in creating the tree but already in creating the sfz 
 datastrcuture: it sometimes contains an empty region -> debug  SfzInstrument::setFromSFZ
 Oh! the problem seems to be not the comment but the initial empty line after the comment. If I
 remove the comment but leave a blank line in, the datastructure somehow creates a "ghost" region
 -> write a unit test that exposes this behavior!

-[fixed?] Change cutoff, then volume then cutoff again ...somewhere, there's still an update 
 missing SfzTreeViewNode.data.playbackSetting still has the old value. 

-When loading Saw_1forAllBiFilterEnv.sfz, we hit an assert in SfzCodeBook::modRoutingToString. It
 seems to be because when building the tree-view and we encounter a mod-routing fileg_depth that 
 goes to two filters, we add two tree-nodes. When building the tree, We need to somehow figure out,
 if we have already a node for a given mod-routing. The problem is that a single mod-routin opcode
 in the sfz can produce two or more mod-connections in the sfz datastructure...

-When switching sfz-files while playing, access violations occur. I think, I need to acquire locks
 in all the GUI functions ...or at least in parseCodeEditorContent
-When saving an sfz file, it seems liek the internal file list is not updated: switching through
 the presets with the forward/backward buttons seems to not load the newly saved file, as if it's 
 not there. Also, the file-widget is not updated, i.e. the new filename is not shown. Maybe these 
 two things have the same underlying reason?
 Maybe in FileManager::openSavingDialog, we need to call updateFileList after successful saving?
 -> check the behavior with xml files? Do we have the same problem there? Nope! There it works as
 it should! -> Check what steps it does and which of these are missing the case of sfz saving!
 hmm - OK: StateFileManager::saveStateToXmlFile actively calls updateFileList() towards the end, 
 so maybe we should aslo do that in SfzPlayer::saveToFile(). Or: move the call into 
 FileManager::openSavingDialog. Also check, if the overwriting code in 
 StateFileManager::saveStateToXmlFile is still needed. I think File::replaceText already implements
 the same functionality (create temp-file, then delete old file, then rename tempfile), so we may 
 not have to do that ourselves.
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
 rsSamplerEngine::handleNoteOn there's code that should do the rollback already - maybe it's buggy?

Performance issues:
-When the Play page is open and we move a controller (via hardware), the CPU load goes up to 12%.
 Strangely, when the Edit page is open, the CPU load does not go up as much - only to ~ 2%. It 
 seems to be caused by the slider updates in SamplerEditor::timerCallback. See comment there. In 
 debug mode, the CPU goes up to 25% with the slider updates active

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
 -The sliders on the "Play" page are shown in order of their controller number. Maybe try to show
  them in the order in which they appear in the <control> section or maybe in the order in which 
  they are assigned to controller DSP objects inside the engine (via the controlN_index opcodes). 
  But actually, via the order under <control> makes more sense from a user's point of view.
  It would also be nice if the initial value could become the default value available via 
  ctrl-click. 
 -Make a GUI page that shows all available patch parameters at once. Maybe call it the "Tweak"
  page. Or: maybe try to make the TreeView nodes themselves "draggable". Maybe we need to subclass
  TreeView to make TreeViewWithDraggableNodes or TreeViewWithNodeWidgets. Or: when the user holds
  shift while hovering over the tree, display an overlaid slider that can be directly manipulated.
  ...ok - this is in the works....
 -OR: allow the user to define custom GUIs either in the sfz itself (maybe under a a <gui> header)
  or in the .xml file
 -Maybe make use of:
    https://sfzformat.com/headers/control
    https://sfzformat.com/opcodes/label_ccN
    https://sfzformat.com/opcodes/set_ccN
  and use MIDI controllers to allow the user to control the patch. The controllers would become 
  available as modulation sources. Thereby, one controller could control multiple opcodes at once
  (like several cutoff frequencies - for different regions). Also, a single opcode could respond
  to multiple controllers (not sure, if that's useful though - maybe for implementing thinsg like
  reso-by-cutoff - resonance would be controlled by the resonance controller but also by cutoff.
  We would realize a many-to-many mapping of sliders to opcodes rather than a one-to-one mapping.
  Maybe allow the user to define sliders for the midi-controllers. A slider would have a position
  defined by x,y,w,h, a mode (horz, vert, rot), a mapping (linear, exponential, etc.), an 
  (optional) image (png, defining stripes or being rotated by code at runtime).
 -Maybe let the opcode editor show not only a single widget for the selected opcode but also 
  widgets for all "sibling" opcodes/parameters which are defined to be all the other the 
  parameters that apply to the same DSP effect. The widgets should be visible regardless whether 
  or not the opcode is specified in the sfz code. If it isn't we may need to take some special 
  care when searching for it - we won't find it. We should perhaps find the last defined sibling 
  and add the new opcode directly after it in the code (maybe on a new line). Maybe the 
  OpcodeEditor should become an EffectEditor
 -The Opcode-Widgets need to be cleared when a new sfz or xml is loaded
 -we need a TreeView that represents the SFZs structure: the root/top-level node represents the 
  whole instrument, subnodes represnet groups and subsubsnodes represent regions. Each of these 
  levels can have opcodes as leafs.
  -The leafs should be wigdges appropriate for the given opcode: sliders, comboboxes, textfields
  -It should be possible to add/remove nodes from the TreeView and edit the data of the leafs
  -The changes should be reflected in the code editor, too. And vice versa, changes to the code 
   should be immediately reflected in the TreeView
  -Maybe the CodeEditor and TreeView should communicate via a Mediator, as I do in Liberty. Let's
   have classes SfzCodeEditor, SfzTreeView, SfzEditingMediator. Maybe later we will also find other
   representations such as a piano-style region editor like in Kontakt.
 -After parsing the code, the sfz-field should get dirtified, maybe also already after editing the 
  code without parsing? After saving, it should be makred as clean again (i think, this should 
  happen automatically? Maybe next to the "Parse" button, there could be a Revert button that 
  reverts to ..yeah - to what? the last saved or the last valid parsed code? Both could be useful.
 -Maybe when the text is unparsed, indicate it in the TreeView - maybe make it empty? Or show some
  text like: "SFZ code was edited. Click "Parse" to update engine"
 -Maybe rename the "Edit" page to "Code" and let the actual Edit page show the tree (like in code) 
  on the left and then have also an area for keymapping, proper sample editing, envelope editors,
  filter editors, LFO editors, etc. Maybe there could also be a "View" page that somehwo shows
  the patch-structure in a more user friendly way...with or without editing capabilities. A proper
  sample editor, could have views for: waveform, spectrogram, spectrum (of whole file), scalogram,
  phasogram, views can swicth between L//R and M/S mode
 -FileManager::markFileAsClean just sets a boolean flag but doesn't rigger any callbacks, so the 
  sfz-widget does not acquire the "dirty" asterisk. Maybe let FileManager have a second callback
  activeFileWasDirtified or activeFileIsOutOfDate, activeFileBecameDirty or similar which the
  load/save widget then overrides
 -Maybe show the compatibility of the patch with sfz 1.0, 2.0, ARIA, RS-MET, etc. Maybe the code 
  editor could somehow be configured to highlight opcodes that require a higher standard. This is
  convenient for authoring sfz files. Maybe the engine could be configured to ignore opcodes of
  ceratin standards. Maybe also show some more status info in additional fields, like number of 
  groups, regions, samples, configuration (i.e. busMode or not)
 -SFZ editor functionality:
  -Dirtify the xml-state when the user loads a new sfz
  -Check dirtification of the .sfz while when the docmument has been edited
  -It would be nice, if the right-click menu could have options:
   -Add Opcode... presenting an (organized) list of opcodes that the user can insert
   -Create Slider... Cretaes a parameter slider for the opcode that is currently selected (there's
    no such selection feature yet, though - another thing to do)
 -Have some routable XY-Pad controllers
 -show some data about the loaded patch: number of samples, regions, groups, filters, equalizers, 
  waveshapers etc.
 -maybe display also, how many filters, eqs, etc. are allocated and how many are currently in use
 -[Done] In the left or right section, show a tree-view of the instrument with levels <global>, 
  <group>, <region> where <global> may be implicit, i.e. the top-level without an opening "folder"
  widget
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
  Maybe later add more views:
  -Play: performance parameters, patch meta-data, metering
  -Code: code editor with tree-view and slider
  -Edit: Visual key/vel mapping editor with piano roll, sample-edit window etc. (maybe call it Map)
  -View: 
  -Maybe some sort of destructive sample editor window ("Sample-Edit") with waveform-view, 
   spectrogram, spectrum, visual representation of sinusoidal model data, etc.

 
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
display the number of DSP modules...if that makes sense - but I think it doesn't because it depends
on the number of played notes...but maybe we can somehow compute the maximum number that could ever 
be needed? That could also help with allocating the objects...like MaxNumFilters, 
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
