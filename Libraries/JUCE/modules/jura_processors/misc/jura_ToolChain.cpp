AudioModule* AudioModuleFactory::createModule(const juce::String& type, CriticalSection *lock, 
  ModulationManager* modMan, MetaParameterManager* metaMan)
{
  // it is important that the strings used her match the moduleTypeName of the respective module
  // otherwise preset recall for ToolChain will not work

  if(type == "None")         return new DummyModule( lock);
  if(type == "DebugAudioModule")        return new DebugAudioModule( lock);

  // modulators:
  if(type == "BreakpointModulator") return new BreakpointModulatorAudioModule(lock);

  // analysis:
  if(type == "Scope")    return new PhaseScope(              lock);
  //if(type == "Scope2")  return new PhaseScope2( lock);
  if(type == "MultiAnalyzer") return new MultiAnalyzerAudioModule(lock);
  if(type == "TrackMeter")    return new TrackMeterAudioModule(   lock);
  if(type == "MidiMonitor")   return new MidiMonitorAudioModule(  lock);

  // generators:
  if(type == "Oscillator3D")   return new RotationOscillatorAudioModule(lock, metaMan, modMan);
  if(type == "RayBouncer")     return new RayBouncerAudioModule(lock);
  //if(type == "WaveOscillator") return new OscillatorStereoAudioModule(lock);
  //if(type == "FourOscSection") return new FourOscSectionAudioModule(lock);

  // effects:
  if(type == "Enveloper")        return new Enveloper(                  lock);
  if(type == "FuncShaper")       return new FuncShaperAudioModule(      lock, metaMan, modMan);
  if(type == "NodeShaper")       return new NodeShaperAudioModule(      lock);

  if(type == "AlgoVerb")         return new AlgoVerbAudioModule(        lock);
  if(type == "EchoLab")          return new EchoLabAudioModule(         lock);
  if(type == "PingPongEcho")     return new PingPongEchoAudioModule(    lock);
  if(type == "StereoDelay")      return new StereoDelayAudioModule(     lock);
  if(type == "PitchShifter")     return new PitchShifterAudioModule(    lock);
  if(type == "Quadrifex")        return new QuadrifexAudioModule(       lock, metaMan, modMan);
  if(type == "Moduluxury")       return new ModuluxuryAudioModule(      lock);
  if(type == "ChannelMatrix2x2") return new ChannelMatrix2x2AudioModule(lock);
  if(type == "DspWorkbench")     return new DspWorkbenchAudioModule(    lock);

  // filters:
  if(type == "Equalizer")       return new EqualizerAudioModule(      lock);
  if(type == "Ladder")          return new Ladder(                    lock, metaMan, modMan);
  if(type == "PhasorFilter")    return new PhasorFilter(              lock);
  if(type == "EngineersFilter") return new EngineersFilterAudioModule(lock);
  if(type == "CrossOver")       return new CrossOverAudioModule(      lock);

  // dynamics:
  if(type == "Limiter")         return new LimiterAudioModule(   lock);
  if(type == "MultiComp")        return new MultiCompAudioModule(lock, metaMan, modMan);

  // instruments:
  if(type == "AcidDevil")     return new AciDevilAudioModule(     lock);
  if(type == "Straightliner") return new StraightlinerAudioModule(lock);
  if(type == "NewSynth")      return new NewSynthAudioModule(     lock);
  if(type == "MagicCarpet")   return new MagicCarpetAudioModule(  lock);
  if(type == "SimpleSampler") return new SimpleSamplerAudioModule(lock);
  if(type == "KeyShot")       return new KeyShotAudioModule(      lock);
  if(type == "Quadriga")      return new QuadrigaAudioModule(     lock);
  if(type == "Workhorse")     return new WorkhorseAudioModule(    lock);
#ifdef _MSC_VER
  if(type == "Liberty")       return new LibertyAudioModule(      lock);
#endif

  jassertfalse;   // unknown module type requested
  //return nullptr;
  return new DummyModule(lock); // to avoid a crash when a user messes up an xml file
}

juce::String AudioModuleFactory::getModuleType(AudioModule *m)
{
  return m->getModuleTypeName();
    // replaces code below, needs testing, then old code can be deleted
    // ...actually, this function itself is superfluous now
}

//=================================================================================================

AudioModuleSelector::AudioModuleSelector() : RComboBox("ModuleSelector")
{
  // populate the tree:

  RTreeViewNode *node;
  int i = 1;           //  the index is actually not used, but we need it as dummy
  popUpMenu->addTreeNodeItem(new RTreeViewNode("None",     i++));
#if JUCE_DEBUG
  popUpMenu->addTreeNodeItem(new RTreeViewNode("DebugAudioModule",    i++));
#endif

  node = new RTreeViewNode("Instruments", -1, "Instruments");
  node->addChildNode(new RTreeViewNode("AcidDevil",       i++));
  node->addChildNode(new RTreeViewNode("Straightliner",  i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);

  node = new RTreeViewNode("Effects", -1, "Effects");
  //node->addChildNode(new RTreeViewNode("Enveloper",     i++));
  node->addChildNode(new RTreeViewNode("FuncShaper",    i++));
  node->addChildNode(new RTreeViewNode("NodeShaper",    i++));
  //node->addChildNode(new RTreeViewNode("StereoDelay",   i++)); // include in Quadrifex
  //node->addChildNode(new RTreeViewNode("PitchShifter",  i++)); // include in Quadrifex
  node->addChildNode(new RTreeViewNode("EchoLab",       i++));
  node->addChildNode(new RTreeViewNode("Quadrifex",     i++));
  node->addChildNode(new RTreeViewNode("PingPongEcho",  i++));
  //node->addChildNode(new RTreeViewNode("AlgoVerb",      i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);

  /*
  node = new RTreeViewNode("Sources", -1, "Sources");
  node->addChildNode(new RTreeViewNode("Oscillator3D",   i++));
  node->addChildNode(new RTreeViewNode("RayBouncer",     i++));
  node->addChildNode(new RTreeViewNode("WaveOscillator",  i++));  // 
  node->addChildNode(new RTreeViewNode("FourOscSection",  i++));
  //node->addChildNode(new RTreeViewNode("NoiseGenerator",  i++));
  //node->addChildNode(new RTreeViewNode("SamplePlayer",    i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);
  */

  node = new RTreeViewNode("Filters", -1, "Filters");
  node->addChildNode(new RTreeViewNode("Ladder",          i++));
  node->addChildNode(new RTreeViewNode("Equalizer",       i++));
  node->addChildNode(new RTreeViewNode("EngineersFilter", i++));
  //node->addChildNode(new RTreeViewNode("PhasorFilter",    i++));
  //node->addChildNode(new RTreeViewNode("CrossOver",       i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);

  node = new RTreeViewNode("Dynamics", -1, "Dynamics");
  //node->addChildNode(new RTreeViewNode("MultiComp",     i++));
  node->addChildNode(new RTreeViewNode("Limiter",     i++));
  popUpMenu->addTreeNodeItem(node);

  node = new RTreeViewNode("Modulators", -1, "Modulators");
  node->addChildNode(new RTreeViewNode("BreakpointModulator",  i++));
  //node->addChildNode(new RTreeViewNode("LowFrequencyOscillator",  i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);

  node = new RTreeViewNode("Analyzers", -1, "Analyzers");
  node->addChildNode(new RTreeViewNode("Scope",    i++));
  //node->addChildNode(new RTreeViewNode("PhaseScope2",   i++));
  node->addChildNode(new RTreeViewNode("MultiAnalyzer", i++));
  //node->addChildNode(new RTreeViewNode("TrackMeter",    i++));
  node->addChildNode(new RTreeViewNode("MidiMonitor",   i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);

  bool showUnfinishedModules = true;
  if(showUnfinishedModules)
  {
    node = new RTreeViewNode("UnderConstruction", -1, "UnderConstruction");
#ifdef _MSC_VER
    node->addChildNode(new RTreeViewNode("Liberty", i++));     // not yet available on gcc
#endif
    node->addChildNode(new RTreeViewNode("NewSynth",  i++));
    //node->addChildNode(new RTreeViewNode("MagicCarpet",    i++));
    node->addChildNode(new RTreeViewNode("SimpleSampler",  i++));
    //node->addChildNode(new RTreeViewNode("KeyShot",        i++));
    //node->addChildNode(new RTreeViewNode("Quadriga",       i++));
    //node->addChildNode(new RTreeViewNode("Workhorse",      i++));
    node->setOpen(false);
    popUpMenu->addTreeNodeItem(node);
  }

  //setSize(300, 300); // has no effect
}

// The current release version includes:
// Instruments:
//   AcidDevil
//   Straightliner
// Effects:
//   FuncShaper
//   EchoLab
// Filters:
//   Equalizer
//   EngineersFilter
// Dynamics:
//   Limiter
// Analyzers:
//   Scope
//   MultiAnalyzer
//   MidiMonitor

void AudioModuleSelector::drawHighlighted(bool shouldBeHighlighted)
{
  highlighted = shouldBeHighlighted;
  repaint();
}

void AudioModuleSelector::paint(Graphics& g)
{
  if(!highlighted)
    RComboBox::paint(g);
  else
  {
    g.fillAll(getHandleColour()); // only difference to baseclass version - maybe refactor
    g.setColour(getOutlineColour());
    g.drawRect(0, 0, getWidth(), getHeight(), 2);
    int x = 4;
    int y = getHeight()/2 - font->getFontAscent()/2;
    drawBitmapFontText(g, x, y, getSelectedItemText(), font, getTextColour());
  }
}

//=================================================================================================

ToolChain::ToolChain(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse) 
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse/*, &modManager*/) // passing modManager causes access violation (not yet constructed)?
  , modManager(lockToUse) // maybe pass the metaManagerToUse to this constructor call
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("ToolChain");
  modManager.setMetaParameterManager(metaManagerToUse);
  setModulationManager(&modManager);
  //createDebugModSourcesAndTargets(); // for debugging the mod-system
  addEmptySlot();
}

ToolChain::~ToolChain()
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++)
    delete modules[i];
}

void ToolChain::addEmptySlot()
{
  addModule("None");
}

bool ToolChain::addModule(const juce::String& type)
{
  ScopedLock scopedLock(*lock);
  AudioModule *m = AudioModuleFactory::createModule(type, lock, &modManager, metaParamManager); // todo: pass the metaParamManager too
  if(m){
    addModule(m);
    return true; 
  }
  return false;
}

void ToolChain::addModule(AudioModule* m)
{
  ScopedLock scopedLock(*lock);
  jassert(m != nullptr);
  setupManagers(m);
  append(modules, m);
  //m->setModuleName("Slot" + String(size(modules)) + "-" + type);
  m->setModuleName("Slot" + String(size(modules)) + "-" + m->getModuleTypeName());
  addToModulatorsIfApplicable(m);
  sendAudioModuleWasAddedNotification(m, size(modules)-1);
}

void ToolChain::deleteModule(int index)
{
  ScopedLock scopedLock(*lock);
  jassert(index >= 0 && index < size(modules)); // index out of range
  if(activeSlot == index)
    activeSlot--;
  sendAudioModuleWillBeDeletedNotification(modules[index], index);
  removeFromModulatorsIfApplicable(modules[index]);
  delete modules[index];
  remove(modules, index);
}

void ToolChain::deleteLastModule()
{
  ScopedLock scopedLock(*lock);
  deleteModule(size(modules) - 1);
}

void ToolChain::replaceModule(int index, const juce::String& type)
{
  ScopedLock scopedLock(*lock);
  jassert(index >= 0 && index < size(modules)); // index out of range
  if(!isModuleOfType(index, type)){              // replace only, if new type is different
    AudioModule* oldModule = modules[index];
    AudioModule* newModule = AudioModuleFactory::createModule(type, lock, &modManager, metaParamManager);
    newModule->setModuleName("Slot" + String(index+1) + "-" + type);

    setupManagers(newModule);
    //newModule->setSmoothingManager(smoothingManager);
    //newModule->setMetaParameterManager(metaParamManager); // superfluous now?

    //newModule->loadDefaultPreset(); // later: either load default preset or recall a stored state
    newModule->setSampleRate(sampleRate);
    modules[index] = newModule;
    removeFromModulatorsIfApplicable(oldModule);
    addToModulatorsIfApplicable(newModule);
    sendAudioModuleWasReplacedNotification(oldModule, newModule, index);
    delete oldModule;
    activeSlot = index;
    ensureOneEmptySlotAtEnd();
  }
}

bool ToolChain::isModuleOfType(int index, const juce::String& type)
{
  ScopedLock scopedLock(*lock);
  jassert(index >= 0 && index < size(modules)); // index out of range
  return type == AudioModuleFactory::getModuleType(modules[index]);
}

AudioModule* ToolChain::getModuleAt(int index)
{
  if(index < 0 || index >= size(modules))  // no assert, this is supposed to happen
    return nullptr;
  return modules[index];
}

void ToolChain::ensureOneEmptySlotAtEnd()
{
  ScopedLock scopedLock(*lock);

  // if the last module/slot is not the "empty" dummy module, add another slot:
  if(!isModuleOfType(size(modules)-1, "None"))
    addEmptySlot();

  // remove superfluous empty slots at end:
  while(modules.size() > 1 && isModuleOfType(size(modules)-1, "None")
                           && isModuleOfType(size(modules)-2, "None"))
  { // if the last two slots are empty, remove the last
    deleteLastModule();
  }
}

void ToolChain::addToolChainObserver(ToolChainObserver *observerToAdd)
{
  ScopedLock scopedLock(*lock);
  appendIfNotAlreadyThere(observers, observerToAdd);
}

void ToolChain::removeToolChainObserver(ToolChainObserver *observerToRemove)
{
  ScopedLock scopedLock(*lock);
  removeFirstOccurrence(observers, observerToRemove);
}

void ToolChain::sendAudioModuleWasAddedNotification(AudioModule *module, int index)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(observers); i++)
    observers[i]->audioModuleWasAdded(this, module, index);
}

void ToolChain::sendAudioModuleWillBeDeletedNotification(AudioModule *module, int index)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(observers); i++)
    observers[i]->audioModuleWillBeDeleted(this, module, index);
}

void ToolChain::sendAudioModuleWasReplacedNotification(AudioModule *oldModule,
  AudioModule *newModule, int index)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(observers); i++)
    observers[i]->audioModuleWasReplaced(this, oldModule, newModule, index);
}

// overrides:

AudioModuleEditor* ToolChain::createEditor()
{
  return new ToolChainEditor(this);
}

void ToolChain::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  if(numChannels != 2)
    return;

  //ScopedLock scopedLock(*lock); already held by wrapping plugin

  bool needsSmoothing  = smoothingManager->needsSmoothing();
  bool needsModulation = modManager.getNumConnections() > 0;

  if( !needsSmoothing && !needsModulation )
  {
    // in case of no modulations, we can use a faster loop
    for(int i = 0; i < size(modules); i++)
      modules[i]->processBlock(inOutBuffer, numChannels, numSamples);
  }
  else
  {
    // we have to iterate through all the samples and for each sample update all the modulators and
    // then compute a sample-frame from each non-modulator module
    for(int n = 0; n < numSamples; n++)
    {
      if(needsSmoothing)
      {
        //smoothingManager->updateSmoothedValues(); // nope, the lock is not the problem
        smoothingManager->updateSmoothedValuesNoLock();
      }
      if(needsModulation)
        modManager.applyModulationsNoLock();
      for(int i = 0; i < size(modules); i++)
      {
        //// for debug:
        //double left  = inOutBuffer[0][n];
        //double right = inOutBuffer[1][n];

        modules[i]->processStereoFrame(&inOutBuffer[0][n], &inOutBuffer[1][n]);
        // AudioModules that are subclasses of ModulationSource have not overriden this function.
        // That means, they inherit the empty baseclass method and do nothing in this call.
      }
    }
  }
}

void ToolChain::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*lock);
  sampleRate = newSampleRate;
  for(int i = 0; i < size(modules); i++)
    modules[i]->setSampleRate(sampleRate);
}

void ToolChain::handleMidiMessage(MidiMessage message)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->handleMidiMessage(message); }

  // todo: maybe let different slots receive MIDI on different channels
  // and/or don't override the noteOn/etc. functions here but rather let the MIDI events also
  // pass through the modules in series. most modules just pass them through, but we can also
  // have MIDI effects such as appregiators and sequencers which modify the sequence and pass
  // the modified sequence to the next module - we could have an appregiator in front of a
  // synth, for example

  // all synthesizer modules should pass through the incoming audio and add their own signal
  // (unless the use it inside for their own signal processing) -> this allows for layering
}

/*
// obsolete, now that we have overriden handleMidiMessage - test if the modulation system works, if
// so, these can eventually be deleted:
void ToolChain::noteOn(int noteNumber, int velocity)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->noteOn(noteNumber, velocity);
  }
}

void ToolChain::noteOff(int noteNumber)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->noteOff(noteNumber);
  }
}

void ToolChain::setMidiController(int controllerNumber, float controllerValue)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->setMidiController(controllerNumber, controllerValue);
  }
}

void ToolChain::setPitchBend(int pitchBendValue)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->setPitchBend(pitchBendValue);
  }
}
*/

void ToolChain::reset()
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++)
    modules[i]->reset();
}

XmlElement* ToolChain::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*lock);
  XmlElement *xml = AudioModule::getStateAsXml(stateName, markAsClean); 
  xml->setAttribute("ActiveSlot", activeSlot+1);
  for(int i = 0; i < size(modules); i++){
    juce::String typeString = AudioModuleFactory::getModuleType(modules[i]);
    XmlElement *child = new XmlElement("Slot");
    child->setAttribute("Type", typeString);
    //child->setAttribute("Bypass", isSlotBypassed(i)); // add later
    child->addChildElement(modules[i]->getStateAsXml(typeString, markAsClean));
    xml->addChildElement(child);
  }

  xml->addChildElement(modManager.getStateAsXml()); 
  // do this also in ModulatableAudioModule...wait no - the mod-settings are only stored in the 
  // top-level module, maybe we should have a ModManagerAudioModule as baseclass which contains the
  // ModulationManager. then, instead of calling xml = AudioModule::getStateAsXml(stateName, markAsClean); 
  // we would call it from the ModManagerAudioModule and this would then already include the mod
  // settings

  return xml;
}

void ToolChain::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  bool markAsClean)
{
  ScopedLock scopedLock(*lock);

  bool smoothingIsByassed = smoothingManager->isSmoothingBypassed();
  smoothingManager->setBypassSmoothing(true);

  //AudioModule::setStateFromXml(xmlState, stateName, markAsClean); 
  // the Chainer has no (global) parameters of its own, nor any static child-modules, so we don't
  // need to call the baseclass method here (it wouldn't do anything)

  recallSlotsFromXml(xmlState, markAsClean);
  recallModulationsFromXml(xmlState); // this should also recall meta-assignments for mod-depths

  // new:
  //recallMidiMappingFromXml(xmlState);
  //recallMetaMappingFromXml(xmlState); // there
  recallMetaValuesFromXml(xmlState);

  smoothingManager->setBypassSmoothing(smoothingIsByassed); // restore old bypassstate
}

void ToolChain::recallSlotsFromXml(const XmlElement &xmlState, bool markAsClean)
{
  activeSlot = -1;  // i think, that's not necessary - should be already -1
  int tmpActiveSlot = xmlState.getIntAttribute("ActiveSlot", 1) - 1;
  clearModulesArray();
  int i = 0;
  forEachXmlChildElementWithTagName(xmlState, slotState, "Slot"){
    juce::String type = slotState->getStringAttribute("Type");
    if(i == tmpActiveSlot)        // hack: we set it before adding the module, so the editor
      activeSlot = tmpActiveSlot; // retrieves the correct value in the moduleAdded callback

    //// old:
    //if(addModule(type))
    //{
    //  XmlElement *moduleState = slotState->getChildElement(0);
    //  modules[i]->setStateFromXml(*moduleState, "", markAsClean);
    //  i++;
    //}
    //// it may be better to first create the module, then recall its state and then add it (instead
    //// of recalling the state after adding) because when loading a preset, the editor my show
    //// a module in the initial state instead of the actual state - check with EchoLab - load
    //// ToolChain presets - the EchoLab editro state is wrong after recall

    // new:
    AudioModule *m = AudioModuleFactory::createModule(type, lock, &modManager, metaParamManager);
    if(m != nullptr)
    {
      //m->setMetaParameterManager(metaParamManager);      // so, the meta-mapping gets recalled
      setupManagers(m);
      XmlElement *moduleState = slotState->getChildElement(0);
      m->setStateFromXml(*moduleState, "", markAsClean); // set the state of the module before...
      addModule(m);                                      // ...adding it, so the newly created
      i++;                                               // editor has correct initial state
    }
  }
}

// move to ModulatableAudioModule:
void ToolChain::recallModulationsFromXml(const XmlElement &xmlState)
{
  modManager.removeAllConnections();
  XmlElement* modXml = xmlState.getChildByName("Modulations");
  if(modXml != nullptr)
    modManager.setStateFromXml(*modXml);  // recall modulation settings
}

void ToolChain::setupManagers(AudioModule* m)
{
  m->setSmoothingManager(smoothingManager);
  m->setMetaParameterManager(metaParamManager); 
  ModulatableAudioModule* mm = dynamic_cast<ModulatableAudioModule*>(m);
  if(mm)
    mm->setModulationManager(&modManager);
}

void ToolChain::addToModulatorsIfApplicable(AudioModule* module)
{
  ModulationSource* ms = dynamic_cast<ModulationSource*> (module);
  if(ms != nullptr)
  {
    assignModulationSourceName(ms);
    modManager.registerModulationSource(ms);
  }
}

void ToolChain::removeFromModulatorsIfApplicable(AudioModule* module)
{
  ModulationSource* ms = dynamic_cast<ModulationSource*> (module);
  if(ms != nullptr)
    modManager.deRegisterModulationSource(ms);
}

void ToolChain::assignModulationSourceName(ModulationSource* source)
{
  juce::String name;

  //// old - sources may be named in different index order after state recall (not good):
  //if(dynamic_cast<BreakpointModulatorAudioModule*> (source))
  //  name = "BM ";
  //// else if...
  //name += String(numRegisteredSourcesOfType(source) + 1);
  //source->setModulationSourceName(name);

  // new:
  AudioModule* am = dynamic_cast<AudioModule*> (source);
  if(am != nullptr)
  {
    int slotIndex = find(modules, am);
    jassert(slotIndex >= 0); // something is wrong - the source is not in the modules array
    name = "Slot" + String(slotIndex+1) + String("-") + AudioModuleFactory::getModuleType(am);
    source->setModulationSourceName(name);
  }
}

void ToolChain::clearModulesArray()
{
  ScopedLock scopedLock(*lock);
  while(size(modules) > 0)
    deleteLastModule();
}

void ToolChain::createDebugModSourcesAndTargets()
{
  // This code was for debugging only and can eventually be thrown away...


  // I think, this is the minimal code that triggers the bug that is apparent in Elan's 
  // SpiralGenerator:

  jura::BreakpointModulatorAudioModule* env = 
    new jura::BreakpointModulatorAudioModule(lock);
  env->setModuleName("Envelope1");
  addChildAudioModule(env);

  modManager.registerModulationSource(env);

  ModulatableParameter* p = 
    new ModulatableParameter("Gain", -60.0, -20.0, 0.0, Parameter::LINEAR, 0.01);
  addObservedParameter(p);
  //p->setValueChangeCallback<ToolChain>(this, &ToolChain::setGain);
  // try to use a lambda function like Elan does

  // 2 bugs:
  // 1: p has no modManager pointer (i.e. null)
  //    -> solved by calling setModulationManager in constructor
  // 2: on destruction, we hit ModulatableParameter.cpp, line 322:
  //    jassert(contains(availableSources, source)); // source was never registered
  //    -> solved by de-registering all sources (and targets) in dtor of ModulationManager

  // ...oh...the modManager's metaManager pointer is still null here - why?
  // -> solved by calling setMetaParameterManager in the ctor

  // int dummy = 0;
}

//=================================================================================================

ToolChainEditor::ToolChainEditor(jura::ToolChain *moduleChainToEdit)
  : AudioModuleEditor(moduleChainToEdit)
{
  ScopedLock scopedLock(*lock);
  chain = moduleChainToEdit;
  setHeadlinePosition(TOP_LEFT);
  numHueOffsets = 2;
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  updateEditorArray();
  updateSelectorArray();
  updateActiveEditor();
  chain->addToolChainObserver(this);
  addChangeListener(this); // we listen to ourselves for deferred destruction of selectors

  //setWidgetAppearance(ColourScheme::DARK_ON_BRIGHT); // bug: doesn't affect state-widgets
}

ToolChainEditor::~ToolChainEditor()
{
  ScopedLock scopedLock(*lock);
  chain->removeToolChainObserver(this);
  clearEditorArray();
}

AudioModuleEditor* ToolChainEditor::getEditorForSlot(int index)
{
  ScopedLock scopedLock(*lock);
  if(size(editors) == 0 || index < 0)            // may happen during xml state recall
    return nullptr;
  jassert(index >= 0 && index < size(editors)); // index out of range
  if(editors[index] == nullptr)
    editors[index] = chain->modules[index]->createEditor();
  return editors[index];
}

void ToolChainEditor::replaceModule(int index, const juce::String& type)
{
  ScopedLock scopedLock(*lock);
  jassert(index >= 0 && index < size(editors));  // index out of range
  if(!chain->isModuleOfType(index, type)){
    chain->replaceModule(index, type);            // will call audioModuleWillBeDeleted
    //AudioModule* m = chain->getModuleAt(index);   // can be 0, if dummy module was placed at end
    updateEditorArray();
    index = chain->activeSlot;
    editors[index] = getEditorForSlot(index);
    updateActiveEditor();
    scheduleSelectorArrayUpdate();                // deferred call to updateSelectorArray
                                                  // may be superfluous now
  }
}

void ToolChainEditor::updateSelectorArray()
{
  ScopedLock scopedLock(*lock);
  int numModules   = size(chain->modules);
  int numSelectors = size(selectors);
  AudioModuleSelector *s;

  // remove superfluous selectors:
  while(numSelectors > numModules){
    s = selectors[numSelectors-1];
    removeWidget(s, true, true);
    remove(selectors, numSelectors-1);
    numSelectors--;
  }

  // add required selectors:
  while(numModules > numSelectors){
    s = new AudioModuleSelector();
    s->setInterceptsMouseClicks(false, false); // we handle them 1st and possibly pass them through
    s->selectItemFromText(
      AudioModuleFactory::getModuleType(chain->modules[numSelectors]), false);
    s->registerComboBoxObserver(this);
    s->setDescriptionField(infoField); // somehow, this doesn't work
    addWidget(s);
    append(selectors, s);
    numSelectors++;
  }

  //// old:
  //// without it, the selector for 1st slot is wrong after preset loading:
  //if(numSelectors > 0 )
  //  selectors[0]->selectItemFromText(AudioModuleFactory::getModuleType(chain->modules[0]), false);

  // let selectors reflect the selected module type and selector for highlight active slot
  for(int i = 0; i < numSelectors; i++)
    selectors[i]->selectItemFromText(AudioModuleFactory::getModuleType(chain->modules[i]), false);
  updateActiveSelector();
}

void ToolChainEditor::updateEditorArray()
{
  ScopedLock scopedLock(*lock);
  int numModules = size(chain->modules);
  int numEditors = size(editors);
  AudioModuleEditor *e;

  // remove superfluous editors:
  while(numEditors > numModules){
    e = editors[numEditors-1];
    if(e == activeEditor)
      activeEditor = nullptr;
    delete e;
    remove(editors, numEditors-1);
    numEditors--;
  }

  // add placeholders for required selectors:
  while(numModules > numEditors){
    append(editors, (AudioModuleEditor*)nullptr);
    numEditors++;
  }
}

void ToolChainEditor::updateActiveEditor()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor* tmpEditor = getEditorForActiveSlot();
  if(tmpEditor != activeEditor){
    removeChildEditor(activeEditor, false);
    addChildEditor(tmpEditor);
    activeEditor = tmpEditor;
    int w = max(240, activeEditor->getWidth());
    int h = max(180, activeEditor->getHeight());
    activeEditor->setBounds(leftColumnWidth, 0, w, h);
    setSize(w + leftColumnWidth, h + bottomRowHeight);
    resized(); // setSize will call resized only if the size actually changes but we need to make
               // sure that it always gets called to arrange the selectors
  }
}

void ToolChainEditor::updateActiveSelector()
{
  for(int i = 0; i < size(selectors); i++){
    if(i == chain->activeSlot)
      selectors[i]->drawHighlighted(true);
    else
      selectors[i]->drawHighlighted(false);
  }
}

void ToolChainEditor::mouseDown(const MouseEvent &e)
{
  //ScopedLock scopedLock(*lock); // blocks audio when popup is open
  bool wasHandled = false;
  int i = chain->activeSlot;
  juce::Rectangle<int> rect = selectors[i]->getBounds();
  if(rect.contains(e.x, e.y)){
    // click was on active slot selector - pass event through:
    selectors[i]->mouseDown(e.getEventRelativeTo(selectors[i]));
    wasHandled = true;
  }
  else{
    for(i = 0; i < size(selectors); i++){
      rect = selectors[i]->getBounds();
      if(rect.contains(e.x, e.y)){
        // click was on inactive slot selector - activate:
        chain->activeSlot = i;
        updateActiveEditor();
        updateActiveSelector();
        repaint();
        wasHandled = true;
      }
    }
  }
  if(wasHandled == false)
    AudioModuleEditor::mouseDown(e);
}

void ToolChainEditor::resized()
{
  ScopedLock scopedLock(*lock);

  Editor::resized();

  int x, y, w, h, dy, margin;
  margin = 4;
  x = margin;
  y = getHeadlineBottom() + margin;
  w = leftColumnWidth - x - margin;
  h = 16;
  stateWidgetSet->setBounds(x, y, w, 32);

  // arrange selectors:
  y  = getPresetSectionBottom() + margin;
  dy = h-2;
  for(int i = 0; i < size(selectors); i++){
    selectors[i]->setBounds(x, y, w, h);
    y += dy;
  }

  // set up bounds of the editor for the active module:
  if(activeEditor != nullptr){
    y = 0;
    x = leftColumnWidth;
    w = getWidth()  - leftColumnWidth;
    h = getHeight() - bottomRowHeight;
    activeEditor->setBounds(x, y, w, h);
  }

  // set up bounds for inherited widgets:
  x = margin;
  y = getHeight() - bottomRowHeight;
  w = getWidth() - x - margin;
  infoField->setBounds(x, y, w, bottomRowHeight);

  x = getWidth() - 108;
  w = getWidth() - x - margin;
  webLink->setBounds(x, y+3, w, bottomRowHeight);

  int buttonWidth = 40;
  x = leftColumnWidth - buttonWidth - margin;
  y = margin;
  setupButton->setBounds(x, y, buttonWidth, 16);

  // If this is a ToolChain wrapped into an AudioPlugIn, we want to resize the whole parent
  // window as well:
  Component *parent =	getParentComponent();
  if(dynamic_cast<AudioPluginEditor*>(parent))
    parent->setSize(getWidth(), getHeight());
}

void ToolChainEditor::paintOverChildren(Graphics& g)
{
  // highlight active slot by drawing a rectangle around it:
  ScopedLock scopedLock(*lock);
  if(size(selectors) == 0)   // occurs during state recall
    return;

  g.setColour(Colour::fromFloatRGBA(0.8125f, 0.8125f, 0.8125f, 1.f));// maybe switch depending on widget color-scheme (dark-on-bright vs bright-on-dark)
  juce::Rectangle<int> rect = selectors[chain->activeSlot]->getBounds();
  g.drawRect(rect, 2);  // 2nd param: thickness
}

void ToolChainEditor::rComboBoxChanged(RComboBox* box)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(selectors); i++){
    if(box == selectors[i]){
      replaceModule(i, box->getSelectedItemText());
    }
  }
}

void ToolChainEditor::changeListenerCallback(ChangeBroadcaster *source)
{
  ScopedLock scopedLock(*lock);
  if(source == this)
  {
    updateSelectorArray();
    resized();  // to arrange selectors
  }
  else
    AudioModuleEditor::changeListenerCallback(source);
}

void ToolChainEditor::audioModuleWasAdded(ToolChain *chain,
  AudioModule *module, int index)
{
  ScopedLock scopedLock(*lock);
  updateEditorArray();
  updateActiveEditor();
  updateSelectorArray();
}

void ToolChainEditor::audioModuleWillBeDeleted(ToolChain *chain,
  AudioModule *module, int index)
{
  ScopedLock scopedLock(*lock);
  deleteEditor(index);
  scheduleSelectorArrayUpdate();
}

void ToolChainEditor::audioModuleWasReplaced(ToolChain *chain,
  AudioModule *oldModule, AudioModule *newModule, int index)
{
  ScopedLock scopedLock(*lock);
  deleteEditor(index);
}

void ToolChainEditor::scheduleSelectorArrayUpdate()
{
  sendChangeMessage();
  // we will receive the message ourselves which causes a call to updateSelectorArray()
}

void ToolChainEditor::deleteEditor(int index)
{
  ScopedLock scopedLock(*lock);
  jassert(index >= 0 && index < size(editors)); // index out of range
  if(activeEditor == editors[index])
    activeEditor = nullptr;
  removeChildColourSchemeComponent(editors[index], true);
  editors[index] = nullptr;
}

void ToolChainEditor::clearEditorArray()
{
  ScopedLock scopedLock(*lock);
  activeEditor = nullptr;
  for(int i = 0; i < size(editors); i++)
    removeChildColourSchemeComponent(editors[i], true);
  editors.clear();
}
