class JUCE_API DummyModule : public jura::AudioModule
{
public:
  DummyModule(CriticalSection *lockToUse) : AudioModule(lockToUse) 
  {
    setModuleTypeName("None");
  }
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override 
  {
    //// for debug:
    //std::vector<double> left(numSamples), right(numSamples);
    //for(int n = 0; n < numSamples; n++)
    //{
    //  left[n]  = inOutBuffer[0][n];
    //  right[n] = inOutBuffer[1][n];
    //}
    //int dummy = 0;
  }
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(DummyModule)
};

//=================================================================================================

ToolChain::ToolChain(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse) 
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse/*, &modManager*/) // passing modManager causes access violation (not yet constructed)?
  /*, modManager(lockToUse)*/, moduleFactory(lockToUse) // maybe pass the metaManagerToUse to this constructor call
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("ToolChain");
  modManager = new ModulationManagerPoly(lockToUse);
  modManager->setMetaParameterManager(metaManagerToUse);
  setModulationManager(modManager);
  modManager->setVoiceManager(&voiceManager);
  voiceSignals.resize(2 * voiceManager.getMaxNumVoices()); // maybe move into allocateVoiceResources later
  voiceManager.setVoiceSignalBuffer(&voiceSignals[0]);
  createMidiModSources();
  //createDebugModSourcesAndTargets(); // for debugging the mod-system
  populateModuleFactory();
  addEmptySlot();
}

ToolChain::~ToolChain()
{
  ScopedLock scopedLock(*lock);

  // Trying to fix the crash of Elan's SeToolChain on destruction:
  modManager->deRegisterAllTargets();
  modManager->deRegisterAllSources();
  // Yep - that seems to fix it. I think the problem was this: when our modManager member goes out
  // of scope, it calls these 2 functions in its destructor - and in these functions the pointers
  // to the source/target modules are referenced (for nulling their modManager pointer) - but at
  // this point, the source/target modules were already deleted (in the loop below, which is 
  // executed before the modManager goes out of scope). But i wonder, why i don't have the 
  // same crash in my ToolChain. hmm....

  for(int i = 0; i < size(modules); i++)
    delete modules[i];

  modManager->setVoiceManager(nullptr); // Dunno if required but better safe than sorry. If the 
  // voiceManager is destructed before the modManager, we might otherwise have a dangling pointer
  // in the modManager for a brief moment during destruction. That may be inconsequential but who 
  // knows... But maybe we can make sure that the modManager is destructed first by declaration 
  // order -> figure out

  deleteMidiModSources();

  delete modManager;
}

void ToolChain::addEmptySlot()
{
  addModule("None");
}

bool ToolChain::addModule(const juce::String& type)
{
  ScopedLock scopedLock(*lock);
  //AudioModule *m = AudioModuleFactory::createModule(type, lock, &modManager, metaParamManager); // todo: pass the metaParamManager too
  AudioModule *m = moduleFactory.createModule(type); // maybe we need to set up the managers?
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

  // handle case when 0th module was deleted (needs test):
  if(activeSlot == -1 && size(modules) > 0 )
    activeSlot = 0;
  // ...in ToolChain, it's impossible to have a completely empty module list, but subclasses (like 
  // Elan's) may want to allow this, so we need to support it
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
    //AudioModule* newModule = AudioModuleFactory::createModule(type, lock, &modManager, metaParamManager);
    AudioModule* newModule =  moduleFactory.createModule(type); // maybe we need to set up the managers?
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
  return type == modules[index]->getModuleTypeName();
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

AudioModuleEditor* ToolChain::createEditor(int type)
{
  return new ToolChainEditor(this);
}

void ToolChain::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  if(numChannels != 2) return;
  //ScopedLock scopedLock(*lock); // lock already held by the wrapping plugin
  bool needsSmoothing  = smoothingManager->needsSmoothing();
  bool needsModulation = modManager->getNumConnections() > 0;
  bool needsVoiceKill  = voiceManager.needsVoiceKillCheck();

  if( !needsSmoothing && !needsModulation && !needsVoiceKill )
    for(size_t i = 0; i < modules.size(); i++)
      modules[i]->processBlock(inOutBuffer, numChannels, numSamples);
  else {
    // we have to iterate through all the samples and for each sample, update all the modulators 
    // and then compute a sample-frame from each non-modulator module:
    for(int n = 0; n < numSamples; n++) 
    {
      if(needsSmoothing)   
        smoothingManager->updateSmoothedValuesNoLock();
      if(needsModulation)  
        modManager->applyModulationsNoLock();
      for(size_t i = 0; i < modules.size(); i++)
        modules[i]->processStereoFrame(&inOutBuffer[0][n], &inOutBuffer[1][n]);
        // AudioModules that are subclasses of ModulationSource have not overriden this function.
        // That means, they inherit the empty baseclass method and do nothing in this call.

      if(needsVoiceKill)
      {
        voiceManager.findAndKillFinishedVoices();
        needsVoiceKill = voiceManager.needsVoiceKillCheck(); // condition may have changed
      }
      // should we also equip rsVoiceManager with a pointer to a mutex and generally lock it and
      // provide a "NoLock" version of these functions? or is it safe to use without?

    }
  }

  // todo: all instrument and source modules should pass through the incoming audio and add their 
  // own signal (unless the use it inside for their own signal processing) -> allows for layering
}

void ToolChain::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*lock);
  sampleRate = newSampleRate;
  voiceManager.setSampleRate(sampleRate);
  for(int i = 0; i < size(modules); i++)
    modules[i]->setSampleRate(sampleRate);
}

void ToolChain::handleMidiMessage(MidiMessage message)
{
  ScopedLock scopedLock(*lock);

  if(message.getChannel() != 1) return;
  // We currently respond only to messages on channel 1. This is a provision for later allowing
  // different slots respond to different channels. I guess it may break behavior, if the user
  // has a saved project which sends on another channel and a later version responds differently
  // ...but maybe not...but better safe than sorry - we now force the user to use channel 1 for 
  // better project compatibility with future versions. Maybe for each slot, we should have a 
  // set of flags 1..16 and responds to messages on all channels which have the flag set - use 
  // class RAPT::rsFlages16


  //rsMidiMessageHandler::MidiHandleInfo info;   // filled out by the voiceManager
  //voiceManager.handleMidiMessage(message, &info);
  // maybe the voiceManager should return some information, specifically, which voice was assigned
  // and we may have to pass this information to the child modules in the llop below - maybe we 
  // need to introduce a new callback handleMidiMessage(const MidiMessage&, int voice)

  int voice = voiceManager.handleMidiMessageReturnVoice(message);
  if(voice < 0) return;
  // If no voice was allocated or used for the event, we don't pass it on to the child modules to 
  // free them from having to handle this condition themselves. So the child modules can expect 
  // that they always receive a valid voice index in their per voice callbacks. This simplifies 
  // their implementations of which we will have many, so it's worth it.

  for(int i = 0; i < size(modules); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->handleMidiMessageForVoice(message, voice); } 

  //voiceManager.handleMidiMessage(message);
  //for(int i = 0; i < size(modules); i++){
  //  AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
  //  if(m != nullptr)
  //    m->handleMidiMessage(message); }   // todo: pass info
  //    // the child modules inquire the voiceIndex from the info

  // todo: maybe let different slots receive MIDI on different channels, maybe 
  // AudioModuleWithMidiIn should have a means to filter midi messages based on their channel and
  // respond only, if the message passes the filter


  // and/or don't override the noteOn/etc. functions here but rather let the MIDI events also
  // pass through the modules in series. most modules just pass them through, but we can also
  // have MIDI effects such as appregiators and sequencers which modify the sequence and pass
  // the modified sequence to the next module - we could have an appregiator in front of a
  // synth, for example
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
    juce::String typeString = modules[i]->getModuleTypeName();
    XmlElement *child = new XmlElement("Slot");
    child->setAttribute("Type", typeString);
    //child->setAttribute("Bypass", isSlotBypassed(i)); // add later
    child->addChildElement(modules[i]->getStateAsXml(typeString, markAsClean));
    xml->addChildElement(child);
  }

  xml->addChildElement(modManager->getStateAsXml()); 
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

  /*
  // just for test:
  saveXmlToFile(xmlState, File("E:/TempData/xmlOld.xml"), false);

  String xmlStr = xmlState.createDocument("");
  File xmlFile("E:/TempData/xmlStr.xml");
  if( xmlFile.existsAsFile() )
    xmlFile.deleteFile();
  xmlFile.create();
  xmlFile.appendText(xmlStr);

  //// doesn't work:
  //XmlElement newXml = XmlElement(xmlStr);
  //saveXmlToFile(newXml, File("E:/TempData/xmlNew.xml"), false);

  XmlElement* newXml = stringToXml(xmlStr);
  saveXmlToFile(*newXml, File("E:/TempData/xmlNew.xml"), false);
  delete newXml;
  */

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
    //// of recalling the state after adding) because when loading a preset, the editor may show
    //// a module in the initial state instead of the actual state - check with EchoLab - load
    //// ToolChain presets - the EchoLab editor state is wrong after recall

    // new:
    //AudioModule *m = AudioModuleFactory::createModule(type, lock, &modManager, metaParamManager);
    AudioModule *m = moduleFactory.createModule(type);
    if(m != nullptr) {
      m->setSampleRate(sampleRate);
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
  modManager->removeAllConnections();
  XmlElement* modXml = xmlState.getChildByName("Modulations");
  if(modXml != nullptr)
    modManager->setStateFromXml(*modXml);  // recall modulation settings
}

void ToolChain::setupManagers(AudioModule* m)
{
  m->setSmoothingManager(smoothingManager);
  m->setMetaParameterManager(metaParamManager); 
  ModulatableAudioModule* mm = dynamic_cast<ModulatableAudioModule*>(m);
  if(mm)
    mm->setModulationManager(modManager);
  AudioModulePoly* pm = dynamic_cast<AudioModulePoly*> (m);
  if(pm != nullptr)
  {
    pm->setVoiceManager(&voiceManager);
    pm->setVoiceSignalBuffer(&voiceSignals[0]);
  }
}

/*
void ToolChain::allocateVoiceResources(rosic::rsVoiceManager* voiceManager)
{
  voiceSignals.resize(2*voiceManager->getMaxNumVoices());
  for(size_t i = 0; i < modules.size(); i++) 
  {
    AudioModulePoly* pm = dynamic_cast<AudioModulePoly*> (modules[i]);
    if(pm != nullptr)
    {
      pm->allocateVoiceResources(&voiceManager);
      pm->setVoiceSignalBuffer(&voiceSignals[0]);
    }
  }
}
*/

void ToolChain::addToModulatorsIfApplicable(AudioModule* module)
{
  ModulationSource* ms = dynamic_cast<ModulationSource*> (module);
  if(ms != nullptr)
  {
    assignModulationSourceName(ms);
    modManager->registerModulationSource(ms);
  }
}

void ToolChain::removeFromModulatorsIfApplicable(AudioModule* module)
{
  ModulationSource* ms = dynamic_cast<ModulationSource*> (module);
  if(ms != nullptr)
    modManager->deRegisterModulationSource(ms);
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
    name = "Slot" + String(slotIndex+1) + String("-") + am->getModuleTypeName();
    source->setModulationSourceName(name);
  }
}

void ToolChain::clearModulesArray()
{
  ScopedLock scopedLock(*lock);
  while(size(modules) > 0)
    deleteLastModule();
}

void ToolChain::createMidiModSources()
{
  constantModulator = 
    new rsConstantOneModulatorModulePoly(lock, metaParamManager, modManager);
  notePitchModulator = 
    new rsNotePitchModulatorModulePoly(lock, metaParamManager, modManager);
  noteFreqModulator = 
    new rsNoteFreqModulatorModulePoly(lock, metaParamManager, modManager);
  noteVelocityModulator = 
    new rsNoteVelocityModulatorModulePoly(lock, metaParamManager, modManager);

  constantModulator->setVoiceManager(&voiceManager);
  notePitchModulator->setVoiceManager(&voiceManager);
  noteFreqModulator->setVoiceManager(&voiceManager);
  noteVelocityModulator->setVoiceManager(&voiceManager);

  modManager->registerModulationSource(constantModulator);
  modManager->registerModulationSource(notePitchModulator);
  modManager->registerModulationSource(noteFreqModulator);
  modManager->registerModulationSource(noteVelocityModulator);
}

void ToolChain::deleteMidiModSources()
{
  delete constantModulator;
  delete notePitchModulator;
  delete noteFreqModulator;
  delete noteVelocityModulator;
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

  modManager->registerModulationSource(env);

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

void ToolChain::populateModuleFactory()
{
  //ScopedLock scopedLock(*lock);

  typedef CriticalSection* CS;
  typedef AudioModule* AM;
  CS cs = lock;
  AudioModuleFactory& f = moduleFactory;

  // todo: maybe allow subcategories, maybe pass the category string as last as last argument to 
  // registerModuleType, default to the empty string (meaning it should appear at root level)
  // and add a subcategory string parameter, following the same convention)

  juce::String s = "";
  f.registerModuleType([](CS cs)->AM { return new DummyModule(cs); },      s, "None");
#ifdef JUCE_DEBUG
  f.registerModuleType([](CS cs)->AM { return new DebugAudioModule(cs); }, s, "DebugAudioModule");
#endif

  s = "Sources";
  // no sources available yet, todo: move TriSawOsc, etc. here, when done, implement QuadSource

  s = "Filters";
  f.registerModuleType([](CS cs)->AM { return new EqualizerAudioModule(cs);       }, s, "Equalizer");
  f.registerModuleType([](CS cs)->AM { return new Ladder(cs);                     }, s, "Ladder");
  f.registerModuleType([](CS cs)->AM { return new EngineersFilterAudioModule(cs); }, s, "EngineersFilter");

  s = "Modulators";
  f.registerModuleType([](CS cs)->AM { return new BreakpointModulatorAudioModule(cs); }, s, "BreakpointModulator");


  s = "Dynamics";
  f.registerModuleType([](CS cs)->AM { return new LimiterAudioModule(cs);   }, s, "Limiter");

  s = "Effects";
  f.registerModuleType([](CS cs)->AM { return new FuncShaperAudioModule(cs);   }, s, "FuncShaper");
  f.registerModuleType([](CS cs)->AM { return new EchoLabAudioModule(cs);      }, s, "EchoLab");
  f.registerModuleType([](CS cs)->AM { return new StereoDelayAudioModule(cs);  }, s, "StereoDelay");

  s = "Analysis";
  f.registerModuleType([](CS cs)->AM { return new PhaseScope(cs); },               s, "Scope");
  f.registerModuleType([](CS cs)->AM { return new MultiAnalyzerAudioModule(cs); }, s, "MultiAnalyzer");
  f.registerModuleType([](CS cs)->AM { return new TrackMeterAudioModule(cs); },    s, "TrackMeter");
  f.registerModuleType([](CS cs)->AM { return new MidiMonitorAudioModule(cs); },   s, "MidiMonitor");

  s = "Instruments";
  f.registerModuleType([](CS cs)->AM { return new AciDevilAudioModule(cs);      }, s, "AcidDevil");
  f.registerModuleType([](CS cs)->AM { return new StraightlinerAudioModule(cs); }, s, "Straightliner");


  // the large collection of modules that are still unfinished (wrap this into an ifdef (or simple if)):
  s = "Under Construction";

  // Sources:
  f.registerModuleType([](CS cs)->AM { return new SineOscAudioModule(cs);  },           s, "SineOscillator");
  f.registerModuleType([](CS cs)->AM { return new SineOscAudioModulePoly(cs);  },       s, "SineOscillatorPoly");
  f.registerModuleType([](CS cs)->AM { return new TriSawOscModule(cs);  },              s, "TriSawOscillator");
  f.registerModuleType([](CS cs)->AM { return new EllipseOscillatorAudioModule(cs);  }, s, "EllipseOscillator");
  f.registerModuleType([](CS cs)->AM { return new RotationOscillatorAudioModule(cs); }, s, "Oscillator3D");
  f.registerModuleType([](CS cs)->AM { return new RayBouncerAudioModule(cs);         }, s, "RayBouncer");
  f.registerModuleType([](CS cs)->AM { return new Snowflake(cs);         },             s, "Snowflake");
  f.registerModuleType([](CS cs)->AM { return new WaveOscModule(cs);   }, s, "WaveOscillator");
  // DualWaveOsc, WaveScanningOsc
  //f.registerModuleType([](CS cs)->AM { return new FourOscSectionAudioModule(cs);     }, s, "FourOscSection");

  // Filters:
  //f.registerModuleType([](CS cs)->AM { return new PhasorFilter(cs);               }, s, "PhasorFilter");
  //f.registerModuleType([](CS cs)->AM { return new CrossOverAudioModule(cs);       }, s, "CrossOver");

  // Modulators:
  f.registerModuleType([](CS cs)->AM { return new AttackDecayEnvelopeModulePoly(cs); },      s, "EnvelopeAD");
  f.registerModuleType([](CS cs)->AM { return new AttackDecayEnvelopeModule(cs); },      s, "AttackDecayEnvelope");
  //f.registerModuleType([](CS cs)->AM { return new AttackDecayEnvelopeModulePoly(cs); },      s, "AttackDecayEnvelope");
  f.registerModuleType([](CS cs)->AM { return new TriSawModulatorModule(cs); },          s, "TriSawModulator");
  // rename to TriSawLFO



  // Effects:
  //f.registerModuleType([](CS cs)->AM { return new NodeShaperAudioModule(cs);   }, s, "NodeShaper");
  //f.registerModuleType([](CS cs)->AM { return new AlgoVerbAudioModule(cs);     }, s, "AlgoVerb");
  //f.registerModuleType([](CS cs)->AM { return new PingPongEchoAudioModule(cs); }, s, "PingPongEcho");
  //f.registerModuleType([](CS cs)->AM { return new PitchShifterAudioModule(cs); }, s, "PitchShifter");
  //f.registerModuleType([](CS cs)->AM { return new DspWorkbenchAudioModule(cs); }, s, "DspWorkbench");
  f.registerModuleType([](CS cs)->AM { return new MultiBandEffect(cs); }, s, "MultiBandEffect");

  // Instruments:
  f.registerModuleType([](CS cs)->AM { return new ModalSynthAudioModule(cs);   }, s, "ModalSynth");
  f.registerModuleType([](CS cs)->AM { return new QuadrifexAudioModule(cs);    }, s, "Quadrifex");
  //f.registerModuleType([](CS cs)->AM { return new NewSynthAudioModule(cs);     }, s, "NewSynth");
  //f.registerModuleType([](CS cs)->AM { return new MagicCarpetAudioModule(cs);   }, s, "MagicCarpet");
  f.registerModuleType([](CS cs)->AM { return new SamplePlayerAudioModule(cs);    }, s, "SamplePlayer");
  f.registerModuleType([](CS cs)->AM { return new BlepOscArrayModule(cs);    }, s, "BlepOscArray");
  //f.registerModuleType([](CS cs)->AM { return new SimpleSamplerAudioModule(cs); }, s, "SimpleSampler");
  //f.registerModuleType([](CS cs)->AM { return new KeyShotAudioModule(cs);       }, s, "KeyShot");
  //f.registerModuleType([](CS cs)->AM { return new QuadrigaAudioModule(cs);      }, s, "Quadriga");
  //f.registerModuleType([](CS cs)->AM { return new WorkhorseAudioModule(cs);     }, s, "Workhorse");
#ifdef _MSC_VER
  f.registerModuleType([](CS cs)->AM { return new LibertyAudioModule(cs);       }, s, "Liberty");
#endif
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
    s = new AudioModuleSelector(&chain->moduleFactory);
    s->setInterceptsMouseClicks(false, false); // we handle them 1st and possibly pass them through
    s->selectItemFromText(chain->modules[numSelectors]->getModuleTypeName(), false);
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
    selectors[i]->selectItemFromText(chain->modules[i]->getModuleTypeName(), false);
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
  // todo: maybe put a number in front of the selector - it's easier to see what module is in which
  // slot when the slots are numbered - especially if there are several slots with the same type
  // of module - in the modualtion setup, we may see Slot7-TriSawModulator, Slot8-TriSawModulator, 
  // etc. and to figure out which is which, we may have to actually count slots

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

  g.setColour(Colour::fromFloatRGBA(0.8125f, 0.8125f, 0.8125f, 1.f));
  // maybe switch depending on widget color-scheme (dark-on-bright vs bright-on-dark)

  jassert(chain->activeSlot < selectors.size()); 
  // ..i once had a weird crash - not sure, if that was out of range

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

/*

Bugs:
-Load patch Acid2 and then _TestSineOscPoly3 -> crash
 -when the patch is loaded, the module will call moduleWillBeDeleted on the editor
 -the editor will call scheduleSelectorArrayUpdate
 -paintOverChildren is called - the selectorArray does not seem to be updated yet - but i guess,
  it should be - seems like updateSelectorArray is called too late (asynchronously)
  -maybe in audioModuleWillBeDeleted, we should call a function deleteSelector(index) similar to
   deleteEditor(index) instead of scheduleSelectorArrayUpdate
-If a limiter (or gain) is placed after a poly-sine osc and the modulator is placed after the 
 limiter, the limiter has no effect
-load a fresh SineOscillatorStereo and play a not: silence. we only get sound if we wire (for 
 example) NoteFrequency to Frequency
-load patch _TestSineOscPoly2 and connect the EnvelopeAD to the frequency - the envelope does
 have effect, but attack/decay parameters in the dsp object do not reflect the settings. moving 
 the sliders for attack/decay also has no effect. in _TestSineOscPoly3 the do have effect. the
 difference is that they are wired to the constant1 modulator (with 0 depth)
-Poly modulators (and generators?) do not update their parameters when the slider is moved but no 
 modulator is connected
-the order of the modulation connections is not always the same as we wired them up
-when moving slider during playing, it sometimes seems to hang and glitch and does not recover
 ...waiting and then playing more notes just produces audio fragments


ToDo:
-FuncShaper: have a DryDelay parameter - mostly for compensating the latency of the AA filter but 
 may be an interestinf effect as well

*/