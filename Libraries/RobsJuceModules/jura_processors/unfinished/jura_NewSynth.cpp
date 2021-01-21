// construction/destruction:

NewSynthAudioModule::NewSynthAudioModule(CriticalSection *lockToUse) 
  : AudioModuleWithMidiIn(lockToUse), sourceFactory(lockToUse), filterFactory(lockToUse)
  , modulatorFactory(lockToUse), modManager(lockToUse)
{
  setModuleTypeName("NewSynth");
  addChildAudioModule(sourceModule = new QuadSourceAudioModule(lock/*, &synthCore.source*/));
  addChildAudioModule(filterModule = new DualFilterAudioModule(lock/*, &synthCore.filter*/));
  addChildAudioModule(modulatorsModule = new PolyModulatorsAudioModule(lock));

  sourceModule->setModuleFactory(&sourceFactory);
  filterModule->setModuleFactory(&filterFactory);
  modulatorsModule->setModuleFactory(&modulatorFactory);

  populateModuleFactories();
  //createParameters();
}

AudioModuleEditor* NewSynthAudioModule::createEditor(int type)
{
  return new jura::NewSynthEditor(lock, this); // get rid of passing the lock
}

void NewSynthAudioModule::populateModuleFactories()
{
  typedef CriticalSection* CS;
  typedef AudioModule* AM;
  CS cs = lock;
  juce::String s;

  // Later, we should use the polyphonic versions of the modules

  AudioModuleFactory& sf = sourceFactory;
  s = "";
  sf.registerModuleType([](CS cs)->AM { return new DummyModule(cs); },      s, "None");
  sf.registerModuleType([](CS cs)->AM { return new EllipseOscillatorAudioModule(cs);  }, s, "EllipseOscillator");

  AudioModuleFactory& ff = filterFactory;
  s = "";
  ff.registerModuleType([](CS cs)->AM { return new DummyModule(cs); }, s, "None");
  ff.registerModuleType([](CS cs)->AM { return new Ladder(cs);      }, s, "Ladder");

  AudioModuleFactory& mf = modulatorFactory;
  s = "";
  mf.registerModuleType([](CS cs)->AM { return new DummyModule(cs); }, s, "None");
  mf.registerModuleType([](CS cs)->AM { return new BreakpointModulatorAudioModule(cs); }, s, "BreakpointModulator");
}

//=================================================================================================

NewSynthEditor::NewSynthEditor(CriticalSection* lockToUse, NewSynthAudioModule* synthToEdit)
  : AudioModuleEditor(lockToUse), synthModule(synthToEdit)
{
  ScopedLock scopedLock(*lock);
  addChildEditor(sourceEditor = new QuadSourceEditor(lock, synthModule->sourceModule));
  addChildEditor(filterEditor = new DualFilterEditor(lock, synthModule->filterModule));
  addChildEditor(modulatorsEditor = new PolyModulatorsEditor(lock, synthModule->modulatorsModule));
  setSize(720, 720+26); // 26 = bottom of preset section + 4
}

void NewSynthEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();

  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;

  h = h/3;
  int dy = h-2;
  sourceEditor->setBounds(x, y, w, h); y += dy;
  filterEditor->setBounds(x, y, w, h); y += dy;
  h = getHeight() - y;
  modulatorsEditor->setBounds(x, y, w, h);
}


/*
ToDo:
implement modules:

Sources:
-TriSawOsc

Filters:
-LadderFilter

Modulators:
-AttackDecayEnv


*/