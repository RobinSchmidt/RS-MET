/*
void AudioModuleCreator::populateFactoryForToolChain(AudioModuleFactory* f, CriticalSection* cs)
{
  typedef CriticalSection* CS;
  typedef AudioModule* AM;

  juce::String s = "";
  f->registerModuleType([](CS cs)->AM { return new DummyModule(cs); },      s, "DummyModule");
  f->registerModuleType([](CS cs)->AM { return new DebugAudioModule(cs); }, s, "DebugAudioModule");

  s = "Analysis";
  f->registerModuleType([](CS cs)->AM { return new PhaseScope(cs); },               s, "Scope");
  f->registerModuleType([](CS cs)->AM { return new MultiAnalyzerAudioModule(cs); }, s, "MultiAnalyzer");
  f->registerModuleType([](CS cs)->AM { return new TrackMeterAudioModule(cs); },    s, "TrackMeter");
  f->registerModuleType([](CS cs)->AM { return new MidiMonitorAudioModule(cs); },   s, "MidiMonitor");

  // maybe this can be made a member of ToolChain ...populateModuleFactory


}
*/