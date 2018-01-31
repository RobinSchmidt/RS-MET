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
  if(type == "EllipseOscillator") return new EllipseOscillatorAudioModule(lock);
  if(type == "Oscillator3D")      return new RotationOscillatorAudioModule(lock, metaMan, modMan);
  if(type == "RayBouncer")        return new RayBouncerAudioModule(lock);
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

void AudioModuleFactory::registerModuleType(std::function<AudioModule*(CriticalSection*)> creatorFunction, 
  const juce::String& category, const juce::String& typeName)
{
#ifdef JUCE_DEBUG
  AudioModule* m = creatorFunction(lock);
  jassert(typeName == m->getModuleTypeName()); // the name under which the module is registered
  delete m;                                    // must match the module's type name (for recall)
#endif

  if(typeName == String::empty)
  {
    AudioModule* m = creatorFunction(lock);
    moduleInfos.push_back(AudioModuleInfo(m->getModuleTypeName(), creatorFunction, category));
    delete m; 
  }
  else
    moduleInfos.push_back(AudioModuleInfo(typeName, creatorFunction, category));
}

AudioModule* AudioModuleFactory::createModule(const juce::String& type)
{
  for(size_t i = 0; i < moduleInfos.size(); i++)
  {
    if(moduleInfos[i].type == type)
      return moduleInfos[i].createInstance(lock);
  }
  jassertfalse;    // unknown module type requested
  return nullptr;
}