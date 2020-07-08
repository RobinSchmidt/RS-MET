
//-------------------------------------------------------------------------------------------------
// construction/destruction:

//QuadrifexAudioModule::QuadrifexAudioModule(CriticalSection *newPlugInLock,
//  rosic::Quadrifex *quadrifexToWrap)
//: AudioModule(newPlugInLock)
QuadrifexAudioModule::QuadrifexAudioModule(CriticalSection* lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse, 
  rosic::Quadrifex *quadrifexToWrap)
  : ModulatableAudioModule(lockToUse, metaManagerToUse, modManagerToUse) 
{
  ScopedLock scopedLock(*lock);

  //jassert(quadrifexToWrap != NULL); // you must pass a valid rosic-object to the constructor
  if( quadrifexToWrap != nullptr )
    wrappedQuadrifex = quadrifexToWrap;
  else
  {
    wrappedQuadrifex = new rosic::Quadrifex;
    wrappedQuadrifexIsOwned = true;
  }

  editor     = NULL; // !!! we keep a pointer to the editor???!!! get rid of that!!!
  setModuleTypeName("Quadrifex");

  matrixModule = new RoutingMatrixAudioModule(lock, &wrappedQuadrifex->mixMatrix);
  matrixModule->setModuleName(juce::String(("RoutingMatrix5x5")));
  addChildAudioModule(matrixModule);

  for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
  {
    // create a bypass-module for each slot:
    rosic::BypassModule* bypassCoreModule =
      static_cast<rosic::BypassModule*> (wrappedQuadrifex->getEffectModule(i)); // this dynamic_cast causes bugs in the release version
    jura::BypassAudioModule *audioModule = new jura::BypassAudioModule(lock,
      bypassCoreModule);
    effectModules[i] = audioModule;
    addChildAudioModule(effectModules[i]);

    // allocate memory to store the states internally:
    bitCrusherStates[i]              = new XmlElement(juce::String(("BitCrusher")));
    chorusStates[i]                  = new XmlElement(juce::String(("Chorus")));
    combBankStates[i]                = new XmlElement(juce::String(("CombBank")));
    combResonatorStates[i]           = new XmlElement(juce::String(("CombResonator")));
    combStereoizerStates[i]          = new XmlElement(juce::String(("CombStereoizer")));
    compressorStates[i]              = new XmlElement(juce::String(("Compressor")));
    dualTwoPoleFilterStates[i]       = new XmlElement(juce::String(("DualTwoPoleFilter")));
    equalizerStates[i]               = new XmlElement(juce::String(("Equalizer")));
    expanderStates[i]                = new XmlElement(juce::String(("Expander")));
    flangerStates[i]                 = new XmlElement(juce::String(("Flanger")));
    formantShifterStates[i]          = new XmlElement(juce::String(("FormantShifter")));
    fourPoleFilterStates[i]          = new XmlElement(juce::String(("FourPoleFilter")));
    frequencyShifterStates[i]        = new XmlElement(juce::String(("FrequencyShifter")));
    harmonicsStates[i]               = new XmlElement(juce::String(("Harmonics")));
    ladderFilterStates[i]            = new XmlElement(juce::String(("LadderFilter")));
    limiterStates[i]                 = new XmlElement(juce::String(("Limiter")));
    modulatedAllpassStates[i]        = new XmlElement(juce::String(("ModulatedAllpass")));
    noiseGateStates[i]               = new XmlElement(juce::String(("NoiseGate")));
    noisifierStates[i]               = new XmlElement(juce::String(("Noisifier")));
    phaserStates[i]                  = new XmlElement(juce::String(("Phaser")));
    phaseStereoizerStates[i]         = new XmlElement(juce::String(("PhaseStereoizer")));
    pingPongEchoStates[i]            = new XmlElement(juce::String(("PingPongEcho")));
    pitchShifterStates[i]            = new XmlElement(juce::String(("PitchShifter")));
    reverbStates[i]                  = new XmlElement(juce::String(("Reverb")));
    ringModulatorStates[i]           = new XmlElement(juce::String(("RingModulator")));
    simpleDelayStates[i]             = new XmlElement(juce::String(("SimpleDelay")));
    sineOscillatorStates[i]          = new XmlElement(juce::String(("SineOscillator")));
    singleSidebandModulatorStates[i] = new XmlElement(juce::String(("SingleSidebandModulator")));
    slewRateLimiterStates[i]         = new XmlElement(juce::String(("SlewRateLimiter")));
    slopeFilterStates[i]             = new XmlElement(juce::String(("SlopeFilter")));
    stereoPanStates[i]               = new XmlElement(juce::String(("StereoPan")));
    stereoWidthStates[i]             = new XmlElement(juce::String(("StereoWidth")));
    tremoloStates[i]                 = new XmlElement(juce::String(("Tremolo")));
    twoPoleFilterStates[i]           = new XmlElement(juce::String(("TwoPoleFilter")));
    vibratoStates[i]                 = new XmlElement(juce::String(("Vibrato")));
    wahWahStates[i]                  = new XmlElement(juce::String(("WahWah")));
    waveShaperStates[i]              = new XmlElement(juce::String(("WaveShaper")));
  }

  initializeAutomatableParameters();
}

QuadrifexAudioModule::~QuadrifexAudioModule()
{
  ScopedLock scopedLock(*lock);

  if(wrappedQuadrifexIsOwned)
    delete wrappedQuadrifex;

  // todo: avoid these deletions, using juce::ScopedPointer or std::unique_ptr - we need also get
  // rid of the deletes in setEffectAlgorithm
  for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
  {
    delete bitCrusherStates[i];
    delete chorusStates[i];
    delete combBankStates[i];
    delete combResonatorStates[i];
    delete combStereoizerStates[i];
    delete compressorStates[i];
    delete dualTwoPoleFilterStates[i];
    delete equalizerStates[i];
    delete expanderStates[i];
    delete flangerStates[i];
    delete formantShifterStates[i];
    delete fourPoleFilterStates[i];
    delete frequencyShifterStates[i];
    delete harmonicsStates[i];
    delete ladderFilterStates[i];
    delete limiterStates[i];
    delete modulatedAllpassStates[i];
    delete noiseGateStates[i];
    delete noisifierStates[i];
    delete phaserStates[i];
    delete phaseStereoizerStates[i];
    delete pingPongEchoStates[i];
    delete pitchShifterStates[i];
    delete reverbStates[i];
    delete ringModulatorStates[i];
    delete simpleDelayStates[i];
    delete sineOscillatorStates[i];
    delete singleSidebandModulatorStates[i];
    delete slewRateLimiterStates[i];
    delete slopeFilterStates[i];
    delete stereoPanStates[i];
    delete stereoWidthStates[i];
    delete tremoloStates[i];
    delete twoPoleFilterStates[i];
    delete vibratoStates[i];
    delete wahWahStates[i];
    delete waveShaperStates[i];
  }
}

AudioModuleEditor* QuadrifexAudioModule::createEditor(int type)
{
  return new QuadrifexModuleEditor(lock, this);
}

//-------------------------------------------------------------------------------------------------
// setup:

void QuadrifexAudioModule::setEditor(QuadrifexModuleEditor *newEditor)
{
  ScopedLock scopedLock(*lock);
  editor = newEditor;
}

void QuadrifexAudioModule::setEffectAlgorithm(int slotIndex, int newAlgorithmIndex)
{
  ScopedLock scopedLock(*lock);

  if( wrappedQuadrifex == NULL )
    return;

  if( slotIndex < 0 || slotIndex >= rosic::Quadrifex::numEffectSlots )
    return;

  // store the state of the old effect to be replaced:
  // this code is horrible - can we replace this maybe with something based on an associative
  // array, like std::map?

  int oldAlgorithmIndex = wrappedQuadrifex->getEffectAlgorithmIndex(slotIndex);
  switch( oldAlgorithmIndex )
  {
  case rosic::Quadrifex::BIT_CRUSHER:
    {
      jura::BitCrusherAudioModule *audioModule =
        static_cast<jura::BitCrusherAudioModule*> (effectModules[slotIndex]);
      delete bitCrusherStates[slotIndex];
      bitCrusherStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::CHORUS:
    {
      jura::ChorusAudioModule *audioModule =
        static_cast<jura::ChorusAudioModule*> (effectModules[slotIndex]);
      delete chorusStates[slotIndex];
      chorusStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::COMB_BANK:
    {
      jura::CombBankAudioModule *audioModule =
        static_cast<jura::CombBankAudioModule*> (effectModules[slotIndex]);
      delete combBankStates[slotIndex];
      combBankStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::COMB_RESONATOR:
    {
      jura::CombResonatorAudioModule *audioModule =
        static_cast<jura::CombResonatorAudioModule*> (effectModules[slotIndex]);
      delete combResonatorStates[slotIndex];
      combResonatorStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::COMB_STEREOIZER:
    {
      jura::CombStereoizerAudioModule *audioModule =
        static_cast<jura::CombStereoizerAudioModule*> (effectModules[slotIndex]);
      delete combStereoizerStates[slotIndex];
      combStereoizerStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::COMPRESSOR:
    {
      jura::CompressorAudioModule *audioModule =
        static_cast<jura::CompressorAudioModule*> (effectModules[slotIndex]);
      delete compressorStates[slotIndex];
      compressorStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::DUAL_TWO_POLE_FILTER:
    {
      jura::DualTwoPoleFilterAudioModule *audioModule =
        static_cast<jura::DualTwoPoleFilterAudioModule*> (effectModules[slotIndex]);
      delete dualTwoPoleFilterStates[slotIndex];
      dualTwoPoleFilterStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::EQUALIZER:
    {
      jura::EqualizerAudioModule *audioModule =
        static_cast<jura::EqualizerAudioModule*> (effectModules[slotIndex]);
      delete equalizerStates[slotIndex];
      equalizerStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::EXPANDER:
    {
      jura::ExpanderAudioModule *audioModule =
        static_cast<jura::ExpanderAudioModule*> (effectModules[slotIndex]);
      delete expanderStates[slotIndex];
      expanderStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::FLANGER:
    {
      jura::FlangerAudioModule *audioModule =
        static_cast<jura::FlangerAudioModule*> (effectModules[slotIndex]);
      delete flangerStates[slotIndex];
      flangerStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::FORMANT_SHIFTER:
    {
      jura::FormantShifterAudioModule *audioModule =
        static_cast<jura::FormantShifterAudioModule*> (effectModules[slotIndex]);
      delete formantShifterStates[slotIndex];
      formantShifterStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::FOUR_POLE_FILTER:
    {
      jura::FourPoleFilterAudioModule *audioModule =
        static_cast<jura::FourPoleFilterAudioModule*> (effectModules[slotIndex]);
      delete fourPoleFilterStates[slotIndex];
      fourPoleFilterStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::FREQUENCY_SHIFTER:
    {
      jura::FrequencyShifterAudioModule *audioModule =
        static_cast<jura::FrequencyShifterAudioModule*> (effectModules[slotIndex]);
      delete frequencyShifterStates[slotIndex];
      frequencyShifterStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::HARMONICS:
    {
      jura::HarmonicsAudioModule *audioModule =
        static_cast<jura::HarmonicsAudioModule*> (effectModules[slotIndex]);
      delete harmonicsStates[slotIndex];
      harmonicsStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::LADDER_FILTER:
    {
      jura::LadderFilterAudioModule *audioModule =
        static_cast<jura::LadderFilterAudioModule*> (effectModules[slotIndex]);
      delete ladderFilterStates[slotIndex];
      ladderFilterStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::LIMITER:
    {
      jura::LimiterAudioModule *audioModule =
        static_cast<jura::LimiterAudioModule*> (effectModules[slotIndex]);
      delete limiterStates[slotIndex];
      limiterStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::MODULATED_ALLPASS:
    {
      jura::ModulatedAllpassAudioModule *audioModule =
        static_cast<jura::ModulatedAllpassAudioModule*> (effectModules[slotIndex]);
      delete modulatedAllpassStates[slotIndex];
      modulatedAllpassStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::NOISE_GATE:
    {
      jura::NoiseGateAudioModule *audioModule =
        static_cast<jura::NoiseGateAudioModule*> (effectModules[slotIndex]);
      delete noiseGateStates[slotIndex];
      noiseGateStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::NOISIFIER:
    {
      jura::NoisifierAudioModule *audioModule =
        static_cast<jura::NoisifierAudioModule*> (effectModules[slotIndex]);
      delete noisifierStates[slotIndex];
      noisifierStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::PHASER:
    {
      jura::PhaserAudioModule *audioModule =
        static_cast<jura::PhaserAudioModule*> (effectModules[slotIndex]);
      delete phaserStates[slotIndex];
      phaserStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::PHASE_STEREOIZER:
    {
      jura::PhaseStereoizerAudioModule *audioModule =
        static_cast<jura::PhaseStereoizerAudioModule*> (effectModules[slotIndex]);
      delete phaseStereoizerStates[slotIndex];
      phaseStereoizerStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::PINGPONG_ECHO:
    {
      jura::PingPongEchoAudioModule *audioModule =
        static_cast<jura::PingPongEchoAudioModule*> (effectModules[slotIndex]);
      delete pingPongEchoStates[slotIndex];
      pingPongEchoStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::PITCH_SHIFTER:
    {
      jura::PitchShifterAudioModule *audioModule =
        static_cast<jura::PitchShifterAudioModule*> (effectModules[slotIndex]);
      delete pitchShifterStates[slotIndex];
      pitchShifterStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::REVERB:
    {
      jura::ReverbAudioModule *audioModule =
        static_cast<jura::ReverbAudioModule*> (effectModules[slotIndex]);
      delete reverbStates[slotIndex];
      reverbStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::RINGMODULATOR:
    {
      jura::RingModulatorAudioModule *audioModule =
        static_cast<jura::RingModulatorAudioModule*> (effectModules[slotIndex]);
      delete ringModulatorStates[slotIndex];
      ringModulatorStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::SIMPLE_DELAY:
    {
      jura::SimpleDelayAudioModule *audioModule =
        static_cast<jura::SimpleDelayAudioModule*> (effectModules[slotIndex]);
      delete simpleDelayStates[slotIndex];
      simpleDelayStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::SINE_OSCILLATOR:
    {
      jura::SineOscillatorAudioModule *audioModule =
        static_cast<jura::SineOscillatorAudioModule*> (effectModules[slotIndex]);
      delete sineOscillatorStates[slotIndex];
      sineOscillatorStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::SLOPE_FILTER:
    {
      jura::SlopeFilterAudioModule *audioModule =
        static_cast<jura::SlopeFilterAudioModule*> (effectModules[slotIndex]);
      delete slopeFilterStates[slotIndex];
      slopeFilterStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::SSB_MODULATOR:
    {
      jura::SingleSidebandModulatorAudioModule *audioModule =
        static_cast<jura::SingleSidebandModulatorAudioModule*> (effectModules[slotIndex]);
      delete singleSidebandModulatorStates[slotIndex];
      singleSidebandModulatorStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::SLEWRATE_LIMITER:
    {
      jura::SlewRateLimiterAudioModule *audioModule =
        static_cast<jura::SlewRateLimiterAudioModule*> (effectModules[slotIndex]);
      delete slewRateLimiterStates[slotIndex];
      slewRateLimiterStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::STEREO_PAN:
    {
      jura::StereoPanAudioModule *audioModule =
        static_cast<jura::StereoPanAudioModule*> (effectModules[slotIndex]);
      delete stereoPanStates[slotIndex];
      stereoPanStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::STEREO_WIDTH:
    {
      jura::StereoWidthAudioModule *audioModule =
        static_cast<jura::StereoWidthAudioModule*> (effectModules[slotIndex]);
      delete stereoWidthStates[slotIndex];
      stereoWidthStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::TREMOLO:
    {
      jura::TremoloAudioModule *audioModule =
        static_cast<jura::TremoloAudioModule*> (effectModules[slotIndex]);
      delete tremoloStates[slotIndex];
      tremoloStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::TWO_POLE_FILTER:
    {
      jura::TwoPoleFilterAudioModule *audioModule =
        static_cast<jura::TwoPoleFilterAudioModule*> (effectModules[slotIndex]);
      delete twoPoleFilterStates[slotIndex];
      twoPoleFilterStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::VIBRATO:
    {
      jura::VibratoAudioModule *audioModule =
        static_cast<jura::VibratoAudioModule*> (effectModules[slotIndex]);
      delete vibratoStates[slotIndex];
      vibratoStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::WAH_WAH:
    {
      jura::WahWahAudioModule *audioModule =
        static_cast<jura::WahWahAudioModule*> (effectModules[slotIndex]);
      delete wahWahStates[slotIndex];
      wahWahStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;
  case rosic::Quadrifex::WAVESHAPER:
    {
      jura::WaveShaperAudioModule *audioModule =
        static_cast<jura::WaveShaperAudioModule*> (effectModules[slotIndex]);
      delete waveShaperStates[slotIndex];
      waveShaperStates[slotIndex] = audioModule->getStateAsXml(juce::String(), false);
    } break;


    //....................tbc...............

  }

  // delete the old effect (and its (sub)editor, if present):
  if( editor != NULL )
    editor->removeChildEditorInSlot(slotIndex);
  effectModules[slotIndex]->removeAllStateWatchers();


  // hangs on plugging in a Reverb...
  removeChildAudioModule(effectModules[slotIndex], true);
  //


  wrappedQuadrifex->setEffectAlgorithm(slotIndex, newAlgorithmIndex);

  // create the new effect and restore its state from any previous use:
  switch( newAlgorithmIndex )
  {
  case rosic::Quadrifex::BIT_CRUSHER:
    {
      rosic::BitCrusherModule *core =
        static_cast<rosic::BitCrusherModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::BitCrusherAudioModule *audioModule = new jura::BitCrusherAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("BitCrusher")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*bitCrusherStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::CHORUS:
    {
      rosic::ChorusModule *core =
        static_cast<rosic::ChorusModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::ChorusAudioModule *audioModule = new jura::ChorusAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Chorus")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*chorusStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::COMB_BANK:
    {
      rosic::CombBankModule *core =
        static_cast<rosic::CombBankModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::CombBankAudioModule *audioModule = new jura::CombBankAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("CombBank")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*combBankStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::COMB_RESONATOR:
    {
      rosic::CombResonatorStereoModule *core =
        static_cast<rosic::CombResonatorStereoModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::CombResonatorAudioModule *audioModule = new jura::CombResonatorAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("CombResonator")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*combResonatorStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::COMB_STEREOIZER:
    {
      rosic::CombStereoizerModule *core =
        static_cast<rosic::CombStereoizerModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::CombStereoizerAudioModule *audioModule = new jura::CombStereoizerAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("CombStereoizer")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*combStereoizerStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::COMPRESSOR:
    {
      rosic::SoftKneeCompressorModule *core =
        static_cast<rosic::SoftKneeCompressorModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::CompressorAudioModule *audioModule = new jura::CompressorAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Compressor")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*compressorStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::DUAL_TWO_POLE_FILTER:
    {
      rosic::DualTwoPoleFilterModule *core =
        static_cast<rosic::DualTwoPoleFilterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::DualTwoPoleFilterAudioModule *audioModule = new jura::DualTwoPoleFilterAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("DualTwoPoleFilter")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*dualTwoPoleFilterStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::EQUALIZER:
    {
      rosic::EqualizerModule *core =
        static_cast<rosic::EqualizerModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::EqualizerAudioModule *audioModule = new jura::EqualizerAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Equalizer")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*equalizerStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::EXPANDER:
    {
      rosic::SoftKneeExpanderModule *core =
        static_cast<rosic::SoftKneeExpanderModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::ExpanderAudioModule *audioModule = new jura::ExpanderAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Expander")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*expanderStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::FLANGER:
    {
      rosic::FlangerModule *core =
        static_cast<rosic::FlangerModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::FlangerAudioModule *audioModule = new jura::FlangerAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Flanger")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*flangerStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::FORMANT_SHIFTER:
    {
      rosic::FormantShifterModule *core =
        static_cast<rosic::FormantShifterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::FormantShifterAudioModule *audioModule = new jura::FormantShifterAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("FormantShifter")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*formantShifterStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::FOUR_POLE_FILTER:
    {
      rosic::FourPoleFilterModule *core =
        static_cast<rosic::FourPoleFilterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::FourPoleFilterAudioModule *audioModule = new jura::FourPoleFilterAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("FourPoleFilter")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*fourPoleFilterStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::FREQUENCY_SHIFTER:
    {
      rosic::FrequencyShifterStereoModule *core =
        static_cast<rosic::FrequencyShifterStereoModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::FrequencyShifterAudioModule *audioModule = new jura::FrequencyShifterAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("FrequencyShifter")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*frequencyShifterStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::HARMONICS:
    {
      rosic::HarmonicsModule *core =
        static_cast<rosic::HarmonicsModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::HarmonicsAudioModule *audioModule = new jura::HarmonicsAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Harmonics")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*harmonicsStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::LADDER_FILTER:
    {
      rosic::LadderFilterModule *core =
        static_cast<rosic::LadderFilterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::LadderFilterAudioModule *audioModule = new jura::LadderFilterAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("LadderFilter")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*ladderFilterStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::LIMITER:
    {
      rosic::LimiterModule *core =
        static_cast<rosic::LimiterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::LimiterAudioModule *audioModule = new jura::LimiterAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Limiter")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*limiterStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::MODULATED_ALLPASS:
    {
      rosic::ModulatedAllpassModule *core =
        static_cast<rosic::ModulatedAllpassModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::ModulatedAllpassAudioModule *audioModule = new jura::ModulatedAllpassAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("ModulatedAllpass")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*modulatedAllpassStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::NOISE_GATE:
    {
      rosic::NoiseGateModule *core =
        static_cast<rosic::NoiseGateModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::NoiseGateAudioModule *audioModule = new jura::NoiseGateAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("NoiseGate")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*noiseGateStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::NOISIFIER:
    {
      rosic::NoisifierModule *core =
        static_cast<rosic::NoisifierModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::NoisifierAudioModule *audioModule = new jura::NoisifierAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Noisifier")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*noisifierStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::PHASER:
    {
      rosic::PhaserModule *core =
        static_cast<rosic::PhaserModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::PhaserAudioModule *audioModule = new jura::PhaserAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Phaser")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*phaserStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::PHASE_STEREOIZER:
    {
      rosic::PhaseStereoizerModule *core =
        static_cast<rosic::PhaseStereoizerModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::PhaseStereoizerAudioModule *audioModule = new jura::PhaseStereoizerAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("PhaseStereoizer")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*phaseStereoizerStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::PINGPONG_ECHO:
    {
      rosic::PingPongEchoModule *core =
        static_cast<rosic::PingPongEchoModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::PingPongEchoAudioModule *audioModule = new jura::PingPongEchoAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("PingPongEcho")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*pingPongEchoStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::PITCH_SHIFTER:
    {
      rosic::PitchShifterModule *core =
        static_cast<rosic::PitchShifterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::PitchShifterAudioModule *audioModule = new jura::PitchShifterAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("PitchShifter")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*pitchShifterStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::REVERB:
    {
      rosic::rsReverbModule *core =
        static_cast<rosic::rsReverbModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::ReverbAudioModule *audioModule = new jura::ReverbAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Reverb")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*reverbStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::RINGMODULATOR:
    {
      rosic::RingModulatorModule *core =
        static_cast<rosic::RingModulatorModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::RingModulatorAudioModule *audioModule = new jura::RingModulatorAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("RingModulator")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*ringModulatorStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::SIMPLE_DELAY:
    {
      rosic::SimpleDelayModule *core =
        static_cast<rosic::SimpleDelayModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::SimpleDelayAudioModule *audioModule = new jura::SimpleDelayAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("SimpleDelay")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*simpleDelayStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::SINE_OSCILLATOR:
    {
      rosic::SineOscillatorModule *core =
        static_cast<rosic::SineOscillatorModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::SineOscillatorAudioModule *audioModule = new jura::SineOscillatorAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("SineOscillator")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*sineOscillatorStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::SLOPE_FILTER:
    {
      rosic::SlopeFilterModule *core =
        static_cast<rosic::SlopeFilterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::SlopeFilterAudioModule *audioModule = new jura::SlopeFilterAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("SlopeFilter")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*slopeFilterStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::SSB_MODULATOR:
    {
      rosic::SingleSidebandModulatorModule *core =
        static_cast<rosic::SingleSidebandModulatorModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::SingleSidebandModulatorAudioModule *audioModule = new jura::SingleSidebandModulatorAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("SingleSidebandModulator")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*singleSidebandModulatorStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::SLEWRATE_LIMITER:
    {
      rosic::SlewRateLimiterStereoModule *core =
        static_cast<rosic::SlewRateLimiterStereoModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::SlewRateLimiterAudioModule *audioModule = new jura::SlewRateLimiterAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("SlewRateLimiter")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*slewRateLimiterStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::STEREO_PAN:
    {
      rosic::StereoPanModule *core =
        static_cast<rosic::StereoPanModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::StereoPanAudioModule *audioModule = new jura::StereoPanAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("StereoPan")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*stereoPanStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::STEREO_WIDTH:
    {
      rosic::StereoWidthModule *core =
        static_cast<rosic::StereoWidthModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::StereoWidthAudioModule *audioModule = new jura::StereoWidthAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("StereoWidth")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*stereoWidthStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::TREMOLO:
    {
      rosic::TremoloModule *core =
        static_cast<rosic::TremoloModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::TremoloAudioModule *audioModule = new jura::TremoloAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Tremolo")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*tremoloStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::TWO_POLE_FILTER:
    {
      rosic::TwoPoleFilterModule *core =
        static_cast<rosic::TwoPoleFilterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::TwoPoleFilterAudioModule *audioModule = new jura::TwoPoleFilterAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("TwoPoleFilter")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*twoPoleFilterStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::VIBRATO:
    {
      rosic::VibratoModule *core =
        static_cast<rosic::VibratoModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::VibratoAudioModule *audioModule = new jura::VibratoAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("Vibrato")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*vibratoStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::WAH_WAH:
    {
      rosic::WahWahModule *core =
        static_cast<rosic::WahWahModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::WahWahAudioModule *audioModule = new jura::WahWahAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("WahWah")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*wahWahStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;
  case rosic::Quadrifex::WAVESHAPER:
    {
      rosic::WaveShaperModule *core =
        static_cast<rosic::WaveShaperModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::WaveShaperAudioModule *audioModule = new jura::WaveShaperAudioModule(lock, core);
      audioModule->setModuleName(juce::String(("WaveShaper")) + juce::String(slotIndex+1));
      addChildAudioModule(audioModule);
      audioModule->setStateFromXml(*waveShaperStates[slotIndex], juce::String(), true);
      effectModules[slotIndex] = audioModule;
    } break;


    // special trivial cases:
  case rosic::Quadrifex::MUTE:
    {
      rosic::MuteModule *core =
        static_cast<rosic::MuteModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::MuteAudioModule *audioModule = new jura::MuteAudioModule(lock, core);
      addChildAudioModule(audioModule);
      effectModules[slotIndex] = audioModule;
    } break;
  default: // bypass by default (i.e. value out of range)
    {
      rosic::BypassModule *core =
        static_cast<rosic::BypassModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      jura::BypassAudioModule *audioModule = new jura::BypassAudioModule(lock, core);
      addChildAudioModule(audioModule);
      effectModules[slotIndex] = audioModule;
    } break;
  }

  // let the editor (if present) create an appropriate child-editor:
  if( editor != NULL )
    editor->createEditorForSlot(slotIndex, newAlgorithmIndex);
}

//-------------------------------------------------------------------------------------------------
// automation:

void QuadrifexAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedLock scopedLock(*lock);

  if( wrappedQuadrifex == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedQuadrifex->setDryWet(  value); break;
  case   1: wrappedQuadrifex->setWetLevel(value); break;
  case   2: triggerInterval = value;
  } // end of switch( parameterIndex )

}

//-------------------------------------------------------------------------------------------------
// state saving and recall:

XmlElement* QuadrifexAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*lock);

  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  if( wrappedQuadrifex != NULL )
  {
    // store the routing and slot-effect assignments:
    xmlState->setAttribute(("Routing"),
      slotRoutingIndexToString(wrappedQuadrifex->getSlotRouting()));
    for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
    {
      xmlState->setAttribute(juce::String(("Slot"))+juce::String(i+1),
        effectAlgorithmIndexToString(wrappedQuadrifex->getEffectAlgorithmIndex(i)) );
    }
  }
  return xmlState;
}

void QuadrifexAudioModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
                                               bool markAsClean)
{
  ScopedLock scopedLock(*lock);

  if( wrappedQuadrifex != NULL )
  {
    // recall the routing and slot-effect assignments:
    wrappedQuadrifex->setSlotRouting(
      stringToSlotRoutingIndex( xmlState.getStringAttribute(("Routing"), ("1>2>3>4")) ) );
    for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
    {
      setEffectAlgorithm(i, stringToEffectAlgorithmIndex(
        xmlState.getStringAttribute( juce::String(("Slot"))+juce::String(i+1), ("Mute"))));
    }
  }
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void QuadrifexAudioModule::initializeAutomatableParameters()
{
  ScopedLock scopedLock(*lock);


  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  juce::Array<double> defaultValues;

  // this pointer will be used to temporarily store the addresses of the created Parameter-objects:
  AutomatableParameter* p;

  // #00:
  p = new AutomatableParameter(lock, "DryWet", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p);

  // #01:
  p = new AutomatableParameter(lock, "WetLevel", -36.0, 6.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #02:
  p = new AutomatableParameter(lock, "TriggerInterval", 0.0, 64.0, 1.0, 8.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

juce::String QuadrifexAudioModule::slotRoutingIndexToString(int index)
{
  ScopedLock scopedLock(*lock);

  // this is kinda stupid - we should use an associative container like std::map

  switch( index )
  {
  case rosic::Quadrifex::R_BYPASS:            return juce::String(("Bypass"));
  case rosic::Quadrifex::R_1TO2TO3TO4:        return juce::String(("1>2>3>4"));
  case rosic::Quadrifex::R_1TO2TO3_PLUS4:     return juce::String(("(1>2>3)+4"));
  case rosic::Quadrifex::R_1TO2_PLUS_3TO4:    return juce::String(("(1>2)+(3>4)"));
  case rosic::Quadrifex::R_1PLUS2PLUS3PLUS4:  return juce::String(("1+2+3+4"));
  case rosic::Quadrifex::R_1PLUS2PLUS3_TO_4:  return juce::String(("(1+2+3)>4"));
  case rosic::Quadrifex::R_1_TO_2PLUS3_TO_4:  return juce::String(("1>(2+3)>4"));
  case rosic::Quadrifex::R_1PLUS2_TO_3TO4:    return juce::String(("(1+2)>3>4"));
  case rosic::Quadrifex::R_1TO2_TO_3PLUS4:    return juce::String(("1>2>(3+4)"));
  case rosic::Quadrifex::R_1PLUS2_TO_3PLUS4:  return juce::String(("(1+2)>(3+4)"));
  case rosic::Quadrifex::R_1_TO_2PLUS3PLUS4:  return juce::String(("1>(2+3+4)"));
  case rosic::Quadrifex::MATRIX:              return juce::String(("Matrix"));
  default:                                    return juce::String(("Bypass"));
  }
}

int QuadrifexAudioModule::stringToSlotRoutingIndex(const juce::String &routingString)
{
  ScopedLock scopedLock(*lock);

  if( routingString == juce::String(("Bypass"))   )       return rosic::Quadrifex::R_BYPASS;
  if( routingString == juce::String(("1>2>3>4"))  )       return rosic::Quadrifex::R_1TO2TO3TO4;
  if( routingString == juce::String(("(1>2>3)+4"))  )     return rosic::Quadrifex::R_1TO2TO3_PLUS4;
  if( routingString == juce::String(("(1>2)+(3>4)"))  )   return rosic::Quadrifex::R_1TO2_PLUS_3TO4;
  if( routingString == juce::String(("1+2+3+4"))  )       return rosic::Quadrifex::R_1PLUS2PLUS3PLUS4;
  if( routingString == juce::String(("(1+2+3)>4"))  )     return rosic::Quadrifex::R_1PLUS2PLUS3_TO_4;
  if( routingString == juce::String(("1>(2+3)>4"))  )     return rosic::Quadrifex::R_1_TO_2PLUS3_TO_4;
  if( routingString == juce::String(("(1+2)>3>4"))  )     return rosic::Quadrifex::R_1PLUS2_TO_3TO4;
  if( routingString == juce::String(("1>2>(3+4)"))  )     return rosic::Quadrifex::R_1TO2_TO_3PLUS4;
  if( routingString == juce::String(("(1+2)>(3+4)"))  )   return rosic::Quadrifex::R_1PLUS2_TO_3PLUS4;
  if( routingString == juce::String(("1>(2+3+4)"))  )     return rosic::Quadrifex::R_1_TO_2PLUS3PLUS4;
  if( routingString == juce::String(("Matrix"))  )        return rosic::Quadrifex::MATRIX;

  return rosic::Quadrifex::R_BYPASS;
}

juce::String QuadrifexAudioModule::effectAlgorithmIndexToString(int index)
{
  ScopedLock scopedLock(*lock);

  // \todo use std::map

  switch( index )
  {
  case rosic::Quadrifex::MUTE:                 return juce::String(("Mute"));
  case rosic::Quadrifex::BYPASS:               return juce::String(("Bypass"));
  case rosic::Quadrifex::BIT_CRUSHER:          return juce::String(("BitCrusher"));
  case rosic::Quadrifex::HARMONICS:            return juce::String(("Harmonics"));
  case rosic::Quadrifex::CHORUS:               return juce::String(("Chorus"));
  case rosic::Quadrifex::COMB_BANK:            return juce::String(("CombBank"));
  case rosic::Quadrifex::COMB_RESONATOR:       return juce::String(("CombResonator"));
  case rosic::Quadrifex::COMB_STEREOIZER:      return juce::String(("CombStereoizer"));
  case rosic::Quadrifex::COMPRESSOR:           return juce::String(("Compressor"));
  //case rosic::Quadrifex::COMP_SHAPER:          return juce::String(("CompShaper"));
  case rosic::Quadrifex::DUAL_TWO_POLE_FILTER: return juce::String(("DualTwoPoleFilter"));
  case rosic::Quadrifex::EQUALIZER:            return juce::String(("Equalizer"));
  case rosic::Quadrifex::EXPANDER:             return juce::String(("Expander"));
  case rosic::Quadrifex::FLANGER:              return juce::String(("Flanger"));
  case rosic::Quadrifex::FORMANT_SHIFTER:      return juce::String(("FormantShifter"));
  case rosic::Quadrifex::FOUR_POLE_FILTER:     return juce::String(("FourPoleFilter"));
  case rosic::Quadrifex::FREQUENCY_SHIFTER:    return juce::String(("FrequencyShifter"));
  case rosic::Quadrifex::LADDER_FILTER:        return juce::String(("LadderFilter"));
  case rosic::Quadrifex::LIMITER:              return juce::String(("Limiter"));
  case rosic::Quadrifex::MODULATED_ALLPASS:    return juce::String(("ModulatedAllpass"));
  case rosic::Quadrifex::NOISE_GATE:           return juce::String(("NoiseGate"));
  case rosic::Quadrifex::NOISIFIER:            return juce::String(("Noisifier"));
  case rosic::Quadrifex::PHASER:               return juce::String(("Phaser"));
  case rosic::Quadrifex::PHASE_STEREOIZER:     return juce::String(("PhaseStereoizer"));
  case rosic::Quadrifex::PINGPONG_ECHO:        return juce::String(("PingPongEcho"));
  case rosic::Quadrifex::PITCH_SHIFTER:        return juce::String(("PitchShifter"));
  case rosic::Quadrifex::REVERB:               return juce::String(("Reverb"));
  case rosic::Quadrifex::RINGMODULATOR:        return juce::String(("RingModulator"));
  case rosic::Quadrifex::SIMPLE_DELAY:         return juce::String(("SimpleDelay"));
  case rosic::Quadrifex::SINE_OSCILLATOR:      return juce::String(("SineOscillator"));
  case rosic::Quadrifex::SSB_MODULATOR:        return juce::String(("SingleSidebandModulator"));
  case rosic::Quadrifex::SLEWRATE_LIMITER:     return juce::String(("SlewRateLimiter"));
  case rosic::Quadrifex::SLOPE_FILTER:         return juce::String(("SlopeFilter"));
  case rosic::Quadrifex::STEREO_PAN:           return juce::String(("StereoPan"));
  case rosic::Quadrifex::STEREO_WIDTH:         return juce::String(("StereoWidth"));
  case rosic::Quadrifex::TREMOLO:              return juce::String(("Tremolo"));
  case rosic::Quadrifex::TWO_POLE_FILTER:      return juce::String(("TwoPoleFilter"));
  case rosic::Quadrifex::VIBRATO:              return juce::String(("Vibrato"));
  case rosic::Quadrifex::WAH_WAH:              return juce::String(("WahWah"));
  case rosic::Quadrifex::WAVESHAPER:           return juce::String(("WaveShaper"));

  default:                                     return juce::String(("Mute"));
  }
}

int QuadrifexAudioModule::stringToEffectAlgorithmIndex(const juce::String &algoString)
{
  ScopedLock scopedLock(*lock);

  if( algoString == juce::String(("Mute"))   )                   return rosic::Quadrifex::MUTE;
  if( algoString == juce::String(("Bypass")) )                   return rosic::Quadrifex::BYPASS;
  if( algoString == juce::String(("BitCrusher")) )               return rosic::Quadrifex::BIT_CRUSHER;
  if( algoString == juce::String(("Harmonics")) )                return rosic::Quadrifex::HARMONICS;
  if( algoString == juce::String(("Chorus")) )                   return rosic::Quadrifex::CHORUS;
  if( algoString == juce::String(("CombBank")) )                 return rosic::Quadrifex::COMB_BANK;
  if( algoString == juce::String(("CombResonator")) )            return rosic::Quadrifex::COMB_RESONATOR;
  if( algoString == juce::String(("CombStereoizer")) )           return rosic::Quadrifex::COMB_STEREOIZER;
  if( algoString == juce::String(("Compressor")) )               return rosic::Quadrifex::COMPRESSOR;
  //if( algoString == juce::String(("CompShaper")) )               return rosic::Quadrifex::COMP_SHAPER;
  if( algoString == juce::String(("DualTwoPoleFilter")) )        return rosic::Quadrifex::DUAL_TWO_POLE_FILTER;
  if( algoString == juce::String(("Equalizer")) )                return rosic::Quadrifex::EQUALIZER;
  if( algoString == juce::String(("Expander")) )                 return rosic::Quadrifex::EXPANDER;
  if( algoString == juce::String(("Flanger")) )                  return rosic::Quadrifex::FLANGER;
  if( algoString == juce::String(("FormantShifter")) )           return rosic::Quadrifex::FORMANT_SHIFTER;
  if( algoString == juce::String(("FourPoleFilter")) )           return rosic::Quadrifex::FOUR_POLE_FILTER;
  if( algoString == juce::String(("FrequencyShifter")) )         return rosic::Quadrifex::FREQUENCY_SHIFTER;
  if( algoString == juce::String(("LadderFilter")) )             return rosic::Quadrifex::LADDER_FILTER;
  if( algoString == juce::String(("Limiter")) )                  return rosic::Quadrifex::LIMITER;
  if( algoString == juce::String(("ModulatedAllpass")) )         return rosic::Quadrifex::MODULATED_ALLPASS;
  if( algoString == juce::String(("NoiseGate")) )                return rosic::Quadrifex::NOISE_GATE;
  if( algoString == juce::String(("Noisifier")) )                return rosic::Quadrifex::NOISIFIER;
  if( algoString == juce::String(("Phaser")) )                   return rosic::Quadrifex::PHASER;
  if( algoString == juce::String(("PhaseStereoizer")) )          return rosic::Quadrifex::PHASE_STEREOIZER;
  if( algoString == juce::String(("PingPongEcho")) )             return rosic::Quadrifex::PINGPONG_ECHO;
  if( algoString == juce::String(("PitchShifter")) )             return rosic::Quadrifex::PITCH_SHIFTER;
  if( algoString == juce::String(("Reverb")) )                   return rosic::Quadrifex::REVERB;
  if( algoString == juce::String(("RingModulator")) )            return rosic::Quadrifex::RINGMODULATOR;
  if( algoString == juce::String(("SimpleDelay")) )              return rosic::Quadrifex::SIMPLE_DELAY;
  if( algoString == juce::String(("SineOscillator")) )           return rosic::Quadrifex::SINE_OSCILLATOR;
  if( algoString == juce::String(("SingleSidebandModulator")) )  return rosic::Quadrifex::SSB_MODULATOR;
  if( algoString == juce::String(("SlewRateLimiter")) )          return rosic::Quadrifex::SLEWRATE_LIMITER;
  if( algoString == juce::String(("SlopeFilter")) )              return rosic::Quadrifex::SLOPE_FILTER;
  if( algoString == juce::String(("StereoPan")) )                return rosic::Quadrifex::STEREO_PAN;
  if( algoString == juce::String(("StereoWidth")) )              return rosic::Quadrifex::STEREO_WIDTH;
  if( algoString == juce::String(("Tremolo")) )                  return rosic::Quadrifex::TREMOLO;
  if( algoString == juce::String(("TwoPoleFilter")) )            return rosic::Quadrifex::TWO_POLE_FILTER;
  if( algoString == juce::String(("Vibrato")) )                  return rosic::Quadrifex::VIBRATO;
  if( algoString == juce::String(("WahWah")) )                   return rosic::Quadrifex::WAH_WAH;
  if( algoString == juce::String(("WaveShaper")) )               return rosic::Quadrifex::WAVESHAPER;

  return rosic::Quadrifex::MUTE;
}

//=================================================================================================


QuadrifexRoutingDiagram::QuadrifexRoutingDiagram(CriticalSection *newPlugInLock) //: RectangleComponent()
{
  plugInLock = newPlugInLock;
  ScopedLock scopedLock(*plugInLock);

  slotRouting = rosic::Quadrifex::R_1TO2TO3TO4;
  for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++ )
    oldAlgorithmIndices[i] = algorithmIndices[i] = rosic::Quadrifex::BYPASS;

  rsPlot::setAxisPositionX(rsPlotSettings::INVISIBLE);
  rsPlot::setAxisPositionY(rsPlotSettings::INVISIBLE);
  rsPlot::setAxisValuesPositionX(rsPlotSettings::INVISIBLE);
  rsPlot::setAxisValuesPositionY(rsPlotSettings::INVISIBLE);
}

QuadrifexRoutingDiagram::~QuadrifexRoutingDiagram()
{
  ScopedLock scopedLock(*plugInLock);
}

//-------------------------------------------------------------------------------------------------
// setup:

void QuadrifexRoutingDiagram::setAlgorithmIndex(int slotIndex, int newAlgorithmIndex)
{
  ScopedLock scopedLock(*plugInLock);
  if( slotIndex >= 0 && slotIndex < rosic::Quadrifex::numEffectSlots )
  {
    if( newAlgorithmIndex >= rosic::Quadrifex::MUTE &&
      newAlgorithmIndex < rosic::Quadrifex::NUM_EFFECT_ALGORITHMS )
    {
      algorithmIndices[slotIndex] = newAlgorithmIndex;
    }
    else
    {
      algorithmIndices[slotIndex] = rosic::Quadrifex::BYPASS;
      jassertfalse; // newAlgorithmIndex out of range
    }
  }
  else
    jassertfalse; // slotIndex out of range
}

void QuadrifexRoutingDiagram::setSlotRouting(int newRouting)
{
  ScopedLock scopedLock(*plugInLock);
  slotRouting = newRouting;
  repaint();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void QuadrifexRoutingDiagram::mouseEnter(const MouseEvent &e)
{

  //RWidget::mouseEnter(e);
}

void QuadrifexRoutingDiagram::mouseExit(const MouseEvent &e)
{
  //RWidget::mouseExit(e);
}

void QuadrifexRoutingDiagram::mouseDown(const MouseEvent &e)
{
  ScopedLock scopedLock(*plugInLock);
  int slotIndex = getSlotIndexAtPixelPosition(e.x, e.y);
  if( slotIndex != -1 )
  {
    if( e.mods.isLeftButtonDown() )
    {
      // preliminary - handle it the same way as right-clicks:
      QuadrifexModuleEditor *qme = static_cast<QuadrifexModuleEditor*> (getParentComponent());
      qme->openEffectSelectionMenuForSlot(slotIndex, Point<int>(getX()+e.x, getY()+e.y));

      // toggle between mute, bypass and some effect:
      if(algorithmIndices[slotIndex] == rosic::Quadrifex::MUTE)
        algorithmIndices[slotIndex] = rosic::Quadrifex::BYPASS;
      else if(algorithmIndices[slotIndex] == rosic::Quadrifex::BYPASS)
        algorithmIndices[slotIndex] = oldAlgorithmIndices[slotIndex];
      else {
        oldAlgorithmIndices[slotIndex] = algorithmIndices[slotIndex];
        algorithmIndices[slotIndex] = rosic::Quadrifex::MUTE; }
      repaint();
    }
    else if( e.mods.isRightButtonDown() )
    {
      QuadrifexModuleEditor *qme = static_cast<QuadrifexModuleEditor*> (getParentComponent());
      qme->openEffectSelectionMenuForSlot(slotIndex, Point<int>(getX()+e.x, getY()+e.y));

      /*
      EffectSelectionPopup menu(this);
      menu.show(false, ROwnedPopUpComponent::BELOW);
      int dummy = 0;
      //DEBUG_BREAK;
      */

      // code below needs to be updated to take into account that RPopUpMenu has become the new baseclass of EffectSelectionPopup
      // instead of RPopUpMenuOld
      // old:
      /*
      EffectSelectionPopup menu;
      int algoIndex = menu.show()-1;
      if( algoIndex > -1 )
      algorithmIndices[slotIndex] = algoIndex;
      */
    }
    sendChangeMessage();
  }
}

void QuadrifexRoutingDiagram::mouseDrag(const juce::MouseEvent &e)
{

}

//-------------------------------------------------------------------------------------------------
// intenal helper functions:

int QuadrifexRoutingDiagram::getSlotIndexAtPixelPosition(int x, int y)
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < rosic::Quadrifex::numEffectSlots; i++ )
  {
    if( slotBoxes[i].contains(x, y) == true )
      return i;
  }
  return -1;
}

//-------------------------------------------------------------------------------------------------
// drawing:

void QuadrifexRoutingDiagram::paint(Graphics &g)
{
  ScopedLock scopedLock(*plugInLock);

  //RectangleComponent::paint(g);

  rsPlot::paint(g);
  drawRoutingDiagram(g);
}

void QuadrifexRoutingDiagram::drawRoutingDiagram(Graphics &g)
{
  ScopedLock scopedLock(*plugInLock);
  //int routing = quadrifexToEdit->getSlotRouting();

  //g.setColour(Colours::red); // test

  g.setColour(plotColourScheme.getCurveColour(0));

  // for nice drawing, me may have to make sure that w, h are divisible by 4...
  float w  = (float) getWidth();
  float w2 = w/2;
  float h  = (float) getHeight();
  float h2 = h/2;
  float x  = 4;
  float y  = h2;
  float s  = 20;        // side length of blocks
  float s2 = s/2;
  float s4 = s/4;
  float cx = 0;
  float cy = 0;
  float d  = sqrt((float)s2);  // diagonal distance of arrowhead from circle for adders with arrows at 45°

                               //Colour colour = Colours::black;

  if( slotRouting == rosic::Quadrifex::R_1TO2TO3TO4 )
  {
    x = w2-4*s-s2;

    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 0, algorithmIndices[0]);
    x += s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 1, algorithmIndices[1]);
    x += s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 2, algorithmIndices[2]);
    x += s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 3, algorithmIndices[3]);
    x += s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
  }
  else if( slotRouting == rosic::Quadrifex::R_1TO2TO3_PLUS4 )
  {
    x = w2-4.5f*s-2;

    g.drawLine(x, y, x+s, y, 2.f);
    x += s;
    g.drawLine(x, y-s, x, y+s, 2.f);
    y -= s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 0, algorithmIndices[0]);
    x += s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2,     s, s, 1, algorithmIndices[1]);
    drawSlotBox(g, x, y-s2+2*s, s, s, 3, algorithmIndices[3]);
    x += s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 2, algorithmIndices[2]);
    x += s;
    g.drawLine(x, y, x+s, y, 2.f);
    x += s;

    cx = x;
    cy = h2;
    g.drawArrow(Line<float>(x, y,     x, cy-s4), 2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x, cy+s4), 2.f, 6.f, 6.f);
    y = h2;
    drawBlockDiagramPlus(g, x-s4, y-s4, s2, s2, 2.f);

    y = h2+s;
    x = slotBoxes[0].getX()-s;
    g.drawArrow(Line<float>(x, y, x+3*s, y), 2.f, 6.f, 6.f);
    x += 4*s;
    g.drawLine(x, y, x+3*s, y, 2.f);

    g.drawArrow(Line<float>(cx+s4, cy, cx+s4+s, cy), 2.f, 6.f, 6.f);
  }
  else if( slotRouting == rosic::Quadrifex::R_1TO2_PLUS_3TO4 )
  {
    x = w2-3.5f*s;

    g.drawLine(x, y, x+s, y, 2.f);
    x += s;
    g.drawLine(x, y-s, x, y+s, 2.f);
    y -= s;
    g.drawArrow(Line<float>(x, y,     x+s, y),     2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x+s, y+2*s), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2,     s, s, 0, algorithmIndices[0]);
    drawSlotBox(g, x, y-s2+2*s, s, s, 2, algorithmIndices[2]);
    x += s;
    g.drawArrow(Line<float>(x, y,     x+s, y),     2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x+s, y+2*s), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2,     s, s, 1, algorithmIndices[1]);
    drawSlotBox(g, x, y-s2+2*s, s, s, 3, algorithmIndices[3]);
    x += s;
    g.drawLine(x, y,     x+s, y,     2.f);
    g.drawLine(x, y+2*s, x+s, y+2*s, 2.f);
    x += s;

    cx = x;
    cy = h2;
    g.drawArrow(Line<float>(x, y,     x, cy-s4), 2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x, cy+s4), 2.f, 6.f, 6.f);
    y = h2;
    drawBlockDiagramPlus(g, x-s4, y-s4, s2, s2, 2.f);

    g.drawArrow(Line<float>(cx+s4, cy, cx+s4+s, cy), 2.f, 6.f, 6.f);
  }
  else if( slotRouting == rosic::Quadrifex::R_1PLUS2PLUS3PLUS4 )
  {
    x = w2 - (3*s);

    g.drawLine(x, y, x+s, y, 2.f);
    x += s;
    g.drawLine(x, y-3*s, x, y+3*s, 2.f); // vertical line
    y = h2-3*s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    y += 2*s;

    x += s;
    y = h2-3*s;
    cx = x+2*s;
    cy = h2;
    drawSlotBox(g, x, y-s2, s, s, 0, algorithmIndices[0]);
    g.drawLine(x+s, y, x+2*s, y, 2.f);
    g.drawArrow(Line<float>(x+2*s, y, x+2*s, cy-s4), 2.f, 6.f, 6.f);
    y += 2*s;
    drawSlotBox(g, x, y-s2, s, s, 1, algorithmIndices[1]);
    y += 2*s;
    drawSlotBox(g, x, y-s2, s, s, 2, algorithmIndices[2]);
    y += 2*s;
    drawSlotBox(g, x, y-s2, s, s, 3, algorithmIndices[3]);
    g.drawLine(x+s, y, x+2*s, y, 2.f);
    g.drawArrow(Line<float>(x+2*s, y, x+2*s, cy+s4), 2.f, 6.f, 6.f);
    y += 2*s;

    x  = (float) slotBoxes[1].getRight();
    y  = (float) slotBoxes[1].getY()+s2;
    g.drawArrow(Line<float>(x, y, cx-d, cy-d), 2.f, 6.f, 6.f);

    x  = (float) slotBoxes[2].getRight();
    y  = (float) slotBoxes[2].getY()+s2;
    g.drawArrow(Line<float>(x, y, cx-d, cy+d), 2.f, 6.f, 6.f);

    x = slotBoxes[0].getRight()+s;
    y = h2;
    drawBlockDiagramPlus(g, x-s4, y-s4, s2, s2, 2.f);

    x = cx+s4;
    g.drawArrow(Line<float>(x, cy, x+s, cy), 2.f, 6.f, 6.f);
  }
  else if( slotRouting == rosic::Quadrifex::R_1PLUS2PLUS3_TO_4 )
  {
    x = w2 - (3.5f*s);

    g.drawLine(x, y, x+s, y, 2.f);
    x += s;
    g.drawLine(x, y-2*s, x, y+2*s, 2.f); // vertical line
    y = h2-2*s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);

    x += s;
    y = h2-2*s;
    drawSlotBox(g, x, y-s2, s, s, 0, algorithmIndices[0]);
    y += 2*s;
    drawSlotBox(g, x, y-s2, s, s, 1, algorithmIndices[1]);
    y += 2*s;
    drawSlotBox(g, x, y-s2, s, s, 2, algorithmIndices[2]);

    cx = x+2*s;
    cy = h2;
    drawBlockDiagramPlus(g, cx-s4, cy-s4, s2, s2, 2.f);

    x += s;
    y = h2-2*s;
    g.drawLine( x,   y, x+s,   y,        2.f);
    g.drawArrow(Line<float>(x+s, y, x+s,   y+2*s-s4), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawArrow(Line<float>(x, y, x+s-s4, y), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawLine(x, y, x+s, y, 2.f);
    g.drawArrow(Line<float>(x+s, y, x+s,   y-2*s+s4), 2.f, 6.f, 6.f);

    x += s;
    y = h2;
    g.drawArrow(Line<float>(x+s4, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 3, algorithmIndices[3]);
    x += s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
  }
  else if( slotRouting == rosic::Quadrifex::R_1_TO_2PLUS3_TO_4 )
  {
    x = w2 - (4.5f*s);

    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 0, algorithmIndices[0]);
    x += s;
    g.drawLine(Line<float>(x, y, x+s, y), 2.f);
    x += s;
    g.drawLine(x, y-s, x, y+s, 2.f); // vertical line
    y = h2-s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);

    x += s;
    drawSlotBox(g, x, y-s2, s, s, 2, algorithmIndices[2]);
    y -= 2*s;
    drawSlotBox(g, x, y-s2, s, s, 1, algorithmIndices[1]);

    x += s;
    g.drawLine(Line<float>(x, y,     x+s, y),     2.f);
    g.drawLine(Line<float>(x, y+2*s, x+s, y+2*s), 2.f);
    x += s;
    g.drawArrow(Line<float>(x, y,     x,   y+s-s4),   2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x,   y+s+s4), 2.f, 6.f, 6.f);

    cx = x;
    cy = h2;
    drawBlockDiagramPlus(g, cx-s4, cy-s4, s2, s2, 2.f);

    y = h2;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 3, algorithmIndices[3]);
    x += s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
  }
  else if( slotRouting == rosic::Quadrifex::R_1PLUS2_TO_3TO4 )
  {
    x = w2-4.5f*s-2;

    g.drawLine(x, y, x+s, y, 2.f);
    x += s;
    g.drawLine(x, y-s, x, y+s, 2.f);
    y -= s;
    g.drawArrow(Line<float>(x, y,     x+s, y),     2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x+s, y+2*s), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2,     s, s, 0, algorithmIndices[0]);
    drawSlotBox(g, x, y-s2+2*s, s, s, 1, algorithmIndices[1]);

    x += s;
    g.drawLine(x, y,     x+s, y,     2.f);
    g.drawLine(x, y+2*s, x+s, y+2*s, 2.f);
    x += s;
    cx = x;
    cy = h2;
    g.drawArrow(Line<float>(x, y,     x, cy-s4), 2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x, cy+s4), 2.f, 6.f, 6.f);
    y = h2;
    drawBlockDiagramPlus(g, x-s4, y-s4, s2, s2, 2.f);

    g.drawArrow(Line<float>(x, y,     x+s, y),     2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 2, algorithmIndices[2]);
    x += s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 3, algorithmIndices[3]);
    x += s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
  }
  else if( slotRouting == rosic::Quadrifex::R_1TO2_TO_3PLUS4 )
  {
    x = w2-4.5f*s-2;

    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 0, algorithmIndices[0]);
    x += s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 1, algorithmIndices[1]);
    x += s;
    g.drawLine(Line<float>(x, y, x+s, y), 2.f);

    x += s;
    g.drawLine(x, y-s, x, y+s, 2.f);
    y -= s;
    g.drawArrow(Line<float>(x, y,     x+s, y),     2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x+s, y+2*s), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2,     s, s, 2, algorithmIndices[2]);
    drawSlotBox(g, x, y-s2+2*s, s, s, 3, algorithmIndices[3]);

    x += s;
    g.drawLine(Line<float>(x, y,     x+s, y),     2.f);
    g.drawLine(Line<float>(x, y+2*s, x+s, y+2*s), 2.f);
    x += s;
    cx = x;
    cy = h2;
    g.drawArrow(Line<float>(x, y,     x, cy-s4), 2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x, cy+s4), 2.f, 6.f, 6.f);
    y = h2;
    drawBlockDiagramPlus(g, x-s4, y-s4, s2, s2, 2.f);

    g.drawArrow(Line<float>(cx+s4, cy, cx+s4+s, cy), 2.f, 6.f, 6.f);
  }
  else if( slotRouting == rosic::Quadrifex::R_1PLUS2_TO_3PLUS4 )
  {
    x = w2-4*s-s2;

    g.drawLine(x, y, x+s, y, 2.f);
    x += s;
    g.drawLine(x, y-s, x, y+s, 2.f);
    y -= s;
    g.drawArrow(Line<float>(x, y,     x+s, y),     2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x+s, y+2*s), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2,     s, s, 0, algorithmIndices[0]);
    drawSlotBox(g, x, y-s2+2*s, s, s, 1, algorithmIndices[1]);

    x += s;
    g.drawLine(Line<float>(x, y,     x+s, y),     2.f);
    g.drawLine(Line<float>(x, y+2*s, x+s, y+2*s), 2.f);
    x += s;
    cx = x;
    cy = h2;
    g.drawArrow(Line<float>(x, y,     x, cy-s4), 2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x, cy+s4), 2.f, 6.f, 6.f);
    y = h2;
    drawBlockDiagramPlus(g, x-s4, y-s4, s2, s2, 2.f);

    g.drawLine(x, y, x+s, y, 2.f);
    x += s;
    g.drawLine(x, y-s, x, y+s, 2.f);
    y -= s;
    g.drawArrow(Line<float>(x, y,     x+s, y),     2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x+s, y+2*s), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2,     s, s, 2, algorithmIndices[2]);
    drawSlotBox(g, x, y-s2+2*s, s, s, 3, algorithmIndices[3]);

    x += s;
    g.drawLine(Line<float>(x, y,     x+s, y),     2.f);
    g.drawLine(Line<float>(x, y+2*s, x+s, y+2*s), 2.f);
    x += s;
    cx = x;
    cy = h2;
    g.drawArrow(Line<float>(x, y,     x, cy-s4), 2.f, 6.f, 6.f);
    g.drawArrow(Line<float>(x, y+2*s, x, cy+s4), 2.f, 6.f, 6.f);
    y = h2;
    drawBlockDiagramPlus(g, x-s4, y-s4, s2, s2, 2.f);

    g.drawArrow(Line<float>(cx+s4, cy, cx+s4+s, cy), 2.f, 6.f, 6.f);
  }
  else if( slotRouting == rosic::Quadrifex::R_1_TO_2PLUS3PLUS4 )
  {
    x = w2 - (4.0f*s);

    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    x += s;
    drawSlotBox(g, x, y-s2, s, s, 0, algorithmIndices[0]);
    x += s;

    g.drawLine(x, y, x+s, y, 2.f);
    x += s;
    g.drawLine(x, y-2*s, x, y+2*s, 2.f); // vertical line
    y = h2-2*s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawArrow(Line<float>(x, y, x+s, y), 2.f, 6.f, 6.f);

    x += s;
    y = h2-2*s;
    drawSlotBox(g, x, y-s2, s, s, 1, algorithmIndices[1]);
    y += 2*s;
    drawSlotBox(g, x, y-s2, s, s, 2, algorithmIndices[2]);
    y += 2*s;
    drawSlotBox(g, x, y-s2, s, s, 3, algorithmIndices[3]);

    cx = x+2*s;
    cy = h2;
    drawBlockDiagramPlus(g, cx-s4, cy-s4, s2, s2, 2.f);

    x += s;
    y = h2-2*s;
    g.drawLine( x,   y, x+s,   y,        2.f);
    g.drawArrow(Line<float>(x+s, y, x+s,   y+2*s-s4), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawArrow(Line<float>(x, y, x+s-s4, y), 2.f, 6.f, 6.f);
    y += 2*s;
    g.drawLine(x, y, x+s, y, 2.f);
    g.drawArrow(Line<float>(x+s, y, x+s,   y-2*s+s4), 2.f, 6.f, 6.f);
    x += s;
    y = h2;
    g.drawArrow(Line<float>(x+s4, y, x+2*s-s4, y), 2.f, 6.f, 6.f);
  }
}

void QuadrifexRoutingDiagram::drawSlotBox(Graphics &g, float x, float y, float w, float h,
  int slotIndex, int algoIndex)
{
  ScopedLock scopedLock(*plugInLock);

  //Colour oldColour = g.getCurrentColour();

  const BitmapFontRoundedBoldA10D0* font = &BitmapFontRoundedBoldA10D0::instance;

  slotBoxes[slotIndex].setBounds((int)x, (int)y, (int)w, (int)h);

  //Colour diagramColour = Colours::black;
  Colour diagramColour = plotColourScheme.getCurveColour(0);

  if( algoIndex == rosic::Quadrifex::MUTE )
  {
    g.setColour(diagramColour);
    g.drawRect(x, y, w, h, 2.f);

    g.setColour(Colours::red.darker(0.5f));
    g.drawLine(x,   y, x+w, y+h, 2.f);
    g.drawLine(x+w, y, x,   y+h, 2.f);

    drawBitmapFontText(g, (int)(x+0.5f*w), (int)(y+0.5f*h), juce::String(slotIndex+1),
      font, diagramColour, -1,
      Justification::centred);
  }
  else if( algoIndex == rosic::Quadrifex::BYPASS )
  {
    g.setColour(diagramColour.withMultipliedAlpha(0.25f));
    g.drawRect(x, y, w, h, 2.f);

    drawBitmapFontText(g, (int)(x+0.5f*w), (int)(y+0.5f*h), juce::String(slotIndex+1),
      font, diagramColour.withMultipliedAlpha(0.5f), -1,
      Justification::centred);

    g.setColour(diagramColour);
    g.drawArrow(Line<float>(x, y+0.5f*h, x+w, y+0.5f*h), 2.f, 6.f, 6.f);
  }
  else
  {
    g.setColour(diagramColour);
    g.drawRect(x, y, w, h, 2.f);

    drawBitmapFontText(g, (int)(x+0.5f*w), (int)(y+0.5f*h), juce::String(slotIndex+1),
      font, diagramColour, -1, Justification::centred);
  }

  //g.setColour(oldColour);
}

//=================================================================================================

QuadrifexModuleEditor::QuadrifexModuleEditor(CriticalSection *newPlugInLock,
  QuadrifexAudioModule* newQuadrifexAudioModule) : AudioModuleEditor(newQuadrifexAudioModule)
{
  ScopedLock scopedLock(*lock);
  // if we don't acquire it here, it hangs on opening the GUI ...why?


  setHeadlineStyle(MAIN_HEADLINE);
  jassert(newQuadrifexAudioModule != NULL ); // you must pass a valid module here
  quadrifexModuleToEdit = newQuadrifexAudioModule;

  // remember for toggling:
  for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++ )
  {
    if( quadrifexModuleToEdit == NULL && quadrifexModuleToEdit->wrappedQuadrifex == NULL )
      oldAlgorithmIndices[i] = quadrifexModuleToEdit->wrappedQuadrifex->getEffectAlgorithmIndex(i);
    else
      oldAlgorithmIndices[i] = rosic::Quadrifex::MUTE;
  }

  addWidget( routingLabel = new RTextField( juce::String(("Routing:"))) );
  routingLabel->setDescription(juce::String(("Choose the routing of the 4 effect slots")));
  routingLabel->setDescriptionField(infoField);

  addWidget( routingComboBox = new RComboBox(juce::String(("RoutingComboBox"))) );
  routingComboBox->setDescription(routingLabel->getDescription());
  routingComboBox->setDescriptionField(infoField);
  routingComboBox->registerComboBoxObserver(this);
  routingComboBox->addItem(Quadrifex::R_BYPASS,           ("Bypass")      );
  routingComboBox->addItem(Quadrifex::R_1TO2TO3TO4,       ("1>2>3>4")     );
  routingComboBox->addItem(Quadrifex::R_1TO2TO3_PLUS4,    ("(1>2>3)+4")   );
  routingComboBox->addItem(Quadrifex::R_1TO2_PLUS_3TO4,   ("(1>2)+(3>4)") );
  routingComboBox->addItem(Quadrifex::R_1PLUS2PLUS3PLUS4, ("1+2+3+4")     );
  routingComboBox->addItem(Quadrifex::R_1PLUS2PLUS3_TO_4, ("(1+2+3)>4")   );
  routingComboBox->addItem(Quadrifex::R_1_TO_2PLUS3_TO_4, ("1>(2+3)>4")   );
  routingComboBox->addItem(Quadrifex::R_1PLUS2_TO_3TO4,   ("(1+2)>3>4")   );
  routingComboBox->addItem(Quadrifex::R_1TO2_TO_3PLUS4,   ("1>2>(3+4)")   );
  routingComboBox->addItem(Quadrifex::R_1PLUS2_TO_3PLUS4, ("(1+2)>(3+4)") );
  routingComboBox->addItem(Quadrifex::R_1_TO_2PLUS3PLUS4, ("1>(2+3+4)")   );
  routingComboBox->addItem(Quadrifex::MATRIX,             ("Matrix")      );

  routingDiagram = new QuadrifexRoutingDiagram(lock);
  routingDiagram->addChangeListener(this);
  addPlot(routingDiagram, true, true);
  //addAndMakeVisible(routingDiagram);

  matrixEditor = new RoutingMatrixModuleEditor(lock, quadrifexModuleToEdit->matrixModule);
  matrixEditor->setHeadlineStyle(AudioModuleEditor::NO_HEADLINE);
  addChildEditor(matrixEditor);

  addWidget( dryWetSlider = new RSlider (("DryWetSlider")) );
  dryWetSlider->setSliderName(juce::String(("Dry/Wet")));
  dryWetSlider->assignParameter( quadrifexModuleToEdit->getParameterByName("DryWet") );
  dryWetSlider->setDefaultValue(0.5);
  dryWetSlider->setDescription( juce::String(("Ratio between dry and wet signal")) );
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&ratioToString0);

  addWidget( wetLevelSlider = new RSlider (("WetLevelSlider")) );
  wetLevelSlider->setSliderName(juce::String(("Wet Level")));
  wetLevelSlider->assignParameter( quadrifexModuleToEdit->getParameterByName("WetLevel") );
  wetLevelSlider->setDescription( juce::String(("Level of the wet signal in dB")) );
  wetLevelSlider->setDescriptionField(infoField);
  wetLevelSlider->setStringConversionFunction(&jura::decibelsToStringWithUnit2);

  // initialize the sub-editor pointers:
  for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
    moduleEditors[i] = NULL;

  // create the popup menu to select an effect algorithm:
  effectSelectionPopup = new EffectSelectionPopup(this);
  effectSelectionPopup->registerPopUpMenuObserver(this);
  slotForWhichMenuIsOpen = -1;



  // attach to the underlying audiomodule to be edited:
  quadrifexModuleToEdit->setEditor(this);

  initializeColourScheme();
  updateWidgetsAccordingToState();

  setSize(840, 540);
}

QuadrifexModuleEditor::~QuadrifexModuleEditor()
{
  ScopedLock scopedLock(*lock);


  delete effectSelectionPopup;

  // detach from the underlying audiomodule to be edited:
  quadrifexModuleToEdit->setEditor(NULL);
}

//-------------------------------------------------------------------------------------------------
// setup:

void QuadrifexModuleEditor::initializeColourScheme()
{

}

//-------------------------------------------------------------------------------------------------
// callbacks:

/*
void QuadrifexModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
ScopedLock scopedLock(*lock);

if( quadrifexModuleToEdit == NULL )
return;
if( quadrifexModuleToEdit->wrappedQuadrifex == NULL )
return;

//...

quadrifexModuleToEdit->markStateAsDirty();
}
*/

void QuadrifexModuleEditor::rComboBoxChanged(jura::RComboBox *rComboBoxThatHasChanged)
{
  ScopedLock scopedLock(*lock);
  if( quadrifexModuleToEdit == NULL )
    return;
  if( quadrifexModuleToEdit->wrappedQuadrifex == NULL )
    return;

  rosic::Quadrifex* core = quadrifexModuleToEdit->wrappedQuadrifex;

  if( rComboBoxThatHasChanged == routingComboBox )
  {
    int id = routingComboBox->getSelectedItemIdentifier();
    routingDiagram->setSlotRouting(id);
    core->setSlotRouting(id);
    if( id ==  rosic::Quadrifex::MATRIX )
      matrixEditor->setVisible(true);
    else
      matrixEditor->setVisible(false);
  }

  if( quadrifexModuleToEdit != NULL )
    quadrifexModuleToEdit->markStateAsDirty();
}

void QuadrifexModuleEditor::rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged)
{
  ScopedLock scopedLock(*lock);
  if( quadrifexModuleToEdit == NULL )
    return;
  if( quadrifexModuleToEdit->wrappedQuadrifex == NULL )
    return;

  if( menuThatHasChanged == effectSelectionPopup )
  {
    RTreeViewNode* selectedNode = effectSelectionPopup->getSelectedItem();
    if( selectedNode != NULL )
    {
      int algoIndex = selectedNode->getNodeIdentifier();
      if( algoIndex > -1 )
        setEffectAlgorithm(slotForWhichMenuIsOpen, algoIndex);
    }
  }

  if( quadrifexModuleToEdit != NULL )
    quadrifexModuleToEdit->markStateAsDirty();
}

void QuadrifexModuleEditor::rSliderValueChanged(RSlider *rSliderThatHasChanged)
{
  ScopedLock scopedLock(*lock);
  if( quadrifexModuleToEdit != NULL )
    quadrifexModuleToEdit->markStateAsDirty();
}


void QuadrifexModuleEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  ScopedLock scopedLock(*lock);
  if( objectThatHasChanged == routingDiagram )
  {
    for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
      setEffectAlgorithm(i, routingDiagram->getAlgorithmIndex(i));
    updateWidgetsAccordingToState();
  }
  else
    AudioModuleEditor::changeListenerCallback(objectThatHasChanged);
}

void QuadrifexModuleEditor::mouseDown(const MouseEvent &e)
{
  ScopedLock scopedLock(*lock);
  int slotIndex = getSlotIndexAtPixelPosition(e.x, e.y);
  if( slotIndex != -1 )
  {
    if( e.mods.isRightButtonDown() )
      openEffectSelectionMenuForSlot(slotIndex, Point<int>(e.x, e.y));
    else if( e.mods.isLeftButtonDown() )
    {
      /*
      // on left-click, toggle between mute, bypass and the most recent actual effect - seems not yet to work:
      rosic::Quadrifex *core = quadrifexModuleToEdit->wrappedQuadrifex;
      if( core->getEffectAlgorithmIndex(slotIndex) == rosic::Quadrifex::MUTE )
      core->setEffectAlgorithm(slotIndex, rosic::Quadrifex::BYPASS);
      else if( core->getEffectAlgorithmIndex(slotIndex) == rosic::Quadrifex::BYPASS )
      core->setEffectAlgorithm(slotIndex, oldAlgorithmIndices[slotIndex]);
      else
      {
      oldAlgorithmIndices[slotIndex] = core->getEffectAlgorithmIndex(slotIndex);
      core->setEffectAlgorithm(slotIndex, rosic::Quadrifex::MUTE);
      }
      updateWidgetsAccordingToState();
      */
    }
  }
}

void QuadrifexModuleEditor::updateWidgetsAccordingToState()
{
  ScopedLock scopedLock(*lock);
  if( quadrifexModuleToEdit == NULL )
    return;
  if( quadrifexModuleToEdit->wrappedQuadrifex == NULL )
    return;

  AudioModuleEditor::updateWidgetsAccordingToState();

  int routingIndex = quadrifexModuleToEdit->wrappedQuadrifex->getSlotRouting();
  routingComboBox->selectItemByIndex(routingIndex, false, false);
  routingDiagram->setSlotRouting(routingIndex);
  if( routingIndex ==  rosic::Quadrifex::MATRIX )
    matrixEditor->setVisible(true);
  else
    matrixEditor->setVisible(false);
  for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
  {
    routingDiagram->setAlgorithmIndex(i,
      quadrifexModuleToEdit->wrappedQuadrifex->getEffectAlgorithmIndex(i));
  }

  // delete old and create new editors:
  removeChildEditors();
  for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
    createEditorForSlot(i, quadrifexModuleToEdit->wrappedQuadrifex->getEffectAlgorithmIndex(i));
}

void QuadrifexModuleEditor::paint(Graphics &g)
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::paint(g);

  fillRectWithBilinearGradient(g, globalRectangle, editorColourScheme.topLeft,
    editorColourScheme.topRight, editorColourScheme.bottomLeft, editorColourScheme.bottomRight);

  g.setColour(editorColourScheme.outline);
  g.drawRect(globalRectangle, 2);

  int x = globalRectangle.getX(); // + middleRectangle.getWidth()/2;
  int y = globalRectangle.getY();

  //drawBitmapFontText(g, x+4, y+4, juce::String(("Global Settings")), 
  //  &BitmapFontRoundedBoldA16D0::instance, editorColourScheme.headline);

  if( quadrifexModuleToEdit == NULL )
    return;
  if( quadrifexModuleToEdit->wrappedQuadrifex == NULL )
    return;

  rosic::Quadrifex* core = quadrifexModuleToEdit->wrappedQuadrifex;
  for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
  {
    fillRectWithBilinearGradient(g, slotRectangles[i], editorColourScheme.topLeft,
      editorColourScheme.topRight, editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
    g.drawRect(slotRectangles[i], 2);

    x = slotRectangles[i].getX();
    y = slotRectangles[i].getY();
    juce::String headlineString = juce::String(i+1) + juce::String((" - ")) +
      quadrifexModuleToEdit->effectAlgorithmIndexToString(core->getEffectAlgorithmIndex(i));

    drawBitmapFontText(g, x+4, y+4, headlineString, &BitmapFontRoundedBoldA16D0::instance,
      editorColourScheme.headline);
  }
}

void QuadrifexModuleEditor::resized()
{
  ScopedLock scopedLock(*lock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getHeadlineBottom();
  int w = getWidth();
  int h = getHeight();

  int leftWidth = 200 + w%2;
  h = infoField->getY()-y-2;
  globalRectangle.setBounds(0, y+4, leftWidth, h);

  x = globalRectangle.getX();
  y = globalRectangle.getY();

  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  stateWidgetSet->setBounds(x+4, y+4, globalRectangle.getWidth()-8, 40);
  y = stateWidgetSet->getBottom();

  //int w2 = 100;
  routingLabel->setBounds(x+4, y+4, 52, 16);
  x = routingLabel->getRight();
  w = globalRectangle.getRight()-x;
  routingComboBox->setBounds(x, y+4, w-4, 16);
  y = routingComboBox->getBottom() - RWidget::outlineThickness;
  routingDiagram->setBounds(globalRectangle.getX()+4, y, leftWidth-8, leftWidth-8);
  matrixEditor->setBounds(routingDiagram->getBounds());

  /*
  int w3 = globalRectangle.getWidth()-routingComboBox->getRight();
  */

  x = globalRectangle.getX();
  y = routingDiagram->getBottom();
  w = globalRectangle.getWidth();
  dryWetSlider->setBounds(x+4, y+8, w-8, 16);

  // setup the sizes of the child editors:
  x = globalRectangle.getRight()-RWidget::outlineThickness;
  y = globalRectangle.getY();
  //w = (getWidth()-x-RWidget::outlineThickness)/2;
  w = (getWidth()-x)/2 + 1;
  h = globalRectangle.getHeight()/2 + RWidget::outlineThickness/2;
  for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
  {
    int x2 = x;
    if( i == 1 || i == 3 )
      x2 += w-RWidget::outlineThickness;

    int y2 = y;
    if( i >=2 )
      y2 += h-RWidget::outlineThickness;

    slotRectangles[i].setBounds(x2, y2,    w, h);
    if( moduleEditors[i] != NULL )
      moduleEditors[i]->setBounds(x2, y2+24, w, h-24);
  }
}

//-------------------------------------------------------------------------------------------------
// intenal helper functions:

void QuadrifexModuleEditor::openEffectSelectionMenuForSlot(int slotIndex, juce::Point<int> menuPosition)
{
  ScopedLock scopedLock(*lock);
  if( slotIndex < 0 || slotIndex > 3 )
    return;

  // set up the menu to highlight the currently selected algorithm:
  int selectedIndex = quadrifexModuleToEdit->wrappedQuadrifex->getEffectAlgorithmIndex(slotIndex);
  effectSelectionPopup->selectItemByIdentifier(selectedIndex, false);
  //effectSelectionPopup->openNodeOfSelectedItem();  // \todo: implement openNodeOfSelectedItem in RTreeView

  // keep track of the slot for which the menu is open:
  slotForWhichMenuIsOpen = slotIndex;

  // show the menu:
  Point<int> thisPosistion = getScreenPosition();
  int x = thisPosistion.getX() + menuPosition.getX();
  int y = thisPosistion.getY() + menuPosition.getY();
  //effectSelectionPopup->showAt(false, x, y, 200, 400);
  effectSelectionPopup->showAt(true, x, y, 200, 400); // modal-state must be true, otherwise it's not shown
}

int QuadrifexModuleEditor::getSlotIndexAtPixelPosition(int x, int y)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < rosic::Quadrifex::numEffectSlots; i++ )
  {
    if( slotRectangles[i].contains(x, y) == true )
      return i;
  }
  return -1;
}

void QuadrifexModuleEditor::setEffectAlgorithm(int slotIndex, int newAlgorithmIndex)
{
  ScopedLock scopedLock(*lock);
  if( slotIndex < 0 || slotIndex >= rosic::Quadrifex::numEffectSlots )
  {
    jassertfalse;
    return;
  }

  if( quadrifexModuleToEdit == NULL )
    return;
  else if( quadrifexModuleToEdit->wrappedQuadrifex->getEffectAlgorithmIndex(slotIndex)
    == newAlgorithmIndex )
  {
    return; // nothing to do
  }

  quadrifexModuleToEdit->setEffectAlgorithm(slotIndex, newAlgorithmIndex);
  // the quadrifexModuleToEdit has a pointer to 'this' editor and will take care of updating
  // the GUI (deleting and creating the appropriate child-editor)

  //algorithmIndices[slotIndex] = newAlgorithmIndex;  // maybe later
}

/*
void QuadrifexModuleEditor::createEditorsIfNeeded()
{
ScopedLock scopedLock(*lock);
for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
createEditorForSlotIfNeeded(i);
}

void QuadrifexModuleEditor::createEditorForSlotIfNeeded(int slotIndex)
{

}
*/

void QuadrifexModuleEditor::removeChildEditors()
{
  ScopedLock scopedLock(*lock);
  for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
    removeChildEditorInSlot(i);
}

void QuadrifexModuleEditor::removeChildEditorInSlot(int slotIndex)
{
  ScopedLock scopedLock(*lock);
  int i = slotIndex;
  if( slotIndex >= 0 && slotIndex < rosic::Quadrifex::numEffectSlots
    && moduleEditors[i] != NULL  )
  {
    moduleEditors[i]->invalidateModulePointer();  // so it doesn't dereference it in the destructor
    removeChildEditor(moduleEditors[i], true);    // deletes the object also
    moduleEditors[i] = NULL;
  }
}

void QuadrifexModuleEditor::createEditorForSlot(int slotIndex, int algorithmIndex)
{
  ScopedLock scopedLock(*lock);
  if( quadrifexModuleToEdit == NULL )
    return;
  if( quadrifexModuleToEdit->wrappedQuadrifex == NULL )
    return;

  jassert( moduleEditors[slotIndex] == NULL );
  // you should delete the old editor and null the pointer before creating a new one

  switch( algorithmIndex )
  {
  case rosic::Quadrifex::BIT_CRUSHER:
  {
    jura::BitCrusherAudioModule *audioModule = static_cast<jura::BitCrusherAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new BitCrusherModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::CHORUS:
  {
    jura::ChorusAudioModule *audioModule = static_cast<jura::ChorusAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new ChorusModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::COMB_BANK:
  {
    jura::CombBankAudioModule *audioModule = static_cast<jura::CombBankAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new CombBankModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::COMB_RESONATOR:
  {
    jura::CombResonatorAudioModule *audioModule = static_cast<jura::CombResonatorAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new CombResonatorModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::COMB_STEREOIZER:
  {
    jura::CombStereoizerAudioModule *audioModule = static_cast<jura::CombStereoizerAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new CombStereoizerModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::COMPRESSOR:
  {
    jura::CompressorAudioModule *audioModule = static_cast<jura::CompressorAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new CompressorModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::DUAL_TWO_POLE_FILTER:
  {
    jura::DualTwoPoleFilterAudioModule *audioModule = static_cast<jura::DualTwoPoleFilterAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new DualTwoPoleFilterModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::EQUALIZER:
  {
    jura::EqualizerAudioModule *audioModule = static_cast<jura::EqualizerAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    EqualizerModuleEditor *editor = new EqualizerModuleEditor(lock, audioModule);
    editor->setLayout(EqualizerModuleEditor::SLIDERS_BELOW);
    editor->setUseShortSliderNames(true);
    editor->setUseSmallComboBox(true);
    moduleEditors[slotIndex] = editor;
  } break;
  case rosic::Quadrifex::EXPANDER:
  {
    jura::ExpanderAudioModule *audioModule = static_cast<jura::ExpanderAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    ExpanderModuleEditor *editor = new ExpanderModuleEditor(lock, audioModule);
    moduleEditors[slotIndex] = editor;
  } break;
  case rosic::Quadrifex::FLANGER:
  {
    jura::FlangerAudioModule *audioModule = static_cast<jura::FlangerAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new FlangerModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::FORMANT_SHIFTER:
  {
    jura::FormantShifterAudioModule *audioModule = static_cast<jura::FormantShifterAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new FormantShifterModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::FOUR_POLE_FILTER:
  {
    jura::FourPoleFilterAudioModule *audioModule = static_cast<jura::FourPoleFilterAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new FourPoleFilterModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::FREQUENCY_SHIFTER:
  {
    jura::FrequencyShifterAudioModule *audioModule = static_cast<jura::FrequencyShifterAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new FrequencyShifterModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::HARMONICS:
  {
    jura::HarmonicsAudioModule *audioModule = static_cast<jura::HarmonicsAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new HarmonicsModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::LADDER_FILTER:
  {
    jura::LadderFilterAudioModule *audioModule = static_cast<jura::LadderFilterAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new LadderFilterModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::LIMITER:
  {
    jura::LimiterAudioModule *audioModule = static_cast<jura::LimiterAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new LimiterModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::MODULATED_ALLPASS:
  {
    jura::ModulatedAllpassAudioModule *audioModule = static_cast<jura::ModulatedAllpassAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new ModulatedAllpassModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::NOISE_GATE:
  {
    jura::NoiseGateAudioModule *audioModule = static_cast<jura::NoiseGateAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new NoiseGateModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::NOISIFIER:
  {
    jura::NoisifierAudioModule *audioModule = static_cast<jura::NoisifierAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new NoisifierModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::PHASER:
  {
    jura::PhaserAudioModule *audioModule = static_cast<jura::PhaserAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new PhaserModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::PHASE_STEREOIZER:
  {
    jura::PhaseStereoizerAudioModule *audioModule = static_cast<jura::PhaseStereoizerAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new PhaseStereoizerModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::PINGPONG_ECHO:
  {
    jura::PingPongEchoAudioModule *audioModule = static_cast<jura::PingPongEchoAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new PingPongEchoModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::PITCH_SHIFTER:
  {
    jura::PitchShifterAudioModule *audioModule = static_cast<jura::PitchShifterAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new PitchShifterModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::REVERB:
  {
    jura::ReverbAudioModule *audioModule = static_cast<jura::ReverbAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new ReverbModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::RINGMODULATOR:
  {
    jura::RingModulatorAudioModule *audioModule = static_cast<jura::RingModulatorAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new RingModulatorModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::SIMPLE_DELAY:
  {
    jura::SimpleDelayAudioModule *audioModule = static_cast<jura::SimpleDelayAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new SimpleDelayModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::SINE_OSCILLATOR:
  {
    jura::SineOscillatorAudioModule *audioModule = static_cast<jura::SineOscillatorAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new SineOscillatorModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::SSB_MODULATOR:
  {
    jura::SingleSidebandModulatorAudioModule *audioModule = static_cast<jura::SingleSidebandModulatorAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new SingleSidebandModulatorModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::SLEWRATE_LIMITER:
  {
    jura::SlewRateLimiterAudioModule *audioModule = static_cast<jura::SlewRateLimiterAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new SlewRateLimiterModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::SLOPE_FILTER:
  {
    jura::SlopeFilterAudioModule *audioModule = static_cast<jura::SlopeFilterAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new SlopeFilterModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::STEREO_PAN:
  {
    jura::StereoPanAudioModule *audioModule = static_cast<jura::StereoPanAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new StereoPanModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::STEREO_WIDTH:
  {
    jura::StereoWidthAudioModule *audioModule = static_cast<jura::StereoWidthAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new StereoWidthModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::TREMOLO:
  {
    jura::TremoloAudioModule *audioModule = static_cast<jura::TremoloAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new TremoloModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::TWO_POLE_FILTER:
  {
    jura::TwoPoleFilterAudioModule *audioModule = static_cast<jura::TwoPoleFilterAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new TwoPoleFilterModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::VIBRATO:
  {
    jura::VibratoAudioModule *audioModule = static_cast<jura::VibratoAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new VibratoModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::WAH_WAH:
  {
    jura::WahWahAudioModule *audioModule = static_cast<jura::WahWahAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new WahWahModuleEditor(lock, audioModule);
  } break;
  case rosic::Quadrifex::WAVESHAPER:
  {
    jura::WaveShaperAudioModule *audioModule = static_cast<jura::WaveShaperAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new WaveShaperModuleEditor(lock, audioModule);
  } break;


  // ------> INSERT NEW CASE-MARK HERE WHEN ADDING A NEW ALGORITHM <------


  case rosic::Quadrifex::MUTE:
  {
    jura::MuteAudioModule *audioModule = static_cast<jura::MuteAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new MuteModuleEditor(lock, audioModule);
  } break;
  default:
  {
    jura::BypassAudioModule *audioModule = static_cast<jura::BypassAudioModule*>
      (quadrifexModuleToEdit->getEffectAudioModule(slotIndex));
    moduleEditors[slotIndex] = new BypassModuleEditor(lock, audioModule);
  }
  }

  routingDiagram->setAlgorithmIndex(slotIndex, algorithmIndex);
  moduleEditors[slotIndex]->setHeadlineStyle(Editor::NO_HEADLINE);
  moduleEditors[slotIndex]->setLinkPosition(AudioModuleEditor::INVISIBLE);
  moduleEditors[slotIndex]->setDescriptionField(infoField, true);
  addChildEditor(moduleEditors[slotIndex]);
  moduleEditors[slotIndex]->copyColourSettingsFrom(this);
  resized();  // to adjust the bounds of the new editor
  repaint();  // to redraw to headline
  moduleEditors[slotIndex]->updateWidgetsAccordingToState();
  setupPopupEditors(slotIndex);

  if( setupDialog != NULL )
    setupDialog->toFront(false);

}

void QuadrifexModuleEditor::setupPopupEditors(int slotIndex)
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor *editor = moduleEditors[slotIndex];
  if( editor == NULL )
    return;

  ModulationEffectModuleEditor* modulationEditor
    = dynamic_cast<ModulationEffectModuleEditor*> (editor);
  if( modulationEditor != NULL )
  {

    //Rectangle r = modulationEditor->lfoEditor->editButton->getBounds();


    int x, y;
    switch( slotIndex )
    {
    case 0: { x = -240; y = 16;   } break;
    case 1: { x = -240; y = 16;   } break;
    case 2: { x = -240; y = -200; } break;
    case 3: { x = -240; y = -200; } break;
    }
    int w = 400;
    int h = 200;


    //modulationEditor->lfoEditor->setPopupEditorBounds(...)
    modulationEditor->setLfoPopUpEditorBounds(x, y, w, h);
    //int dummy = 0;
  }

}

/*
Bugs:
-set the routing to 1+2+3+4 and mute effect 4 - some parts of the graphics turn red that shouldn't

*/
