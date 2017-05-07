//#include "rosof_QuadrifexAudioModule.h"
#include "rosof_QuadrifexModuleEditor.h"
using namespace rosof;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

QuadrifexAudioModule::QuadrifexAudioModule(CriticalSection *newPlugInLock, rosic::Quadrifex *quadrifexToWrap)
: AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(quadrifexToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedQuadrifex = quadrifexToWrap;
  editor           = NULL;
  moduleName       = juce::String(T("Quadrifex"));
  setActiveDirectory(getApplicationDirectory() + juce::String(T("/QuadrifexPresets")) );

  matrixModule = new RoutingMatrixAudioModule(plugInLock, &wrappedQuadrifex->mixMatrix);
  matrixModule->setModuleName(juce::String(T("RoutingMatrix5x5")));
  addChildAudioModule(matrixModule);

  for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
  {
    // create a bypass-module for each slot:
    rosic::BypassModule* bypassCoreModule = 
      static_cast<rosic::BypassModule*> (wrappedQuadrifex->getEffectModule(i)); // this dynamic_cast causes bugs in the release version
    rosof::BypassAudioModule *audioModule = new rosof::BypassAudioModule(plugInLock, bypassCoreModule);
    effectModules[i] = audioModule;
    addChildAudioModule(effectModules[i]);

    // allocate memory to store the states internally:
    bitCrusherStates[i]              = new XmlElement(juce::String(T("BitCrusher")));
    chorusStates[i]                  = new XmlElement(juce::String(T("Chorus")));
    combBankStates[i]                = new XmlElement(juce::String(T("CombBank")));
    combResonatorStates[i]           = new XmlElement(juce::String(T("CombResonator")));
    combStereoizerStates[i]          = new XmlElement(juce::String(T("CombStereoizer")));
    compressorStates[i]              = new XmlElement(juce::String(T("Compressor")));
    dualTwoPoleFilterStates[i]       = new XmlElement(juce::String(T("DualTwoPoleFilter")));
    equalizerStates[i]               = new XmlElement(juce::String(T("Equalizer")));
    expanderStates[i]                = new XmlElement(juce::String(T("Expander")));
    flangerStates[i]                 = new XmlElement(juce::String(T("Flanger")));
    formantShifterStates[i]          = new XmlElement(juce::String(T("FormantShifter")));
    fourPoleFilterStates[i]          = new XmlElement(juce::String(T("FourPoleFilter")));
    frequencyShifterStates[i]        = new XmlElement(juce::String(T("FrequencyShifter")));
    harmonicsStates[i]               = new XmlElement(juce::String(T("Harmonics")));
    ladderFilterStates[i]            = new XmlElement(juce::String(T("LadderFilter")));
    limiterStates[i]                 = new XmlElement(juce::String(T("Limiter")));
    modulatedAllpassStates[i]        = new XmlElement(juce::String(T("ModulatedAllpass")));
    noiseGateStates[i]               = new XmlElement(juce::String(T("NoiseGate")));
    noisifierStates[i]               = new XmlElement(juce::String(T("Noisifier")));
    phaserStates[i]                  = new XmlElement(juce::String(T("Phaser")));
    phaseStereoizerStates[i]         = new XmlElement(juce::String(T("PhaseStereoizer")));
    pingPongEchoStates[i]            = new XmlElement(juce::String(T("PingPongEcho")));
    pitchShifterStates[i]            = new XmlElement(juce::String(T("PitchShifter")));
    reverbStates[i]                  = new XmlElement(juce::String(T("Reverb")));
    ringModulatorStates[i]           = new XmlElement(juce::String(T("RingModulator")));
    simpleDelayStates[i]             = new XmlElement(juce::String(T("SimpleDelay")));
    sineOscillatorStates[i]          = new XmlElement(juce::String(T("SineOscillator")));
    singleSidebandModulatorStates[i] = new XmlElement(juce::String(T("SingleSidebandModulator")));
    slewRateLimiterStates[i]         = new XmlElement(juce::String(T("SlewRateLimiter")));
    slopeFilterStates[i]             = new XmlElement(juce::String(T("SlopeFilter")));
    stereoPanStates[i]               = new XmlElement(juce::String(T("StereoPan")));
    stereoWidthStates[i]             = new XmlElement(juce::String(T("StereoWidth")));
    tremoloStates[i]                 = new XmlElement(juce::String(T("Tremolo")));
    twoPoleFilterStates[i]           = new XmlElement(juce::String(T("TwoPoleFilter")));
    vibratoStates[i]                 = new XmlElement(juce::String(T("Vibrato")));
    wahWahStates[i]                  = new XmlElement(juce::String(T("WahWah")));
    waveShaperStates[i]              = new XmlElement(juce::String(T("WaveShaper")));
  }

  initializeAutomatableParameters();
}

QuadrifexAudioModule::~QuadrifexAudioModule()
{
  ScopedLock scopedLock(*plugInLock);

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

//-------------------------------------------------------------------------------------------------
// setup:

void QuadrifexAudioModule::setEditor(QuadrifexModuleEditor *newEditor)
{
  ScopedLock scopedLock(*plugInLock);
  editor = newEditor;
}

void QuadrifexAudioModule::setEffectAlgorithm(int slotIndex, int newAlgorithmIndex)
{
  ScopedLock scopedLock(*plugInLock);

  if( wrappedQuadrifex == NULL )
    return;

  if( slotIndex < 0 || slotIndex >= rosic::Quadrifex::numEffectSlots )
    return;

  // store the state of the old effect to be replaced:
  int oldAlgorithmIndex = wrappedQuadrifex->getEffectAlgorithmIndex(slotIndex);
  switch( oldAlgorithmIndex )
  {
  case rosic::Quadrifex::BIT_CRUSHER: 
    {
      rosof::BitCrusherAudioModule *audioModule = 
        static_cast<rosof::BitCrusherAudioModule*> (effectModules[slotIndex]);
      delete bitCrusherStates[slotIndex];
      bitCrusherStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::CHORUS: 
    {
      rosof::ChorusAudioModule *audioModule = 
        static_cast<rosof::ChorusAudioModule*> (effectModules[slotIndex]);
      delete chorusStates[slotIndex];
      chorusStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::COMB_BANK: 
    {
      rosof::CombBankAudioModule *audioModule = 
        static_cast<rosof::CombBankAudioModule*> (effectModules[slotIndex]);
      delete combBankStates[slotIndex];
      combBankStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::COMB_RESONATOR: 
    {
      rosof::CombResonatorAudioModule *audioModule = 
        static_cast<rosof::CombResonatorAudioModule*> (effectModules[slotIndex]);
      delete combResonatorStates[slotIndex];
      combResonatorStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::COMB_STEREOIZER: 
    {
      rosof::CombStereoizerAudioModule *audioModule = 
        static_cast<rosof::CombStereoizerAudioModule*> (effectModules[slotIndex]);
      delete combStereoizerStates[slotIndex];
      combStereoizerStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::COMPRESSOR: 
    {
      rosof::CompressorAudioModule *audioModule = 
        static_cast<rosof::CompressorAudioModule*> (effectModules[slotIndex]);
      delete compressorStates[slotIndex];
      compressorStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::DUAL_TWO_POLE_FILTER: 
    {
      rosof::DualTwoPoleFilterAudioModule *audioModule = 
        static_cast<rosof::DualTwoPoleFilterAudioModule*> (effectModules[slotIndex]);
      delete dualTwoPoleFilterStates[slotIndex];
      dualTwoPoleFilterStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::EQUALIZER: 
    {
      rosof::EqualizerAudioModule *audioModule = 
        static_cast<rosof::EqualizerAudioModule*> (effectModules[slotIndex]);
      delete equalizerStates[slotIndex];
      equalizerStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::EXPANDER: 
    {
      rosof::ExpanderAudioModule *audioModule = 
        static_cast<rosof::ExpanderAudioModule*> (effectModules[slotIndex]);
      delete expanderStates[slotIndex];
      expanderStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::FLANGER: 
    {
      rosof::FlangerAudioModule *audioModule = 
        static_cast<rosof::FlangerAudioModule*> (effectModules[slotIndex]);
      delete flangerStates[slotIndex];
      flangerStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::FORMANT_SHIFTER: 
    {
      rosof::FormantShifterAudioModule *audioModule = 
        static_cast<rosof::FormantShifterAudioModule*> (effectModules[slotIndex]);
      delete formantShifterStates[slotIndex];
      formantShifterStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::FOUR_POLE_FILTER: 
    {
      rosof::FourPoleFilterAudioModule *audioModule = 
        static_cast<rosof::FourPoleFilterAudioModule*> (effectModules[slotIndex]);
      delete fourPoleFilterStates[slotIndex];
      fourPoleFilterStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::FREQUENCY_SHIFTER: 
    {
      rosof::FrequencyShifterAudioModule *audioModule = 
        static_cast<rosof::FrequencyShifterAudioModule*> (effectModules[slotIndex]);
      delete frequencyShifterStates[slotIndex];
      frequencyShifterStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::HARMONICS: 
    {
      rosof::HarmonicsAudioModule *audioModule = 
        static_cast<rosof::HarmonicsAudioModule*> (effectModules[slotIndex]);
      delete harmonicsStates[slotIndex];
      harmonicsStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::LADDER_FILTER: 
    {
      rosof::LadderFilterAudioModule *audioModule = 
        static_cast<rosof::LadderFilterAudioModule*> (effectModules[slotIndex]);
      delete ladderFilterStates[slotIndex];
      ladderFilterStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::LIMITER: 
    {
      rosof::LimiterAudioModule *audioModule = 
        static_cast<rosof::LimiterAudioModule*> (effectModules[slotIndex]);
      delete limiterStates[slotIndex];
      limiterStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::MODULATED_ALLPASS: 
    {
      rosof::ModulatedAllpassAudioModule *audioModule = 
        static_cast<rosof::ModulatedAllpassAudioModule*> (effectModules[slotIndex]);
      delete modulatedAllpassStates[slotIndex];
      modulatedAllpassStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::NOISE_GATE: 
    {
      rosof::NoiseGateAudioModule *audioModule = 
        static_cast<rosof::NoiseGateAudioModule*> (effectModules[slotIndex]);
      delete noiseGateStates[slotIndex];
      noiseGateStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::NOISIFIER: 
    {
      rosof::NoisifierAudioModule *audioModule = 
        static_cast<rosof::NoisifierAudioModule*> (effectModules[slotIndex]);
      delete noisifierStates[slotIndex];
      noisifierStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::PHASER: 
    {
      rosof::PhaserAudioModule *audioModule = 
        static_cast<rosof::PhaserAudioModule*> (effectModules[slotIndex]);
      delete phaserStates[slotIndex];
      phaserStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::PHASE_STEREOIZER: 
    {
      rosof::PhaseStereoizerAudioModule *audioModule = 
        static_cast<rosof::PhaseStereoizerAudioModule*> (effectModules[slotIndex]);
      delete phaseStereoizerStates[slotIndex];
      phaseStereoizerStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::PINGPONG_ECHO: 
    {
      rosof::PingPongEchoAudioModule *audioModule = 
        static_cast<rosof::PingPongEchoAudioModule*> (effectModules[slotIndex]);
      delete pingPongEchoStates[slotIndex];
      pingPongEchoStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::PITCH_SHIFTER: 
    {
      rosof::PitchShifterAudioModule *audioModule = 
        static_cast<rosof::PitchShifterAudioModule*> (effectModules[slotIndex]);
      delete pitchShifterStates[slotIndex];
      pitchShifterStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::REVERB: 
    {
      rosof::ReverbAudioModule *audioModule = 
        static_cast<rosof::ReverbAudioModule*> (effectModules[slotIndex]);
      delete reverbStates[slotIndex];
      reverbStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::RINGMODULATOR: 
    {
      rosof::RingModulatorAudioModule *audioModule = 
        static_cast<rosof::RingModulatorAudioModule*> (effectModules[slotIndex]);
      delete ringModulatorStates[slotIndex];
      ringModulatorStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::SIMPLE_DELAY: 
    {
      rosof::SimpleDelayAudioModule *audioModule = 
        static_cast<rosof::SimpleDelayAudioModule*> (effectModules[slotIndex]);
      delete simpleDelayStates[slotIndex];
      simpleDelayStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::SINE_OSCILLATOR: 
    {
      rosof::SineOscillatorAudioModule *audioModule = 
        static_cast<rosof::SineOscillatorAudioModule*> (effectModules[slotIndex]);
      delete sineOscillatorStates[slotIndex];
      sineOscillatorStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::SLOPE_FILTER: 
    {
      rosof::SlopeFilterAudioModule *audioModule = 
        static_cast<rosof::SlopeFilterAudioModule*> (effectModules[slotIndex]);
      delete slopeFilterStates[slotIndex];
      slopeFilterStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::SSB_MODULATOR: 
    {
      rosof::SingleSidebandModulatorAudioModule *audioModule = 
        static_cast<rosof::SingleSidebandModulatorAudioModule*> (effectModules[slotIndex]);
      delete singleSidebandModulatorStates[slotIndex];
      singleSidebandModulatorStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::SLEWRATE_LIMITER: 
    {
      rosof::SlewRateLimiterAudioModule *audioModule = 
        static_cast<rosof::SlewRateLimiterAudioModule*> (effectModules[slotIndex]);
      delete slewRateLimiterStates[slotIndex];
      slewRateLimiterStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::STEREO_PAN: 
    {
      rosof::StereoPanAudioModule *audioModule = 
        static_cast<rosof::StereoPanAudioModule*> (effectModules[slotIndex]);
      delete stereoPanStates[slotIndex];
      stereoPanStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::STEREO_WIDTH: 
    {
      rosof::StereoWidthAudioModule *audioModule = 
        static_cast<rosof::StereoWidthAudioModule*> (effectModules[slotIndex]);
      delete stereoWidthStates[slotIndex];
      stereoWidthStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::TREMOLO: 
    {
      rosof::TremoloAudioModule *audioModule = 
        static_cast<rosof::TremoloAudioModule*> (effectModules[slotIndex]);
      delete tremoloStates[slotIndex];
      tremoloStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::TWO_POLE_FILTER: 
    {
      rosof::TwoPoleFilterAudioModule *audioModule = 
        static_cast<rosof::TwoPoleFilterAudioModule*> (effectModules[slotIndex]);
      delete twoPoleFilterStates[slotIndex];
      twoPoleFilterStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::VIBRATO: 
    {
      rosof::VibratoAudioModule *audioModule = 
        static_cast<rosof::VibratoAudioModule*> (effectModules[slotIndex]);
      delete vibratoStates[slotIndex];
      vibratoStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::WAH_WAH: 
    {
      rosof::WahWahAudioModule *audioModule = 
        static_cast<rosof::WahWahAudioModule*> (effectModules[slotIndex]);
      delete wahWahStates[slotIndex];
      wahWahStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;
  case rosic::Quadrifex::WAVESHAPER: 
    {
      rosof::WaveShaperAudioModule *audioModule = 
        static_cast<rosof::WaveShaperAudioModule*> (effectModules[slotIndex]);
      delete waveShaperStates[slotIndex];
      waveShaperStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
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
      rosof::BitCrusherAudioModule *audioModule = new rosof::BitCrusherAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("BitCrusher")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*bitCrusherStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::CHORUS: 
    {
      rosic::ChorusModule *core = 
        static_cast<rosic::ChorusModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::ChorusAudioModule *audioModule = new rosof::ChorusAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Chorus")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*chorusStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::COMB_BANK: 
    {
      rosic::CombBankModule *core = 
        static_cast<rosic::CombBankModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::CombBankAudioModule *audioModule = new rosof::CombBankAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("CombBank")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*combBankStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::COMB_RESONATOR: 
    {
      rosic::CombResonatorStereoModule *core = 
        static_cast<rosic::CombResonatorStereoModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::CombResonatorAudioModule *audioModule = new rosof::CombResonatorAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("CombResonator")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*combResonatorStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::COMB_STEREOIZER: 
    {
      rosic::CombStereoizerModule *core = 
        static_cast<rosic::CombStereoizerModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::CombStereoizerAudioModule *audioModule = new rosof::CombStereoizerAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("CombStereoizer")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*combStereoizerStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::COMPRESSOR: 
    {
      rosic::SoftKneeCompressorModule *core = 
        static_cast<rosic::SoftKneeCompressorModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::CompressorAudioModule *audioModule = new rosof::CompressorAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Compressor")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*compressorStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::DUAL_TWO_POLE_FILTER: 
    {
      rosic::DualTwoPoleFilterModule *core = 
        static_cast<rosic::DualTwoPoleFilterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::DualTwoPoleFilterAudioModule *audioModule = new rosof::DualTwoPoleFilterAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("DualTwoPoleFilter")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*dualTwoPoleFilterStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::EQUALIZER: 
    {
      rosic::EqualizerModule *core = 
        static_cast<rosic::EqualizerModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::EqualizerAudioModule *audioModule = new rosof::EqualizerAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Equalizer")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*equalizerStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::EXPANDER: 
    {
      rosic::SoftKneeExpanderModule *core = 
        static_cast<rosic::SoftKneeExpanderModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::ExpanderAudioModule *audioModule = new rosof::ExpanderAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Expander")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*expanderStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::FLANGER: 
    {
      rosic::FlangerModule *core = 
        static_cast<rosic::FlangerModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::FlangerAudioModule *audioModule = new rosof::FlangerAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Flanger")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*flangerStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::FORMANT_SHIFTER: 
    {
      rosic::FormantShifterModule *core = 
        static_cast<rosic::FormantShifterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::FormantShifterAudioModule *audioModule = new rosof::FormantShifterAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("FormantShifter")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*formantShifterStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::FOUR_POLE_FILTER: 
    {
      rosic::FourPoleFilterModule *core = 
        static_cast<rosic::FourPoleFilterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::FourPoleFilterAudioModule *audioModule = new rosof::FourPoleFilterAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("FourPoleFilter")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*fourPoleFilterStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::FREQUENCY_SHIFTER: 
    {
      rosic::FrequencyShifterStereoModule *core = 
        static_cast<rosic::FrequencyShifterStereoModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::FrequencyShifterAudioModule *audioModule = new rosof::FrequencyShifterAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("FrequencyShifter")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*frequencyShifterStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::HARMONICS: 
    {
      rosic::HarmonicsModule *core = 
        static_cast<rosic::HarmonicsModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::HarmonicsAudioModule *audioModule = new rosof::HarmonicsAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Harmonics")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*harmonicsStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::LADDER_FILTER: 
    {
      rosic::LadderFilterModule *core = 
        static_cast<rosic::LadderFilterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::LadderFilterAudioModule *audioModule = new rosof::LadderFilterAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("LadderFilter")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*ladderFilterStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::LIMITER: 
    {
      rosic::LimiterModule *core = 
        static_cast<rosic::LimiterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::LimiterAudioModule *audioModule = new rosof::LimiterAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Limiter")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*limiterStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::MODULATED_ALLPASS: 
    {
      rosic::ModulatedAllpassModule *core = 
        static_cast<rosic::ModulatedAllpassModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::ModulatedAllpassAudioModule *audioModule = new rosof::ModulatedAllpassAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("ModulatedAllpass")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*modulatedAllpassStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::NOISE_GATE: 
    {
      rosic::NoiseGateModule *core = 
        static_cast<rosic::NoiseGateModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::NoiseGateAudioModule *audioModule = new rosof::NoiseGateAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("NoiseGate")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*noiseGateStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::NOISIFIER: 
    {
      rosic::NoisifierModule *core = 
        static_cast<rosic::NoisifierModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::NoisifierAudioModule *audioModule = new rosof::NoisifierAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Noisifier")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*noisifierStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::PHASER: 
    {
      rosic::PhaserModule *core = 
        static_cast<rosic::PhaserModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::PhaserAudioModule *audioModule = new rosof::PhaserAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Phaser")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*phaserStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::PHASE_STEREOIZER: 
    {
      rosic::PhaseStereoizerModule *core = 
        static_cast<rosic::PhaseStereoizerModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::PhaseStereoizerAudioModule *audioModule = new rosof::PhaseStereoizerAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("PhaseStereoizer")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*phaseStereoizerStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::PINGPONG_ECHO: 
    {
      rosic::PingPongEchoModule *core = 
        static_cast<rosic::PingPongEchoModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::PingPongEchoAudioModule *audioModule = new rosof::PingPongEchoAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("PingPongEcho")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*pingPongEchoStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::PITCH_SHIFTER: 
    {
      rosic::PitchShifterModule *core = 
        static_cast<rosic::PitchShifterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::PitchShifterAudioModule *audioModule = new rosof::PitchShifterAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("PitchShifter")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*pitchShifterStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::REVERB: 
    {
      rosic::ReverbModule *core = 
        static_cast<rosic::ReverbModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::ReverbAudioModule *audioModule = new rosof::ReverbAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Reverb")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*reverbStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::RINGMODULATOR: 
    {
      rosic::RingModulatorModule *core = 
        static_cast<rosic::RingModulatorModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::RingModulatorAudioModule *audioModule = new rosof::RingModulatorAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("RingModulator")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*ringModulatorStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::SIMPLE_DELAY: 
    {
      rosic::SimpleDelayModule *core = 
        static_cast<rosic::SimpleDelayModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::SimpleDelayAudioModule *audioModule = new rosof::SimpleDelayAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("SimpleDelay")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*simpleDelayStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::SINE_OSCILLATOR: 
    {
      rosic::SineOscillatorModule *core = 
        static_cast<rosic::SineOscillatorModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::SineOscillatorAudioModule *audioModule = new rosof::SineOscillatorAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("SineOscillator")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*sineOscillatorStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::SLOPE_FILTER: 
    {
      rosic::SlopeFilterModule *core = 
        static_cast<rosic::SlopeFilterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::SlopeFilterAudioModule *audioModule = new rosof::SlopeFilterAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("SlopeFilter")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*slopeFilterStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::SSB_MODULATOR: 
    {
      rosic::SingleSidebandModulatorModule *core = 
        static_cast<rosic::SingleSidebandModulatorModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::SingleSidebandModulatorAudioModule *audioModule = new rosof::SingleSidebandModulatorAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("SingleSidebandModulator")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*singleSidebandModulatorStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::SLEWRATE_LIMITER: 
    {
      rosic::SlewRateLimiterStereoModule *core = 
        static_cast<rosic::SlewRateLimiterStereoModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::SlewRateLimiterAudioModule *audioModule = new rosof::SlewRateLimiterAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("SlewRateLimiter")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*slewRateLimiterStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::STEREO_PAN: 
    {
      rosic::StereoPanModule *core = 
        static_cast<rosic::StereoPanModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::StereoPanAudioModule *audioModule = new rosof::StereoPanAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("StereoPan")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*stereoPanStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::STEREO_WIDTH: 
    {
      rosic::StereoWidthModule *core = 
        static_cast<rosic::StereoWidthModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::StereoWidthAudioModule *audioModule = new rosof::StereoWidthAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("StereoWidth")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*stereoWidthStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::TREMOLO: 
    {
      rosic::TremoloModule *core = 
        static_cast<rosic::TremoloModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::TremoloAudioModule *audioModule = new rosof::TremoloAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Tremolo")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*tremoloStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::TWO_POLE_FILTER: 
    {
      rosic::TwoPoleFilterModule *core = 
        static_cast<rosic::TwoPoleFilterModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::TwoPoleFilterAudioModule *audioModule = new rosof::TwoPoleFilterAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("TwoPoleFilter")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*twoPoleFilterStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::VIBRATO: 
    {
      rosic::VibratoModule *core = 
        static_cast<rosic::VibratoModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::VibratoAudioModule *audioModule = new rosof::VibratoAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("Vibrato")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*vibratoStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::WAH_WAH: 
    {
      rosic::WahWahModule *core = 
        static_cast<rosic::WahWahModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::WahWahAudioModule *audioModule = new rosof::WahWahAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("WahWah")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*wahWahStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  case rosic::Quadrifex::WAVESHAPER: 
    {
      rosic::WaveShaperModule *core = 
        static_cast<rosic::WaveShaperModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::WaveShaperAudioModule *audioModule = new rosof::WaveShaperAudioModule(plugInLock, core);
      audioModule->setModuleName(juce::String(T("WaveShaper")) + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*waveShaperStates[slotIndex], juce::String::empty, true);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;


    // special trivial cases:
  case rosic::Quadrifex::MUTE: 
    {
      rosic::MuteModule *core = 
        static_cast<rosic::MuteModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::MuteAudioModule *audioModule = new rosof::MuteAudioModule(plugInLock, core);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  default: // bypass by default (i.e. value out of range)
    {
      rosic::BypassModule *core = 
        static_cast<rosic::BypassModule*> (wrappedQuadrifex->getEffectModule(slotIndex));
      rosof::BypassAudioModule *audioModule = new rosof::BypassAudioModule(plugInLock, core);
      effectModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
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
  ScopedLock scopedLock(*plugInLock);

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
  ScopedLock scopedLock(*plugInLock);

  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  if( wrappedQuadrifex != NULL )
  {
    // store the routing and slot-effect assignments:
    xmlState->setAttribute(T("Routing"), 
      slotRoutingIndexToString(wrappedQuadrifex->getSlotRouting()));
    for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
    {
      xmlState->setAttribute(juce::String(T("Slot"))+juce::String(i+1), 
        effectAlgorithmIndexToString(wrappedQuadrifex->getEffectAlgorithmIndex(i)) );
    }
  }
  return xmlState;
}

void QuadrifexAudioModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
                                               bool markAsClean)
{
  ScopedLock scopedLock(*plugInLock);

  if( wrappedQuadrifex != NULL )
  {
    // recall the routing and slot-effect assignments:
    wrappedQuadrifex->setSlotRouting(
      stringToSlotRoutingIndex( xmlState.getStringAttribute(T("Routing"), T("1>2>3>4")) ) );
    for(int i=0; i<rosic::Quadrifex::numEffectSlots; i++)
    {
      setEffectAlgorithm(i, stringToEffectAlgorithmIndex( 
        xmlState.getStringAttribute( juce::String(T("Slot"))+juce::String(i+1), T("Mute"))));
    }
  }
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void QuadrifexAudioModule::initializeAutomatableParameters()
{
  ScopedLock scopedLock(*plugInLock);


  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  juce::Array<double> defaultValues;

  // this pointer will be used to temporarily store the addresses of the created Parameter-objects:
  AutomatableParameter* p;

  // #00:
  p = new AutomatableParameter(plugInLock, "DryWet", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR); 
  addObservedParameter(p);

  // #01:
  p = new AutomatableParameter(plugInLock, "WetLevel", -36.0, 6.0, 0.01, 0.0, Parameter::LINEAR); 
  addObservedParameter(p);

  // #02:
  p = new AutomatableParameter(plugInLock, "TriggerInterval", 0.0, 64.0, 1.0, 8.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) observedParameters.size(); i++ )
    parameterChanged(observedParameters[i]);
}

juce::String QuadrifexAudioModule::slotRoutingIndexToString(int index)
{
  ScopedLock scopedLock(*plugInLock);

  switch( index )
  {
  case rosic::Quadrifex::R_BYPASS:            return juce::String(T("Bypass"));
  case rosic::Quadrifex::R_1TO2TO3TO4:        return juce::String(T("1>2>3>4"));
  case rosic::Quadrifex::R_1TO2TO3_PLUS4:     return juce::String(T("(1>2>3)+4"));
  case rosic::Quadrifex::R_1TO2_PLUS_3TO4:    return juce::String(T("(1>2)+(3>4)"));
  case rosic::Quadrifex::R_1PLUS2PLUS3PLUS4:  return juce::String(T("1+2+3+4"));
  case rosic::Quadrifex::R_1PLUS2PLUS3_TO_4:  return juce::String(T("(1+2+3)>4"));
  case rosic::Quadrifex::R_1_TO_2PLUS3_TO_4:  return juce::String(T("1>(2+3)>4"));
  case rosic::Quadrifex::R_1PLUS2_TO_3TO4:    return juce::String(T("(1+2)>3>4"));
  case rosic::Quadrifex::R_1TO2_TO_3PLUS4:    return juce::String(T("1>2>(3+4)"));
  case rosic::Quadrifex::R_1PLUS2_TO_3PLUS4:  return juce::String(T("(1+2)>(3+4)"));
  case rosic::Quadrifex::R_1_TO_2PLUS3PLUS4:  return juce::String(T("1>(2+3+4)"));
  case rosic::Quadrifex::MATRIX:              return juce::String(T("Matrix"));
  default:                                    return juce::String(T("Bypass"));
  }
}

int QuadrifexAudioModule::stringToSlotRoutingIndex(const juce::String &routingString)
{
  ScopedLock scopedLock(*plugInLock);

  if( routingString == juce::String(T("Bypass"))   )       return rosic::Quadrifex::R_BYPASS;
  if( routingString == juce::String(T("1>2>3>4"))  )       return rosic::Quadrifex::R_1TO2TO3TO4;
  if( routingString == juce::String(T("(1>2>3)+4"))  )     return rosic::Quadrifex::R_1TO2TO3_PLUS4;
  if( routingString == juce::String(T("(1>2)+(3>4)"))  )   return rosic::Quadrifex::R_1TO2_PLUS_3TO4;
  if( routingString == juce::String(T("1+2+3+4"))  )       return rosic::Quadrifex::R_1PLUS2PLUS3PLUS4;
  if( routingString == juce::String(T("(1+2+3)>4"))  )     return rosic::Quadrifex::R_1PLUS2PLUS3_TO_4;
  if( routingString == juce::String(T("1>(2+3)>4"))  )     return rosic::Quadrifex::R_1_TO_2PLUS3_TO_4;
  if( routingString == juce::String(T("(1+2)>3>4"))  )     return rosic::Quadrifex::R_1PLUS2_TO_3TO4;
  if( routingString == juce::String(T("1>2>(3+4)"))  )     return rosic::Quadrifex::R_1TO2_TO_3PLUS4;
  if( routingString == juce::String(T("(1+2)>(3+4)"))  )   return rosic::Quadrifex::R_1PLUS2_TO_3PLUS4;
  if( routingString == juce::String(T("1>(2+3+4)"))  )     return rosic::Quadrifex::R_1_TO_2PLUS3PLUS4;
  if( routingString == juce::String(T("Matrix"))  )        return rosic::Quadrifex::MATRIX;

  return rosic::Quadrifex::R_BYPASS;
}

juce::String QuadrifexAudioModule::effectAlgorithmIndexToString(int index)
{
  ScopedLock scopedLock(*plugInLock);

  switch( index )
  {
  case rosic::Quadrifex::MUTE:                 return juce::String(T("Mute"));
  case rosic::Quadrifex::BYPASS:               return juce::String(T("Bypass"));
  case rosic::Quadrifex::BIT_CRUSHER:          return juce::String(T("BitCrusher"));
  case rosic::Quadrifex::HARMONICS:            return juce::String(T("Harmonics"));
  case rosic::Quadrifex::CHORUS:               return juce::String(T("Chorus"));
  case rosic::Quadrifex::COMB_BANK:            return juce::String(T("CombBank"));
  case rosic::Quadrifex::COMB_RESONATOR:       return juce::String(T("CombResonator"));
  case rosic::Quadrifex::COMB_STEREOIZER:      return juce::String(T("CombStereoizer"));
  case rosic::Quadrifex::COMPRESSOR:           return juce::String(T("Compressor"));
  //case rosic::Quadrifex::COMP_SHAPER:          return juce::String(T("CompShaper"));
  case rosic::Quadrifex::DUAL_TWO_POLE_FILTER: return juce::String(T("DualTwoPoleFilter"));
  case rosic::Quadrifex::EQUALIZER:            return juce::String(T("Equalizer"));
  case rosic::Quadrifex::EXPANDER:             return juce::String(T("Expander"));
  case rosic::Quadrifex::FLANGER:              return juce::String(T("Flanger"));
  case rosic::Quadrifex::FORMANT_SHIFTER:      return juce::String(T("FormantShifter"));
  case rosic::Quadrifex::FOUR_POLE_FILTER:     return juce::String(T("FourPoleFilter"));
  case rosic::Quadrifex::FREQUENCY_SHIFTER:    return juce::String(T("FrequencyShifter"));
  case rosic::Quadrifex::LADDER_FILTER:        return juce::String(T("LadderFilter"));
  case rosic::Quadrifex::LIMITER:              return juce::String(T("Limiter"));
  case rosic::Quadrifex::MODULATED_ALLPASS:    return juce::String(T("ModulatedAllpass"));
  case rosic::Quadrifex::NOISE_GATE:           return juce::String(T("NoiseGate"));
  case rosic::Quadrifex::NOISIFIER:            return juce::String(T("Noisifier"));
  case rosic::Quadrifex::PHASER:               return juce::String(T("Phaser"));
  case rosic::Quadrifex::PHASE_STEREOIZER:     return juce::String(T("PhaseStereoizer"));
  case rosic::Quadrifex::PINGPONG_ECHO:        return juce::String(T("PingPongEcho")); 
  case rosic::Quadrifex::PITCH_SHIFTER:        return juce::String(T("PitchShifter")); 
  case rosic::Quadrifex::REVERB:               return juce::String(T("Reverb"));
  case rosic::Quadrifex::RINGMODULATOR:        return juce::String(T("RingModulator"));
  case rosic::Quadrifex::SIMPLE_DELAY:         return juce::String(T("SimpleDelay"));
  case rosic::Quadrifex::SINE_OSCILLATOR:      return juce::String(T("SineOscillator"));
  case rosic::Quadrifex::SSB_MODULATOR:        return juce::String(T("SingleSidebandModulator"));
  case rosic::Quadrifex::SLEWRATE_LIMITER:     return juce::String(T("SlewRateLimiter"));
  case rosic::Quadrifex::SLOPE_FILTER:         return juce::String(T("SlopeFilter"));
  case rosic::Quadrifex::STEREO_PAN:           return juce::String(T("StereoPan"));
  case rosic::Quadrifex::STEREO_WIDTH:         return juce::String(T("StereoWidth"));
  case rosic::Quadrifex::TREMOLO:              return juce::String(T("Tremolo"));
  case rosic::Quadrifex::TWO_POLE_FILTER:      return juce::String(T("TwoPoleFilter"));
  case rosic::Quadrifex::VIBRATO:              return juce::String(T("Vibrato"));
  case rosic::Quadrifex::WAH_WAH:              return juce::String(T("WahWah"));
  case rosic::Quadrifex::WAVESHAPER:           return juce::String(T("WaveShaper"));

  default:                                     return juce::String(T("Mute"));
  }
}

int QuadrifexAudioModule::stringToEffectAlgorithmIndex(const juce::String &algoString)
{
  ScopedLock scopedLock(*plugInLock);

  if( algoString == juce::String(T("Mute"))   )                   return rosic::Quadrifex::MUTE;
  if( algoString == juce::String(T("Bypass")) )                   return rosic::Quadrifex::BYPASS;
  if( algoString == juce::String(T("BitCrusher")) )               return rosic::Quadrifex::BIT_CRUSHER;
  if( algoString == juce::String(T("Harmonics")) )                return rosic::Quadrifex::HARMONICS;
  if( algoString == juce::String(T("Chorus")) )                   return rosic::Quadrifex::CHORUS;
  if( algoString == juce::String(T("CombBank")) )                 return rosic::Quadrifex::COMB_BANK;
  if( algoString == juce::String(T("CombResonator")) )            return rosic::Quadrifex::COMB_RESONATOR;
  if( algoString == juce::String(T("CombStereoizer")) )           return rosic::Quadrifex::COMB_STEREOIZER;
  if( algoString == juce::String(T("Compressor")) )               return rosic::Quadrifex::COMPRESSOR;
  //if( algoString == juce::String(T("CompShaper")) )               return rosic::Quadrifex::COMP_SHAPER;
  if( algoString == juce::String(T("DualTwoPoleFilter")) )        return rosic::Quadrifex::DUAL_TWO_POLE_FILTER;
  if( algoString == juce::String(T("Equalizer")) )                return rosic::Quadrifex::EQUALIZER;
  if( algoString == juce::String(T("Expander")) )                 return rosic::Quadrifex::EXPANDER;
  if( algoString == juce::String(T("Flanger")) )                  return rosic::Quadrifex::FLANGER;
  if( algoString == juce::String(T("FormantShifter")) )           return rosic::Quadrifex::FORMANT_SHIFTER;
  if( algoString == juce::String(T("FourPoleFilter")) )           return rosic::Quadrifex::FOUR_POLE_FILTER;
  if( algoString == juce::String(T("FrequencyShifter")) )         return rosic::Quadrifex::FREQUENCY_SHIFTER;
  if( algoString == juce::String(T("LadderFilter")) )             return rosic::Quadrifex::LADDER_FILTER;
  if( algoString == juce::String(T("Limiter")) )                  return rosic::Quadrifex::LIMITER;
  if( algoString == juce::String(T("ModulatedAllpass")) )         return rosic::Quadrifex::MODULATED_ALLPASS;
  if( algoString == juce::String(T("NoiseGate")) )                return rosic::Quadrifex::NOISE_GATE;
  if( algoString == juce::String(T("Noisifier")) )                return rosic::Quadrifex::NOISIFIER;
  if( algoString == juce::String(T("Phaser")) )                   return rosic::Quadrifex::PHASER;
  if( algoString == juce::String(T("PhaseStereoizer")) )          return rosic::Quadrifex::PHASE_STEREOIZER;
  if( algoString == juce::String(T("PingPongEcho")) )             return rosic::Quadrifex::PINGPONG_ECHO;
  if( algoString == juce::String(T("PitchShifter")) )             return rosic::Quadrifex::PITCH_SHIFTER;
  if( algoString == juce::String(T("Reverb")) )                   return rosic::Quadrifex::REVERB;
  if( algoString == juce::String(T("RingModulator")) )            return rosic::Quadrifex::RINGMODULATOR;
  if( algoString == juce::String(T("SimpleDelay")) )              return rosic::Quadrifex::SIMPLE_DELAY;
  if( algoString == juce::String(T("SineOscillator")) )           return rosic::Quadrifex::SINE_OSCILLATOR;
  if( algoString == juce::String(T("SingleSidebandModulator")) )  return rosic::Quadrifex::SSB_MODULATOR;
  if( algoString == juce::String(T("SlewRateLimiter")) )          return rosic::Quadrifex::SLEWRATE_LIMITER;
  if( algoString == juce::String(T("SlopeFilter")) )              return rosic::Quadrifex::SLOPE_FILTER;
  if( algoString == juce::String(T("StereoPan")) )                return rosic::Quadrifex::STEREO_PAN;
  if( algoString == juce::String(T("StereoWidth")) )              return rosic::Quadrifex::STEREO_WIDTH;
  if( algoString == juce::String(T("Tremolo")) )                  return rosic::Quadrifex::TREMOLO;
  if( algoString == juce::String(T("TwoPoleFilter")) )            return rosic::Quadrifex::TWO_POLE_FILTER;
  if( algoString == juce::String(T("Vibrato")) )                  return rosic::Quadrifex::VIBRATO;
  if( algoString == juce::String(T("WahWah")) )                   return rosic::Quadrifex::WAH_WAH;
  if( algoString == juce::String(T("WaveShaper")) )               return rosic::Quadrifex::WAVESHAPER;

  return rosic::Quadrifex::MUTE;
}