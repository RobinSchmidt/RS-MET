#ifndef rosof_VariousModulesAndEditors_h
#define rosof_VariousModulesAndEditors_h
// rename to VariousEffects, override createEditor everywhere, let Quadrifex use an
// AudioModuleFactory and get rid of that ugly, messy code there related to selecting
// which module to plug in

/** This file defines a bunch of smaller AudioModules and their editors.

mmm... lots of boilerplate code here - is it perhaps possible to avoid some of it by using 
templates (maybe with explicit specializations)?
-maybe we can avoid the verbose constructors of the editors by having a means of auto-generating 
 widgets -> Parameter class needs description- and stringConversionFunction member (refactored 
 from RWidget/RSlider)
-maybe introduce some class ParameterGroup which corresponds to the labels that we have in some 
 classes (like curveLabel, timeLabel, ...) */

//=================================================================================================
// Trivial 'Effects':

//-------------------------------------------------------------------------------------------------
// Bypass:

class BypassAudioModule : public AudioModule
{
public:
  BypassAudioModule(CriticalSection *newPlugInLock, rosic::BypassModule *newBypassToWrap = nullptr)
    : AudioModule(newPlugInLock) {}
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override {}
  juce_UseDebuggingNewOperator;
};

class BypassModuleEditor : public AudioModuleEditor
{
public:
  BypassModuleEditor(CriticalSection *newPlugInLock, BypassAudioModule* newBypassAudioModule)
    :AudioModuleEditor(newBypassAudioModule)
  { }
  juce_UseDebuggingNewOperator;
};

//-----------------------------------------------------------------------------------------------
// Mute:

class MuteAudioModule : public AudioModule
{
public:
  MuteAudioModule(CriticalSection *newPlugInLock, rosic::MuteModule *newMuteToWrap) 
    : AudioModule(newPlugInLock)  {}
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    for(int i = 0; i < numChannels; i++)
      fillWithZeros(inOutBuffer[i], numSamples);
  }
  //virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  //{
  //  *inOutL = *inOutR = 0.0;
  //}
  juce_UseDebuggingNewOperator;
};

class MuteModuleEditor : public AudioModuleEditor
{
public:
  MuteModuleEditor(CriticalSection *newPlugInLock, MuteAudioModule* newMuteAudioModule)
    :AudioModuleEditor(newMuteAudioModule)
  { }
  juce_UseDebuggingNewOperator;
};

//===============================================================================================
// Distortion Effects:

//-----------------------------------------------------------------------------------------------
// BitCrusher:

class BitCrusherAudioModule : public ModulatableAudioModule
{
public:
  BitCrusherAudioModule(CriticalSection *newPlugInLock, rosic::BitCrusher *newBitCrusherToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::BitCrusher *wrappedBitCrusher;
};

class BitCrusherModuleEditor : public AudioModuleEditor
{
public:
  BitCrusherModuleEditor(CriticalSection *newPlugInLock, 
    BitCrusherAudioModule* newBitCrusherAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *decimationSlider, *quantizationSlider, *amountSlider;
};

//-----------------------------------------------------------------------------------------------
// Harmonics:

class HarmonicsAudioModule : public ModulatableAudioModule
{
public:
  HarmonicsAudioModule(CriticalSection *newPlugInLock, rosic::Harmonics *newHarmonicsToWrap);
  virtual void parameterChanged(Parameter* parameterThatHasChanged); // maybe to be deprecated
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::Harmonics *wrappedHarmonics;
};

class HarmonicsModuleEditor : public AudioModuleEditor
{
public:
  HarmonicsModuleEditor(CriticalSection *newPlugInLock, 
    HarmonicsAudioModule* newHarmonicsAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  juce::Rectangle<int> globalRect, harmonicsRect;
  RTextField *globalLabel, *harmonicsLabel, *inFilterLabel, *outFilterLabel;
  rsModulatableSlider *driveSlider, *dryWetSlider, *inHighpassSlider, *inLowpassSlider,
    *outHighpassSlider, *outLowpassSlider;
  RSlider    *h02Slider, *h03Slider, *h04Slider, *h05Slider, *h06Slider, *h07Slider, *h08Slider,
    *h09Slider, *h10Slider, *h11Slider, *h12Slider;

};

//-----------------------------------------------------------------------------------------------
// ModulatedAllpass:

class ModulatedAllpassAudioModule : public ModulatableAudioModule
{
public:
  ModulatedAllpassAudioModule(CriticalSection *newPlugInLock, rosic::ModulatedAllpass *newModulatedAllpassToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::ModulatedAllpass *wrappedModulatedAllpass;
};

class ModulatedAllpassModuleEditor : public AudioModuleEditor
{
public:
  ModulatedAllpassModuleEditor(CriticalSection *newPlugInLock, ModulatedAllpassAudioModule* newModulatedAllpassAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *factorSlider, *offsetSlider;
};

//-----------------------------------------------------------------------------------------------
// SlewRateLimiter:

class SlewRateLimiterAudioModule : public ModulatableAudioModule
{
public:
  SlewRateLimiterAudioModule(CriticalSection *newPlugInLock, rosic::SlewRateLimiterStereo *newSlewRateLimiterToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::SlewRateLimiterStereo *wrappedSlewRateLimiter;
};

class SlewRateLimiterModuleEditor : public AudioModuleEditor
{
public:
  SlewRateLimiterModuleEditor(CriticalSection *newPlugInLock, SlewRateLimiterAudioModule* newSlewRateLimiterAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *attackSlider, *releaseSlider;
};

//-----------------------------------------------------------------------------------------------
// WaveShaper:

class WaveShaperAudioModule : public ModulatableAudioModule
{
  friend class WaveShaperModuleEditor;
public:
  WaveShaperAudioModule(CriticalSection *newPlugInLock, rosic::WaveShaper *newWaveShaperToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::WaveShaper *wrappedWaveShaper;
};

class WaveShaperModuleEditor : public AudioModuleEditor, public RSliderListener, public RComboBoxObserver
{
public:
  WaveShaperModuleEditor(CriticalSection *newPlugInLock, WaveShaperAudioModule* newWaveShaperAudioModule);
  virtual ~WaveShaperModuleEditor();
  virtual void resized();
  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void rSliderValueChanged(RSlider *rSliderThatHasChanged);
  virtual void updateWidgetsAccordingToState();
  virtual void updateWidgetEnablement();
  virtual void updatePlot();
  juce_UseDebuggingNewOperator;
protected:
  WaveShaperAudioModule *waveShaperModuleToEdit;
  RComboBox *curveComboBox;
  rsModulatableSlider *driveSlider, *dcSlider, *amountSlider, *outputLevelSlider;
  RSlider  *oversamplingSlider, *slopeSlider, *interceptSlider;
  rsDataPlot *plot;
  double *xValues, *yValues;
  int    numValues;
};

//-----------------------------------------------------------------------------------------------
// CompShaper:

class CompShaperAudioModule : public ModulatableAudioModule
{
public:
  CompShaperAudioModule(CriticalSection *newPlugInLock, rosic::CompShaper *newCompShaperToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::CompShaper *wrappedCompShaper;
};

class CompShaperModuleEditor : public AudioModuleEditor
{
public:
  CompShaperModuleEditor(CriticalSection *newPlugInLock, CompShaperAudioModule* newCompShaperAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  juce::Rectangle<int> curveParametersRect, timeParametersRect, otherParametersRect;
  RTextField *curveLabel, *timeLabel, *othersLabel;
  rsModulatableSlider *driveSlider, *outLevelSlider, *amountSlider, *thresholdSlider, *ratioSlider, 
    *kneeSlider;
  RButton    *clipButton;
};

//===============================================================================================
// Dynamics:

//-----------------------------------------------------------------------------------------------
// Compressor:

class CompressorAudioModule : public ModulatableAudioModule // make baseclass DynamicsAudioModule
{
public:
  CompressorAudioModule(CriticalSection *newPlugInLock, 
    rosic::SoftKneeCompressor *newCompressorToWrap = nullptr);
  virtual ~CompressorAudioModule();

  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void reset() override;
  AudioModuleEditor* createEditor(int type) override;

protected:
  virtual void createStaticParameters();
  rosic::SoftKneeCompressor *wrappedCompressor;
  bool wrappedCompressorIsOwned = false;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(CompressorAudioModule)
};

class CompressorModuleEditor : public AudioModuleEditor
{
public:
  CompressorModuleEditor(CriticalSection *newPlugInLock, CompressorAudioModule* newCompressorAudioModule);
  virtual void resized() override;
protected:
  juce::Rectangle<int> curveParametersRect, timeParametersRect, otherParametersRect;
  RTextField *curveLabel, *timeLabel, *othersLabel;
  rsModulatableSlider *attackSlider, *releaseSlider, *lookAheadSlider, *inLevelSlider, 
    *outLevelSlider, *dryWetSlider, *thresholdSlider, *ratioSlider, *kneeSlider;
  RButton *autoGainButton, *limitButton, *antiAliasButton;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(CompressorModuleEditor)
};

//-----------------------------------------------------------------------------------------------
// Expander:

class ExpanderAudioModule : public ModulatableAudioModule
{
public:
  ExpanderAudioModule(CriticalSection *newPlugInLock, rosic::SoftKneeExpander *newExpanderToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::SoftKneeExpander *wrappedExpander;
};

class ExpanderModuleEditor : public AudioModuleEditor
{
public:
  ExpanderModuleEditor(CriticalSection *newPlugInLock, ExpanderAudioModule* newExpanderAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  juce::Rectangle<int> curveParametersRect, timeParametersRect, otherParametersRect;
  RTextField *curveLabel, *timeLabel, *othersLabel;
  rsModulatableSlider *attackSlider, *releaseSlider, *lookAheadSlider, *inLevelSlider, 
    *outLevelSlider, *dryWetSlider, *thresholdSlider, *ratioSlider, *kneeSlider;
  RButton    *gateButton;
};

//-----------------------------------------------------------------------------------------------
// Limiter:

class LimiterAudioModule : public ModulatableAudioModule
{
public:
  LimiterAudioModule(CriticalSection *newPlugInLock, rosic::Limiter* newLimiterToWrap = nullptr);
  virtual ~LimiterAudioModule() { if(wrappedLimiterIsOwned) delete wrappedLimiter; }
  virtual AudioModuleEditor* createEditor(int type) override;
  virtual void processStereoFrame(double *left, double *right) override 
  { 
    wrappedLimiter->getSampleFrameStereo(left, right); 
  }
protected:
  virtual void createParameters();
  rosic::Limiter *wrappedLimiter;
  bool wrappedLimiterIsOwned = false;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LimiterAudioModule)
};

class LimiterModuleEditor : public AudioModuleEditor
{
public:
  LimiterModuleEditor(CriticalSection *newPlugInLock, LimiterAudioModule* newLimiterAudioModule);
  virtual void resized();
protected:
  virtual void createWidgets();
  juce::Rectangle<int> curveParametersRect, timeParametersRect, otherParametersRect;
  RTextField *curveLabel, *timeLabel, *othersLabel;
  rsModulatableSlider *attackSlider, *releaseSlider, *lookAheadSlider, *inLevelSlider, 
    *outLevelSlider, *dryWetSlider, *limitSlider;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LimiterModuleEditor)
};

//-----------------------------------------------------------------------------------------------
// NoiseGate:

class NoiseGateAudioModule : public ModulatableAudioModule
{
public:
  NoiseGateAudioModule(CriticalSection *newPlugInLock, rosic::NoiseGate *newNoiseGateToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::NoiseGate *wrappedNoiseGate;
};

class NoiseGateModuleEditor : public AudioModuleEditor
{
public:
  NoiseGateModuleEditor(CriticalSection *newPlugInLock, NoiseGateAudioModule* newNoiseGateAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  juce::Rectangle<int> curveParametersRect, timeParametersRect, otherParametersRect;
  RTextField *curveLabel, *timeLabel, *othersLabel;
  rsModulatableSlider *attackSlider, *holdSlider, *releaseSlider, *lookAheadSlider, *inLevelSlider, 
    *outLevelSlider, *dryWetSlider, *thresholdSlider, *hysteresisSlider;
};


//===============================================================================================
// Filters:

//-----------------------------------------------------------------------------------------------
// CombBank:

class CombBankAudioModule : public ModulatableAudioModule
{
public:
  CombBankAudioModule(CriticalSection *newPlugInLock, rosic::CombBank *newCombBankToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::CombBank *wrappedCombBank;
};

class CombBankModuleEditor : public AudioModuleEditor
{
public:
  CombBankModuleEditor(CriticalSection *newPlugInLock, CombBankAudioModule* newCombBankAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  juce::Rectangle<int> toneParametersRect, decayParametersRect; //, otherParametersRect;
  RTextField *toneLabel, *decayLabel; //, *othersLabel;
  rsModulatableSlider *dryWetSlider, *levelSlider, *frequencySlider, *detuneSlider, *pan1Slider, 
    *pan2Slider, *decayTimeSlider,*highDecayScaleSlider, *lowDecayScaleSlider, *highFreqSlider, 
    *lowFreqSlider;
  RButton *oddOnlyButton;
};

//-----------------------------------------------------------------------------------------------
// CombResonator:

class CombResonatorAudioModule : public ModulatableAudioModule
{
public:
  CombResonatorAudioModule(CriticalSection *newPlugInLock, rosic::CombResonatorStereo *newCombResonatorToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::CombResonatorStereo *wrappedCombResonator;
};

class CombResonatorModuleEditor : public AudioModuleEditor
{
public:
  CombResonatorModuleEditor(CriticalSection *newPlugInLock, CombResonatorAudioModule* newCombResonatorAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  juce::Rectangle<int> toneParametersRect, decayParametersRect; //, otherParametersRect;
  RTextField *toneLabel, *decayLabel; //, *othersLabel;
  rsModulatableSlider *dryWetSlider, *levelSlider, *frequencySlider, *detuneSlider, *pan1Slider, 
    *pan2Slider, *decayTimeSlider, *highDecayScaleSlider, *lowDecayScaleSlider, *highFreqSlider, 
    *lowFreqSlider;
  RButton    *oddOnlyButton;
};

//-----------------------------------------------------------------------------------------------
// DualTwoPoleFilter:

class DualTwoPoleFilterAudioModule : public ModulatableAudioModule
{
  friend class DualTwoPoleFilterModuleEditor;
public:
  DualTwoPoleFilterAudioModule(CriticalSection *newPlugInLock, rosic::DualTwoPoleFilter *newDualTwoPoleFilterToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::DualTwoPoleFilter *wrappedDualTwoPoleFilter;
};

class DualTwoPoleFilterModuleEditor : public AudioModuleEditor, public RComboBoxObserver
{
public:
  DualTwoPoleFilterModuleEditor(CriticalSection *newPlugInLock, DualTwoPoleFilterAudioModule* newDualTwoPoleFilterAudioModule);
  virtual void resized();
  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void updateWidgetEnablement();
  juce_UseDebuggingNewOperator;
protected:
  DualTwoPoleFilterAudioModule *twoPoleFilterModuleToEdit;
  juce::Rectangle<int> globalRect, filter1Rect, filter2Rect;
  RTextField *globalLabel, *filter1Label, *filter2Label;
  RComboBox  *modeComboBox1, *modeComboBox2;
  rsModulatableSlider *frequencySlider1, *gainSlider1, *bandwidthSlider1, *frequencySlider2, 
    *gainSlider2, *bandwidthSlider2, *serialParallelBlendSlider, *frequencyScaleSlider, 
    *gainScaleSlider, *bandwidthScaleSlider;
};

//-----------------------------------------------------------------------------------------------
// FourPoleFilter:

class FourPoleFilterAudioModule : public ModulatableAudioModule
{
  friend class FourPoleFilterModuleEditor;
public:
  FourPoleFilterAudioModule(CriticalSection *newPlugInLock, rosic::FourPoleFilter *newFourPoleFilterToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::FourPoleFilter *wrappedFourPoleFilter;
};

class FourPoleFilterModuleEditor : public AudioModuleEditor, public RComboBoxObserver
{
public:
  FourPoleFilterModuleEditor(CriticalSection *newPlugInLock, FourPoleFilterAudioModule* newFourPoleFilterAudioModule);
  virtual void resized();
  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void updateWidgetEnablement();
  juce_UseDebuggingNewOperator;
protected:
  FourPoleFilterAudioModule *fourPoleFilterModuleToEdit;
  //RComboBox *modeComboBox;
  FourPoleFilterModeComboBox *modeComboBox;
  rsModulatableSlider *frequencySlider, *gainSlider, *bandwidthSlider;
};

//-----------------------------------------------------------------------------------------------
// LadderFilter:

class LadderFilterAudioModule : public ModulatableAudioModule
{
  friend class LadderFilterModuleEditor;
public:
  LadderFilterAudioModule(CriticalSection *newPlugInLock, rosic::LadderFilterOld *newLadderFilterToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::LadderFilterOld *wrappedLadderFilter;
};

class LadderFilterModuleEditor : public AudioModuleEditor, public RComboBoxObserver
{
public:
  LadderFilterModuleEditor(CriticalSection *newPlugInLock, LadderFilterAudioModule* newLadderFilterAudioModule);
  virtual void resized();
  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void updateWidgetEnablement();
  juce_UseDebuggingNewOperator;
protected:
  LadderFilterAudioModule *ladderFilterModuleToEdit;
  RComboBox *modeComboBox;
  rsModulatableSlider *frequencySlider, *resonanceSlider, *makeUpSlider, *driveSlider, *orderSlider,
    *morphSlider;
};

//-----------------------------------------------------------------------------------------------
// SlopeFilter:

class SlopeFilterAudioModule : public ModulatableAudioModule
{
  friend class SlopeFilterModuleEditor;
public:
  SlopeFilterAudioModule(CriticalSection *newPlugInLock, rosic::SlopeFilter *newSlopeFilterToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::SlopeFilter *wrappedSlopeFilter;
};

class SlopeFilterModuleEditor : public AudioModuleEditor
{
public:
  SlopeFilterModuleEditor(CriticalSection *newPlugInLock, SlopeFilterAudioModule* newSlopeFilterAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  SlopeFilterAudioModule *slopeFilterModuleToEdit;
  rsModulatableSlider *slopeSlider;
};

//-----------------------------------------------------------------------------------------------
// TwoPoleFilter:

class TwoPoleFilterAudioModule : public ModulatableAudioModule
{
  friend class TwoPoleFilterModuleEditor;
public:
  TwoPoleFilterAudioModule(CriticalSection *newPlugInLock, rosic::TwoPoleFilter *newTwoPoleFilterToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::TwoPoleFilter *wrappedTwoPoleFilter;
};

class TwoPoleFilterModuleEditor : public AudioModuleEditor, public RComboBoxObserver
{
public:
  TwoPoleFilterModuleEditor(CriticalSection *newPlugInLock, TwoPoleFilterAudioModule* newTwoPoleFilterAudioModule);
  virtual void resized();
  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void updateWidgetEnablement();
  juce_UseDebuggingNewOperator;
protected:
  TwoPoleFilterAudioModule *twoPoleFilterModuleToEdit;
  RComboBox *modeComboBox;
  rsModulatableSlider *frequencySlider, *gainSlider, *bandwidthSlider, *radiusSlider;
};

//===============================================================================================
// Delay Effects:

//-----------------------------------------------------------------------------------------------
// PingPongEcho:

class PingPongEchoAudioModule : public ModulatableAudioModule
{
public:
  PingPongEchoAudioModule(CriticalSection *lockTouse, rosic::PingPongEcho *echoToWrap = nullptr);
  virtual ~PingPongEchoAudioModule() { if(wrappedEchoIsOwned) delete wrappedEcho; }
  virtual void setBeatsPerMinute(double newBpm)
  {
    ScopedLock scopedLock(*lock);
    wrappedEcho->setTempoInBPM(newBpm);
  }
  virtual AudioModuleEditor* createEditor(int type) override;
  virtual void processStereoFrame(double *left, double *right) override 
  { 
    wrappedEcho->getSampleFrameStereo(left, right); 
  }
protected:
  virtual void createParameters();
  rosic::PingPongEcho *wrappedEcho;
  bool wrappedEchoIsOwned = false;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PingPongEchoAudioModule)
};

class PingPongEchoModuleEditor : public AudioModuleEditor
{
public:
  PingPongEchoModuleEditor(CriticalSection *newPlugInLock, PingPongEchoAudioModule* newPingPongEchoAudioModule);
  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void resized();
protected:
  rsModulatableSlider *delayTimeSlider, *dryWetSlider, *feedbackSlider, *panSlider, *highDampSlider,
    *lowDampSlider;
  RButton *pingPongButton, *tempoSyncButton, *trueStereoButton;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PingPongEchoModuleEditor)
};

//-----------------------------------------------------------------------------------------------
// Reverb:

class ReverbAudioModule : public ModulatableAudioModule
{
public:
  ReverbAudioModule(CriticalSection *newPlugInLock, rosic::rsReverb *newReverbToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::rsReverb *wrappedReverb;
};

class ReverbModuleEditor : public AudioModuleEditor
{
public:
  ReverbModuleEditor(CriticalSection *newPlugInLock, ReverbAudioModule* newReverbAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *dryWetSlider, *firstEchoSlider, *preDelaySlider, *decayTimeSlider,
    *highDecayScaleSlider, *lowDecayScaleSlider, *highFreqSlider, *lowFreqSlider;
  RButton *pinkButton, *stereoSwapButton;
};

//-----------------------------------------------------------------------------------------------
// SimpleDelay:

class SimpleDelayAudioModule : public ModulatableAudioModule
{
public:
  SimpleDelayAudioModule(CriticalSection *newPlugInLock, rosic::FractionalDelayLineStereo *newSimpleDelayToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::FractionalDelayLineStereo *wrappedSimpleDelay;
};

class SimpleDelayModuleEditor : public AudioModuleEditor
{
public:
  SimpleDelayModuleEditor(CriticalSection *newPlugInLock, SimpleDelayAudioModule* newSimpleDelayAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *delaySlider;
};


//===============================================================================================
// Modulation Effects:

//-----------------------------------------------------------------------------------------------
// ModulationEffect baseclass:

class ModulationEffectAudioModule : public ModulatableAudioModule //, public AudioFileManager //, public StateFileManager
{
  friend class ModulationEffectModuleEditor;
public:
  ModulationEffectAudioModule(CriticalSection *newPlugInLock, rosic::ModulationEffect *newModulationEffectToWrap);
  juce_UseDebuggingNewOperator;
protected:
  LowFrequencyOscillatorAudioModule *lfoModule;
  rosic::ModulationEffect *wrappedModulationEffect;
};

class ModulationEffectModuleEditor : public AudioModuleEditor
{
public:
  ModulationEffectModuleEditor(CriticalSection *newPlugInLock, ModulationEffectAudioModule* newModulationEffectAudioModule);
  virtual void setLfoPopUpEditorBounds(int x, int y, int w, int h)
  {
    ScopedLock scopedLock(*lock);
    lfoEditor->setPopUpEditorBounds(x, y, w, h);
  }
  virtual void resized();
  virtual void updateWidgetsAccordingToState();
  juce_UseDebuggingNewOperator;
protected:
  LowFrequencyOscillatorEditor *lfoEditor;
  ModulationEffectAudioModule *modulationEffectModuleToEdit;
  juce::Rectangle<int> lfoRect, effectRect;
  RTextField *lfoLabel, *effectLabel;
};

//-----------------------------------------------------------------------------------------------
// Flanger:

class FlangerAudioModule : public ModulationEffectAudioModule
{
  friend class FlangerModuleEditor;
public:
  FlangerAudioModule(CriticalSection *newPlugInLock, rosic::Flanger *newFlangerToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::Flanger *wrappedFlanger;
};

class FlangerModuleEditor : public ModulationEffectModuleEditor
{
public:
  FlangerModuleEditor(CriticalSection *newPlugInLock, FlangerAudioModule* newFlangerAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  FlangerAudioModule *flangerModuleToEdit;
  rsModulatableSlider *depthSlider, *dryWetSlider, *frequencySlider, *feedbackSlider;
  RButton *invertButton;
};

//-----------------------------------------------------------------------------------------------
// Phaser:

class PhaserAudioModule : public ModulationEffectAudioModule
{
  friend class PhaserModuleEditor;
public:
  PhaserAudioModule(CriticalSection *newPlugInLock, rosic::Phaser *newPhaserToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::Phaser *wrappedPhaser;
};

class PhaserModuleEditor : public ModulationEffectModuleEditor, public RComboBoxObserver
{
public:
  PhaserModuleEditor(CriticalSection *newPlugInLock, PhaserAudioModule* newPhaserAudioModule);
  virtual void resized();
  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void updateWidgetEnablement();
  juce_UseDebuggingNewOperator;
protected:
  PhaserAudioModule *phaserModuleToEdit;
  RTextField *filterLabel;
  rsModulatableSlider *depthSlider, *dryWetSlider, *frequencySlider, *qSlider, *feedbackSlider,
    *stagesSlider;
  RComboBox  *modeComboBox;
};

//-----------------------------------------------------------------------------------------------
// Tremolo:

class TremoloAudioModule : public ModulationEffectAudioModule
{
  friend class TremoloModuleEditor;
public:
  TremoloAudioModule(CriticalSection *newPlugInLock, rosic::Tremolo *newTremoloToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::Tremolo *wrappedTremolo;
};

class TremoloModuleEditor : public ModulationEffectModuleEditor
{
public:
  TremoloModuleEditor(CriticalSection *newPlugInLock, TremoloAudioModule* newTremoloAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  TremoloAudioModule *tremoloModuleToEdit;
  rsModulatableSlider *depthSlider;
};

//-----------------------------------------------------------------------------------------------
// Vibrato:

class VibratoAudioModule : public ModulationEffectAudioModule
{
  friend class VibratoModuleEditor;
public:
  VibratoAudioModule(CriticalSection *newPlugInLock, rosic::Vibrato *newVibratoToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::Vibrato *wrappedVibrato;
};

class VibratoModuleEditor : public ModulationEffectModuleEditor
{
public:
  VibratoModuleEditor(CriticalSection *newPlugInLock, VibratoAudioModule* newVibratoAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  VibratoAudioModule *vibratoModuleToEdit;
  rsModulatableSlider *depthSlider, *dryWetSlider;
};

//-----------------------------------------------------------------------------------------------
// WahWah:

class WahWahAudioModule : public ModulationEffectAudioModule
{
  friend class WahWahModuleEditor;
public:
  WahWahAudioModule(CriticalSection *newPlugInLock, rosic::WahWah *newWahWahToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::WahWah *wrappedWahWah;
};

class WahWahModuleEditor : public ModulationEffectModuleEditor, public RComboBoxObserver
{
public:
  WahWahModuleEditor(CriticalSection *newPlugInLock, WahWahAudioModule* newWahWahAudioModule);
  virtual void resized();
  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void updateWidgetEnablement();
  juce_UseDebuggingNewOperator;
protected:
  WahWahAudioModule *wahWahModuleToEdit;
  RTextField *filterLabel;
  rsModulatableSlider *depthSlider, *dryWetSlider, *frequencySlider, *gainSlider, *bandwidthSlider;
  RComboBox  *modeComboBox;
};

//===============================================================================================
// Spectral Effects:

//-----------------------------------------------------------------------------------------------
// FormantShifter:

class FormantShifterAudioModule : public ModulatableAudioModule
{
public:
  FormantShifterAudioModule(CriticalSection *newPlugInLock, rosic::FormantShifterStereo *newFormantShifterToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::FormantShifterStereo *wrappedFormantShifter;
};

class FormantShifterModuleEditor : public AudioModuleEditor
{
public:
  FormantShifterModuleEditor(CriticalSection *newPlugInLock, FormantShifterAudioModule* newFormantShifterAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *formantScaleSlider, *formantOffsetSlider, *dryWetSlider;
  // blockSizeComboBox/Slider
};


//===============================================================================================
// Other Effects:

//-----------------------------------------------------------------------------------------------
// Chorus:

class ChorusAudioModule : public ModulatableAudioModule
{
public:
  ChorusAudioModule(CriticalSection *newPlugInLock, rosic::Chorus *newChorusToWrap);
  virtual void parameterChanged(Parameter* parameterThatHasChanged); // remnant because of two-parametric callbacks
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::Chorus *wrappedChorus;
};

class ChorusModuleEditor : public AudioModuleEditor
{
public:
  ChorusModuleEditor(CriticalSection *newPlugInLock, ChorusAudioModule* newChorusAudioModule);
  virtual void resized();
  // todo: update per-voice slider enablement on switchin voices on/off
  juce_UseDebuggingNewOperator;
protected:
  juce::Rectangle<int> globalRect, voice1Rect, voice2Rect, voice3Rect, voice4Rect;
  RTextField *globalLabel;
  RButton    *voice1Button, *voice2Button, *voice3Button, *voice4Button;
  rsModulatableSlider *delaySlider, *cycleLengthSlider, *depthSlider, *globalFeedbackSlider, 
    *crossMixSlider, *feedback2Slider, *dryWetSlider, *stereoPhaseSlider;
  RSlider    *voice1DelaySlider, *voice1DepthSlider, *voice1AmpSlider,
    *voice2DelaySlider, *voice2DepthSlider, *voice2AmpSlider,
    *voice3DelaySlider, *voice3DepthSlider, *voice3AmpSlider,
    *voice4DelaySlider, *voice4DepthSlider, *voice4AmpSlider;
};

//-----------------------------------------------------------------------------------------------
// FrequencyShifter:

class FrequencyShifterAudioModule : public ModulatableAudioModule
{
public:
  FrequencyShifterAudioModule(CriticalSection *newPlugInLock, rosic::FrequencyShifterStereo *newFrequencyShifterToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::FrequencyShifterStereo *wrappedFrequencyShifter;
};

class FrequencyShifterModuleEditor : public AudioModuleEditor
{
public:
  FrequencyShifterModuleEditor(CriticalSection *newPlugInLock, FrequencyShifterAudioModule* newFrequencyShifterAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *shiftSlider, *feedbackSlider, *stereoOffsetSlider, *dryWetSlider, *midSideSlider;
};

//-----------------------------------------------------------------------------------------------
// PhaseStereoizer:

class PhaseStereoizerAudioModule : public ModulatableAudioModule
{
public:
  PhaseStereoizerAudioModule(CriticalSection *newPlugInLock, rosic::PhaseStereoizer *newPhaseStereoizerToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::PhaseStereoizer *wrappedPhaseStereoizer;
};

class PhaseStereoizerModuleEditor : public AudioModuleEditor
{
public:
  PhaseStereoizerModuleEditor(CriticalSection *newPlugInLock, PhaseStereoizerAudioModule* newPhaseStereoizerAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *phaseOffsetSlider, *dryWetRatioSlider, *sideLowpassSlider, *sideHighpassSlider,
    *midSideRatioSlider, *gainSlider;
  RButton *channelSwapButton;
};

//-----------------------------------------------------------------------------------------------
// RingModulator:

class RingModulatorAudioModule : public ModulatableAudioModule
{
public:
  RingModulatorAudioModule(CriticalSection *newPlugInLock, rosic::RingModulatorStereo *newRingModulatorToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::RingModulatorStereo *wrappedRingModulator;
};

class RingModulatorModuleEditor : public AudioModuleEditor
{
public:
  RingModulatorModuleEditor(CriticalSection *newPlugInLock, RingModulatorAudioModule* newRingModulatorAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *frequencySlider, *feedbackSlider, *stereoOffsetSlider, *dryWetSlider;
  RButton *antiAliasButton;
};

//-----------------------------------------------------------------------------------------------
// SingleSidebandModulator:

class SingleSidebandModulatorAudioModule : public ModulatableAudioModule
{
public:
  SingleSidebandModulatorAudioModule(CriticalSection *newPlugInLock,
    rosic::SingleSidebandModulatorStereo *newSingleSidebandModulatorToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::SingleSidebandModulatorStereo *wrappedSingleSidebandModulator;
};

class SingleSidebandModulatorModuleEditor : public AudioModuleEditor
{
public:
  SingleSidebandModulatorModuleEditor(CriticalSection *newPlugInLock,
    SingleSidebandModulatorAudioModule* newSingleSidebandModulatorAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *frequencySlider, *feedbackSlider, *stereoOffsetSlider, 
    *upperSidebandLevelSlider, *lowerSidebandLevelSlider, *dryWetSlider;
  RButton *antiAliasButton;
};

//-----------------------------------------------------------------------------------------------
// StereoPan:

class StereoPanAudioModule : public ModulatableAudioModule
{
  friend class StereoPanModuleEditor;
public:
  StereoPanAudioModule(CriticalSection *newPlugInLock, rosic::StereoPan *newStereoPanToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::StereoPan *wrappedStereoPan;
};

class StereoPanModuleEditor : public AudioModuleEditor //, public RSliderListener, public RComboBoxObserver
{
public:
  StereoPanModuleEditor(CriticalSection *newPlugInLock, StereoPanAudioModule* newStereoPanAudioModule);
  virtual ~StereoPanModuleEditor();
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  StereoPanAudioModule *stereoPanModuleToEdit;
  RTextField           *panLawLabel;
  RComboBox            *panLawComboBox;
  rsModulatableSlider    *panSlider, *gainSlider;
  StereoPanPlotEditor  *plot;
  double *xValues, *yValues;
  int    numValues;
};

//-----------------------------------------------------------------------------------------------
// StereoWidth:

class StereoWidthAudioModule : public ModulatableAudioModule
{
public:
  StereoWidthAudioModule(CriticalSection *newPlugInLock, rosic::StereoWidth *newStereoWidthToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::StereoWidth *wrappedStereoWidth;
};

class StereoWidthModuleEditor : public AudioModuleEditor
{
public:
  StereoWidthModuleEditor(CriticalSection *newPlugInLock, StereoWidthAudioModule* newStereoWidthAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *midSideRatioSlider, *gainSlider; // maybe include a correlation meter?
  RButton *monoButton;
};

//===============================================================================================
// Signal Generators:

//-----------------------------------------------------------------------------------------------
// SineOscillator:

class SineOscillatorAudioModule : public ModulatableAudioModule
{
public:
  SineOscillatorAudioModule(CriticalSection *newPlugInLock, rosic::SineOscillator *newSineOscillatorToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::SineOscillator *wrappedSineOscillator;
};

class SineOscillatorModuleEditor : public AudioModuleEditor
{
public:
  SineOscillatorModuleEditor(CriticalSection *newPlugInLock, SineOscillatorAudioModule* newSineOscillatorAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *frequencySlider, *levelSlider;
};

//-----------------------------------------------------------------------------------------------
// Noisifier:

class NoisifierAudioModule : public ModulatableAudioModule
{
public:
  NoisifierAudioModule(CriticalSection *newPlugInLock, rosic::Noisifier *newNoisifierToWrap);
  juce_UseDebuggingNewOperator;
protected:
  virtual void createStaticParameters();
  rosic::Noisifier *wrappedNoisifier;
};

class NoisifierModuleEditor : public AudioModuleEditor
{
public:
  NoisifierModuleEditor(CriticalSection *newPlugInLock, NoisifierAudioModule* newNoisifierAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  rsModulatableSlider *passLevelSlider, *noiseLevelSlider, *spectralSlopeSlider, *lowestFreqSlider,
    *highestFreqSlider;
};

#endif 
