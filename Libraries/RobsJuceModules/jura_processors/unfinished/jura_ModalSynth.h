#pragma once

/** A module for framework debugging puposes */

class JUCE_API ModalSynthAudioModule : public AudioModuleWithMidiIn
{

public:

  ModalSynthAudioModule(CriticalSection *lockToUse);

  virtual void createParameters();

  // overriden from AudioModule baseclass:
  AudioModuleEditor* createEditor(int type) override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;
  virtual void noteOn( int key, int vel) override { core.noteOn(key, vel); }
  virtual void noteOff(int key)          override { core.noteOn(key, 0); }


protected:

  void populateFreqRatioProfileParam(Parameter* p);


  // global parameters
  Parameter *level, *levelByKey, *levelByVel;
  Parameter *detune, *detuneByKey, *detuneByVel;

  // frequency parameters:
  Parameter *maxNumModes;
  Parameter *ratioProfileTopLeft, *ratioProfileTopRight, *ratioProfileBottomLeft, 
    *ratioProfileBottomRight;
  Parameter *freqRatiosX, *freqRatiosY; // should also have key/vel scaling

  // amplitude-spectrum and envelope parameters:
  Parameter *ampSlope, *ampSlopeByKey, *ampSlopeByVel;
  Parameter *attack, *attackByRatio, *attackByKey, *attackByVel;
  Parameter *decay,  *decayByRatio,  *decayByKey,  *decayByVel;
  Parameter *attackScale, *decayScale, *freqDelta, *phaseDelta;
  Parameter *blend,  *blendByKey,  *blendByVel;

  // actually, we not have to store the pointers in member variables if ywe don't need to refer
  // to them - LowestMode/HighestMode

  // phase-spectrum parameters:





  rosic::rsModalSynth core;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModalSynthAudioModule)
};

//=================================================================================================

/** Editor for the modal synthesizer. */

class JUCE_API ModalSynthEditor : public AudioModuleEditor
{

public:

  ModalSynthEditor(jura::ModalSynthAudioModule *newModalSynthModuleToEdit);
  virtual void resized() override;

protected:

  void createWidgets();

  ModalSynthAudioModule* modalModule;

  // global parameters:
  RSlider *sldLevel, *sldLevelByKey, *sldLevelByVel;
  RSlider *sldDetune, *sldDetuneByKey, *sldDetuneByVel;

  // freq-ratio parameters:
  RComboBox *boxTopLeftRatios, *boxTopRightRatios, *boxBottomLeftRatios, *boxBottomRightRatios;
  rsVectorPad *xyPadRatios;
  RSlider *sldLowestMode, *sldHighestMode;
  RSlider /* *sldMaxNumModes, */ *sldRatiosX, *sldRatiosY;

  // magnitude spectrum parameters:
  RSlider *sldAmpSlope, *sldAmpSlopeByKey, *sldAmpSlopeByVel;
  RSlider *sldAttack, *sldAttackByRatio, *sldAttackByKey, *sldAttackByVel;
  RSlider *sldDecay,  *sldDecayByRatio,  *sldDecayByKey,  *sldDecayByVel;
  RSlider *sldAttackScale, *sldDecayScale, *sldFreqDelta, *sldPhaseDelta;
  RSlider *sldBlend,  *sldBlendByKey,  *sldBlendByVel;

  //rsSlider *leftSlider, *rightSlider, *smoothSlider, *testSlider;

  //rsButton *invertButton;  
  // use ModulatabelSlider, etc later

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModalSynthEditor)
};