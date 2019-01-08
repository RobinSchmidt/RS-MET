#pragma once

/** A module for framework debugging puposes */

class JUCE_API ModalSynthAudioModule : public AudioModuleWithMidiIn
{

public:

  ModalSynthAudioModule(CriticalSection *lockToUse);

  virtual void createParameters();

  // overriden from AudioModule baseclass:
  //AudioModuleEditor* createEditor(int type) override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;
  virtual void noteOn( int key, int vel) override { core.noteOn(key, vel); }
  virtual void noteOff(int key)          override { core.noteOn(key, 0); }


protected:

  void populateFreqRatioProfileParam(Parameter* p);

  Parameter *ratioProfileTopLeft, *ratioProfileTopRight, *ratioProfileBottomLeft, 
    *ratioProfileBottomRight;

  Parameter *freqRatiosX, *freqRatiosY;

  Parameter *maxNumModes;

  ParameterWithKeyVelScaling *attack, *decay;

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

  rsVectorPad *xyPad;

  //rsSlider *leftSlider, *rightSlider, *smoothSlider, *testSlider;
  //rsComboBox *modeComboBox;
  //rsButton *invertButton;  
  // use ModulatabelSlider, etc later

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModalSynthEditor)
};