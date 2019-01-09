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


  // frequency parameters:
  Parameter *maxNumModes;
  Parameter *ratioProfileTopLeft, *ratioProfileTopRight, *ratioProfileBottomLeft, 
    *ratioProfileBottomRight;
  Parameter *freqRatiosX, *freqRatiosY; // should also have key/vel scaling

  // have a general "Tune" parameter (maybe also with ByKey, ByVel options) because the percieved 
  // frequency may not coincide with the lowest frequency - and when ratios are key/vel dependent,
  // that may also change as function of key/vel

  // amp and envelope parameters:
  Parameter *amp,    *ampByRatio,    *ampByKey,    *ampByVel;
  Parameter *attack, *attackByRatio, *attackByKey, *attackByVel;
  Parameter *decay,  *decayByRatio,  *decayByKey,  *decayByVel;


  // phase parameters:





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


  RSlider   *sldMaxNumModes;
  RComboBox *boxTopLeftRatios, *boxTopRightRatios, *boxBottomLeftRatios, *boxBottomRightRatios;
  rsVectorPad *xyPadRatios;



  //rsSlider *leftSlider, *rightSlider, *smoothSlider, *testSlider;

  //rsButton *invertButton;  
  // use ModulatabelSlider, etc later

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModalSynthEditor)
};