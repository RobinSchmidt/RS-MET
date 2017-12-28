#ifndef jura_DebugModule_h
#define jura_DebugModule_h

/** A module for framework debugging puposes */

class JUCE_API DebugAudioModule : public jura::AudioModuleWithMidiIn
{

public:

  DebugAudioModule(CriticalSection *lockToUse);
    
  virtual void createParameters();

  // overriden from AudioModule baseclass:
  AudioModuleEditor* createEditor() override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;
  virtual void setMidiController(int controllerNumber, float controllerValue) override;

  // callback target functions:
  void setLeftValue(double newValue);
  void setRightValue(double newValue);
  void setSmoothingTime(double newTime);

protected:


  static const int numValues = 2;
  double values[numValues] = { 0, 0 };

  MetaControlledParameter *leftParam, *rightParam, *smoothParam;

  EqualizerAudioModule* eqModule == nullptr;
  // We use an equalizer to see if dealing with child-modules and works well. Also, the eq
  // has a dynamic number of parameters, so we can check that stuff with eq, too.


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(DebugAudioModule)
};

//=================================================================================================

/** Editor for the debugging AudioModule. */

class JUCE_API DebugModuleEditor : public AudioModuleEditor
{

public:

  DebugModuleEditor(jura::DebugAudioModule *newDebugModuleToEdit);
  virtual void resized() override;

protected:

  void createWidgets();

  DebugAudioModule *debugModule;

  rsVectorPad *xyPad;

  AutomatableSlider *leftSlider, *rightSlider, *smoothSlider;
  //AutomatableComboBox *modeComboBox;
  //AutomatableButton *invertButton;  
   // use ModulatabelSlider, etc later

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(DebugModuleEditor)
};


#endif 