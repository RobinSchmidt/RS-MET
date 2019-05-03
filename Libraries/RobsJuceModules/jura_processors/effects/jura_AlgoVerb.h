#ifndef jura_AlgoVerbAudioModule_h
#define jura_AlgoVerbAudioModule_h


/** This class wraps rosic::AlgoVerb into a rosof::AudioModule to facilitate its use as plugIn. */

class AlgoVerbAudioModule : public AudioModule
{

  friend class AlgoVerbModuleEditor;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  AlgoVerbAudioModule(CriticalSection *newPlugInLock, rosic::AlgoVerb *algoVerbToWrap);

  AlgoVerbAudioModule(CriticalSection *newPlugInLock);

  virtual ~AlgoVerbAudioModule();

  AudioModuleEditor* createEditor(int type) override;
  // still commneted out/preliminary because the AudioModule doesn't have all the parameters that
  // the custom editor expects - fix that

  //-----------------------------------------------------------------------------------------------
  // automation and state management:

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate) override
  {
    wrappedAlgoVerb->setSampleRate((float)newSampleRate);
  }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedAlgoVerb->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    for(int n = 0; n < numSamples; n++)
      wrappedAlgoVerb->getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
  }

  //---------------------------------------------------------------------------------------------
  // event processing:

  virtual void reset() override
  {
    wrappedAlgoVerb->reset();
  }

protected:

  void initializeAutomatableParameters();

  rosic::AlgoVerb *wrappedAlgoVerb;
  bool wrappedAlgoVerbIsOwned = false;


  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class AlgoVerbModuleEditor : public AudioModuleEditor, public RComboBoxObserver
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  AlgoVerbModuleEditor(CriticalSection *newPlugInLock, 
    AlgoVerbAudioModule* newAlgoVerbAudioModule);

  //---------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();



protected:

  /** Updates the plots on the GUI. */
  virtual void updatePlots();

  AlgoVerbAudioModule *algoVerbModuleToEdit;

  juce::Rectangle<int> globalRect, earlyRect, lateRect;

  RTextField *globalLabel, *earlyLabel, *lateLabel;

  RSlider   *dryWetSlider, *lateLevelSlider, *referenceDelayTimeSlider, *densitySlider, 
    *diffusionSlider, *latePreDelaySlider, *decayTimeSlider, *lowDecayScaleSlider, 
    *lowCrossFreqSlider, *highDecayScaleSlider, *highCrossFreqSlider;

  RComboBox *feedbackMatrixComboBox, *injectionVectorComboBox, *outputVectorComboBox;

  RButton   *pingButton, *earlyPingButton, *latePingButton, *allpassModeButton, 
    *wetPinkingButton;

  /*
  RSlider *coarseSlider, *fineSlider, *grainLengthInMillisecondsSlider, *grainLengthInCyclesSlider,
  *grainLengthInBeatsSlider, *feedbackSlider, *dryWetSlider;
  RComboBox *grainLengthUnitComboBox;
  RButton *invertButton, *reverseButton, *antiAliasButton; // *formantPreserveButton, *monoButton;
  */

  WaveformDisplay *impulseResponsePlot;


  juce_UseDebuggingNewOperator;
};


#endif 
