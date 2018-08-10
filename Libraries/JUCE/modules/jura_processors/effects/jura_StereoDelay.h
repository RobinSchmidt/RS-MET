#ifndef jura_StereoDelay_h
#define jura_StereoDelay_h

class StereoDelayAudioModule : public AudioModule
{

  friend class StereoDelayModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  StereoDelayAudioModule(CriticalSection *newPlugInLock, 
    rosic::StereoDelay *stereoDelayToWrap = nullptr);

  virtual ~StereoDelayAudioModule();

  AudioModuleEditor* createEditor(int type) override;

  //---------------------------------------------------------------------------------------------
  // automation and state management:

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate) override
  {
    wrappedStereoDelay->setSampleRate(newSampleRate);
  }

  virtual void setBeatsPerMinute(double newBpm) override
  {
    wrappedStereoDelay->setBeatsPerMinute(newBpm);
  }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    for(int n = 0; n < numSamples; n++)
      wrappedStereoDelay->getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
  }

  virtual void processStereoFrame(double *left, double *right) override
  {
    wrappedStereoDelay->getSampleFrameStereo(left, right);
  }

  //---------------------------------------------------------------------------------------------
  // event processing:

  virtual void reset() override
  {
    wrappedStereoDelay->reset();
  }


protected:

  void initializeAutomatableParameters();

  rosic::StereoDelay *wrappedStereoDelay;
  bool wrappedStereoDelayIsOwned = false;


  juce_UseDebuggingNewOperator;
};

//=================================================================================================


class StereoDelayModuleEditor : public AudioModuleEditor, public RComboBoxObserver
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  StereoDelayModuleEditor(CriticalSection *newPlugInLock, 
    StereoDelayAudioModule* newStereoDelayAudioModule);
  /**< Constructor. */

  //virtual ~StereoDelayModuleEditor();
  /**< Destructor. */

  //---------------------------------------------------------------------------------------------
  // callbacks:

  //virtual void rButtonClicked(RButton *buttonThatWasClicked);
  /**< Implements the purely virtual rButtonClicked()-method of the ButtonListener base-class. */

  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);

  //---------------------------------------------------------------------------------------------
  // others:

  /** Overrides paint(). */   
  virtual void paint(Graphics &g);

  /** Overrides resized(). */    
  virtual void resized();


protected:

  virtual void updateWidgetsAccordingToState();
  /**< Overrides the method inherited from RPolyphonicInstrumentEditor. */

  StereoDelayAudioModule *stereoDelayAudioModule;
  /**< This is the actual plugin engine which does all the dsp and automation handling. */

  juce::Rectangle<int> leftRectangle, rightRectangle;
  /**< Some rectangles to subdivide the GUI. */

  // global widgets:
  RSlider *dryWetSlider, *cutoffScaleSlider;

  // left delayline widgets:
  RTextField *delayLineLabelL, *inputLabelL, *diffusorLabelL, *filterLabelL, 
    *feedbackLabelL, *outputLabelL;
  RTextField *delayLineLabelR, *inputLabelR, *diffusorLabelR, *filterLabelR, 
    *feedbackLabelR, *outputLabelR;

  RSyncIntervalComboBox *delayComboBoxL;
  RSyncIntervalComboBox *delayComboBoxR;

  RSlider *delayScaleSliderL; 
  RSlider *inputSliderL2L, *inputSliderR2L;
  RSlider *diffusorTimeSliderL, *diffusorAmountSliderL;
  RSlider *lowpassSliderL, *highpassSliderL; 
  RSlider *feedbackSliderL, *crossFeedbackSliderL;
  RSlider *outputSliderL2L, *outputSliderL2R, *outDelaySliderL;
  RSlider *delayScaleSliderR; 
  RSlider *inputSliderL2R, *inputSliderR2R;
  RSlider *diffusorTimeSliderR, *diffusorAmountSliderR;
  RSlider *lowpassSliderR, *highpassSliderR; 
  RSlider *feedbackSliderR, *crossFeedbackSliderR;
  RSlider *outputSliderR2L, *outputSliderR2R, *outDelaySliderR;

  /*
  RSlider *coarseSlider, *fineSlider, *grainLengthInMillisecondsSlider, 
  *grainLengthInCyclesSlider, *grainLengthInBeatsSlider, *feedbackSlider, *dryWetSlider;
  RLabel    *grainLengthUnitLabel;
  RComboBox *grainLengthUnitComboBox;
  RButton *invertButton, *reverseButton, *formantPreserveButton, *antiAliasButton;
  */

  juce_UseDebuggingNewOperator;
};


#endif 
