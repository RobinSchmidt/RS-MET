#ifndef jura_DspWorkbench_h
#define jura_DspWorkbench_h

class DspWorkbenchAudioModule : public AudioModule
{

  friend class DspWorkbenchModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  DspWorkbenchAudioModule(CriticalSection *newPlugInLock, 
    rosic::DspWorkbench *dspWorkbenchToWrap = nullptr);

  virtual ~DspWorkbenchAudioModule();

  AudioModuleEditor* createEditor(int type) override;



  //---------------------------------------------------------------------------------------------
  // automation and state management:

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  /** Sets up the algorithm to be used as a juse::juce::String object. */
  virtual bool setAlgorithmString(const juce::String& newAlgorithmString);

  /** Returns the algorithm as juse::juce::String object */
  virtual juce::String getAlgorithmString();

  virtual void setSampleRate(double newSampleRate) override
  {
    wrappedDspWorkbench->setSampleRate(newSampleRate);
  }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inL, double* inR, double* outL, double* outR)
  {
    wrappedDspWorkbench->getSampleFrameStereo(inL, inR, outL, outR);
  }

protected:

  void initializeAutomatableParameters();

  /** Strips off the comments from a program text written in SPSL (signal processing scripting
  language). */
  juce::String stripOffComments(const juce::String& inputText);

  juce::String algorithmWithComments, algorithmWithoutComments;

  rosic::DspWorkbench *wrappedDspWorkbench;
  bool wrappedDspWorkbenchIsOwned = false;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class DspWorkbenchModuleEditor : public AudioModuleEditor, public RTextEntryFieldObserver, 
  public RSliderListener, public RTextEditorListener
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  DspWorkbenchModuleEditor(CriticalSection *newPlugInLock, DspWorkbenchAudioModule* newDspWorkbenchAudioModule);

  /** Destructor. */
  //virtual ~DspWorkbenchModuleEditor();

  //---------------------------------------------------------------------------------------------
  // setup:

  //---------------------------------------------------------------------------------------------
  // callbacks:

  /** Implements the purely virtual rButtonClicked()-method of the ButtonListener base-class. */
  virtual void rButtonClicked(RButton *buttonThatWasClicked) override;

  /** Implements the purely virtual rButtonClicked()-method of the ComboBoxListener base-class. */
  virtual void comboBoxChanged(ComboBox *comboBoxThatHasChanged);

  /** Implements the purely virtual labelTextChanged()-method of the LabelListener base-class. */
  //virtual void rLabelTextChanged(RLabel *rLabelThatHasChanged);

  /** Implements the purely virtual rSliderValueChanged-method of the RSliderListener 
  base-class. */
  virtual void rSliderValueChanged(RSlider *sliderThatHasChanged) override;

  /** Implements the purely virtual rTextFieldChanged-method of the RTextFieldListener 
  base-class. */
  virtual void rTextFieldChanged(RTextField *textFieldThatHasChanged);

  /** Overrides paint(). */   
  virtual void paint(Graphics &g) override;

  /** Overrides resized(). */    
  virtual void resized() override;

  //virtual void textEditorTextChanged(TextEditor &editor);
  //virtual void textEditorReturnKeyPressed(TextEditor &editor);
  //virtual void textEditorEscapeKeyPressed(TextEditor &editor);
  //virtual void textEditorFocusLost(TextEditor &editor);

  virtual void rTextEditorTextChanged(RTextEditor& editor) override;
  virtual void rTextEditorReturnKeyPressed(RTextEditor& editor) override;
  virtual void rTextEditorEscapeKeyPressed(RTextEditor& editor) override;
  virtual void rTextEditorFocusLost(RTextEditor& editor) override;


  virtual void textChanged(RTextEntryField *rTextEntryFieldThatHasChanged) override {}

  //---------------------------------------------------------------------------------------------
  // others:

  // some widgets need to be made public in order to let them be set from outside in repsonse to
  // automation and MIDI-events:
  //RButton* noteOnButton;
  //RButton* noteOffButton;



protected:

  /** Overrides the method inherited from RPolyphonicInstrumentEditor. */
  virtual void updateWidgetsAccordingToState() override;

  /** Shows or hides certain modulation-target specific wdgets according to the currently 
  selected modulation-target. */
  virtual void showOrHideTargetSpecificWidgets();

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  DspWorkbenchAudioModule *dspWorkbenchAudioModule;


  // some rectangles to subdivide the GUI:
  juce::Rectangle<int> topLeftRectangle, topRightRectangle, bottomLeftRectangle, 
    bottomRightRectangle;

  // the main code editor:
  RTextEditor *codeEditor;

  // the embedded sub-editors:
  BreakpointModulatorEditor *breakpointModulatorEditor;
  //SampleModulatorEditor          *sampleModulatorEditor;
  MultiModeFilterModuleEditor     *filterEditor;

  // some global widgets:
  RSlider *oversamplingSlider;

  // the widgets for the user parameters:
  //RLabel*     parameterLabels[rosic::DspWorkbench::numParameters];
  RTextEntryField* parameterLabels[rosic::DspWorkbench::numParameters];
  RSlider*         parameterSliders[rosic::DspWorkbench::numParameters];
  RButton*         parameterExpButtons[rosic::DspWorkbench::numParameters];
  RTextEntryField* parameterMinFields[rosic::DspWorkbench::numParameters];
  RTextEntryField* parameterMaxFields[rosic::DspWorkbench::numParameters];
  RTextEntryField* parameterDefaultFields[rosic::DspWorkbench::numParameters];
  RTextEntryField* parameterNameFields[rosic::DspWorkbench::numParameters];

  /*
  // the labels for the user parameters:
  RLabel *par01Label, par02Label;

  // the sliders for the user parameters:
  RSlider *par01Slider, *par02Slider;

  // the buttons for the exponential scaling of the user parameters:
  RButton *par01ExpButton, par02ExpButton;

  // the textfields for the minimum, maximum and default values and names of the user parameters:
  RTextField *par01MinField, *par01MaxField, *par01DefaultField, *par01NameField,
  *par02MinField, *par02MaxField, *par02DefaultField, *par02NameField;
  */

  juce_UseDebuggingNewOperator;
};



#endif 