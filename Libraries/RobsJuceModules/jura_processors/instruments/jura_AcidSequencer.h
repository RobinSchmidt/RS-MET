#ifndef jura_AcidSequencer_h
#define jura_AcidSequencer_h

//#include "../../../rosic/sequencing/rosic_AcidSequencer.h"
//using namespace rosic;
//
//#include "../rosof_AudioModule.h"


/** This class wraps rosic::AcidSequencer into a rosof::AudioModule. */

class AcidSequencerAudioModule : public AudioModule
{

  friend class AcidSequencerModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  AcidSequencerAudioModule(CriticalSection *newPlugInLock, rosic::AcidSequencer *acidSequencerToWrap);

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  /*
  virtual void setSampleRate(double newSampleRate)
  {
    if( wrappedAcidSequencer != NULL )
      wrappedAcidSequencer->setSampleRate(newSampleRate);
  }
  */

  //---------------------------------------------------------------------------------------------
  // automation and state management:

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;


  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    // empty, only to satisfy compiler
  }

  //---------------------------------------------------------------------------------------------
  // others:

protected:

  void createParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::AcidSequencer *wrappedAcidSequencer;

  juce_UseDebuggingNewOperator;
};


//=================================================================================================

class AcidPatternEditor : public Component // public RWidget ...maybe needs to inherit from ColouSchemeComponent?
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  AcidPatternEditor(rosic::AcidSequencer *sequencerToEdit);

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Sets the pattern that is to be edited. */
  void setPatternToEdit(rosic::AcidPattern *newPatternToEdit);

  //virtual void setColourScheme(const WidgetColourScheme& newColourScheme) override;

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the step-number (starting at 0) that corresponds to the given x-coordinate. 
  Returns -1, if no step is at the given x-coordinate (i.e. the x-coordinate is left from the 
  actual sequencer). */
  int getStepAt(float xCoordinate);

  /** Returns the note-number (0...1) that corresponds to the given y-coordinate. 
  Returns -1, if no note is at the given y-coordinate. (i.e. the y-coordinate is above the actual
  sequencer). */
  int getNoteAt(float yCoordinate);

  /** Returns the key-number on the virtual keyboard at the given position if any, -1 if none. */
  int getVirtualKeyAt(float x, float y);

  /** Returns true if given y-coordinate is in the gate-row, false otherwise). */
  bool isInGateRow(float yCoordinate);

  /** Returns true if given y-coordinate is in the accent-row, false otherwise). */
  bool isInAccentRow(float yCoordinate);

  /** Returns true if given y-coordinate is in the slide-row, false otherwise). */
  bool isInSlideRow(float yCoordinate);

  /** Returns true if given y-coordinate is in the octave-row, false otherwise). */
  bool isInOctaveRow(float yCoordinate);

  /** Returns true if given x-coordinate is in the keyboard column, false otherwise). */
  bool isInKeyboardColumn(float xCoordinate);

  /** Returns the height of one row in the top lane (for accent, slide, etc.) */
  int getTopLaneHeight() const { return (int) topLaneHeight; } // why is this a fdloat?

  //---------------------------------------------------------------------------------------------
  // callbacks:

  virtual void mouseDown(const MouseEvent &e) override;
  virtual void paint(Graphics &g) override;



protected:

  rosic::AcidPattern   *patternToEdit;
  rosic::AcidSequencer *sequencerToEdit;

  Colour whiteKeyColour, blackKeyColour, backgroundColourWhiteKey, backgroundColourBlackKey,
    handleColor, textColour, lineColour;

  float rowHeight, columnWidth, keyLength, topLaneHeight;

  juce_UseDebuggingNewOperator;
};


//===============================================================================================
// class AcidSequencerModuleEditor:

class AcidSequencerModuleEditor : public AudioModuleEditor, public RSliderListener
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  AcidSequencerModuleEditor(CriticalSection *newPlugInLock, AcidSequencerAudioModule* newAcidSequencerAudioModule);

  //---------------------------------------------------------------------------------------------
  // setup:

  //---------------------------------------------------------------------------------------------
  // callbacks:

  /** Overrides rButtonClicked() for triggering various actions and update the display 
  accordingly. */
  virtual void rButtonClicked(RButton *buttonThatWasClicked);

  /** Overrides rSliderValueChanged to update the sequencer display when the steLength slider 
  changes. */  
  virtual void rSliderValueChanged(RSlider *rSliderThatHasChanged);

  /** Overrides resized(). */    
  virtual void paint(Graphics &g);

  /** Overrides resized(). */    
  virtual void resized();

  /** Overrides the method inherited from AudioModuleEditor. */
  virtual void updateWidgetsAccordingToState();

protected:

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  AcidSequencerAudioModule *acidSequencerModuleToEdit;
  AcidPatternEditor        *patternEditor;

  RTextField *modeLabel, *shiftLabel;  // shiftLabel not use anymore
  RComboBox  *modeBox;
  RSlider    *stepLengthSlider;
  RButton    *shiftLeftButton, *shiftRightButton;
  RButton    *lockButton;  // prevents the sequencer data from being changed and preset-change

  // new:
  RButton    *shiftAccentsLeftButton, *shiftAccentsRightButton;
  RButton    *shiftSlidesLeftButton,  *shiftSlidesRightButton;
  RButton    *shiftNotesLeftButton,   *shiftNotesRightButton;
  RButton    *shiftOctavesLeftButton, *shiftOctavesRightButton;

  RButton    *reverseAllButton;
  RButton    *reverseAccentsButton;
  RButton    *reverseSlidesButton;
  RButton    *reverseNotesButton;
  RButton    *reverseOctavesButton;
  RButton    *swapAccentsSlidesButton;
  RButton    *xorAccentsSlidesButton;
  RButton    *xorSlidesAccentsButton;

  RButton    *invertAccentsButton;
  RButton    *invertSlidesButton;
  RButton    *invertOctavesButton;


  juce_UseDebuggingNewOperator;
};

#endif
