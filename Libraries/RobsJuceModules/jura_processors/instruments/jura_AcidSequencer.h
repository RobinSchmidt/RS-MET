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

class AcidPatternEditor : public jura::ColourSchemeComponent // old: : public Component 
{

  // maybe we should derive from jura::RWidget instead. That may be more economic

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Constructor. */
  AcidPatternEditor(rosic::AcidSequencer *sequencerToEdit);

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the pattern that is to be edited. */
  void setPatternToEdit(rosic::AcidPattern *newPatternToEdit);

  //virtual void setColourScheme(const WidgetColourScheme& newColourScheme) override;

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

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
  int getTopLaneHeight() const { return (int) topLaneHeight; } // why is this a float?



  juce::Colour getColorBackground()
  {
    return widgetColourScheme.background;
    //return juce::Colours::black; 
  }
  // preliminary - maybe use the average/flattened background color fo the plots

  /** Returns the color to be used to draw the white keys in the little keyboard at the left. */
  //juce::Colour getColorWhiteKeys() { return whiteKeyColour; }
  juce::Colour getColorWhiteKeys() { return widgetColourScheme.handle.brighter(0.75f); }
  // maybe use a whitened widget handle color

  /** Returns the color to be used to draw the black keys in the little keyboard at the left. */
  //juce::Colour getColorBlackKeys() { return blackKeyColour; }
  juce::Colour getColorBlackKeys() 
  { 
    //return Colours::red;  // debug
    return widgetColourScheme.handle.darker(0.75f); 
  }
  // maybe use a darkened widget handle color

  /** Returns the color to be used to draw the lanes for the white keys. */
  //juce::Colour getColorWhiteLanes() { return backgroundColourWhiteKey; }
  juce::Colour getColorWhiteLanes() { return widgetColourScheme.background; }
  // use plotColourScheme.getFlatBackgroundColour()


  /** Returns the color to be used to draw the lanes for the black keys. */
  //juce::Colour getColorBlackLanes() { return backgroundColourBlackKey; }
  juce::Colour getColorBlackLanes() { return widgetColourScheme.background.brighter(0.15f); }
  // use same color as for white key lanes with some additional layover gray using brighter works
  // only for dark-on-bright, i think

  //juce::Colour getColorHandles() { return handleColor; }
  //juce::Colour getColorHandles() { return widgetColourScheme.handle; }
  juce::Colour getColorHandles() { return widgetColourScheme.weakHighlight; }
  // use widgetColourScheme.handle

  //juce::Colour getColorText() { return editorColourScheme.text; }
  juce::Colour getColorText() { return widgetColourScheme.text; }
  // use widgetColourScheme.text

  juce::Colour getColorLines() { return widgetColourScheme.outline; }

  //juce::Colour getColorLines() { return plotColourScheme.coarseGrid; }
  // use wdigetColourScheme.coarseGrid
  // widgetColourScheme.outline


  //---------------------------------------------------------------------------------------------
  // \name Callback overrides :

  virtual void mouseDown(const MouseEvent &e) override;
  virtual void paint(Graphics &g) override;



protected:

  rosic::AcidPattern   *patternToEdit;
  rosic::AcidSequencer *sequencerToEdit;


  float rowHeight, columnWidth, keyLength, topLaneHeight;
  // Why are these floats? Can't we make them integers?
  // rename keyLength to keySize..hmm..no

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

  RButton    *lockButton;  // prevents the sequencer data from being changed and preset-change

  // Sequence manipulators:
  RButton    *shiftLeftButton, *shiftRightButton;
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

  // use RClickButton instead


  juce_UseDebuggingNewOperator;
};

#endif
