#ifndef jura_Snowflake_h
#define jura_Snowflake_h

class JUCE_API Snowflake : public jura::AudioModuleWithMidiIn, public ChangeBroadcaster
{

public:

  Snowflake(CriticalSection *lockToUse);

  virtual void createParameters();

  // overriden from AudioModule baseclass:
  AudioModuleEditor* createEditor() override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void reset() override;
  virtual void noteOn(int noteNumber, int velocity) override;

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;
  virtual XmlElement convertXmlStateIfNecessary(const XmlElement& xmlState) override;

  // override set/getXml to store strings for rules and seed

  void setResetRatio1( double newRatio)  { core.setResetRatio( 0, newRatio);  }
  void setResetOffset1(double newOffset) { core.setResetOffset(0, newOffset); }
  void setResetRatio2( double newRatio)  { core.setResetRatio( 1, newRatio);  }
  void setResetOffset2(double newOffset) { core.setResetOffset(1, newOffset); }

  void setReverseRatio1( double newRatio)  { core.setReverseRatio( 0, newRatio);  }
  void setReverseOffset1(double newOffset) { core.setReverseOffset(0, newOffset); }


  /** Sets the axiom (i.e. seed") for the L-system. */
  void setAxiom(const juce::String& newAxiom);

  /** Sets the string containing the L-system rules and returns whether or not this was sucessful. 
  It may fail, if the string is malformed in which case and empty set of rules will be used by 
  default. */
  bool setRules(const juce::String& newRules);

  /** Returns true, if the given string with L-system rules is properly formed. */
  bool validateRuleString(const juce::String& newRules);

  juce::String getAxiom() const { return axiom; }
  juce::String getRules() const { return rules; }
  int getNumTurtleLines() { return core.getNumTurtleLines(); }

protected:

  rosic::Snowflake core;

  juce::String axiom, rules;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Snowflake)
};

//=================================================================================================

class JUCE_API SnowflakePlot2D
{

};

class JUCE_API SnowflakePlot1D
{

};

//=================================================================================================

/** Editor for the Snowflake AudioModule.
todo:
-add plots for:
 -initiator curve (axiom)
 -generator curves (rules)
 -resulting curve in 2D
 -resulting x,y signals
OR:
 -add a ScopeDisplay that shows the axiom and rules when no note is played and (optionally) shows 
  the realtime output when a note is played. scope widgets can be overlaid over the display when
  the mouse is over the display
  -subdivide the scope screen into a grid of squares, 2x2, 3x3, 4x4 - whatever is necessary to show
   the axiom plus all rules - use one little square for each rule, overlay the plot with the string
   of the rule, like A=AF+FA-BF+ ..or whatever the rule is (use semi-transparent text)
*/

class JUCE_API SnowflakeEditor : public AudioModuleEditor, public RTextEditorListener, 
  public RSliderListener
{

public:

  SnowflakeEditor(jura::Snowflake *snowFlakeToEdit);

  void createWidgets();

  virtual void resized() override;
  virtual void rTextEditorTextChanged(RTextEditor& editor) override;
  virtual void rSliderValueChanged(RSlider* slider) override;
  virtual void updateWidgetsAccordingToState() override;

  void updateNumTurtleLinesField();

  //virtual void changeListenerCallback(ChangeBroadcaster* source) override;

protected:

  jura::Snowflake *snowflakeModule;

  RTextField  *axiomLabel,  *rulesLabel;
  RTextEditor *axiomEditor, *rulesEditor;
  RTextField  *numLinesLabel;  

  RSlider *sliderIterations, *sliderAngle, *sliderPhase, *sliderSkew, *sliderAmplitude, 
    *sliderRotation;

  RSlider *sliderResetRatio1,   *sliderResetOffset1,   *sliderResetRatio2,   *sliderResetOffset2;
  RSlider *sliderReverseRatio1, *sliderReverseOffset1, *sliderReverseRatio2, *sliderReverseOffset2;

  RSlider *sliderFreqScaler; // maybe use this just for development and remove later
  RButton *buttonAntiAlias, *buttonUseTable;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SnowflakeEditor)
};


#endif