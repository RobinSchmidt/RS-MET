#ifndef jura_Snowflake_h
#define jura_Snowflake_h

class JUCE_API Snowflake : public jura::AudioModuleWithMidiIn
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

  // override set/getXml to store strings for rules and seed

  /** Sets the seed (i.e. "axiom") for the L-system. */
  void setSeed(const juce::String& newSeed);

  /** Sets the string containing the L-system rules and returns whether or not this was sucessful. 
  It may fail, if the string is malformed in which case and empty set of rules will be used by 
  default. */
  bool setRules(const juce::String& newRules);

  /** Returns true, if the given string with L-system rules is properly formed. */
  bool validateRuleString(const juce::String& newRules);

protected:

  rosic::Snowflake core;

  juce::String rules, seed;

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
*/

class JUCE_API SnowflakeEditor : public AudioModuleEditor
{

public:

  SnowflakeEditor(jura::Snowflake *snowFlakeToEdit);

  void createWidgets();

  virtual void resized() override;

protected:

  jura::Snowflake *snowflakeModule;

  RSlider *sliderIterations, *sliderAngle;

  RTextField  *seedLabel,  *rulesLabel;
  RTextEditor *seedEditor, *rulesEditor;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SnowflakeEditor)
};


#endif