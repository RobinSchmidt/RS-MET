#ifndef jura_Snowflake_h
#define jura_Snowflake_h

class JUCE_API Snowflake : public jura::AudioModuleWithMidiIn
{

public:

  Snowflake(CriticalSection *lockToUse);

  virtual void createParameters();

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void reset() override;
  virtual void noteOn(int noteNumber, int velocity) override;

  // override set/getXml to store strings for rules and seed

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

class JUCE_API SnowflakeEditor : public AudioModuleEditor
{

public:

  SnowflakeEditor(jura::Snowflake *snowFlakeToEdit);
  virtual void resized() override;

protected:

  jura::Snowflake *snowflakeModule;

  //AutomatableSlider *frequencySlider;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SnowflakeEditor)
};


#endif