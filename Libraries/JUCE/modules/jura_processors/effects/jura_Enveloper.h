#ifndef jura_Enveloper_h
#define jura_Enveloper_h
  
//=================================================================================================

/** Uses the RAPT::BreakpointModulator for an amplitude modulation/enveloping effect. The main 
purpose of this class is actually to test the wrapper and editor classes for the breakpoint
modulator in a minimal context. */

class JUCE_API Enveloper : public jura::AudioModuleWithMidiIn
{

public:

  Enveloper(CriticalSection *lockToUse);

  // overriden from AudioModule baseclass:
  AudioModuleEditor *createEditor() override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void noteOn(int noteNumber, int velocity) override;
  virtual void noteOff(int noteNumber) override;
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) override;
  //virtual void reset() override;

protected:

  rosic::BreakpointModulator envGen;
  //jura::BreakpointModulatorAudioModule envGenWrapper; // rename to envGenModule
  jura::BreakpointModulatorAudioModule* envGenWrapper;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Enveloper)
};

#endif 