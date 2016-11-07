#ifndef jura_Enveloper_h
#define jura_Enveloper_h
  
//=================================================================================================

/** Uses the RAPT::BreakpointModulator for an amplitude modulation/enveloping effect. The main 
purpose of this class is actually to test the wrapper and editor classes for the brakpoint
modulator in a minimal context. */

class JUCE_API Enveloper : public jura::AudioModuleWithMidiIn
{

public:

  Enveloper(CriticalSection *lockToUse);

  // overriden from AudioModule baseclass:
  AudioModuleEditor *createEditor() override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  //virtual void reset() override;

  // we need to override noteOn and noteOff, too

protected:

  RAPT::rsBreakpointModulator envGen;
  jura::BreakpointModulatorAudioModule envGenWrapper; // rename to envGenModule

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Enveloper)
};

////=================================================================================================
//
///** This is the GUI editor class for the jura::Enveloper audio module. 
//\todo: maybe, we don't need this - we may let the Enveloper's createEditor method return an object
//of class jura::BreakpointModulatorEditor. */
//
//class JUCE_API EnveloperEditor : public AudioModuleEditor
//{
//
//public:
//
//  EnveloperEditor(jura::Enveloper *newEnveloperToEdit);
//  virtual void resized() override;
//
//protected:
//
//  Enveloper *enveloperToEdit;
//  BreakpointModulatorEditor *envelopeEditor;
//
//  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(EnveloperEditor)
//};

#endif 