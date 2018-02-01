#ifndef jura_NewSynth_h
#define jura_NewSynth_h

class JUCE_API NewSynthAudioModule : public jura::AudioModuleWithMidiIn
{

public:

  /** Constructor. */
  NewSynthAudioModule(CriticalSection *lockToUse);


  virtual void processBlock(double **inOut, int numChannels, int numSamples) override
  {
    for(int n = 0; n < numSamples; n++)
      synthCore.getSampleFrameStereo(&inOut[0][n], &inOut[1][n]);
  }
  virtual void processStereoFrame(double *left, double *right) override
  {
    synthCore.getSampleFrameStereo(left, right);
  }

protected:

  rosic::rsNewSynth synthCore;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NewSynthAudioModule)
};

//=================================================================================================

class JUCE_API NewSynthEditor : public jura::AudioModuleEditor
{

public:

protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NewSynthEditor)
};

#endif
