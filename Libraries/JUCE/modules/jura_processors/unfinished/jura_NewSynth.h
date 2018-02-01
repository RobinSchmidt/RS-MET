#ifndef jura_NewSynth_h
#define jura_NewSynth_h

class JUCE_API NewSynthAudioModule : public jura::AudioModuleWithMidiIn
{

public:

  /** Constructor. */
  NewSynthAudioModule(CriticalSection *lockToUse);

  AudioModuleEditor* createEditor() override;

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

  // child modules:
  QuadSourceAudioModule* sourceModule;
  DualFilterAudioModule* filterModule;
  PolyModulatorsAudioModule* modulatorsModule;

  // dsp core:
  rosic::rsNewSynth synthCore;

  friend class NewSynthEditor;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NewSynthAudioModule)
};

//=================================================================================================

class JUCE_API NewSynthEditor : public jura::AudioModuleEditor
{

public:

  NewSynthEditor(CriticalSection* lockToUse, NewSynthAudioModule* synthToEdit);


  virtual void resized() override;


protected:

  // child editors:
  QuadSourceEditor*     sourceEditor;
  DualFilterEditor*     filterEditor;
  PolyModulatorsEditor* modulatorsEditor;

  // edited module
  NewSynthAudioModule* synthModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NewSynthEditor)
};

#endif
