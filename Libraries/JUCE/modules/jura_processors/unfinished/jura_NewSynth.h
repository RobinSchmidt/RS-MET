#ifndef jura_NewSynth_h
#define jura_NewSynth_h

class JUCE_API NewSynthAudioModule : public jura::AudioModuleWithMidiIn
{

public:

  /** Constructor. */
  NewSynthAudioModule(CriticalSection *lockToUse);

  AudioModuleEditor* createEditor(int type) override;

  virtual void processBlock(double **inOut, int numChannels, int numSamples) override
  {
    sourceModule->processBlock(inOut, numChannels, numSamples);
    filterModule->processBlock(inOut, numChannels, numSamples);
    // todo: apply modulations and/or smoothing

    //for(int n = 0; n < numSamples; n++)
    //  synthCore.getSampleFrameStereo(&inOut[0][n], &inOut[1][n]);
  }
  virtual void processStereoFrame(double *left, double *right) override
  {
    sourceModule->processStereoFrame(left, right);
    filterModule->processStereoFrame(left, right);
    //synthCore.getSampleFrameStereo(left, right);
  }

protected:

  /** Populates the factory objects that create the actual modules. */
  void populateModuleFactories();

  // child modules:
  QuadSourceAudioModule* sourceModule;
  DualFilterAudioModule* filterModule;
  PolyModulatorsAudioModule* modulatorsModule;

  // factory objects to create sources, filters and modulators:
  AudioModuleFactory sourceFactory, filterFactory, modulatorFactory;

   
  ModulationManager modManager; // use a PolyModulationManager


  // dsp core:
  //rosic::rsNewSynth synthCore; //...hmm...maybe get rid of this

  friend class NewSynthEditor;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NewSynthAudioModule)
};

//=================================================================================================

class JUCE_API NewSynthEditor : public jura::AudioModuleEditor
{

public:

  NewSynthEditor(CriticalSection* lockToUse, NewSynthAudioModule* synthToEdit);


  virtual void resized() override;

  // maybe have the option for a tabbed gui, showing one of the 3 editors at a time - better for 
  // smaller screens


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
